#### Libraries ####
library(pointr)
library(parallel)
library(viridis)
library(gridExtra)
library(ggplot2)
library(jsonlite)
library(ggnewscale)
library(scales) # color ??
# library(beeswarm) # similar to dotplot
# library(ggridges) # nice density representation
# mixOmics #required

#### Common ggplot theme ####
# ggplot <- function(...) ggplot2::ggplot(...) + scale_color_gradient2(low = muted('blue'), 
#                                                                      high = muted('red')) +
#         scale_fill_gradient2(low = muted('blue'), high = muted('red'))
theme_set(ggplot2::theme_classic() + ggplot2::theme(axis.title = element_text(size  = 15, color = 'black'),
                                                    axis.text = element_text(size = 15, color = 'black'),
                                                    legend.text = element_text(size =15, color = 'black'),
                                                    legend.title = element_text(size = 15, color = 'black')))

Path2File <- normalizePath(paste0(gsub('generic.R', '', paste0(sys.frames()[[1]]$ofile)), '/../data/hex.csv'))
if(!file.exists(Path2File)){
        warning("File hex.csv does not exist. Trying to download from repository.")
        download.file("https://github.com/Polfe94/AnTracks/tree/main/data/hex.csv",
                      Path2File)
} 
hex <- read.csv(Path2File)
rm(Path2File)

#### Functions +++ Computation ####

#' Initializes data object to a node-matrix, simulation or coordinates class
#' 
#' @param obj A data.frame containing data about X, Y coordinates, time, and others
#' @param refcoords A data.frame containing the coordinates of the vertices of the hexagonal lattice
#' @param class A character, one of: "nodes", "simulation", "coords". Class to convert data to.
#' @param ... Additional arguments to be passed, such as already computed segments, experiment date...
#' @TO.BE.ADDED: AVAILABLE ARGUMENTS (date, segments, )
#' @return Object of the class specified in \code{class}
init <- function(obj, refcoords, class, ...){
        if(!class %in% c('nodes', 'simulation', 'coords')){
                stop('only "nodes", "simulation" or "coords" are currently available classes')
        }
        tmp <- list(data = as.data.frame(obj), refcoords = refcoords)
        args <- list(...)
        if(length(args) && length(names(args)) && length(args) == length(names(args))){
                n <- names(args)
                for(i in seq_along(args)){
                        tmp[[n[i]]] <- args[[i]]
                }
        }
        
        # assign radius to compute neighbours
        class(tmp) <- class
        if(class == 'nodes' || class == 'coords'){
                r <- 51
        } else {
                r <- 1.01 
        }
        tmp$r <- r
        tmp
}

#' Computes the pairwise euclidean distance between two sets of coordinates
#' 
#' @param A A matrix of dim(N, 2)
#' @param B A matrix of dim(N, 2)
#' @param ret.vec A boolean indicating whether or not a vector should be returned
#' @return A numeric vector (or matrix, if ret.vec = FALSE) of distances (in whatever unit)
#' with length = nrow(A)*nrow(B) (or matrix with rows = nrow(A) and columns = nrow(B))
pdist <- function(A, B, ret.vec = TRUE){
     an = apply(A, 1, function(rvec) crossprod(rvec,rvec))
     bn = apply(B, 1, function(rvec) crossprod(rvec,rvec))
     
     m = nrow(A)
     n = nrow(B)
     
     tmp = matrix(rep(an, n), nrow=m) 
     tmp = tmp +  matrix(rep(bn, m), nrow=m, byrow=TRUE)
     
     d <- rbind(sqrt(tmp - 2 * tcrossprod(A,B)))
     
     if(ret.vec == TRUE){
          as.numeric(d)
     } else {
          d
     }
}

get_node <- function(x, y, ref){
        xy <- t(c(x, y))
        which.min(pdist(xy, as.matrix(ref)))
}

#' Normalizes a numeric vector to fit in a given range
#' 
#' @param x A numeric vector to be normalized
#' @param a The lower boundary of the range
#' @param b The upper boundary of the range
#' @return The normalized version of x
norm_range <- function(x, a = -1, b = 1){
        (b - a) * (x - min(x)) / (max(x) - min(x)) + a
}

#' Rotates the coordinates to a certain angle
#' 
#' @param x A numeric (x) coordinate
#' @param y A numeric (y) coordinate
#' @param theta A number representing an angle
#' @return A data.frame with the rotated coordinates
rotate.theta <- function(x, y, theta = pi/2){
     
     if(length(x) != length(y)){
          stop("x and y need to be the same length!")
     }
     
     m <- t(vapply(seq_along(x), function(i){
          x1 <- round(x[i]*cos(theta)-y[i]*sin(theta), 2)
          y1 <- round(x[i]*sin(theta)+y[i]*cos(theta), 2)
          c(x = x1, y = y1)
     }, numeric(2)))
     
     as.data.frame(m)
}

get_nest <- function(region = 'top'){
        if(region == 'top'){
                c(1020, 1015)
        } else if (grepl('bot', region)){
                c(990, 990)
        } else {
                rbind(c(1020, 1015), c(990, 990))
        }
}

closest_node <- function(obj){
        if(!'coords' %in% class(obj)){
                return(obj)
        }

        xy <- cbind(obj$data$Xmm, obj$data$Ymm)
        # h <- as.matrix(obj$refcoords)
        idx <- obj$refcoords[, 'y'] > 1000
        n <- obj$refcoords$node[idx]
        h <- as.matrix(obj$refcoords[idx, c('x', 'y')])
        # idx <- vapply(seq_len(nrow(xy)), function(i){
        #         which.min(pdist(h, t(xy[i, ])))
        # }, integer(1))
        idx <- vapply(seq_len(nrow(xy)), function(i){
                n[which.min(pdist(h, t(xy[i, ])))]
        }, integer(1))
        obj$data$node <- idx
        obj
}

#' Checks if a list of points is inside a radius
#' 
#' @param coords A numeric of length 2 or data frame with 2 columns (x, y)
#' @param center A numeric indicating  center of the circle
#' @param r A number indicating the size of the radius (in whatever unit)
#' @return A vector of TRUE / FALSE of length equal to coords
inRadius <- function(coords, center, r){
        if(length(center) != 2){
                stop('center needs to be an x, y coordinate')
        }
        
        if(is.numeric(coords) && length(coords) == 2){
                warning('coords is not a data.frame; treating coords as a single x, y coordinate')
                coords <- t(coords)
        }
        ((coords[, 1] - center[1]) ^ 2 + (coords[, 2] - center[2]) ^ 2) < (r ^ 2)

}

#' Computes the "origin - destination" (neighbouring) segments in the hexagonal lattice
#' 
#' @param obj An object of whatever class (coords, lattice or simulations)
#' @return A data.frame with the origin coordinates and the destination (neighbouring) coordinates
compute_segments <- function(obj){
     refcoords <- obj$refcoords
     r <- obj$r
     # xy <- as.matrix(refcoords[, c('x', 'y')])
     xy <- as.matrix(refcoords[refcoords[, 'y'] > 1000, c('x', 'y')])
     df <- c()

     for(i in 1:nrow(xy)){
          idx <- pdist(t(xy[i, ]), xy) < r
          idx[i] <- FALSE
          
          df <- rbind(df, data.frame(x = rep(xy[i, 1], sum(idx)),
                                     y = rep(xy[i, 2], sum(idx)),
                                     xend = xy[idx, 1], yend = xy[idx, 2],
                                     o = rep(rownames(xy)[i], sum(idx)),
                                     d = rownames(xy)[idx]))
     }
     df
}

aggregate_time <- function(t, tau = 10){
     sq <- seq(min(t), max(t), tau)
     
     idx <- lapply(2:length(sq), function(i){
          which(t >= sq[i-1] & t < sq[i])
     })
     
     names(idx) <- sq[-1]
     
     idx <- idx[lapply(idx, length) > 0]
     
     idx
}

interaction_matrix <- function(obj){
        if(!'coords' %in% class(obj)){
                stop('Object must be of class "coords"')
        }
        if(!'node' %in% colnames(obj$data)){
                stop('Object must have "node" computed')
        }
        n <- obj$data$node
        t <- obj$data$Frame[obj$data$Crossings > 0]
        m <- matrix(data = 0, ncol = length(unique(n)),
                    nrow = max(obj$data$Frame))
        dimnames(m) <- list(seq_len(nrow(m)), unique(n))
        
        for(i in seq_along(t)){
                idx <- which(colnames(m) == n[i])
                m[t[i], idx] <- m[t[i], idx] + 1
                # m[t[i], idx] <- 1
        }
        
        m
}

interactions_in_trail <- function(obj, ...){
        if(!'trails' %in% names(obj)){
                obj$trails <- food_trails(obj, ...)
        }
        if(!'IiT' %in% names(obj)){
                m <- interaction_matrix(obj)
                df<- data.frame(t = as.integer(rownames(m)), 
                                N = rowSums(m),
                                inTrail = rowSums(m[, as.character(obj$trails)]),
                                outTrail = rowSums(m[, !colnames(m) %in% as.character(obj$trails)]))
                mlt_df <- reshape2::melt(df, id.vars = 't')
        } else {
                mlt_df <- obj$IiT
        }
        mlt_df
}


activity_in_trail <- function(obj, ...){
        if(!'trails' %in% names(obj)){
                obj$trails <- food_trails(obj, ...)
        }
        if(!'AiT' %in% names(obj)){
                m <- coords2matrix(obj, data = 0)
                df<- data.frame(t = as.integer(rownames(m)), 
                                N = rowSums(m),
                                inTrail = rowSums(m[, as.character(obj$trails)]),
                                outTrail = rowSums(m[, !colnames(m) %in% as.character(obj$trails)]))
                mlt_df <- reshape2::melt(df, id.vars = 't')
        } else {
                mlt_df <- obj$AiT
        }
        # ptr('tmp', obj$food[1:2])
        # name <- deparse(substitute(obj))
        # if(grepl('[[i]]', name)){
        #         name <- gsub('i', eval(i), name)
        # }
        # name <- paste0(name, '["food"][1:2]')
        # print(name)
        # print(as.character(name))
        # x <- list(data = mlt_df, food = parse(text = paste0("ptr('tmp', ", name, ")")))
        mlt_df
        # x <- list(data = mlt_df, food = obj$food)
        # class(x) <- 'AiT'
        # x
}

transform_nodes <- function(x){
        r <- numeric(length(x))
        idx <- 1
        for(i in x){
                r[idx] <- hex$node[hex$old_node == i]
                idx <- idx + 1
        }
        r
}

food_trails <- function(obj, method = 'post-exploration', prob = 0.90){
        
        if(!'food' %in% names(obj)){
                stop('Compute food times before searching for food trails!')
        } 
        if(method == 'in-recruitment'){
                intervals <- range(do.call('rbind', obj$food)$t)
        } else if(method == 'post-recruitment'){
                intervals <- c(max(do.call('rbind', obj$food)$t), max(obj$data$Frame))
        } else if(method == 'post-exploration'){
                intervals <- c(min(do.call('rbind', obj$food)$t), max(obj$data$Frame))
        } else if(method == 'mid-recruitment') {
                intervals <- c(as.integer(mean(do.call('rbind', obj$food)$t)), max(obj$data$Frame))
        } else {
                stop('Currently available methods are "in-recruitment" and "post-recruitment"')
        }
        
        m <- coords2matrix(obj, data = 0)
        m <- m[intervals[1]:intervals[2], ]
        
        x <- colSums(m)
        x[x < quantile(x, prob = prob)] <- 0
        y <- as.integer(names(x))
        if(1241 %in% y){
                y <- transform_nodes(y)
                names(x) <- y
        }
        y[x > 0]
}


cluster_lengths <- function(x, tau = 10){
     
     idx <- aggregate_time(x$t, tau = tau)
     output <- vapply(idx, function(i){
          result <- connectivity(x[i, c('idx')])
     }, numeric(1))
     
     data.frame(t = c(0, as.numeric(names(idx))),
                k = c(0, output))
} 

#### Functions +++ Visualization ####

#' Draws the hexagonal lattice
#' 
#' @param obj An object of whatever class (coords, lattice, simulation)
#' @param add Either NULL or a ggplot object
#' @param z Either NULL or a numeric vector with a Z value to produce a heatmap
#' @param ... Additional parameters to be passed to ggplot (color, size, alpha ...)
#' @return A ggplot object drawing the hexagonal lattice on top of whatever (if any)
#' layer is passed to the argument add. The colour of the hexagons will be based on
#' the values provided in \code{z} (if any, black otherwise).
draw_hexagons <- function(obj, add = NULL, z = NULL, ...){
     if(!'segments' %in% names(obj)){
          obj$segments <- compute_segments(obj)
     }
     
     if(!is.null(add)){
          pl <- add
     } else {
          pl <- ggplot()
     }
     
     if(!is.null(z)){
          pl <- pl + geom_segment(data = obj$segments, 
                                  aes(x = x, xend = xend, y = y, yend = yend, color = z), ...)
     } else {
          pl <- pl + geom_segment(data = obj$segments, aes(x = x, xend = xend, y = y, yend = yend), ...)
     }
     pl + xlab('') + ylab('') + theme(axis.text = element_blank(), axis.ticks = element_blank(),
                                        axis.line = element_blank())
}

#' Draws the food patches
#' 
#' @param obj An object of whatever class (coords, lattice, simulation)
#' @param add Either NULL or a ggplot object
#' @param ... Additional parameters to be passed to ggplot (color, size, alpha ...)
#' @return A ggplot object drawing the food patches on top of whatever (if any)
#' layer is passed to the argument add.
draw_FoodPatches <- function(obj, add = NULL, ...){
        if(is.null(add)){
                add <- ggplot()
        }
        
        if(all('coords' %in% sapply(obj, class))){
                for(x in seq_along(obj)){
                        tmp <- obj[[x]]
                        for(i in seq_along(tmp$food)){
                                add <- add + geom_polygon(data = tmp$food[[i]], aes(x, y), ...) +
                                        xlab('') + ylab('')
                        }
                }
        } else {
                for(i in seq_along(obj$food)){
                        add <- add + geom_polygon(data = obj$food[[i]], aes(x, y), ...) +
                                xlab('') + ylab('')
                }
        }
        add
}

plot.IiT <- function(obj){
        df <- obj$data
        
        ylim <- c(0, 0, rep(max(df$value), 2))
        xlim <- range(do.call('rbind', obj$food)$t)[c(1:2, 2:1)]/120
        
        ggplot(data = df, aes(t/120, value, color = variable))+
                geom_line() + scale_y_continuous('', breaks = seq(0, max(ylim), 0.5)) +
                scale_color_viridis_d('', labels = c('Total Interactions', 
                                                     'Interactions in trail', 'Interactions out trail'))+
                scale_x_continuous('Time (min)', breaks = seq(0, 180, 15))+
                geom_polygon(data = data.frame(x = xlim,
                                               y = ylim), aes(x, y),
                             fill = muted('green'), color = 'green', alpha = 0.15, linetype = 2, size = 1)+
                guides(color = guide_legend(override.aes = list(size = 2)))

        
}

plot_trails <- function(obj){
        draw_hexagons(obj, add = draw_FoodPatches(obj, fill = 'grey50') +
                              geom_point(data = data.frame(x = hex$x[obj$trails],
                                                           y = hex$y[obj$trails]),
                                         aes(x, y), fill = muted('blue'), size = 4, shape = 21))
}

plot_IiT <- function(obj){
        df <- obj$IiT
        tmp <- obj$food
        # eval(obj$food)
        
        
        ylim <- c(0, 0, rep(max(df$value), 2))
        xlim <- range(do.call('rbind', tmp)$t)[c(1:2, 2:1)]/120
        
        ggplot(data = df, aes(t/120, value, color = variable))+
                geom_line() + scale_y_continuous('', breaks = seq(0, max(ylim), 1)) +
                scale_color_viridis_d('', labels = c('Total activity', 'Activity in trail', 'Activity out trail'))+
                scale_x_continuous('Time (min)', breaks = seq(0, 180, 15))+
                geom_polygon(data = data.frame(x = xlim,
                                               y = ylim), aes(x, y),
                             fill = muted('green'), color = 'green', alpha = 0.15, linetype = 2, size = 1)+
                guides(color = guide_legend(override.aes = list(size = 2)))
        
}

plot_interactions <- function(obj){
        m <- interaction_matrix(obj)
        nodes <- transform_nodes(as.integer(colnames(m)))
        y <- colSums(m)
        idx <- y > 0
        draw_hexagons(obj, add = draw_FoodPatches(obj, fill = 'grey50') +
                              geom_point(data = data.frame(x = hex$x[nodes[idx]],
                                                           y = hex$y[nodes[idx]],
                                                           z = y[idx]),
                                         aes(x, y, fill = z), size = 4, shape = 21, show.legend = F))+
                geom_point(data = as.data.frame(matrix(get_nest(), ncol = 2, dimnames = list(1, c('x', 'y')))),
                           aes(x, y), shape = 17, size = 2.5)
}

plot_AiT <- function(obj){
        df <- obj$AiT
        tmp <- obj$food
        # eval(obj$food)
        
        
        ylim <- c(0, 0, rep(max(df$value), 2))
        xlim <- range(do.call('rbind', tmp)$t)[c(1:2, 2:1)]/120

        ggplot(data = df, aes(t/120, value, color = variable))+
                geom_line() + scale_y_continuous('', breaks = seq(0, max(ylim), 5)) +
                scale_color_viridis_d('', labels = c('Total activity', 'Activity in trail', 'Activity out trail'))+
                scale_x_continuous('Time (min)', breaks = seq(0, 180, 15))+
                geom_polygon(data = data.frame(x = xlim,
                                               y = ylim), aes(x, y),
                             fill = muted('green'), color = 'green', alpha = 0.15, linetype = 2, size = 1)+
                guides(color = guide_legend(override.aes = list(size = 2)))
        
}

geom_circle <- function(center, r, npoints = 100, ...){
        sq <- seq(0, 2*pi, length.out = npoints)
        if(is.data.frame(center)){
                
                if(nrow(center) == 1){
                        x <- center[, 1] + r * cos(sq)
                        y <- center[, 2] + r * sin(sq)
                        geom_path(data = data.frame(x = x, y = y), aes(x, y), ...)
                } else {
                        x <- c()
                        y <- c()
                        g <- c()
                        for(i in seq_len(nrow(center))){
                                x <- c(x, center[i, 1] + r * cos(sq))
                                y <- c(y, center[i, 2] + r * sin(sq))
                                g <- c(g, rep(i, npoints))
                        }
                        geom_path(data = data.frame(x = x, y = y, g = g), aes(x, y, group = factor(g)), ...)
                }
                
        } else if (is.atomic(center) && length(center) == 2){
                
                geom_path(data = data.frame(x = center[1] + r * cos(sq),
                                            y = center[2] + r * sin(sq)),
                          aes(x, y, group = factor(g)), ...)
        } else {
                stop('Please use a data.frame or an atomic vector of length to for the center')
        }
        
}

#### Methods ####
.plot <- function(obj, type){
        UseMethod('plot', type)
}

get_N <- function(obj){
     UseMethod('get_N')
}

get_I <- function(obj){
        UseMethod('get_I')
}

connectivity <- function(obj, tau){
     UseMethod('connectivity')
}

node_idx <- function(obj, row){
        UseMethod('node_idx')
}

local_cov <- function(obj){
        UseMethod('local_cov')
}

pairwise_cov <- function(obj, t){
        UseMethod('pairwise_cov')
}

mutual_info <- function(obj, t, nodes, subset){
        UseMethod('mutual_info')
}

get_neighbors <- function(obj, row){
        UseMethod('get_neighbors')
}

foodpatches <- function(obj){
        UseMethod('foodpatches')
}

get_foodPatches <- function(obj){
        UseMethod('get_foodPatches')
}

food_detection <- function(obj, r){
        UseMethod('food_detection')
}

coords2matrix <- function(obj, data){
        UseMethod('coords2matrix')
}

# .foodTrails <- function(obj, prob){
#         UseMethod('food_trails')
# }

#### +++ Extract meta-information from experiments +++ ####

#' Extracts additional information from experiments (such as patch location)
#' 
#' @param spreadsheet A data.frame containing the experiment information
#' @param expdate A character containing the experiment date
#' @return A list of length two, with a data.frame containing the
#' experiment information in the specified dates, both for green and red colony
get_metainfo <- function(expdate, spreadsheet = exp_spreadsheet){

        exptime <- vapply(expdate, function(i){
                ifelse(grepl('M', i), 'M', 'T')
        }, character(1), USE.NAMES = FALSE)

        d <- format(strptime(expdate, '%Y%m%d'), format = '%d/%m/%Y')
        
        idx <- vapply(seq_along(expdate), function(i){
                which(spreadsheet$Date == d[i] & spreadsheet$MATI.TARDA == exptime[i])
        }, integer(1))

        output <- list(G = spreadsheet[idx, c(1:10, 18:24, 35:46)], # Green colony
                       R = spreadsheet[idx, c(1:17, 26:34, 44:46)]) # Red colony
        output
}

#' Extracts the food location from a metainfo object
#' 
#' @param metainfo A data.frame with the metainfo of a single experiment
#' @return A data.frame with the food patch x and y coordinates
get_food_from_metainfo <- function(metainfo){
        
        patch2df <- function(patch_coord){
                s <- gsub('[()]', '', patch_coord)
                s <- strsplit(s, ',')
                data.frame(x = as.integer(s[[1]][1]), 
                           y = as.integer(s[[1]][2]))
        }
        
        df <- data.frame(colony = c('G', 'G', 'R', 'R'),
                         patch = c(1, 2, 1, 2))
        
        m = list(GP1 = metainfo$G$P1cord..x.y..2,
                 GP2 = metainfo$G$P2cord,
                 RP1 = metainfo$R$P1cord..x.y.,
                 RP2 = metainfo$R$P1cord..x.y..1)
        
        cbind(df, do.call('rbind', 
                          lapply(m, patch2df)))
        
}

foodpatches.top <- function(obj){
        # hexagon dimension
        dy <- 50
        dx <- 86.62
        ref <- 80 # node index top left
        
        start <- as.numeric(hex[ref, c('x', 'y')])
        start[2] <- start[2] - dy
        pos <- start

        move <- c(-dx / 2, -1.5*dy)
        
        for(i in seq_len(obj[1] -1)){
                pos <- pos + move
                move[1] <- move[1] * -1
        }
        
        for(i in seq_len(obj[2]-1)){
                pos[1] <- pos[1] + dx
        }
        
        pos
}

foodpatches.bot <- function(obj){
        # hexagon dimension
        dy <- 50
        dx <- 86.62
        ref <- 1161 # node index bottom right
        
        start <- as.numeric(hex[ref, c('x', 'y')])
        start[2] <- start[2] + dy
        pos <- start
        
        move <- c(dx / 2, 1.5*dy)
        
        for(i in seq_len(obj[1] -1)){
                pos <- pos + move
                move[1] <- move[1] * -1
        }
        
        for(i in seq_len(obj[2]-1)){
                pos[1] <- pos[1] - dx
        }
        
        pos
}

set_foodpatches <- function(patches){
        
        p <- lapply(seq_len(nrow(patches)), function(i){
                x <- patches[i, ]
                c <- ifelse(grepl('G', x$colony), 'top', 'bot')
                x <- c(x$x, x$y)
                class(x) <- c
                foodpatches(x)
        })
        
        p <- lapply(p, function(i){
                k <- inRadius(hex, i, r = 51)
                h <- hex[k, ]
                h[chull(h), ]
        })
        
        names(p) <- paste(patches$colony, 'P', patches$patch, sep = '')
        p
}
