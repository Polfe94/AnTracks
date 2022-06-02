#### Libraries ####
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
ggplot <- function(...) ggplot2::ggplot(...) + scale_color_gradient2(low = muted('blue'), 
                                                                     high = muted('red')) +
        scale_fill_gradient2(low = muted('blue'), high = muted('red'))
theme_set(ggplot2::theme_classic() + ggplot2::theme(axis.title = element_text(size  = 15, color = 'black'),
                                                    axis.text = element_text(size = 15, color = 'black'),
                                                    legend.text = element_text(size =15, color = 'black'),
                                                    legend.title = element_text(size = 15, color = 'black')))


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

closest_node <- function(obj){
        if(!'coords' %in% class(obj)){
                return(obj)
        }

        xy <- cbind(obj$data$Xmm, obj$data$Ymm)
        h <- as.matrix(obj$refcoords)
        idx <- vapply(seq_len(nrow(xy)), function(i){
                which.min(pdist(h, t(xy[i, ])))
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
     xy <- as.matrix(refcoords[, c('x', 'y')])
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

#' Transforms an object of class coords to an object of class nodes
#' 
#' @param obj An object of class coords
#' @return An object of class nodes
coords2matrix <- function(obj){
        
        if(!'coords' %in% class(obj)){
                stop('Object must be of class "coords"')
        }
        if(!'node' %in% colnames(obj$data)){
                stop('Object must have "node" computed')
        }
        n <- obj$data$node
        t <- obj$data$Frame
        m <- matrix(data = -1, ncol = length(unique(n)), nrow = max(t))
        dimnames(m) <- list(seq_len(nrow(m)), unique(n))
        
        for(i in seq_along(t)){
                idx <- which(colnames(m) == n[i])
                m[t[i], idx] <- 1
        }
        class(m) <- 'nodes'
        m
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

#### Methods ####
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

mutual_info <- function(obj, t){
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
        
        start <- as.numeric(hex[ref, ])
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
        
        start <- as.numeric(hex[ref, ])
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
