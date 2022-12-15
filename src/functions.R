## --- LIBRARIES --- ##
library(stringr)
library(viridis)
library(gridExtra)
library(ggpubr)
library(data.table)
library(reshape2)
library(ggplot2)
library(jsonlite)
library(ggnewscale)
library(scales) 


## --- SET GGPLOT THEME --- ##
theme_set(ggplot2::theme_classic() + ggplot2::theme(axis.title = element_text(size  = 15, color = 'black'),
                                                    axis.text = element_text(size = 15, color = 'black'),
                                                    legend.text = element_text(size =15, color = 'black'),
                                                    legend.title = element_text(size = 15, color = 'black')))

#### +++ GENERIC FUNCTIONS +++ ####
compute_edges <- function(refcoords = hex[hex$y > 1000, ], r = 51){
     xy <- as.matrix(refcoords[, c('x', 'y')])
     
     d <- pdist(xy, xy)
     idx <- which(d < r & d > 0, arr.ind = TRUE)
     
     edges <- data.frame(x = xy[idx[, 1], 1], y = xy[idx[, 1], 2],
                                xend = xy[idx[, 2], 1], yend = xy[idx[, 2], 2],
                                o = refcoords$node[idx[, 1]], d = refcoords$node[idx[, 2]])
     data.table(edges)
}

#' Returns the closest node to the provided coordinate
#' 
#' @param ... one of [1] 'data.frame' with 2 dimensions, [2] 'numeric' of length 2, [3] two 'numeric' of length 1
#' @return numeric vector of length equal to the number of coordinates passed to the function
get_node <- function(...){
     l <- list(...)
     if('xy' %in% names(l)){
             xy <- as.matrix(l$xy[, c('x', 'y')])
             n <- l$xy$node
             l$xy <- NULL
     } else {
             xy <- as.matrix(hex[, c('x', 'y')])
             n <- hex$node
     }

     if(length(l) == 2){
          idx <- which.min(pdist(t(do.call('c', l)), xy))
     } else if(length(l) == 1){
          x <- l[[1]]
          if(is.atomic(x)){
               idx <- which.min(pdist(t(x), xy))
          } else {
               if(dim(x)[2] == 2){
                       # idx <- apply(fastPdist2(as.matrix(x), xy), 1, which.min) 
                    idx <- apply(pdist(as.matrix(x), xy), 1, which.min)
               } else {
                    stop('Object must have X and Y coordinates')
               }
          }
          
     } else {
          stop('Could not read coordinates')
     }
     n[idx]
}

#' Computes the pairwise euclidean distance between two sets of coordinates
#' 
#' @param A A matrix of dim(N, 2)
#' @param B A matrix of dim(N, 2)
#' @return A matrix of distances of dimensions = nrow(A) X nrow(B)
#' with length = nrow(A)*nrow(B) (or matrix with rows = nrow(A) and columns = nrow(B))
Rcpp::cppFunction(
        '
NumericMatrix pdist(NumericMatrix A, NumericMatrix B) {
     int m = A.nrow(), 
          n = B.nrow(),
          k = A.ncol();
     arma::mat Ar = arma::mat(A.begin(), m, k, false); 
     arma::mat Br = arma::mat(B.begin(), n, k, false); 
     
     arma::colvec An =  sum(square(Ar),1);
     arma::colvec Bn =  sum(square(Br),1);
     
     arma::mat C = -2 * (Ar * Br.t());
     C.each_col() += An;
     C.each_row() += Bn.t();
     
     return wrap(sqrt(C)); 
}', c('Rcpp', 'RcppArmadillo'))

nest_influence <- function(r = 101){
        idx <- inRadius(hex[, c('x', 'y')], as.numeric(hex[663, c('x', 'y')]), r)
        hex$node[idx]
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
rotate <- function(x, y, theta = pi/2){
     
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

foodpatches <- function(obj){
     UseMethod('foodpatches')
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
          k <- inRadius(hex[, c('x', 'y')], i, r = 51)
          h <- hex[k, c('x', 'y')]
          h[chull(h), ]
     })
     
     names(p) <- paste(patches$colony, 'P', patches$patch, sep = '')
     p
}

get_foodPatches <- function(date){
     if(is.character(date) && str_count(date, '[0-9]') == 8 && str_count(date, '[M, T]') == 1){
          m <- get_metainfo(date)
          p <- get_food_from_metainfo(m)
          p <- set_foodpatches(p)
          p
     } else {
          stop('Invalid experiment date')
     }
}

#' Transforms a simple matrix M(Frame, node) to a filled array A(Frame, node, value)
#' 
#' @param DT a data.table containing (ONLY) a Frame and a node column
#' @param data a numeric to complete missing cases with
#' @return a three column matrix
M2M <- function(DT, data, maxt = NULL){
        if(is.null(maxt)){
                maxt <- max(DT$Frame)
        }
        if(is.data.frame(DT) && !is.data.table(DT)){
                DT <- data.table(DT)
        }
     M <- DT[CJ(Frame = seq_len(maxt), node = node, unique = T), on=. (Frame, node)]
     M$N[is.na(M$N)] <- data
     M
}

#' Helper function: transforms a filled array (see function @code M2M) to a subsetted matrix
#' 
#' @param M a data.table as produced by function @code M2M
#' @param t a numeric of length n specifying the frames to subset the matrix
#' @return a subsetted (and decasted) matrix of M(i, j) with i = length(t) and j = number of nodes
M2sbst <- function(M, t){
     sbst <- M[Frame %in% t, ]
     sbst <- data.table::dcast(sbst, Frame ~ node, value.var = 'N')
     sbst$Frame <- NULL
     rownames(sbst) <- t
     sbst
}

get_neighbors <- function(nodes, ...){
        l <- list(...)
        if('edges' %in% names(l)){
                edges <- l$edges
        } 
        if(!exists('edges', envir = .GlobalEnv)){
                if('refcoords' %in% names(l)){
                        edges <<- compute_edges(refcoords)
                } else {
                        edges <<- compute_edges()
                }
        }
        as.integer(edges$d[edges$o %in% nodes])
}

optimal_path <- function(start, target, refcoords = hex[hex$y > 1000, ], r = 51){
        
        edges <- compute_edges(refcoords, r)
        
        x1 <- as.matrix(refcoords[refcoords$node == start, c('x', 'y')])
        x2 <- as.matrix(refcoords[refcoords$node == target, c('x', 'y')])
        
        path <- start
        pos <- start
        
        while(pos != target){
                
                available_positions <- get_neighbors(pos, edges = edges)
                distances <- pdist(as.matrix(refcoords[refcoords$node %in% available_positions, c('x', 'y')]), x2)
                idx <- which.min(distances)
                pos <- available_positions[idx]
                path <- c(path , pos)
                
        }
        unique(path)
}

possible_paths <- function(start, target, refcoords = hex[hex$y > 1000, ], r = 51){
        
        edges <- compute_edges(refcoords, r)
        
        x1 <- as.matrix(refcoords[refcoords$node == start, c('x', 'y')])
        x2 <- as.matrix(refcoords[refcoords$node == target, c('x', 'y')])
        
        d <- pdist(x1, x2)
        
        paths <- data.table(node = start, t = 1)
        
        current_origin <- start
        current_distances <- d
        future_origins <- c()
        future_distances <- c()
        
        i <- 2
        iters <- 0
        
        while(TRUE){
                
                for(pos in seq_along(current_origin)){
                        available_positions <- get_neighbors(current_origin[pos], edges = edges)
                        distances <- pdist(as.matrix(refcoords[refcoords$node %in% available_positions, c('x', 'y')]), x2)
                        idx <- which(distances < current_distances[pos])
                        
                        if(length(idx)){
                                future_origins <- c(future_origins, available_positions[idx])
                                future_distances <- c(future_distances, distances[idx])
                        }
                        
                }
                if(length(future_origins)){
                        paths <- rbind(paths, data.table(node = future_origins, t = i))
                        current_origin <- future_origins
                        current_distances <- future_distances
                        future_origins <- c()
                        future_distances <- c()
                        i <- i + 1
                } else {
                        break
                }
        }
        unique(paths$node)
}

moving_average <- function(x, t, overlap = 0){
        
        slide <- t %/% 2
        if(overlap > slide){
                warning('Overlap is bigger than window size. Capping overlap to window size.')
                overlap <- slide
        }
        x0 <- 1 + slide
        xn <- length(x) - slide
        
        sq <- seq(x0, xn, (slide - overlap + 1))
        v <- vapply(sq, function(i){
                mean(x[seq(i-slide, i+slide)])
        }, numeric(1))
        
        # fill in start and finish positions
        if(overlap == slide){
                c(x[1:slide], v, x[(xn+1):length(x)])
        } else {
                v
        }

}
