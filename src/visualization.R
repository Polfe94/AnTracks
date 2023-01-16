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
draw_hexagons <- function(edges = NULL, add = NULL, z = NULL, ...){
        if(is.null(edges)){
                edges <- compute_edges(hex[hex$y > 1000, ])
        }

     if(!is.null(add)){
          pl <- add
     } else {
          pl <- ggplot()
     }
     
     if(!is.null(z)){
          pl <- pl + geom_segment(data = edges, 
                                  aes(x = x, xend = xend, y = y, yend = yend, color = z), ...)
     } else {
          pl <- pl + geom_segment(data = edges, aes(x = x, xend = xend, y = y, yend = yend), ...)
     }
     pl + xlab('') + ylab('') + theme(axis.text = element_blank(), axis.ticks = element_blank(),
                                      axis.line = element_blank())
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

geom_food <- function(t,  
                      rectangle_params = list(fill = muted('green'),
                                              color = 'green', 
                                              alpha = 0.15, linetype = 2, size = 1),
                      complete = FALSE, ylim = c(-Inf, Inf)){
        
        ylim <- sort(ylim)
        trange <- range(t)
        
        rectangle_params$data <- data.frame(x = c(trange[1], trange[1], trange[2], trange[2]),
                                            y = c(ylim[1], ylim[2], ylim[2], ylim[1]))
        rectangle_params$mapping <- aes(x, y)
        
        rectangle_geom <- do.call('geom_polygon', args = rectangle_params)
        
        if(complete){
                segment_params <- rectangle_params
                segment_params$fill <- NULL
                segment_params$alpha <- 1
                segment_params$size <- 0.5
                segment_params$data <- data.frame(x = t[t > trange[1] & t < trange[2]],
                                                  y = ylim[1], yend = ylim[2])
                
                segment_params$mapping <- aes(x = x, xend = x, y = y, yend = yend)
                
                food_geom <- do.call('geom_segment', args = segment_params)
                
                return(list(rectangle_geom, food_geom))
        }
        
        rectangle_geom
}
