node_idx.nodes <- function(obj, row){
        as.integer(rownames(obj$data[row, ]))
}

food_trails.nodes <- function(obj, method = 'post-recruitment', prob = 0.85){
        if(method == 'in-recruitment'){
                obj$intervals <- range(do.call('rbind', obj$food)$t)
        } else if(method == 'post-recruitment'){
                obj$intervals <- c(max(do.call('rbind', obj$food)$t), max(obj$data$Frame))
        } else {
                stop('Currently available methods are "in-recruitment" and "post-recruitment"')
        }
        
        m <- obj$data
        m <- m[obj$intervals[1]:obj$intervals[2], ]
        x <- colSums(m)
        x[x < quantile(x, prob = prob)] <- 0
        y <- as.integer(names(x))
        if(1241 %in% y){
                y <- transform_nodes(y)
                names(x) <- y
        }
        y[x > 0]
}

mutual_info.nodes <- function(obj, t, nodes = NULL, subset = FALSE){
     m <- obj$data[t, ]
     
     z <- numeric(length(det[[1]]$segments$o))
     s <- det[[1]]$segments[, c('o', 'd')]
     
     existing_segments <- s[s[, 'o'] %in% colnames(m) & s[, 'd'] %in% colnames(m), ]
         
     if(!is.null(nodes)){
             existing_segments <- existing_segments[existing_segments[, 'o'] %in% nodes &
                                                            existing_segments[, 'd'] %in% nodes, ]
     } 
     
     indices <- c()
     
     for(i in 1:nrow(existing_segments)){
          
          idx <- which(s$o == existing_segments[i, 1] & s$d == existing_segments[i, 2])
          z[idx] <- infotheo::mutinformation(m[, existing_segments[i, 1]], m[, existing_segments[i, 2]])
          indices <- c(indices, idx)
          
     }
     if(subset){
             z <- z[indices]
     }
     z
}
