node_idx.nodes <- function(obj, row){
        as.integer(rownames(obj$data[row, ]))
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
