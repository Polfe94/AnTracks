## Class setter
setClass('simulation', slots = list(data = 'data.frame',
                                    params = 'list',
                                    refcoords = 'data.frame'))

## Class initializer
init.simulation <- function(obj, ...){
     tmp <- list(data = obj)
     args <- list(...)
     if(length(args) && length(names(args)) && length(args) == length(names(args))){
          n <- names(args)
          for(i in seq_along(args)){
               tmp[[n[i]]] <- args[[i]]
          }
     }
     tmp$r <- 1.01
     tmp
}

connectivity <- function(pos){
     
     if(length(pos)>0){
          pos <- unique(pos)
          k_length <- c()
          branch <- c()
          
          while(length(pos)){
               
               current_path <- pos[length(pos)]
               pos <- pos[-length(pos)]
               
               while(length(current_path)){
                    
                    target <- current_path[length(current_path)]
                    neighbors <- get_neighbors(target, hex[, 5:6], r = 1.1)
                    idx <- which(neighbors %in% pos)
                    branch <- c(branch, neighbors[idx])
                    while(length(branch)){
                         
                         current_path <- c(current_path, branch[length(branch)])
                         branch <- branch[-length(branch)]
                         if(current_path[length(current_path)] %in% pos){
                              pos <- pos[-which(pos == current_path[length(current_path)])]
                         } else {
                              next
                         }
                         target <- current_path[length(current_path)]
                         neighbors <- get_neighbors(target, hex[, 5:6], r = 1.1)
                         idx <- which(neighbors %in% pos)
                         branch <- c(branch, neighbors[idx])
                    }
                    k_length <- c(k_length, length(current_path))
                    current_path <- c()
                    branch <- c()
               }
          }
          return(mean(k_length))
     } else {
          0
     }
}


