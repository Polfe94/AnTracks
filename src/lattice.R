## Class setter
setClass('lattice', slots = list(data = 'matrix', 
                                 exp.condition = 'character', 
                                 is.average = 'logical',
                                 refcoords = 'data.frame'))

## Class initializer
init.lattice <- function(obj, ...){
     tmp <- list(data = obj)
     args <- list(...)
     if(length(args) && length(names(args)) && length(args) == length(names(args))){
          n <- names(args)
          for(i in seq_along(args)){
               tmp[[n[i]]] <- args[[i]]
          }
     }
     tmp
}