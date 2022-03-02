## Class setter
setClass('coords', slots = list(data = 'data.frame',
                                exp.condition = 'character',
                                is.average = 'logical',
                                refcoords = 'data.frame'))

## Class initializer


class(obj) <- 'lattice'
obj <- new('coords', data = obj, ...)
obj
setMethod('nest_boundaries', 'coords', function(obj){
     obj@nest_boundaries <- ...
     obj
})
nest_boundaries <- function(){
     tnest <- closest.node(1000,1000)
     bnest <- closest.node(980, 980)
     tnest <- calculate_tiers(tnest, table = 'htop', lim = 6, subset = 250)
     bnest <- calculate_tiers(bnest, table = 'hbot', lim = 6, subset = 250)
     return(list(green = tnest[,1:2], red = bnest[,1:2]))
}


coords2idx <- function(df, par = F, threads = 4L){
     if(par){
          cl <- makeCluster(threads)
          clusterExport(cl, varlist = c('closest.node', 'command', 'min_move', 'pdist', 'hex90'))
          x <- parLapply(cl, 1:nrow(df), function(i) closest.node(df[i,], get.idx = T))
          stopCluster(cl)
          as.integer(x)
     } else {
          sapply(1:nrow(df), function(i) closest.node(df[i,], get.idx = T))
     }
}


# Returns a list of the possible movements from the perspective of the input node.
#' @param currentNode vector of length 2 (a coordinate X, Y) as returned by \code{closest.node}
#' @return named list of length 1, 2 or 3, depending on the amount of neighbors
#' @example command(closest.node(40, 40))
command <- function(currentNode){
     r <- 55 # radius to consider available neighbours
     
     commands <- vector('list',4)
     
     names(commands) <- c('north','east','south','west')
     commands[[1]] <- which(hex90[,1] == currentNode[,1] &
                                 hex90[,2] > currentNode[,2] &
                                 hex90[,2] <= currentNode[,2]+r)
     commands[[2]] <- which(hex90[,1] > currentNode[,1] &
                                 hex90[,1] <= currentNode[,1]+r &
                                 hex90[,2] > currentNode[,2] -r &
                                 hex90[,2] < currentNode[,2] +r)
     commands[[3]] <- which(hex90[,1] == currentNode[,1] &
                                 hex90[,2] < currentNode[,2] &
                                 hex90[,2] >= currentNode[,2]-r)
     commands[[4]] <- which(hex90[,1] < currentNode[,1] &
                                 hex90[,1] >= currentNode[,1]-r &
                                 hex90[,2] > currentNode[,2] - r&
                                 hex90[,2] <= currentNode[,2] + r)
     
     
     # take only available positions
     commands <- commands[lapply(commands, length)>0]
     return(commands)
}
