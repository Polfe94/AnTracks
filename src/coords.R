node_idx.coords <- function(obj, row){
        m <- cbind(obj$data[row, 'Xmm'],obj$data[row, 'Ymm'])
        h <- cbind(obj$refcoords[, 1], obj$refcoords[, 2])
        d <- pdist(h, m, ret.vec = FALSE)
        idx <- apply(d, 2, which.min)
        as.integer(rownames(obj$refcoords)[idx])
}

get_neighbors.coords <- function(obj, row){
        
        if(!'segments' %in% names(obj)){
                obj$segments <- compute_segments(obj)
        }
        if(sum(grepl('node', colnames(obj$data))) == 0){
                nodes <- node_idx(obj, row)
        } else {
                nodes <- obj$data$node[row]
        }
        
        as.integer(obj$segments$d[obj$segments$o %in% nodes])
}


local_cov.coords <- function(obj, t = c(50, 5)){
        # initialize time vectors
        tvec <- obj$data$Frame
        t0 <- seq(0, max(obj$data$Frame), t[2])
        t1 <- t0 + t[1]
        sq <- lapply(seq_along(t0), function(i){
                which(tvec > t0[i] & tvec <= t1[i])
        })
        
        # initialize a vector of 0s (covariance)
        x <- numeric(length(sq))
        names(sq) <- seq_along(sq)
        
        sq <- sq[lapply(sq, length) > 0]
        idx <- as.integer(names(sq))
        
        # calculate covariance for local neighbors
        y <- vapply(seq_along(sq), function(i){
                data <- obj$data[sq[[i]], c('Frame', 'node')]
                s <- obj$segments[obj$segments$o %in% data$node, c('o', 'd')]
                neighbors <- get_neighbors(obj, sq[[i]])
                data <- rbind(data, cbind(Frame = NA, node = neighbors))
                s <- s[s$d %in% data$node, ]
                t <- table(data)
                t[t == 0] <- -1
                z <- cov(t)
                z <- vapply(seq_len(nrow(s)), function(k){
                        z[rownames(z) == s[k, 1] , colnames(z) == s[k, 2]]
                }, numeric(1)) 
                mean(z, na.rm = T)
        }, numeric(1))
        
        # replace 0s with actual mean covariance
        x[idx] <- y
        x[!is.finite(x)] <- 0
        x
}

pairwise_cov.coords <- function(obj, t){
        m <- coords2matrix(obj)
        m <- m[t, ]
        
        mc <- cov(m)
        z <- numeric(length(det[[1]]$segments$o))
        s <- det[[1]]$segments[, c('o', 'd')]
        n <- colnames(mc)
        for(i in n){
                idx <- which(s$o == i)
                sb <- s[idx, ]
                tmp <- numeric(nrow(sb))
                for(x in seq_len(nrow(sb))){
                        r <- mc[n == sb$o[x], n == sb$d[x]]
                        if(length(r)){
                                tmp[x] <- r
                        }
                }
                if(length(idx) != length(tmp)){
                        print(i)
                }
                z[idx] <- tmp
        }
        z
}

get_N.coords <- function(obj){
        N <- numeric(length(seq(0, max(obj$data$Frame))))
        t <- unique(obj$data$Frame)
        n <- vapply(t, function(k){
                sum(obj$data$Frame == k)
        }, numeric(1))
        
        N[t] <- n
        N
}

get_I.coords <- function(obj){
        I <- numeric(length(seq(0, max(obj$data$Frame))))
        t <- unique(obj$data$Frame[obj$data$Crossings != 0])
        i <- vapply(t, function(k){
                sum(obj$data$Frame == k & obj$data$Crossings > 0)
        }, numeric(1))
        
        I[t] <- i
        I
}

get_foodPatches.coords <- function(obj){
        m <- get_metainfo(obj$date)
        p <- get_food_from_metainfo(m)
        p <- set_foodpatches(p)
        p
}

food_detection.coords <- function(obj, r = 7){
        if(!'food' %in% names(obj)){
                obj$food <- get_foodPatches(obj)
        }
        food <- do.call('rbind', obj$food[1:2]) ## only top
        
        if('t' %in% colnames(food)){
                return(food$t)
        }
        xy <- obj$data[order(obj$data$Frame), c('Xmm', 'Ymm', 'Frame')]
        
        apply(food, 1, function(i){
                idx <- which(
                        (xy$Xmm - i[1]) ^ 2 + (xy$Ymm - i[2]) ^ 2 < (r ^ 2)
                )
                min(xy$Frame[idx])
        })
}

cluster_lengths <- function(x, tau = 10){
        
        idx <- aggregate_time(x$t, tau = tau)
        output <- vapply(idx, function(i){
                result <- connectivity(x[i, c('idx')])
        }, numeric(1))
        
        data.frame(t = c(0, as.numeric(names(idx))),
                   k = c(0, output))
}

connectivity.coords <- function(obj, tau){
        
        data <- obj$data
        
        if(!'node' %in% colnames(data)){
                stop('Prior to connectivity, object must have node computed')
        }
        
        result <- numeric(max(data$Frame))
        
        L <- lapply(unique(data$Frame), function(i){
                idx <- which(data$Frame %in% i:(i+tau))
                list(Frame = i, Pos = data$node[idx])
        })
        
        for(i in seq_along(L)){
                
                pos <- unique(L[[i]]$Pos)
                k_length <- c()
                branch <- c()
                
                while(length(pos)){
                        
                        current_path <- pos[length(pos)]
                        pos <- pos[-length(pos)]
                        
                        while(length(current_path)){
                                
                                target <- current_path[length(current_path)]
                                neighbors <- as.integer(obj$segments$d[obj$segments$o %in% target])
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
                                        neighbors <- as.integer(obj$segments$d[obj$segments$o %in% target])
                                        idx <- which(neighbors %in% pos)
                                        branch <- c(branch, neighbors[idx])
                                }
                                k_length <- c(k_length, length(current_path))
                                current_path <- c()
                                branch <- c()
                        }
                }
                result[L[[i]]$Frame] <- mean(k_length)
        }
        
        return(result)
}

# ETAfood <- function(obj, r = 3.5){
#         food <- do.call('rbind', obj$food)
#         xy <- obj$data[order(obj$data$Frame), c('Xmm', 'Ymm', 'Frame')]
#         
#         f <- apply(food, 1, function(i){
#                 idx <- which(
#                         (xy$Xmm - i[1]) ^ 2 + (xy$Ymm - i[2]) ^ 2 < (r ^ 2)
#                 )
#                 min(xy$Frame[idx])
#         })
# }
# 
# 
# ## gives an estimated time of arrival to the food vertices (returns a dataframe)
# ## note that the function assumes the time to be ordered
# ## need to provide the patch configuration (foodpatches = list of points (vertices of the foodpatches))
# ETAfood <- function(foodpatches, exp, radius = 3.5){
#         exp <- exp[order(exp$Time_sec),]
#         coordinates <- do.call('rbind', foodpatches)
#         positions <- vector(mode = 'list', length = nrow(coordinates))
#         ## take all left positions in the specified radius near the food vertices
#         positions <- apply(
#                 coordinates, 1, function(i) exp[which(exp$Xmm >= i[1] - radius &
#                                                               exp$Xmm < i[1] + radius &
#                                                               exp$Ymm >= i[2] - radius &
#                                                               exp$Ymm < i[2] + radius),c('Time_sec','Xmm', 'Ymm')]
#         )
#         names(positions)
#         positions <- positions[lapply(positions, nrow)>0]
#         times <- do.call('rbind', lapply(
#                 positions, function(i) i[1,]
#         ))
#         return(times)
# }

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
