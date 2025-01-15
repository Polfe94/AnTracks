## --- LOAD GENERIC FUNCTIONS --- ##
Path2File <- normalizePath(paste0(gsub('Experiment.R', '', paste0(sys.frames()[[1]]$ofile)), '/functions.R'))
source(Path2File)

Path2File <- normalizePath(paste0(gsub('Experiment.R', '', paste0(sys.frames()[[1]]$ofile)), '/visualization.R'))
source(Path2File)

## --- LOAD HEXAGONAL REFERENCE COORDINATES & SPREADSHEET REFERENCE --- ##
Path2File <- normalizePath(paste0(gsub('Experiment.R', '', paste0(sys.frames()[[1]]$ofile)), '/../data/hex.csv'))
if(!file.exists(Path2File)){
        warning("File hex.csv does not exist. Trying to download from repository.")
        download.file("https://github.com/Polfe94/AnTracks/tree/main/data/hex.csv",
                      Path2File)
} 
hex <- read.csv(Path2File)
edges <- compute_edges()


Path2File <- normalizePath(paste0(gsub('Experiment.R', '', paste0(sys.frames()[[1]]$ofile)), '/../data/Experiments_Spreadsheet_ANTS2018.csv'))
if(!file.exists(Path2File)){
        warning("Could not load Experiments_Spreadsheet_ANTS2018.csv: file does not exist.")
} else {
        exp_spreadsheet <- read.csv2(Path2File)
} 
rm(Path2File)

#### +++ CLASS DEFINITION +++ ####
setClass('Experiment', representation(
        data = 'data.frame',
        type = 'character',
        matrix = 'data.frame',
        interactions = 'data.frame',
        date = 'character', 
        food = 'list', 
        N = 'numeric',
        I = 'numeric',
        connectivity = 'list',
        trails = 'integer',
        nest = 'integer', 
        AiT = 'data.frame',
        IiT = 'data.frame',
        .__classVersion__ = 'character'
), prototype = list(
        nest = nest_influence(r = 101), 
        .__classVersion__ = '1.0'
))

setGeneric('updateObject', function(.Object){
        standardGeneric('updateObject')
})

setMethod('updateObject', 'Experiment', function(.Object){
        if(.hasSlot(.Object, '.__classVersion__') && .Object@.__classVersion__ == '1.0'){
                return(.Object)
        } else {
                attrs <- getSlots('Experiment')
                slots <- list(Class = 'Experiment')
                for(i in names(attrs)){
                        if(.hasSlot(.Object, i)){
                                x <- slot(.Object, i)
                                if(attrs[i] %in% class(x)){
                                        slots[[i]] <- x
                                }
                        }
                }
                .Object <- do.call('new', args = slots)
                .Object
        }
})

#### +++ GENERIC FUNCTIONS +++ ####
setGeneric('compute_nodes', function(.Object){
        standardGeneric('compute_nodes')
})

setGeneric('activity_matrix', function(.Object, binary = TRUE){
        standardGeneric('activity_matrix')
})

setGeneric('interaction_matrix', function(.Object){
        standardGeneric('interaction_matrix')
})

setGeneric('get_matrix', function(.Object, ...){
        standardGeneric('get_matrix')
})

setGeneric('get_N', function(.Object){
        standardGeneric('get_N')
})

setGeneric('get_I', function(.Object){
        standardGeneric('get_I')
})

# setGeneric('compute_connectivity', signature = c('Experiment'),
#            def = function(.Object, t, edges){
#         standardGeneric('compute_connectivity')
# })

setGeneric('mutual_info', function(.Object, t, edges, nodes = NULL){
        standardGeneric('mutual_info')
})

setGeneric('get_fp', function(.Object){
        standardGeneric('get_fp')
})

setGeneric('food_detection', function(.Object){
        standardGeneric('food_detection')
})

setGeneric('food_trails', function(.Object){
        standardGeneric('food_trails')
})

setGeneric('optimal_food_trails', function(.Object){
        standardGeneric('optimal_food_trails')
})

setGeneric('ait', function(.Object, include_nest = FALSE){
        standardGeneric('ait')
})

setGeneric('iit', function(.Object, include_nest = FALSE){
        standardGeneric('iit')
})

setGeneric('draw_FoodPatches', function(.Object, add = NULL, ...){
        standardGeneric('draw_FoodPatches')
})

setGeneric('plot_trails', function(.Object){
        standardGeneric('plot_trails')
})

setGeneric('plot_AiT', function(.Object, ...){
        standardGeneric('plot_AiT')
})

setGeneric('plot_IiT', function(.Object, norm = FALSE, ...){
        standardGeneric('plot_IiT')
})

setGeneric('get_eff', function(.Object){
        standardGeneric('get_eff')
})

# setGeneric('move_stats', function(.Object, scouts = 'LR', maxt = 2400, maxd = 500, min_pos = 10){
#         standardGeneric('move_stats')
# })

setGeneric('move_stats', function(.Object, ...){
        standardGeneric('move_stats')
})

setGeneric('get_straightness', function(.Object, ...){
        standardGeneric('get_straightness')
})

setGeneric('class_ids', function(.Object, dists = list(det = 200, sto = NA, nf = 350)){
        standardGeneric('class_ids')
})

#### +++ CLASS METHODS +++ ####
setMethod('get_N', 'Experiment', function(.Object){
        dt <- data.table(.Object@data)
        l <- max(dt$Frame)
        t <- seq_len(l)
        N <- numeric(l)
        tmp <- dt[, .N, by = Frame]
        N[tmp$Frame] <- tmp$N
        .Object@N <- N
        .Object
})

setMethod('get_I', 'Experiment', function(.Object){
        dt <- data.table(.Object@data)
        l <- max(dt$Frame)
        t <- seq_len(l)
        I <- numeric(l)
        tmp <- dt[Crossings > 0, .N, by = Frame]
        I[tmp$Frame] <- tmp$N
        .Object@I <- I
        .Object
})


setMethod('compute_nodes', signature = 'Experiment', function(.Object){
        if(!'node' %in% colnames(.Object@data)){
                
                if(.hasSlot(.Object, 'refcoords')){
                        .Object@data$node <- get_node(.Object@data[, c('Xmm', 'Ymm')], xy = .Object@refcoords) 
                } else {
                        .Object@data$node <- get_node(.Object@data[, c('Xmm', 'Ymm')]) 
                }

                
                # clean memmory
                do.call('gc', list(FALSE))
        } 
        .Object
})

setMethod('activity_matrix', 'Experiment', function(.Object, binary = TRUE){

        .Object <- compute_nodes(.Object)
        dt <- data.table(.Object@data[, c('Frame', 'node')])
        if(binary){
                m <- dt[, as.integer(.N > 0), by = c('node', 'Frame')] 
                colnames(m)[3] <- 'N'
        } else {
                m <- dt[, .N, by = c('node', 'Frame')]
                
        }
        .Object@matrix <- m[, c('Frame', 'node', 'N')]
        .Object
})

setMethod('interaction_matrix', 'Experiment', function(.Object){
        if(nrow(.Object@interactions) == 0){
                .Object <- compute_nodes(.Object)
                dt <- data.table(.Object@data[, c('Frame', 'node', 'Crossings')])
                m <- unique(dt[Crossings > 0, ])
                m$N <- 1
                .Object@interactions <- m[, c('Frame', 'node', 'N')]
        }
        .Object
})

setMethod('get_matrix', 'Experiment', function(.Object, ...){
        l <- list(...)
        n <- names(l)
        if(!'data' %in% n){
                data <- -1 
        } else {
                data <- l$data
        }
        if(!'t' %in% n){
                t <- seq_len(max(.Object@data$Frame))
        } else {
                t <- l$t
        }
        if(!'maxt' %in% n){
                maxt <- max(.Object@data$Frame)
        } else {
                maxt <- l$maxt
        }
        if('type' %in% n && grepl('inter', l$type)){
                
                .Object <- interaction_matrix(.Object)
                M <- .Object@interactions
        } else {
                .Object <- activity_matrix(.Object)
                M <- .Object@matrix
        }
        
        m <- M2M(M, data = data, maxt = maxt)
        M2sbst(m, t)
})

# setMethod('activity_matrix', 'Experiment', function(.Object, data = -1){
#      .Object <- compute_nodes(.Object)
#      n <- .Object@data$node
#      t <- .Object@data$Frame
#      m <- matrix(data = data, ncol = length(unique(n)), nrow = max(t))
#      dimnames(m) <- list(seq_len(nrow(m)), unique(n))
#      
#      for(i in seq_along(t)){
#           idx <- which(colnames(m) == n[i])
#           m[t[i], idx] <- 1
#      }
#      .Object@matrix <- m
#      .Object
# })

setMethod('get_fp', 'Experiment', function(.Object){
        if(!length(.Object@food)){
                d <- .Object@date
                if(length(d)){
                        food <- get_foodPatches(d)
                        .Object@food <- food
                } else {
                        stop('Invalid experiment date')    
                }
        } 
        .Object
})

setMethod('food_detection', 'Experiment', function(.Object){
        if(!length(.Object@food) | length(.Object@food) == 0){
                .Object <- get_fp(.Object)
        }
        if('t' %in% colnames(.Object@food[[1]])){
                return(.Object)
        }
        
        food <- data.table(do.call('rbind', .Object@food[1:2]))[, node := mapply(get_node, x, y)]
        food$GP <- c(rep(1, 6), rep(2, 6))
        xy <- data.table(.Object@data)[node %in% food$node, lapply(.SD, min), .SDcols = c('Frame'), by = node]
        mrgd <- merge(food, xy)
        df <- data.frame(x = mrgd$x, y = mrgd$y, t = mrgd$Frame, GP = mrgd$GP)
        food <- list(GP1 = df[df$GP == 1, c('x', 'y', 't')], GP2 = df[df$GP == 2, c('x', 'y', 't')])
        .Object@food <- food
        .Object
})

setMethod('mutual_info', 'Experiment', function(.Object, t, edges = edges, nodes = NULL){
        .Object <- activity_matrix(.Object)
        m <- M2M(.Object@matrix, data = -1)
        m <- M2sbst(m, t)
        if(!exists('edges', envir = .GlobalEnv)){
                edges <<- compute_edges()
        }
        z <- numeric(nrow(edges))
        if(!is.null(nodes)){
                existing_edges <- edges[o %in% nodes & d %in% nodes]
        } else {
                existing_edges <- edges[o %in% as.integer(colnames(m)) & d %in% as.integer(colnames(m))]
        }
        
        for(i in seq_len(nrow(existing_edges))){
                idx <- edges[o == existing_edges[i, o] & d == existing_edges[i, d], which = TRUE]
                n <- as.character(existing_edges[i, c('o', 'd')])
                z[idx] <- infotheo::mutinformation(m[, ..n])[1, 2]
        }
        z
})



setMethod('compute_connectivity', signature = c('Experiment', 'Simulation'),
          function(.Object, t, edges){
                  
                  .Object@results[['k']] <- .Object@data[, .(k = lapply(list(node), function(i){
                          mean(connectivity(unlist(i), edges))
                  })), by = .(Frame = ((Frame -1) %/% t) *t)]
                  
                  .Object
                  
          })

#' Find food trails (according to: time from first food finding until experiment end + 90 percentile activity)
#' 
#' NOTE!!! There are OTHER criteria that can be apply to calculate food trails;
#' These are not covered in this function for purposes of consistence and simplicity
#' 
#' @param .Object an object of class Experiment
#' @return the object with the food trails in the "trails" slot
setMethod('food_trails', 'Experiment', function(.Object){
        .Object <- activity_matrix(.Object, binary = FALSE)
        .Object <- food_detection(.Object)
        
        m <- .Object@matrix
        t <- c(min(do.call('rbind', .Object@food)$t), max(.Object@data$Frame))
        m <- m[Frame %in% t[1]:t[2], lapply(.SD, sum), .SDcols = 'N', by = node]
        
        .Object@trails <- m[N > quantile(N, prob = 0.9)]$node
        .Object
})

setMethod('optimal_food_trails', 'Experiment', function(.Object){
        targets <- c(get_node(.Object@food[[1]][, 1:2]), get_node(.Object@food[[2]][, 1:2]))
        origins <- nest_influence(r = 101)
        # origins <- origins[origins %in% hex$node[hex$y > 1030]]
        
        nodes <- c()
        for(i in seq_along(origins)){
                for(j in seq_along(targets)){
                        nodes <- c(nodes, optimal_path(origins[i], targets[j]))
                }
        }
        nest <- nest_influence(r = 199)
        .Object@trails <- unique(c(nodes, nest[nest %in% hex$node[hex$y > 1030]]))
        .Object
})

setMethod('ait', 'Experiment', function(.Object, include_nest = FALSE){
        if(length(.Object@trails) == 0){
                .Object@trails <- food_trails(.Object)
        }

        m <- as.data.frame(get_matrix(.Object, data = 0))
        if(include_nest){
                
                trails <- .Object@trails[!.Object@trails %in% .Object@nest]
                df <- data.frame(t = as.integer(rownames(m)),
                                 N = rowSums(m),
                                 inTrail = rowSums(m[, colnames(m) %in% as.character(trails)]),
                                 nest = rowSums(m[, as.character(.Object@nest)]),
                                 outTrail = rowSums(m[, !colnames(m) %in% as.character(.Object@trails)]))
        } else {
                df <- data.frame(t = as.integer(rownames(m)),
                                 N = rowSums(m),
                                 inTrail = rowSums(m[, as.character(.Object@trails)]),
                                 outTrail = rowSums(m[, !colnames(m) %in% as.character(.Object@trails)]))
        }
        
        mlt_df <- melt(df, id.vars = 't')
        .Object@AiT <- mlt_df

        .Object
})

setMethod('iit', 'Experiment', function(.Object, include_nest = FALSE){
        if(length(.Object@trails) == 0){
                .Object@trails <- food_trails(.Object)
        }

        m <- as.data.frame(get_matrix(.Object, data = 0, t = 1:max(.Object@data$Frame), type = 'interaction'))
        if(include_nest){
                trails <- as.character(.Object@trails[!.Object@trails %in% .Object@nest])
                df <- data.frame(t = as.integer(rownames(m)),
                                 N = rowSums(m),
                                 inTrail = rowSums(m[, colnames(m) %in% as.character(trails)]),
                                 nest = rowSums(m[, colnames(m) %in% as.character(.Object@nest)]),
                                 outTrail = rowSums(m[, !colnames(m) %in% as.character(.Object@trails)]))
        } else {
                df <- data.frame(t = as.integer(rownames(m)),
                                 N = rowSums(m),
                                 inTrail = rowSums(m[, colnames(m) %in% as.character(.Object@trails)]),
                                 outTrail = rowSums(m[, !colnames(m) %in% as.character(.Object@trails)]))
        }

        mlt_df <- melt(df, id.vars = 't')
        .Object@IiT <- mlt_df

        .Object
})


setMethod('get_eff', 'Experiment', function(.Object){
        food <- range(rbindlist(.Object@food)[['t']])
        dt <- setDT(.Object@data)
        filterdt <- dt[Frame <= food[2], .(id = N_ind), by = 'Frame']
        M <- dcast(data = filterdt, Frame ~ id, value.var = 'id')
        .tmp <- rbindlist(lapply(2:ncol(M), function(i){
                idx <- which(!is.na(M[[i]]))
                data.frame(id = unique(M[[i]][idx]),Frame = M[[1]][idx], d=c(1, diff(idx)))
        }))
        .tmp[d > 1, 'd'] <- 0
        .tmp[, interval := cumsum(c(TRUE, diff(d) != 0)), by = id]
        .tmp[, interval := ifelse(interval %% 2 == 0, interval + 1, interval), by = id]
        .tmp_result <- .tmp[, .(start_frame = min(Frame), end_frame = max(Frame)), by = .(id, interval)]
        
        .tp1 <- .tmp_result[start_frame <= food[1]][, end_frame := ifelse(end_frame < food[1], end_frame, food[1])]
        tp1_value <- sum(.tp1[, end_frame] - .tp1[, start_frame])/2
        .tp2 <- .tmp_result[start_frame <= food[2] & end_frame >= food[1]][, start_frame := ifelse(start_frame < food[1], food[1], start_frame)]
        .tp2 <- .tp2[, end_frame := ifelse(end_frame < food[2], end_frame, food[2])]
        tp2_value <- sum(.tp2[, end_frame] - .tp2[, start_frame])/2
        data.table(tp1 = tp1_value, tp2 = tp2_value)
})

setMethod('move_stats', 'Experiment', function(.Object, ...){
        require(trajr)
        .Object <- class_ids(.Object, ...)
        data <- .Object@data[type != 'UNKNOWN']
        ids <- data[, unique(N_ind)]
        
        rbindlist(lapply(seq_along(ids), function(i){
                sbst <- data[N_ind == ids[i]]
                xy <- TrajFromCoords(sbst[, c('Xmm', 'Ymm', 'Frame')], fps = 2)
                
                # slight filter
                if(nrow(xy) > 10){
                        d <- which.max(as.numeric(pdist(as.matrix(sbst[, c('Xmm', 'Ymm')]), 
                                                        as.matrix(hex[hex$node == 634, c('x', 'y')]))))
                        df <- data.frame(v = mean(Mod(TrajVelocity(xy)), na.rm = TRUE),
                                         acc = mean(Mod(TrajAcceleration(xy)), na.rm = TRUE),
                                         # l = TrajLength(xy),
                                         d = TrajDistance(xy, endIndex = d),
                                         t = TrajDuration(xy))
                } else {
                        # df <- data.frame(v = NA, acc = NA, l = NA, d = NA, t = NA)
                        df <- data.frame(v = NA, acc = NA, d = NA, t = NA)
                }
                df[['type']] <- sbst[, unique(type)]
                df[d > 0]
        }), idcol = 'ind')
})

setMethod('get_straightness', 'Experiment', function(.Object, ...){
        require(trajr)
        .Object <- class_ids(.Object, ...)
        data <- .Object@data[type != 'UNKNOWN']
        ids <- data[, unique(N_ind)]
        
        rbindlist(lapply(seq_along(ids), function(i){
                sbst <- data[N_ind == ids[i]]
                xy <- TrajFromCoords(sbst[, c('Xmm', 'Ymm', 'Frame')], fps = 2)
                
                # slight filter
                if(nrow(xy) > 10){
                        Ls <- vapply(2:nrow(xy), function(x) TrajLength(xy, endIndex = x), numeric(1))
                        maxL <- max(Ls)
                        if(maxL >= 100){
                                idxs <- vapply(seq(100, max(Ls), 100), function(x) which.min(abs(Ls - x)), numeric(1))
                        } else {
                                idxs <- nrow(xy)
                        }
                        straightness <- vapply(idxs, function(l) TrajStraightness(xy[seq_len(l)]), numeric(1))
                        data.frame(l = Ls[idxs], s = straightness)
                } else {
                        df <- data.frame(l = NA, s = NA)
                }
                df[['type']] <- sbst[, unique(type)]
        }), idcol = 'ind')
})


# setMethod('move_stats', 'Experiment', function(.Object, scouts = 'LR', maxt = 2400, maxd = 500, min_pos = 10){
#         data <- setDT(.Object@data)[Frame <= maxt]
#         dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
#         set(data, j = 'd', value = as.numeric(dmatrix))
#         inds <- unique(data[d > maxd, N_ind])
#         if(scouts == 'SR'){
#                 inds <- unique(data[!N_ind %in% inds, N_ind])
#                 # idx <- data[, .(idx = which.min(abs((d - maxd)))), by = 'N_ind'][['idx']]
#         } else {
#                 # idx <- data[, .(idx = which.max((d - maxd)> 0)), by = 'N_ind'][['idx']]
#         }
#         
#         data <- data[N_ind %in% inds]
#         
#         result <- as.data.frame(t(vapply(seq_along(inds), function(ii){
#                 sbst <- data[N_ind == inds[ii]]
#                 xy <- TrajFromCoords(sbst[, c('Xmm', 'Ymm', 'Frame')], fps = 2)
#                 if(nrow(xy) > min_pos){
#                         c(TrajStraightness(xy), TrajDistance(xy), TrajLength(xy), 
#                           mean(as.numeric(Mod(TrajVelocity(xy))), na.rm = TRUE), 
#                           mean(as.numeric(Mod(TrajAcceleration(xy))), na.rm = TRUE),
#                           nrow(xy) / 120, nrow(sbst)/120)
#                 } else {rep(0, 7)}
#                 
#         }, numeric(7))))
#         colnames(result) <- c('Straightness', 'Diffusion', 'Distance', 'Mean_v', 'Mean_acc', 'Time', 'Total_time')
#         result[result[['Time']] > 0, ]
# })

# setMethod('class_ids',
#           'Experiment', # dists are the breakpoint between the two <distance> scales
#           function(.Object, dists = list(det = 200, sto = NA, nf = 350)){
#         
#         if('type' %in% colnames(.Object@data)){
#                 return(.Object)
#         }
#         
#         setDT(.Object@data)
#         .Object@data[['type']] <- 'UNKNOWN'
#         
#         # init max time and distance to classify ants
#         maxd <- dists[[.Object@type]]
#         if(length(.Object@food)){
#                 maxt <- min(rbindlist(.Object@food)[['t']])
#         } else {
#                 maxt <- .Object@data[, max(Frame)]		
#         }
#         
#         # subset data
#         data <- setDT(.Object@data)[Frame <= maxt]
#         data[['d']] <- as.numeric(pdist(as.matrix(data[, c('Xmm', 'Ymm')]), 
#                                         as.matrix(hex[hex$node == 634, c('x', 'y')])))
#         
#         ids <- data[, unique(N_ind)]
#         scouts <- unique(data[d > maxd, N_ind])
#         recruits <- ids[!ids %in% scouts]
#         
#         .Object@data[N_ind %in% scouts, 'type'] <- 'Scout'
#         .Object@data[N_ind %in% recruits, 'type'] <- 'Recruit'
#         
#         .Object
# })

setMethod('class_ids',
          'Experiment', # dists are the breakpoint between the two <distance> scales
          function(.Object, dists = list(det = 200, sto = NA, nf = 350)){
                  
                  if('type' %in% colnames(.Object@data)){
                          return(.Object)
                  }
                  
                  setDT(.Object@data)
                  .Object@data[['type']] <- 'UNKNOWN'
                  
                  # init max time and distance to classify ants
                  maxd <- dists[[.Object@type]]
                  if(length(.Object@food)){
                          maxt <- min(rbindlist(.Object@food)[['t']])
                  } else {
                          maxt <- .Object@data[, max(Frame)]		
                  }
                  
                  # subset data
                  ids_filtered <- setDT(.Object@data)[Frame <= maxt, .N, by='N_ind'][N > 10, N_ind]
                  
                  data <- .Object@data[N_ind %in% ids_filtered]
                  
                  data[['d']] <- as.numeric(pdist(as.matrix(data[, c('Xmm', 'Ymm')]), 
                                                  as.matrix(hex[hex$node == 634, c('x', 'y')])))
                  
                  ids <- data[, .(d = mean(d)), by = 'N_ind']
                  scouts <- unique(ids[d > maxd, N_ind])
                  recruits <- unique(ids[d <= maxd, N_ind])
                  
                  .Object@data[N_ind %in% scouts, 'type'] <- 'Scout'
                  .Object@data[N_ind %in% recruits, 'type'] <- 'Recruit'
                  
                  .Object
          })



#### +++ VISUALIZATION METHODS +++ ####
#' Draws the food patches
#' 
#' @param .Object An object of "Experiment" class
#' @param add Either NULL or a ggplot object
#' @param ... Additional parameters to be passed to ggplot (color, size, alpha ...)
#' @return A ggplot object drawing the food patches on top of whatever (if any)
#' layer is passed to the argument add.
setMethod('draw_FoodPatches', 'Experiment', function(.Object, add = NULL, ...){
        if(is.null(add)){
                add <- ggplot()
        }
        .Object <- get_fp(.Object)
        
        for(i in seq_along(.Object@food[1:2])){
                f <- .Object@food[[i]][chull(.Object@food[[i]][, c('x', 'y')]), ]
                add <- add + geom_polygon(data = f, aes(x, y), ...) +
                        xlab('') + ylab('')
                
        }
        add
})

setMethod('plot_trails', 'Experiment', function(.Object){
        draw_hexagons(add = draw_FoodPatches(.Object, fill = 'grey50') +
                              geom_point(data = data.frame(x = hex$x[.Object@trails],
                                                           y = hex$y[.Object@trails]),
                                         aes(x, y), fill = muted('blue'), size = 4, shape = 21))
})

setMethod('plot_AiT', 'Experiment', function(.Object, ...){
        
        df <- .Object@AiT
        food <- .Object@food
        
        ylim <- c(0, 0, rep(max(df$value), 2))
        t <- do.call('rbind', food)$t / 120
        xlim <- range(t[c(1:2, 2:1)])

        if(length(levels(df$variable)) == 4){
                pl <- ggplot(data = df, aes(t/120, value, color = variable))+
                        geom_line() + scale_y_continuous('', breaks = seq(0, max(ylim), 5)) +
                        scale_color_viridis_d('', labels = c('Total activity', 'Activity in trail',
                                                             'Activity in nest',
                                                             'Activity out trail'))+
                        scale_x_continuous('Time (min)', breaks = seq(0, 180, 15))
        } else {
                pl <- ggplot(data = df, aes(t/120, value, color = variable))+
                        geom_line() + scale_y_continuous('', breaks = seq(0, max(ylim), 5)) +
                        scale_color_viridis_d('', labels = c('Total activity', 'Activity in trail',
                                                             'Activity out trail'))+
                        scale_x_continuous('Time (min)', breaks = seq(0, 180, 15))
        }

        pl + geom_food(t = t, ylim = range(ylim), ...) +
                guides(color = guide_legend(override.aes = list(size = 2)))
})

setMethod('plot_IiT', 'Experiment', function(.Object, norm = FALSE, ...){
        
        df <- .Object@IiT
        if(norm){
                for(i in unique(df$variable)){
                        y <- cumsum(df$value[df$variable == i])
                        df$value[df$variable == i] <- norm_range(y, 0, 1)
                }
                # ybrks <- 1/5
        } else {
                for(i in unique(df$variable)){
                        df$value[df$variable == i] <- cumsum(df$value[df$variable == i])
                }
                # ybrks <- max(df$value) %/% 8
        }
        
        food <- .Object@food

        ylim <- c(0, 0, rep(max(df$value), 2))
        t <- do.call('rbind', food)$t / 120
        xlim <- range(t[c(1:2, 2:1)])
        
        
        if(length(levels(df$variable)) == 4){
                pl <- ggplot(data = df, aes(t/120, value, color = variable))+
                        geom_line(size = 1.2) + scale_y_continuous('', breaks = seq(0, max(ylim), length.out = 6)) +
                        scale_color_viridis_d('', labels = c('Total interactions', 'Interactions in trail',
                                                             'Interactions in nest',
                                                             'Interactions out trail'),
                                              end = 0.85)+
                        scale_x_continuous('Time (min)', breaks = seq(0, 180, 15))
        } else {
                pl <- ggplot(data = df, aes(t/120, value, color = variable))+
                        geom_line(size = 1.2) + scale_y_continuous('', breaks = seq(0, max(ylim), length.out = 9)) +
                        scale_color_viridis_d('', labels = c('Total interactions', 'Interactions in trail', 
                                                             'Interactions out trail'),
                                              end = 0.85)+
                        scale_x_continuous('Time (min)', breaks = seq(0, 180, 15))
        }
        
        pl + geom_food(t = t, ylim = range(ylim), ...) +
                guides(color = guide_legend(override.aes = list(size = 2)))
})


