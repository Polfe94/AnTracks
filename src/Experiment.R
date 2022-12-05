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
        refcoords = 'data.frame',
        r = 'numeric',
        # edges = 'data.frame', 
        # local_cov = 'numeric',
        date = 'character', 
        food = 'list', 
        N = 'numeric',
        I = 'numeric',
        connectivity = 'numeric',
        trails = 'numeric',
        AiT = 'data.frame',
        IiT = 'data.frame'
), prototype = list(
        r = 51,
        refcoords = hex[hex$y > 1000, ]
))

#### +++ GENERIC FUNCTIONS +++ ####
setGeneric('compute_nodes', function(obj){
        standardGeneric('compute_nodes')
})

setGeneric('activity_matrix', function(obj, binary = TRUE){
        standardGeneric('activity_matrix')
})

setGeneric('interaction_matrix', function(obj){
        standardGeneric('interaction_matrix')
})

setGeneric('get_matrix', function(obj, ...){
        standardGeneric('get_matrix')
})

setGeneric('get_N', function(obj){
        standardGeneric('get_N')
})

setGeneric('get_I', function(obj){
        standardGeneric('get_I')
})

setGeneric('mutual_info', function(obj, t, nodes = NULL){
        standardGeneric('mutual_info')
})

setGeneric('get_fp', function(obj){
        standardGeneric('get_fp')
})

setGeneric('food_detection', function(obj){
        standardGeneric('food_detection')
})

setGeneric('food_trails', function(obj){
        standardGeneric('food_trails')
})

setGeneric('ait', function(obj){
        standardGeneric('ait')
})

setGeneric('iit', function(obj){
        standardGeneric('iit')
})

setGeneric('draw_FoodPatches', function(obj, add = NULL, ...){
        standardGeneric('draw_FoodPatches')
})

setGeneric('plot_trails', function(obj){
        standardGeneric('plot_trails')
})

setGeneric('plot_AiT', function(obj, rectangle_params = list(fill = muted('green'),
                                                             color = 'green', 
                                                             alpha = 0.15, linetype = 2, size = 1)){
        standardGeneric('plot_AiT')
})

setGeneric('plot_IiT', function(obj, rectangle_params = list(fill = muted('green'),
                                                             color = 'green', 
                                                             alpha = 0.15, linetype = 2, size = 1)){
        standardGeneric('plot_IiT')
})

#### +++ CLASS METHODS +++ ####
setMethod('get_N', 'Experiment', function(obj){
        dt <- data.table(obj@data)
        l <- max(dt$Frame)
        t <- seq_len(l)
        N <- numeric(l)
        tmp <- dt[, .N, by = Frame]
        N[tmp$Frame] <- tmp$N
        obj@N <- N
        obj
})

setMethod('get_I', 'Experiment', function(obj){
        dt <- data.table(obj@data)
        l <- max(dt$Frame)
        t <- seq_len(l)
        I <- numeric(l)
        tmp <- dt[Crossings > 0, .N, by = Frame]
        I[tmp$Frame] <- tmp$N
        obj@I <- I
        obj
})

setMethod('compute_nodes', signature = 'Experiment', function(obj){
        if(!'node' %in% colnames(obj@data)){
                
                obj@data$node <- get_node(obj@data[, c('Xmm', 'Ymm')], xy = obj@refcoords)
                
                # clean memmory
                do.call('gc', list(FALSE))
        } 
        obj
})

setMethod('activity_matrix', 'Experiment', function(obj, binary = TRUE){

        obj <- compute_nodes(obj)
        dt <- data.table(obj@data[, c('Frame', 'node')])
        if(binary){
                m <- dt[, as.integer(.N > 0), by = c('node', 'Frame')] 
                colnames(m)[3] <- 'N'
        } else {
                m <- dt[, .N, by = c('node', 'Frame')]
                
        }
        obj@matrix <- m[, c('Frame', 'node', 'N')]
        obj
})

setMethod('interaction_matrix', 'Experiment', function(obj){
        if(nrow(obj@interactions) == 0){
                obj <- compute_nodes(obj)
                dt <- data.table(obj@data[, c('Frame', 'node', 'Crossings')])
                m <- unique(dt[Crossings > 0, ])
                m$N <- 1
                obj@interactions <- m[, c('Frame', 'node', 'N')]
        }
        obj
})

setMethod('get_matrix', 'Experiment', function(obj, ...){
        l <- list(...)
        n <- names(l)
        if(!'data' %in% n){
                data <- -1 
        } else {
                data <- l$data
        }
        if(!'t' %in% n){
                t <- seq_len(max(obj@data$Frame))
        } else {
                t <- l$t
        }
        if(!'maxt' %in% n){
                maxt <- max(obj@data$Frame)
        } else {
                maxt <- l$maxt
        }
        if('type' %in% n && grepl('inter', l$type)){
                
                obj <- interaction_matrix(obj)
                M <- obj@interactions
        } else {
                obj <- activity_matrix(obj)
                M <- obj@matrix
        }
        
        m <- M2M(M, data = data, maxt = maxt)
        M2sbst(m, t)
})

# setMethod('activity_matrix', 'Experiment', function(obj, data = -1){
#      obj <- compute_nodes(obj)
#      n <- obj@data$node
#      t <- obj@data$Frame
#      m <- matrix(data = data, ncol = length(unique(n)), nrow = max(t))
#      dimnames(m) <- list(seq_len(nrow(m)), unique(n))
#      
#      for(i in seq_along(t)){
#           idx <- which(colnames(m) == n[i])
#           m[t[i], idx] <- 1
#      }
#      obj@matrix <- m
#      obj
# })

setMethod('get_fp', 'Experiment', function(obj){
        if(!length(obj@food)){
                d <- obj@date
                if(length(d)){
                        food <- get_foodPatches(d)
                        obj@food <- food
                } else {
                        stop('Invalid experiment date')    
                }
        } 
        obj
})

setMethod('food_detection', 'Experiment', function(obj){
        if(!length(obj@food) | length(obj@food) == 0){
                obj <- get_fp(obj)
        }
        if('t' %in% colnames(obj@food[[1]])){
                return(obj)
        }
        
        food <- data.table(do.call('rbind', obj@food[1:2]))[, node := mapply(get_node, x, y)]
        food$GP <- c(rep(1, 6), rep(2, 6))
        xy <- data.table(obj@data)[node %in% food$node, lapply(.SD, min), .SDcols = c('Frame'), by = node]
        mrgd <- merge(food, xy)
        df <- data.frame(x = mrgd$x, y = mrgd$y, t = mrgd$Frame, GP = mrgd$GP)
        food <- list(GP1 = df[df$GP == 1, c('x', 'y', 't')], GP2 = df[df$GP == 2, c('x', 'y', 't')])
        obj@food <- food
        obj
})

setMethod('mutual_info', 'Experiment', function(obj, t, nodes = NULL){
        obj <- activity_matrix(obj)
        m <- M2M(obj@matrix, data = -1)
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

#' Find food trails (according to: time from first food finding until experiment end + 90 percentile activity)
#' 
#' NOTE!!! There are OTHER criteria that can be apply to calculate food trails;
#' These are not covered in this function for purposes of consistence and simplicity
#' 
#' @param obj an object of class Experiment
#' @return the object with the food trails in the "trails" slot
setMethod('food_trails', 'Experiment', function(obj){
        obj <- activity_matrix(obj, binary = FALSE)
        obj <- food_detection(obj)
        
        m <- obj@matrix
        t <- c(min(do.call('rbind', obj@food)$t), max(obj@data$Frame))
        m <- m[Frame %in% t[1]:t[2], lapply(.SD, sum), .SDcols = 'N', by = node]
        
        obj@trails <- m[N > quantile(N, prob = 0.9)]$node
        obj
})

setMethod('ait', 'Experiment', function(obj){
        if(length(obj@trails) == 0){
                obj@trails <- food_trails(obj)
        }
        if(nrow(obj@AiT) == 0){
                m <- as.data.frame(get_matrix(obj, data = 0))
                df <- data.frame(t = as.integer(rownames(m)),
                                 N = rowSums(m),
                                 inTrail = rowSums(m[, as.character(obj@trails)]),
                                 outTrail = rowSums(m[, !colnames(m) %in% as.character(obj@trails)]))
                mlt_df <- melt(df, id.vars = 't')
                obj@AiT <- mlt_df
        }
        obj
})

setMethod('iit', 'Experiment', function(obj){
        if(length(obj@trails) == 0){
                obj@trails <- food_trails(obj)
        }
        if(nrow(obj@IiT) == 0){
                m <- as.data.frame(get_matrix(obj, data = 0, t = 1:max(obj@data$Frame), type = 'interaction'))
                df <- data.frame(t = as.integer(rownames(m)),
                                 N = rowSums(m),
                                 inTrail = rowSums(m[, colnames(m) %in% as.character(obj@trails)]),
                                 outTrail = rowSums(m[, !colnames(m) %in% as.character(obj@trails)]))
                mlt_df <- melt(df, id.vars = 't')
                obj@IiT <- mlt_df
        }
        obj
})

#### +++ VISUALIZATION METHODS +++ ####
#' Draws the food patches
#' 
#' @param obj An object of "Experiment" class
#' @param add Either NULL or a ggplot object
#' @param ... Additional parameters to be passed to ggplot (color, size, alpha ...)
#' @return A ggplot object drawing the food patches on top of whatever (if any)
#' layer is passed to the argument add.
setMethod('draw_FoodPatches', 'Experiment', function(obj, add = NULL, ...){
        if(is.null(add)){
                add <- ggplot()
        }
        obj <- get_fp(obj)
        
        for(i in seq_along(obj@food[1:2])){
                f <- obj@food[[i]][chull(obj@food[[i]][, c('x', 'y')]), ]
                add <- add + geom_polygon(data = f, aes(x, y), ...) +
                        xlab('') + ylab('')
                
        }
        add
})

setMethod('plot_trails', 'Experiment', function(obj){
        draw_hexagons(add = draw_FoodPatches(obj, fill = 'grey50') +
                              geom_point(data = data.frame(x = hex$x[obj@trails],
                                                           y = hex$y[obj@trails]),
                                         aes(x, y), fill = muted('blue'), size = 4, shape = 21))
})

setMethod('plot_AiT', 'Experiment', function(obj, 
                                             rectangle_params = list(fill = muted('green'),
                                             color = 'green', alpha = 0.15, linetype = 2, size = 1)){
        obj <- ait(obj)
        
        df <- obj@AiT
        food <- obj@food
        
        
        
        ylim <- c(0, 0, rep(max(df$value), 2))
        xlim <- range(do.call('rbind', food)$t)[c(1:2, 2:1)]/120

        pl <- ggplot(data = df, aes(t/120, value, color = variable))+
                geom_line() + scale_y_continuous('', breaks = seq(0, max(ylim), 5)) +
                scale_color_viridis_d('', labels = c('Total activity', 'Activity in trail', 'Activity out trail'))+
                scale_x_continuous('Time (min)', breaks = seq(0, 180, 15))

        rectangle_params$data <- data.frame(x = xlim, y = ylim)
        rectangle_params$mapping <- aes(x, y)
        
        pl <- pl + do.call('geom_polygon', args = rectangle_params) +
                guides(color = guide_legend(override.aes = list(size = 2)))
        
        pl
})

setMethod('plot_IiT', 'Experiment', function(obj, 
                                             rectangle_params = list(fill = muted('green'),
                                                                     color = 'green', alpha = 0.15,
                                                                     linetype = 2, size = 1)){
        obj <- iit(obj)
        
        df <- obj@IiT
        for(i in unique(df$variable)){
                df$value[df$variable == i] <- cumsum(df$value[df$variable == i])
        }
        food <- obj@food
        
        ylim <- c(0, 0, rep(max(df$value), 2))
        xlim <- range(do.call('rbind', food)$t)[c(1:2, 2:1)]/120
        
        pl <- ggplot(data = df, aes(t/120, value, color = variable))+
                geom_line(size = 1.2) + scale_y_continuous('', breaks = seq(0, max(ylim), max(ylim) %/% 8)) +
                scale_color_viridis_d('', labels = c('Total interactions', 'Interactions in trail', 
                                                     'Interactions out trail'),
                                      end = 0.85)+
                scale_x_continuous('Time (min)', breaks = seq(0, 180, 15))
        
        rectangle_params$data <- data.frame(x = xlim, y = ylim)
        rectangle_params$mapping <- aes(x, y)
        
        pl <- pl + do.call('geom_polygon', args = rectangle_params) +
                guides(color = guide_legend(override.aes = list(size = 2)))
        
        pl
})
