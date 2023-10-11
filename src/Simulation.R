## --- LOAD GENERIC FUNCTIONS --- ##
Path2File <- normalizePath(paste0(gsub('Simulation.R', '', paste0(sys.frames()[[1]]$ofile)), '/functions.R'))
source(Path2File)

Path2File <- normalizePath(paste0(gsub('Simulation.R', '', paste0(sys.frames()[[1]]$ofile)), '/visualization.R'))
source(Path2File)

## --- LOAD HEXAGONAL REFERENCE COORDINATES & SPREADSHEET REFERENCE --- ##
Path2File <- normalizePath(paste0(gsub('Simulation.R', '', paste0(sys.frames()[[1]]$ofile)), '../data/hex_sims.csv'))
if(!file.exists(Path2File)){
	warning("File hex_sims.csv does not exist. Trying to download from repository.")
	download.file("https://github.com/Polfe94/AnTracks/tree/main/data/hex_sims.csv",
		      Path2File)
} 
hex_sim <- setDT(read.csv(Path2File, sep = ',', dec = '.', stringsAsFactors = F))
sim_edges <- compute_edges(refcoords = hex_sim, r = 1.01)
rm(Path2File)

#### +++ CLASS DEFINITION +++ ####
setClass('Simulation', representation(
	data = 'data.frame',
	results = 'list',
	type = 'character',
	# food = 'data.frame',
	.__classVersion__ = 'character'
), prototype = list(
	# food = data.frame(node = c("(6, 33)", "(6, 34)", "(7, 34)", # patch 1
	#                            "(7, 33)", "(7, 32)", "(6, 32)",
	#                            "(6, 11)", "(6, 12)", "(7, 12)", # patch 2
	#                            "(7, 11)", "(7, 10)", "(6, 10)")),
	.__classVersion__ = '1.0'
))

setGeneric('updateObject', function(.Object){
	standardGeneric('updateObject')
})

setMethod('updateObject', 'Simulation', function(.Object){
	if(.hasSlot(.Object, '.__classVersion__') && .Object@.__classVersion__ == '1.0'){
		return(.Object)
	} else {
		attrs <- getSlots('Simulation')
		slots <- list(Class = 'Simulation')
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

# is called when the new object is created
setMethod('initialize', 'Simulation', function(.Object, ...){
	.Object <- callNextMethod()
	setDT(.Object@data)
	.Object@data[['node']] <- .Object@data[, .(node = lapply(pos, parse_nodes)), by = 'Frame'][['node']]
	# .Object@data[['Frame']] <- NULL # don't know why the column "Frame" gets duplicated...
	.Object
})

setGeneric('compute_connectivity', function(.Object, t, edges){
	standardGeneric('compute_connectivity')
})

# setGeneric('process_data', function(.Object, frames = 21600,
# 				    cols = c('N', 'I', 'SiOut', 'pos', 'k')){
# 	standardGeneric('process_data')
# })

setGeneric('activity_matrix', function(.Object, binary = TRUE){
	standardGeneric('activity_matrix')
})

setGeneric('get_matrix', function(.Object, ...){
	standardGeneric('get_matrix')
})

setGeneric('mutual_info', function(.Object, t, edges, nodes = NULL){
	standardGeneric('mutual_info')
})

#### +++ CLASS METHODS +++ ####
setMethod('compute_connectivity', signature = c('Experiment', 'Simulation'),
	  function(.Object, t, edges){
	
	.Object@results[['k']] <- .Object@data[, .(k = lapply(list(node), function(i){
		mean(connectivity(unlist(i), edges))
	})), by = .(Frame = ((Frame -1) %/% t) *t)]

	.Object

})

# setMethod('process_data', signature = 'Simulation',
# 	  function(.Object, frames = 21600, cols = c('N', 'I', 'SiOut', 'pos', 'k')){
# 	data <- .Object@origin_data
# 	setDT(data)
# 	if('T' %in% colnames(data)){
# 		data[['Frame']] <- cut(data[['T']], seq(0, frames / 2, 0.5), labels = F)
# 		data[['T']] <- NULL
# 
# 	}
# 	data <- merge(data.table(Frame = seq_len(frames)), data,
# 		   all.x = TRUE, by = 'Frame')
# 	
# 	dt <- data.table(Frame = data[['Frame']])
# 	for(i in cols){
# 		if(i %in% colnames(data)){
# 			if(is.numeric(data[[i]])){
# 				if(is.na(data[[i]][1])){
# 					data[[i]][1] <- 0
# 				}
# 				dt[[i]] <- fillVec(data[[i]])				
# 			} else {
# 				if(i == 'pos'){
# 					data[[i]][1] <- 'nest'
# 				} else {
# 					if(is.na(data[[i]][1])){
# 						data[[i]][1] <- '0'
# 					}
# 				}
# 
# 				dt[[i]] <- fillCharVec(data[[i]])
# 			}
# 
# 		}
# 	}
# 	.Object@data <- dt
# 	.Object
# })



setMethod('activity_matrix', 'Simulation', function(.Object, binary = TRUE){
	
	if(!'node' %in% colnames(.Object@data)){
		.Object@data[, node := lapply(.SD, function(i){
			x <- parse_nodes(i)
			list(x[!is.na(x)])
		}), .SDcols = 'pos', by = 'Frame']
	}

	dt <- .Object@data[, .(node = unlist(node)), by = .(Frame)]
	dt <- dt[!is.na(node)]
	
	if(binary){
		m <- dt[, .(N = as.integer(.N > 0)), by = c('node', 'Frame')] 
	} else {
		m <- dt[, .N, by = c('node', 'Frame')]
		
	}
	.Object@matrix <- m[, c('Frame', 'node', 'N')]
	setDT(.Object@matrix)
	.Object
})

setMethod('mutual_info', 'Simulation', function(.Object, t, edges = sim_edges, nodes = NULL){
	.Object <- activity_matrix(.Object)
	m <- M2M(.Object@matrix, data = -1)
	m <- M2sbst(m, t)
	z <- numeric(nrow(edges))
	if(!is.null(nodes)){
		existing_edges <- edges[o %in% nodes & d %in% nodes]
	} else {
		existing_edges <- edges[o %in% colnames(m) & d %in% colnames(m)]
	}
	
	for(i in seq_len(nrow(existing_edges))){
		idx <- edges[o == existing_edges[i, o] & d == existing_edges[i, d], which = TRUE]
		n <- as.character(existing_edges[i, c('o', 'd')])
		z[idx] <- infotheo::mutinformation(m[, ..n])[1, 2]
	}
	z
})
