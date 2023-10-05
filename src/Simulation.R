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
hex_sim <- read.csv(Path2File, sep = ',', dec = '.', stringsAsFactors = F)
rm(Path2File)

#### +++ CLASS DEFINITION +++ ####
setClass('Simulation', representation(
	origin_data = 'data.frame',
	data = 'data.frame',
	results = 'list',
	type = 'character',
	refcoords = 'data.frame',
	r = 'numeric',
	food = 'data.frame',
	.__classVersion__ = 'character'
), prototype = list(
	r = 1.01,
	refcoords = hex_sim,
	food = data.frame(node = c("(6, 33)", "(6, 34)", "(7, 34)", # patch 1
	                           "(7, 33)", "(7, 32)", "(6, 32)",
	                           "(6, 11)", "(6, 12)", "(7, 12)", # patch 2
	                           "(7, 11)", "(7, 10)", "(6, 10)")),
	.__classVersion__ = '1.0'
))

setGeneric('updateObject', function(obj){
	standardGeneric('updateObject')
})

setMethod('updateObject', 'Simulation', function(obj){
	if(.hasSlot(obj, '.__classVersion__') && obj@.__classVersion__ == '1.0'){
		return(obj)
	} else {
		attrs <- getSlots('Simulation')
		slots <- list(Class = 'Simulation')
		for(i in names(attrs)){
			if(.hasSlot(obj, i)){
				x <- slot(obj, i)
				if(attrs[i] %in% class(x)){
					slots[[i]] <- x
				}
			}
		}
		obj <- do.call('new', args = slots)
		obj
	}
})

#### +++ GENERIC FUNCTIONS +++ ####
setGeneric('compute_connectivity', function(obj, t, ...){
	standardGeneric('compute_connectivity')
})

setGeneric('process_data', function(obj, frames = 21600,
				    cols = c('N', 'I', 'SiOut', 'pos', 'k')){
	standardGeneric('process_data')
})

#### +++ CLASS METHODS +++ ####
setMethod('compute_connectivity', signature = 'Simulation', function(obj, t, ...){

	if(!is.data.table(obj@data)){
		setDT(obj@data)
	}
	result <- obj@data[, .(N= mean(N), k = as.numeric(lapply(list(pos), function(i){
		x <- parse_nodes(i)
		mean(connectivity(x, t, ...))
	}))),
	     by = .(Frame = ((Frame - 1) %/% t)*t)]
	
	obj@results[['Frame']] <- result[['Frame']]
	obj@results[['N']] <- result[['N']]
	obj@results[['k']] <- result[['k']]
	obj

})

setMethod('process_data', signature = 'Simulation',
	  function(obj, frames = 21600, cols = c('N', 'I', 'SiOut', 'pos', 'k')){
	data <- obj@origin_data
	setDT(data)
	if('T' %in% colnames(data)){
		data[['Frame']] <- cut(data[['T']], seq(0, frames / 2, 0.5), labels = F)
		data[['T']] <- NULL

	}
	data <- merge(data.table(Frame = seq_len(frames)), data,
		   all.x = TRUE, by = 'Frame')
	
	dt <- data.table(Frame = data[['Frame']])
	for(i in cols){
		if(i %in% colnames(data)){
			if(is.numeric(data[[i]])){
				if(is.na(data[[i]][1])){
					data[[i]][1] <- 0
				}
				dt[[i]] <- fillVec(data[[i]])				
			} else {
				if(i == 'pos'){
					data[[i]][1] <- 'nest'
				} else {
					if(is.na(data[[i]][1])){
						data[[i]][1] <- '0'
					}
				}

				dt[[i]] <- fillCharVec(data[[i]])
			}

		}
	}
	obj@data <- dt
	obj
})
