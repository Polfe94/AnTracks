source('G:/research/2022/AnTracks/src/Experiment.R')


nms_det <- gsub('json', '', names(det))
dates <- strsplit(exp_spreadsheet$Date[exp_spreadsheet$CONDITION == 'STOCH'], '/')
nms_sto <- paste0(vapply(dates, function(i) paste0(i[[3]], i[[2]], i[[1]]), character(1)), c('M', 'T'))

Ldet <- lapply(1:10, function(i){
        x <- new('Experiment', data = det[[i]]$data, refcoords = hex[hex$y > 1000, ], 
                 date = nms_det[i], type = 'det')
        x@data$node <- NULL
        x
})
Lsto <- lapply(1:10, function(i){
        x <- new('Experiment', data = sto[[i]]$data, refcoords = hex[hex$y > 1000, ], 
                 date = nms_sto[i], type = 'sto')
        x@data$node <- NULL
        x
})

det <- Ldet
sto <- Lsto
rm(Ldet, Lsto)


det <- lapply(det, function(i){
        i@trails <- integer()
        refcoords <- hex[hex$y > 1000, ]
        i@trails <- refcoords$node[inRadius(refcoords[, c('x', 'y')], c(1020, 1015), 400)]
        i
})

det <- lapply(det, function(i){
        i@AiT <- data.frame()
        i@IiT <- data.frame()
        m <- as.data.frame(get_matrix(i, data = 0))
        df <- data.frame(t = as.integer(rownames(m)),
                         N = rowSums(m),
                         inTrail = rowSums(m[, colnames(m) %in% as.character(i@trails)]),
                         outTrail = rowSums(m[, !colnames(m) %in% as.character(i@trails)]))
        mlt_df <- melt(df, id.vars = 't')
        i@AiT <- mlt_df
        i <- iit(i)
        i
})

sto <- lapply(sto, function(i){
        i@trails <- integer()
        refcoords <- hex[hex$y > 1000, ]
        i@trails <- refcoords$node[inRadius(refcoords[, c('x', 'y')], c(1020, 1015), 400)]
        i
})

sto <- lapply(sto, function(i){
        i@AiT <- data.frame()
        i@IiT <- data.frame()
        m <- as.data.frame(get_matrix(i, data = 0))
        df <- data.frame(t = as.integer(rownames(m)),
                         N = rowSums(m),
                         inTrail = rowSums(m[, colnames(m) %in% as.character(i@trails)]),
                         outTrail = rowSums(m[, !colnames(m) %in% as.character(i@trails)]))
        mlt_df <- melt(df, id.vars = 't')
        i@AiT <- mlt_df
        i <- iit(i)
        i
})

