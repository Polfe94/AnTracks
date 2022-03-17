source('G:/research/MoveAnts/code/config.R')
source('C:/Users/POL/Desktop/gants_current/code/dev/dev_functions.R')

hex <- hex90
colnames(hex) <- c('x', 'y')
h <- hex90[hex90$rotY > 980, ]
colnames(h) <- c('x', 'y')
source('G:/research/2022/AnTracks/src/generic.R')
source('G:/research/2022/AnTracks/src/coords.R')
det <- lapply(expL$det, function(i) init(i, refcoords = h,
                                      class = 'coords'))

d <- exp_spreadsheet[exp_spreadsheet$CONDITION == 'DET', c('Date', 'MATI.TARDA')][1:10, ]
d <- paste(format(strptime(d[, 1], '%d/%m/%Y'), format = '%Y%m%d'), d[, 2], sep='')

det <- lapply(seq_along(det), function(i){
        det[[i]]$date <- d[i]
        det[[i]]
})

s <- compute_segments(det[[1]])
det <- lapply(det, function(i){
        c <- class(i)
        x <- append(i, list(segments = s))
        class(x) <- c
        x
        })


# visualization
draw_hexagons(det[[1]])
food <- get.exp_foodLocation(det[1])[[1]][1:2]
draw_foodpatches(food, q = draw_hexagons(det[[1]]), a = 0.8, col = 'black', fill = 'grey20')

food <- get.exp_foodLocation(det)

det <- lapply(seq_along(det), function(i){
        
        det[[i]]$food <- food[[i]][1:2]
        det[[i]]
})


det <- lapply(det, function(i){
        p = node_idx(i, seq_len(nrow(i$data)))
        i$data$node <- p
        i
})

tt <- Sys.time()
det <- lapply(det, function(i){
        c <- class(i)
        z <- local_cov(i)
        i$local_cov <- z
        class(i) <- c
        i
})
Sys.time() - tt
beepr::beep(3)


tt <- Sys.time()
det <- lapply(det, function(i){
        i$N <- get_N(i)
        i
})
Sys.time() -tt
beepr::beep(3)


tt <- Sys.time()
det <- lapply(det, function(i){
        i$I <- get_I(i)
        i
})
Sys.time() -tt
beepr::beep(3)
