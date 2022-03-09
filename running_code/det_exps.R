source('G:/research/MoveAnts/code/config.R')
source('C:/Users/POL/Desktop/gants_current/code/dev/dev_functions.R')
# expL <- lapply(expL$det, function(i) i[i$Ymm > 1000, ])
h <- hex90[hex90$rotY > 980, ]
colnames(h) <- c('x', 'y')
source('G:/research/2022/AnTracks/src/generic.R')
expL <- lapply(expL$det, function(i) init(i, refcoords = h,
                                      class = 'coords'))
s <- compute_segments(expL[[1]])
expL <- lapply(expL, function(i){
        c <- class(i)
        x <- append(i, list(segments = s))
        class(x) <- c
        x
        })


# visualization
draw_hexagons(expL[[1]])
food <- get.exp_foodLocation(expL[1])[[1]][1:2]
draw_foodpatches(food, q = draw_hexagons(expL[[1]]), a = 0.8, col = 'black', fill = 'grey20')


expL <- lapply(expL, function(i){
        p = node_idx(i, seq_len(nrow(i$data)))
        i$data$node <- p
        i
})

tt <- Sys.time()
expL <- lapply(expL, function(i){
        c <- class(i)
        z <- local_cov(i)
        i$local_cov <- z
        class(i) <- c
        i
})
Sys.time() - tt
beepr::beep(3)
