
source('G:/research/2022/AnTracks/src/generic.R')
source('G:/research/2022/AnTracks/src/coords.R')
source('G:/research/2022/AnTracks/src/nodes.R')
load('G:/research/2022/AnTracks/results/det_coords.RData')
load('G:/research/2022/AnTracks/results/sto_coords.RData')

if(file.exists("G:/research/2022/AnTracks/data/Det_I.dat")){
        d <- read.table(file = "G:/research/2022/AnTracks/data/Det_I.dat", header = F, sep = "")
} else {
        if(file.exists("~/research/2022/AnTracks/data/Det_I.dat")){
                d <- read.table(file = "~/research/2022/AnTracks/data/Det_I.dat", header = F, sep = "")
        } else {
                stop("Could not find data !")
        }
}

coord_correspondence <- cbind(d[, 1:2])
d[, 2] <- d[, 2]  + 1000
xy <- as.matrix(hex)
node <- sapply(1:nrow(d), function(i){which.min(pdist(
        as.matrix(d[i, 1:2]), xy
))})
coord_correspondence$node <- node

coord_correspondence <- cbind(coord_correspondence, xy[coord_correspondence$node, ])
colnames(coord_correspondence) <- c('jx', 'jy', 'node', 'px', 'py')
write.table(coord_correspondence, 
            file = "G:/research/2022/AnTracks/data/coord_correspondence.csv", sep = ',', dec = '.')

# ggplot(data = coord_correspondence) + geom_point(aes(jx, jy+ 1000)) + geom_point(aes(px, py), color = 'red')

det_food <- rbind(data.frame(node = vapply(1:6, function(i){
        get_node(det[[1]]$food[[1]][i, 1], det[[1]]$food[[1]][i, 2], hex)
        }, numeric(1)), patch = "p1"), 
        data.frame(node = vapply(1:6, function(i){
                get_node(det[[1]]$food[[2]][i, 1], det[[1]]$food[[2]][i, 2], hex)
        }, numeric(1)), patch = "p2"))
write.table(det_food, file = "G:/research/2022/AnTracks/data/det_food.csv", sep = ',', dec = '.')


sto_food <- c()
for(x in seq_along(sto)){
        sto_food <- rbind(sto_food, rbind(
                data.frame(node = vapply(1:6, function(i){
                        get_node(sto[[x]]$food[[1]][i, 1], sto[[x]]$food[[1]][i, 2], hex)
                }, numeric(1)), patch = "p1", exp = sto[[x]]$date), 
                data.frame(node = vapply(1:6, function(i){
                        get_node(sto[[x]]$food[[2]][i, 1], sto[[x]]$food[[2]][i, 2], hex)
                }, numeric(1)), patch = "p2", exp = sto[[x]]$date)
        ))
}
write.table(sto_food, file = "G:/research/2022/AnTracks/data/sto_food.csv", sep = ',', dec = '.')


## interaction matrices
sbst <- as.character(which(hex$y > 1000))
det_I <- lapply(seq_along(det), function(i){
        interaction_matrix(det[[i]])
})
for(i in seq_along(det_I)){
        det_I[[i]] <- det_I[[i]][, colnames(det_I[[i]]) %in% sbst]
        write.table(det_I[[i]], 
                    file = paste0("G:/research/2022/AnTracks/results/Dani_Campos/", 
                                  det[[i]]$date,".csv"), sep = ',', dec = '.')
}

sto_I <- lapply(seq_along(sto), function(i){
        interaction_matrix(sto[[i]])
})
for(i in seq_along(sto_I)){
        sto_I[[i]] <- sto_I[[i]][, colnames(sto_I[[i]]) %in% sbst]
        write.table(sto_I[[i]], 
                    file = paste0("G:/research/2022/AnTracks/results/Dani_Campos/", 
                                  sto[[i]]$date,".csv"), sep = ',', dec = '.')
}
