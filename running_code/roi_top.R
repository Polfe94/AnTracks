if(!exists('det')){
        if(exists('sto')){
                det <- list(list(segments = sto[[1]]$segments))
        } else {
                stop('Could not find a valid "segments" object.')
        }
}

sgm <- det[[1]]$segments
food <- det[[1]]$food
centers <- rbind(food[[1]][5, 1:2], food[[2]][5, 1:2])
centers$y <- centers$y - 50

nest <- get_nest()
nest[2] <- nest[2] + 100
centers <- rbind(centers, nest)

idx_1 <- inRadius(coords = sgm[, 1:2], center = c(centers[1, ], recursive = T), r = 150)
idx_2 <- inRadius(coords = sgm[, 1:2], center = c(centers[2, ], recursive = T), r = 150)
idx_3 <- inRadius(coords = sgm[, 1:2], center = nest, r = 150)

regions <- sgm[idx_1 | idx_2 | idx_3, c('x', 'y')]
regions$region <- c(rep('food_patch_1', sum(idx_1)),rep('food_patch_2', sum(idx_2)),rep('nest', sum(idx_3)))

if(exists('do_plot') && do_plot == TRUE){
     print(
          suppressMessages(
                  draw_FoodPatches(det[[1]], add = draw_hexagons(det[[1]]), fill = 'grey50')+
                          geom_point(data = sgm[idx_1 | idx_2 | idx_3, ], aes(x, y))+
                          geom_circle(centers, 150, color = 'blue') + 
                          geom_point(data = data.frame(x = nest[1], y = nest[2]-100),
                                     aes(x, y), color = 'red', size = 3)
                  
          ))
}

rm(list = c(ls()[grepl('idx_', ls())], 'sgm', 'food', 'centers', 'nest'))
