source('G:/research/2022/AnTracks/src/generic.R')
source('G:/research/2022/AnTracks/src/coords.R')
source('G:/research/2022/AnTracks/src/nodes.R')
load('G:/research/2022/AnTracks/results/det_coords.RData')
load('G:/research/2022/AnTracks/results/sto_coords.RData')

nest <- get_nest()


### +++ STATIONARY VIEW

# [[DETERMINIST]]
example_det <- coords2matrix(det[[8]])
example_det <- example_det[, colnames(example_det) %in% rownames(hex[hex$y > 1000, ])]
example_det[example_det == -1] <- 0

# avg_det <- apply(example_det, 2, mean)
avg_det <- apply(example_det, 2, sum)
nodes_det <- cbind(hex[colnames(example_det), ] , node = colnames(example_det))
n_det <- sum(example_det)

# [[STOCHASTIC]]
example_sto <- coords2matrix(sto[[8]])
example_sto <- example_sto[, colnames(example_sto) %in% rownames(hex[hex$y > 1000, ])]
example_sto[example_sto == -1] <- 0

# avg_sto <- apply(example_sto, 2, mean)
avg_sto <- apply(example_sto, 2, sum)
nodes_sto <- cbind(hex[colnames(example_sto), ] , node = colnames(example_sto))
n_sto <- sum(example_sto)



results <- list(det = c(), sto = c())
sq <- seq(0, 750, 25)

for(i in seq_along(sq)){
        inbound_det <- inRadius(coords = nodes_det[, 1:2], center = nest, r = sq[i])
        inb_det <- sum(avg_det[inbound_det])
        inbound_sto <- inRadius(coords = nodes_sto[, 1:2], center = nest, r = sq[i])
        inb_sto <- sum(avg_sto[inbound_sto])
        
        results[['det']] <- rbind(results[['det']], 
                                  data.frame(# inB = mean(avg_det[inbound_det]),
                                          inB = inb_det,
                                          outB = n_det - inb_det,
                                          # outB = mean(avg_det[!inbound_det]),
                                          r = sq[i]))
        
        results[['sto']] <- rbind(results[['sto']], 
                                  data.frame(# inB = mean(avg_sto[inbound_sto]),
                                          # outB = mean(avg_sto[!inbound_sto]),
                                          inB = inb_sto,
                                          outB = n_sto - inb_sto,
                                          r = sq[i]))
}

# 450 < r < 475
mltdet <- reshape2::melt(results$det, id.vars = 'r')
pl1 <- ggplot(data = mltdet, aes(r, value, color = variable)) +
        geom_vline(xintercept = 455, linetype = 'dashed')+
        geom_hline(yintercept = 141500, linetype = 'dashed')+
        geom_path(size = 2)+

        scale_x_continuous('Radius size (mm)', breaks = seq(0, 1000, 75))+
        scale_y_continuous('Activity', breaks = seq(0, 250000, 50000))+
        scale_color_viridis_d('', begin = 0.2, end = 0.9, labels = c('Inside radius', "Outside radius"))+
        ggtitle('Information about food location')



# 575 < r < 600
mltsto <- reshape2::melt(results$sto, id.vars = 'r')
pl2 <- ggplot(data = mltsto, aes(r, value, color = variable)) +
        geom_vline(xintercept = 575, linetype = 'dashed')+
        geom_hline(yintercept = 160750, linetype = 'dashed')+
        geom_path(size = 2)+
        scale_x_continuous('Radius size (mm)', breaks = seq(0, 1000, 75))+
        scale_y_continuous('Activity', breaks = seq(0, 350000, 50000))+
        scale_color_viridis_d('', begin = 0.2, end = 0.9, labels = c('Inside radius', "Outside radius"))+
        ggtitle('No nformation about location')

ggpubr::ggarrange(pl1, pl2, ncol = 2, common.legend = T, legend = 'bottom')



### +++ DYNAMIC VIEW

results_dynamic <- list(det = c(), sto = c())
sq <- seq(1, min(nrow(example_det), nrow(example_sto)), 50)

for(i in seq_along(sq)){
        inbound_det <- inRadius(coords = nodes_det[, 1:2], center = nest, r = 600)
        inb_det = sum(example_det[1:sq[i], inbound_det])
        n_det <- sum(example_det[1:sq[i], ])
        inbound_sto <- inRadius(coords = nodes_sto[, 1:2], center = nest, r = 600)
        inb_sto = sum(example_sto[1:sq[i], inbound_sto])
        n_sto <- sum(example_sto[1:sq[i], ])
        results_dynamic[['det']] <- rbind(results_dynamic[['det']], 
                                          data.frame(# inB = mean(avg_det[inbound_det]),
                                                  inB = inb_det,
                                                  outB = n_det - inb_det,
                                                  # outB = mean(avg_det[!inbound_det]),
                                                  t = sq[i]))
        
        results_dynamic[['sto']] <- rbind(results_dynamic[['sto']], 
                                          data.frame(# inB = mean(avg_sto[inbound_sto]),
                                                  # outB = mean(avg_sto[!inbound_sto]),
                                                  inB = inb_sto,
                                                  outB = n_sto - inb_sto,
                                                  t = sq[i]))
}

# 450 < r < 475
mltdet <- reshape2::melt(results_dynamic$det, id.vars = 't')
ggplot(data = mltdet, aes(t, value, color = variable)) + geom_path()
ggplot(data = results_dynamic$det, aes(t, inB - outB)) + geom_path()

# 575 < r < 600
mltsto <- reshape2::melt(results_dynamic$sto, id.vars = 't')
ggplot(data = mltsto, aes(t, value, color = variable)) + geom_path()
ggplot(data = results_dynamic$sto, aes(t, inB - outB)) + geom_path()


### +++ FOR DIFERENT R

diff_radius <- lapply(seq(100, 750, 150), function(x){
        det <- c()
        sto <- c()
        for(i in seq_along(sq)){
                inbound_det <- inRadius(coords = nodes_det[, 1:2], center = nest, r = x)
                inb_det = sum(example_det[1:sq[i], inbound_det])
                n_det <- sum(example_det[1:sq[i], ])
                inbound_sto <- inRadius(coords = nodes_sto[, 1:2], center = nest, r = x)
                inb_sto = sum(example_sto[1:sq[i], inbound_sto])
                n_sto <- sum(example_sto[1:sq[i], ])
                det <- rbind(det, data.frame(dif = 2 * inb_det - n_det,
                                                          t = sq[i]))
                sto <- rbind(sto, data.frame(dif = 2 * inb_sto - n_sto,
                                             t = sq[i]))
                
        }
        list(det = det, sto = sto, r = x)
})

df_summary <- do.call('rbind', lapply(diff_radius, function(i){
        cbind.data.frame(t = i[[1]][, 2], det = i[[1]][, 1], sto = i[[2]][, 1], r = i[[3]])}))
df_summary <- do.call('rbind', lapply(diff_radius, function(i){
        cbind.data.frame(t = i[[1]][, 2], det = i[[1]][, 1], sto = i[[2]][, 1], r = i[[3]])}))
mlt_summary <- reshape2::melt(df_summary, id.vars = c('t', 'r'))

ggplot(data = df_summary, aes(t, det)) + geom_path() + facet_wrap(~r)
ggplot(data = df_summary, aes(t, sto)) + geom_path() + facet_wrap(~r)

ggplot(data = mlt_summary, aes(t/120, value, color = variable)) + geom_path(size = 1.5) +
        ylab('Activity inside radius - Activity outside radius')+
        facet_wrap(~r, scales = 'free')+
        scale_color_viridis_d('Information', labels = c('Yes', 'No'), begin = 0.2, end = 0.9)+
        scale_x_continuous('Time (min)')



inset <- draw_hexagons(det[[1]])
inset + 
        geom_point(data = hex[hex$y > 1000, ], aes(x, y), shape = 21, fill = 'white')+
        geom_circle(data.frame(x = nest[1], y = nest[2]), r = 100, color = 'blue', size = 1.5)+
        geom_circle(data.frame(x = nest[1], y = nest[2]), r = 250, color = 'blue', size = 1.5)+
        geom_circle(data.frame(x = nest[1], y = nest[2]), r = 400, color = 'blue', size = 1.5)+
        geom_circle(data.frame(x = nest[1], y = nest[2]), r = 550, color = 'blue', size = 1.5)+
        geom_circle(data.frame(x = nest[1], y = nest[2]), r = 700, color = 'blue', size = 1.5)+
        scale_y_continuous(limits = c(990, 2000))


ggplot(data = data.frame(t = seq(0, 21600, 150)/ 120, I = sapply(MI_det150, mean)), aes(t, I))+
        geom_point()+
        scale_x_continuous('Time (min)', breaks = seq(0, 185, 15))+ ylab('Neighbouring Mutual Information')+
        geom_polygon(data = data.frame(x = c(min(do.call('rbind', det[[9]]$food)$t) / 120, 
                                       max(do.call('rbind', det[[9]]$food)$t) / 120,
                                       max(do.call('rbind', det[[9]]$food)$t) / 120,
                                      min(do.call('rbind', det[[9]]$food)$t) / 120),
                                      y = c(0, 0, 0.0043, 0.0043)),aes(x, y),
                     fill = 'lightgreen', alpha = 0.2, color = 'lightgreen', linetype = 'dashed')+
        geom_line()
        


draw_hexagons(det[[1]], z = apply(do.call('rbind',MI_det150[1:4]), 2, mean), size = 5,
              add = draw_hexagons(det[[1]], size = 6), show.legend = F)+
        geom_point(data = hex[hex$y > 985, ], aes(x, y), size = 5, fill = 'white', color = 'black', shape = 21)+
        scale_color_gradient2(low = muted('blue'),high = muted('red'))


draw_hexagons(det[[1]], z = apply(do.call('rbind',MI_det150[5:34]), 2, mean), size = 5,
              add = draw_hexagons(det[[1]], size = 6), show.legend = F)+
        geom_point(data = hex[hex$y > 985, ], aes(x, y), size = 5, fill = 'white', color = 'black', shape = 21)+
        scale_color_gradient2(low = muted('blue'),high = muted('red'))


draw_hexagons(det[[1]], z = apply(do.call('rbind',MI_det150[35:145]), 2, mean), size = 5,
              add = draw_hexagons(det[[1]], size = 6), show.legend = F)+
        geom_point(data = hex[hex$y > 985, ], aes(x, y), size = 5, fill = 'white', color = 'black', shape = 21)+
        scale_color_gradient2(low = muted('blue'),high = muted('red'))
