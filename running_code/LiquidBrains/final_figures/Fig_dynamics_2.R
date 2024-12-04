#### LIBRARIES, DATA AND GENERIC FUNCTIONS ####
library(dtw)
library(dtwclust)
library(latex2exp)
library(arrow)
# library(infotheo)
source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')

load('~/research/gits/AnTracks/data/det.RData')


## EXPERIMENTS
foods_det <- rbindlist(lapply(det, function(i){
        a <- rbindlist(i@food)[['t']]
        data.table(mint = min(a)/2, maxt = max(a)/2)
}))[, lapply(.SD, median)] # <--- median !

det_avg <- rbindlist(lapply(det, function(i){
        setDT(i@data)
        x <- merge(i@data[, .N, by = 'Frame'], data.table(Frame = 1:21600), all = TRUE, by = 'Frame')
        if(i@date != det[[7]]@date){
                x[is.na(N), 'N'] <- 0		
        }
        x
}))[, .(N = mean(N, na.rm = T)), by = 'Frame']

det_N <- rbindlist(lapply(det, function(i){
        setDT(i@data)
        x <- merge(i@data[, .N, by = 'Frame'], data.table(Frame = 1:21600), all = TRUE, by = 'Frame')
        if(i@date != det[[7]]@date){
                x[is.na(N), 'N'] <- 0		
        }
        x[, .(N = movingAverage(N, t = 60)[1:695], Frame = seq(31, 21545, 31))]
}), idcol = TRUE)


process_data <- function(path, 
                         vars = c('Frame', 'N')){
        ref <- data.table(Frame = 1:21600)
        t <- movingAverage(ref[['Frame']], 60, 0)
        
        files <- list.files(path)
        files <- files[!grepl('food', files) & !grepl('position', files) &
                               !grepl('data', files) & !grepl('keys', files)]
        
        a <- lapply(files, function(i){
                
                y <- merge(ref, data.table(read_parquet(paste0(path, i)))[, ..vars], all = TRUE, by = 'Frame')
                y[['N']] <- fillVec(y[['N']])
                n <- movingAverage(y[['N']], 60L, 0L)
                z <- zscore(n)
                data.table(Frame = t, N = n, Z = z)
                
        })
        names(a) <- files
        a
}

path <- '~/research/gits/AutomatAnts/results/rho_0.4_eps_1/'
files <- list.files(path)
files <- unique(unlist(regmatches(files, gregexpr('rho_0.401_epsilon_1.001_\\d{1,2}', files))))


#### SPATIAL DYNAMICS
food <- data.table(read_parquet(paste0(path, files[1], '_food.parquet')))[, node := revert_node(node)]
food <- as.data.frame(merge(food, hex_sim, by = 'node', sort = FALSE)[, c('x', 'y')])
food <- list(GP1 = food[1:6, ], GP2 = food[7:12, ])

pos_data <- rbindlist(lapply(files, function(i){ ## <--- takes a couple of minutes and a bit of RAM
        do.call('gc', args = list(verbose = FALSE))
        data.table(read_parquet(paste0(path, i, '_data.parquet')))[, .(node = parse_nodes(pos))][, .N, by = 'node']
}))

sum_pos_data <- merge(hex_sim, pos_data[, .(N = sum(N)), by = 'node'], by = 'node')
c <- draw_nodes(sum_pos_data, 
                add = draw_hexagons(sim_edges, linewidth = 2, color = 'black',
                                    add = draw_hexagons(sim_edges, linewidth = 2, color = 'black', 
                                                        add = geom_foodpatches(food = food, fill = 'mediumpurple', alpha = 0.9))), 
                z = rank(sum_pos_data[['N']]), size = 3.54, show.legend = FALSE) + 
        theme(aspect.ratio = 0.5 , plot.margin = margin(t = -5, b = -70), plot.subtitle = element_text(size = 14))+
        scale_fill_viridis(option = 'plasma') +
        ggtitle('','Model')

pos_det <- merge(hex, rbindlist(lapply(det, function(i){
        i@data[, .N, by = 'node']
}), idcol = TRUE)[, .(N = sum(N)), by = 'node'], by = 'node')

c2 <- draw_nodes(pos_det, 
                 add = draw_hexagons(edges, linewidth = 2, color = 'black',
                                     add = draw_hexagons(edges, linewidth = 2, color = 'black', 
                                                         add = geom_foodpatches(fill = 'mediumpurple', alpha = 0.9))), 
                 z = rank(pos_det[['N']]), size = 3.5) + 
        scale_fill_viridis('Density',option = 'plasma', breaks = c(1, 620),
                           labels = c('Low', 'High')) +
        theme(legend.position = 'bottom', aspect.ratio = 0.5,
              plot.margin = margin(t = 0, b = -90), plot.subtitle = element_text(size = 14))+
        guides(fill = guide_colorbar(title.position = 'top', barwidth = 15, title.hjust = 0.5))+
        ggtitle('', 'Experimental')


#### POPULATION DYNAMICS ####

data <- process_data(path)

full_data <- rbindlist(lapply(files, function(i){
        data.table(read_parquet(paste0(path, i,'.parquet')))
}), idcol = TRUE)[order(Frame)]
full_data_avg <- full_data[Frame < 175*120, .(N = mean(N)), by = 'Frame'] ## filtered

foods <- rbindlist(lapply(files, function(i){
        food <- data.table(read_parquet(paste0(path, i,'_food.parquet')))[['t']]
        data <- data.table(read_parquet(paste0(path, i,'.parquet')))
        mint <- data[N > 0, min(Frame)] /2
        r <- range(food) - mint
        data.table(mint = r[1], maxt = r[2])
        
}), idcol = TRUE)[is.finite(mint),lapply(.SD, median), .SDcols = c('mint', 'maxt')]


data_peak <- rbindlist(data, idcol = TRUE)

# comparison of population dynamics
a <- ggplot(data = data_peak[Frame < 175*120], aes(Frame, N)) +
        geom_path(aes(group = .id), color = 'grey80', alpha = 0.1, linewidth = 0.8)+
        geom_path(data = det_N[Frame < 175*120], aes(group = .id), color = 'grey30', alpha = 0.5, linewidth = 0.8)+
        
        geom_path(data = full_data_avg, color = 'black', linewidth = 3, alpha = 0.5)+
        geom_path(data = full_data_avg, aes(color = 'Simulations'), linewidth = 2, alpha = 0.5)+
        geom_path(data = det_avg[Frame < 175*120], color = 'black', linewidth = 3, alpha = 0.5)+
        geom_path(data = det_avg[Frame < 175*120], aes(color = 'Experiments'), linewidth = 2, alpha = 0.5)+
        scale_x_continuous('Time (min)', breaks = seq(0, 150, 50)*120, labels = seq(0, 150, 50))+
        scale_y_continuous('N (number of ants in arena)')+
        scale_color_manual('',values = c('grey30', 'grey80'), labels = c('Experimental', 'Model'))+
        geom_vline(xintercept = unlist(foods)*2, linetype = 2, linewidth = 1)+
        geom_vline(xintercept = unlist(foods_det)*2, linetype = 3, linewidth = 1)+
        theme(legend.position = c(0.8, 0.85),
              legend.background = element_rect(fill = 'white', color = 'black'),
              legend.title = element_blank())+
        geom_text(data = data.frame(x = c(120, 51*120), y = 40, label = c('TP1', 'TP2')), 
                 aes(x, y, label = label), size = 10)+
        ggtitle('')


png('~/research/gits/AnTracks/plots/liquid_brains/Fig_dynamics.png', 8000, 4000, res = 600)
ggarrange(a + theme(aspect.ratio = 0.6),
          ggarrange(c2, c, ncol = 1, common.legend = TRUE, legend = 'bottom'), widths = c(0.7, 0.3),
          labels = c('A', 'B'), font.label = list(face = 'plain', size = 22), heights = c(0.6, 0.4))
# grid.arrange(a + theme(aspect.ratio = 0.6), 
#           ggarrange(c2, c, ncol = 1, common.legend = TRUE, legend = 'bottom'), widths = c(0.7, 0.3), 
#           labels = c('A', 'B'), font.label = list(face = 'plain', size = 22), heights = c(0.6, 0.4))
dev.off()

