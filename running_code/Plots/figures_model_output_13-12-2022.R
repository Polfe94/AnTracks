# library(data.table)
# library(ggpubr)
# library(ggplot2)
source('~/research/2022/ANTS/AnTracks/src/Experiment.R')
path <- '~/research/2022/ANTS/AnTracks/plots/xerrada_13_desembre/'
load('~/research/2022/ANTS/AnTracks/data/det.RData')
files <- list.files(path)

files_N <- files[grepl('N.csv', files)]
files_Z <- files[grepl('xyz.csv', files)]

N <- vector('list', length(files_N))
Z <- vector('list', length(files_Z))

fillfun <- function(x, t){
     i <- 2
     time <- t[i]
     dt <- diff(t)
     
     N <- numeric(max(t))
     
     while(time < max(t)){
          N[time:(time + dt[i] - 1)] <- x[i]
          
          time <- time + dt[i]
          i <- i + 1
     }
     N[length(N)] <- x[length(x)]
     N
}

for(i in seq_along(files_N)){
     N[[i]] <- unique(read.csv(paste0(path, files_N[i]))[,-1])
     N[[i]] <- data.table(N = fillfun(N[[i]]$N, N[[i]]$`T`), t = seq(max(N[[i]]$`T`)))
     Z[[i]] <- data.frame(read.csv(paste0(path, files_Z[i]))[, -1])
}

## GENERAL PLOT TO SEE EACH INDIVIDUAL SIMULATION
plot_list <- lapply(N, function(i){
     ggplot(data = i, aes(t / 120, N)) + geom_line()
})
ggarrange(plotlist = plot_list)

## SOME FILTERING OF THE SIMULATIONS
# N <- N[lapply(N, function(i) max(i$N)) > 25]
N <- N[lapply(N, function(i) max(i$N)) > 20]
# N <- N[lapply(N, function(i) i$t[which.max(i$N)] < 10000) == T]

## AVERAGE DENSITY IN NODES [SIMULATION]
z <- do.call('rbind', lapply(Z, function(i) i$z))
z <- apply(z, 2, mean)
xyz <- data.frame(x = Z[[1]]$x, y = Z[[1]]$y, z = z, node = 1:length(Z[[1]]$x))

## SIMULATION PLOT [DENSITY]
edges <- compute_edges(xyz[xyz$y > 0, ], r = 1.01)
sim_xyz <- draw_hexagons(edges, add = ggplot(data = xyz[xyz$y > 0, ], aes(x, y)) + 
                   geom_polygon(data = data.frame(x = c(-28.58, -29.44, -29.44, -28.58, -27.71, -27.71),
                                                  y = c(9, 9.5, 10.5, 11, 10.5, 9.5)), aes(x, y), fill = 'grey20')+
                   geom_polygon(data = data.frame(x = c(-9.53, -10.39, -10.39, -9.53, -8.66, -8.66),
                                                  y = c(9, 9.5, 10.5, 11, 10.5, 9.5)), aes(x, y), fill = 'grey20')+
                   geom_point(size = 5, aes(color = z), show.legend = F) + scale_color_viridis() , color = 'grey20')

## AVERAGE DENSITY IN NODES [DETERMINIST EXPERIMENTS]
det_avg_xyz <- lapply(det, function(i){
     M <- get_matrix(i, data = 0)
     node_act <- colSums(M)
     node_act[node_act > quantile(node_act, probs = 0.95)] <- quantile(node_act, probs = 0.95)
     data.frame(x = hex[hex$node %in% as.integer(colnames(M)), 'x'],
     y = hex[hex$node %in% as.integer(colnames(M)), 'y'],
     z = node_act)
})
det_xyz <- data.table(do.call('rbind', det_avg_xyz))
det_xyz <- det_xyz[, lapply(.SD, mean), .SDcols = 'z', by = c('x', 'y')]

## INDIVIDUAL EXPERIMENT
# M <- get_matrix(det[[10]], data = 0)
# node_act <- colSums(M)
# node_act[node_act > quantile(node_act, probs = 0.95)] <- quantile(node_act, probs = 0.95)

## EXPERIMENTS PLOT [DENSITY]
exp_xyz <- draw_hexagons(add = draw_FoodPatches(det[[1]])+ 
                              geom_point(data = det_xyz, aes(x, y, color = z), size = 5, show.legend = F) +
                              scale_color_viridis()) 

# exp_xyz <- draw_hexagons(add = draw_FoodPatches(det[[2]])+ 
#                    geom_point(data = rbind(data.frame(x = hex[hex$node %in% as.integer(colnames(M)), 'x'],
#                                                       y = hex[hex$node %in% as.integer(colnames(M)), 'y'],
#                                                       z = node_act),
#                                            data.frame(x = hex[hex$y > 1000 & !hex$node %in% as.integer(colnames(M)), 'x'],
#                                                       y = hex[hex$y > 1000 & !hex$node %in% as.integer(colnames(M)), 'y'],
#                                                       z = 0)), aes(x, y, color = z), size = 5, show.legend = F) +
#                    scale_color_viridis()) 

ggarrange(sim_xyz, exp_xyz, ncol = 1, labels = c('Simulation', 'Experiment'), vjust = 1.3)
# 
# ggplot(data = xyz, aes(x, y)) + 
#      geom_polygon(data = data.frame(x = c(-28.58, -29.44, -29.44, -28.58, -27.71, -27.71),
#                                     y = c(9, 9.5, 10.5, 11, 10.5, 9.5)), aes(x, y), fill = 'grey20')+
#      geom_polygon(data = data.frame(x = c(-9.53, -10.39, -10.39, -9.53, -8.66, -8.66),
#                                     y = c(9, 9.5, 10.5, 11, 10.5, 9.5)), aes(x, y), fill = 'grey20')+
#      geom_point(size = 5, aes(color = z), show.legend = F) + scale_color_viridis() + 
#      theme_void()

## TEMPORAL DYNAMICS [SIMULATION SUBSET]
sbst_sim <- c(3, 4, 5, 9, 12, 13, 14, 16, 27)
sbst_sim <- c(27, 9, 13, 14, 3)
# sbst_sim <- c(6, 8, 9, 29)
sbst_exp <- c(1, 2, 5, 8, 9)

list2plot <- vector('list', length = length(sbst_exp) + length(sbst_sim))
for(i in seq(length(sbst_exp) + length(sbst_sim))){
        if(i <= length(sbst_exp)){

                list2plot[[i]] <- data.frame(t = seq(11, 21589, 6), 
                                             N = moving_average(det_N[sbst_exp[i], 1:21599], 20, 5),
                                             exp = paste('Experiment', i), cond = 'Experiment')

        } else {
                k <- i - length(sbst_exp)
                list2plot[[i]] <- data.frame(t = seq(11, 21589, 6), 
                                             N = moving_average(N[[sbst_sim[k]]]$N[1:21599], 20, 5),
                                             exp = paste('Simulation', k), cond = 'Simulation')

        }
        
}
dfsubsetted <- do.call('rbind', list2plot)
ggplot(data = dfsubsetted, aes(t/120, N,color = cond, size = cond), alpha = 0.8) + geom_line() +
        scale_size_manual(values = c(0.75, 1.25), guide = 'none')+
        facet_wrap(~ exp, ncol = 5, scales = 'free_y')+
        scale_color_manual('',values = c('grey10', 'grey60'))+
        scale_x_continuous('Time (min)', breaks = seq(0, 180, 30))+
        scale_y_continuous('Activity', breaks = seq(0, 40, 7.5))+
        theme_bw() + 
        theme(axis.title = element_text(size = 22),
              axis.text = element_text(size = 18),
              legend.text = element_text(size = 22),
              strip.text = element_text(size = 22))+
        guides(color = guide_legend(override.aes = list(size = 3)))


example_df <- rbind(
        # SIMULATION
        data.frame(t = 1:21599, N = moving_average(N[[9]]$N, 20, 10), exp = 'Simulation'),
        # EXPERIMENT
        data.frame(t = 1:21599, N = moving_average(det_N[7, 1:21599], 20, 10), exp = 'Experiment')
)


## TEMPORAL DYNAMICS [SIMULATION AVERAGE]
Nt <- data.table(do.call('rbind', N))
Nt <- Nt[, c(sd = lapply(.SD, sd), N = lapply(.SD, mean)), by = `t`]
colnames(Nt) <- c('t', 'sd', 'N')
Nt$exp <- 'Simulations average'

## TEMPORAL DYNAMICS [EXPERIMENTS AVERAGE]
det_N <- do.call('rbind', lapply(det, function(i){
     # ifelse(length(i@N) > 21600, T, F)
     if(length(i@N) > 21600){
          i@N[1:21600]
     } else {
          NULL
     }
}))
avg_N <- data.frame(t = 1:21600, N = apply(det_N, 2, mean),
                    sd = apply(det_N, 2, sd), 
                    exp = 'Experiments average')
avg_df <- rbind(Nt, avg_N)

# Average
ggplot(data = avg_df[avg_df$t < 21550], aes(t / 120, N)) + 
     geom_ribbon(aes(ymin = N - sd, ymax = N + sd, fill = exp), alpha = 0.8, show.legend = F)+
     # geom_line(aes(y = N - sd), color = 'grey10')+
     # geom_line(aes(y = N + sd), color = 'grey10')+
     geom_line(size = 1, aes(color = factor(exp))) +
     scale_color_manual('', values = c('grey25', 'grey40'))+
     scale_fill_manual(values = c('grey60', 'grey80'))+
     scale_x_continuous('Time (min)', breaks = seq(0, 180, 15))+
     scale_y_continuous('Activity', breaks = seq(0, 40, 7.5))+
     theme_bw() + 
     theme(axis.title = element_text(size = 25),
           axis.text = element_text(size = 22),
           legend.text = element_text(size = 22))+
     guides(color = guide_legend(override.aes = list(size = 3)))


## INDIVIDUAL EXPERIMENT COMPARISON 
example_df <- rbind(
     # SIMULATION
     data.frame(t = 1:21599, N = moving_average(N[[9]]$N, 20, 10), exp = 'Simulation'),
     # EXPERIMENT
     data.frame(t = 1:21599, N = moving_average(det_N[7, 1:21599], 20, 10), exp = 'Experiment')
     )

# Example
ggplot(data = example_df, aes(t / 120, N)) + 
     geom_line(size = 1, aes(color = factor(exp))) +
     scale_color_manual('', values = c('black', 'grey50'))+
     
     scale_x_continuous('Time (min)', breaks = seq(0, 180, 15))+
     scale_y_continuous('Activity', breaks = seq(0, 40, 7.5))+
     theme_bw() + 
     theme(axis.title = element_text(size = 25),
           axis.text = element_text(size = 22),
           legend.text = element_text(size = 22))+
     guides(color = guide_legend(override.aes = list(size = 3)))





## TESTS ##
# ggplot(data = Nt, aes(t / 120, N)) + 
#      geom_ribbon(aes(ymin = N - sd, ymax = N + sd), fill = 'grey80')+
#      geom_line(aes(y = N - sd), color = 'grey25')+
#      geom_line(aes(y = N + sd), color = 'grey25')+
#      geom_line(size = 1) + 
#      scale_x_continuous('Time (min)', breaks = seq(0, 180, 15))+
#      geom_line(data = data.frame(x = seq_along(avg_N) / 120,
#                                  y =avg_N), aes(x, y), color = 'red')+
#      scale_y_continuous('Activity', breaks = seq(0, 30, 5))+
#      theme_bw() + 
#      theme(axis.title = element_text(size = 25),
#            axis.)

# ggplot(data = Nt, aes(t / 120, N)) + 
#      geom_ribbon(aes(ymin = N - sd, ymax = N + sd), fill = 'grey80')+
#      geom_line(aes(y = N - sd), color = 'grey25')+
#      geom_line(aes(y = N + sd), color = 'grey25')+
#      geom_line(size = 1) + 
#      scale_x_continuous('Time (min)', breaks = seq(0, 180, 15))+
#      geom_line(data = data.frame(x = seq_along(det[[4]]@N) / 120,
#                                  y = det[[4]]@N), aes(x, y), color = 'red')
# 
# ggplot(data = Nt, aes(t / 120, N)) + 
#      geom_ribbon(aes(ymin = N - sd, ymax = N + sd), fill = 'grey80')+
#      geom_line(aes(y = N - sd), color = 'grey25')+
#      geom_line(aes(y = N + sd), color = 'grey25')+
#      geom_line(size = 1) + 
#      scale_x_continuous('Time (min)', breaks = seq(0, 180, 15))+
#      geom_line(data = data.frame(x = seq_along(det[[8]]@N) / 120,
#                                  y = det[[8]]@N), aes(x, y), color = 'red')
# 
# 
# ggplot(data = Nt, aes(t / 120, N)) + 
#      geom_ribbon(aes(ymin = N - sd, ymax = N + sd), fill = 'grey80')+
#      geom_line(aes(y = N - sd), color = 'grey25')+
#      geom_line(aes(y = N + sd), color = 'grey25')+
#      geom_line(size = 1) + 
#      scale_x_continuous('Time (min)', breaks = seq(0, 180, 15))+
#      geom_line(data = data.frame(x = seq_along(det[[6]]@N) / 120,
#                                  y = det[[6]]@N), aes(x, y), color = 'red')
