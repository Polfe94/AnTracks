exp_condition <- c('')

source('~/research/gAnts/code/config.R')
source('~/gants_current/code/dev/dev_functions.R')

hex <- hex90
colnames(hex) <- c('x', 'y')
h <- hex90[hex90$rotY > 980, ]
colnames(h) <- c('x', 'y')
source('~/research/2022/ANTS/AnTracks/src/generic.R')
source('~/research/2022/ANTS/AnTracks/src/coords.R')

load('/home/polfer/research/2022/ANTS/AnTracks/results/sto_coords.RData')
load('/home/polfer/research/2022/ANTS/AnTracks/results/det_coords.RData')

rm(list = ls()[grepl('json', ls())], visual, colony, current_dir, dirL, exp_condition, exps, conditionL)


#######                       #######
## +++ DETERMINIST EXPERIMENTS +++ ##
#######                       #######

det_phase <- c()

for(i in seq_along(det)){
        x <- food_detection(det[[i]])
        det[[i]]$food$GP1$t <- x[1:6]
        det[[i]]$food$GP2$t <- x[7:12]
        det_phase <- rbind(det_phase, c(min(x, na.rm = T), max(x, na.rm = T)))
}

for(i in seq_along(det)){
        print(paste0('Iter ',i))
        det[[i]]$mutual_info_p1 <- mutual_info(det[[i]], 1:det_phase[i, 1])
        det[[i]]$mutual_info_p2 <- mutual_info(det[[i]], (det_phase[i, 1] + 1):det_phase[i, 2])
        det[[i]]$mutual_info_p3 <- mutual_info(det[[i]], (det_phase[i, 2]+1):max(det[[i]]$data$Frame))
}

dp1 <- c()
dp2 <- c()
dp3 <- c()

for(i in seq_along(det)){
        dp1 <- rbind(dp1, det[[i]]$mutual_info_p1)
        dp2 <- rbind(dp2, det[[i]]$mutual_info_p2)
        dp3 <- rbind(dp3, det[[i]]$mutual_info_p3)
}


plot_dp1 <- draw_hexagons(det[[1]], z = colMeans(dp1, na.rm = T), size = 1, 
                          add = draw_hexagons(det[[1]], size = 1.05,
                                              add = draw_FoodPatches(det[[1]], fill = 'grey70')))+
        geom_point(data = hex[634, ], aes(x, y-50), size = 2, color = 'grey40', shape = 17)+
        theme_void()+
        theme(aspect.ratio = 0.5, legend.position = 'none', plot.title = element_text(hjust = 0.5),
              axis.title.y = element_text(angle = 90, margin = unit(c(0, 0, 0, 0), 'pt'),
                                          size = 9),
              plot.margin = unit(c(0, 0, -140,0), 'pt')) + ggtitle('DET')+ ylab('Exploration')+
        scale_color_gradient(low = "white", high = "#4C0000")

plot_dp2 <- draw_hexagons(det[[1]], z = colMeans(dp2, na.rm = T), size = 1, 
                          add = draw_hexagons(det[[1]], size = 1.05,
                                              add = draw_FoodPatches(det[[1]], fill = 'grey70')))+
        geom_point(data = hex[634, ], aes(x, y-50), size = 2, color = 'grey40', shape = 17)+
        theme_void()+
        theme(aspect.ratio = 0.5, legend.position = 'none', plot.title = element_text(hjust = 0.5),
              axis.title.y = element_text(angle = 90, margin = unit(c(0, 0, 0, 0), 'pt'),
                                          size = 9),
              plot.margin = unit(c(-70, 0, -70,0), 'pt')) + ggtitle('')+ ylab('Recruitment')+
        scale_color_gradient(low = "white", high = "#4C0000")

plot_dp3 <- draw_hexagons(det[[1]], z = colMeans(dp3, na.rm = T), size = 1, 
              add = draw_hexagons(det[[1]], size = 1.05,
                                  add = draw_FoodPatches(det[[1]], fill = 'grey70')))+
        geom_point(data = hex[634, ], aes(x, y-50), size = 2, color = 'grey40', shape = 17)+
        theme_void()+
        theme(aspect.ratio = 0.5, legend.position = 'none', plot.title = element_text(hjust = 0.5),
              axis.title.y = element_text(angle = 90, margin = unit(c(0, 0, 0, 0), 'pt'),
                                          size = 9),
              plot.margin = unit(c(-140, 0, 0,0), 'pt')) + ggtitle('')+ ylab('Post-recruitment')+
        scale_color_gradient(low = "white", high = "#4C0000")


#######                      #######
## +++ STOCHASTIC EXPERIMENTS +++ ##
#######                      #######

d <- exp_spreadsheet[exp_spreadsheet$CONDITION == 'STOCH', c('Date', 'MATI.TARDA')][1:10, ]
d <- paste(format(strptime(d[, 1], '%d/%m/%Y'), format = '%Y%m%d'), d[, 2], sep='')

for(i in seq_along(sto)){
        sto[[i]]$date <- d[i]
}

sto_phase <- c()

for(i in seq_along(sto)){
        food <- get_foodPatches(sto[[i]])
        sto[[i]]$food <- food[1:2]
        x <- food_detection(sto[[i]])
        sto[[i]]$food$GP1$t <- x[1:6]
        sto[[i]]$food$GP2$t <- x[7:12]
        sto_phase <- rbind(sto_phase, c(min(x, na.rm = T), max(x, na.rm = T)))
}

for(i in seq_along(sto)){
        print(paste0('Iter ',i))
        sto[[i]]$mutual_info_p1 <- mutual_info(sto[[i]], 1:sto_phase[i, 1])
        sto[[i]]$mutual_info_p2 <- mutual_info(sto[[i]], (sto_phase[i, 1] + 1):sto_phase[i, 2])
        sto[[i]]$mutual_info_p3 <- mutual_info(sto[[i]], (sto_phase[i, 2]+1):max(sto[[i]]$data$Frame))
}

sp1 <- c()
sp2 <- c()
sp3 <- c()

for(i in seq_along(sto)){
        sp1 <- rbind(sp1, sto[[i]]$mutual_info_p1)
        sp2 <- rbind(sp2, sto[[i]]$mutual_info_p2)
        sp3 <- rbind(sp3, sto[[i]]$mutual_info_p3)
}

plot_sp1 <- draw_hexagons(det[[1]], z = colMeans(sp1[-2, ], na.rm = T), size = 1, 
                          add = draw_hexagons(det[[1]], size = 1.05,
                                              add = draw_FoodPatches(sto[-2], fill = 'grey70')))+
        geom_point(data = hex[634, ], aes(x, y-50), size = 2, color = 'grey40', shape = 17)+
        theme_void()+
        theme(aspect.ratio = 0.5, legend.position = 'none', plot.title = element_text(hjust = 0.5),
              axis.title.y = element_text(angle = 90, margin = unit(c(0, 0, 0, 0), 'pt'),
                                          size = 9),
              plot.margin = unit(c(0, 0, -140,0), 'pt')) + ggtitle('STO')+ ylab('')+
        scale_color_gradient(low = "white", high = "#4C0000")

plot_sp2 <- draw_hexagons(det[[1]], z = colMeans(sp2[-2, ], na.rm = T), size = 1, 
                          add = draw_hexagons(det[[1]], size = 1.05,
                                              add = draw_FoodPatches(sto[-2], fill = 'grey70')))+
        geom_point(data = hex[634, ], aes(x, y-50), size = 2, color = 'grey40', shape = 17)+
        theme_void()+
        theme(aspect.ratio = 0.5, legend.position = 'none', plot.title = element_text(hjust = 0.5),
              axis.title.y = element_text(angle = 90, margin = unit(c(0, 0, 0, 0), 'pt'),
                                          size = 9),
              plot.margin = unit(c(-70, 0, -70,0), 'pt')) + ggtitle('')+ ylab('')+
        scale_color_gradient(low = "white", high = "#4C0000")

plot_sp3 <- draw_hexagons(det[[1]], z = colMeans(sp3[-2, ], na.rm = T), size = 1, 
                          add = draw_hexagons(det[[1]], size = 1.05,
                                              add = draw_FoodPatches(sto[-2], fill = 'grey70')))+
        geom_point(data = hex[634, ], aes(x, y-50), size = 2, color = 'grey40', shape = 17)+
        theme_void()+
        theme(aspect.ratio = 0.5, legend.position = 'none', plot.title = element_text(hjust = 0.5),
              axis.title.y = element_text(angle = 90, margin = unit(c(0, 0, 0, 0), 'pt'),
                                          size = 9),
              plot.margin = unit(c(-140, 0, 0,0), 'pt')) + ggtitle('')+ ylab('')+
        scale_color_gradient(low = "white", high = "#4C0000")


png(filename = '~/research/2022/ANTS/Figs_spinGlasses/fig_sencera.png', 4000, 6000, res = 900)
grid.arrange(plot_dp1, plot_sp1, plot_dp2,
             plot_sp2, plot_dp3, plot_sp3, ncol = 2, nrow = 3)
dev.off()
