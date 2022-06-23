exp_condition <- c('')

source('G:/research/MoveAnts/code/config.R')
source('C:/Users/POL/Desktop/gants_current/code/dev/dev_functions.R')

hex <- hex90
colnames(hex) <- c('x', 'y')
h <- hex90[hex90$rotY > 980, ]
colnames(h) <- c('x', 'y')
source('G:/research/2022/AnTracks/src/generic.R')
source('G:/research/2022/AnTracks/src/coords.R')

rm(list = ls()[grepl('json', ls())], visual, colony, current_dir, dirL, exp_condition, exps, conditionL)

load('G:/research/2022/AnTracks/results/sto_coords.RData')
load('G:/research/2022/AnTracks/results/det_coords.RData')


############## ++++ DETERMINIST EXPERIMENTS ++++ ############## 

#### PHASE 1 ####
dp1 <- init(read.table('G:/research/2022/Figs_SpinGlasses/Det_I.dat', 
                       col.names = c('Xmm', 'Ymm', 'Zsim', 'Zexp')), h, 'coords', segments = det[[1]]$segments)
dp1$data$Ymm <- dp1$data$Ymm + 995
dp1 <- closest_node(dp1)
df_dp1 <-  data.frame(x = h[dp1$data$node, 1], y = h[dp1$data$node, 2], z = dp1$data$Zexp, z1 = dp1$data$Zsim)


d1a <- draw_FoodPatches(det[[1]], fill = 'grey25', add = draw_hexagons(dp1, color = 'grey85',
                                                                add = ggplot()+ geom_point(data = df_dp1, aes(x,y,color = z),
                                                                                           size = 5, show.legend = F)+
                                                                        scale_color_gradient(low = 'white', high = "#4C0000")+
                                                                        theme(aspect.ratio = 0.5)) )+
        theme(plot.margin = unit(c(-20,10,-260,-20), 'pt'))
d1b <- draw_FoodPatches(det[[1]], fill = 'grey25', add = draw_hexagons(dp1, color = 'grey85',
                                                                add = ggplot()+ geom_point(data = df_dp1, aes(x,y,color = z1),
                                                                                           size = 5, show.legend = F)+
                                                                        scale_color_gradient(low = 'white', high = "#4C0000")+
                                                                        theme(aspect.ratio = 0.5)) )+
        theme(plot.margin = unit(c(-260,10,0,-20), 'pt'))


#### PHASE 2 ####
dp2 <- init(read.table('G:/research/2022/Figs_SpinGlasses/Det_II.dat', 
                       col.names = c('Xmm', 'Ymm', 'Zsim', 'Zexp')), h, 'coords', segments = det[[1]]$segments)
dp2$data$Ymm <- dp2$data$Ymm + 995
dp2 <- closest_node(dp2)
df_dp2 <-  data.frame(x = h[dp2$data$node, 1], y = h[dp2$data$node, 2], z = dp2$data$Zexp, z1 = dp2$data$Zsim)

d2a <- draw_FoodPatches(det[[1]], fill = 'grey25', add = draw_hexagons(dp2, color = 'grey85', add = ggplot()+ geom_point(data = df_dp2, aes(x,y,color = z),
                                                                                                                  size = 5, show.legend = F)+
                                                                        scale_color_gradient(low = 'white', high = "#4C0000")+
                                                                        theme(aspect.ratio = 0.5)) )+
        theme(plot.margin = unit(c(-20,10,-260,-20), 'pt'))
d2b <- draw_FoodPatches(det[[1]], fill = 'grey25', add = draw_hexagons(dp2, color = 'grey85', add = ggplot()+ geom_point(data = df_dp2, aes(x,y,color = z1),
                                                                                                                  size = 5, show.legend = F)+
                                                                        scale_color_gradient(low = 'white', high = "#4C0000")+
                                                                        theme(aspect.ratio = 0.5)) )+
        theme(plot.margin = unit(c(-260,10,0,-20), 'pt'))


#### PHASE 3 ####
dp3 <- init(read.table('G:/research/2022/Figs_SpinGlasses/Det_III.dat', 
                       col.names = c('Xmm', 'Ymm', 'Zsim', 'Zexp')), h, 'coords', segments = det[[1]]$segments)
dp3$data$Ymm <- dp3$data$Ymm + 995
dp3 <- closest_node(dp3)
df_dp3 <-  data.frame(x = h[dp3$data$node, 1], y = h[dp3$data$node, 2], z = dp3$data$Zexp, z1 = dp3$data$Zsim)

d3a <- draw_FoodPatches(det[[1]], fill = 'grey25', add = draw_hexagons(dp3, color = 'grey85', add = ggplot()+ 
                                                                        geom_point(data = df_dp3, aes(x,y,color = z),
                                                                                   size = 5, show.legend = F)+
                                                                        scale_color_gradient(low = 'white', high = "#4C0000")+
                                                                        theme(aspect.ratio = 0.5)) )+
        theme(plot.margin = unit(c(-20,10,-260,-20), 'pt'))
d3b <- draw_FoodPatches(det[[1]], fill = 'grey25', add = draw_hexagons(dp3, color = 'grey85', add = ggplot()+ 
                                                                        geom_point(data = df_dp3, aes(x,y,color = z1),
                                                                                   size = 5, show.legend = F)+
                                                                        scale_color_gradient(low = 'white', high = "#4C0000")+
                                                                        theme(aspect.ratio = 0.5)) )+
        theme(plot.margin = unit(c(-260,10,0,-20), 'pt'))

png(filename = 'G:/research/2022/Figs_SpinGlasses/test_global_det.png', 6000, 4000, res = 300)
ggarrange(d1a, d2a, d3a, d1b, d2b, d3b, nrow = 2, ncol = 3, legend = 'none')
dev.off()


############## ++++ STOCHASTIC EXPERIMENTS ++++ ############## 

#### PHASE 1 ####
sp1 <- init(read.table('G:/research/2022/Figs_SpinGlasses/Sto_I.dat', 
                       col.names = c('Xmm', 'Ymm', 'Zsim', 'Zexp')), h, 'coords', segments = det[[1]]$segments)
sp1$data$Ymm <- sp1$data$Ymm + 995
sp1 <- closest_node(sp1)
df_sp1 <-  data.frame(x = h[sp1$data$node, 1], y = h[sp1$data$node, 2], z = sp1$data$Zexp, z1 = sp1$data$Zsim)


s1a <- draw_FoodPatches(sto, fill = 'grey25', add = draw_hexagons(sp1, color = 'grey85',
                                                                add = ggplot()+ geom_point(data = df_sp1, aes(x,y,color = z),
                                                                                           size = 5, show.legend = F)+
                                                                        scale_color_gradient(low = 'white', high = "#4C0000")+
                                                                        theme(aspect.ratio = 0.5)) )+
        theme(plot.margin = unit(c(-20,10,-260,-20), 'pt'))
s1b <- draw_FoodPatches(sto, fill = 'grey25', add = draw_hexagons(sp1, color = 'grey85',
                                                                add = ggplot()+ geom_point(data = df_sp1, aes(x,y,color = z1),
                                                                                           size = 5, show.legend = F)+
                                                                        scale_color_gradient(low = 'white', high = "#4C0000")+
                                                                        theme(aspect.ratio = 0.5)) )+
        theme(plot.margin = unit(c(-260,10,0,-20), 'pt'))


#### PHASE 2 ####
sp2 <- init(read.table('G:/research/2022/Figs_SpinGlasses/Sto_II.dat', 
                       col.names = c('Xmm', 'Ymm', 'Zsim', 'Zexp')), h, 'coords', segments = det[[1]]$segments)
sp2$data$Ymm <- sp2$data$Ymm + 995
sp2 <- closest_node(sp2)
df_sp2 <-  data.frame(x = h[sp2$data$node, 1], y = h[sp2$data$node, 2], z = sp2$data$Zexp, z1 = sp2$data$Zsim)

s2a <- draw_FoodPatches(sto, fill = 'grey25', add = draw_hexagons(sp2, color = 'grey85', add = ggplot()+ geom_point(data = df_sp2, aes(x,y,color = z),
                                                                                                                  size = 5, show.legend = F)+
                                                                        scale_color_gradient(low = 'white', high = "#4C0000")+
                                                                        theme(aspect.ratio = 0.5)) )+
        theme(plot.margin = unit(c(-20,10,-260,-20), 'pt'))
s2b <- draw_FoodPatches(sto, fill = 'grey25', add = draw_hexagons(sp2, color = 'grey85', add = ggplot()+ geom_point(data = df_sp2, aes(x,y,color = z1),
                                                                                                                  size = 5, show.legend = F)+
                                                                        scale_color_gradient(low = 'white', high = "#4C0000")+
                                                                        theme(aspect.ratio = 0.5)) )+
        theme(plot.margin = unit(c(-260,10,0,-20), 'pt'))


#### PHASE 3 ####
sp3 <- init(read.table('G:/research/2022/Figs_SpinGlasses/Sto_III.dat', 
                       col.names = c('Xmm', 'Ymm', 'Zsim', 'Zexp')), h, 'coords', segments = det[[1]]$segments)
sp3$data$Ymm <- sp3$data$Ymm + 995
sp3 <- closest_node(sp3)
df_sp3 <-  data.frame(x = h[sp3$data$node, 1], y = h[sp3$data$node, 2], z = sp3$data$Zexp, z1 = sp3$data$Zsim)

s3a <- draw_FoodPatches(sto, fill = 'grey25', add = draw_hexagons(sp3, color = 'grey85', add = ggplot()+ 
                                                                        geom_point(data = df_sp3, aes(x,y,color = z),
                                                                                   size = 5, show.legend = F)+
                                                                        scale_color_gradient(low = 'white', high = "#4C0000")+
                                                                        theme(aspect.ratio = 0.5)) )+
        theme(plot.margin = unit(c(-20,10,-260,-20), 'pt'))
s3b <- draw_FoodPatches(sto, fill = 'grey25', add = draw_hexagons(sp3, color = 'grey85', add = ggplot()+ 
                                                                        geom_point(data = df_sp3, aes(x,y,color = z1),
                                                                                   size = 5, show.legend = F)+
                                                                        scale_color_gradient(low = 'white', high = "#4C0000")+
                                                                        theme(aspect.ratio = 0.5)) )+
        theme(plot.margin = unit(c(-260,10,0,-20), 'pt'))

png(filename = 'G:/research/2022/Figs_SpinGlasses/test_global_sto.png', 6000, 4000, res = 300)
ggarrange(s1a, s2a, s3a, s1b, s2b, s3b, nrow = 2, ncol = 3, legend = 'none')
dev.off()
