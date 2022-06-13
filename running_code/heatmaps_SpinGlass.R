exp_condition <- c('')

source('G:/research/MoveAnts/code/config.R')
source('C:/Users/POL/Desktop/gants_current/code/dev/dev_functions.R')

hex <- hex90
colnames(hex) <- c('x', 'y')
h <- hex90[hex90$rotY > 980, ]
colnames(h) <- c('x', 'y')
source('G:/research/2022/AnTracks/src/generic.R')
source('G:/research/2022/AnTracks/src/coords.R')

load('G:/research/2022/AnTracks/results/sto_coords.RData')
load('G:/research/2022/AnTracks/results/det_coords.RData')

rm(list = ls()[grepl('json', ls())], visual, colony, current_dir, dirL, exp_condition, exps, conditionL)

library(scico)

## read data
# det
dp1 <- init(read.table('G:/research/2022/Figs_SpinGlasses/Det_I.dat', 
                       col.names = c('Xmm', 'Ymm', 'Zexp', 'Zsim')), h, 'coords', segments = det[[1]]$segments)
dp1$data$Ymm <- dp1$data$Ymm + 995
dp1 <- closest_node(dp1)
df_dp1 <-  data.frame(x = h[dp1$data$node, 1], y = h[dp1$data$node, 2], z = dp1$data$Zexp)
draw_hexagons(dp1, color = 'grey85', add = ggplot()+ geom_point(data = df_dp1, aes(x,y,color = z),
                                                                          size = 5, alpha = 0.9)+
                      scale_color_gradient(low = 'white', high = muted('red'))+
                      theme(aspect.ratio = 0.5)) 
# df_dp1$z[df_dp1$z < 0] <- norm_range(df_dp1$z[df_dp1$z < 0], -1, -.Machine$double.eps)
# df_dp1$z[df_dp1$z > 0] <- norm_range(df_dp1$z[df_dp1$z > 0], .Machine$double.eps, 1)
# draw_hexagons(dp1, size = 2) + geom_point(data = df_dp1, aes(x,y,color = z), size = 5) + 
#         scale_color_gradient2(limits = c(-1, 1), low = muted('blue'), high = muted('red'))


dp2 <- init(read.table('G:/research/2022/Figs_SpinGlasses/Det_II.dat', 
                       col.names = c('Xmm', 'Ymm', 'Zexp', 'Zsim')), h, 'coords', segments = det[[1]]$segments)
dp2$data$Ymm <- dp2$data$Ymm + 995
dp2 <- closest_node(dp2)
df_dp2 <-  data.frame(x = h[dp2$data$node, 1], y = h[dp2$data$node, 2], z = dp2$data$Zexp)
draw_hexagons(dp2, color = 'grey85', add = ggplot()+ geom_point(data = df_dp2, aes(x,y,color = z),
                                                                size = 5, alpha = 0.9)+
                      scale_color_gradient(low = 'white', high = muted('red'))+
                      theme(aspect.ratio = 0.5))
# df_dp2$z <- norm_range(df_dp2$z, -1, 0)
# df_dp2$z[df_dp2$z < 0] <- norm_range(df_dp2$z[df_dp2$z < 0], -1, -.Machine$double.eps)
# df_dp2$z[df_dp2$z > 0] <- norm_range(df_dp2$z[df_dp2$z > 0], .Machine$double.eps, 1)
# draw_hexagons(dp2, size = 2, color = 'grey80') + geom_point(data = df_dp2, aes(x,y,color = z),
#                                                             size = 7.5, alpha = 0.7)+
#         scale_color_viridis(option = 'plasma')
# draw_hexagons(dp2, size = 2) + geom_point(data = df_dp2, aes(x,y,color = z), size = 5) +
#         scale_color_gradient2(limits = c(-1, 0), low = muted('blue'), high = muted('red'), midpoint = -0.5)


dp3 <- init(read.table('G:/research/2022/Figs_SpinGlasses/Det_III.dat', 
                       col.names = c('Xmm', 'Ymm', 'Zexp', 'Zsim')), h, 'coords', segments = det[[1]]$segments)
dp3$data$Ymm <- dp3$data$Ymm + 995
dp3 <- closest_node(dp3)
df_dp3 <-  data.frame(x = h[dp3$data$node, 1], y = h[dp3$data$node, 2], z = dp3$data$Zexp)
draw_hexagons(dp3, color = 'grey85', add = ggplot()+ geom_point(data = df_dp3, aes(x,y,color = z),
                                                                size = 5, alpha = 0.9)+
                      scale_color_gradient(low = 'white', high = "#4C0000")+
                      theme(aspect.ratio = 0.5))
# draw_hexagons(dp3, size = 2, color = 'grey80') + geom_point(data = df_dp3, aes(x,y,color = z),
#                                                             size = 7.5, alpha = 0.7)+
#         scale_color_viridis(option = 'plasma')
# df_dp3$z[df_dp3$z < 0] <- norm_range(df_dp3$z[df_dp3$z < 0], -1, -.Machine$double.eps)
# df_dp3$z[df_dp3$z > 0] <- norm_range(df_dp3$z[df_dp3$z > 0], .Machine$double.eps, 1)
# draw_hexagons(dp3, size = 2) + geom_point(data = df_dp3, aes(x,y,color = z), size = 5) + 
#         scale_color_gradient2(limits = c(-1, 1), low = muted('blue'), high = muted('red'))


draw_hexagons(dp1) + geom_point(data = data.frame(x = h[dp1$data$node, 1],
                                                  y = h[dp1$data$node, 2],
                                                  z = dp1$data$Zexp), aes(x, y, color = z))

# sto
sp1 <- read.table('G:/research/2022/Figs_SpinGlasses/Sto_I.dat', 
                  col.names = c('Xmm', 'Ymm', 'Zexp', 'Zsim'))
sp2 <- read.table('G:/research/2022/Figs_SpinGlasses/Sto_II.dat', 
                  col.names = c('Xmm', 'Ymm', 'Zexp', 'Zsim'))
sp3 <- read.table('G:/research/2022/Figs_SpinGlasses/Sto_III.dat', 
                  col.names = c('Xmm', 'Ymm', 'Zexp', 'Zsim'))


draw_hexagons(det[[1]])+ geom_point(data = dp1, aes(V1-2.5, V2+995, color = V3))

## DET ##


ggplot(data = dp1, aes(V1, V2, color = V3)) + geom_point()


ggplot(data = dp2, aes(V1, V2, color = V3)) + geom_point()


ggplot(data = dp3, aes(V1, V2, color = V3)) + geom_point()

## STO ##


ggplot(data = sp1, aes(V1, V2, color = V3)) + geom_point()


ggplot(data = sp2, aes(V1, V2, color = V3)) + geom_point()


ggplot(data = sp3, aes(V1, V2, color = V3)) + geom_point()

