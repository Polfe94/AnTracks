f <- function(x, t){
     i <- 1 + t/2
     xf <- length(x) - t/2
     y <- x
     
     while(i <= xf){
          y[i] <- mean(x[(i-t/2):(i+t/2)], na.rm = T)
          i <- i+1
     }
     y
}

#### +++ DETERMINIST EXPERIMENTS +++ ####

df_det <- c()
for(i in seq_along(det)){
     
     df_det <- rbind(df_det, data.frame(t = seq_along(det[[i]]$I),
                                      i = f(det[[i]]$I, 50),
                                      n = f(det[[i]]$N, 50),
                                      exp = i))
}
det_times <- cbind.data.frame(x = c(rep(det_phase[, 1],2), rep(det_phase[, 2], 2)),
                              y = c(rep(0, 10), rep(Inf, 20), rep(0, 10)),
                              exp = rep(1:10, 2))

# overlapping interactions (all exps in the same plot)
ggplot(data = df_det, aes(t / 120, i, color = factor(exp))) + geom_line(show.legend = F)+
     scale_color_viridis_d() + scale_x_continuous('Time (min)', breaks = seq(0, 180, 15), limits = c(0, 180))+
     geom_polygon(data = data.frame(x = c(rep(mean(det_phase[, 1]), 2), rep(mean(det_phase[, 2]),2)),
                                    y = c(0, Inf, Inf, 0)), fill = viridis(10)[7], 
                  color = viridis(10)[7], alpha = 0.2, aes(x/120, y), linetype = 2)+
     scale_y_continuous('Interactions', breaks = seq(0, 3, length.out = 7))

# faceted interactions and N by experiment
ggplot(data = df_det, aes(t / 120, norm_range(i, a = 0, b = 1))) + geom_line(show.legend = F, color = viridis(5)[3])+
     geom_line(data = df_det, aes(t / 120, norm_range(n, a = 0, b = 1)), color = 'black', show.legend = F)+
     geom_polygon(data = det_times, aes(x/120, y), fill = viridis(10)[7], 
                  color = viridis(10)[7], alpha = 0.2, show.legend = F, linetype = 2)+
     facet_wrap(~factor(exp))+ theme(strip.text = element_text(size = 15))+
      scale_y_continuous('Interactions', breaks = seq(0, 3/max(df_det$i), length.out = 5), 
                         labels = seq(0, 3, length.out = 5), 
                         sec.axis = sec_axis(trans = function(x) x, breaks = seq(0, 45/max(df_det$n), length.out = 10),
                                             labels = seq(0, 45, length.out = 10))) +
     scale_x_continuous('Time (min)', breaks = seq(0, 180, 25), limits = c(0, 180))


#### +++ STOCHASTIC EXPERIMENTS +++ ####
df_sto <- c()
for(i in seq_along(sto)){
     
     df_sto <- rbind(df_sto, data.frame(t = seq_along(sto[[i]]$I),
                                        i = f(sto[[i]]$I, 50),
                                        n = f(sto[[i]]$N, 50),
                                        exp = i))
}
sto_times <- cbind.data.frame(x = c(rep(sto_phase[, 1],2), rep(sto_phase[, 2], 2)),
                              y = c(rep(0, 10), rep(Inf, 20), rep(0, 10)),
                              exp = rep(1:10, 2))

# overlapping interactions (all exps in the same plot)
ggplot(data = df_sto, aes(t / 120, i, color = factor(exp))) + geom_line(show.legend = F)+
     scale_color_viridis_d() + scale_x_continuous('Time (min)', breaks = seq(0, 180, 15), limits = c(0, 180))+
     geom_polygon(data = data.frame(x = c(rep(mean(sto_phase[, 1]), 2), rep(mean(sto_phase[, 2]),2)),
                                    y = c(0, Inf, Inf, 0)), fill = viridis(10)[7], 
                  color = viridis(10)[7], alpha = 0.2, aes(x/120, y), linetype = 2)+
     scale_y_continuous('Interactions', breaks = seq(0, 2.5, length.out = 6))

# faceted interactions and N by experiment
ggplot(data = df_sto, aes(t / 120, norm_range(i, a = 0, b = 1))) + geom_line(show.legend = F, color = viridis(5)[3])+
     geom_line(data = df_sto, aes(t / 120, norm_range(n, a = 0, b = 1)), color = 'black', show.legend = F)+
     geom_polygon(data = sto_times, aes(x/120, y), fill = viridis(10)[7], 
                  color = viridis(10)[7], alpha = 0.2, show.legend = F, linetype = 2)+
     facet_wrap(~factor(exp))+ theme(strip.text = element_text(size = 15))+
     scale_y_continuous('Interactions', breaks = seq(0, 2.5/max(df_sto$i), length.out = 6), 
                        labels = seq(0, 2.5, length.out = 6), 
                        sec.axis = sec_axis(trans = function(x) x, breaks = seq(0, 55/max(df_sto$n), length.out = 12),
                                            labels = seq(0, 55, length.out = 12))) +
     scale_x_continuous('Time (min)', breaks = seq(0, 180, 25), limits = c(0, 180))
