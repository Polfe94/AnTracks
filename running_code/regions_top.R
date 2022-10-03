
# regions <- list(tl = det[[1]]$segments[det[[1]]$segments$x > 0 & det[[1]]$segments$x < 2000/3+20 &
#                                             det[[1]]$segments$y > 1500 & det[[1]]$segments$y < 2000, ],
#                 tm = det[[1]]$segments[det[[1]]$segments$x > 2000/3+20 & det[[1]]$segments$x < 2*2000/3 &
#                                             det[[1]]$segments$y > 1500 & det[[1]]$segments$y < 2000, ],
#                 tr = det[[1]]$segments[det[[1]]$segments$x > 2*2000/3 & det[[1]]$segments$x < 2000 &
#                                             det[[1]]$segments$y > 1500 & det[[1]]$segments$y < 2000, ],
#                 br = det[[1]]$segments[det[[1]]$segments$x > 2*2000/3 & det[[1]]$segments$x < 2000 &
#                                             det[[1]]$segments$y > 1000 & det[[1]]$segments$y < 1500, ],
#                 bm = det[[1]]$segments[det[[1]]$segments$x > 2000/3+20 & det[[1]]$segments$x < 2*2000/3 &
#                                             det[[1]]$segments$y > 1000 & det[[1]]$segments$y < 1500, ],
#                 bl = det[[1]]$segments[det[[1]]$segments$x > 0 & det[[1]]$segments$x < 2000/3+20 &
#                                             det[[1]]$segments$y > 1000 & det[[1]]$segments$y < 1500, ])


regions <- rbind(cbind.data.frame(det[[1]]$segments[det[[1]]$segments$x > 0 & det[[1]]$segments$x < 2000/3+20 &
                                               det[[1]]$segments$y > 1500 & det[[1]]$segments$y < 2000, ],
                            region = 'tl'),
                cbind.data.frame(det[[1]]$segments[det[[1]]$segments$x > 2000/3+20 & det[[1]]$segments$x < 2*2000/3 &
                                               det[[1]]$segments$y > 1500 & det[[1]]$segments$y < 2000, ],
                                 region = 'tm'),
                cbind.data.frame(det[[1]]$segments[det[[1]]$segments$x > 2*2000/3 & det[[1]]$segments$x < 2000 &
                                               det[[1]]$segments$y > 1500 & det[[1]]$segments$y < 2000, ],
                                 region = 'tr'),
                cbind.data.frame(det[[1]]$segments[det[[1]]$segments$x > 2*2000/3 & det[[1]]$segments$x < 2000 &
                                               det[[1]]$segments$y > 1000 & det[[1]]$segments$y < 1500, ],
                                 region = 'br'),
                cbind.data.frame(det[[1]]$segments[det[[1]]$segments$x > 2000/3+20 & det[[1]]$segments$x < 2*2000/3 &
                                               det[[1]]$segments$y > 1000 & det[[1]]$segments$y < 1500, ],
                                 region = 'bm'),
                cbind.data.frame(det[[1]]$segments[det[[1]]$segments$x > 0 & det[[1]]$segments$x < 2000/3+20 &
                                               det[[1]]$segments$y > 1000 & det[[1]]$segments$y < 1500, ],
                                 region = 'bl'))


if(exists('do_plot') && do_plot == TRUE){
     print(
        suppressMessages(
        draw_hexagons(det[[1]]) + geom_polygon(data = data.frame(x = c(0, 0, 2000, 2000), 
                                                                 y = c(1000, 2000, 2000, 1000)), aes(x,y), 
                                               fill = NA, color = 'red')+
                geom_line(data = data.frame(x = c(0, 2000+20), y = c(1500, 1500)), aes(x, y), color = 'red')+
                geom_line(data = data.frame(x = c(20+2000/3, 0+2000/3), y = c(2000, 1000)), aes(x, y), color = 'red')+
                geom_line(data = data.frame(x =c(2000/3+2000/3, 2000/3+2000/3), y = c(2000, 1000)), aes(x, y), color = 'red')+
                geom_point(data = regions, aes(x, y, color = factor(region)), show.legend = F)+
                scale_color_viridis_d()
        ))
}

# draw_hexagons(det[[1]]) + geom_polygon(data = data.frame(x = c(0, 0, 2000, 2000), 
#                                                          y = c(1000, 2000, 2000, 1000)), aes(x,y), 
#                                        fill = NA, color = 'red')+
#      geom_line(data = data.frame(x = c(0, 2000), y = c(1500, 1500)), aes(x, y), color = 'red')+
#      geom_line(data = data.frame(x = c(0+2000/3, 0+2000/3), y = c(2000, 1000)), aes(x, y), color = 'red')+
#      geom_line(data = data.frame(x =c(2000/3+2000/3, 2000/3+2000/3), y = c(2000, 1000)), aes(x, y), color = 'red')