source('/home/polfer/research/2022/ANTS/AnTracks/src/Experiment.R')
path <- '/home/polfer/Documents/Tracker_LEOV/ANTS/data/'
exps <- list.files(path)
load('~/research/2022/ANTS/AnTracks/data/det.RData')
load('~/research/2022/ANTS/AnTracks/data/sto.RData')
load('~/research/2022/ANTS/AnTracks/data/nf.RData')


#### +++ DETERMINIST EXPS +++ ####



idx <- exp_spreadsheet$CONDITION == 'DET'
det_exps <- vector('list', sum(idx))
names(det_exps) <- paste0(format(as.POSIXct(exp_spreadsheet$Date, format = '%d/%m/%Y'), '%Y%m%d'),
                          exp_spreadsheet$MATI.TARDA)[idx]

for(exp in names(det_exps)){
        jsonexp <- read_json(paste0(path, exp, '/iTracks.json'))
        det_exps[[exp]] <- alpha(jsonexp)
}

det_ccfs <- lapply(det_exps, function(i){
        ccf(moving_average(i$x, 60, 30), moving_average(i$y, 60, 30), lag.max = 1000)
})

avg_det_ccf <- colMeans(do.call('rbind', lapply(det_ccfs, function(i){
        i$acf
})))
plot(-1000:1000, avg_det_ccf)
plot(0:1000, avg_det_ccf[1000:2000], type = 'l')

t0 <- Sys.time()
DET_results <- lapply(1:1200, function(i){
     det_exps <- vector('list', sum(idx))
     names(det_exps) <- paste0(format(as.POSIXct(exp_spreadsheet$Date, format = '%d/%m/%Y'), '%Y%m%d'),
                               exp_spreadsheet$MATI.TARDA)[idx]
     
     for(exp in names(det_exps)){
          jsonexp <- read_json(paste0(path, exp, '/iTracks.json'))
          det_exps[[exp]] <- alpha(jsonexp, i)
     }
     
     r <- vapply(seq_along(det_exps), function(z){
             x <- moving_average(det_exps[[z]]$x, 60, 30)
             y <- moving_average(det_exps[[z]]$y, 60, 30)
             
             x[x >= 1] <- 0
             y[y >= 1] <- 0
             cor(x, y)
     }, numeric(1))
     c(mean = mean(r), sd = sd(r))
})
Sys.time() - t0



# DET_CORS <- lapply(DET_results, function(i){
# 
#      
#      r <- vapply(seq_along(i), function(z){
#           x <- moving_average(i[[z]]$x, 60, 30)
#           y <- moving_average(i[[z]]$y, 60, 30)
#           
#           x[x >= 1] <- 0
#           y[y >= 1] <- 0
#           cor(x, y)
#      }, numeric(1))
#      
#      c(mean = mean(r), sd = sd(r))
# 
# })





# df <- do.call('rbind.data.frame', DET_CORS)
df <- do.call('rbind.data.frame', DET_results)
colnames(df) <- c('mean', 'sd')
df$t <- 1:nrow(df)
df <- rbind(df, c(mean = 0.85039122, sd = 0.05570747, t = 0))
plot(0:1000, avg_det_ccf[1000:2000], type = 'l', ylim = c(0, 1))
lines(df$t[1:600], df$mean[1:600], col = 'red')

ggplot(data = df, aes(t, mean)) + 
     geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.3, color = 'black')+ 
     geom_point(size = 3)+ geom_line()+
     scale_x_continuous('Temporal buffer (frames)', breaks = seq(0, 600, length.out = 8))+
     scale_y_continuous('Average correlation', breaks = seq(0, 1, length.out = 6))+
     ggtitle('Determinist') + theme(plot.title = element_text(hjust = 0.5, size = 16))
# fitdistr(df$mean, 'exponential')
# fitted_exp <- rexp(1000, rate = 2.2567071)
# fitted_exp <- fitted_exp[fitted_exp <= max(df$mean) & fitted_exp >= min(df$mean)]
# fitted_exp <- sort(fitted_exp, decreasing = T)
# tx <- seq_along(fitted_exp) / length(fitted_exp) * 60
# plot(df$t, df$mean)
# lines(tx, fitted_exp, col = 'red')
# 
# ggplot(data= df, aes(t, log(mean))) + geom_point()

#### +++ STOCHASTIC EXPS +++ ####

idx <- exp_spreadsheet$CONDITION == 'STOCH'
sto_exps <- vector('list', sum(idx))
names(sto_exps) <- paste0(format(as.POSIXct(exp_spreadsheet$Date, format = '%d/%m/%Y'), '%Y%m%d'),
                          exp_spreadsheet$MATI.TARDA)[idx]


# t0 <- Sys.time()
# STO_results <- lapply(1:600, function(i){
#      sto_exps <- vector('list', sum(idx))
#      names(sto_exps) <- paste0(format(as.POSIXct(exp_spreadsheet$Date, format = '%d/%m/%Y'), '%Y%m%d'),
#                                exp_spreadsheet$MATI.TARDA)[idx]
#      
#      for(exp in names(sto_exps)){
#           jsonexp <- read_json(paste0(path, exp, '/iTracks.json'))
#           sto_exps[[exp]] <- alpha(jsonexp, i)
#      }
#      sto_exps
# })
# Sys.time() - t0

t0 <- Sys.time()
STO_results <- lapply(1:1200, function(i){
        sto_exps <- vector('list', sum(idx))
        names(sto_exps) <- paste0(format(as.POSIXct(exp_spreadsheet$Date, format = '%d/%m/%Y'), '%Y%m%d'),
                                  exp_spreadsheet$MATI.TARDA)[idx]
        
        for(exp in names(sto_exps)){
                jsonexp <- read_json(paste0(path, exp, '/iTracks.json'))
                sto_exps[[exp]] <- alpha(jsonexp, i)
        }
        
        r <- vapply(seq_along(sto_exps), function(z){
                x <- moving_average(sto_exps[[z]]$x, 60, 30)
                y <- moving_average(sto_exps[[z]]$y, 60, 30)
                
                x[x >= 1] <- 0
                y[y >= 1] <- 0
                cor(x, y)
        }, numeric(1))
        c(mean = mean(r), sd = sd(r))
})
Sys.time() - t0

# STO_CORS <- lapply(STO_results, function(i){
#      
#      
#      r <- vapply(seq_along(i), function(z){
#           x <- moving_average(i[[z]]$x, 60, 30)
#           y <- moving_average(i[[z]]$y, 60, 30)
#           
#           x[x >= 1] <- 0
#           y[y >= 1] <- 0
#           cor(x, y)
#      }, numeric(1))
#      
#      c(mean = mean(r), sd = sd(r))
#      
# })

df <- do.call('rbind.data.frame', STO_results)
colnames(df) <- c('mean', 'sd')
df$t <- 1:nrow(df)
ggplot(data = df, aes(t, mean)) + 
     geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.3, color = 'black')+ 
     geom_point(size = 3)+ geom_line()+
     scale_x_continuous('Temporal buffer (frames)', breaks = seq(0, 60, 10))+
     scale_y_continuous('Average correlation', breaks = seq(0, 1, length.out = 6))+
     ggtitle('Stochastic') + theme(plot.title = element_text(hjust = 0.5, size = 16))
#### +++ NO FOOD EXPS +++ ####

idx <- exp_spreadsheet$CONDITION == 'NO FOOD'
nf_exps <- vector('list', sum(idx))
names(nf_exps) <- paste0(format(as.POSIXct(exp_spreadsheet$Date, format = '%d/%m/%Y'), '%Y%m%d'),
                          exp_spreadsheet$MATI.TARDA)[idx]


# t0 <- Sys.time()
# NF_results <- lapply(1:600, function(i){
#      nf_exps <- vector('list', sum(idx))
#      names(nf_exps) <- paste0(format(as.POSIXct(exp_spreadsheet$Date, format = '%d/%m/%Y'), '%Y%m%d'),
#                                exp_spreadsheet$MATI.TARDA)[idx]
#      
#      for(exp in names(nf_exps)){
#           jsonexp <- read_json(paste0(path, exp, '/iTracks.json'))
#           nf_exps[[exp]] <- alpha(jsonexp, i)
#      }
#      nf_exps
# })
# Sys.time() - t0

t0 <- Sys.time()
NF_results <- lapply(1:1200, function(i){
        nf_exps <- vector('list', sum(idx))
        names(nf_exps) <- paste0(format(as.POSIXct(exp_spreadsheet$Date, format = '%d/%m/%Y'), '%Y%m%d'),
                                  exp_spreadsheet$MATI.TARDA)[idx]
        
        for(exp in names(nf_exps)){
                jsonexp <- read_json(paste0(path, exp, '/iTracks.json'))
                nf_exps[[exp]] <- alpha(jsonexp, i)
        }
        
        r <- vapply(seq_along(nf_exps), function(z){
                x <- moving_average(nf_exps[[z]]$x, 60, 30)
                y <- moving_average(nf_exps[[z]]$y, 60, 30)
                
                x[x >= 1] <- 0
                y[y >= 1] <- 0
                cor(x, y)
        }, numeric(1))
        c(mean = mean(r), sd = sd(r))
})
Sys.time() - t0

# NF_CORS <- lapply(NF_results, function(i){
#      
#      
#      r <- vapply(seq_along(i), function(z){
#           x <- moving_average(i[[z]]$x, 60, 30)
#           y <- moving_average(i[[z]]$y, 60, 30)
#           
#           x[x >= 1] <- 0
#           y[y >= 1] <- 0
#           cor(x, y)
#      }, numeric(1))
#      
#      c(mean = mean(r), sd = sd(r))
#      
# })


df <- do.call('rbind.data.frame', NF_results)
colnames(df) <- c('mean', 'sd')
df$t <- 1:nrow(df)
ggplot(data = df, aes(t, mean)) + 
     geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.3, color = 'black')+ 
     geom_point(size = 3)+ geom_line()+
     scale_x_continuous('Temporal buffer (frames)', breaks = seq(0, 60, 10))+
     scale_y_continuous('Average correlation', breaks = seq(0, 1, length.out = 6))+
     ggtitle('No food') + theme(plot.title = element_text(hjust = 0.5, size = 16))

save.image('/home/polfer/research/2022/ANTS/AnTracks/results/WS_det_sto_nf_DptEnt.RData')
