source('/home/polfer/research/2022/ANTS/AnTracks/src/Experiment.R')
path <- '/home/polfer/Documents/Tracker_LEOV/ANTS/data/'
exps <- list.files(path)
load('~/research/2022/ANTS/AnTracks/data/det.RData')
load('~/research/2022/ANTS/AnTracks/data/sto.RData')

## Calculate nest entry and departure rates over time ##
alpha <- function(jsonexp, min_time = 0, min_length = 0){
     x <- numeric(21600)
     y <- numeric(21600)
     
     for(id in seq_along(jsonexp)){
          
          l <- length(jsonexp[[id]])
          mint <- jsonexp[[id]][[1]][[1]]
          maxt <- jsonexp[[id]][[l]][[1]]
          
          if(l > min_length && maxt > (mint + min_time)){
               
               if(maxt > 21600){
                    maxt <- 21600
               }
               if(mint > 21600){
                    mint <- 21600
               }
               x[mint] <- x[mint] + 1
               y[maxt] <- y[maxt] + 1
          }
          
     }
     
     df <- data.frame(x = x, y = y, 
                      xsmth = moving_average(x, t = 60, overlap = 30),
                      ysmth = moving_average(y, t = 60, overlap = 30))
     df$xsmth[df$xsmth >= 1] <- 0
     df$ysmth[df$ysmth >= 1] <- 0
     
     df
     
}

#### ------ =+=+= SOFT CRITERIA =+=+= ------ ####

#### +++ DETERMINIST EXPS +++ ####

idx <- exp_spreadsheet$CONDITION == 'DET'
det_exps <- vector('list', sum(idx))
names(det_exps) <- paste0(format(as.POSIXct(exp_spreadsheet$Date, format = '%d/%m/%Y'), '%Y%m%d'),
                          exp_spreadsheet$MATI.TARDA)[idx]
for(exp in names(det_exps)){
     jsonexp <- read_json(paste0(path, exp, '/iTracks.json'))
     det_exps[[exp]] <- alpha(jsonexp, 30)
}

det_plots <- lapply(det_exps, function(i){
     ggplot(data = data.frame(t = 1:21600, x = i$xsmth), aes(t, x)) + geom_path()
})
ggarrange(plotlist = det_plots, labels = names(det_exps))

det_plots <- lapply(det_exps, function(i){
     ggplot(data = data.frame(t = 1:21600, y = i$ysmth), aes(t, y)) + geom_path()
})
ggarrange(plotlist = det_plots, labels = names(det_exps))

det_corplots <- lapply(det_exps, function(i){
     c <- cor(i$xsmth, i$ysmth)
     ggplot(data = data.frame(x = i$xsmth, y = i$ysmth), aes(x, y)) + geom_point()+
          annotate('text', label = paste0('R = ', round(c, 2)),
                   x = max(i$xsmth) / 5, y = 3* max(i$ysmth)/4, size = 6)+
          geom_smooth(method = 'lm', formula = y ~ x)
})
ggarrange(plotlist = det_corplots, labels = names(det_exps))

det_diff <- lapply(det_exps, function(i){
     ggplot(data = data.frame(x = i$xsmth, y = i$ysmth, t = 1:21600), aes(t, x - y)) + geom_path()
})
ggarrange(plotlist = det_diff, labels = names(det_exps))

#### +++ STOCHASTIC EXPS +++ ####

idx <- exp_spreadsheet$CONDITION == 'STOCH'
sto_exps <- vector('list', sum(idx))
names(sto_exps) <- paste0(format(as.POSIXct(exp_spreadsheet$Date, format = '%d/%m/%Y'), '%Y%m%d'),
                          exp_spreadsheet$MATI.TARDA)[idx]

for(exp in names(sto_exps)){
     jsonexp <- read_json(paste0(path, exp, '/iTracks.json'))
     sto_exps[[exp]] <- alpha(jsonexp)
}

sto_plots <- lapply(sto_exps, function(i){
     ggplot(data = data.frame(t = 1:21600, x = i$xsmth), aes(t, x)) + geom_path()
})
ggarrange(plotlist = sto_plots, labels = names(sto_exps))

sto_plots <- lapply(sto_exps, function(i){
     ggplot(data = data.frame(t = 1:21600, y = i$ysmth), aes(t, y)) + geom_path()
})
ggarrange(plotlist = sto_plots, labels = names(sto_exps))

sto_corplots <- lapply(sto_exps, function(i){
     c <- cor(i$xsmth, i$ysmth)
     ggplot(data = data.frame(x = i$xsmth, y = i$ysmth), aes(x, y)) + geom_point()+
          annotate('text', label = paste0('R = ', round(c, 2)),
                   x = max(i$xsmth) / 5, y = 3* max(i$ysmth)/4, size = 6)+
          geom_smooth(method = 'lm', formula = y ~ x)
})
ggarrange(plotlist = sto_corplots, labels = names(sto_exps))

sto_diff <- lapply(sto_exps, function(i){
     ggplot(data = data.frame(x = i$xsmth, y = i$ysmth, t = 1:21600), aes(t, x - y)) + geom_path()
})
ggarrange(plotlist = sto_diff, labels = names(sto_exps))

#### +++ CONTROL (NO FOOD) EXPS +++ ####
idx <- exp_spreadsheet$CONDITION == 'NO FOOD'
nf_exps <- vector('list', sum(idx))
names(nf_exps) <- paste0(format(as.POSIXct(exp_spreadsheet$Date, format = '%d/%m/%Y'), '%Y%m%d'),
                          exp_spreadsheet$MATI.TARDA)[idx]

for(exp in names(nf_exps)){
     jsonexp <- read_json(paste0(path, exp, '/iTracks.json'))
     nf_exps[[exp]] <- alpha(jsonexp)
}

nf_plots <- lapply(nf_exps, function(i){
     ggplot(data = data.frame(t = 1:21600, x = i$xsmth), aes(t, x)) + geom_path()
})
ggarrange(plotlist = nf_plots, labels = names(nf_exps))

nf_plots <- lapply(nf_exps, function(i){
     ggplot(data = data.frame(t = 1:21600, y = i$ysmth), aes(t, y)) + geom_path()
})
ggarrange(plotlist = nf_plots, labels = names(nf_exps))

nf_corplots <- lapply(nf_exps, function(i){
     c <- cor(i$xsmth, i$ysmth)
     ggplot(data = data.frame(x = i$xsmth, y = i$ysmth), aes(x, y)) + geom_point()+
          annotate('text', label = paste0('R = ', round(c, 2)),
                   x = max(i$xsmth) / 5, y = 3* max(i$ysmth)/4, size = 6)+
          geom_smooth(method = 'lm', formula = y ~ x)
})
ggarrange(plotlist = nf_corplots, labels = names(nf_exps))

nf_diff <- lapply(nf_exps, function(i){
     ggplot(data = data.frame(x = i$xsmth, y = i$ysmth, t = 1:21600), aes(t, x - y)) + geom_path()
})
ggarrange(plotlist = nf_diff, labels = names(nf_exps))



#### ------ =+=+= HARD CRITERIA =+=+= ------ ####

#### +++ DETERMINIST EXPS +++ ####

idx <- exp_spreadsheet$CONDITION == 'DET'
det_exps <- vector('list', sum(idx))
names(det_exps) <- paste0(format(as.POSIXct(exp_spreadsheet$Date, format = '%d/%m/%Y'), '%Y%m%d'),
                          exp_spreadsheet$MATI.TARDA)[idx]
for(exp in names(det_exps)){
     jsonexp <- read_json(paste0(path, exp, '/iTracks.json'))
     det_exps[[exp]] <- alpha(jsonexp, min_time = 30, min_length = 10)
}

det_plots <- lapply(det_exps, function(i){
     ggplot(data = data.frame(t = 1:21600, x = i$xsmth), aes(t, x)) + geom_path()
})
ggarrange(plotlist = det_plots, labels = names(det_exps))

det_plots <- lapply(det_exps, function(i){
     ggplot(data = data.frame(t = 1:21600, y = i$ysmth), aes(t, y)) + geom_path()
})
ggarrange(plotlist = det_plots, labels = names(det_exps))

det_corplots <- lapply(det_exps, function(i){
     c <- cor(i$xsmth, i$ysmth)
     ggplot(data = data.frame(x = i$xsmth, y = i$ysmth), aes(x, y)) + geom_point()+
          annotate('text', label = paste0('R = ', round(c, 2)),
                   x = max(i$xsmth) / 5, y = 3* max(i$ysmth)/4, size = 6)+
          geom_smooth(method = 'lm', formula = y ~ x)
})
ggarrange(plotlist = det_corplots, labels = names(det_exps))

det_diff <- lapply(seq_along(det_exps), function(i){
     ggplot(data = data.frame(x = det_exps[[i]]$xsmth, y = det_exps[[i]]$ysmth, t = 1:21600),
            aes(t, x - y)) + geom_path()+
          geom_food(do.call('rbind', det[[i]]@food)$t, complete = T)
})
ggarrange(plotlist = det_diff, labels = names(det_exps))

#### +++ STOCHASTIC EXPS +++ ####

idx <- exp_spreadsheet$CONDITION == 'STOCH'
sto_exps <- vector('list', sum(idx))
names(sto_exps) <- paste0(format(as.POSIXct(exp_spreadsheet$Date, format = '%d/%m/%Y'), '%Y%m%d'),
                          exp_spreadsheet$MATI.TARDA)[idx]

for(exp in names(sto_exps)){
     jsonexp <- read_json(paste0(path, exp, '/iTracks.json'))
     sto_exps[[exp]] <- alpha(jsonexp, min_time = 30, min_length = 10)
}

sto_plots <- lapply(sto_exps, function(i){
     ggplot(data = data.frame(t = 1:21600, x = i$xsmth), aes(t, x)) + geom_path()
})
ggarrange(plotlist = sto_plots, labels = names(sto_exps))

sto_plots <- lapply(sto_exps, function(i){
     ggplot(data = data.frame(t = 1:21600, y = i$ysmth), aes(t, y)) + geom_path()
})
ggarrange(plotlist = sto_plots, labels = names(sto_exps))

sto_corplots <- lapply(sto_exps, function(i){
     c <- cor(i$xsmth, i$ysmth)
     ggplot(data = data.frame(x = i$xsmth, y = i$ysmth), aes(x, y)) + geom_point()+
          annotate('text', label = paste0('R = ', round(c, 2)),
                   x = max(i$xsmth) / 5, y = 3* max(i$ysmth)/4, size = 6)+
          geom_smooth(method = 'lm', formula = y ~ x)
})
ggarrange(plotlist = sto_corplots, labels = names(sto_exps))

sto_diff <- lapply(seq_along(sto_exps), function(i){
     ggplot(data = data.frame(x = sto_exps[[i]]$xsmth, y = sto_exps[[i]]$ysmth, t = 1:21600),
            aes(t, x - y)) + geom_path()+
          geom_food(do.call('rbind', sto[[i]]@food)$t, complete = T)
})
ggarrange(plotlist = sto_diff, labels = names(sto_exps))

#### +++ CONTROL (NO FOOD) EXPS +++ ####
idx <- exp_spreadsheet$CONDITION == 'NO FOOD'
nf_exps <- vector('list', sum(idx))
names(nf_exps) <- paste0(format(as.POSIXct(exp_spreadsheet$Date, format = '%d/%m/%Y'), '%Y%m%d'),
                         exp_spreadsheet$MATI.TARDA)[idx]

for(exp in names(nf_exps)){
     jsonexp <- read_json(paste0(path, exp, '/iTracks.json'))
     nf_exps[[exp]] <- alpha(jsonexp)
}

nf_plots <- lapply(nf_exps, function(i){
     ggplot(data = data.frame(t = 1:21600, x = i$xsmth), aes(t, x)) + geom_path()
})
ggarrange(plotlist = nf_plots, labels = names(nf_exps))

nf_plots <- lapply(nf_exps, function(i){
     ggplot(data = data.frame(t = 1:21600, y = i$ysmth), aes(t, y)) + geom_path()
})
ggarrange(plotlist = nf_plots, labels = names(nf_exps))

nf_corplots <- lapply(nf_exps, function(i){
     c <- cor(i$xsmth, i$ysmth)
     ggplot(data = data.frame(x = i$xsmth, y = i$ysmth), aes(x, y)) + geom_point()+
          annotate('text', label = paste0('R = ', round(c, 2)),
                   x = max(i$xsmth) / 5, y = 3* max(i$ysmth)/4, size = 6)+
          geom_smooth(method = 'lm', formula = y ~ x)
})
ggarrange(plotlist = nf_corplots, labels = names(nf_exps))

nf_diff <- lapply(nf_exps, function(i){
     ggplot(data = data.frame(x = i$xsmth, y = i$ysmth, t = 1:21600), aes(t, x - y)) + geom_path()
})
ggarrange(plotlist = nf_diff, labels = names(nf_exps))