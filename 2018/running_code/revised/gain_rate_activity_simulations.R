library(ggplot2)

#### +++ CLASS DEFINITION +++ ####
setClass('Model', representation(
  agents = 'list',
  params = 'list',
  default_params = 'list',
  .__classVersion__ = 'character'
), prototype = list(
  params = list(),
  default_params = list(
  Jij = data.frame(aa = 1, ab = 1, ac = 1,
                    ba = 1, bb = 1, bc = 1,
                    ca = 1, cb = 1, cc = 1),
  g = 0.8,
  rate = 10,
  threshold = 0.5,
  max_ints = 4
  ),
  .__classVersion__ = '1.0'
))

setGeneric('init_params', function(obj){
  standardGeneric('init_params')
})
setMethod('init_params', 'Model', function(obj){
  params_left <- names(obj@default_params)[!names(obj@default_params) %in% names(obj@params)]
  for(i in params_left){
    obj@params[[i]] <- obj@default_params[[i]]
  }
  obj
})

setGeneric('run_serio', function(obj, tmax = 250){
  standardGeneric('run_serio')
})

setMethod('run_serio', 'Model', function(obj, tmax = 250){
  
  obj <- init_params(obj)
  obj@agents <- lapply(obj@agents, update, t = 0)
  
  
  rate <- obj@params$rate
  threshold <- obj@params$threshold
  individuals <- seq_along(obj@agents)
  max_ints <- obj@params$max_ints
  
  rates <- sapply(obj@agents, function(i) abs(i@Si))
  ratesum <- sum(rates)
  ratenorm <- rates / ratesum
  
  time <- 0
  
  out <- c()
  
  
  while(time < tmax){
    sampled_t <- rexp(1, ratesum)
    ind <- sample(seq_along(rates), size = 1, prob = ratenorm)
    I <- obj@agents[[ind]]
    if(ind %in% out){
      if(length(out) > 4){
        ints <- sample(out[out != ind], size = sample(1:4, size = 1), replace = F)
      } else if(length(out) > 1){
        ints <- sample(out[out != ind], size = sample(seq_len(length(out)-1), size = 1), replace = F)
      } else {
        ints <- c()
      }
      I <- activity(I, obj@agents[ints])
    } else if (I@Si > 1e-3){
      out <- c(out, ind)
    } else {
      I <- activity(I, obj@agents[sample(seq_along(rates)[-out], size = sample(1:4, size = 1), replace = F)])
    }

    I <- update(I, t = time)
    obj@agents[[ind]] <- I
    time <- time + sampled_t
  }
  obj
})

setGeneric('run', function(obj, iters = 250){
  standardGeneric('run')
})

setMethod('run', 'Model', function(obj, iters = 250){
  
  obj <- init_params(obj)
  obj@agents <- lapply(obj@agents, update, t = 0)
  
  
  rate <- obj@params$rate
  threshold <- obj@params$threshold
  individuals <- seq_along(obj@agents)
  max_ints <- obj@params$max_ints
  
  for(i in 1:iters){
    ord <- sample(individuals, size = length(individuals), replace = F)
    for(ii in ord){
      I <- obj@agents[[ii]]
      if(rexp(1, rate = rate) > threshold){
        nints <- sample(seq_along(max_ints), 1)
        ints <- sample(individuals[individuals != ii], size = nints, replace = F)
        I <- activity(I, obj@agents[ints])
      } else {
        I <- activity(I, c())
      }
      I <- update(I, t = i)
      obj@agents[[ii]] <- I
    }
  }
  obj
})

setGeneric('manipulated_run', function(obj, g = NULL, Si = NULL, rate = NULL, t = 100, iters = 250){
  standardGeneric('manipulated_run')
})

setMethod('manipulated_run', 'Model', function(obj, g = NULL, Si = NULL, rate = NULL, t = 100, iters = 250){
  
  obj <- init_params(obj)
  obj@agents <- lapply(obj@agents, update, t = 0)
  
  rate_post <- rate
  rate <- obj@params$rate
  threshold <- obj@params$threshold
  individuals <- seq_along(obj@agents)
  max_ints <- obj@params$max_ints
  
  for(i in 1:t){
    ord <- sample(individuals, size = length(individuals), replace = F)
    for(ii in ord){
      I <- obj@agents[[ii]]
      if(rexp(1, rate = rate) > threshold){
        nints <- sample(seq_along(max_ints), 1)
        ints <- sample(individuals[individuals != ii], size = nints, replace = F)
        I <- activity(I, obj@agents[ints])
      } else {
        I <- activity(I, c())
      }
      I <- update(I, t = i)
      obj@agents[[ii]] <- I
    }
  }
  if(!is.null(g) && is.numeric(g)){
    for(a in seq_along(obj@agents)){
      obj@agents[[a]]@g <- g
    }
  }
  
  if(!is.null(Si) && is.numeric(Si)){
    for(a in seq_along(obj@agents)){
      obj@agents[[a]]@Si <- Si
    }
  }
  
  if(!is.null(rate_post) && is.numeric(rate_post)){
    rate <- rate_post
  }
  
  for(i in (t+1):iters){
    ord <- sample(individuals, size = length(individuals), replace = F)
    for(ii in ord){
      I <- obj@agents[[ii]]
      if(rexp(1, rate = rate) > threshold){
        nints <- sample(seq_along(max_ints), 1)
        ints <- sample(individuals[individuals != ii], size = nints, replace = F)
        I <- activity(I, obj@agents[ints])
      } else {
        I <- activity(I, c())
      }
      I <- update(I, t = i)
      obj@agents[[ii]] <- I
    }
  }
  obj
})

setGeneric('plot', function(obj, id = 1){
  standardGeneric('plot')
})

setMethod('plot', 'Model', function(obj, id = 1){
  print(obj@params)
  
  df <- obj@agents[[id]]@df
  ggplot(data = df, aes(t, Si))+ geom_path() +
    theme_bw() + 
    scale_x_continuous('Time', breaks = seq(0, max(df$t), length.out = 6))+
    scale_y_continuous('Activity (Si)')+
    theme(axis.title = element_text(size  = 15, color = 'black'),
          axis.text = element_text(size = 15, color = 'black'),
          legend.text = element_text(size =15, color = 'black'),
          legend.title = element_text(size = 15, color = 'black'))
})

setGeneric('global_plot', function(obj){
  standardGeneric('global_plot')
})

setMethod('global_plot', 'Model', function(obj){
  print(obj@params)
  df <- do.call('rbind', lapply(obj@agents, function(i){df <- i@df; df$id <- i@id; df}))
  
  ggplot(data = df, aes(t, Si, group = id))+ geom_path() +
    theme_bw() + 
    scale_x_continuous('Time', breaks = seq(0, max(df$t), length.out = 6))+
    scale_y_continuous('Activity (Si)')+
    theme(axis.title = element_text(size  = 15, color = 'black'),
          axis.text = element_text(size = 15, color = 'black'),
          legend.text = element_text(size =15, color = 'black'),
          legend.title = element_text(size = 15, color = 'black'))
})

#### +++ CLASS DEFINITION +++ ####
setClass('Ant', representation(
  id = 'numeric',
  Si = 'numeric',
  state = 'character',
  g = 'numeric',
  df = 'data.frame',
  .__classVersion__ = 'character'
), prototype = list(
  Si = runif(1, -1, 1),
  state = 'a', 
  g = runif(1),
  .__classVersion__ = '1.0'
))

setGeneric('activity', function(obj, neighbors){
  standardGeneric('activity')
})

setMethod('activity', 'Ant', function(obj, neighbors){
  z <- sum(vapply(neighbors, function(i){
    Jij[[paste0(obj@state, i@state)]] * i@Si
  }, numeric(1)))
  obj@Si <- tanh(obj@g * (z + obj@Si))
  obj@g <- mean(vapply(neighbors, function(i) i@g, numeric(1)), obj@g)
  obj
})


# setMethod('activity', 'Ant', function(obj, neighbors){
#   z <- sum(vapply(neighbors, function(i){
#     Jij[[paste0(obj@state, i@state)]] * i@Si
#   }, numeric(1)))
#   obj@Si <- tanh(obj@g * (z + obj@Si))
#   obj
# })

setGeneric('init', function(obj){
  standardGeneric('init')
})

setMethod('init', 'Ant', function(obj){
  obj@Si = runif(1, -1, 1)
  # obj@g = runif(1)
  obj
})

setGeneric('update', function(obj, t){
  standardGeneric('update')
})

setMethod('update', 'Ant', function(obj, t){
  obj@df <- rbind(obj@df, c(obj@Si, t))
  colnames(obj@df) <- c('Si', 't')
  obj
})


comp_z <- function(rate, gain){
  
  m <- matrix(nrow = length(rate), ncol = length(gain))
  col <- 1
  row <- 1

  for(g in gain){
    for(r in rate){
      M <- new('Model', params = list(g = g, rate = r))
      M@agents <- lapply(1:25, function(i){
        new('Ant', g = M@params$g, Si = runif(1,  0, 1), id = i)
      })
      M <- run(M)
      v <- as.numeric(sapply(M@agents, function(i) i@df$Si))
      m[row, col] <- mean(v)
      row <- row + 1
    }
    row <- 1
    col <- col + 1
  }
  # cbind(0, m)
  m
}

gain <- seq(0, 1, 0.05)
rate <- seq(1/20, 4, length.out = 20)
Si <- comp_z(rate = rate, gain = gain)

persp(rate, gain, Si, col = 'lightblue', shade = 0.5, ticktype = 'detailed', theta = 60, phi = 45)

t0 <- Sys.time()
gain <- seq(0, 1.5, 0.01)
rate <- seq(1/50, 10, length.out = 151)
Si <- comp_z(rate = rate, gain = gain)
Sys.time() - t0

jet.colors <- viridis::viridis
zfacet <- scan$Si[-1, -1] + scan$Si[-1, -ncol(scan$Si)] +
  scan$Si[-nrow(scan$Si), -1] + scan$Si[-nrow(scan$Si), -ncol(scan$Si)]
facetcol <- cut(zfacet, 100)

persp(rate, gain, Si, col = 'lightblue', shade = 0.1, ticktype = 'detailed', theta = 60, phi = 45)
persp(rate, gain, Si, col = viridis::viridis(100)[facetcol],
      ticktype = 'detailed', theta = 300, phi = 45, lwd = 0.1)


#################### 

M <- new('Model', params = list(g = 0.75, rate = 1))
M@agents <- lapply(1:100, function(i){
  new('Ant', g = rnorm(1, 0.5, 0.1), Si = runif(1, 0, 1), id = i)
})
M <- run(M)

plot(M, 3)
global_plot(M)


result <- c()
for(i in seq(0, 1, 0.05)){
  M <- new('Model', params = list(g = i, rate = 1))
  M@agents <- lapply(1:50, function(ii){
    new('Ant', g = M@params$g, Si = runif(1, 0, 1), id = ii)
  })
  M <- run(M)
  # v <- as.numeric(sapply(M@agents, function(i) i@df$Si))
  # result <- rbind.data.frame(result, c(mean(v), sd(v), i))
  tmp <- do.call('rbind', lapply(M@agents, function(i) i@df$Si))
  result <- rbind.data.frame(result, cbind(colMeans(tmp), apply(tmp, 2, sd), seq(0, 250), i))

}
# colnames(result) <- c('mean', 'sd', 'g')
colnames(result) <- c('mean', 'sd', 't', 'g')

ggplot(data = result, aes(g, mean)) + geom_ribbon(aes(x = g, ymin = mean - sd, ymax = mean + sd),
                                                  color = 'black', fill = 'grey80', alpha = 0.3)+
  geom_point(size = 2.5)+ 
  scale_x_continuous('Gain (g)', breaks = seq(0, 1, length.out = 5))+
  scale_y_continuous('Average activity (<Si>)', breaks = seq(0, 1, length.out = 5), limits = c(NA, NA))+
  theme_bw() + 
  theme(axis.title = element_text(size  = 15, color = 'black'),
        axis.text = element_text(size = 15, color = 'black'),
        legend.text = element_text(size =15, color = 'black'),
        legend.title = element_text(size = 15, color = 'black'))

ggplot(data = result, aes(t, mean, color = g, group = g)) + 
  geom_path(size = 2.5)+ 
  scale_x_continuous('Time (steps)')+
  scale_y_continuous('Average activity (<Si>)', breaks = seq(0, 1, length.out = 5), limits = c(NA, NA))+
  scale_color_viridis_c('Gain (g)')+
  theme_bw() + 
  theme(axis.title = element_text(size  = 15, color = 'black'),
        axis.text = element_text(size = 15, color = 'black'),
        legend.text = element_text(size =15, color = 'black'),
        legend.title = element_text(size = 15, color = 'black'))



print()



M <- new('Model', params = list(g = 1, rate = 25))
M@agents <- lapply(1:100, function(i){
  new('Ant', g = rnorm(1, 0.5, 0.2), Si = runif(1, 0, 1), id = i)
})
M <- manipulated_run(M, rate = 1, t = 100, iters = 500)

plot(M, 3)
global_plot(M)


result <- c()
for(i in seq(0, 1, 0.05)){
  M <- new('Model', params = list(g = i, rate = 1))
  M@agents <- lapply(1:50, function(ii){
    new('Ant', g = M@params$g, Si = runif(1, 0, 1), id = ii)
  })
  M <- run(M)
  # v <- as.numeric(sapply(M@agents, function(i) i@df$Si))
  # result <- rbind.data.frame(result, c(mean(v), sd(v), i))
  tmp <- do.call('rbind', lapply(M@agents, function(i) i@df$Si))
  result <- rbind.data.frame(result, cbind(colMeans(tmp), apply(tmp, 2, sd), seq(0, 250), i))
  
}
# colnames(result) <- c('mean', 'sd', 'g')
colnames(result) <- c('mean', 'sd', 't', 'g')

ggplot(data = result, aes(g, mean)) + geom_ribbon(aes(x = g, ymin = mean - sd, ymax = mean + sd),
                                                  color = 'black', fill = 'grey80', alpha = 0.3)+
  geom_point(size = 2.5)+ 
  scale_x_continuous('Gain (g)', breaks = seq(0, 1, length.out = 5))+
  scale_y_continuous('Average activity (<Si>)', breaks = seq(0, 1, length.out = 5), limits = c(NA, NA))+
  theme_bw() + 
  theme(axis.title = element_text(size  = 15, color = 'black'),
        axis.text = element_text(size = 15, color = 'black'),
        legend.text = element_text(size =15, color = 'black'),
        legend.title = element_text(size = 15, color = 'black'))





mfp <- function(N, R = 50, L = 2000 * 1000, v = 12){
  lambda <- L / (2 * N * R)
  lambda / v 
}

mfp <- function(N, R = 50, L = 1766 * 50, v = 12){
  lambda <- L / (2 * N * R)
  lambda * v 
}

plot.default(as.numeric(sapply(1:20, mfp)))





M <- new('Model', params = list(g = 1, rate = 25))
M@agents <- lapply(1:100, function(i){
  new('Ant', g = rnorm(1, 0.5, 0.2), Si = runif(1, 0, 1), id = i)
})
M <- run_serio(M, 2000)

global_plot(M)

# 
# A <- new('Ant', Si = runif(1, -1, 1), g = 0.25)
# B <- new('Ant', Si = runif(1, -1, 1), g = 0.25)
# C <- new('Ant', Si = runif(1, -1, 1), g = 0.25)
# D <- new('Ant', Si = runif(1, -1, 1), g = 0.25)
# E <- new('Ant', Si = runif(1, -1, 1), g = 0.25)
# 
# rate <- 1/20
# individuals <- c('A', 'B', 'C', 'D', 'E')
# 
# for(i in 1:250){
#   ord <- sample(individuals, size = 5, replace = F)
#   for(s in ord){
#     I <- get(s)
#     if(rexp(1, rate = rate) > 0.5){
#       nints <- sample(1:4, 1)
#       ints <- sample(individuals[individuals != s], size = nints, replace = F)
#       I <- activity(I, lapply(ints, get))
#     } else {
#       I <- activity(I, c())
#     }
#     I <- update(I, t = i)
#     assign(s, I)
#   }
# }
# 
# data_rate1_g0.8 <- rbind(A@df, B@df, C@df, D@df, E@df)
# data_rate1_g0.8$Ant <- as.character(sapply(individuals, function(i){rep(i, 10000)}))
# ggplot(data = data_rate1_g0.8, aes(t, Si)) + facet_wrap(~ Ant) + geom_path()
# 
# data_rate5_g0.8 <- rbind(A@df, B@df, C@df, D@df, E@df)
# data_rate5_g0.8$Ant <- as.character(sapply(individuals, function(i){rep(i, 10000)}))
# ggplot(data = data_rate5_g0.8, aes(t, Si)) + facet_wrap(~ Ant) + geom_path()
# 
# data_rate1_g0.5 <- rbind(A@df, B@df, C@df, D@df, E@df)
# data_rate1_g0.5$Ant <- as.character(sapply(individuals, function(i){rep(i, 10000)}))
# ggplot(data = data_rate1_g0.5, aes(t, Si)) + facet_wrap(~ Ant) + geom_path()
# 
# data_rate0.5_g0.25 <- rbind(A@df, B@df, C@df, D@df, E@df)
# data_rate0.5_g0.25$Ant <- as.character(sapply(individuals, function(i){rep(i, 1000)}))
# ggplot(data = data_rate0.5_g0.25, aes(t, Si)) + facet_wrap(~ Ant) + geom_path()
# 
# data_rate0.2_g0.25 <- rbind(A@df, B@df, C@df, D@df, E@df)
# data_rate0.2_g0.25$Ant <- as.character(sapply(individuals, function(i){rep(i, 1000)}))
# ggplot(data = data_rate0.2_g0.25, aes(t, Si)) + facet_wrap(~ Ant) + geom_path()
# 
# data_rate0.05_g0.25 <- rbind(A@df, B@df, C@df, D@df, E@df)
# data_rate0.05_g0.25$Ant <- as.character(sapply(individuals, function(i){rep(i, 1000)}))
# ggplot(data = data_rate0.05_g0.25, aes(t, Si)) + facet_wrap(~ Ant) + geom_path()
