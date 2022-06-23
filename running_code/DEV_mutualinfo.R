t0 <- Sys.time()
mi <- mutinformation(as.data.frame(m[1:1000, ]))
Sys.time()- t0

mi_func <- function(X){
        s0 <- c(-1, -1)
        s1 <- c(-1, 1)
        s2 <- c(1, -1)
        s3 <- c(1, 1)
        
        
}

md1 <- coords2matrix(det[[1]])

prs <- t(combn(colnames(md1), 2))
head(prs)

det_phase <- t(vapply(det, function(i){
        x <- food_detection(i)
        c(min(x, na.rm = T), max(x, na.rm = T))
}, numeric(2)))

for(i in seq_along(det)){
        det[[i]]$mutual_info_p1 <- mutual_info(det[[i]], 1:det_phase[i, 1])
        det[[i]]$mutual_info_p2 <- mutual_info(det[[i]], (det_phase[i, 1] + 1):det_phase[i, 2])
        det[[i]]$mutual_info_p3 <- mutual_info(det[[i]], (det_phase[i, 2]+1):max(det[[i]]$data$Frame))
}






d <- exp_spreadsheet[exp_spreadsheet$CONDITION == 'STOCH', c('Date', 'MATI.TARDA')][1:10, ]
d <- paste(format(strptime(d[, 1], '%d/%m/%Y'), format = '%Y%m%d'), d[, 2], sep='')

for(i in seq_along(sto)){
        sto[[i]]$date <- d[i]
}

for(i in seq_along(sto)){
        food <- get_foodPatches(sto[[i]])
        sto[[i]]$food <- food[1:2]
}

sto_phase <- t(vapply(sto, function(i){
        x <- food_detection(i)
        c(min(x, na.rm = T), max(x, na.rm = T))
}, numeric(2)))

for(i in seq_along(sto)){
        sto[[i]]$mutual_info_p1 <- mutual_info(sto[[i]], 1:sto_phase[i, 1])
        sto[[i]]$mutual_info_p2 <- mutual_info(sto[[i]], (sto_phase[i, 1] + 1):sto_phase[i, 2])
        sto[[i]]$mutual_info_p3 <- mutual_info(sto[[i]], (sto_phase[i, 2]+1):max(sto[[i]]$data$Frame))
}

######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################

exp_condition <- c('')

source('G:/research/MoveAnts/code/config.R')
source('C:/Users/POL/Desktop/gants_current/code/dev/dev_functions.R')

hex <- hex90
colnames(hex) <- c('x', 'y')
h <- hex90[hex90$rotY > 980, ]
colnames(h) <- c('x', 'y')
source('G:/research/2022/AnTracks/src/generic.R')
source('G:/research/2022/AnTracks/src/coords.R')

rm(list = ls()[grepl('json', ls())], visual, colony, current_dir, exp_condition, exps, conditionL)


load('G:/research/2022/AnTracks/results/det_coords.RData')


mi_func <- function(X){
        dims <- dim(X)
        m <- dims[1]
        n <- dims[2]
        
        M <- matrix(0, nrow = n, ncol = n, dimnames = list(colnames(X), colnames(X)))
        
        entropy <- function(x){x * log(x)}
        
        # combinations <- combn(1:)
        
        for(i in 1:(n-1)){
                for(j in 2:n){
                        
                        px <- c(sum(X[, i] == -1), sum(X[, i] == 1)) / m
                        py <- c(sum(X[, j] == -1), sum(X[, j] == 1)) /m 
                        
                        if(any(c(px, py) == 0)){
                                
                                I <- 0
                                
                        } else {
                                
                                p <- c(
                                        p0 = sum(X[, i] == -1 & X[, j] == -1) / m,
                                        p1 = sum(X[, i] == -1 & X[, j] == 1) / m,
                                        p2 = sum(X[, i] == 1 & X[, j] == -1) / m,
                                        p3 = sum(X[, i] == 1 & X[, j] == 1) / m
                                )
                                
                                if(any(p == 0)){
                                        
                                        I <- 0
                                } else {
                                        
                                        Hx <- -sum(vapply(px, entropy, numeric(1)))
                                        Hy <- -sum(vapply(py, entropy, numeric(1)))
                                        H <- -sum(vapply(p, entropy, numeric(1)))
                                        
                                        
                                        I <- Hx + Hy - H
                                }
                                
                        }
                        M[i, j] <- I
                        
                }
        }
        M
        
}


z <- numeric(length(det[[1]]$segments$o))
s <- det[[1]]$segments[, c('o', 'd')]
n <- colnames(I)
for(i in n){
        idx <- which(s$o == i)
        sb <- s[idx, ]
        tmp <- numeric(nrow(sb))
        for(x in seq_len(nrow(sb))){
                r <- I[n == sb$o[x], n == sb$d[x]]
                if(length(r)){
                        tmp[x] <- r
                }
        }
        if(length(idx) != length(tmp)){
                print(i)
        }
        z[idx] <- tmp
}
z
