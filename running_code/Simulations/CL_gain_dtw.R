library(dtw)
library(data.table)

source('~/research/AnTracks/src/simulations.R')
load('~/research/AnTracks/data/det_ref_seqs.RData')

dt <- read_gainsOutput('~/research/AutomatAnts/results/gains/')

aligned_dt <- align_sequences(dt, sim_type = 'gain')
gain_dtw <- pdtw(aligned_dt, det_NS, align = FALSE, only.average = FALSE)

save(gain_dtw, file = '~/research/AutomatAnts/results/gains/outputs/gain_dtw.RData')