library(dtw)
library(data.table)

source('~/research/AnTracks/src/simulations.R')
load('~/research/AnTracks/data/det_ref_seqs.RData')

dt <- read_FilteredOutput('~/research/AutomatAnts/results/pheromone/', filter = 'pheromone')

aligned_dt <- align_sequences(dt, sim_type = 'pheromone')
pheromone_dtw <- mass_dtw(aligned_dt, det_NS, align = FALSE, only.average = FALSE)

save(pheromone_dtw, file = '~/research/AutomatAnts/results/pheromone/outputs/pheromone_dtw.RData')