node_idx.nodes <- function(obj, row){
        as.integer(rownames(obj$data[row, ]))
}