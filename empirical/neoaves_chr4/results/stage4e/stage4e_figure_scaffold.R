#!/usr/bin/env Rscript
required <- c("T0","NONCHR4","T_PNAS","T_STRUCT_PRIMARY","T_BLOCK_500K")
tree_dir <- "empirical/neoaves_chr4/results/stage4e/trees"
missing <- required[!file.exists(file.path(tree_dir, paste0(required, ".nwk")))]
if (length(missing)) {
  stop(paste("PENDING_CLUSTER_OUTPUT:", paste(missing, collapse=", ")))
}
message("Stage 4E figure scaffolding is ready; implement final panels after tree validation.")
