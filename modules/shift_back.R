args <- commandArgs(trailingOnly = TRUE)

.libPaths=(args[3])
#library(dplyr)
library(stringr)
#library(tidyr)

shift_back = function(x) {
	if (x < 8570) {
          return(x + 8000)
        } else {
          return (x - 8569)
        }
      }

#control_region_shifted = read.table("control_region_shifted.tsv", header=T)
control_region_shifted = read.table(args[1], header=T)
shifted_back = sapply(control_region_shifted[,"pos"], shift_back)
control_region_shifted[,"pos"] = shifted_back

beginning = subset(control_region_shifted, control_region_shifted[,'pos']<8000)
end = subset(control_region_shifted, control_region_shifted[,'pos']>8000)

#non_control_region = read.table("non_control_region.tsv", header=T)
non_control_region = read.table(args[2], header=T)
combined_table = rbind(beginning, non_control_region, end)
write.table(combined_table, "per_base_coverage.tsv", row.names=F, col.names=T, quote=F, sep="\t")
