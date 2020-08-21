#####################################################################################################################################################
############################################################ Bipartite Network ######################################################################
#####################################################################################################################################################
library(igraph)
setwd("/home/frankaylward/bipartite_network")
viral_matrix <- read.table(file="vog_table.txt", header=T, row.names=1, check.names = F)
#viral_subset <- viral_matrix[,colSums(viral_matrix) > 10] 
viral_subset2 <- viral_subset[rowSums(viral_matrix) > 1,] #only include VOGs present in at least 2 genomes

set.seed=123
igr <- graph.incidence(viral_subset2, weighted = T)
l <- layout.fruchterman.reingold(igr,niter=5000)

features <- read.table("features.txt", header=T, sep="\t", check.names = F)
V(igr)$color <- as.character(features$color[match(V(igr)$name, features$node)])

features$size <- features$size / 100 # scale the size so that nodes aren't too big. 
V(igr)$size <- as.numeric(features$size[match(V(igr)$name, features$node)])
V(igr)$size[!V(igr)$type] <- 0.5
#V(igr)$size <- as.numeric(V(igr)$size)
#V(igr)$size[V(igr)$type] <- 4

plot.igraph(igr,vertex.label=NA, layout=l, edge.color="grey90", vertex.frame.color="grey30", edge.width=0.2)

