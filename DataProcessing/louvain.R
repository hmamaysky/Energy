### 	Function: louvain algorithm for clustering the cosine similarity matrix
### 	Input: cosine similarity matrix, Output: Clusters 

install.packages('igraph')

rm(list = ls())
setwd("/Users/Roya/Desktop/cosine/louvain")

#read the cosine file
df <- read.table(file="cosine.csv", header = TRUE,  row.names = 1, sep = ',')
colnames(df) <- rownames(df)
names <- rownames(df)

#iterate over different shuffling of words to get the best modularity
s <- 0
s.names <- names

for (i in 1:1000){
names.shuffle <- sample(names)
df.shuffle <- df[names.shuffle,names.shuffle] 
m <- as.matrix(df.shuffle)
mode(m) <- "numeric"
library(igraph)
ig <- graph.adjacency(m, mode="undirected", weighted=TRUE,diag = FALSE)
c1 = cluster_louvain(ig)
if (modularity(c1) > s) {s = modularity(c1)
s.names = names.shuffle}
}


#get the comunities for the best permutation of words
df.final <- df[s.names,s.names] 
m <- as.matrix(df.final)
mode(m) <- "numeric"
ig <- graph.adjacency(m, mode="undirected", weighted=TRUE,diag = FALSE)
c1 = cluster_louvain(ig)

#see modularity and sizes of clusters
print(modularity(c1))
print(sizes(c1))

#generate communities in different .csv files
for (i in 1:length(c1)){
  write.table(c1[i],paste("group", i, ".csv", sep = ""), quote = FALSE, row.names = FALSE)
}


