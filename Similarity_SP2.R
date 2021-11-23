###### This script is to compute the similarity for BRDs with SP=2 or SP=3
library(igraph)
### Firstly to try with small data set, only in BRD-BRD interaction, ATAD2 friends with SP=2
ATAD2_SP2 <- dist.BRD[,colnames(dist.BRD)=="ATAD2"]
ATAD2_SP2_df <- as.data.frame(ATAD2_SP2[ATAD2_SP2==2])
ATAD2_SP2_list <- toString(rownames(ATAD2_SP2_df))
### Then using data frame to do 
BRD_SP2_df <- data.frame(matrix(ncol = 2, nrow = 42))
colnames(BRD_SP2_df) <- c("BRD", "Neighbors")


for (i in 1:42) {
  BRD_SP2_df$BRD[i]= colnames(dist.BRD)[i]
  a <- dist.BRD[,i]
  aa <- as.data.frame(a[a==2])
  BRD_SP2_df$Neighbors[i] <- toString(rownames(aa))
}

aa=unlist(strsplit(BRD_SP2_df$Neighbors[1], ","))
bb=unlist(strsplit(BRD_SP2_df$Neighbors[2], ","))
cc=intersect(aa,bb); dd=union(aa, bb)
length(cc)

BRD_SP2_similarity <- data.frame(matrix(ncol = 42, nrow = 42))
rownames(BRD_SP2_similarity) = colnames(dist.BRD); colnames(BRD_SP2_similarity) = colnames(dist.BRD)

for(i in 1:42) {
  for(j in 1:42) {
    aa=unlist(strsplit(BRD_SP2_df$Neighbors[i], ","))
    bb=unlist(strsplit(BRD_SP2_df$Neighbors[j], ","))
    cc=intersect(aa,bb)
    dd=union(aa, bb)
    BRD_SP2_similarity[i,j] = length(cc)/length(dd)
  }
}


###### Do it using whole human data #####
E(Human_igraph_unweighted)$weight; E(Human_igraph_unweighted)$weight = 1
dist_human_BRD <- distances(Human_igraph_unweighted, v = V(Human_igraph_unweighted), 
                            to = V(Human_igraph_unweighted)[is.element(V(Human_igraph_unweighted)$name, BRD_list$V1)], mode = "all", weights = NA, algorithm = "automatic")

colnames(dist_human_BRD) ##!!!

ALL_BRD_SP2_df <- data.frame(matrix(ncol = 2, nrow = 42))
colnames(ALL_BRD_SP2_df) <- c("BRD", "Neighbors")


for (i in 1:42) {
  ALL_BRD_SP2_df$BRD[i]= colnames(dist_human_BRD)[i]
  a <- dist_human_BRD[,i]
  aa <- as.data.frame(a[a==2])
  ALL_BRD_SP2_df$Neighbors[i] <- toString(rownames(aa))
}

ALL_BRD_SP2_similarity <- data.frame(matrix(ncol = 42, nrow = 42))
rownames(ALL_BRD_SP2_similarity) = colnames(dist_human_BRD); colnames(ALL_BRD_SP2_similarity) = colnames(dist_human_BRD)

for(i in 1:42) {
  for(j in 1:42) {
    aa=unlist(strsplit(ALL_BRD_SP2_df$Neighbors[i], ","))
    bb=unlist(strsplit(ALL_BRD_SP2_df$Neighbors[j], ","))
    cc=intersect(aa,bb)
    dd=union(aa, bb)
    ALL_BRD_SP2_similarity[i,j] = length(cc)/length(dd)
  }
}

### Plot
col.order <- as.character(BRD_level)
row.order <- as.character(BRD_level)
ALL_BRD_SP2_similarity_ordered <- ALL_BRD_SP2_similarity
ALL_BRD_SP2_similarity_ordered <- as.matrix(ALL_BRD_SP2_similarity_ordered)
ALL_BRD_SP2_similarity_ordered <- ALL_BRD_SP2_similarity_ordered[row.order,col.order]
# ALL_BRD_SP2_similarity_ordered[upper.tri(ALL_BRD_SP2_similarity_ordered)] <- NA

dev.off()
graphics.off()

col <- colorRampPalette(c("ivory2", "lightgoldenrod3", "lightsteelblue","olivedrab1", "springgreen")) #2edit
breaks <- c(0, 0.1, 0.2, 0.5, 0.9, 1)
library(gplots)
gc()
heatmap.2(ALL_BRD_SP2_similarity_ordered, Colv=NULL, Rowv= NULL, col=col, 
          breaks=breaks,
          trace= "none", denscol="black", keysize=1.5, 
          RowSideColors = rep(FuncGroup_summary$color, times = FuncGroup_summary$Num_BRDs), main = "Similarity Heatmap of Bromodomain-containing Proteins(SP=2)", xlab = "BRDs")




hist(ALL_BRD_SP2_similarity_ordered); median(ALL_BRD_SP2_similarity_ordered, na.rm=T)

###### SP=3
ALL_BRD_SP3_df <- data.frame(matrix(ncol = 2, nrow = 42))
colnames(ALL_BRD_SP3_df) <- c("BRD", "Neighbors")


for (i in 1:42) {
  ALL_BRD_SP3_df$BRD[i]= colnames(dist_human_BRD)[i]
  a <- dist_human_BRD[,i]
  aa <- as.data.frame(a[a==3])
  ALL_BRD_SP3_df$Neighbors[i] <- toString(rownames(aa))
}

ALL_BRD_SP3_similarity <- data.frame(matrix(ncol = 42, nrow = 42))
rownames(ALL_BRD_SP3_similarity) = colnames(dist_human_BRD); colnames(ALL_BRD_SP3_similarity) = colnames(dist_human_BRD)

for(i in 1:42) {
  for(j in 1:42) {
    aa=unlist(strsplit(ALL_BRD_SP3_df$Neighbors[i], ","))
    bb=unlist(strsplit(ALL_BRD_SP3_df$Neighbors[j], ","))
    cc=intersect(aa,bb)
    dd=union(aa, bb)
    ALL_BRD_SP3_similarity[i,j] = length(cc)/length(dd)
  }
}

### Plot
col.order <- as.character(BRD_level)
row.order <- as.character(BRD_level)
ALL_BRD_SP3_similarity_ordered <- ALL_BRD_SP3_similarity
ALL_BRD_SP3_similarity_ordered <- as.matrix(ALL_BRD_SP3_similarity_ordered)
ALL_BRD_SP3_similarity_ordered <- ALL_BRD_SP3_similarity_ordered[row.order,col.order]
ALL_BRD_SP3_similarity_ordered[upper.tri(ALL_BRD_SP3_similarity_ordered)] <- NA


col <- colorRampPalette(c("ivory2", "lightgoldenrod3", "lightsteelblue","olivedrab1", "springgreen")) #2edit
breaks <- c(0, 0.1, 0.2, 0.5, 0.9, 1)

graphics.off()
gc()
heatmap.2(ALL_BRD_SP3_similarity_ordered, Colv=NULL, Rowv= NULL, col=col, 
          breaks=breaks,
          trace= "none", denscol="black", keysize=1.5, 
          RowSideColors = rep(FuncGroup_summary$color, times = FuncGroup_summary$Num_BRDs), main = "Similarity Heatmap of Bromodomain-containing Proteins(SP=3)", xlab = "BRDs")




hist(ALL_BRD_SP3_similarity_ordered); median(ALL_BRD_SP3_similarity_ordered, na.rm=T)

mean(BRD_simi, na.rm=T); median(BRD_simi, na.rm=T)
mean(ALL_BRD_SP2_similarity_ordered, na.rm=T); median(ALL_BRD_SP2_similarity_ordered, na.rm=T)
mean(ALL_BRD_SP3_similarity_ordered, na.rm=T); median(ALL_BRD_SP3_similarity_ordered, na.rm=T)

### trial
sp1 <- as.vector(BRD_simi)
sp2 <- as.vector(ALL_BRD_SP2_similarity_ordered)
sp3 <- as.vector(ALL_BRD_SP3_similarity_ordered)

cc <- as.data.frame(cbind(sp1, sp2, sp3))
similarity_tidy <- melt(cc)
library(tidyr)
similarity_tidy <- similarity_tidy %>% drop_na()

### plot
library(ggplot2)
similarity_tidy$variable <- as.factor(similarity_tidy$variable)
p<-ggplot(similarity_tidy, aes(x=variable, y=value, fill=variable)) +
  geom_boxplot()+scale_fill_brewer(palette="Dark2")
p+labs(title="Plot of similarities for different shortest paths",x="Shortest Paths", y = "Similarity") + 
  stat_summary(fun="mean", color="red", shape=15)

### save objects 
ALL_BRD_SP2_similarity <- as.matrix(ALL_BRD_SP2_similarity)
ALL_BRD_SP3_similarity <- as.matrix(ALL_BRD_SP3_similarity)
save(ALL_BRD_SP2_similarity, ALL_BRD_SP3_similarity, 
     file = "~/BRD_interaction/BRD_relationship/global_BRD_similarity_SP2&3.Rdata")