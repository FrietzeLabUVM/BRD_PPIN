load("~/BRD_interaction/Interactions_2021/BRD_Interaction_Data.RData")
load("~/BRD_interaction/Interactions_2021/HomoPhysical_Interactions.Rdata")

Human_Complete_Interaction_all <- rbind(comp_HomoPhysical, add_on[,1:2])

#check
table(nzchar(Human_Complete_Interaction_all$Interactor_A)); table(nzchar(Human_Complete_Interaction_all$Interactor_B)) 
Human_Complete_Interaction <- Human_Complete_Interaction_all[nzchar(Human_Complete_Interaction_all$Interactor_A)==TRUE & 
                             nzchar(Human_Complete_Interaction_all$Interactor_B)==TRUE,]
table(is.na(Human_Complete_Interaction$Interactor_A))
table(is.na(Human_Complete_Interaction$Interactor_B))
library(tidyr)
Human_Complete_Interaction <- Human_Complete_Interaction %>% drop_na() #1026053

Human_interactor_list <- unique(as.list(rbind(as.character(Human_Complete_Interaction$Interactor_A), 
                                              as.character(Human_Complete_Interaction$Interactor_B)))) #19843

### Construct objects ###
## igraph
library(dplyr)
Human_comp_edge <- Human_Complete_Interaction %>% mutate(type="Human_Physical", weight=1)
colnames(Human_comp_edge)[1:2] <- c("from","to")
Human_comp_nodes <- as.list(Human_interactor_list)
library(data.table)
Human_comp_node=as.data.frame(t(setDT(Human_comp_nodes, keep.rownames=FALSE))); rownames(Human_comp_node) = NULL
colnames(Human_comp_node) = "id"
# check
table(nzchar(Human_comp_edge$from)); table(nzchar(Human_comp_edge$to))
Human_comp_node$id <- as.character(Human_comp_node$id)
table(nzchar(Human_comp_node$id))

Human_comp_node_copy <- Human_comp_node %>% drop_na() ###19843
Human_comp_edge_copy <- Human_comp_edge %>% drop_na() ###1026053
save(Human_Complete_Interaction, Human_interactor_list, Human_comp_node_copy, Human_comp_edge_copy, 
     file= "/slipstream/home/conggao/BRD_interaction/Interactions_2021/Human_physical_preparation.RData")

library(igraph)
Human_igraph_all <- igraph::graph_from_data_frame(d=Human_comp_edge_copy, vertices=Human_comp_node_copy, directed=F) 
class(Human_igraph_all)
Human_igraph_unweighted <- igraph::simplify(Human_igraph_all, remove.multiple = T, remove.loops = T, edge.attr.comb=list(weight="sum","ignore"))
Human_unweighted_edges <- igraph::as_data_frame(Human_igraph_unweighted, what="edges") ##559183
Human_unweighted_nodes <- igraph::as_data_frame(Human_igraph_unweighted, what="vertices") ##19843


BRD_class <- read.table("~/BRD_interaction/BRD_classBET.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
Human_unweighted_nodes$type = "non-BRDs"
Human_unweighted_nodes$type[is.element(Human_unweighted_nodes$name, BRD_class$non_BET)] = "non_BET"
Human_unweighted_nodes$type[is.element(Human_unweighted_nodes$name, BRD_class$BET)] = "BET"#...
table(Human_unweighted_nodes$type)
Human_unweighted_nodes$BRD = 0; Human_unweighted_nodes$BRD[Human_unweighted_nodes$type != "non-BRDs"] =1 # for stat analysis

save(Human_igraph_unweighted, Human_unweighted_nodes, Human_unweighted_edges, 
     file="/slipstream/home/conggao/BRD_interaction/Interactions_2021/Human_Physical_Interaction.RData")

#### Check whether the degree is correct
BRD4_links <- Human_unweighted_edges[Human_unweighted_edges$from=="BRD4" | Human_unweighted_edges$to=="BRD4",]
BRD4_interactor_list <- unique(as.list(rbind(as.character(BRD4_links$from), 
                                              as.character(BRD4_links$to)))) #1400, one is BRD4, correct!!



### Graph Analysis
Human_unweighted_nodes$degree <- igraph::degree(Human_igraph_unweighted, loops = FALSE, normalized = FALSE) 
# to run
hist(Human_unweighted_nodes$degree, breaks=1:3000, main="Histogram of node degree", xlab="Degree")
deg.dist <- degree_distribution(Human_igraph_unweighted, cumulative=T, loops = FALSE)
plot( x=0:max(Human_unweighted_nodes$degree), y=1-deg.dist, pch=19, cex=1.2, col="orange", 
      xlab="Degree", ylab="Cumulative Frequency")


# when finish this and BRD-BRD analysis, need to update the Human_Physical_Interaction.RData

### Get BRD ####
BRD_list <- read.table("~/BRD_interaction/BRD_list.txt", stringsAsFactors = FALSE, quote = "\t", header = FALSE)
BRD_comp_edge <- Human_unweighted_edges[is.element(Human_unweighted_edges$from, BRD_list$V1) | is.element(Human_unweighted_edges$to, BRD_list$V1),] #8612
BRD_interactor_list <- unique(as.list(rbind(as.character(BRD_comp_edge$from), as.character(BRD_comp_edge$to)))) #4054
BRD_comp_edge_SL <- Human_unweighted_edges[is.element(Human_unweighted_edges$from, BRD_interactor_list) & 
                                             is.element(Human_unweighted_edges$to, BRD_interactor_list),] #192785
BRD_comp_edge_SL_unique <- unique(BRD_comp_edge_SL[,1:2]) #192785

BRD_comp_nodes <- as.list(BRD_interactor_list)
library(data.table)
BRD_comp_node=as.data.frame(t(setDT(BRD_comp_nodes, keep.rownames=FALSE))); rownames(BRD_comp_node) = NULL
colnames(BRD_comp_node) = "id"

#Complete BRD igraph
BRD_igraph_all <- igraph::graph_from_data_frame(d=BRD_comp_edge, vertices=BRD_comp_node, directed=F) 
class(BRD_igraph_all)
BRD_igraph_unweighted <- igraph::simplify(BRD_igraph_all, remove.multiple = T, remove.loops = T, edge.attr.comb=list(weight="sum","ignore"))
BRD_igraph_complete <- igraph::graph_from_data_frame(d=BRD_comp_edge_SL, vertices=BRD_comp_node, directed=F) 
BRD_comp_unweighted <- igraph::simplify(BRD_igraph_complete, remove.multiple = T, remove.loops = T, edge.attr.comb=list(weight="sum","ignore"))


save(BRD_comp_edge, BRD_interactor_list, BRD_comp_edge_SL_unique, BRD_igraph_unweighted, BRD_igraph_complete, BRD_comp_unweighted,
     file="/slipstream/home/conggao/BRD_interaction/Interactions_2021/BRD_interaction_20210509.RData")



