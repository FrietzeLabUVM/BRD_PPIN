#### updated Graph Analysis on the whole Human Interaction Network #####
rm(list=ls())
load("/slipstream/home/conggao/BRD_interaction/Interactions_2021/Human_Physical_Interaction.RData")
library(igraph)
## check
#Human_igraph_unweighted_copy <- igraph::simplify(Human_igraph_unweighted, remove.multiple = T, remove.loops = T, edge.attr.comb=list(weight="sum","ignore")) # passed


# general info
edge_num = gsize(Human_igraph_unweighted) #of edges 559183
node_num = gorder(Human_igraph_unweighted) #of nodes 19843

# Note that edge weights are used by default, unless set to NA.
# K-core decomposition 
Human_unweighted_nodes$Kcore <- coreness(Human_igraph_unweighted)
# Clustering coeffcients
Human_unweighted_nodes$clust_coef<- igraph::transitivity(Human_igraph_unweighted, type = "local", vids = NULL, weights = NA, isolates = "NaN")  
#Centrality Measurements
Human_unweighted_nodes$closeness <- igraph::closeness(Human_igraph_unweighted, vids = V(Human_igraph_unweighted), weights = NA, normalized = FALSE) 
Human_unweighted_nodes$eigen_centrality <- igraph::eigen_centrality(Human_igraph_unweighted, directed=F, weights=NA)$vector
Human_unweighted_nodes$betweenness <- igraph::betweenness(Human_igraph_unweighted, directed=F, weights=NA) 
Human_unweighted_edges$edge_betweenness <- igraph::edge_betweenness(Human_igraph_unweighted, directed=F, weights=NA) 
Human_unweighted_nodes$hub_score <- hub_score(Human_igraph_unweighted, weights=NA)$vector
eccentricity <- eccentricity(Human_igraph_unweighted, vids = V(Human_igraph_unweighted), mode = "all")
Human_unweighted_nodes$eccentricity <- eccentricity
save(Human_igraph_unweighted, Human_unweighted_nodes, Human_unweighted_edges, file="/slipstream/home/conggao/BRD_interaction/Interactions_2021/Human_Physical_Interaction.RData")                                                                                                                                                                                                                          
save(BRD_comp_unique, BRD_comp_edge_SL_unique, BRD_interactor_list, Human_unweighted_nodes, file = "/slipstream/home/conggao/BRD_interaction/Interactions_2021/BRD_Interactions_Analysis_20210902.RData")

####################
graph_clustcoef <- igraph::transitivity(Human_igraph_unweighted, type = "global", vids = NULL, weights = NA, isolates ="NaN")
graph_clustcoef2 <- igraph::transitivity(Human_igraph_unweighted, type = "localaverage", vids = NULL, weights = NA, isolates ="NaN")

edge_density <- edge_density(Human_igraph_unweighted, loops=F)

# need to pay attention if the graph is not connected
diameter <- diameter(Human_igraph_unweighted, directed = F, unconnected = TRUE, weights = NULL)# quite slow
diameter2 <- diameter(Human_igraph_unweighted, directed = F, unconnected = TRUE, weights = NA)# quite slow, correct one

diam <- get_diameter(Human_igraph_unweighted, directed=F); diam <- as.vector(diam) # returns a path with the actual diameter. If there are many shortest paths of the length of the diameter, then it returns the first one found.
radius <- igraph::radius(Human_igraph_unweighted) # quite slow

centr_degree <- centr_degree(Human_igraph_unweighted, normalized=T) #!to understand
centr_clo <- centr_clo(Human_igraph_unweighted, normalized=T)
centr_eigen <- centr_eigen(Human_igraph_unweighted, directed=F, normalized=T) # 
centr_betw <- centr_betw(Human_igraph_unweighted, directed=F, normalized=T) #


# 1 strong component
vertex_connectivit <- vertex_connectivity(Human_igraph_unweighted, source = NULL, target = NULL,
                                          checks = TRUE) #0
is_connected(Human_igraph_unweighted) #false
# component_distribution(Human_igraph_unweighted, cumulative = FALSE, mul.size = FALSE); table(component_distribution(Human_igraph_unweighted, cumulative = FALSE, mul.size = FALSE)) # then plot?
stong_component <- components(Human_igraph_unweighted, mode = "strong")

# distance, whether not-connected
mean_dist <- mean_distance(Human_igraph_unweighted, directed=F)
# dis_matirx = distances(BRD_igraph, weights=NA)# 

save(mean_dist, stong_component, vertex_connectivit, diameter, diam, radius, 
     centr_degree, centr_clo, centr_betw, centr_eigen, edge_density, graph_clustcoef, graph_clustcoef2,
     edge_num, node_num, file = "/slipstream/home/conggao/BRD_interaction/Interactions_2021/Human_Interaction_Graph.RData")

### update on 05/20/2021, analysis on degree
# to run
hist(Human_unweighted_nodes$degree, main="Human Interaction Network", xlab="Degree")
deg.dist <- degree_distribution(Human_igraph_unweighted, cumulative=T, loops = FALSE)
plot( x=0:max(Human_unweighted_nodes$degree), y=1-deg.dist, pch=19, cex=1.2, col="orange", 
      xlab="Degree", ylab="Cumulative Frequency")
gdf = data.frame(Human_unweighted_nodes$degree); colnames(gdf)="Degree"
ggplot(gdf, aes(x=Degree)) + 
  geom_histogram(binwidth=6, color="darkred", fill="red")
ggplot(gdf, aes(Degree)) +
  stat_ecdf(geom = "point", pad = FALSE, color="darkred")+labs(x="Degree", y="Cumulative Frequency")

fit_degree <- fit_power_law(Human_unweighted_nodes$degree+1) 
fit_degree$KS.p  #0.9106701; 
#Small p-values (less than 0.05) indicate that the test rejected the hypothesis 
##that the original data could have been drawn from the fitted power-law distribution.
fit_degree$KS.stat #0.02201176
#The test statistic of a Kolmogorov-Smirnov test that compares the fitted distribution with the input vector. 
##Smaller scores denote better fit.


## Specifically for BRD PPIN
hist(Human_unweighted_nodes$degree[is.element(Human_unweighted_nodes$name, BRD_interactor_list)], main=NULL, xlab="Degree", ylab="Frequency", col = "darkred")
df = data.frame(Human_unweighted_nodes$degree[is.element(Human_unweighted_nodes$name, BRD_interactor_list)]); colnames(df)="Degree"
ggplot(df, aes(x=Degree)) + 
  geom_histogram(binwidth=6, color="darkred", fill="red")
ggplot(df, aes(Degree)) +
  stat_ecdf(geom = "point", pad = FALSE, color="darkred")+labs(x="Degree", y="Cumulative Frequency")

fit_degree <- fit_power_law(Human_unweighted_nodes$degree[is.element(Human_unweighted_nodes$name, BRD_interactor_list)]+1) 
fit_degree$KS.p  #0.7308415 
fit_degree$KS.stat #0.02931668


### update on 08/01/2021, type the edges
Human_unweighted_edges$type = ifelse(is.element(Human_unweighted_edges$from, BRD_list$V1) | is.element(Human_unweighted_edges$to, BRD_list$V1), "B","N") # B:BRD_associated; N: not BRD_associated
table(Human_unweighted_edges$type) #8612 & 550571
Human_PPI <- Human_unweighted_edges[, c(1,2,5)]
write.table(Human_PPI, file = "~/BRD_interaction/Interactions_2021/Human_PPI.txt", row.names = F, col.names = T, 
            sep = "\t", quote = F)

#plot(Human_igraph_unweighted)



