#### This script is for the newest BioGrid data collected on Feb 2, 2021 ####
#### BioGrid download tab3 version:4.3.194, modified on January 28th, 2021 ####

### Process ALL BioGrid data#758468 -> 1985662
biogrid_all <- read.table("~/BRD_interaction/Interactions_2021/BIOGRID-ALL-4.3.194.tab3.txt", quote="", sep="\t", fill = TRUE, header = TRUE, stringsAsFactors = FALSE)
# extract column 8,9,13,18,36,37 & rename columns(IntA refers to Interactor_A)
BioGrid_All_21Feb <- biogrid_all[,c(8,9,13,18,36,37)]
colnames(BioGrid_All_21Feb) = c("Interactor_A", "Interactor_B", "Experimental_type", "Throughput", "Organism_IntA", "Organism_IntB")
BioGrid_All21Feb_HomoPhysical <- subset(BioGrid_All_21Feb, Organism_IntA=="Homo sapiens" & Organism_IntB=="Homo sapiens" & Experimental_type=="physical")
# check! #634922
table(BioGrid_All21Feb_HomoPhysical$Experimental_type); table(BioGrid_All21Feb_HomoPhysical$Organism_IntA); table(BioGrid_All21Feb_HomoPhysical$Organism_IntB)
#!! Using this generate BRD interaction data to see !! since the file are not read in completely before
BRD_list <- read.table("~/Useful/BRD_list.txt", sep="\t", header = F, stringsAsFactors = FALSE)
BioGrid_BRD_21Feb <- BioGrid_All21Feb_HomoPhysical[is.element(BioGrid_All21Feb_HomoPhysical$Interactor_A, BRD_list$V1) | is.element(BioGrid_All21Feb_HomoPhysical$Interactor_B, BRD_list$V1),]
#check ##1
listA <- BRD_list$V1[is.element(BRD_list$V1, BioGrid_BRD_21Feb$Interactor_A)==FALSE]
listB <- BRD_list$V1[is.element(BRD_list$V1, BioGrid_BRD_21Feb$Interactor_B)==FALSE]

# ### Process BRD interaction data; 11407 interactions for human physical
# BRD_all <- read.table("~/BRD_interaction/Interactions_2021/new.txt", sep="\t", quote="", fill=TRUE, header = TRUE, stringsAsFactors = FALSE)
# BRD_All_21Feb <- BRD_all[,c(8,9,13,18,36,37)]
# colnames(BRD_All_21Feb) = c("Interactor_A", "Interactor_B", "Experimental_type", "Throughput", "Organism_IntA", "Organism_IntB")
# BRD_All21Feb_HomoPhysical <- subset(BRD_All_21Feb, Organism_IntA=="Homo sapiens" & Organism_IntB=="Homo sapiens" & Experimental_type=="physical")
# # check, listA and listB not overlap
# listA <- BRD_list$V1[is.element(BRD_list$V1, BRD_All21Feb_HomoPhysical$Interactor_A)==FALSE]
# listB <- BRD_list$V1[is.element(BRD_list$V1, BRD_All21Feb_HomoPhysical$Interactor_B)==FALSE]

### Process add_on data; 502 interactions, note quote=""
add_on <- read.table("~/BRD_interaction/add_interactions.txt", sep="\t", quote="", header = TRUE, stringsAsFactors = FALSE)
add_on <- add_on[,1:2]
colnames(add_on) <- c("Interactor_A", "Interactor_B")
add_on$Experimental_type ="physical"
add_on$Throughput="High Throughput add_on"
add_on$Organism_IntA="Homo sapiens"
add_on$Organism_IntB="Homo sapiens"
### to edit
add_on <- subset(add_on, Interactor_B != "LOC100996657" & Interactor_B != "ATP5L2")
add_on$Interactor_B[add_on$Interactor_B == "C20orf20"] = "MRGBP"
add_on$Interactor_B[add_on$Interactor_B == "GTF2H2D"] = "GTF2H2C_2"
add_on$Interactor_B[add_on$Interactor_B == "MKI67IP"] = "NIFK"
add_on$Interactor_B[add_on$Interactor_B == "COBRA1"] = "NELFB"
add_on$Interactor_B[add_on$Interactor_B == "TH1L"] = "NELFCD"
add_on$Interactor_B[add_on$Interactor_B == "WHSC2"] = "NELFA" 

### New database, HIPPIE and STRING, human and physical
## a1. HIPPIE #391410
hippie_all <- read.table("~/BRD_interaction/Interactions_2021/HIPPIE-current.mitab.txt", sep="\t", quote="", fill=TRUE, header = TRUE, stringsAsFactors = FALSE)
table(hippie_all$Taxid.Interactor.A); table(hippie_all$Taxid.Interactor.B)
hippie_All_21Feb <- hippie_all[,17:18]
colnames(hippie_All_21Feb) <- c("Interactor_A", "Interactor_B")
hippie_BRD_21Feb <- hippie_All_21Feb[is.element(hippie_All_21Feb$Interactor_A, BRD_list$V1) | is.element(hippie_All_21Feb$Interactor_B, BRD_list$V1),]
listAA <- BRD_list$V1[is.element(BRD_list$V1, hippie_BRD_21Feb$Interactor_A)==FALSE]
listBB <- BRD_list$V1[is.element(BRD_list$V1, hippie_BRD_21Feb$Interactor_B)==FALSE]

## SAVE BRD interactions
save(add_on, BioGrid_BRD_21Feb, hippie_BRD_21Feb, file = "~/BRD_interaction/Interactions_2021/BRD_Interaction_Data.RData")

# ## a2. STRING, not use STRING
# STRING_all <- read.table("~/BRD_interaction/Interactions_2021/9606.protein.physical.links.detailed.v11.0.txt", sep="", dec = ".", header = TRUE, stringsAsFactors = FALSE)


### Combine BRD interaction with add_on data --> First level BRD interaction data, total 16000
comp_BRD_interaction_pre <- rbind(BioGrid_BRD_21Feb, add_on)
comp_BRD_interaction <- rbind(comp_BRD_interaction_pre[,1:2], hippie_BRD_21Feb) # BRD-associated physical interactions in human 16001-->16000
comp_HomoPhysical <- rbind(BioGrid_All21Feb_HomoPhysical[,1:2], hippie_All_21Feb) # All physical interactions in human, 1026332
#check
table(nzchar(comp_BRD_interaction$Interactor_B)); table(nzchar(comp_BRD_interaction$Interactor_A)) # seems empty string exists in comp_BRD_interaction
comp_BRD_interaction[nzchar(comp_BRD_interaction$Interactor_B)==FALSE,]
comp_BRD_interaction <- comp_BRD_interaction[row.names(comp_BRD_interaction) != "60982", , drop = FALSE]

### Combine whole BRD interactions with BioGrid&HIPPIE data(interactions between non_BRDs) --> Second level(SL) BRD interaction data
whole_interactor_list <- unique(as.list(rbind(as.character(comp_BRD_interaction$Interactor_A), as.character(comp_BRD_interaction$Interactor_B))))
table(nzchar(whole_interactor_list)) #4054
comp_BRD_interaction_SL <- comp_HomoPhysical[is.element(comp_HomoPhysical$Interactor_A, whole_interactor_list) & is.element(comp_HomoPhysical$Interactor_B, whole_interactor_list),] #382283
#check
table(nzchar(comp_BRD_interaction_SL$Interactor_B)); table(nzchar(comp_BRD_interaction_SL$Interactor_A)) #381844

### Save raw data
save(comp_HomoPhysical, comp_BRD_interaction, whole_interactor_list, file="~/BRD_interaction/Interactions_2021/HomoPhysical_Interactions.Rdata")

### Construct objects ###
## igraph
library(dplyr)
BRD_comp_edge <- comp_BRD_interaction_SL %>% mutate(type="Human_Physical", weight=1)
colnames(BRD_comp_edge)[1:2] <- c("from","to")
BRD_comp_nodes <- as.list(whole_interactor_list)
library(data.table)
BRD_comp_node=as.data.frame(t(setDT(BRD_comp_nodes, keep.rownames=FALSE))); rownames(BRD_comp_node) = NULL
colnames(BRD_comp_node) = "id"
# check
table(nzchar(BRD_comp_edge$from)); table(nzchar(BRD_comp_edge$to))

#detach("package:graph", unload=TRUE)
#detach("package:BioNet", unload=TRUE)
library(igraph)
BRD_igraph_all <- igraph::graph_from_data_frame(d=BRD_comp_edge, vertices=BRD_comp_node, directed=F) 
class(BRD_igraph_all)
BRD_igraph_unweighted <- igraph::simplify(BRD_igraph_all, remove.multiple = T, remove.loops = T, edge.attr.comb=list(weight="sum","ignore"))
BRD_unweighted_edges <- igraph::as_data_frame(BRD_igraph_unweighted, what="edges")
BRD_unweighted_nodes <- igraph::as_data_frame(BRD_igraph_unweighted, what="vertices")
#label nodes, get non-BRDs with degree 1
BRD_unweighted_nodes$degree <- igraph::degree(BRD_igraph_unweighted, loops = FALSE, normalized = FALSE) #Error in degree(BRD_igraph_unweighted, loops = FALSE, normalized = FALSE) : 
#unused arguments (loops = FALSE, normalized = FALSE), need to re-run!
BRD_class <- read.table("~/BRD_interaction/BRD_classBET.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
BRD_unweighted_nodes$type = "non-BRDs"
BRD_unweighted_nodes$type[is.element(BRD_unweighted_nodes$name, BRD_class$non_BET)] = "non_BET"
BRD_unweighted_nodes$type[is.element(BRD_unweighted_nodes$name, BRD_class$BET)] = "BET"#...
table(BRD_unweighted_nodes$type)
BRD_unweighted_nodes$BRD = 0; BRD_unweighted_nodes$BRD[is.element(BRD_unweighted_nodes$name, BRD_list$V1)] =1 # for stat analysis

# remove nodes with 2-degree, to generate smaller files for clique
BRD_symp_rm <- as.list(BRD_unweighted_nodes$name[BRD_unweighted_nodes$degree<=1])
BRD_comp_edge_smaller <- BRD_comp_edge[is.element(BRD_comp_edge$from, BRD_symp_rm)==FALSE & is.element(BRD_comp_edge$to, BRD_symp_rm)==FALSE,]
BRD_igraph_smaller <- igraph::graph_from_data_frame(d=BRD_comp_edge_smaller, vertices=NULL, directed=F) 
class(BRD_igraph_smaller)
BRD_igraph_unweighted_smaller <- igraph::simplify(BRD_igraph_smaller, remove.multiple = T, remove.loops = T, edge.attr.comb=list(weight="sum","ignore"))
length(V(BRD_igraph_unweighted_smaller)); length(E(BRD_igraph_unweighted_smaller))
#BRD_symp_rm1 <- as.list(BRD_unweighted_nodes$name[BRD_unweighted_nodes$degree==1 & BRD_unweighted_nodes$type=="non-BRDs"]) #not run yet, stop here.
#BRD_symp_rm2 <- as.list(BRD_unweighted_nodes$name[BRD_unweighted_nodes$degree==2 & BRD_unweighted_nodes$type=="non-BRDs"]) #not run yet, stop here.
##igraph::max_cliques(BRD_igraph_unweighted_smaller, file = "~/BRD_interaction/Interactions_2021/All_MaxCliques.rds") #2run
##igraph::max_cliques(BRD_igraph_unweighted_smaller, subset = V(BRD_igraph_unweighted_smaller)[is.element(V(BRD_igraph_unweighted_smaller)$name, BRD_list$V1)], file = "~/BRD_interaction/Interactions_2021/BRD_MaxCliques.rds") #2try

### Save igraph objects, node, edge and symplified object
save(BRD_igraph_unweighted_smaller, BRD_igraph_unweighted, file="~/BRD_interaction/Interactions_2021/BRD_igraphSL_objects.Rdata")
save(BRD_unweighted_edges, BRD_unweighted_nodes, file="~/BRD_interaction/Interactions_2021/BRD_interaction_analysis.Rdata")
save(BRD_comp_edge, BRD_comp_node, file="~/BRD_interaction/Interactions_2021/igraph_prepare.Rdata")

## BioNet
library(BioNet)
BRD_graphNEL_pre <- BRD_igraph_unweighted;  E(BRD_graphNEL_pre)$weight = 1
BRD_graphNEL <-as_graphnel(BRD_graphNEL_pre)
### save BRD_graphNEL object
saveRDS(BRD_graphNEL, "~/BRD_interaction/Interactions_2021/BRD_graphNEL.rds") 

#### The end of this processing script
rm(list = ls())
