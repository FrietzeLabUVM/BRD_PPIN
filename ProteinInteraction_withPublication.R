# Summarize publication for each BRD

biogrid_all <- read.table("~/BRD_interaction/Interactions_2021/BIOGRID-ALL-4.3.194.tab3.txt", quote="", sep="\t", fill = TRUE, header = TRUE, stringsAsFactors = FALSE)
BioGrid_All <- biogrid_all[,c(8,9,13,18,36,37,15)]
colnames(BioGrid_All) = c("Interactor_A", "Interactor_B", "Experimental_type", "Throughput", "Organism_IntA", "Organism_IntB", "Publication")
BioGrid_All_HomoPhysical <- subset(BioGrid_All, Organism_IntA=="Homo sapiens" & Organism_IntB=="Homo sapiens" & Experimental_type=="physical")

#
hippie_all <- read.table("~/BRD_interaction/Interactions_2021/HIPPIE-current.mitab.txt", sep="\t", quote="", fill=TRUE, header = TRUE, stringsAsFactors = FALSE)
library(data.table)
dt = fread("~/BRD_interaction/Interactions_2021/HIPPIE-current.mitab.txt", sep="\t", dec=".", stringsAsFactors = FALSE)
dt = dt[, c(9,17,18)]
colnames(dt) = c("Publication", "Interactor_A", "Interactor_B")
dt_pubs = dt[, .(pub_ids = tstrsplit(Publication, "\\|")) , .( Interactor_A,  Interactor_B )]
dt_pubs$Publication <- toupper(dt_pubs$pub_ids)
#dt_merged = merge(dt, dt_pubs, by = c("Interactor_A",  "Interactor_B"))
# dt_test = fread("/slipstream/home/conggao/BRD_interaction/Interactions_2021/HIPPIE_head10.txt")
# dt_pubs = dt_test[, .(pub_ids = tstrsplit(Publication.Identifiers, "\\|")) , .( ID.Interactor.A,  ID.Interactor.B )]
# dt_merged = merge(dt_test, dt_pubs, by = c("ID.Interactor.A",  "ID.Interactor.B"))
Interaction_pub <- rbind(dt_pubs[,c(1,2,4)], BioGrid_All_HomoPhysical[,c(1,2,7)])
add_on$Publication = "PUBMED:30554943"
All_Interaction_Pub <- rbind(Interaction_pub, add_on[,c(1,2,7)])
BRD_Interaction_Pub <- All_Interaction_Pub[is.element(All_Interaction_Pub$Interactor_A, BRD_list$V1) | 
                                             is.element(All_Interaction_Pub$Interactor_B, BRD_list$V1),] #18450
BRD_Interaction_Pub_unique <- unique(BRD_Interaction_Pub) #13731

save(All_Interaction_Pub, BRD_Interaction_Pub, BRD_Interaction_Pub_unique, file="./Interaction_Publication_raw.RData") # 2run!!

colnames(BRD_Interaction_Pub_unique) = c("from", "to", "publication")
BRD_Pub_interactor_list <- unique(as.list(rbind(as.character(BRD_Interaction_Pub_unique$from), 
                                                as.character(BRD_Interaction_Pub_unique$to)))) 
### Backup if unique function doesn't work well, but should be fine
BRD_Pub_List2 <- unique(as.list(rbind(as.character(BRD_Interaction_Pub$Interactor_A), 
                                                as.character(BRD_Interaction_Pub$Interactor_B)))) 


BRD_Pub_nodes <- as.list(BRD_Pub_interactor_list)
BRD_Pub_node=as.data.frame(t(setDT(BRD_Pub_nodes, keep.rownames=FALSE))); rownames(BRD_Pub_node) = NULL
colnames(BRD_Pub_node) = "id"

#Complete BRD igraph
BRD_Pub_igraph <- igraph::graph_from_data_frame(d=BRD_Interaction_Pub_unique, vertices=BRD_Pub_node, directed=F) 
### Backup
BRD_Pub_igraph_backup <- igraph::graph_from_data_frame(d=BRD_Interaction_Pub, vertices=BRD_Pub_node, directed=F) 

class(BRD_Pub_igraph)
E(BRD_Pub_igraph)$weight=1
#E(BRD_Pub_igraph_backup)$weight=1

BRD_Pub_edges <- igraph::as_data_frame(BRD_Pub_igraph, what="edges") 
BRD_Pub_nodes <- igraph::as_data_frame(BRD_Pub_igraph, what="vertices") 

### remove.multiple=T will remove all of the duplicated links regardless of publication
BRD_Pub_unweighted <- igraph::simplify(BRD_Pub_igraph, remove.multiple = F, remove.loops = T)
BRD_Pub_uw_edges <- igraph::as_data_frame(BRD_Pub_unweighted, what="edges") 
BRD_Pub_uw_nodes <- igraph::as_data_frame(BRD_Pub_unweighted, what="vertices") 
BRD_Pub_uw_nodes$Degree <- degree(BRD_Pub_unweighted, mode="all") # Degree is the number of publications for each node
#BRD_Ref <- BRD_Pub_uw_nodes[is.element(BRD_Pub_uw_nodes$name, BRD_list$V1),] # it is not correct

save(BRD_Pub_unweighted, BRD_Pub_uw_edges, BRD_Pub_uw_nodes, 
     file = "./Interaction_Publication_igraph.RData")
#saveRDS(BRD_Ref, file = "./BRD_Publication_Stats.Rds") #update later

# update on 05/13/2021
rm(list = ls())
load("~/BRD_interaction/Interactions_2021/Interaction_Publication_igraph.RData")
table(BRD_Pub_uw_edges$weight) #1*13650
BRD_Ref = readRDS("~/BRD_interaction/Interactions_2021/BRD_Publication_Stats.Rds")
BRD_Ref$Pub=0

for(i in 1:42) {
  interaction_df = BRD_Pub_uw_edges[BRD_Pub_uw_edges$from==BRD_Ref$name[i] | 
                                      BRD_Pub_uw_edges$to==BRD_Ref$name[i], 1:3]
  BRD_Ref$Pub[i] <- length(unique(interaction_df$publication))
}
saveRDS(BRD_Ref, file = "~/BRD_interaction/Interactions_2021/BRD_Publication_Stats.Rds")

BRD_pub_degree <- merge(BRD_Publication_Stats, Human_unweighted_nodes[, c(1,4,6)], by="name", all.x=T)
reg1 <- lm(BRD_pub_degree$Pub~BRD_pub_degree$degree)
plot(BRD_pub_degree$degree, BRD_pub_degree$Pub, pch=19); abline(reg1, col="red", lwd=2); abline(a=0,b=1, col="darkblue", lwd=2) 
identify(BRD_pub_degree$degree, BRD_pub_degree$Pub, labels=BRD_pub_degree$name) 

### need to swap x and y
reg11 <- lm(BRD_pub_degree$degree~BRD_pub_degree$Pub)
summary(reg11)
cor.test(x=BRD_pub_degree$Pub, y=BRD_pub_degree$degree, method = c("pearson")) #0.4955036, p-value = 0.0008482
plot(BRD_pub_degree$Pub, BRD_pub_degree$degree, pch=19); abline(reg11, col="red", lwd=2); abline(a=0,b=1, col="darkblue", lwd=2) 
identify(BRD_pub_degree$Pub, BRD_pub_degree$degree, labels=BRD_pub_degree$name)
table(I(BRD_pub_degree$Pub >56))

### 2021/10/17
reg2 <- lm(BRD_pub_degree$clust_coef*100~BRD_pub_degree$Pub)
plot(BRD_pub_degree$Pub, BRD_pub_degree$clust_coef*100, pch=19); abline(reg2, col="red", lwd=2); abline(a=0,b=1, col="darkblue", lwd=2) 
identify(BRD_pub_degree$Pub, BRD_pub_degree$clust_coef*100, labels=BRD_pub_degree$name) 

plot(BRD_pub_degree$degree, BRD_pub_degree$clust_coef*100, pch=19)

### test
interaction_dff = BRD_Pub_uw_edges[BRD_Pub_uw_edges$from=="BRD4" | 
                                    BRD_Pub_uw_edges$to=="BRD4", 1:3]
npub <- length(unique(interaction_dff$publication))
