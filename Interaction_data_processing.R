#### BRD Interaction Processing ####
### fetch the BioGRID dataset; download file and save to working directory ###
getwd(); setwd("~/BRD_interaction")
url <- "https://downloads.thebiogrid.org/File/BioGRID/Release-Archive/BIOGRID-4.2.191/BIOGRID-ALL-4.2.191.tab3.zip"
destfile <- "BioGRID_all_2020Oct.tab.zip"
download.file(url, destfile, method = "wget")
unzip(zipfile = "./BioGRID_all_2020Oct.tab.zip", exdir = "./BioGRID_all_2020Oct") #!not working
### Retry with all organism data downloaded to lacal and transfer the one for human to the server ###
### url: https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-4.2.191/BIOGRID-ORGANISM-4.2.191.tab3.zip ###

### read and processing the tab file ###
human_all <- read.table("./BIOGRID-ORGANISM-Homo_sapiens-4.2.191.tab3.txt", sep = "\t", fill = TRUE) # since a # at the header, so no header here, need to refer to the local file
table(human_all$V13)
table(human_all$V36)
table(human_all$V37)
Interaction_h2h <- human_all[human_all$V36=="Homo sapiens" & human_all$V37=="Homo sapiens",] #extract the sub-data where interactor A and B are both for human
table(Interaction_h2h$V13)
table(Interaction_h2h$V36)
table(Interaction_h2h$V37)
h2h_all <- subset(Interaction_h2h, select=c(V8, V9, V13, V36, V37)) #extract information data includes the names of Interactors, expreiment type and the organisms of interactors(for further check)

library(data.table) #change column names
setnames(h2h_all, "V8", "Interactor_A")
setnames(h2h_all, "V9", "Interactor_B")
setnames(h2h_all, "V13", "Experiment_type")
setnames(h2h_all, "V36", "InterA_organism")
setnames(h2h_all, "V37", "InterB_organism")

interactions_addon <- read.table("./add_interactions.txt", sep = "\t", header = TRUE)
setnames(interactions_addon, "Bait", "Interactor_A")
setnames(interactions_addon, "PreyGene", "Interactor_B")
interactions_addon$BFDR = "physical"
setnames(interactions_addon, "BFDR", "Experiment_type")
interactions_addon <- subset(interactions_addon, Interactor_B != "LOC100996657" & Interactor_B != "ATP5L2")
interactions_addon$Interactor_B <- unfactor(interactions_addon$Interactor_B)
interactions_addon$Interactor_B[interactions_addon$Interactor_B == "C20orf20"] = "MRGBP"

interactions_addon$Interactor_B[interactions_addon$Interactor_B == "GTF2H2D"] = "GTF2H2C_2"
interactions_addon$Interactor_B[interactions_addon$Interactor_B == "MKI67IP"] = "NIFK"
interactions_addon$Interactor_B[interactions_addon$Interactor_B == "COBRA1"] = "NELFB"
interactions_addon$Interactor_B[interactions_addon$Interactor_B == "TH1L"] = "NELFCD"
interactions_addon$Interactor_B[interactions_addon$Interactor_B == "WHSC2"] = "NELFA"

### Extract BRD interactions ###
BRD_list <- read.table("~/Useful/BRD_list.txt", stringsAsFactors = FALSE) #BRD_list$V1
#typeof(BRD_list$V1)
#[1] "integer"
#typeof(h2h_all$Interactor_A)
#[1] "integer"   #same data type, good to go.
BRD_all <- h2h_all[is.element(h2h_all$Interactor_A, BRD_list$V1) | is.element(h2h_all$Interactor_B, BRD_list$V1),]
#TRIM66_test <- BRD_all[BRD_all$Interactor_A=="TRIM66" | BRD_all$Interactor_B=="TRIM66",]
BRD_physical <- BRD_all[BRD_all$Experiment_type=="physical", 1:3]

#BRD_genetic <- subset(BRD_all, Experiment_type=="genetic", select=c(Interactor_A, Interactor_B, Experiment_type))
#write.table(BRD_physical, "./BRD_physical.tab", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

BRD_interactions <- rbind(BRD_physical, interactions_addon) # 7519, only interactions including BRDs
whole_interactor_list <- unique(as.list(rbind(as.character(BRD_interactions$Interactor_A), as.character(BRD_interactions$Interactor_B))))
Update_BRD_All <- h2h_all[is.element(h2h_all$Interactor_A, whole_interactor_list) & is.element(h2h_all$Interactor_B, whole_interactor_list),]
Update_BRD_physical <- Update_BRD_All[Update_BRD_All$Experiment_type=="physical", 1:3] #108471 interactions, BRD and between non-BRDs(ones interacting with BRDs)
#write.table(BRD_interactions, "./BRD_only_interactions.sif", quote = FALSE, sep = "\t")
#write.table(Update_BRD_physical, "./BRD_all_interactions.sif", quote = FALSE, sep = "\t")


SL_interaction <- h2h_all[is.element(h2h_all$Interactor_A, SL_list) | is.element(h2h_all$Interactor_B, SL_list),]
interactor_list <- as.list(rbind(as.character(SL_interaction$Interactor_A), as.character(SL_interaction$Interactor_B)))
interactor_list <- unique(interactor_list)
SL_interaction_all <- h2h_all[is.element(h2h_all$Interactor_A, interactor_list) & is.element(h2h_all$Interactor_B, interactor_list),]
write.table(SL_interaction_all, "./SL_interaction_all.sif", sep = "\t", row.names = FALSE, quote = FALSE)
### construct input file for igraph ###
#install.packages("igraph")


