rm(list=ls())
load("~/BRD_interaction/Interactions_2021/BRD_Interaction_Data.RData")

table(add_on$Interactor_A)
Interaction_summary <- as.data.frame(table(add_on$Interactor_A))

Interactome_list=c("BRD2", "BRD3", "BRD4", "BRD9")
#BioGrid_BRD_21Feb_short <- BioGrid_BRD_21Feb[,1:2]
Database_collection <- as.data.frame(rbind(BioGrid_BRD_21Feb_short, hippie_BRD_21Feb))

BRDS_collection <- Database_collection[is.element(Database_collection$Interactor_A, Interactome_list) | is.element(Database_collection$Interactor_B, Interactome_list),]

# BRD2
BRD2_Int_addon <- add_on$Interactor_B[add_on$Interactor_A=="BRD2"]
BRD2_Int_collection <- BRDS_collection[BRDS_collection$Interactor_A=="BRD2" | BRDS_collection$Interactor_B=="BRD2",]
BRD2_Int_list <- unique(as.list(rbind(as.character(BRD2_Int_collection$Interactor_A), as.character(BRD2_Int_collection$Interactor_B))))
length(BRD2_Int_list)-1; length(setdiff(BRD2_Int_addon, BRD2_Int_list))
Interaction_summary$Int_collection[1]=230; Interaction_summary$setdiff[[1]]=70

# BRD3
BRD3_Int_addon <- add_on$Interactor_B[add_on$Interactor_A=="BRD3"]
BRD3_Int_collection <- BRDS_collection[BRDS_collection$Interactor_A=="BRD3" | BRDS_collection$Interactor_B=="BRD3",]
BRD3_Int_list <- unique(as.list(rbind(as.character(BRD3_Int_collection$Interactor_A), as.character(BRD3_Int_collection$Interactor_B))))
Interaction_summary$Int_collection[2]=length(BRD3_Int_list)-1
Interaction_summary$setdiff[2]=length(setdiff(BRD3_Int_addon, BRD3_Int_list))
  
# BRD4
BRD4_Int_addon <- add_on$Interactor_B[add_on$Interactor_A=="BRD4"]
BRD4_Int_collection <- BRDS_collection[BRDS_collection$Interactor_A=="BRD4" | BRDS_collection$Interactor_B=="BRD4",]
BRD4_Int_list <- unique(as.list(rbind(as.character(BRD4_Int_collection$Interactor_A), as.character(BRD4_Int_collection$Interactor_B))))
Interaction_summary$Int_collection[3]=length(BRD4_Int_list)-1
Interaction_summary$setdiff[3]=length(setdiff(BRD4_Int_addon, BRD4_Int_list))

# BRD9
BRD9_Int_addon <- add_on$Interactor_B[add_on$Interactor_A=="BRD9"]
BRD9_Int_collection <- BRDS_collection[BRDS_collection$Interactor_A=="BRD9" | BRDS_collection$Interactor_B=="BRD9",]
BRD9_Int_list <- unique(as.list(rbind(as.character(BRD9_Int_collection$Interactor_A), as.character(BRD9_Int_collection$Interactor_B))))
Interaction_summary$Int_collection[4]=length(BRD9_Int_list)-1
Interaction_summary$setdiff[4]=length(setdiff(BRD9_Int_addon, BRD9_Int_list))

#### Plot
Interaction_tidy <- subset(Interaction_summary, select=c(BRDS, Int_collection, setdiff))
Interaction_tidy <- melt(Interaction_tidy)

Interaction_tidy$Category=factor(Interaction_tidy$Category, levels = c("setdiff", "Int_collection"))
ggplot(data=Interaction_tidy, aes(x=BRDS, y=Interactions, fill=Category)) +
  geom_bar(stat="identity")

library(plyr)
Interaction_cumsum <- ddply(Interaction_tidy, "BRDS",
                   transform, label_ypos=cumsum(Interactions))

ggplot(data=Interaction_cumsum, aes(x=BRDS, y=Interactions, fill=Category)) +
  geom_bar(stat="identity")+
  geom_text(aes(y=label_ypos, label=Interactions), vjust=1.6, 
            color="white", size=3.5)+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()