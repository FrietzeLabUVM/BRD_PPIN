getwd()
## [1] "/slipstream/home/conggao/BRD_interaction/Interactions_2021"
non_BRD_interactors <- readRDS("./non-BRD_Interactors.rds")
#load("./BRD_FL.RData") # to have BRD_FL_edges

############# Updated on 06/02/2021
BRD_FL_interaction <- read.table("./BRD_links_FL.txt", sep="\t", stringsAsFactors = FALSE)
BRD_nonBRD_interaction <- BRD_FL_interaction[is.element(BRD_FL_interaction$V1, non_BRD_interactors$Interactor_Name) | 
                                               is.element(BRD_FL_interaction$V2, non_BRD_interactors$Interactor_Name),1:2]#8491
                                   
BRD_nonBRD_interaction$non_BRD <- ifelse(is.element(BRD_nonBRD_interaction$V1, non_BRD_interactors$Interactor_Name), BRD_nonBRD_interaction$V1, BRD_nonBRD_interaction$V2)
BRD_nonBRD_interaction$BRD <- ifelse(is.element(BRD_nonBRD_interaction$V1, non_BRD_interactors$Interactor_Name), BRD_nonBRD_interaction$V2, BRD_nonBRD_interaction$V1)

# BRD_FL_interaction <- BRD_FL_edges[is.element(BRD_FL_edges$from, non_BRD_interactors$Interactor_Name) 
#                                    | is.element(BRD_FL_edges$to, non_BRD_interactors$Interactor_Name),1:2]
# BRD_FL_interaction$non_BRD <- ifelse(is.element(BRD_FL_interaction$from, non_BRD_interactors$Interactor_Name), BRD_FL_interaction$from, BRD_FL_interaction$to)
# BRD_FL_interaction$BRD <- ifelse(is.element(BRD_FL_interaction$from, non_BRD_interactors$Interactor_Name), BRD_FL_interaction$to, BRD_FL_interaction$from)

table(BRD_nonBRD_interaction$BRD)
BRD_list <- read.table("../BRD_list.txt", sep="\t", header = FALSE)
BRD_BET <- read.table("../BRD_classBET.txt", sep="\t", header = TRUE)
save(BRD_list, BRD_BET, file = "./BRD_Type.Rdata")
table(is.element(BRD_nonBRD_interaction$BRD, BRD_list$V1)) # TRUE 8491
table(is.element(BRD_nonBRD_interaction$non_BRD, BRD_list$V1))  #FALSE 8491
table(is.element(non_BRD_interactors$Interactor_Name, BRD_nonBRD_interaction$non_BRD)) #TRUE 4012

BRD_Int_nonBRD <- BRD_nonBRD_interaction[,c(4,3)] #8612
saveRDS(BRD_Int_nonBRD, file = "./BRD_Int_nonBRD_0604.rds")
