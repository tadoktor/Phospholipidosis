

library(collapsibleTree)
#TG Gates Pathologies

TG_Gates_pathologies<- read.csv("file:///C:/Users/tatyana/Documents/Projects/Advance project/open_tggates_pathology.csv",sep=",", header=TRUE)
PC<- read.csv2("file:///C:/Users/tatyana/Documents/projects/KIT_phospholipidosis/Phospholipidosis_positive_list.csv")
NC<- read.csv2("file:///C:/Users/tatyana/Documents/projects/KIT_phospholipidosis/Phospholipidosis_negative_list.csv")



library(dplyr)
PC <- as.data.frame(apply(PC,2,toupper))
NC<- as.data.frame(apply(NC,2,toupper))

TG_Gates_pathologies<-as.data.frame(apply(TG_Gates_pathologies, 2, toupper))
Pathologies_unique_chem<-TG_Gates_pathologies %>% distinct(COMPOUND_NAME)


TG_Gates_biochem_PC<-merge (x= TG_Gates_pathologies, y=PC, by.x="COMPOUND_NAME", by.y= "Compound.name")
TG_Gates_biochem_NC<-merge (x= TG_Gates_pathologies, y=NC, by.x="COMPOUND_NAME", by.y= "Compound.name")



Liver_only_PC<- filter(TG_Gates_biochem_PC, ORGAN == "LIVER")
Liver_only_NC<- filter(TG_Gates_biochem_NC, ORGAN == "LIVER")


Liver_only_unique_PC<-Liver_only_PC %>% distinct(COMPOUND_NAME)
Liver_only_unique_NC<-Liver_only_NC %>% distinct(COMPOUND_NAME)



PC_Grade_type_severe<- filter(Liver_only_PC,GRADE_TYPE == "SEVERE")
PC_Grade_type_moderate<- filter(Liver_only_PC,GRADE_TYPE == "MODERATE")
PC_Grade_type_slight<- filter(Liver_only_PC,GRADE_TYPE == "SLIGHT")
PC_Grade_type_minimal<- filter(Liver_only_PC,GRADE_TYPE == "MINIMAL")



NC_Grade_type_severe<- filter(Liver_only_NC,GRADE_TYPE == "SEVERE")
NC_Grade_type_moderate<- filter(Liver_only_NC,GRADE_TYPE == "MODERATE")
NC_Grade_type_slight<- filter(Liver_only_NC,GRADE_TYPE == "SLIGHT")
NC_Grade_type_minimal<- filter(Liver_only_NC,GRADE_TYPE == "MINIMAL")

#SEVERE

PC_sub_severe<- subset(PC_Grade_type_severe, select=c("SACRIFICE_PERIOD","DOSE_LEVEL" , "FINDING_TYPE", "TOPOGRAPHY_TYPE"))
NC_sub_severe<- subset(NC_Grade_type_severe, select=c("SACRIFICE_PERIOD","DOSE_LEVEL" , "FINDING_TYPE", "TOPOGRAPHY_TYPE"))

PC_sub_severe$new<-paste(PC_sub_severe$TOPOGRAPHY_TYPE,PC_sub_severe$FINDING_TYPE,sep="_")
NC_sub_severe$new<-paste(NC_sub_severe$TOPOGRAPHY_TYPE,NC_sub_severe$FINDING_TYPE,sep="_")

PC_severe_list<-unique(subset(PC_sub_severe, select=c("new")))
NC_severe_list<-unique(subset(NC_sub_severe, select=c("new")))

#MODERATE

PC_sub_moderate<- subset(PC_Grade_type_moderate, select=c("SACRIFICE_PERIOD","DOSE_LEVEL" , "FINDING_TYPE", "TOPOGRAPHY_TYPE"))
NC_sub_moderate<- subset(NC_Grade_type_moderate, select=c("SACRIFICE_PERIOD","DOSE_LEVEL" , "FINDING_TYPE", "TOPOGRAPHY_TYPE"))

PC_sub_moderate$new<-paste(PC_sub_moderate$TOPOGRAPHY_TYPE,PC_sub_moderate$FINDING_TYPE,sep="_")
NC_sub_moderate$new<-paste(NC_sub_moderate$TOPOGRAPHY_TYPE,NC_sub_moderate$FINDING_TYPE,sep="_")

PC_moderate_list<-unique(subset(PC_sub_moderate, select=c("new")))
NC_moderate_list<-unique(subset(NC_sub_moderate, select=c("new")))



#SLIGHT
PC_sub_slight<- subset(PC_Grade_type_slight, select=c("SACRIFICE_PERIOD","DOSE_LEVEL" , "FINDING_TYPE", "TOPOGRAPHY_TYPE"))
NC_sub_slight<- subset(NC_Grade_type_slight, select=c("SACRIFICE_PERIOD","DOSE_LEVEL" , "FINDING_TYPE", "TOPOGRAPHY_TYPE"))

PC_sub_slight$new<-paste(PC_sub_slight$TOPOGRAPHY_TYPE,PC_sub_slight$FINDING_TYPE,sep="_")
NC_sub_slight$new<-paste(NC_sub_slight$TOPOGRAPHY_TYPE,NC_sub_slight$FINDING_TYPE,sep="_")

PC_slight_list<-unique(subset(PC_sub_slight, select=c("new")))
NC_slight_list<-unique(subset(NC_sub_slight, select=c("new")))




#MINIMAL
PC_sub_minimal<- subset(PC_Grade_type_minimal, select=c("SACRIFICE_PERIOD","DOSE_LEVEL" , "FINDING_TYPE", "TOPOGRAPHY_TYPE"))
NC_sub_minimal<- subset(NC_Grade_type_minimal, select=c("SACRIFICE_PERIOD","DOSE_LEVEL" , "FINDING_TYPE", "TOPOGRAPHY_TYPE"))

PC_sub_minimal$new<-paste(PC_sub_minimal$TOPOGRAPHY_TYPE,PC_sub_minimal$FINDING_TYPE,sep="_")
NC_sub_minimal$new<-paste(NC_sub_minimal$TOPOGRAPHY_TYPE,NC_sub_minimal$FINDING_TYPE,sep="_")

PC_minimal_list<-unique(subset(PC_sub_minimal, select=c("new")))
NC_minimal_list<-unique(subset(NC_sub_minimal, select=c("new")))

#Pathologies in common
Severe_common<-intersect(PC_severe_list, NC_severe_list)
Moderate_common<-intersect(PC_moderate_list, NC_moderate_list)
Slight_common<-intersect(PC_slight_list, NC_slight_list)
Minimal_common<-intersect(PC_minimal_list, NC_minimal_list)

PC_severe_specific_1<-setdiff(PC_severe_list, NC_severe_list)
PC_moderate_specific_1<-setdiff(PC_moderate_list, NC_moderate_list)



PC_common_pathologies<-data.frame (unique(rbind(PC_severe_specific_1, PC_moderate_specific_1)))
PC_Pathology<-PC_common_pathologies %>% separate(new, c("a", "b"), extra = "merge", fill = "left")
collapsibleTree(df=PC_Pathology, c("a", "b"),  fill = "lightgreen")

NC_severe_specific_1<-setdiff(NC_severe_list, PC_severe_list)
NC_moderate_specific_1<-setdiff(NC_moderate_list, PC_moderate_list)

NC_common_pathologies<-data.frame (unique(rbind(NC_severe_specific_1, NC_moderate_specific_1)))
library(tidyr)

NC_Pathology<-NC_common_pathologies %>% separate(new, c("a", "b"), extra = "merge", fill = "left")

collapsibleTree(df=NC_Pathology, c("a", "b"),  fill = "lightgreen")



