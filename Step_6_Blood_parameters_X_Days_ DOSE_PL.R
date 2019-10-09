

Root_location <- "C:\\Users\\tatyana\\Documents\\Projects\\Advance project"

#if (!require(ggplot2)) {
   # install.packages("ggplot2")}
    library(ggplot2)


#if (!require(plyr)) {
    #install.packages("plyr")}
    library(plyr)


#if (!require(dplyr)) {
   # install.packages("dplyr")
    library(dplyr)


#if (!require(tidyr)) {install.packages("tidyr")}
    library(tidyr)


#if (!require(tidyverse)) {install.packages("tidyverse")}
    library(tidyverse)


#if (!require(reshape2)) { install.packages("reshape2")}
    library(reshape2)


analyze_dataset <- function(df, organ, exp_time) {
    # function takes as argument data frame df and creates a summary data frame for a given organ and exposure time
    
    # select a subset for a specific organ and exposure time
    df_sub <- subset(df, ( (organ_id == as.character(organ)) & (exposure_time == as.character(exp_time)) ))

    # add a column for unique sample names
    df_sub$sample_name <- paste(df_sub$compound_name, df_sub$exposure_time, df_sub$dose, df_sub$dose_level)

    # ensure that numeric columns are numeric
    ilow = 22  # lowest index of numeric columns
    ihigh = 66  # highest index of numeric columns
    for (i in c(ilow : ihigh)) {
        df_sub[i] <- as.numeric(as.character(df_sub[, i]))
    }

    
    # initiate summary data frame for unique sample names
    summary <- data.frame("sample_name" = unique(df_sub$sample_name))

    # add columns "compound_name", "dose_level", "dose", "exposure_time", "organ_id"
    for (smp in unique(df_sub$sample_name)) {
        summary[summary$sample_name == smp, "compound_name"] <- (df_sub[df_sub$sample_name == smp, "compound_name"])[1]
        summary[summary$sample_name == smp, "dose_level"] <- (df_sub[df_sub$sample_name == smp, "dose_level"])[1]
        summary[summary$sample_name == smp, "dose"] <- (df_sub[df_sub$sample_name == smp, "dose"])[1]
        summary[summary$sample_name == smp, "exposure_time"] <- (df_sub[df_sub$sample_name == smp, "exposure_time"])[1]
        summary[summary$sample_name == smp, "organ_id"] <- (df_sub[df_sub$sample_name == smp, "organ_id"])[1]
    }

    
    
  
    
    
    # calculation of mean, median and sd
    # loop through numeric columns
    for (i in c(ilow : ihigh)) {
        # loop through unique sample names
        for (smp in unique(df_sub$sample_name)) {
            # select replicates with the same sample name
            x <- subset(df_sub, sample_name == smp)[,i]
            
            name <- paste(names(df_sub[i]), '.mean', sep="")  # set the new column name
            summary[summary$sample_name == smp, name] <- mean(x)  # calculate mean of the replicates
        
            name <- paste(names(df_sub[i]), '.median', sep="")  # set the new column name
            summary[summary$sample_name == smp, name] <- median(x)  # calculate median of the replicates
            
            name <- paste(names(df_sub[i]), '.sd', sep="")  # set the new column name
            summary[summary$sample_name == smp, name] <- sd(x)  # calculate sd of the replicates
        }
    }
    

    # calculation of p-values
    # loop through numeric columns
    for (i in c(ilow : ihigh)) {
        # loop through unique compound names
        for (comp in unique(summary$compound_name)) {
            # select the CONTROL replicates of the compound. These are saved in data frame x and used in t-test.
            x <- df_sub[df_sub$compound_name == comp & df_sub$dose_level == "CONTROL", i]

            # loop through all dose levels for the compound. These are saved in data frame y and used in t-test.
            for (dl in c("CONTROL", "LOW", "MIDDLE", "HIGH")) {
                y <- df_sub[df_sub$compound_name == comp & df_sub$dose_level == dl, i]
                
                # calculation of p-values
                # try-catch section for t-test. If t-test returns any kind of error, then assign NA to the p-value.
                obj <- try(tt <- t.test(x,y), silent=TRUE)
                if (is(obj, "try-error")) {
                    pval <- NA
                }
                else {
                    pval <- broom::tidy(tt)$p.value
                }
                
                name <- paste(names(df_sub[i]), '.pval', sep="")  # set the new column name
                summary[summary$compound_name == comp & summary$dose_level == dl, name] <- pval  # add p-value to the summary data frame
                
                
            }
        }
    }
    
    
    # calculation of adjusted p-values
    # loop through sample names in the summary data frame. For each sample name collect p-values and compute adjusted p-values from them.
    for (smp in summary$sample_name) {
        # select all columns with ".pval"
        pvalue_columns <- grep(".pval", colnames(summary), value=TRUE)

        # create a vector of new column names (replace "pval" by "padj")
        padj_columns <- sub("pval", "padj", pvalue_columns)
        
        pval <- summary[summary$sample_name == smp, pvalue_columns]  # vector of p-values
        padj <- p.adjust(pval, method = "bonferroni")  # vector of adjusted p-values
        names(padj) <- padj_columns  # add names of the padj_columns
        
        # update summary data frame
        for (col in padj_columns) {
            summary[summary$sample_name == smp, col] <- padj[col] 
        }
    }
    
    return(summary)
}

TGB<- read.table("file:///C:/Users/tatyana/Documents/Projects/Advance project/rat_in_vivo.tsv",sep="\t", header=TRUE)
PC<- read.csv2("file:///C:/Users/tatyana/Documents/projects/KIT_phospholipidosis/Phospholipidosis_positive_list.csv")
NC<- read.csv2("file:///C:/Users/tatyana/Documents/projects/KIT_phospholipidosis/Phospholipidosis_negative_list.csv")



TGB <- as.data.frame(apply(TGB, 2, toupper))
PC <- as.data.frame(apply(PC, 2, toupper))
NC <- as.data.frame(apply(NC, 2, toupper))

TGB_PC <- merge(x=TGB, y=PC, by.x="compound_name", by.y="Compound.name")
TGB_NC <- merge(x=TGB, y=NC, by.x="compound_name", by.y="Compound.name")


TGB_PC_unique<-as.data.frame(unique(TGB_PC$compound_name))

TGB_NC_unique<-as.data.frame(unique(TGB_NC$compound_name))






# orgnas = LIVER, KIDNEY
# exposure times = 3 HR, 6 HR, 9 HR, 24 HR, 4 DAY, 8 DAY, 15 DAY, 29 DAY

organs = c("LIVER")
exp_times = c("29 DAY", "15 DAY", "8 DAY",'3 HR', '6 HR', '9 HR', '24 HR', '4 DAY')

# initiate empty data frames for summary data
SUMMARY_PC = data.frame()
SUMMARY_NC = data.frame()


# run the analysis for each organ/exposure time
for (organ in organs) {
    for (exp_time in exp_times) {
        SUMMARY_PC <- rbind(SUMMARY_PC, analyze_dataset(TGB_PC, organ, exp_time))
        SUMMARY_NC <- rbind(SUMMARY_NC, analyze_dataset(TGB_NC, organ, exp_time))
    
    }
}

organ <- "LIVER"
exp_time <- "29 DAY"

#Subset the p-values only

#Final_GC_adj<- SUMMARY_GC[,c(1,187:228)]

CountMySignificant <- function(x) {
    count = 0
    for(i in 1:length(x)) {
      if (!(is.na(x[i]))) {
        if(x[i] < 0.05) {
          count = count + 1
        }
      }
      
    }
    return(count)
  }


dose_level <- "HIGH"
exp_time <- "29 DAY"

pval_columns <- grep(".pval", colnames(SUMMARY_PC), value = TRUE)


Active_compounds = data.frame()

for (col in pval_columns) {
  df <- SUMMARY_PC[SUMMARY_PC$dose_level == dose_level & SUMMARY_PC$exposure_time == exp_time, pval_columns]
  Active_compounds["PC", col] <- CountMySignificant(df[[col]])
  
  df <- SUMMARY_NC[SUMMARY_NC$dose_level == dose_level & SUMMARY_NC$exposure_time == exp_time, pval_columns]
  Active_compounds["NC", col] <- CountMySignificant(df[[col]])
  
  
}

Counts<-t(Active_compounds)
Counts<-as.data.frame(Counts)

#Fold changes

mean_columns <- grep(".mean", colnames(SUMMARY_PC), value = TRUE)

for (col in mean_columns) {
  for (cmp in unique(SUMMARY_PC$compound_name)) {
    for (exp_time in unique(SUMMARY_PC$exposure_time)) {
      # control parameter value
      x_control <- SUMMARY_PC[SUMMARY_PC$compound_name == cmp & SUMMARY_PC$exposure_time == exp_time & SUMMARY_PC$dose_level == "CONTROL", col]
      
      for (dl in c("LOW", "MIDDLE", "HIGH")) {
        # treatment parameter value
        x_treatment <- SUMMARY_PC[SUMMARY_PC$compound_name == cmp & SUMMARY_PC$exposure_time == exp_time & SUMMARY_PC$dose_level == dl, col]
        
        # new column name
        new_col <- sub(".mean", ".fmean", col)
        
        # add factor to the new column in SUMMARY_GC
        SUMMARY_PC[SUMMARY_PC$compound_name == cmp & SUMMARY_PC$exposure_time == exp_time & SUMMARY_PC$dose_level == dl, new_col] <- x_treatment/x_control
      }
    }
  }
}




mean_columns <- grep(".mean", colnames(SUMMARY_NC), value = TRUE)

for (col in mean_columns) {
  for (cmp in unique(SUMMARY_NC$compound_name)) {
    for (exp_time in unique(SUMMARY_NC$exposure_time)) {
      # control parameter value
      x_control <- SUMMARY_NC[SUMMARY_NC$compound_name == cmp & SUMMARY_NC$exposure_time == exp_time & SUMMARY_NC$dose_level == "CONTROL", col]
      
      for (dl in c("LOW", "MIDDLE", "HIGH")) {
        # treatment parameter value
        x_treatment <- SUMMARY_NC[SUMMARY_NC$compound_name == cmp & SUMMARY_NC$exposure_time == exp_time & SUMMARY_NC$dose_level == dl, col]
        
        # new column name
        new_col <- sub(".mean", ".fmean", col)
        
        # add factor to the new column in SUMMARY_GC
        SUMMARY_NC[SUMMARY_NC$compound_name == cmp & SUMMARY_NC$exposure_time == exp_time & SUMMARY_NC$dose_level == dl, new_col] <- x_treatment/x_control
      }
    }
  }
}





 # Select only the parameters significantly deregulated in min 3 of the chemicals
Counts_PC<- Counts[Counts$PC > 3,]
Counts_NC<- Counts[Counts$NC > 3,]

Counts_PC$names <- rownames(Counts_PC)
Counts_NC$names <- rownames(Counts_NC)



Counts_PC<-subset(Counts_PC, select=c("names"))
Counts_NC<-subset(Counts_NC, select=c("names"))


#Identification of toxic group- specific parameters



Toxic_specific<-setdiff(Counts_PC, Counts_NC)



##Replacement of pval with mean in order to match properly later
Toxic_specific<-data.frame(Toxic_specific)
data<- lapply(Toxic_specific, gsub, pattern = "pval", replacement = "fmean", fixed = TRUE)
data<-as.data.frame (data)


#Filtering of the GC group for different conditions
##Toxic Compounds
Heatmap_PC<-filter(SUMMARY_PC, dose_level== "HIGH", exposure_time== "4 DAY")
Heatmap_PC <- Heatmap_PC[,-1]
rownames(Heatmap_PC) <- Heatmap_PC[,1]
Heatmap_PC<-t(Heatmap_PC)
Heatmap_PC<-as.data.frame(Heatmap_PC)
Heatmap_PC$names <- rownames(Heatmap_PC)

##Non-toxic
Heatmap_NC<-filter(SUMMARY_NC, dose_level== "HIGH",  exposure_time== "4 DAY")
Heatmap_NC <- Heatmap_NC[,-1]
rownames(Heatmap_NC) <- Heatmap_NC[,1]
Heatmap_NC<-t(Heatmap_NC)
Heatmap_NC<-as.data.frame(Heatmap_NC)
Heatmap_NC$names <- rownames(Heatmap_NC)


#merge with the mean fold changes
Final_toxic_specific<-merge(data,Heatmap_PC, by="names")
Final_toxic_specific<-merge(Final_toxic_specific,Heatmap_NC, by="names")

#Final_NGTX_specific<- lapply(Final_NGTX_specific, gsub, pattern = ".fmean", replacement = " ", fixed = TRUE)

#colnames(Final_toxic_specific)<- c("parameter","AMI","CHL", "CHP","CLO","FLU","GEN","HAL",
 #                                 "IMI","KETO","PER","PRO","TAM",
  #                                "THI","ACE", "CAF","CAR", "DIA"
  #                                ,"DIC","DOX","FEN","FLU",
  #                                "FUR","GEM","ISO"
   #                           )

m <- (apply(Final_toxic_specific[, -1],2,as.numeric))
rownames(m) = as.character( Final_toxic_specific[,1])

library(RColorBrewer)
library("devtools")
source("C://users/tatyana/Documents/Projects/Advance project/heamap.3.R")

cola <- (c("red", "red",  "red", "red",  "red", "red",
           "red", "red", "red","red","red","red","red",
            "green", "green","green","green","green",
           "green","green","green",
           "green","green","green"))
cola <- as.matrix(cola)
library(gplots)

heatmap.3(m, col= redgreen(75), scale = "none", ColSideColors=cola, 
          key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol = 0.8, 
          keysize = 1, sepwidth = c(0.05, 0.05))


#Explore all parameters deregulated at least in 3 conditions

Parameters<- rbind(Counts_PC, Counts_NC)

Parameters<- lapply(Parameters, gsub, pattern = "pval", replacement = "fmean", fixed = TRUE)
Parameters<-as.data.frame (Parameters)


#Add the fold changes to all
Final<-merge(Parameters,Heatmap_PC, by="names")
Final<-unique (merge(Final,Heatmap_NC, by="names"))






m1 <- (apply(Final[, -1],2,as.numeric))
rownames(m1) = as.character( Final[,1])




cola <- as.matrix(cola)
library(gplots)

heatmap.3(m1, col= redgreen(20), scale = "none", ColSideColors=cola, 
          key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol = 0.8, 
          keysize = 1, sepwidth = c(1,2,3, 4,5,6))


SUMMARY_PC$group <- "PC"
SUMMARY_NC$group <- "NC"



TOTAL <- rbind(SUMMARY_PC, SUMMARY_NC)

TOTAL[TOTAL$exposure_time %in% exp_time & TOTAL$dose_level == dl, ]$Ht.fmean


exp_time <- c("29 DAY", "15 DAY", "8 DAY", '4 DAY')
dl <- "HIGH"

TOTAL$exposure_time_factor = factor(TOTAL$exposure_time, levels=c('4 DAY', '8 DAY','15 DAY','29 DAY'))


# estimate limits
min <- min(TOTAL[TOTAL$exposure_time %in% exp_time & TOTAL$dose_level == dl, ]$Ht.fmean)
max <- max(TOTAL[TOTAL$exposure_time %in% exp_time & TOTAL$dose_level == dl, ]$Ht.fmean)
print(paste("Minimum limit should be", log2(min), "or less"))
print(paste("Minimum limit should be", log2(max), "or more"))

ggplot(TOTAL[TOTAL$exposure_time %in% exp_time & TOTAL$dose_level == dl, ], aes(y=log(Ht.fmean, 2), x=compound_name)) +
  scale_y_continuous(name="log2 FC Ht", limits=c(- 2, 2), breaks=c(-2,0,2), labels=c("-2", "0", "2")) +
  geom_bar(stat="identity", colour="white", aes(fill=group)) + 
  coord_flip() +
  facet_grid(rows=vars(group), cols=vars(exposure_time_factor), scale="free_y") +
  scale_fill_brewer(palette= "Dark2") +
  theme(axis.title.y=element_blank())

ggsave("Ret_29Days_MIDDLE.tiff")









