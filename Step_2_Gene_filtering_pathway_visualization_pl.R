#Author: Tatyana Doktorova and Noffisat Oki
#Version: 12_12_2018
# The script filters the most important genes and facilitates their visulaisation



library(igraph)
library (dplyr)
library(VennDiagram)
library(treemap)
library(data.tree)
library(data.table)

library(DiagrammeR)
library(collapsibleTree)
library(RColorBrewer)

setwd("C:\\Users\\tatyana\\Documents\\Projects\\KIT_Phospholipidosis")

#Upload disease relevant files

TC_perc <- read.csv("file:///C:/Users/tatyana/Documents/Projects/KIT_Phospholipidosis/TC_perc_lm.txt",sep = '\t', header=TRUE)
CTD_perc <- read.csv("file:///C:/Users/tatyana/Documents/Projects/KIT_Phospholipidosis/CTD_perc_lm.txt",sep = '\t', header=TRUE)
CTDTC_perc <- read.csv("file:///C:/Users/tatyana/Documents/Projects/KIT_Phospholipidosis/CTDTC_perc_lm.txt",sep = '\t', header=TRUE)

#Select the unique categories of interest

TC_perc_sub<- subset(TC_perc, select=c("Parent.Name","lhs" , "Major", "Minor"))
CTD_perc_sub<- subset(CTD_perc, select=c("Parent.Name","lhs" , "Major", "Minor"))
CTDTC_perc_sub<- subset(CTDTC_perc, select=c("Parent.Name","lhs",  "Major", "Minor"))

TC_perc_parent = unique(TC_perc_sub)
CTD_perc_parent = unique(CTD_perc_sub)
CTDTC_perc_parent = unique(CTDTC_perc_sub)

# Indentify unique genes only 

TC_gene<- unique(subset(TC_perc, select=c("lhs"  )))
CTDTC_gene<- unique(subset(CTDTC_perc, select=c("lhs"  )))
CTD_gene<- unique(subset(CTD_perc, select=c("lhs"  )))

TC_perc_gene = TC_perc_sub %>% distinct(lhs)
CTD_perc_gene = CTD_perc_sub %>% distinct(lhs)
CTDTC_perc_gene = CTDTC_perc_sub %>% distinct(lhs)

#Look for overlapping genes

TC_CTDTC_CTD_overlap<- Reduce(intersect, list(TC_perc_gene$lhs,CTD_perc_gene$lhs, CTDTC_perc_gene$lhs ))
TC_CTDTC_overlap<-Reduce(intersect, list(TC_perc_gene$lhs, CTDTC_perc_gene$lhs ))
TC_CTD_overlap<-Reduce(intersect, list(TC_perc_gene$lhs, CTD_perc_gene$lhs ))
CTDTC_CTD_overlap<-Reduce(intersect, list(CTDTC_perc_gene$lhs, CTD_perc_gene$lhs ))

#Merge all unique genes
all<-unique(rbind(TC_gene,CTDTC_gene, CTD_gene))

#Overlaps as dataframes
TC_CTDTC_CTD_overlap_df <- as.data.frame(TC_CTDTC_CTD_overlap)
TC_CTDTC_overlap_df <- as.data.frame(TC_CTDTC_overlap)
TC_CTD_overlap_df <- as.data.frame(TC_CTD_overlap)
CTDTC_CTD_overlap_df <- as.data.frame(CTDTC_CTD_overlap)


#Merge the genes plus respective pathway info together
merge_all<-unique(rbindlist(list(CTD_perc_parent,TC_perc_parent,CTDTC_perc_parent)))
write.csv(merge_all, file="HCC_Toxcast_CTD_lp.csv")

ALL = merge_all %>% distinct(lhs)
#ifelse(merge_all$Parent.Name==merge_all$Minor,1)

gene_major_subset<-unique(data.frame(merge_all$"lhs",merge_all$"Major"))



HCC<-data.frame(merge_all)

#tree diagram (interactive)


collapsibleTree(df=HCC, c("Major", "Minor","Parent.Name", "lhs"),  fill = "lightgreen")
             

collapsibleTree(df=merge_all, c("Major",  "lhs"),  fill = "lightpink")
                
                
collapsibleTree(df=merge_all, c("lhs","Parent.Name","Minor","Major"),  fill = "green")


gene_subset<-unique(data.frame(merge_all$"lhs"))

major_subset<-unique(data.frame(merge_all$"Major"))

minor_subset<-unique(data.frame(merge_all$"Minor"))


GENE<- filter(merge_all,lhs == "ACTB")
collapsibleTree(df=GENE, c("lhs","Major","Minor","Parent.Name"),  fill = "green")

GENE1<- filter(merge_all,lhs == "XRCC1")
collapsibleTree(df=GENE1, c("lhs","Major","Minor","Parent.Name"),  fill = "green")

#igraph 
#Gene major subset visualisation

gene_major_subset <- unique(gene_major_subset[,2:1])
colnames(gene_major_subset) <- c("V1", "V2")
g <- graph.data.frame(gene_major_subset)
E(g)$curved <- 0
E(g)$label <- rep(1, nrow(gene_major_subset))

g <- simplify(g, remove.multiple = T, remove.loops = T)
plot.igraph(g, vertex.size=0, edge.arrow.size=0 ,
            layout=-layout.reingold.tilford(g)[,2:1])
l <- layout_with_kk(g)

tkplot(g, canvas.width = 1500, canvas.height = 1000,edge.arrow.size=.2, edge.color="orange",
       vertex.color="orange", vertex.frame.color="#ffffff",
       vertex.label=V(g)$Major, vertex.label.font=6, vertex.label.color="black", edge.curved=.1, layout=l*2.0, dim=3)
tkplot(g, canvas.width = 1500, canvas.height = 750, vertex.color="lightpink", vertex=10, size=8, vertex.size=5)

library (rgl)
rglplot(g, canvas.width = 1500, canvas.height = 1000,edge.arrow.size=.2, edge.color="orange",
        vertex.color="orange", vertex.frame.color="#ffffff",
        vertex.label=V(g)$Major, vertex.label.font=6, vertex.label.color="black", edge.curved=.1, layout=l*2.0)



