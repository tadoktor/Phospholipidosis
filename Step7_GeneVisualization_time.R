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

TOTAL_all<-read.csv("file:///C:/Users/tatyana/Documents/projects/Advance project/TOTAL.csv",sep=",", header=TRUE)
TOTAL<- filter (TOTAL_all, Gene== "Cited2", doseLevel== "High")

TOTAL$hour <- "h"
                       
exp_time <- c("696")
dl <- "High"
TOTAL[TOTAL$timepointHr %in% exp_time & TOTAL$doseLevel == dl, ]$Zfp36

TOTAL$exposure_time_factor = factor(TOTAL$timepointHr, levels=c(  "696"))


# estimate limits
min <- min(TOTAL[TOTAL$timepointHr%in% exp_time & TOTAL$doseLevel == dl, ]$FC)
print (min)
max <- max(TOTAL[TOTAL$timepointHr %in% exp_time & TOTAL$doseLevel == dl, ]$FC)
print (max)

ggplot(TOTAL[TOTAL$timepointHr %in% exp_time & TOTAL$doseLevel == dl, ], aes(y=FC, x=Compound)) +
  scale_y_continuous(name="log2 Zfp36", limits=c(-2,2), breaks=c(-2,0,2), labels=c("-2", "0", "2")) +
  geom_bar(stat="identity", colour="white", aes(fill=group)) + 
  coord_flip() +
  facet_grid(rows=vars(group), cols=vars(exposure_time_factor), scale="free_y") +
  scale_fill_brewer(palette="Set2") +
  theme(axis.title.y=element_blank())

ggsave("Gpnmb_Middle.png")
