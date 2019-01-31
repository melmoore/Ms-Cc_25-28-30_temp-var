#Calculating survival for all Ms+Cc_temp-var treatments

#load libraries

library(readr)
library(plyr)
library(ggplot2)
library(Rmisc)
library(dplyr)
library(tidyr)
library(reshape2)


#------------------------

#load data

tv_raw <- read_csv("data files/25-28-30_RAW_final.csv")
View(tv_raw)


#--------------------

#Make binary column for if individual died

#make date.died.j NAs==0 for ifelse statement
tv_raw$date.died.j[is.na(tv_raw$date.died.j)]<-0

#ifelse statement to create binary "died" column
tv_raw$died<-ifelse(tv_raw$date.died.j>0, 1, 0)

#----------------------------

#count the TOTAL number in each treatment, and the % survival in each treatment
  ##this includes field bugs, many of whom died--may need to remove

tv_raw<-subset(tv_raw, pop=="lab")

#create a table with the total number in each treatment
tot.table<-table(tv_raw$temp.avg, tv_raw$temp.var, tv_raw$treatment)
View(tot.table)

#create a table with the number that died in each treatment
d.table<-table(tv_raw$temp.avg, tv_raw$temp.var, tv_raw$treatment, tv_raw$died)
View(d.table)

#rename table columns to more useful names

tot.table<-as.data.frame(tot.table)
tot.table <- tot.table %>% rename(temp.avg=Var1, temp.var=Var2, treat=Var3, treat.tot=Freq)

d.table<-as.data.frame(d.table)
d.table <- d.table %>% rename(temp.avg=Var1, temp.var=Var2, treat=Var3, died=Var4, smn=Freq)

#to calcuate number that survived in each treatment, subset d.table to only died==0
ds.table<-subset(d.table, died==0)

#add total number from tot.table to ds.table
ds.table$treat.tot<-tot.table$treat.tot

#calculate percent survived in each treatment
ds.table$perc.surv<-ds.table$smn/ds.table$treat.tot


