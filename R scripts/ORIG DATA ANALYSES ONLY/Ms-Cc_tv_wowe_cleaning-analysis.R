# Ms Cc Temp Var WOWE data cleaning and analysis script

#load libraries

library(Rmisc)
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)

#load data
wowe <- read_csv("data files/Ms_Cc_TV_30-10_WOWE_raw_data.csv",
                 col_types = cols(temp.avg = col_factor(levels = c("25","28", "30")), 
                                  temp.var = col_factor(levels = c("0", "5", "10")), 
                                  treat = col_factor(levels = c("control","para"))))

tv<-read_csv("data files/25-28-30_tv-final_clean.csv",
             col_types = cols(temp.avg = col_factor(levels = c("25","28", "30")), 
                              temp.var = col_factor(levels = c("0", "5", "10")), 
                              treatment = col_factor(levels = c("control","para"))))


#-------------------------

#remove the info about what temp they were put in (extra list in columns), put in separate object 

info_list <- wowe[, 21:23]

wowe <- wowe[, -(21:23)]


#--------------------------

#COMBINE FULL DATA AND WOWE TEMP DATA FRAMES

#getting full data for WOWE individuals from tv by subsetting to the 30+/-10 treatment

#subset to only 30+/-10 treatment
full <- subset(tv, temp.avg==30 & temp.var==10)

#subset to only parasitized individuals
full <- subset(full, treatment=="para")


#wowe is 51 rows, full is 57. See which ones aren't in wowe
full$match <- full$bug.id %in% wowe$bug.id


#I think most of the ones not in wowe died before I started doing this
#check with table
table(full$ttend, full$match)


#rename full$treatment to full$treat so it matches wowe
full <- rename(full, treat=treatment)

#join data frames
full_wowe <- full_join(wowe, full, by=c("bug.id", "treat", "temp.avg", "temp.var"))

#-------------------------

#EXPLORE DATA

#See how many WOWEs were put in each temp treatment after 10 days (I think?), not sure what after.10 means

#assuming WOWEs that were not in wowe temp data frame were just kept at 30+/-10--changing after.10 column
#to reflect this
full_wowe$after.10[is.na(full_wowe$after.10)]<-0
full_wowe$after.10 <- as.character(full_wowe$after.10)
full_wowe$after.10 <- ifelse(full_wowe$after.10=="0", "30.10", full_wowe$after.10)

#tabulate the number in each after.10 group
table(full_wowe$after.10)

#remove individuals that had wasp emergence
full_wowe$date.em.j[is.na(full_wowe$date.em.j)] <- 0
full_wowe <- subset(full_wowe, date.em.j==0)


#find mean and summary statistics for ttend for each group
ttend_sum <- summarySE(full_wowe, measurevar = "ttend",
                       groupvars = "after.10",
                       na.rm = TRUE)
ttend_sum


#see range of ttend for those at 25 and those at 30.10
fw_hot <- subset(full_wowe, after.10=="30.10")
fw_cool <- subset(full_wowe, after.10=="25.00")

range(fw_hot$ttend)
range(fw_cool$ttend)


#---------------------

#do same calculations as above, but with field individuals removed

#remove field individuals
full_wowe_nof <- subset(full_wowe, pop=="lab")


#assuming WOWEs that were not in wowe temp data frame were just kept at 30+/-10--changing after.10 column
#to reflect this
full_wowe_nof$after.10[is.na(full_wowe_nof$after.10)]<-0
full_wowe_nof$after.10 <- as.character(full_wowe_nof$after.10)
full_wowe_nof$after.10 <- ifelse(full_wowe_nof$after.10=="0", "30.10", full_wowe_nof$after.10)

#tabulate the number in each after.10 group
table(full_wowe_nof$after.10)

#remove individuals that had wasp emergence
full_wowe_nof$date.em.j[is.na(full_wowe_nof$date.em.j)] <- 0
full_wowe_nof <- subset(full_wowe_nof, date.em.j==0)


#find mean and summary statistics for ttend for each group
ttend_sum <- summarySE(full_wowe_nof, measurevar = "ttend",
                       groupvars = "after.10",
                       na.rm = TRUE)
ttend_sum


#see range of ttend for those at 25 and those at 30.10
fw_hot <- subset(full_wowe_nof, after.10=="30.10")
fw_cool <- subset(full_wowe_nof, after.10=="25.00")

range(fw_hot$ttend)
range(fw_cool$ttend)



#----------

#see how long individuals were at 25C or 30.10 after deciding to move them to diff temp treats

#change date.mass1 into jdate

##Converts x into julian date
j.date<-function(x){
  strptime(x, "%m/%d")$yday+1
}


full_wowe_nof$date.mass.1 <- j.date(full_wowe_nof$date.mass.1)

full_wowe_nof$time_in_trt <- full_wowe_nof$date.died.j - full_wowe_nof$date.mass.1


time_sum <- summarySE(full_wowe_nof, measurevar = "time_in_trt",
                      groupvars = "after.10",
                      na.rm = TRUE)
time_sum


#-------------------

#plotting life span data

#make a frequency density plot to show life span for each group (with field)

ttend_plot <- ggplot(full_wowe, aes(x=ttend, group=after.10, color=after.10))
ttend_plot + geom_density()



#make a frequency density plot to show life span for each group (without field)

ttend_plot_nof <- ggplot(full_wowe_nof, aes(x=ttend, group=after.10, color=after.10))
ttend_plot_nof + geom_density()


#--------------------

#comparing 30 const data from cpt (reg diet) and tv (tobacco diet)

cpt_30 <- subset(cpt, temp==30)
tv_30 <- subset(tv, temp.avg==30 & temp.var==0)
tv_30 <- subset(tv_30, pop=="lab")
tv_30 <- drop_na(tv_30, mass.end)

tvmass_sum <- summarySE(tv_30, measurevar = "mass.end", 
                        groupvars = "treatment",
                        na.rm = TRUE)
tvmass_sum


cpt_npmss_sum <- summarySE(cpt_30, measurevar = "mass.wander",
                           groupvars = "treatment",
                           na.rm = TRUE)
cpt_npmss_sum



cpt_pmss_sum <- summarySE(cpt_30, measurevar = "mass.befem",
                          groupvars = "treatment",
                          na.rm = TRUE)
cpt_pmss_sum




tv_30$prop_em <- tv_30$num.em / tv_30$load
cpt_30$prop_em <- cpt_30$num.em / cpt_30$load


tv_wspsurv_sum <- summarySE(tv_30, measurevar = "prop_em",
                            groupvars = "treatment",
                            na.rm = TRUE)
tv_wspsurv_sum



cpt_wspsurv_sum <- summarySE(cpt_30, measurevar = "prop_em",
                             groupvars = "treatment",
                             na.rm = TRUE)
cpt_wspsurv_sum

