#combining Ms Cc temp var experiment (2017) data, with the 30C replication data (2020)


#load libraries
library(scales)
library(readr)
library(nlme)
library(lme4)
library(lmerTest)
library(ggplot2)
library(car)
library(tidyr)
library(mgcv)
library(dplyr)




#load data (using wide data sets that include dead individuals, to compare mortality; long data sets are clean)
tv <- read_csv("data files/25-28-30_tv-final_ed.csv")

tvr <- read_csv("data files/Ms-Cc_tv-rep-30_ed.csv")


#----------------------------------

#look at column names to see if they match

cn_tv <- colnames(tv)
cn_tvr <- colnames(tvr)

View(cn_tv)
View(cn_tvr)

#remove field individuals from tv
tv <- subset(tv, pop=="lab")

#rename columns in tv data set to match tvr data set
tv <- tv %>% rename(tt3 = timeto3, tt4 = timeto4, tt5 = timeto5, ps.ecl = tot.surv, 
                    ps.em = tot.elsurv, ttem.h = ttem, ttem.w = waspdev.int, ttecl = waspdev.tot)


#remove columns about field collection information from tv (have removed all field individuals), remove extraneous columns that I don't use in analyses
tv <- tv %>% select(-c("pop", "date.coll.j", "coll.loc", "date.laid.j", "int3", "int4", "int5", 
                       "use", "stsp.llsurv", "stsp.pupsurv", "tot.llsurv", "end"))


#create died and died.bf5 columns for tv 
tv$date.died.j[is.na(tv$date.died.j)]<-0
tv$date.5.j[is.na(tv$date.5.j)] <- 0
tv$died <- ifelse(tv$date.died.j > 0, 1, 0)
tv$died.bf5 <- ifelse(tv$date.died.j > 0 & tv$date.5.j==0, 1, 0)


#add date.cull and mass.cull, left_out and ttcull columns to tv (so it matches tvr)--values will all be NA
tv$date.cull.j <- NA
tv$mass.cull <- NA
tv$left_out <- NA
tv$ttcull <- NA

#create a column indicating if an individual is categorised as a mongo or not. (in 30.10 para treat, has date died, but did not die before 5th)
tv$mongo <- ifelse(tv$treatment=="para" & tv$temp.avg==30 & tv$temp.var==10 & tv$died==1 & tv$died.bf5==0, 1, 0)

#create an end.class column, classifying if that individual wandered, died as mongo (calling dead_cull to 
#keep it clear with tvr data set), had em, or died before 5th and was not a mongo

tv$date.wand.j[is.na(tv$date.wand.j)]<-0
tv$date.em.j[is.na(tv$date.em.j)]<-0

tv$end.class <- ifelse(tv$date.wand.j > 0, "wand",
                       ifelse(tv$date.em.j > 0, "em",
                              ifelse(tv$mongo==1, "dead_cull", "dead")))


#create column indicating which experiment each row of data is from--here, "orig"
tv$expt <- "orig"




#remove p5 data from tvr data set, dates that are not julian day and inst.6 column
tvr <- tvr %>% select(-c("date.p5.1", "mass.p5.1", 
                         "date.p5.2", "mass.p5.2",
                         "date.p5.3", "mass.p5.3"))

#remove date columns that are not julian day
tvr <- tvr %>% select(-c("date.3", "date.4", "date.5", "date.cull", "date.wand", "date.em", "date.ecl", 
                         "date.died", "date.hatch", "date.ovp", "date.p5.1.j", "date.p5.2.j", "date.p5.3.j", "inst.6"))


#rename num column as bug.id, and convert to character
tvr <- tvr %>% rename("bug.id" = "num")
tvr$bug.id <- as.character(tvr$bug.id)

#create a diet column, that indicates all individuals fed on colony diet
tvr$diet <- "col"


#create mongo column that correctly categorizes if an individual was a mongo (treat para, date cull or died after 5th)
tvr$date.cull.j[is.na(tvr$date.cull.j)]<-0
tvr$mongo <- ifelse(tvr$date.cull.j > 0, 1, 
                    ifelse(tvr$treatment=="para" & tvr$temp.var==10 & tvr$died==1 & tvr$died.bf5==0, 1, 0))
  
#remove notes columns  
tvr$notes <- NULL
tvr$diss.notes <- NULL

#create column indicating which experiment each row of data is from--here, "repl"
tvr$expt <- "repl"


#make sure column names match
cn_tv2 <- colnames(tv)
cn_tvr2 <- colnames(tvr)

cn_tv2 %in% cn_tvr2


#bind into one data frame
tv_all <- bind_rows(tv, tvr)


#create clean version, create sorting column "keep"
tv_all$keep <- ifelse(tv_all$died==0, 1,
                      ifelse(tv_all$died==1 & tv_all$died.bf5==0 & tv_all$mongo==1, 1, 0))

#subset to only those with 1 in keep column
tv_all_cl <- subset(tv_all, keep==1)


#--------------

#create a long dataframe from combined, clean data sets

tv_almass <- gather(tv_all_cl, instar, mass, mass.3, mass.4, mass.5, mass.end)
tv_almass$instar <- gsub("mass.", "", tv_almass$instar)

tv_alage <- gather(tv_all_cl, instar, age, tt3, tt4, tt5, ttend)
tv_alage$instar <- gsub("tt", "", tv_alage$instar)
tv_alage <- select(tv_alage, bug.id, temp.avg, temp.var, treatment, instar, age)

tv_allng <- merge(tv_almass, tv_alage, by=c("bug.id", "temp.avg", "temp.var", "treatment", "instar"))
View(tv_allng)


#write to csv
write.csv(tv_all, "Ms-Cc_tv-orig-rep_comb_wdead.csv", row.names = FALSE)
write.csv(tv_all_cl, "Ms-Cc_tv-orig-rep_comb_cl.csv", row.names = FALSE)
write.csv(tv_allng, "Ms-Cc_tv-orig-rep_comb_lng.csv", row.names = FALSE)

