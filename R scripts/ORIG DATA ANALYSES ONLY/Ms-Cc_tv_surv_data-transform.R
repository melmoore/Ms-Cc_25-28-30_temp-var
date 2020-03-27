#Ms Cc temp var survival analysis for time and mass at end of dev (wandering or wasp emergence)
  ##DATA TRANSFORMATION SCRIPT


#load libraries

library(readr)
library(tidyr)
library(dplyr)
library(survival)


#load data

tv<-read_csv("data files/25-28-30_tv-final_clean.csv",
             col_types = cols(temp.avg = col_factor(levels = c("25","28", "30")), 
                              temp.var = col_factor(levels = c("0", "5", "10")), 
                              treatment = col_factor(levels = c("control","para"))))

tv_lng <- read_csv("data files/25-28-30_tv-final_clean_LONG.csv", 
                   col_types = cols(temp.avg = col_factor(levels = c("25", "28", "30")), 
                                    temp.var = col_factor(levels = c("0", "5", "10")), 
                                    treatment = col_factor(levels = c("control", "para"))))


#-------------------------------------------------------

#practice with examples from tutorial

View(cgd0)

newcgd <- tmerge(data1=cgd0[, 1:13], data2=cgd0, id=id, tstop=futime)
newcgd <- tmerge(newcgd, cgd0, id=id, infect = event(etime1))
newcgd <- tmerge(newcgd, cgd0, id=id, infect = event(etime2))
newcgd <- tmerge(newcgd, cgd0, id=id, infect = event(etime3))
newcgd <- tmerge(newcgd, cgd0, id=id, infect = event(etime4))
newcgd <- tmerge(newcgd, cgd0, id=id, infect = event(etime5))
newcgd <- tmerge(newcgd, cgd0, id=id, infect = event(etime6))
newcgd <- tmerge(newcgd, cgd0, id=id, infect = event(etime7))
newcgd <- tmerge(newcgd, newcgd, id, enum=cumtdc(tstart))


View(jasa)

jasa$subject <- 1:nrow(jasa)


tdata <- with(jasa, data.frame(subject = subject,
                               futime= pmax(.5, fu.date - accept.dt),
                               txtime= ifelse(tx.date== fu.date,(tx.date -accept.dt) -.5,
                                              (tx.date - accept.dt)),fustat = fustat))


sdata <- tmerge(jasa, tdata, id=subject,
                death = event(futime, fustat),
                trt   =  tdc(txtime),
                options= list(idname="subject"))



temp <- subset(pbc, id <= 312, select=c(id:sex, stage))
pbc2 <- tmerge(temp, temp, id=id, death = event(time, status))
pbc2 <- tmerge(pbc2, pbcseq, id=id, 
               ascites = tdc(day, ascites),
               bili = tdc(day, bili), 
               albumin = tdc(day, albumin),
               protime = tdc(day, protime), 
               alk.phos = tdc(day, alk.phos))

#------------------------

#attempting to create a start,stop data set for use in survival analysis

#to create a stop start dataset, need a wide format and a long format that have the same unique IDs
  ##problem is wowes, which sort out as dead in the wide dataframe, but not in the long data frame
  ##some 30+/-10 bugs wandered, which are sorted out in the long data frame (because of the wowe column),
  ##but not the wide.

  ##for initial analysis, removing all "para" treatment bugs that wandered, because it is unknown if there
  ##was a failed ovp, rescue effect, human error, or other that caused it.


#sort out dead and para wandered individuals
tv$date.died.j[is.na(tv$date.died.j)]<-0
tv_cl <- subset(tv, date.died.j==0)

tv_cl$mass.wand[is.na(tv_cl$mass.wand)]<-0
tv_cl$para_wand <- ifelse(tv_cl$treatment=="para" & tv_cl$mass.wand>0, 1, 0)
tv_cl <- subset(tv_cl, para_wand==0)
tv_cl$mass.wand[tv_cl$mass.wand==0]<-NA


#make a sorting column for wowes
tv_lng$date.em.j[is.na(tv_lng$date.em.j)]<-0
tv_lng$wowe <- ifelse(tv_lng$temp.avg=="30" & tv_lng$temp.var=="10" & tv_lng$date.em.j=="0"
                      & tv_lng$treatment=="para", 1, 0)
tv_lng$date.em.j[tv_lng$date.em.j==0]<-NA

tv_lngc <- subset(tv_lng, wowe==0)

#making a sorting column for para wanderers
tv_lngc$mass.wand[is.na(tv_lngc$mass.wand)]<-0
tv_lngc$para_wand <- ifelse(tv_lngc$treatment=="para" & tv_lngc$mass.wand>0, 1, 0)
tv_lngc <- subset(tv_lngc, para_wand==0)
tv_lngc$mass.wand[tv_lngc$mass.wand==0]<-NA


#making a sorting/event column to indicate the event of wandering or wasp emergence
tv_lngc$end_evnt <- ifelse(tv_lngc$instar=="end", 1, 0)

tv_cl$date.wand.j[is.na(tv_cl$date.wand.j)]<-0
tv_cl$end_ev <- ifelse(tv_cl$date.wand.j>0, "wand", "em")


tv_sel <- select(tv_cl, bug.id, treatment, temp.avg, temp.var, ttend, end_ev)

test <- tmerge(data1 = tv_sel, data2 = tv_sel, id=bug.id, end=event(ttend, end_ev))
test <- tmerge(data1 = test, data2 = tv_lngc, id=bug.id, instar=tdc(day.age, instar))
test <- tmerge(data1 = test, data2 = tv_lngc, id=bug.id, mass=tdc(day.age, mass))




