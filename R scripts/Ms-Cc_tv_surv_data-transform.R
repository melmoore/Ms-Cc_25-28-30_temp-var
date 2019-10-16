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

#sort out dead individuals
tv$date.died.j[is.na(tv$date.died.j)]<-0
tv_cl <- subset(tv, date.died.j==0)

tv_lng$date.died.j[is.na(tv_lng$date.died.j)]<-0
tv_lngc <- subset(tv_lng, date.died.j==0)

#making a sorting/event column to indicate the event of wandering or wasp emergence
tv_lngc$end_evnt <- ifelse(tv_lngc$instar=="end", 1, 0)


tv_sel <- select(tv_cl, bug.id, treatment, temp.avg, temp.var, ttend)

test <- tmerge(data1 = tv_sel, data2 = tv_lngc, id=bug.id, end = event(end_evnt))
test <- tmerge(data1 = test, data2 = tv_lngc, id=bug.id, instr = tdc(day.age, instar))

#running into problems because of mongos, need to sort better so can remove them from both dataframes

