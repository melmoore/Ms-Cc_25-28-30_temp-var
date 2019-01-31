

library(readr)
library(plyr)
library(ggplot2)
library(Rmisc)
library(dplyr)
library(tidyr)
library(reshape2)

theme_set(theme_classic())

#28 data
hill <- read_csv("~/Manduca expts/Summer+Fall 2017/Hill rt and hs/data files/hill rt hs ed LOAD.csv")
View(hill)


#25 and 30 data
mem<-read_csv("~/Manduca expts/Summer+Fall 2017/25.30 wasp temp var/data files/25-30_cf_final_1-26-18.csv", 
              col_types = cols(temp.avg = col_factor(levels = c("25", "30")), 
                               temp.var = col_factor(levels = c("0", "5", "10"))))
View(mem)


#Renaming hill data columns to match mem colnames

colnames(hill)<-c("bug.id","num","temp.avg","temp.var","treatment","laid",
                  "date.hatch","date.3","mass.3","date.4","mass.4","hs.4",
                  "date.5","mass.5","hs.5","date.wand","mass.wand","date.det",
                  "mass.48em","hs.treat","date.ovp","num.ovp","wasp.gen","suc.ovp","date.em",
                  "instar.em","bled.em","num.em","num.coc","num.fail.spin","date.ecl",
                  "num.ecl","num.ecl.fail","num.unem","load","hs",
                  "date.laid.j","date.3.j","date.4.j","date.5.j","date.wand.j",
                  "date.det.j","date.ovp.j","date.em.j","date.ecl.j","laidto4",
                  "laidto5","laidtowand","waspdev.tot","waspdev.int","stsp.llsurv")

#Creating columns in hill for diet and population data

hill$diet<-"col"
hill$pop<-"lab"
hill$date.coll.j<-NA
hill$coll.loc<-NA
hill$date.died.j<-NA




#Converting mem dates to Julian

mem$date.hatch.j<-strptime(mem$date.hatch, "%m/%d")$yday+1
mem$date.died.j<-strptime(mem$date.died, "%m/%d")$yday+1
mem$date.3.j<-strptime(mem$date.3, "%m/%d")$yday+1
mem$date.4.j<-strptime(mem$date.4, "%m/%d")$yday+1
mem$date.5.j<-strptime(mem$date.5, "%m/%d")$yday+1
mem$date.wand.j<-strptime(mem$date.wand, "%m/%d")$yday+1
mem$date.ovp.j<-strptime(mem$date.ovp, "%m/%d")$yday+1
mem$date.em.j<-strptime(mem$date.em, "%m/%d")$yday+1
mem$date.ecl.j<-strptime(mem$date.ecl, "%m/%d")$yday+1
mem$date.coll.j<-strptime(mem$date.coll, "%m/%d")$yday+1



#Calculating 25 and 30 dev times

mem$timeto3<-mem$date.3.j-mem$date.hatch.j
mem$int3<-mem$date.4.j-mem$date.3.j
mem$timeto4<-mem$date.4.j-mem$date.hatch.j
mem$int4<-mem$date.5.j-mem$date.4.j
mem$timeto5<-mem$date.5.j-mem$date.hatch.j


mem$ttwand<-mem$date.wand.j-mem$date.hatch.j
mem$ttem<-mem$date.em.j-mem$date.hatch.j


mem$waspdev.int<-mem$date.em.j-mem$date.ovp.j
mem$waspdev.tot<-mem$date.ecl.j-mem$date.ovp.j
mem$waspdev.coc<-mem$date.ecl.j-mem$date.em.j


mem$date.laid.j<-NA

#Calculating wasp surv prop

mem$stsp.llsurv<-mem$num.coc/mem$num.em
mem$stsp.pupsurv<-mem$num.ecl/mem$num.coc

mem$tot.elsurv<-mem$num.em/mem$load
mem$tot.llsurv<-mem$num.coc/mem$load
mem$tot.surv<-mem$num.ecl/mem$load  #Same thing as tot.pupsurv



#Making an approximate date.hatch.j column, by adding 1 to the date.laid.j column (talk to Chr about protocol she used)

hill$date.hatch.j<-hill$date.laid.j+1


#Calculating hill dev times, using hatch date instead of laid date:

hill$timeto3<-hill$date.3.j-hill$date.hatch.j
hill$int3<-hill$date.4.j-hill$date.3.j
hill$timeto4<-hill$date.4.j-hill$date.hatch.j
hill$int4<-hill$date.5.j-hill$date.4.j
hill$timeto5<-hill$date.5.j-hill$date.hatch.j


hill$ttwand<-hill$date.wand.j-hill$date.hatch.j
hill$ttem<-hill$date.em.j-hill$date.hatch.j


hill$stsp.pupsurv<-hill$num.ecl/hill$num.coc

hill$tot.elsurv<-hill$num.em/hill$load
hill$tot.llsurv<-hill$num.coc/hill$load
hill$tot.surv<-hill$num.ecl/hill$load



#renaming factor levels in hill:

hill$treatment<-revalue(hill$treatment, c("NP"="control", "P"="para"))


#Removing pupal mass (mass_det [now mass48.em] for wanderers)

hill$mass.48em[is.na(hill$mass.48em)]<-0
hill$date.wand.j[is.na(hill$date.wand.j)]<-0

hill$mass.48em<-ifelse(hill$treatment=="control" | hill$date.wand.j>0,0,hill$mass.48em )

hill$mass.48em[(hill$mass.48em=="0")]<-NA
hill$date.wand.j[(hill$date.wand.j=="0")]<-NA


#subsetting hill to only hs controls (none that were heatshocked)

hill.nhs<-subset(hill, hs=="C")



#Making both data sets have the same number/name of columns:


keepvars<-c("bug.id","temp.avg","temp.var","treatment","mass.3","mass.4","mass.5","mass.wand","mass.48em","num.ovp","suc.ovp","instar.em","bled.em","num.em",
            "num.coc","num.fail.spin","num.ecl","num.unem","load","date.laid.j","date.hatch.j","date.died.j","date.3.j","date.4.j","date.5.j","date.wand.j","date.ovp.j","date.em.j","date.ecl.j",
            "waspdev.tot","waspdev.int","stsp.llsurv","stsp.pupsurv","tot.elsurv","tot.llsurv","tot.surv",
            "timeto3","timeto4","timeto5","int3","int4","ttwand","ttem","date.coll.j","coll.loc","pop","diet")

hill.keep<-hill.nhs[keepvars]
mem.keep<-mem[keepvars]


tv<-rbind(hill.keep,mem.keep)
View(tv)



write.csv(tv,"25-28-30_RAW_final.csv",row.names = FALSE)




#Creating a column "end" that has the julian day of when the data for each caterpillar stopped being recorded
#(date wander, date emergence, or date death)

tv$end<-coalesce(tv$date.died.j,tv$date.wand.j)
tv$end<-coalesce(tv$end,tv$date.em.j)

tv<- tv[!is.na(tv$end),]

#Creating a column for the length of the 5th instar by using "end"--calculates 5th length regardless of treatment

tv$int5<-(tv$end-tv$date.5.j)




#Creating binary columns to eventually sort data (can't sort by death, since 30+/-10 all have death dates, but want their data)


#Subset by treatment

tv.para<-subset(tv,treatment=="para")


#Make NA's == 0, since logical arguments won't run properly with NAs

tv.para$date.5.j[is.na(tv.para$date.5.j)] <- 0
tv.para$date.died.j[is.na(tv.para$date.died.j)] <- 0
tv.para$date.em.j[is.na(tv.para$date.em.j)] <- 0
tv.para$date.wand.j[is.na(tv.para$date.wand.j)] <- 0


#Logical nested ifelse statement that creates a binary column for use for the 30+/-10 para treatment
#First part says give a 1 if temp.avg is 30 and temp.var is 10, otherwise give a 0
#Second part (nested in 30.10 treatment group), give a 1 if date.5>0, and either date.died, date.em or date.wander is
#greater than 0


tv.para$use.30f10<-ifelse(tv.para$temp.avg=="30" & tv.para$temp.var=="10",
                          ifelse(tv.para$date.5.j>0 & tv.para$date.died.j>0 | tv.para$date.em.j>0 |
                                   tv.para$date.wand.j>0,1,0),0)


#Logical ifelse statements for each temp.avg and temp.var group, creating a binary column for each
#If date died == 0 (did not die), put a 1. If not, (cat died, has julian date), put a 0

tv.para$use.28c<-ifelse(tv.para$temp.avg=="28" & tv.para$temp.var=="0" & tv.para$date.died.j=="0",1,0)
tv.para$use.28f<-ifelse(tv.para$temp.avg=="28" & tv.para$temp.var=="10" & tv.para$date.died.j=="0",1,0)
tv.para$use.25c<-ifelse(tv.para$temp.avg=="25" & tv.para$temp.var=="0" & tv.para$date.died.j=="0",1,0)
tv.para$use.25f<-ifelse(tv.para$temp.avg=="25" & tv.para$temp.var=="10" & tv.para$date.died.j=="0",1,0)
tv.para$use.30c<-ifelse(tv.para$temp.avg=="30" & tv.para$temp.var=="0" & tv.para$date.died.j=="0",1,0)
tv.para$use.30f<-ifelse(tv.para$temp.avg=="30" & tv.para$temp.var=="5" & tv.para$date.died.j=="0",1,0)

#creating a final use column to sort all data with: if any of the other use columns have a 1, put a 1 in this column
#if none of the other columns have a 1, put a 0

tv.para$use<-ifelse(tv.para$use.30f10=="1" | tv.para$use.28c=="1" | tv.para$use.28f=="1" | tv.para$use.25c=="1" |
                      tv.para$use.25f=="1" | tv.para$use.30c=="1" | tv.para$use.30f=="1",1,0)


#Comparing relevant columns to double check everything worked

test<-tv.para[,c(1,2,3,22,50,51:57)]


#Subset by treatment

tv.con<-subset(tv,treatment=="control")


#Make NA's == 0, since logical arguments won't run properly with NAs

tv.con$date.died.j[is.na(tv.con$date.died.j)] <- 0


#Logical ifelse statements for each temp.avg and temp.var group, creating a binary column for each
#If date died == 0 (did not die), put a 1. If not, (cat died, has julian date), put a 0

tv.con$use.28c<-ifelse(tv.con$temp.avg=="28" & tv.con$temp.var=="0" & tv.con$date.died.j=="0",1,0)
tv.con$use.28f<-ifelse(tv.con$temp.avg=="28" & tv.con$temp.var=="10" & tv.con$date.died.j=="0",1,0)
tv.con$use.25c<-ifelse(tv.con$temp.avg=="25" & tv.con$temp.var=="0" & tv.con$date.died.j=="0",1,0)
tv.con$use.25f<-ifelse(tv.con$temp.avg=="25" & tv.con$temp.var=="10" & tv.con$date.died.j=="0",1,0)
tv.con$use.30c<-ifelse(tv.con$temp.avg=="30" & tv.con$temp.var=="0" & tv.con$date.died.j=="0",1,0)
tv.con$use.30f<-ifelse(tv.con$temp.avg=="30" & tv.con$temp.var=="5" & tv.con$date.died.j=="0",1,0)
tv.con$use.30f10<-ifelse(tv.con$temp.avg=="30" & tv.con$temp.var=="10" & tv.con$date.died.j=="0",1,0)


#creating a final use column to sort all data with: if any of the other use columns have a 1, put a 1 in this column
#if none of the other columns have a 1, put a 0

tv.con$use<-ifelse(tv.con$use.28c=="1" | tv.con$use.28f=="1" | tv.con$use.25c=="1" | tv.con$use.25f=="1" |
                     tv.con$use.30c=="1" | tv.con$use.30f=="1" | tv.con$use.30f10=="1",1,0)
                     


#rbind the subsetted data sets (with new use columns) back into 1 dataset

tv<-rbind(tv.con,tv.para)                     


#Look at relevant columns to double check things worked

test.tv<-tv[,c(1,2,3,22,24,25,48,50,51,52,53,54,55,56,57)]              


#set extraneous use columns to NULL to get rid of them
                     
tv$use.25c<-NULL
tv$use.25f<-NULL
tv$use.28c<-NULL
tv$use.28f<-NULL
tv$use.30c<-NULL
tv$use.30f<-NULL
tv$use.30f10<-NULL


#Revert 0s back to NAs

tv$date.5.j[(tv$date.5.j=="0")] <- NA
tv$date.died.j[(tv$date.died.j=="0")] <- NA
tv$date.em.j[(tv$date.em.j=="0")] <- NA
tv$date.wand.j[(tv$date.wand.j=="0")]<-NA

#Calculating length of the 4th instar for caterpillars that had wasp emergence in the 4th instar


#Rename old int4 column

tv<-rename(tv,int4.1 = int4)


#Set NAs to 0s so logical statements will work properly

tv$instar.em[is.na(tv$instar.em)] <- 0
tv$date.em.j[is.na(tv$date.em.j)] <- 0
tv$date.4.j[is.na(tv$date.4.j)] <- 0


#Logical statement that creates a 2nd int4 column--if the column instar.em == 4, then subtract the date of emergence
#from the date of 4th instar molt to find the length of the 4th instar. 
#If instar.em does not equal 4, put a 0

tv$int4.2<-ifelse(tv$instar.em=="4",(tv$date.em.j-tv$date.4.j),0)


#Convert 0s back to NAs

tv$int4.2[(tv$int4.2=="0")] <- NA
tv$date.em.j[(tv$date.em.j=="0")] <- NA
tv$date.4.j[(tv$date.4.j=="0")] <- NA
tv$instar.em[(tv$instar.em=="0")]<-NA

#Combine the 2 int4 columns into on int4 column using coalesce

tv$int4<-coalesce(tv$int4.1,tv$int4.2)


#Remove extraneous int4 columns

tv$int4.1<-NULL
tv$int4.2<-NULL



#Making wasp number data (emerged, unemerged, etc) be 0 for 30.10 treatment group, for analyses and figures

#Turn NAs to 0s for relevant columns (used date.wand for para with suc.ovp==0)

tv$num.em[is.na(tv$num.em)] <- 0
tv$num.coc[is.na(tv$num.coc)] <- 0
tv$num.fail.spin[is.na(tv$num.fail.spin)] <- 0
tv$num.ecl[is.na(tv$num.ecl)] <- 0
tv$num.unem[is.na(tv$num.unem)] <- 0
tv$date.wand.j[is.na(tv$date.wand.j)]<-0


#Ifelse statements that gives any cat with control in treatment column (or a date.wand > 0) a wasp metric of 1.5
#Leaves everything else the same (else==tv$num.em, etc)
#Turns anything in a wasp metric column with value of 1.5 into an NA

tv$num.em<-ifelse( tv$treatment=="control" |  tv$date.wand.j>0,1.5, tv$num.em)
tv$num.em[( tv$num.em=="1.5")]<-NA

tv$num.coc<-ifelse( tv$treatment=="control" |  tv$date.wand.j>0,1.5, tv$num.coc)
tv$num.coc[( tv$num.coc=="1.5")]<-NA

tv$num.fail.spin<-ifelse( tv$treatment=="control" |  tv$date.wand.j>0,1.5, tv$num.fail.spin)
tv$num.fail.spin[( tv$num.fail.spin=="1.5")]<-NA

tv$num.ecl<-ifelse( tv$treatment=="control" |  tv$date.wand.j>0,1.5, tv$num.ecl)
tv$num.ecl[( tv$num.ecl=="1.5")]<-NA

tv$num.unem<-ifelse( tv$treatment=="control" |  tv$date.wand.j>0,1.5, tv$num.unem)
tv$num.unem[( tv$num.unem=="1.5")]<-NA

tv$load<-ifelse( tv$treatment=="control" |  tv$date.wand.j>0,1.5, tv$load)
tv$load[( tv$load=="1.5")]<-NA

#Turning date.wand 0s back into NAs

tv$date.wand.j[(tv$date.wand.j=="0")]<-NA


#Recalculating wasp prop surv columns so 30.10 group have 0s instead of NAs

tv$stsp.llsurv<-tv$num.coc/tv$num.em
tv$stsp.pupsurv<-tv$num.ecl/tv$num.coc

tv$tot.elsurv<-tv$num.em/tv$load
tv$tot.llsurv<-tv$num.coc/tv$load
tv$tot.surv<-tv$num.ecl/tv$load


tv$stsp.llsurv[is.nan(tv$stsp.llsurv)] = 0
tv$stsp.pupsurv[is.nan(tv$stsp.pupsurv)] = 0
tv$tot.elsurv[is.nan(tv$tot.elsurv)] = 0
tv$tot.llsurv[is.nan(tv$tot.llsurv)] = 0
tv$tot.surv[is.nan(tv$tot.surv)] = 0


#Making control and suc.ovp==0 individuals have NAs for wasp metrics

#Make date.wand NAs == 0
tv$date.wand.j[is.na(tv$date.wand.j)]<-0 


#Ifelse statements--if treatment is control, or date.wand is greater than 0 (for individuals with unsuccessful ovp),
#put a 1.5. If no, leave as is (tv$num.em, etc)

tv$num.em<-ifelse(tv$treatment=="control" | tv$date.wand.j>0 ,1.5,tv$num.em)
tv$num.coc<-ifelse(tv$treatment=="control" | tv$date.wand.j>0 ,1.5,tv$num.coc)
tv$num.fail.spin<-ifelse(tv$treatment=="control" | tv$date.wand.j>0 ,1.5,tv$num.fail.spin)
tv$num.ecl<-ifelse(tv$treatment=="control" | tv$date.wand.j>0 ,1.5,tv$num.ecl)
tv$num.unem<-ifelse(tv$treatment=="control" | tv$date.wand.j>0 ,1.5,tv$num.unem)


#Convert 1.5s to NAs

tv$num.em[( tv$num.em=="1.5")]<-NA
tv$num.coc[( tv$num.coc=="1.5")]<-NA
tv$num.fail.spin[( tv$num.fail.spin=="1.5")]<-NA
tv$num.ecl[( tv$num.ecl=="1.5")]<-NA
tv$num.unem[( tv$num.unem=="1.5")]<-NA


#Converting date.wand back to NA

tv$date.wand.j[(tv$date.wand.j=="0")]<-NA


#Creating a time to end column for figures and analysis

tv$ttend<-tv$end-tv$date.hatch.j


#Adding a mass.end column, combining mass at wandering and mass at 48 em (may need to get mass at death for 30.10)

tv$mass.48em<-as.numeric(tv$mass.48em)
tv$mass.end<-coalesce(tv$mass.48em,tv$mass.wand)



#Subset whole data frame by the use column--data should now be clean (ish)

tv.cl<-subset(tv,use=="1")

View(tv.cl)


#Removing para caterpillars that have a suc.ovp==0 (were not successfully parasitized)

#Make NAs==0 so ifelse will work
tv.cl$suc.ovp[is.na(tv.cl$suc.ovp)]<-0
tv.cl$mass.wand[is.na(tv.cl$mass.wand)]<-0
tv.cl$date.em.j[is.na(tv.cl$date.em.j)]<-0

#Make an ifelse statement that will put a 1.5 if treatment==control or mass.wand==0 and date.em.j==0, 
#leave else as is

tv.cl$suc.ovp<-ifelse(tv.cl$treatment=="control" | tv.cl$mass.wand==0 & tv.cl$date.em.j==0, 1.5, tv.cl$suc.ovp)


#Remove suc.ovp==0

tv.cl<-subset(tv.cl,suc.ovp!=0)


#Turn suc.ovp that == 1.5 back to NA, return NA's to other relevant columns

tv.cl$suc.ovp[tv.cl$suc.ovp==1.5]<-NA
tv.cl$mass.wand[tv.cl$mass.wand==0]<-NA
tv.cl$date.em.j[tv.cl$date.em.j==0]<-NA


#Write csv's

write.csv(tv,"25-28-30_tv-final_ed.csv",row.names = FALSE)
write.csv(tv.cl,"25-28-30_tv-final_clean.csv",row.names=FALSE)




#Creating a long data set from tv.cl

#Gathering development time

test.long<-gather(tv.cl,instar,day.age,timeto3,timeto4,timeto5,ttend)
View(test.long)

#Whittling down to only necessary columns (for merging later), and renaming instar levels

test.long1<-select(test.long,bug.id,treatment,temp.avg,temp.var,instar,day.age)
test.long1$instar<- gsub("timeto", "",test.long1$instar)
test.long1$instar<- gsub("tt", "",test.long1$instar)

#Gathering mass data

test.long2<-gather(tv.cl,instar,mass,mass.3,mass.4,mass.5,mass.end)
test.long2$instar<-gsub("mass.","",test.long2$instar)

#Merging the two long data sets by id, treatments and instar

tv.cl.long<-merge(test.long1,test.long2,by=c("bug.id","treatment","temp.avg","temp.var","instar"))




write.csv(tv.cl.long,"25-28-30_tv-final_clean_LONG.csv",row.names=FALSE)




