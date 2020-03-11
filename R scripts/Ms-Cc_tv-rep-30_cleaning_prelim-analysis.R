#TVR Replication of mean 30C treatments (spring 2020), cleaning, prelim analyses, adding to full data set

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


#load data

tvr <- read_csv("data files/Ms-Cc_tv_rep_30_no-wowe-diss_3-9-2020.csv", 
                col_types = cols(temp.var = col_factor(levels = c("0", "10")), 
                                 treatment = col_factor(levels = c("control", "para"))))
View(tvr)



#-----------------------

#convert dates to julian day

##Converts x into julian date
j.date<-function(x){
  strptime(x, "%m/%d")$yday+1
}


#Takes all columns that have "date." in the name, and converts contents to Julian day using j.date function. Renames columns (adds a 
##j to end of column name), and binds the out put julian day columns to the original data set

lapj.date<-function(df){
  date.j<-lapply(df[,grep("date.",colnames(df))],j.date)
  date.j<-as.data.frame(date.j)
  colnames(date.j)<-paste(colnames(date.j), "j", sep = ".")
  output.df<-cbind(df,date.j)
  output.df
}



tvr <- lapj.date(tvr)


#create a "died" and a "died.bf5" sorting column

tvr$date.died.j[is.na(tvr$date.died.j)]<-0
tvr$date.5.j[is.na(tvr$date.5.j)] <- 0
tvr$died <- ifelse(tvr$date.died.j > 0, 1, 0)
tvr$died.bf5 <- ifelse(tvr$date.died > 0 & tvr$date.5.j==0, 1, 0)


#------------------------


#calculate development times

#for host/caterpillars
tvr$tt3 <- tvr$date.3.j - tvr$date.hatch.j
tvr$tt4 <- tvr$date.4.j - tvr$date.hatch.j
tvr$tt5 <- tvr$date.5.j - tvr$date.hatch.j
tvr$ttwand <- tvr$date.wand.j - tvr$date.hatch.j
tvr$ttem.h <- tvr$date.em.j - tvr$date.hatch.j
tvr$ttcull <- tvr$date.cull.j - tvr$date.hatch.j


#for wasps
tvr$ttem.w <- tvr$date.em.j - tvr$date.ovp.j
tvr$ttecl <- tvr$date.ecl.j - tvr$date.ovp.j


#Wasp survival metrics
tvr$num.em[is.na(tvr$num.em)] <- 0
tvr$num.ecl[is.na(tvr$num.ecl)] <- 0
tvr$load[is.na(tvr$load)] <- 0

tvr$ps.em <- tvr$num.em / tvr$load
tvr$ps.ecl <- tvr$num.ecl / tvr$load

tvr$ps.em[is.nan(tvr$ps.em)] <- 0
tvr$ps.ecl[is.nan(tvr$ps.ecl)] <- 0


#combine wandering, cull and em into date and mass end columns
tvr$ttend1 <- coalesce(tvr$ttwand, tvr$ttem.h)
tvr$ttend <- coalesce(tvr$ttend1, tvr$ttcull)
tvr$ttend1 <- NULL

tvr$mass.end1 <- coalesce(tvr$mass.wand, tvr$mass.48em)
tvr$mass.end <- coalesce(tvr$mass.cull, tvr$mass.end1)
tvr$mass.end1 <- NULL


#make a class column to indicate stage of outcome (wand, cull, em)
tvr$date.wand.j[is.na(tvr$date.wand.j)] <- 0
tvr$date.em.j[is.na(tvr$date.em.j)] <- 0
tvr$date.cull.j[is.na(tvr$date.cull.j)] <- 0

tvr$end.class <- ifelse(tvr$date.wand.j > 0, "wand",
                        ifelse(tvr$date.em.j > 0, "em",
                               ifelse(tvr$date.cull.j > 0, "cull", "dead")))



#---------------------

#creating long data set with only live individuals

tvr_cl <- subset(tvr, died==0) 

tvr_lng <- tvr_cl %>% gather(instar, mass, mass.3, mass.4, mass.5, mass.end)
tvr_lng$instar <- gsub("mass.", "", tvr_lng$instar)

tvr_age <- tvr_cl %>% gather(instar, age, tt3, tt4, tt5, ttend)
tvr_age$instar <- gsub("tt", "", tvr_age$instar)
tvr_age <- select(tvr_age, num, treatment, temp.var, instar, age)

tvr_lng <- merge(tvr_lng, tvr_age, by=c("num", "treatment", "temp.var", "instar"))



#--------------------

#looking at how many died in each treatment, both before and after 5th

#combo para and control treatments
table(tvr$temp.var, tvr$died)
table(tvr$temp.var, tvr$died.bf5)

#subset into para and control
tvr_p <- subset(tvr, treatment=="para")
tvr_c <- subset(tvr, treatment=="control")


#table of deaths in para treatment--most deaths occurred here, after 5th (WOWEs that died before culling)
table(tvr_p$temp.var, tvr_p$died)
table(tvr_p$temp.var, tvr_p$died.bf5)

#table of deaths in para treatment
table(tvr_c$temp.var, tvr_c$died)
table(tvr_c$temp.var, tvr_c$died.bf5)




#-----------------------


#exploration of WOWE data

#examine 30.10 para treatment--see how many died that don't have masses
tvr_wowe <- subset(tvr_p, temp.var=="10" & died.bf5==0)


#see how many have date.died and date.culled--none
tvr_wowe$date.cull.j[is.na(tvr_wowe$date.cull.j)] <- 0
which(tvr_wowe$date.died.j > 0 & tvr_wowe$date.cull.j >0)

#see if any that had died==1 have notes about mass at death in the notes column
tvr_dw <- subset(tvr_wowe, died==1)
tvr_lw <- subset(tvr_wowe, died==0)

notes <- tvr_dw[, c("num", "died", "notes")]
View(notes)

notes2 <- tvr_lw[, c("num", "died", "notes")]
View(notes2)


#search for the phrase "culled early" in notes for WOWEs with out a date died
grep("culled early", tvr_lw$notes)

#create a column to indicate if a WOWE was culled early
tvr_lw$early_cull <- ifelse(grepl("culled early", tvr_lw$notes)=="TRUE", 1, 0)

check <- tvr_lw[, c("num", "date.p5.1", "date.p5.2", "date.p5.3", "date.cull.j", "early_cull")]
View(check)



#calculate development times for post 5th measurements
tvr_lw$ttp5.1 <- tvr_lw$date.p5.1.j - tvr_lw$date.hatch.j
tvr_lw$ttp5.2 <- tvr_lw$date.p5.2.j - tvr_lw$date.hatch.j
tvr_lw$ttp5.3 <- tvr_lw$date.p5.3.j - tvr_lw$date.hatch.j



#create long data frame of live wowes, for age and mass at each post 5th measurement and at cull
tvr_lw1 <- tvr_lw %>% gather(p5.meas, mass, mass.p5.1, mass.p5.2, mass.p5.3, mass.end)
tvr_lw1$p5.meas <- gsub("mass.", "", tvr_lw1$p5.meas)


tvrlw2 <- tvr_lw %>% gather(p5.meas, age, ttp5.1, ttp5.2, ttp5.3, ttend)
tvrlw2$p5.meas <- gsub("tt", "", tvrlw2$p5.meas)
tvrlw2 <- select(tvrlw2, num, p5.meas, age)


tvr_wlng <- merge(tvr_lw1, tvrlw2, by=c("num", "p5.meas"))


#----------------------

#investigate whether individuals that were left out at wasp emergence had a difference in wasp survival

#subset to only those that should have had emergence (tv = 0)
tvr_pem <- subset(tvr_p, temp.var==0)

#convert NAs to 0 for the left_out column
tvr_pem$left_out[is.na(tvr_pem$left_out)]<-0

#plot num_ecl against load, colored by left_out column
leftout_plot <- ggplot(tvr_pem, aes(x=load, y=num.ecl, color=as.character(left_out)))
leftout_plot + geom_point(
) + geom_smooth(method="lm")


#plot psecl against load, colored by left_out column
leftout_plot2 <- ggplot(tvr_pem, aes(x=load, y=ps.ecl, color=as.character(left_out)))
leftout_plot2 + geom_point(
) + geom_smooth(method="lm")


#quick and dirty linear model of num.ecl
leftout_mod <- lm(num.ecl ~ load * left_out,
                  data=tvr_pem,
                  na.action = na.omit)

anova(leftout_mod)
summary(leftout_mod)



#quick and dirty linear model of ps.ecl
leftout_mod2 <- lm(ps.ecl ~ left_out,
                  data=tvr_pem,
                  na.action = na.omit)

anova(leftout_mod2)
summary(leftout_mod2)


#-------------------------

#write to csv file, with datasets with and without dead individuals

write.csv(tvr, "Ms-Cc_tv-rep-30_ed.csv", row.names = FALSE)
write.csv(tvr_cl, "Ms-Cc_tv-rep-30_clean.csv",row.names = FALSE)
write.csv(tvr_lng, "Ms-Cc_tv-rep-30_cl_lng.csv", row.names = FALSE)

write.csv(tvr_wlng, "Ms-Cc_tv-rep-30_WOWES_lng.csv", row.names = FALSE)
write.csv(tvr_lw, "Ms-Cc_tv-rep-30_WOWES.csv", row.names = FALSE)



