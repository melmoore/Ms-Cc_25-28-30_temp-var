##ANALYSES FOR Ms Cc TEMPERATURE VARIATION MANUSCRIPT
##FINAL VERSIONS--USING REPL 30 DATA ONLY, NO ORIG 30 DATA



#load libraries
library(scales)
library(Rmisc)
library(readr)
library(nlme)
library(lme4)
library(lmerTest)
library(car)
library(tidyr)
library(mgcv)
library(dplyr)



#load data
tvor <- read_csv("data files/Ms-Cc_tv-orig-rep_comb_cl.csv", 
                 col_types = cols(temp.avg = col_factor(levels = c("25", "28", "30"))))


tvor_lng <- read_csv("data files/Ms-Cc_tv-orig-rep_comb_lng.csv",
                     col_types = cols(temp.avg = col_factor(levels = c("25", "28", "30"))))



#make sorting column to remove orig 30 data
tvor$keep <- ifelse(tvor$temp.avg==30 & tvor$expt=="orig", 0, 1)
tvor_lng$keep <- ifelse(tvor_lng$temp.avg==30 & tvor_lng$expt=="orig", 0, 1)

tvor <- subset(tvor, keep==1)
tvor_lng <- subset(tvor_lng, keep==1)


#----------------------

#necessary data cleaning and configuration for plots

#log mass
tvor_lng$log_mss <- log(tvor_lng$mass)


#remove individuals with overly large loads (>300)
tvor$keep_ld <- ifelse(tvor$treatment=="para" & tvor$load > 300, 0, 1)
tvor_lng$keep_ld <- ifelse(tvor_lng$end.class=="em" & tvor_lng$load > 300, 0, 1)

tvor <- subset(tvor, keep_ld==1)
tvor_lng <- subset(tvor_lng, keep_ld==1)

#remove wanderers and WOWEs that wandered
tvor$keep_p <- ifelse(tvor$treatment=="para" & tvor$end.class=="wand", 0, 1)
tvor_lng$keep_p <- ifelse(tvor_lng$treatment=="para" & tvor_lng$end.class=="wand", 0, 1)

tvor <- subset(tvor, keep_p==1)
tvor_lng <- subset(tvor_lng, keep_p==1)


#remove WOWEs in 30C treatment
tvor$keep_p2 <- ifelse(tvor$temp.avg==30 & tvor$temp.var==0 & tvor$treatment=="para" & tvor$end.class!="em", 0, 1)
tvor_lng$keep_p2 <- ifelse(tvor_lng$temp.avg==30 & tvor_lng$temp.var==0 & tvor_lng$treatment=="para" & tvor_lng$end.class!="em", 0, 1)

tvor <- subset(tvor, keep_p2==1)
tvor_lng <- subset(tvor_lng, keep_p2==1)


#make temp.var a factor, removing the now non existant +/-5 treatment
tvor$temp.var <- factor(tvor$temp.var, levels = c(0, 10))
tvor_lng$temp.var <- factor(tvor_lng$temp.var, levels = c(0, 10))


#create tot.died column for analysing wasp survival 
tvor$tot.died <- tvor$load - tvor$num.ecl

#make subsets that have only parasitized treatments
tvor_p <- subset(tvor, treatment=="para")
tvor_lngp <- subset(tvor_lng, treatment=="para")



#------------------------

#ANALYSIS OF HOST MASS USING GAMM MODELS

#subset to only columns in model, remove rows with NAs (so that predicted and fitted values can be added
#to the dataframe easily)
tvor_mass <- select(tvor_lng, bug.id, temp.avg, temp.var, treatment, log_mss, age)
tvor_mass <- na.omit(tvor_mass)

#make bug.id a factor so it will work as a random effect in the GAMM model
tvor_mass$bug.id <- as.factor(tvor_mass$bug.id)


#Full GAMM model 
#Response = log of mass. Smooth by age, with an interaction of parasitization treatment, mean temperature and
#temperature fluctuation. Smooth of individual added to act as a random effect. Parasitization treatment, mean
#temperature and temperature fluctuation used as fixed effects. 
gam_mass_mod<-gam(log_mss ~ s(age, by = interaction(treatment, temp.avg, temp.var, k=20, bs="ts")) 
                  + s(bug.id, bs ="re") + treatment * temp.avg * temp.var,
                  method="ML", data=tvor_mass, na.action = na.omit)


anova(gam_mass_mod)
summary(gam_mass_mod)

gam.check(gam_mass_mod, type = "deviance")



#add predicted and residual values to model data set, plot results
tvor_mass$pred <- predict(gam_mass_mod, level=0)
tvor_mass$resid <- residuals(gam_mass_mod, level=0)


#plot model residuals against age, color by mean temperature, facet by parasitization treatment and fluctuation
md_gammod_ra<-ggplot(tvor_mass, aes(x=age, y=resid, color=temp.avg))
md_gammod_ra+geom_point(size=4, shape=1
)+geom_hline(aes(yintercept=0),
             color="black",
             size=1.5, linetype="dashed"
)+facet_wrap(treatment~temp.var)



md_gammod_fit<-ggplot(tvor_mass, aes(x=age, y=log_mss, group=interaction(bug.id, temp.avg), color=temp.avg))
md_gammod_fit+geom_point(size=3, shape=1
)+geom_line(aes(y=pred, group=interaction(bug.id, temp.avg))
)+facet_wrap(treatment~temp.var)





