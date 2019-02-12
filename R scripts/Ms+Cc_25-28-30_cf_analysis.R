#Moore--Ms+Cc constant and fluctuating temp (25, 28, 30) experiment--ANALYSIS


#load libraries

library(scales)
library(readr)
library(nlme)
library(lme4)
library(lmerTest)
library(ggplot2)



#load data

tv<-read_csv("data files/25-28-30_tv-final_clean.csv",
             col_types = cols(temp.avg = col_factor(levels = c("25","28", "30")), 
                              temp.var = col_factor(levels = c("0", "5", "10")), 
                              treatment = col_factor(levels = c("control","para"))))

View(tv)


tv.long <- read_csv("data files/25-28-30_tv-final_clean_LONG.csv",
                    col_types = cols(temp.avg = col_factor(levels = c("25", "28", "30")), 
                                     temp.var = col_factor(levels = c("0","5", "10")), 
                                     treatment = col_factor(levels = c("control","para"))))

View(tv.long)



#creating a column for log.mass in tv.long

tv.long$log.mass<-log(tv.long$mass)


#---------------------------------------------


#BUILDING A MIXED EFFECTS MODEL FOR CATERPILLAR MASS

##Basing the model on the model I created with James and Joel (KU model) for CxPxT experiment
##age interacts with all terms as a linear and quadratic term, either temp.avg will be a stand alone, 
##to account for the effect of temperature on mass at 3rd (intercept).
             

##James said the model doesn't run with the +/- 5 treatment because there are probably too many variables, can't converge. 
  ##Suggested taking out +/-5, model seems to work without


tv.long.no5<-subset(tv.long,temp.var!=5)


##MODEL 1--no +/- 5 fluctuation
##log.mass is numeric, day.age is integer, treat and temp.avg and temp.var are factors
##random effect is linear day age by individual
  ##temp.avg is included as a fixed effect that doesn't interact with age to account for differing intercepts (mass at 3rd)
  ##due to different rearing temperatures
    ##syntax of (x+y+z)^2 removes the 4 way interactions, leaving only 3 way interactions with the different age terms



lms.mod1<-lme(log.mass~(day.age+I(day.age^2)):(temp.var+treatment+temp.avg)^2+temp.avg,
             random=~day.age|bug.id,
             data=tv.long.no5,
             na.action=na.omit,
             method="ML")

anova(lms.mod1)
summary(lms.mod1)


#Trying model with just 30 temps

tv.long.30<-subset(tv.long, temp.avg=="30")

lms.mod2<-lme(log.mass~(day.age+I(day.age^2)):(temp.var+treatment)^2,
              random=~day.age|bug.id,
              data=tv.long.30,
              na.action=na.omit,
              method="ML")

anova(lms.mod2)
summary(lms.mod2)


#trying model of just para treatments with load

tv.long.para<-subset(tv.long.no5, treatment=="para")

lms.mod2<-lme(log.mass~(day.age+I(day.age^2)):(temp.var+load+temp.avg)^2+temp.avg,
              random=~day.age|bug.id,
              data=tv.long.para,
              na.action=na.omit,
              method="ML",
              control = lmeControl(opt='optim'))

anova(lms.mod2)
summary(lms.mod2)

#---------------------------

#MODELLING MASS AT END OF DEV (WAND OR EM)

tv$lmend<-log(tv$mass.end)

lmend.mod1<-lm(lmend~temp.avg*temp.var*treatment,
               data=tv,
               na.action=na.omit)

anova(lmend.mod1)
summary(lmend.mod1)


#without +/-5

tv.no5<-subset(tv, temp.var!=5)

lmend.mod2<-lm(lmend~temp.avg*temp.var*treatment,
               data=tv.no5,
               na.action=na.omit)

anova(lmend.mod2)
summary(lmend.mod2)


#-------------------------------------------------------------------------------------------------------------------------------


#MODELLING WASP TOTAL SURVIVAL

#subset out +/-5 treatment

tv.no5<-subset(tv,temp.var!=5)
tv.no5$temp.var<-factor(tv.no5$temp.var)

##subset to only parasitized treatment

tv.para<-subset(tv.no5,treatment=="para")


#Making a column for total died (load-num.ecl)

tv.para$tot.died<-tv.para$load-tv.para$num.ecl


#Making mongo tot.died==NA (doesn't make sense to analyze them as 0s)--ask Joel and James about this

#Making tot.died==NA for mongos (load==0)

tv.para$tot.died<-ifelse(tv.para$load=="0", 1.5, tv.para$tot.died)
tv.para$tot.died[tv.para$tot.died=="1.5"]<-NA

tv.para$load[tv.para$load==0]<-NA


#making a glm model of wasp total survival to test for overdispersion
  ##temp.avg==factor, temp.var==factor, load==numeric
    ###overdispersion is high, should add random effect of individual-James says this is fine for now, due to problems with
    ###running the glmer (conversion problems due to multiple issues--see lab notebook)

wtots.mod1<-glm(cbind(num.ecl,tot.died)~temp.avg*temp.var*resc.ld,
                family=quasibinomial,
                data=tv.para,
                na.action = na.omit)

anova(wtots.mod1,test="F")
summary(wtots.mod1)



#rescaling load, tot.died and num.ecl

tv.para$resc.ld<-rescale(tv.para$load,to=c(0,1))
tv.para$man.resc.ld<-tv.para$load/321

tv.para$temp.var.num<-as.numeric(tv.para$temp.var)

tv.para.no10<-subset(tv.para,temp.avg!="30" & temp.var!="10")


#Making a glmer (binomial) model of wasp total survival
  ##temp.avg==factor, temp.var==numeric, resc.load==numeric, mongos treated as NAs until I determine a better way to deal with them
    ###won't run without rescaling load
    ###as of 2.27.18, won't run even when load is rescaled. Don't know what's going on

wtots.mod2<-glmer(cbind(tot.died,num.ecl)~temp.avg*temp.var*resc.ld+(1|bug.id),
                family=binomial,
                data=tv.para.no10,
                na.action = na.omit,
                control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))


anova(wtots.mod2,test="F")
summary(wtots.mod2)


wtots.mod3<-glm(cbind(num.ecl,tot.died)~temp.avg*temp.var*load,
                family=quasibinomial,
                data=tv.para.no10,
                na.action = na.omit)

anova(wtots.mod1,test="F")
summary(wtots.mod1)


#running model with load not rescaled
wtots.mod2a<-glmer(cbind(tot.died,num.ecl)~temp.avg*temp.var+load+(1|bug.id),
                  family=binomial,
                  data=tv.para,
                  na.action = na.omit,
                  control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))


anova(wtots.mod2a,test="F")
summary(wtots.mod2a)



#running model without 30+/-10 treatment
  ##won't run due to: "Error in `contrasts<-`(`*tmp*`, value = contr.funs[1 + isOF[nn]]) : 
  ##contrasts can be applied only to factors with 2 or more levels"

wtots.mod3<-glmer(cbind(tot.died,num.ecl)~temp.avg*temp.var+resc.ld+(1|bug.id),
                  family=binomial,
                  data=tv.para.no10,
                  na.action = na.omit,
                  control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))


anova(wtots.mod3,test="F")
summary(wtots.mod3)



#modelling wasp survival to emergence   This model runs!

tv.para$num.em[tv.para$num.em==0]<-NA
tv.para$num.unem[tv.para$num.unem==0]<-NA

wems.mod1<-glmer(cbind(num.unem,num.em)~temp.avg*temp.var*resc.ld+(1|bug.id),
                 family=binomial,
                 data=tv.para,
                 na.action=na.omit,
                 control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

anova(wems.mod1)
summary(wems.mod1)



#comparing models with only 1 variable against full model to get p values

wems.modTA<-glmer(cbind(num.unem,num.em)~temp.avg+(1|bug.id),
                  family=binomial,
                  data=tv.para,
                  na.action=na.omit,
                  control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

wems.modTV<-glmer(cbind(num.unem,num.em)~temp.var+(1|bug.id),
                    family=binomial,
                    data=tv.para,
                    na.action=na.omit,
                    control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

wems.modL<-glmer(cbind(num.unem,num.em)~resc.ld+(1|bug.id),
                 family=binomial,
                 data=tv.para,
                 na.action=na.omit,
                 control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))


wems.mod.TAV<-glmer(cbind(num.unem,num.em)~(temp.var:temp.avg)+(1|bug.id),
                    family=binomial,
                    data=tv.para,
                    na.action=na.omit,
                    control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))


wems.mod.TAL<-glmer(cbind(num.unem,num.em)~(resc.ld:temp.avg)+(1|bug.id),
                    family=binomial,
                    data=tv.para,
                    na.action=na.omit,
                    control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))


wems.mod.TVL<-glmer(cbind(num.unem,num.em)~(resc.ld:temp.var)+(1|bug.id),
                    family=binomial,
                    data=tv.para,
                    na.action=na.omit,
                    control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

wems.mod.TATVL<-glmer(cbind(num.unem,num.em)~(resc.ld:temp.var:temp.avg)+(1|bug.id),
                      family=binomial,
                      data=tv.para,
                      na.action=na.omit,
                      control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))


#testing temp.avg

anova(wems.mod1, wems.modTA, test="chi")
anova(wems.mod1, wems.modTV, test="chi")
anova(wems.mod1, wems.modL, test="chi")
anova(wems.mod1, wems.mod.TAV, test="chi")
anova(wems.mod1, wems.mod.TAL, test="chi")
anova(wems.mod1, wems.mod.TVL, test="chi")
anova(wems.mod1, wems.mod.TATVL, test="chi")


#---------------------------------------------------------------------------------------

#MODELLING WASP DEVELOPMENT TIME


#Making a linear mixed effect model for wasp internal development time
  ##random effect==random intercept of individual
    ###temp.avg==factor, temp.var==numeric, resc.load==numeric, mongos treated as NAs 

wdint.mod1<-lme(waspdev.int~temp.avg*temp.var*load,
               random=~1|bug.id,
               data=tv.para,
               na.action = na.omit,
               method="ML")

anova(wdint.mod1)
summary(wdint.mod1)


#Testing full vs reduced models
  
#removing 3 way interaction
wdint.mod1a<-lme(waspdev.int~temp.avg+temp.var+load+
                            temp.avg:temp.var+
                            temp.avg:load+
                            temp.var:load,
                random=~1|bug.id,
                data=tv.para,
                na.action = na.omit,
                method="ML")

anova(wdint.mod1a)


anova(wdint.mod1,wdint.mod1a)

#removing temp.avg:load
wdint.mod1b<-lme(waspdev.int~temp.avg+temp.var+load+
                   temp.avg:temp.var+
                   temp.var:load,
                 random=~1|bug.id,
                 data=tv.para,
                 na.action = na.omit,
                 method="ML")

anova(wdint.mod1b)

anova(wdint.mod1,wdint.mod1a,wdint.mod1b)


#removing load
wdint.mod1c<-lme(waspdev.int~temp.avg+temp.var+
                   temp.avg:temp.var+
                   temp.var:load,
                 random=~1|bug.id,
                 data=tv.para,
                 na.action = na.omit,
                 method="ML")

anova(wdint.mod1c)

anova(wdint.mod1, wdint.mod1a, wdint.mod1b, wdint.mod1c)


#removing temp.var:load
wdint.mod1d<-lme(waspdev.int~temp.avg+temp.var+
                   temp.avg:temp.var,
                 random=~1|bug.id,
                 data=tv.para,
                 na.action = na.omit,
                 method="ML")

anova(wdint.mod1d)

anova(wdint.mod1, wdint.mod1a, wdint.mod1b, wdint.mod1c, wdint.mod1d)

#reduced model without load or its interactions seems best


#Making a linear mixed effect model for wasp total development time
  ##random effect==random intercept of individual
    ###temp.avg==factor, temp.var==numeric, resc.load==numeric, mongos treated as NAs

#Won't run because there is no waspdev.tot for 30+/-10: all wasps died after emergence. Ask Joel/James what to do with this

wdtot.mod1<-lme(waspdev.tot~temp.avg*temp.var*resc.ld,
                random=~1|bug.id,
                data=tv.para.no10,
                na.action=na.omit,
                method="ML")

anova(wdtot.mod1)
summary(wdtot.mod1)





