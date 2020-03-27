#Moore--Ms+Cc constant and fluctuating temp (25, 28, 30) experiment--ANALYSIS


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
library(viridis)

#----------------------


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


#Create datasets without the +/-5 treatment 
tv.long.no5<-subset(tv.long,temp.var!=5)
tv.no5<-subset(tv, temp.var!=5)

#create data sets with out field individuals
tv.no5 <- subset(tv.no5, pop=="lab")
tv.long.no5 <- subset(tv.long, pop=="lab")

#create datasets with only parasitized hosts

#long dataframe
tv.long.para<-subset(tv.long.no5, treatment=="para")

#remove parasitized bugs that wandered
tv.long.para$date.wand.j[is.na(tv.long.para$date.wand.j)]<-0
tv.long.para<-subset(tv.long.para, date.wand.j==0)
tv.long.para$date.wand.j[tv.long.para$date.wand.j==0]<-NA

#remove individual with load > 300
tv.long.para <- subset(tv.long.para, load < 300)


#wide dataframe
tv.para<-subset(tv.no5, treatment=="para")

#removing individuals that were parasitized but wandered
tv.para$date.wand.j[is.na(tv.para$date.wand.j)]<-0
tv.para<-subset(tv.para, date.wand.j==0)

#remove individual with load > 300
tv.para <- subset(tv.para, load < 300)


#---------------------------------------------


#BUILDING A MIXED EFFECTS MODEL FOR CATERPILLAR MASS

#Polynomial models

##Basing the model on the model I created with James and Joel (KU model) for CxPxT experiment
##age interacts with all terms as a linear and quadratic term, either temp.avg will be a stand alone, 
##to account for the effect of temperature on mass at 3rd (intercept).
             

##James said the model doesn't run with the +/- 5 treatment because there are probably too many variables, can't converge. 
  ##Suggested taking out +/-5, model seems to work without


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

lms.mod2<-lme(log.mass~(day.age+I(day.age^2)):(temp.var+load+temp.avg)^2+temp.avg,
              random=~day.age|bug.id,
              data=tv.long.para,
              na.action=na.omit,
              method="ML",
              control = lmeControl(opt='optim'))

anova(lms.mod2)
summary(lms.mod2)


#plotting residuals for model with both para treatments

msmd<-tv.long.no5
msmd<-drop_na(msmd, log.mass, day.age)


msmd$pred<-predict(lms.mod1)
msmd$resid<-residuals(lms.mod1)
msmd$pred_f<-predict(lms.mod1, level=0)


msmd_pr_fit<-ggplot(msmd, aes(x=pred, y=resid, color=temp.avg))
msmd_pr_fit+geom_point(shape=1, size=5
)+geom_hline(aes(yintercept=0), 
             color="black", linetype="dashed"
)+facet_wrap(treatment~temp.var)


msmd_pdat_fit<-ggplot(msmd, aes(x=day.age, y=log.mass, color=temp.avg))
msmd_pdat_fit+geom_point(aes(shape=treatment),
                         size=3
)+geom_line(aes(y=pred_f, group=interaction(temp.avg, treatment)),
            size=2
)+facet_wrap(treatment~temp.var)


msmd_ra_fit<-ggplot(msmd, aes(x=day.age, y=resid, color=temp.avg))
msmd_ra_fit+geom_point(shape=1, size=5
)+geom_hline(aes(yintercept=0), 
             color="black", linetype="dashed"
)+facet_wrap(treatment~temp.var)

#residuals for this polynomial model don't look that bad



#plotting residuals for model with para treatment and load

pmd<-tv.long.para
pmd<-drop_na(pmd, log.mass, day.age, load)


pmd$pred<-predict(lms.mod2)
pmd$resid<-residuals(lms.mod2)
pmd$pred_f<-predict(lms.mod2, level=0)


pmd_pr_fit<-ggplot(pmd, aes(x=pred, y=resid, color=temp.avg))
pmd_pr_fit+geom_point(shape=1, size=5
)+geom_hline(aes(yintercept=0), 
             color="black", linetype="dashed"
)+facet_wrap(~temp.var)

#can't get the fixed effects fit lines for some reason--can only get the random effects fit lines
pmd_pdat_fit<-ggplot(pmd, aes(x=day.age, y=log.mass, color=temp.avg))
pmd_pdat_fit+geom_point(size=3
)+geom_line(aes(y=pred_f, group=interaction(temp.avg, bug.id)),
            size=2
)+facet_wrap(~temp.var)


pmd_ra_fit<-ggplot(pmd, aes(x=day.age, y=resid, color=temp.avg))
pmd_ra_fit+geom_point(shape=1, size=5
)+geom_hline(aes(yintercept=0), 
             color="black", linetype="dashed"
)+facet_wrap(~temp.var)


#---------------------------------


#Trying GAMM model for caterpillar growth and dev

#subset to only columns in model, and remove rows with NAs (so that predicted and fitted values can be
#added to the dataframe easily)
tv_mass<-select(tv.long.no5, bug.id, temp.avg, temp.var, treatment, log.mass, day.age)
tv_mass<-na.omit(tv_mass)


#make bug.id a factor so it will work as a random effect in the GAMM model
tv_mass$bug.id<-as.factor(tv_mass$bug.id)


#run a full GAMM model (where the smooth of age is also affected by the interaction of temp and treat)
gam_mass_mod<-gam(log.mass ~ s(day.age, by= interaction(treatment,temp.avg, temp.var, k=20,bs="ts")) 
                  + s(bug.id,bs="re") + treatment * temp.avg * temp.var,
                  method="ML", data=tv_mass, na.action = na.omit)
anova(gam_mass_mod)
summary(gam_mass_mod)

gam.check(gam_mass_mod, type = "deviance")



#Make a null model for model testing, where the smooth of age is not affected by temp.var, temp.avg or treat
gam_mnull_mod<-gam(log.mass ~ s(day.age, k=10,bs="ts") 
                   + s(bug.id,bs="re") + treatment * temp.avg * temp.var,
                   method="ML", data=tv_mass, na.action = na.omit)

anova(gam_mnull_mod)


#run a full GAMM model (where the smooth of age is also affected by the interaction of temp and treat)
  ##leaving knots unspecified, to test against k=10
gam_mass_usk_mod<-gam(log.mass ~ s(day.age, by= interaction(treatment,temp.avg, temp.var, bs="ts")) 
                  + s(bug.id,bs="re") + treatment * temp.avg * temp.var,
                  method="ML", data=tv_mass, na.action = na.omit)
anova(gam_mass_usk_mod)
summary(gam_mass_usk_mod)

gam.check(gam_mass_usk_mod)

#compare full and null models
#full model is much better!
#models with k=10 and k unspecified don't differ
anova(gam_mass_mod, gam_mass_usk_mod, gam_mnull_mod, test="Chisq")
AIC(gam_mass_mod, gam_mass_usk_mod, gam_mnull_mod)



#make columns with predicted and residual values for plotting
tv_mass$pred<-predict(gam_mass_mod, level=0)
tv_mass$resid<-residuals(gam_mass_mod, level=0)



md_gammod_ra<-ggplot(tv_mass, aes(x=day.age, y=resid, color=temp.avg))
md_gammod_ra+geom_point(size=4, shape=1
)+geom_hline(aes(yintercept=0),
             color="black",
             size=1.5, linetype="dashed"
)+facet_wrap(treatment~temp.var)


md_gammod_fit<-ggplot(tv_mass, aes(x=day.age, y=log.mass, color=temp.avg))
md_gammod_fit+geom_point(size=3, shape=1
)+geom_line(aes(y=pred, group=temp.avg)
)+facet_wrap(treatment~temp.var)




#para treatment with load as a predictor

#subset to only columns in model, and remove rows with NAs (so that predicted and fitted values can be
#added to the dataframe easily)
p_mass<-select(tv.long.para, bug.id, temp.avg, temp.var, load, log.mass, day.age)
p_mass<-na.omit(p_mass)

#convert bug.id to factor so it functions properly as a random effect
p_mass$bug.id<-as.factor(p_mass$bug.id)


#run gam mode with age and load within a single smooth, with an interaction of temp.avg and temp.var
##check knots etc for fit
gam_pml_mod<-gam(log.mass ~ s(day.age, load, by=interaction(temp.avg, temp.var, bs="ts")) 
                 + s(bug.id,bs="re") + temp.avg * temp.var, method="ML", data=p_mass, na.action = na.omit)

anova(gam_pml_mod)
summary(gam_pml_mod)

gam.check(gam_pml_mod)

#run a GAMM with age and load as separate smooths--does not have an interaction with temp this way (2 2D surfaces, instead
## of one 3D surface)
gam_pml_nointmod<-gam(log.mass ~ s(day.age, by=interaction(temp.avg, temp.var, k=10, bs="ts")) 
                      + s(load, by=interaction(temp.avg, temp.var, k=10, bs="ts")) + s(bug.id,bs="re") 
                      + temp.avg * temp.var, method="ML", data=p_mass, na.action = na.omit)

anova(gam_pml_nointmod)
summary(gam_pml_nointmod)
#what does plot() plot for a gam object?

gam.check(gam_pml_nointmod)


#make a null model with an intercept (age and load in the same smooth)
gam_pml_null_mod<-gam(log.mass ~ s(day.age, load, bs="ts") 
                 + s(bug.id,bs="re") + temp.avg * temp.var, method="ML", data=p_mass, na.action = na.omit)


#Make a null model with no intercept (age and load in separate smooths)
gam_pml_null_nointmod<-gam(log.mass ~ s(day.age, bs="ts") 
                      + s(load, bs="ts") + s(bug.id,bs="re") 
                      + temp.avg * temp.var, method="ML", data=p_mass, na.action = na.omit)



#compare models with and without interaction of age and load in smooth
#the no interaction mod seems better
anova(gam_pml_mod, gam_pml_nointmod, gam_pml_null_mod, gam_pml_null_nointmod, test="Chisq")
AIC(gam_pml_mod, gam_pml_nointmod, gam_pml_null_mod, gam_pml_null_nointmod)


#Model with age and load in same smooth
#make columns with predicted and residual values for plotting
p_mass$pred<-predict(gam_pml_mod, level=0)
p_mass$resid<-residuals(gam_pml_mod, level=0)

#residuals don't look great--investigate more
pml_gammod_ra<-ggplot(p_mass, aes(x=day.age, y=resid, color=temp.avg))
pml_gammod_ra+geom_point(size=4, shape=1
)+geom_hline(aes(yintercept=0),
             color="black",
             size=1.5, linetype="dashed"
)+facet_wrap(~temp.var)


pml_gammod_fit<-ggplot(p_mass, aes(x=day.age, y=log.mass, color=load))
pml_gammod_fit+geom_point(size=3, shape=1
)+geom_line(aes(y=pred, group=interaction(temp.avg, bug.id))
)+scale_color_viridis(
)+facet_wrap(temp.avg~temp.var)



#Model with age and load in same smooth
#make columns with predicted and residual values for plotting
p_mass$pred_ni<-predict(gam_pml_nointmod, level=0)
p_mass$resid_ni<-residuals(gam_pml_nointmod, level=0)


pmlni_gammod_ra<-ggplot(p_mass, aes(x=day.age, y=resid_ni, color=temp.avg))
pmlni_gammod_ra+geom_point(size=4, shape=1
)+geom_hline(aes(yintercept=0),
             color="black",
             size=1.5, linetype="dashed"
)+facet_wrap(~temp.var)


pmlni_gammod_fit<-ggplot(subset(p_mass, load<300), aes(x=day.age, y=log.mass, color=load))
pmlni_gammod_fit+geom_point(size=4, shape=1
)+geom_line(aes(y=pred_ni, group=interaction(temp.avg, bug.id)),
            size=1
)+scale_color_viridis(begin=.1, end=1,
                      option = "viridis"
)+facet_wrap(temp.avg~temp.var)


#------------------------

#MODELLING WASP TOTAL SURVIVAL--GLMM


#Making a column for total died (load-num.ecl)

tv.para$tot.died<-tv.para$load-tv.para$num.ecl


#rescaling load, tot.died and num.ecl

tv.para$resc.ld<-rescale(tv.para$load,to=c(0,1))

#tv.para$man.resc.ld<-tv.para$load/321
#tv.para$temp.var.num<-as.numeric(tv.para$temp.var)
#tv.para.no10<-subset(tv.para,temp.avg!="30" & temp.var!="10")


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




#Making a glmer (binomial) model of wasp total survival
##temp.avg==factor, temp.var==numeric, resc.load==numeric, mongos treated as NAs until I determine a better way to deal with them
###won't run without rescaling load
###as of 2.27.18, won't run even when load is rescaled. Don't know what's going on

wtots.mod2<-glmer(cbind(tot.died,num.ecl)~temp.avg*temp.var*resc.ld+(1|bug.id),
                  family=binomial,
                  data=tv.para,
                  na.action = na.omit,
                  control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))


anova(wtots.mod2,test="F")
summary(wtots.mod2)








#---------------------------

#MODELLING MASS AT END OF DEV (WAND OR EM)

tv$lmend<-log(tv$mass.end)
tv.no5$lmend<-log(tv.no5$mass.end)
tv.para$lmend<-log(tv.para$mass.end)

lmend.mod1<-lm(lmend~temp.avg*temp.var*treatment,
               data=tv,
               na.action=na.omit)

anova(lmend.mod1)
summary(lmend.mod1)


#model without the +/-5 treatment
lmend.mod2<-lme(lmend~temp.avg*temp.var*treatment,
                random = ~1|bug.id,
               data=tv.no5,
               method="ML",
               na.action=na.omit)

anova(lmend.mod2)
summary(lmend.mod2)


#removing non significant terms for model testing
lmend.mod2_red<-lme(lmend~temp.avg+temp.var+treatment + temp.avg:treatment + temp.var:treatment,
                    random = ~1|bug.id,
                    data=tv.no5,
                    method = "ML",
                    na.action=na.omit)
anova(lmend.mod2_red)


anova(lmend.mod2, lmend.mod2_red)


#Additive model for model testing

lmend.mod2_add<-lme(lmend~temp.avg+temp.var+treatment,
                    random = ~1|bug.id,
                    data=tv.no5,
                    method = "ML",
                    na.action=na.omit)


anova(lmend.mod2, lmend.mod2_add, lmend.mod2_red) #reduced model is best


#looking for only para, with load as fixed effect

lmend.mod3<-lme(lmend~temp.avg*temp.var*load,
               random=~1|bug.id,
               data=tv.para,
               na.action = na.omit,
               method="ML")

anova(lmend.mod3)


#making a reduced model for testing
lmend.mod3_add<-lme(lmend~temp.avg+temp.var+load,
                    random=~1|bug.id,
                    data=tv.para,
                    na.action = na.omit,
                    method="ML")

anova(lmend.mod3_add)

anova(lmend.mod3, lmend.mod3_add)


#plotting mass at end by load

lmend_plot<-ggplot(tv.para, aes(x=load, y=lmend, color=temp.var))
lmend_plot+geom_point(
)+geom_smooth(method="lm"
)+facet_wrap(~temp.avg)


#plotting model residuals

lmend_dat<-tv.no5
lmend_dat<-drop_na(lmend_dat, lmend)

lmend_dat$pred_f<-predict(lmend.mod2, level = 0)
lmend_dat$resid_f<-residuals(lmend.mod2, level = 0)

lmend_pr_plot<-ggplot(lmend_dat, aes(x=pred_f, y=resid_f, color=temp.avg))
lmend_pr_plot+geom_point(shape=1
)+geom_hline(aes(yintercept=0), size=1, linetype="dashed"
)+facet_wrap(treatment ~ temp.var)



plmend<-tv.para
plmend<-drop_na(plmend, lmend, load)

plmend$pred_f<-predict(lmend.mod3, level = 0)
plmend$resid_f<-residuals(lmend.mod3, level=0)

plmend_pr_plot<-ggplot(plmend, aes(x=pred_f, y=resid_f, color=temp.avg))
plmend_pr_plot+geom_point(shape=1
)+geom_hline(aes(yintercept=0), size=1, linetype="dashed"
)+facet_wrap(~temp.var)

plmend_lr_plot<-ggplot(plmend, aes(x=load, y=resid_f, color=temp.avg))
plmend_pr_plot+geom_point(shape=1
)+geom_hline(aes(yintercept=0)
)+facet_wrap(~temp.var)





#-----------------------------------

#Trying Joel's suggestion of analyzing wasp survival by running a linear model with the response variable a logit transformed percentage of 
##survival

tv.para$lgt_surv<-logit(tv.para$tot.surv)

lmlogit_mod<-lme(lgt_surv ~ temp.avg*temp.var*load,
                 random=~1|bug.id,
                 data=tv.para,
                 na.action = na.omit,
                 method="ML")

anova(lmlogit_mod)
summary(lmlogit_mod)


#make a model without load, see which is better

lmlogit_mod_red<-lme(lgt_surv ~ temp.avg*temp.var,
                  random=~1|bug.id,
                  data=tv.para,
                  na.action = na.omit,
                  method="ML")

anova(lmlogit_mod, lmlogit_mod_red)


#reduced model without non significant terms
lmlogit_mod_red2<-lme(lgt_surv ~ temp.avg*temp.var + temp.avg:load,
                     random=~1|bug.id,
                     data=tv.para,
                     na.action = na.omit,
                     method="ML")


anova(lmlogit_mod, lmlogit_mod_red, lmlogit_mod_red2)


#plotting model results

lgt<-tv.para

lgt$pred<-predict(lmlogit_mod)
lgt$resid<-residuals(lmlogit_mod)
lgt$pred_f<-predict(lmlogit_mod, level=0)


lmlgt_pr_fit<-ggplot(lgt, aes(x=pred, y=resid, color=temp.avg))
lmlgt_pr_fit+geom_point(shape=1, size=5
)+geom_hline(aes(yintercept=0), 
             color="black", linetype="dashed"
)+facet_wrap(~temp.var)

lmlgt_lr_fit<-ggplot(lgt, aes(x=load, y=resid, color=temp.avg))
lmlgt_lr_fit+geom_jitter(shape=1, size=5
)+geom_hline(aes(yintercept=0), 
             color="black", linetype="dashed"
)+facet_wrap(~temp.var)


lmlgt_modfit<-ggplot(lgt, aes(x=load, y=lgt_surv, color=temp.avg))
lmlgt_modfit+geom_point(size=4
)+geom_line(aes(y=pred_f, group=temp.avg),
            size=2
)+facet_wrap(~temp.var)



#Plot wasp survival to ecl by load, not transformed

wsurv_ld_plot<-ggplot(tv.para, aes(x=load, y=tot.surv, color=temp.avg))
wsurv_ld_plot+geom_point(size=4
)+geom_smooth(method="lm"
)+facet_wrap(~temp.var)


#plot wasp survival to em by load, not transformed
tv.para$perc_em<-tv.para$num.em / tv.para$load

percem_ld_plot<-ggplot(tv.para, aes(x=load, y=perc_em, color=temp.avg))
percem_ld_plot+geom_point(size=4
)+geom_smooth(method="lm"
)+facet_wrap(~temp.var)


#Linear model of logit transformed perc_em data

tv.para$lgt_emsurv<-logit(tv.para$perc_em)

lmlgt_em_mod<-lme(lgt_emsurv ~ temp.avg*temp.var*load,
                  random=~1|bug.id,
                  data=tv.para,
                  na.action = na.omit,
                  method="ML")

anova(lmlgt_em_mod)


lmlgt_em_nullmod<-lme(lgt_emsurv ~ temp.avg*temp.var,
                      random=~1|bug.id,
                      data=tv.para,
                      na.action = na.omit,
                      method="ML")

anova(lmlgt_em_nullmod)


anova(lmlgt_em_mod, lmlgt_em_nullmod)




#-------------------------------------------------------------------------------------------------------------------------------


#looking at effects of load, by separating data into 25 and 28 mean temps to try and get rid of 30+/-10 effects

#Making a column for total died (load-num.ecl)

tv$tot.died<-tv$load-tv$num.ecl

tv$resc.ld<-rescale(tv$load,to=c(0,1))


tv.low<-subset(tv.para, temp.avg!=30)
tv.high<-subset(tv.para, temp.avg==30)

wsload.mod1<-glmer(cbind(num.ecl, tot.died)~temp.avg*temp.var*resc.ld + (1|bug.id),
                 family = binomial,
                 data=tv.low,
                 na.action=na.omit,
                 control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

anova(wsload.mod1)
summary(wsload.mod1)


wsload.mod2<-glmer(cbind(num.ecl, tot.died)~temp.var*resc.ld + (1|bug.id),
                   family = binomial,
                   data=tv.high,
                   na.action=na.omit,
                   control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

anova(wsload.mod2)
summary(wsload.mod2)



wemsload.mod1<-glmer(cbind(num.em, num.unem)~temp.avg*temp.var*resc.ld + (1|bug.id),
                   family = binomial,
                   data=tv.low,
                   na.action=na.omit,
                   control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

anova(wemsload.mod1)
summary(wemsload.mod1)




wemsload.mod2<-glmer(cbind(num.em, num.unem)~temp.var*resc.ld + (1|bug.id),
                     family = binomial,
                     data=tv.high,
                     na.action=na.omit,
                     control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

anova(wemsload.mod2)
summary(wemsload.mod2)



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
                data=tv.para,
                na.action=na.omit,
                method="ML")

anova(wdtot.mod1)
summary(wdtot.mod1)


#trying a non-mixed effects model

wdtot.mod2<-lm(waspdev.tot~temp.avg*temp.var*resc.ld,
               data=tv.para,
               na.action = na.omit)

anova(wdtot.mod2)

#-----------------------------------


#MODELLING WASP TOTAL SURVIVAL--GLMM


#Making a column for total died (load-num.ecl)

tv.para$tot.died<-tv.para$load-tv.para$num.ecl


#rescaling load, tot.died and num.ecl

tv.para$resc.ld<-rescale(tv.para$load,to=c(0,1))

#tv.para$man.resc.ld<-tv.para$load/321
#tv.para$temp.var.num<-as.numeric(tv.para$temp.var)
#tv.para.no10<-subset(tv.para,temp.avg!="30" & temp.var!="10")


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




#Making a glmer (binomial) model of wasp total survival
##temp.avg==factor, temp.var==numeric, resc.load==numeric, mongos treated as NAs until I determine a better way to deal with them
###won't run without rescaling load
###as of 2.27.18, won't run even when load is rescaled. Don't know what's going on

wtots.mod2<-glmer(cbind(tot.died,num.ecl)~temp.avg*temp.var*resc.ld+(1|bug.id),
                  family=binomial,
                  data=tv.para,
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
