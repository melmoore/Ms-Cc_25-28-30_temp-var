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
library(ggplot2)
library(extrafont)
library(MuMIn)


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


#set plot theme
theme_set(theme_classic())

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


#make a label for the facet wrap panels
fnames <- c("control" = "NP", "para" = "P", "0" = "0", "10" = "10")



#plot model residuals against age, color by mean temperature, facet by parasitization treatment and fluctuation
md_gammod_ra<-ggplot(tvor_mass, aes(x=age, y=resid, color=temp.avg))
md_gammod_ra + geom_point(size=4, shape=1
) + geom_hline(aes(yintercept=0),
             color="black",
             size=1.5, linetype="dashed"
) + scale_color_manual(values=c("#009E73","#E69F00","#000000"),name=c("Avg. Temp. [C]"),
                       breaks=c("25","28","30"),labels=c("25","28","30"),
                       guide=guide_legend(keywidth=3)   
) + labs(x="Age [Days]", y="GAMM Model Residuals"
) + facet_wrap(treatment~temp.var, labeller = as_labeller(fnames)
) + theme(text = element_text(family=("Cambria")),
          strip.background = element_rect(colour="black",linetype = "solid",fill="white",
                                          size = 1),
          strip.text = element_text(size=18),
          axis.line.x=element_line(colour = 'black', size = 1),
          axis.line.y=element_line(colour = 'black', size = 1),
          axis.ticks = element_line(colour = 'black', size = 1),
          axis.ticks.length = unit(2, "mm"),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          legend.background = element_rect(color="black",linetype="solid"),
          legend.text = element_text(size=16),
          legend.title = element_text(size=16))



md_gammod_fit<-ggplot(tvor_mass, aes(x=age, y=log_mss, group=interaction(bug.id, temp.avg), color=temp.avg))
md_gammod_fit+geom_point(size=3, shape=1
) + geom_line(aes(y=pred, group=interaction(bug.id, temp.avg))
) + scale_color_manual(values=c("#009E73","#E69F00","#000000"),name=c("Avg. Temp. [C]"),
                       breaks=c("25","28","30"),labels=c("25","28","30"),
                       guide=guide_legend(keywidth=3)   
) + labs(x="Age [Days]", y="Log(Mass [mg])"
) + facet_wrap(treatment~temp.var, labeller = as_labeller(fnames)
) + theme(text = element_text(family=("Cambria")),
          strip.background = element_rect(colour="black",linetype = "solid",fill="white",
                                          size = 1),
          strip.text = element_text(size=18),
          axis.line.x=element_line(colour = 'black', size = 1),
          axis.line.y=element_line(colour = 'black', size = 1),
          axis.ticks = element_line(colour = 'black', size = 1),
          axis.ticks.length = unit(2, "mm"),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          legend.background = element_rect(color="black",linetype="solid"),
          legend.text = element_text(size=16),
          legend.title = element_text(size=16))






#using geom_smooth instead of geom_line to try and get an average line over temp.avg and not by individual?
#not sure if this is appropriate.
md_gammod_fit<-ggplot(tvor_mass, aes(x=age, y=log_mss, color=temp.avg))
md_gammod_fit + geom_point(size=3, shape=1
) + geom_smooth(aes(y=pred, group=temp.avg)
) + scale_color_manual(values=c("#009E73","#E69F00","#000000"),name=c("Avg. Temp. [C]"),
                       breaks=c("25","28","30"),labels=c("25","28","30"),
                       guide=guide_legend(keywidth=3)   
) + facet_wrap(treatment~temp.var, labeller = as_labeller(fnames)
) + labs(x="Age [Days]", y="Log(Mass [mg])"
) + theme(text = element_text(family=("Cambria")),
          strip.background = element_rect(colour="black",linetype = "solid",fill="white",
                                          size = 1),
          strip.text = element_text(size=18),
          axis.line.x=element_line(colour = 'black', size = 1),
          axis.line.y=element_line(colour = 'black', size = 1),
          axis.ticks = element_line(colour = 'black', size = 1),
          axis.ticks.length = unit(2, "mm"),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          legend.background = element_rect(color="black",linetype="solid"),
          legend.text = element_text(size=16),
          legend.title = element_text(size=16))




#--------------------------


#para treatment with load as a predictor

#subset to only columns in model, and remove rows with NAs (so that predicted and fitted values can be
#added to the dataframe easily)
p_mass<-select(tvor_lngp, bug.id, temp.avg, temp.var, load, log_mss, age)
p_mass<-na.omit(p_mass)

#convert bug.id to factor so it functions properly as a random effect
p_mass$bug.id<-as.factor(p_mass$bug.id)


#run a GAMM with age and load as separate smooths--does not have an interaction with temp this way (2 2D surfaces, instead
## of one 3D surface)
gam_pml_nointmod<-gam(log_mss ~ s(age, by=interaction(temp.avg, temp.var, bs="ts")) 
                      + s(load, by=interaction(temp.avg, temp.var, k=10, bs="ts")) + s(bug.id, bs="re") 
                      + temp.avg * temp.var, method="ML", data=p_mass, na.action = na.omit)

anova(gam_pml_nointmod)
summary(gam_pml_nointmod)

gam.check(gam_pml_nointmod)

#Model with age and load in separate smooths
#make columns with predicted and residual values for plotting
p_mass$pred_ni<-predict(gam_pml_nointmod, level=0)
p_mass$resid_ni<-residuals(gam_pml_nointmod, level=0)

#residuals against age
#residuals much better than model with interaction between age and load
pmlni_gammod_ra<-ggplot(p_mass, aes(x=age, y=resid_ni, color=temp.avg))
pmlni_gammod_ra+geom_point(size=4, shape=1
)+geom_hline(aes(yintercept=0),
             color="black",
             size=1.5, linetype="dashed"
)+facet_wrap(~temp.var)


#residuals against load
#look ok except for all the WOWEs with 0 load
pmlni_gammod_rl<-ggplot(p_mass, aes(x=load, y=resid_ni, color=temp.avg))
pmlni_gammod_rl+geom_point(size=4, shape=1
)+geom_hline(aes(yintercept=0),
             color="black",
             size=1.5, linetype="dashed"
) + scale_color_manual(values=c("#009E73","#E69F00","#000000"),name=c("Avg. Temp. [C]"),
                       breaks=c("25","28","30"),labels=c("25","28","30"),
                       guide=guide_legend(keywidth=3)   
) + labs(x="Total Load", y="GAMM Model Residuals"
) + facet_wrap(~temp.var, labeller = as_labeller(fnames)
) + theme(text = element_text(family=("Cambria")),
          strip.background = element_rect(colour="black",linetype = "solid",fill="white",
                                          size = 1),
          strip.text = element_text(size=18),
          axis.line.x=element_line(colour = 'black', size = 1),
          axis.line.y=element_line(colour = 'black', size = 1),
          axis.ticks = element_line(colour = 'black', size = 1),
          axis.ticks.length = unit(2, "mm"),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          legend.background = element_rect(color="black",linetype="solid"),
          legend.text = element_text(size=16),
          legend.title = element_text(size=16))



#model fit, colored by load
pmlni_gammod_fit<-ggplot(p_mass, aes(x=age, y=log_mss, color=load))
pmlni_gammod_fit+geom_point(size=4, shape=1
) + geom_line(aes(y=pred_ni, group=interaction(temp.avg, bug.id)),
            size=1
) + scale_color_viridis(begin=.1, end=1,
                      option = "viridis",
                      name="Total Load"
) + facet_wrap(temp.avg~temp.var
) + labs(x="Age [Days]", y="Log(Mass [mg]"
) + theme(text = element_text(family=("Cambria")),
          strip.background = element_rect(colour="black",linetype = "solid",fill="white",
                                          size = 1),
          strip.text = element_text(size=18),
          axis.line.x=element_line(colour = 'black', size = 1),
          axis.line.y=element_line(colour = 'black', size = 1),
          axis.ticks = element_line(colour = 'black', size = 1),
          axis.ticks.length = unit(2, "mm"),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          legend.background = element_rect(color="black",linetype="solid"),
          legend.text = element_text(size=16),
          legend.title = element_text(size=16),
          legend.position = c(0.95, 0.8))



#--------------------

#ANALYSES OF WASP SURVIVAL TO EMERGENCE AND ECLOSION

#Remove the 30+/-10 treatment, as no wasps survived, the mechanism is not the same as the other treatments

#create combo column with temp.avg and temp.var
tvor_p <- unite(tvor_p, "tatv", temp.avg, temp.var, remove = FALSE)

tvor_nw <- subset(tvor_p, tatv!="30_10")

#make temp var a character instead of factor
tvor_nw$temp.var <- as.character(tvor_nw$temp.var)


#rescale load 
tvor_nw$resc_ld <- rescale(tvor_nw$load, to=c(0,1))


#Full GLMER model, binomial distribution. Number survived to eclosion (success) vs number died (failure) as response variable
#Mean temperature, temperature fluctuation and rescaled load as fixed effects. Random intercept of individual
ws_nowowe_re_mod1 <- glmer(cbind(num.ecl, tot.died) ~ temp.avg * temp.var * resc_ld + (1|bug.id),
                           family=binomial,
                           data=tvor_nw,
                           na.action=na.omit,
                           control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

anova(ws_nowowe_re_mod1)
summary(ws_nowowe_re_mod1)





#model selection using dredge() 

#dredge requires dataframe with no NAs--subsetting to only columns in the model
tvor_wsmod <- tvor_nw %>% select(bug.id, temp.avg, temp.var, resc_ld, num.ecl, tot.died)
tvor_wsmod <- drop_na(tvor_wsmod)

#model with data frame with no NAs, na.action set to na.fail
ws_nowowe_re_mod <- glmer(cbind(num.ecl, tot.died) ~ temp.avg * temp.var * resc_ld + (1|bug.id),
                          family=binomial,
                          data=tvor_wsmod,
                          na.action=na.fail,
                          control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))


wsmod_dredge <- dredge(ws_nowowe_re_mod)
wsmod_dredge

#best model is one without the interaction between tv and load, and without 3 way interaction


#reduced model from dredge:
ws_nowowe_re_mod_rd <- glmer(cbind(num.ecl, tot.died) ~ temp.avg + temp.var + resc_ld + temp.avg:temp.var +
                           temp.avg:resc_ld + (1|bug.id),
                           family=binomial,
                           data=tvor_nw,
                           na.action=na.omit,
                           control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

anova(ws_nowowe_re_mod_rd)



#models without fixed effects of interest: removing main effect and interactions involving main effect

ws_nowowe_re_mod_noload <- glmer(cbind(num.ecl, tot.died) ~ temp.avg + temp.var + temp.avg:temp.var + (1|bug.id),
                             family=binomial,
                             data=tvor_nw,
                             na.action=na.omit,
                             control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))


ws_nowowe_re_mod_nota <- glmer(cbind(num.ecl, tot.died) ~ temp.var + resc_ld + (1|bug.id),
                               family=binomial,
                               data=tvor_nw,
                               na.action=na.omit,
                               control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))


ws_nowowe_re_mod_notv <- glmer(cbind(num.ecl, tot.died) ~ temp.avg + resc_ld + temp.avg:resc_ld + (1|bug.id),
                               family=binomial,
                               data=tvor_nw,
                               na.action=na.omit,
                               control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))



ws_nowowe_re_mod_notatv <- glmer(cbind(num.ecl, tot.died) ~ temp.avg + temp.var + resc_ld + temp.avg:resc_ld + (1|bug.id),
                             family=binomial,
                             data=tvor_nw,
                             na.action=na.omit,
                             control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))


ws_nowowe_re_mod_notald <- glmer(cbind(num.ecl, tot.died) ~ temp.avg + temp.var + resc_ld + temp.avg:temp.var + (1|bug.id),
                             family=binomial,
                             data=tvor_nw,
                             na.action=na.omit,
                             control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))




anova(ws_nowowe_re_mod_rd, ws_nowowe_re_mod_noload, test="Chisq")
anova(ws_nowowe_re_mod_rd, ws_nowowe_re_mod_nota, test="Chisq")
anova(ws_nowowe_re_mod_rd, ws_nowowe_re_mod_notv, test="Chisq")
anova(ws_nowowe_re_mod_rd, ws_nowowe_re_mod_notatv, test="Chisq")
anova(ws_nowowe_re_mod_rd, ws_nowowe_re_mod_notald, test="Chisq")




#----------------------

#Analysis of wasp survival to emergence--individuals from 30.10 treatment removed

#Full GLMER model, binomial distribution. Number survived to emergence (success) vs number died (failure) as response variable
#Mean temperature, temperature fluctuation and rescaled load as fixed effects. Random intercept of individual
wsem_nw_re_mod1 <- glmer(cbind(num.em, num.unem) ~ temp.avg * temp.var * resc_ld + (1|bug.id),
                           family=binomial,
                           data=tvor_nw,
                           na.action=na.omit,
                           control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

anova(wsem_nw_re_mod1)
summary(wsem_nw_re_mod1)



#model selection using dredge()

#dredge requires dataframe with no NAs--subsetting to only columns in the model
tvor_wsem <- tvor_nw %>% select(bug.id, temp.avg, temp.var, resc_ld, num.em, num.unem)
tvor_wsem <- drop_na(tvor_wsem)


#model with dataset without NAs, na.action set to na.fail
wsem_nw_re_mod <- glmer(cbind(num.em, num.unem) ~ temp.avg * temp.var * resc_ld + (1|bug.id),
                         family=binomial,
                         data=tvor_wsem,
                         na.action=na.fail,
                         control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))


wsem_dredge <- dredge(wsem_nw_re_mod)
wsem_dredge


#reduced model from dredge
wsem_nw_mod_rd <- glmer(cbind(num.em, num.unem) ~ temp.avg + temp.var + resc_ld + temp.avg:temp.var 
                        + (1|bug.id),
                         family=binomial,
                         data=tvor_nw,
                         na.action=na.omit,
                         control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

anova(wsem_nw_mod_rd)
summary(wsem_nw_mod_rd)



#determing effects of fixed effects: dropping main effects and interactions involving them, comparing to dredge model
wsem_nw_mod_nold <- glmer(cbind(num.em, num.unem) ~ temp.avg + temp.var + temp.avg:temp.var + (1|bug.id),
                          family=binomial,
                          data=tvor_nw,
                          na.action=na.omit,
                          control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))


wsem_nw_mod_nota <- glmer(cbind(num.em, num.unem) ~ temp.var + resc_ld + (1|bug.id),
                          family=binomial,
                          data=tvor_nw,
                          na.action=na.omit,
                          control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))



wsem_nw_mod_notv <- glmer(cbind(num.em, num.unem) ~ temp.avg + resc_ld + (1|bug.id),
                          family=binomial,
                          data=tvor_nw,
                          na.action=na.omit,
                          control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))



wsem_nw_mod_notatv <- glmer(cbind(num.em, num.unem) ~ temp.avg + temp.var + resc_ld + (1|bug.id),
                            family=binomial,
                            data=tvor_nw,
                            na.action=na.omit,
                            control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))



#compare to dredge model
anova(wsem_nw_mod_rd, wsem_nw_mod_nold, test="Chisq")
anova(wsem_nw_mod_rd, wsem_nw_mod_nota, test="Chisq")
anova(wsem_nw_mod_rd, wsem_nw_mod_notv, test="Chisq")
anova(wsem_nw_mod_rd, wsem_nw_mod_notatv, test="Chisq")



#-----------------------

#Analysis of wasp development time

#Use the dataset without 30.10 treatment group

#Linear model of development time to emergence from the host; time from oviposition to emergence (in days)
#as response variable, mean temperature (factor), fluctuation (character) and load (numeric) as fixed effects.  
wdev_em_mod1 <- lm(ttem.w ~ temp.avg * temp.var * load,
                    data = tvor_nw,
                    na.action = na.omit)

anova(wdev_em_mod1)



#model selection using dredge() 

#dredge requires dataframe with no NAs--subsetting to only columns in the model
tvor_wdevem_mod <- tvor_nw %>% select(bug.id, temp.avg, temp.var, load, ttem.w)
tvor_wdevem_mod <- drop_na(tvor_wdevem_mod)

#model with data frame with no NAs, na.action set to na.fail
wdev_em_mod <- lm(ttem.w ~ temp.avg * temp.var * load,
                   data = tvor_wdevem_mod,
                   na.action = na.fail)

#dredged models
wdev_em_dredge <- dredge(wdev_em_mod)
wdev_em_dredge



#Reduced model from dredge
wdev_em_mod2 <- lm(ttem.w ~ temp.avg + temp.var + load + temp.avg:temp.var + temp.avg:load,
                   data = tvor_nw,
                   na.action = na.omit)

anova(wdev_em_mod2)




#looking at model fit
dt_em <- tvor_nw

#adding predicted and residual values
dt_em$pred <- predict(wdev_em_mod1)
dt_em$resid <- residuals(wdev_em_mod1)


#plot predicted values against residuals
wdtem_rf_plot <- ggplot(dt_em, aes(x=pred, y=resid, group=temp.avg, color=temp.avg))
wdtem_rf_plot + geom_point(shape=1, size=5
) + geom_hline(aes(yintercept=0), color="black", linetype="dashed", size=2
) + facet_wrap(~temp.var)


#plot residual values against load
wdtem_rld_plot <- ggplot(dt_em, aes(x=load, y=resid, group=temp.avg, color=temp.avg))
wdtem_rld_plot + geom_point(shape=1, size=5
) + geom_hline(aes(yintercept=0), color="black", linetype="dashed", size=2
) + facet_wrap(~temp.var)


#plot fitted values against raw data
wdtem_fld_plot <- ggplot(dt_em, aes(x=load, y=ttem.w, group=temp.avg, color=temp.avg))
wdtem_fld_plot + geom_point(shape=1, size=5
) + geom_line(aes(y=pred, group=temp.avg, color=temp.avg),
              size=2
) + facet_wrap(~temp.var)




#Linear model of development time to eclosion from the host; time from oviposition to eclosion (in days)
#as response variable, mean temperature (factor), fluctuation (character) and load (numeric) as fixed effects.  
wdev_ecl_mod1 <- lm(ttecl ~ temp.avg * temp.var * load,
                   data = tvor_nw,
                   na.action = na.omit)

anova(wdev_ecl_mod1)


#dredge requires dataframe with no NAs--subsetting to only columns in the model
tvor_wdevecl_mod <- tvor_nw %>% select(bug.id, temp.avg, temp.var, load, ttecl)
tvor_wdevecl_mod <- drop_na(tvor_wdevecl_mod)

#model with data frame with no NAs, na.action set to na.fail
wdev_ecl_mod <- lm(ttecl ~ temp.avg * temp.var * load,
                  data = tvor_wdevecl_mod,
                  na.action = na.fail)


wdev_ecl_dredge <- dredge(wdev_ecl_mod)
wdev_ecl_dredge


#reduced model from dredge:
wdev_ecl_mod2 <- lm(ttecl ~ temp.avg + temp.var + load + temp.avg:temp.var,
                    data = tvor_nw,
                    na.action = na.omit)


anova(wdev_ecl_mod2)



#looking at model fit
dt_ecl <- drop_na(tvor_nw, ttecl)

#adding predicted and residual values
dt_ecl$pred <- predict(wdev_ecl_mod1)
dt_ecl$resid <- residuals(wdev_ecl_mod1)


#plot predicted values against residuals
wdtecl_rf_plot <- ggplot(dt_ecl, aes(x=pred, y=resid, group=temp.avg, color=temp.avg))
wdtecl_rf_plot + geom_point(shape=1, size=5
) + geom_hline(aes(yintercept=0), color="black", linetype="dashed", size=2
) + facet_wrap(~temp.var)


#plot residual values against load
wdtecl_rld_plot <- ggplot(dt_ecl, aes(x=load, y=resid, group=temp.avg, color=temp.avg))
wdtecl_rld_plot + geom_point(shape=1, size=5
) + geom_hline(aes(yintercept=0), color="black", linetype="dashed", size=2
) + facet_wrap(~temp.var)



#plot fitted values against raw data
wdtecl_fld_plot <- ggplot(dt_ecl, aes(x=load, y=ttecl, group=temp.avg, color=temp.avg))
wdtecl_fld_plot + geom_point(shape=1, size=5
) + geom_line(aes(y=pred, group=temp.avg, color=temp.avg),
              size=2
) + facet_wrap(~temp.var)




#----------------------

#Analysis of final mass by mean temperature, fluctuation and parasitization treatment

#convert final mass to log scale
tvor$log_mssend <- log(tvor$mass.end)

#Linear mixed effects model, with log of final mass as response, mean temperature, fluctuation, parasitization
#treatment and their interactions as fixed effects. Individual included as a random slope. 
finmass_mod <- lme(log_mssend ~ temp.avg * temp.var * treatment,
                   random = ~ 1|bug.id,
                   data = tvor,
                   method = "ML",
                   na.action = na.omit)

anova(finmass_mod)
summary(finmass_mod)



#---------------------------

#Analysis of final mass for parasitized caterpillars by mean temperature, fluctuation and load

#Convert final mass to log scale for parasitized caterpillars
tvor_p$log_mssend <- log(tvor_p$mass.end)


#Remove 30.10 treatment, as they have no load and it therefor cannot affect final mass
tvor_p_nw <- subset(tvor_p, end.class!="cull")

#Convert fluctuation to character, instead of factor
tvor_p_nw$temp.var <- as.character(tvor_p_nw$temp.var)


#Linear model, with log of final mass as response, mean temperature, fluctuation, load,
#and their interactions as fixed effects. 
finmass_p_mod <- lm(log_mssend ~ temp.avg * temp.var * load,
                   data = tvor_p_nw,
                   na.action = na.omit)

anova(finmass_p_mod)
summary(finmass_p_mod)





