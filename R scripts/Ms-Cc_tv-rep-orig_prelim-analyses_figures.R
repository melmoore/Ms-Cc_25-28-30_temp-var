#Ms Cc temp var experiment, prelim figure and analyses comparing the original data to the replication expt

#load libraries
library(scales)
library(Rmisc)
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
library(cowplot)



#load data

tvor <- read_csv("data files/Ms-Cc_tv-orig-rep_comb_cl.csv", 
                 col_types = cols(temp.avg = col_factor(levels = c("25", "28", "30")), 
                                  temp.var = col_factor(levels = c("0", "5", "10")), 
                                  treatment = col_factor(levels = c("control", "para"))))

tvor_lng <- read_csv("data files/Ms-Cc_tv-orig-rep_comb_lng.csv", 
                     col_types = cols(temp.avg = col_factor(levels = c("25",  "28", "30")),
                                      temp.var = col_factor(levels = c("0", "5", "10")), 
                                      treatment = col_factor(levels = c("control", "para"))))

#take log of caterpillar mass
tvor_lng$log_mss <- log(tvor_lng$mass)

#create tot.died column
tvor$tot.died <- tvor$load - tvor$num.ecl


#----------------------------


#raw data of wasp survival for the original and replication experiments

#subsetting to only 30C treatments to compare wasp survival between orig and repl (remove +/-5)
tvor_30 <- subset(tvor, temp.avg==30 & temp.var!=5)

#plot number emerged by load, color by expt
numem_plot <- ggplot(tvor_30, aes(x=load, y=num.em, color=expt))
numem_plot + geom_point(
) + geom_smooth(method = "lm")



#plot number eclosed by load, color by expt
numecl_plot <- ggplot(tvor_30, aes(x=load, y=num.ecl, color=expt))
numecl_plot + geom_point(
) + geom_smooth(method = "lm")



#plot proportion emerged by load, color by expt--geom_smooth is doing something weird, don't know why
psem_plot <- ggplot(tvor_30, aes(x=load, y=ps.em, group=expt, color=expt))
psem_plot + geom_point(
) + geom_smooth(method = "lm")


#plot proportion eclosedd by load, color by expt--geom_smooth is doing something weird, don't know why
psecl_plot <- ggplot(tvor_30, aes(x=load, y=ps.ecl, group=expt, color=expt))
psecl_plot + geom_point(
)+geom_smooth(method="lm")


#subset to only const and para treatments
tvor_p0 <- subset(tvor_30, temp.var==0 & treatment=="para")


#boxplot of ps.em
psem_boxplot <- ggplot(tvor_p0, aes(x=expt, y=ps.em))
psem_boxplot + geom_boxplot()


#boxplot of ps.ecl
psecl_boxplot <- ggplot(tvor_p0, aes(x=expt, y=ps.ecl))
psecl_boxplot + geom_boxplot()


#----------------------

#quick and dirty analysis of wasp survival between 30C orig and repl data

expt_wtots_mod1 <- glm(cbind(num.ecl, tot.died) ~ expt * load,
                       family = quasibinomial,
                       data = tvor_p0,
                       na.action = na.omit)
anova(expt_wtots_mod1)
summary(expt_wtots_mod1)


tvor_p0$resc.ld <- rescale(tvor_p0$load, to=c(0,1))

expt_wtots_modre <- glmer(cbind(num.ecl, tot.died) ~ expt * load + (1|bug.id),
                          family = binomial,
                          data = tvor_p0,
                          na.action = na.omit)
anova(expt_wtots_modre)
summary(expt_wtots_modre)

#---------------------------


#mean wasp survival for original and replication experiments--comparing orig and repl data
theme_set(theme_classic())

#subset to only parasitized individuals, and removing +/-5 treatment
tvor_p <- subset(tvor, treatment=="para" & temp.var!=5)


#proportion surviving to emergence

psem_sum <- summarySE(tvor_p, measurevar = "ps.em",
                      groupvars = c("temp.avg", "temp.var", "expt"),
                      na.rm = TRUE)
psem_sum


#making temp.avg numeric instead of a factor
psem_sum$temp.avg <- as.numeric(psem_sum$temp.avg)
psem_sum$temp.avg <- ifelse(psem_sum$temp.avg==1, 25,
                            ifelse(psem_sum$temp.avg==2, 28, 30))



#plot of mn ps.em, with numeric temp.avg on the x axis, psem on the y axis, grouped by temp.var
mn_psem_plot <- ggplot(psem_sum, aes(x=temp.avg, y=ps.em, group=interaction(temp.var, expt),
                                     color=temp.var))
mn_psem_plot + geom_point(aes(shape=expt),
                        size=5
) + geom_line(aes(linetype=temp.var),
              size = 1.2
) + geom_errorbar(aes(ymin = ps.em-se, ymax = ps.em+se)) 





#proportion surviving to eclosion

psecl_sum <- summarySE(tvor_p, measurevar = "ps.ecl",
                      groupvars = c("temp.avg", "temp.var", "expt"),
                      na.rm = TRUE)
psecl_sum

#making temp.avg numeric instead of a factor
psecl_sum$temp.avg <- as.numeric(psecl_sum$temp.avg)
psecl_sum$temp.avg <- ifelse(psecl_sum$temp.avg==1, 25,
                             ifelse(psecl_sum$temp.avg==2, 28, 30))


#plot of mn ps.ecl, with numeric temp.avg on the x axis, psecl on the y axis, grouped by temp.var
mn_psecl_plot <- ggplot(psecl_sum, aes(x=temp.avg, y=ps.ecl, group=interaction(temp.var, expt),
                                     color=temp.var))
mn_psecl_plot + geom_point(aes(shape=expt),
                          size=5
) + geom_line(aes(linetype=temp.var),
              size = 1.2
) + geom_errorbar(aes(ymin = ps.ecl-se, ymax = ps.ecl+se)) 




#subset out the orig 30 treatment, so I can get plots of the repl data with the orig other 2 temp avg treatments

#prop em
psem_sum2 <- psem_sum[-c(5,7),]

#plot of mn ps.em, with numeric temp.avg on the x axis, psem on the y axis, grouped by temp.var
#--ONLY repl data for 30C temp.avg
mn_psem_plot2 <- ggplot(psem_sum2, aes(x=temp.avg, y=ps.em, group=temp.var,
                                     color=temp.var))
mn_psem_plot2 <- mn_psem_plot2 + geom_point(size=6
) + geom_line(size = 2
) + geom_errorbar(aes(ymin = ps.em-se, ymax = ps.em+se),
                  width=.5, size=1.2
) + scale_color_manual(values=c("#56B4E9","#D55E00"),name=c("Fluctuation [C]"),
                       breaks=c("0","10"),labels=c("0","10"),
                       guide=guide_legend(keywidth = 2.5)
) + scale_x_continuous(limits=c(24.5,30.5),
                       breaks = c(25, 28, 30)
) + scale_y_continuous(limits = c(0, 0.9),
                       breaks = c(0, 0.2, 0.4, 0.6, 0.8)
) + labs(x="Mean Temperature [C]", y="% Emergence"
) + theme(axis.line.x=element_line(colour = 'black', size = 1),
          axis.line.y=element_line(colour = 'black', size = 1),
          axis.ticks = element_line(colour = 'black', size = 1),
          axis.ticks.length = unit(2, "mm"),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          legend.background = element_rect(color="black",linetype="solid"),
          legend.position = "none")





#prop ecl
psecl_sum2 <- psecl_sum[-c(5,7),]


#plot of mn ps.ecl, with numeric temp.avg on the x axis, psecl on the y axis, grouped by temp.var
#--ONLY repl data for 30C temp.avg treatments
mn_psecl_plot2 <- ggplot(psecl_sum2, aes(x=temp.avg, y=ps.ecl, group=temp.var,
                                       color=temp.var))
mn_psecl_plot2 <- mn_psecl_plot2 + geom_point(size=6
) + geom_line(size = 2
) + geom_errorbar(aes(ymin = ps.ecl-se, ymax = ps.ecl+se),
                  width=.5, size=1.2
) + scale_color_manual(values=c("#56B4E9","#D55E00"),name=c("Fluctuation [C]"),
                         breaks=c("0","10"),labels=c("0","10"),
                         guide=guide_legend(keywidth = 2.5)
) + scale_x_continuous(limits=c(24.5,30.5),
                         breaks = c(25, 28, 30)
) + scale_y_continuous(limits = c(0, 0.75),
                       breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
) + labs(x="Mean Temperature [C]", y="% Eclosion"
) + theme(axis.line.x=element_line(colour = 'black', size = 1),
        axis.line.y=element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = 'black', size = 1),
        axis.ticks.length = unit(2, "mm"),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.background = element_rect(color="black",linetype="solid"),
        legend.position = c(.8, .85))

mn_psecl_plot2



#combine into one figure using cowplot
surv_fig_repl <- plot_grid(mn_psem_plot2, mn_psecl_plot2, labels=c("A", "B"), align = "h")
surv_fig_repl






#plot of mn ps.em, with numeric temp.avg on the x axis, psem on the y axis, grouped by temp.var
#--data from orig expt added as disconnected point

orig <- psem_sum[c(5,7),]

mn_psem_plot3 <- ggplot(psem_sum2, aes(x=temp.avg, y=ps.em, group=temp.var))
mn_psem_plot3 <- mn_psem_plot2 + geom_point(size=6, aes(color=temp.var)
) + geom_line(size = 2, aes(color=temp.var)
) + geom_errorbar(aes(ymin = ps.em-se, ymax = ps.em+se),
                  width=.5, size=1.2
) + geom_point(data=orig, aes(x=temp.avg, y=ps.em),
               color="black", size=6, shape=17
) + geom_errorbar(data=orig, aes(ymin=ps.em - se, ymax=ps.em +se),
                  width=.5, size=1.2, color="black"
) + scale_color_manual(values=c("#56B4E9","#D55E00"),name=c("Fluctuation [C]"),
                       breaks=c("0","10"),labels=c("0","10"),
                       guide=guide_legend(keywidth = 2.5)
) + scale_x_continuous(limits=c(24.5,30.5),
                       breaks = c(25, 28, 30)
) + scale_y_continuous(limits = c(0, 0.9),
                       breaks = c(0, 0.2, 0.4, 0.6, 0.8)
) + labs(x="Mean Temperature [C]", y="% Emergence"
) + theme(axis.line.x=element_line(colour = 'black', size = 1),
          axis.line.y=element_line(colour = 'black', size = 1),
          axis.ticks = element_line(colour = 'black', size = 1),
          axis.ticks.length = unit(2, "mm"),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          legend.background = element_rect(color="black",linetype="solid"),
          legend.position = "none")

mn_psem_plot3





#plot of mn ps.ecl, with numeric temp.avg on the x axis, psecl on the y axis, grouped by temp.var
##--data from orig expt added as disconnected point

orig_ecl <- psecl_sum[c(5,7),]

mn_psecl_plot3 <- ggplot(psecl_sum2, aes(x=temp.avg, y=ps.ecl, group=temp.var))
mn_psecl_plot3 <- mn_psecl_plot2 + geom_point(size=6, aes(color=temp.var)
) + geom_line(size = 2, aes(color=temp.var)
) + geom_errorbar(aes(ymin = ps.ecl-se, ymax = ps.ecl+se),
                  width=.5, size=1.2
) + geom_point(data=orig_ecl, aes(x=temp.avg, y=ps.ecl),
               color="black", shape=17, size = 6
) + geom_errorbar(data=orig_ecl, aes(ymin = ps.ecl-se, ymax = ps.ecl+se),
                  width=.5, size=1.2, color="black"
) + scale_color_manual(values=c("#56B4E9","#D55E00"),name=c("Fluctuation [C]"),
                       breaks=c("0","10"),labels=c("0","10"),
                       guide=guide_legend(keywidth = 2.5)
) + scale_x_continuous(limits=c(24.5,30.5),
                       breaks = c(25, 28, 30)
) + scale_y_continuous(limits = c(0, 0.75),
                       breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
) + labs(x="Mean Temperature [C]", y="% Eclosion"
) + theme(axis.line.x=element_line(colour = 'black', size = 1),
          axis.line.y=element_line(colour = 'black', size = 1),
          axis.ticks = element_line(colour = 'black', size = 1),
          axis.ticks.length = unit(2, "mm"),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          legend.background = element_rect(color="black",linetype="solid"),
          legend.position = c(.8, .85))

mn_psecl_plot3


#combine into one figure with cowplot
surv_fig_rwo <- plot_grid(mn_psem_plot3, mn_psecl_plot3, labels = c("A", "B", align="h"))
surv_fig_rwo




#--------------------


#look at effects of host growth between orig and repl experiment

#subset to only 30C mean temp treatments, removing +/-5
tvor_lng30 <- subset(tvor_lng, temp.avg==30 & temp.var!=5)


#find mean values of host growth and age between orig and repl expts

#log mass
mn_lmss_sum <- summarySE(tvor_lng30, measurevar = "log_mss",
                         groupvars = c("expt", "temp.var", "treatment", "instar"),
                         na.rm = TRUE)
mn_lmss_sum



#age
mn_age_sum <- summarySE(tvor_lng30, measurevar = "age",
                        groupvars = c("expt", "temp.var", "treatment", "instar"),
                        na.rm = TRUE)
mn_age_sum



#add age and age_se columns to log mass sum
mn_lmss_sum$age <- mn_age_sum[, 6]
mn_lmss_sum$age_se <- mn_age_sum[, 8]



#plot mean mass by mean age--color by expt, facet wrap by tempvar and treatment
mn_lma_plot <- ggplot(mn_lmss_sum, aes(x=age, y=log_mss, group=expt, color=expt))
mn_lma_plot + geom_point(aes(shape=expt), size=5
)+geom_line(size=1.2
)+geom_errorbar(aes(ymin=log_mss - se, ymax=log_mss + se),
                width=.5, size=1
)+geom_errorbarh(aes(xmin=age - age_se, xmax=age + age_se),
                 height=.5, size=1
)+facet_wrap(temp.var~treatment)



#plot mean mass by age--color by treatment, facet wrap by temp.var and expt
mn_lma_plot2 <- ggplot(mn_lmss_sum, aes(x=age, y=log_mss, group=treatment, color=treatment))
mn_lma_plot2 + geom_point(aes(shape=treatment), size=5
)+geom_line(size=1.2
)+geom_errorbar(aes(ymin=log_mss - se, ymax=log_mss + se),
                width=.5, size=1
)+geom_errorbarh(aes(xmin=age - age_se, xmax=age + age_se),
                 height=.5, size=1
)+facet_wrap(temp.var~expt)



#--------------------

#plot mn host mass and age orig data for 25 and 28 with repl data for 30


#subset out the orig 30C treatments--create keep sorting column
tvor_lng$keep <- ifelse(tvor_lng$temp.avg==30 & tvor_lng$expt=="orig", 0, 1)
tvor_lng_ro <- subset(tvor_lng, keep==1)


#find mean values of host growth and age 

#log mass
mn_lmssro_sum <- summarySE(tvor_lng_ro, measurevar = "log_mss",
                         groupvars = c("temp.avg", "temp.var", "treatment", "instar"),
                         na.rm = TRUE)
mn_lmssro_sum



#age
mn_agero_sum <- summarySE(tvor_lng_ro, measurevar = "age",
                        groupvars = c("temp.avg", "temp.var", "treatment", "instar"),
                        na.rm = TRUE)
mn_agero_sum



#add age and age_se columns to log mass sum
mn_lmssro_sum$age <- mn_agero_sum[, 6]
mn_lmssro_sum$age_se <- mn_agero_sum[, 8]



#plot mn log mass by mn age, color by temp.var, line type by treatment, facet wrap by temp.avg
mn_lmaro_plot <- ggplot(mn_lmssro_sum, aes(x=age, y=log_mss, group=interaction(temp.var, treatment),
                                         color=temp.var))
mn_lmaro_plot + geom_point(aes(shape=treatment),
                           size=5
) + geom_line(aes(linetype=treatment),
              size=1.7
) + geom_errorbar(aes(ymin=log_mss - se, ymax=log_mss + se),
                  width=1.7, size=1.2
) + geom_errorbarh(aes(xmin=age - age_se, xmax=age + age_se),
                   height=.3, size=1.2
) + scale_color_manual(values=c("#56B4E9","#D55E00"),name=c("Fluctuation [C]"),
                     breaks=c("0","10"),labels=c("0","10"),
                     guide=guide_legend(keywidth = 2.5)
) + scale_linetype_manual(values=c("solid","dashed"),name="Treatment",
                        breaks=c("control","para"),labels=c("Cotrol","Parasitized"),
                        guide=guide_legend(keywidth = 2.5)
) + scale_shape_manual(values = c(16,17),name="Treatment",
                     breaks=c("control","para"),labels=c("Cotrol","Parasitized"),
                     guide=guide_legend(keywidth = 2.5)
) + labs(x="Age [days]",y="Log(Mass) [mg]"
) + facet_wrap(~temp.avg
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
        legend.position = c(.9, .25))


#-------------------

#run a prelim gamm model of host mass and age by temp avg, temp var and treatment

#subset to only columns in model, remove rows with NAs (so that predicted and fitted values can be added
#to the dataframe easily)
tvor_mass <- select(tvor_lng_ro, bug.id, temp.avg, temp.var, treatment, log_mss, age)
tvor_mass <- na.omit(tvor_mass)

#make bug.id a factor so it will work as a random effect in the GAMM model
tvor_mass$bug.id <- as.factor(tvor_mass$bug.id)


#run a full GAMM model (where the smooth of age is also affected by the interaction of temp and treat)
gam_mass_mod<-gam(log_mss ~ s(age, by = interaction(treatment, temp.avg, temp.var, k=20, bs="ts")) 
                  + s(bug.id, bs ="re") + treatment * temp.avg * temp.var,
                  method="ML", data=tvor_mass, na.action = na.omit)

anova(gam_mass_mod)
summary(gam_mass_mod)

gam.check(gam_mass_mod, type = "deviance")




#Make a null model for model testing, where the smooth of age is not affected by temp.var, temp.avg or treat
gam_mnull_mod<-gam(log_mss ~ s(age, k=10, bs="ts") 
                   + s(bug.id, bs="re") + treatment * temp.avg * temp.var,
                   method="ML", data=tvor_mass, na.action = na.omit)

anova(gam_mnull_mod)


anova(gam_mnull_mod, gam_mass_mod)
AIC(gam_mnull_mod, gam_mass_mod)


#add predicted and residual values to model data set, plot results
tvor_mass$pred <- predict(gam_mass_mod, level=0)
tvor_mass$resid <- residuals(gam_mass_mod, level=0)



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


#using geom_smooth instead of geom_line to try and get an average line over temp.avg and not by individual?
#not sure if this is appropriate.
md_gammod_fit<-ggplot(tvor_mass, aes(x=age, y=log_mss, color=temp.avg))
md_gammod_fit+geom_point(size=3, shape=1
)+geom_smooth(aes(y=pred, group=temp.avg)
)+facet_wrap(treatment~temp.var)





#para treatment with load as a predictor

#subset tvor_lng_ro to only para
tvor_lng_rop <- subset(tvor_lng_ro, treatment=="para" & end.class!="wand")

#range of load
range(tvor_lng_rop$load)


#remove individuals with load over 300
tvor_lng_rop <- subset(tvor_lng_rop, load < 300)


#subset to only columns in model, and remove rows with NAs (so that predicted and fitted values can be
#added to the dataframe easily)
p_mass<-select(tvor_lng_rop, bug.id, temp.avg, temp.var, load, log_mss, age)
p_mass<-na.omit(p_mass)

#convert bug.id to factor so it functions properly as a random effect
p_mass$bug.id<-as.factor(p_mass$bug.id)


#run gam mode with age and load within a single smooth, with an interaction of temp.avg and temp.var
##check knots etc for fit
gam_pml_mod<-gam(log_mss ~ s(age, load, by=interaction(temp.avg, temp.var, bs="ts")) 
                 + s(bug.id, bs="re") + temp.avg * temp.var, method="ML", data=p_mass, na.action = na.omit)

anova(gam_pml_mod)
summary(gam_pml_mod)

gam.check(gam_pml_mod)

#run a GAMM with age and load as separate smooths--does not have an interaction with temp this way (2 2D surfaces, instead
## of one 3D surface)
gam_pml_nointmod<-gam(log_mss ~ s(age, by=interaction(temp.avg, temp.var, k=10, bs="ts")) 
                      + s(load, by=interaction(temp.avg, temp.var, k=10, bs="ts")) + s(bug.id, bs="re") 
                      + temp.avg * temp.var, method="ML", data=p_mass, na.action = na.omit)

anova(gam_pml_nointmod)
summary(gam_pml_nointmod)


#make a null model with an intercept (age and load in the same smooth)
gam_pml_null_mod<-gam(log_mss ~ s(age, load, bs="ts") 
                      + s(bug.id,bs="re") + temp.avg * temp.var, 
                      method="ML", data=p_mass, na.action = na.omit)


#Make a null model with no intercept (age and load in separate smooths)
gam_pml_null_nointmod<-gam(log_mss ~ s(age, bs="ts") 
                           + s(load, bs="ts") + s(bug.id,bs="re") 
                           + temp.avg * temp.var, method="ML", data=p_mass, na.action = na.omit)




#compare models with and without interaction of age and load to see which is better (and with null models)
#model without interaction seems best
anova(gam_pml_mod, gam_pml_nointmod, gam_pml_null_mod, gam_pml_null_nointmod)
AIC(gam_pml_mod, gam_pml_nointmod, gam_pml_null_mod, gam_pml_null_nointmod)



#Model with age and load in same smooth
#make columns with predicted and residual values for plotting
p_mass$pred<-predict(gam_pml_mod, level=0)
p_mass$resid<-residuals(gam_pml_mod, level=0)


#residuals against age
#residuals don't look great
pml_gammod_ra<-ggplot(p_mass, aes(x=age, y=resid, color=temp.avg))
pml_gammod_ra+geom_point(size=4, shape=1
)+geom_hline(aes(yintercept=0),
             color="black",
             size=1.5, linetype="dashed"
)+facet_wrap(~temp.var)


#residuals against load
pml_gammod_rl<-ggplot(p_mass, aes(x=load, y=resid, color=temp.avg))
pml_gammod_rl+geom_point(size=4, shape=1
)+geom_hline(aes(yintercept=0),
             color="black",
             size=1.5, linetype="dashed"
)+facet_wrap(~temp.var)


pml_gammod_fit<-ggplot(p_mass, aes(x=age, y=log_mss, color=load))
pml_gammod_fit+geom_point(size=3, shape=1
)+geom_line(aes(y=pred, group=interaction(temp.avg, bug.id))
)+scale_color_viridis(
)+facet_wrap(temp.avg~temp.var)




#Model with age and load in same smooth
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
)+facet_wrap(~temp.var)


#model fit, colored by load
pmlni_gammod_fit<-ggplot(p_mass, aes(x=age, y=log_mss, color=load))
pmlni_gammod_fit+geom_point(size=4, shape=1
)+geom_line(aes(y=pred_ni, group=interaction(temp.avg, bug.id)),
            size=1
)+scale_color_viridis(begin=.1, end=1,
                      option = "viridis"
)+facet_wrap(temp.avg~temp.var)




#--------------------------

#quick and dirty analysis of wasp survival, only looking at constant temps--to see if the repl data 
#influences effects of mean temp on survival at constant temps (so I don't have to worry about 
#WOWEs at 30.10)


#subset to only constant temps and parasitized treatment
tvor_con <- subset(tvor, temp.var==0 & treatment=="para" & end.class=="em")

#make a "keep" column to remove individuals from the orig 30C treatment
tvor_con$keep <- ifelse(tvor_con$temp.avg==30 & tvor_con$expt=="orig", 0, 1)

#subset out the orig 30C data
tvor_con <- subset(tvor_con, keep==1)


#glm model with quasibinomial distribution to look at overdispersion
wtotsurv_mod1<-glm(cbind(num.ecl,tot.died) ~ temp.avg * load,
                        family=quasibinomial,
                        data=tvor_con,
                        na.action = na.omit)

anova(wtotsurv_mod1)
summary(wtotsurv_mod1)


#scale load
tvor_con$resc.ld <- rescale(tvor_con$load, to = c(0,1))

#quite over dispersed, so adding a random effect of individual
wtots_re_mod1 <- glmer(cbind(num.ecl, tot.died) ~ temp.avg * resc.ld +(1|bug.id),
                       family = binomial,
                       data=tvor_con,
                       na.action = na.omit)

anova(wtots_re_mod1)
summary(wtots_re_mod1)



#finding p values of fixed effects

wtots_re_mod_ta <- glmer(cbind(num.ecl, tot.died) ~ temp.avg +(1|bug.id),
                         family = binomial,
                         data=tvor_con,
                         na.action = na.omit)


wtots_re_mod_ld <- glmer(cbind(num.ecl, tot.died) ~ resc.ld +(1|bug.id),
                         family = binomial,
                         data=tvor_con,
                         na.action = na.omit)



wtots_re_mod_int <- glmer(cbind(num.ecl, tot.died) ~ temp.avg:resc.ld +(1|bug.id),
                          family = binomial,
                          data=tvor_con,
                          na.action = na.omit)


#null model
wtots_re_mod_null <- glmer(cbind(num.ecl, tot.died) ~ 1 +(1|bug.id),
                           family = binomial,
                           data=tvor_con,
                           na.action = na.omit)


anova(wtots_re_mod_null, wtots_re_mod1, wtots_re_mod_ta, wtots_re_mod_ld, wtots_re_mod_int, test="Chisq")


#-------------

#plot num.ecl by load for orig 25, 28 data and 30 repl data

#create keep column to sort to only the orig 25, 28 and the repl 30 data
tvor$expt_keep <- ifelse(tvor$temp.avg==30 & tvor$expt=="orig", 0, 1)

#subset to only orig 25, 28 and repl 30 data
tvor_ro <- subset(tvor, expt_keep==1)

#subset to only para treatment that had emergence
tvor_rop <- subset(tvor_ro, treatment=="para" & end.class=="em")


#plot num.ecl by load, color by temp.avg, facet_wrap by temp.var
ecl_load_plot <- ggplot(tvor_rop, aes(x=load, y=num.ecl, group=temp.avg, color=temp.avg))
ecl_load_plot + geom_point(size=4, alpha=.7
) + geom_smooth(method="lm", se=FALSE, size=1.5
) + scale_color_manual(values=c("#009E73","#E69F00","#000000"),name=c("Avg. Temp. [C]"),
                     breaks=c("25","28","30"),labels=c("25","28","30"),
                     guide=guide_legend(keywidth=3)
) + facet_wrap(~temp.var)




#-----------------------------


#analysis of wasp survival for only 25 and 28 mean temp treatments--removing 30 so the 30.10 complete separation
#doesn't cause the model to not work properly

tvor_ncs <- subset(tvor, temp.avg!=30)

tvor_ncs$resc.ld <- rescale(tvor_ncs$load, to = c(0,1))

#set tvor_ncs so that 30 doesn't show up as a factor level, same with temp var 
tvor_ncs$temp.avg <- factor(tvor_ncs$temp.avg, levels = c(25, 28))
tvor_ncs$temp.var <- factor(tvor_ncs$temp.var, levels=c(0,10))



wtots_ncs_mod1 <- glmer(cbind(num.ecl, tot.died) ~ temp.avg * temp.var * resc.ld + (1|bug.id),
                        data = tvor_ncs,
                        family = binomial,
                        na.action = na.omit)

anova(wtots_ncs_mod1, test="Chisq")
summary(wtots_ncs_mod1)



#testing fixed effects

wtots_ncs_mod_ta <- glmer(cbind(num.ecl, tot.died) ~ temp.avg + (1|bug.id),
                          data = tvor_ncs,
                          family = binomial,
                          na.action = na.omit)


wtots_ncs_mod_tv <- glmer(cbind(num.ecl, tot.died) ~ temp.var + (1|bug.id),
                          data = tvor_ncs,
                          family = binomial,
                          na.action = na.omit)


wtots_ncs_mod_ld <- glmer(cbind(num.ecl, tot.died) ~ resc.ld + (1|bug.id),
                          data = tvor_ncs,
                          family = binomial,
                          na.action = na.omit)


wtots_ncs_mod_tatv <- glmer(cbind(num.ecl, tot.died) ~ temp.avg:temp.var + (1|bug.id),
                            data = tvor_ncs,
                            family = binomial,
                            na.action = na.omit)


wtots_ncs_mod_tald <- glmer(cbind(num.ecl, tot.died) ~ temp.avg:resc.ld + (1|bug.id),
                            data = tvor_ncs,
                            family = binomial,
                            na.action = na.omit)


wtots_ncs_mod_tvld <- glmer(cbind(num.ecl, tot.died) ~ temp.var:resc.ld + (1|bug.id),
                            data = tvor_ncs,
                            family = binomial,
                            na.action = na.omit)


wtots_ncs_mod_tatvld <- glmer(cbind(num.ecl, tot.died) ~ temp.avg:temp.var:resc.ld + (1|bug.id),
                              data = tvor_ncs,
                              family = binomial,
                              na.action = na.omit)


wtots_ncs_mod_null <- glmer(cbind(num.ecl, tot.died) ~ 1 + (1|bug.id),
                            data = tvor_ncs,
                            family = binomial,
                            na.action = na.omit)


anova(wtots_ncs_mod_null, wtots_ncs_mod1, wtots_ncs_mod_ta, wtots_ncs_mod_tv, wtots_ncs_mod_ld, wtots_ncs_mod_tatv,
      wtots_ncs_mod_tald, wtots_ncs_mod_tvld, wtots_ncs_mod_tatvld, test="Chisq")


