#comparing 25 and 30 constant treatments between CxPxT and TV experiments


#load libraries
library(scales)
library(readr)
library(plyr)
library(ggplot2)
library(Rmisc)
library(dplyr)
library(tidyr)
library(reshape2)
library(cowplot)
library(viridis)
library(extrafont)
library(lme4)
library(MuMIn)
library(mgcv)


#load data

#CxPxT
cpt <- read_csv("data files/cpt gr wide.csv")


#TV
tvor <- read_csv("data files/Ms-Cc_tv-orig-rep_comb_cl.csv")


#-----------------------------

#CxPxT cleaning

#remove columns not needed for this comparison (all diet and consumption measurements)
cpt <- cpt %>% select(bug.id, treatment, temp, date.hatch, date.died, died, date.ovp, suc.ovp, date.3, date.4,
                      date.5, date.wander, date.em, instar.em, mass.T0, mass.4, mass.5, mass.wander, mass.befem,  
                      num.em, num.coc, num.fail.spin, num.unem, load)

#rename mass.T0 to be mass.3
cpt <- rename(cpt, mass.3 = mass.T0)


#calculate dev times
cpt$tt3 <- cpt$date.3 - cpt$date.hatch
cpt$tt4 <- cpt$date.4 - cpt$date.hatch
cpt$tt5 <- cpt$date.5 - cpt$date.hatch
cpt$ttwand <- cpt$date.wander - cpt$date.hatch
cpt$ttem.h <- cpt$date.em - cpt$date.hatch
cpt$ttem.w <- cpt$date.em - cpt$date.ovp


#calculate wasp survival
cpt$ps.em <- cpt$num.em/cpt$load


#remove dead individuals
cpt <- subset(cpt, died==0)


#remove individuals with loads greater than 300
cpt$keep_ld <- ifelse(cpt$treatment=="para" & cpt$load > 300, 0, 1)
cpt <- subset(cpt, keep_ld==1)


#create subset of only parasitized individuals
cpt_p <- subset(cpt, treatment=="para")


#----------------------

#TV data cleaning

#remove parasitized wanderers and WOWEs that wandered
tvor$keep_p <- ifelse(tvor$treatment=="para" & tvor$end.class=="wand", 0, 1)

tvor <- subset(tvor, keep_p==1)


#remove parasitized individuals with too large loads (>300)
tvor$keep_ld <- ifelse(tvor$treatment=="para" & tvor$load > 300, 0, 1)
tvor <- subset(tvor, keep_ld==1)


#remove WOWEs in 30C treatment
tvor$keep_p2 <- ifelse(tvor$temp.avg==30 & tvor$temp.var==0 & tvor$treatment=="para" & tvor$end.class!="em", 0, 1)
tvor <- subset(tvor, keep_p2==1)


#remove +/- 5 temp fluc treatment
tvor <- subset(tvor, temp.var!=5)


#create subset of parasitized individuals
tvor_p <- subset(tvor, treatment=="para")


#------------------------

#COMBINE DATA SETS

colnames(cpt)
colnames(tvor)



#CPT changes

#make experiment column for cpt to differentiate those rows from tvor, make temp.var column (0),
#make other columns with NAs to match tvor
cpt$expt <- "cpt"
cpt$temp.var <- 0
cpt$date.cull <- NA
cpt$mass.cull <- NA
cpt$ttcul <- NA

#create ttend and mass end columns for cpt
cpt$ttend <- coalesce(cpt$ttwand, cpt$ttem.h)
cpt$mass.end <- coalesce(cpt$mass.wander, cpt$mass.befem)

#making an end.class column for cpt
cpt$date.wander[is.na(cpt$date.wander)]<-0
cpt$date.em[is.na(cpt$date.em)] <- 0

cpt$end.class <- ifelse(cpt$date.wander > 0, "wand",
                        ifelse(cpt$date.em > 0, "em", "unk"))


#remove individuals with an end class of "unknown"
cpt <- subset(cpt, end.class!="unk")


#rename temp, wander and emergence columns
cpt <-  cpt %>% rename(date.wand = date.wander, mass.wand = mass.wander, mass.48em = mass.befem, temp.avg = temp)



#TVOR changes

#remove columns from tvor that are unneeded
tvor <- tvor %>% select(-c(ps.ecl, diet, died.bf5, left_out, mongo, keep, keep2, keep_p, keep_ld, keep_p2, num.ecl, ttecl, 
                           date.ecl.j, bled.em))


#rename date columns to match cpt
tvor <- tvor %>% rename(date.hatch = date.hatch.j, date.3 = date.3.j, date.4 = date.4.j, date.5 = date.5.j,
                        date.wand = date.wand.j, date.em = date.em.j, date.cull = date.cull.j, 
                        date.died = date.died.j, date.ovp = date.ovp.j)



#COMBINE DATA SETS
all_dat <- bind_rows(cpt, tvor)



#-----------------------

#find mean wasp survival for each experiment and temperature

#subset combined dataset to only be parasitized treatments
ad_p <- subset(all_dat, treatment=="para")

#NUMBER EMERGED

ws_num_sum <- summarySE(ad_p, measurevar = "num.em",
                        groupvars = c("temp.avg", "expt"),
                        na.rm = TRUE)
ws_num_sum



#PROP EMERGED

ws_psem_sum <- summarySE(ad_p, measurevar = "ps.em",
                         groupvars = c("temp.avg", "expt"),
                         na.rm = TRUE)

ws_psem_sum



#---------------------

#plot load effects on number emerged, compared by experiment

#remove all temp treatments except 25 and 30 (repl)
ad_psub <- subset(ad_p, temp.avg!=20 | temp.avg!=28)
ad_psub <- subset(ad_psub, temp.var==0)

ad_psub$keep <- ifelse(ad_psub$temp.avg==30 & ad_psub$expt=="orig", 0, 1)
ad_psub <- subset(ad_psub, keep==1)


#rename "repl" and "orig" to tvor, for easier comparison in plots and analyses
ad_psub$expt <- ifelse(ad_psub$expt=="orig" | ad_psub$expt=="repl", "tvor", ad_psub$expt)


#plot of num em by load, colored by expt, faceted by temp
com_wsnmem_plot <- ggplot(ad_psub, aes(x=load, y=num.em, group=expt, color=expt))
com_wsnmem_plot + geom_point(size=6
) + geom_smooth(method = "lm"
) + scale_color_manual(values=c("black", "#E69F00")
) + facet_wrap(~temp.avg)


#----------------

#quick and dirty analysis of wasp survival to emergence between the two experiments

#make temp.avg factor
ad_psub$temp.avg <- factor(ad_psub$temp.avg)


wsem_mod1 <- glm(cbind(num.em, num.unem) ~ temp.avg * load * expt,
                 family = binomial,
                 data = ad_psub,
                 na.action=na.omit)
anova(wsem_mod1)
summary(wsem_mod1)



#model selection using dredge
ad_wsem <- select(ad_psub, bug.id, temp.avg, load, expt, num.em, num.unem)
ad_wsem <- drop_na(ad_wsem)


wsem_mod_fd <- glm(cbind(num.em, num.unem) ~ temp.avg * load * expt,
                   family = binomial,
                   data = ad_wsem,
                   na.action=na.fail)

wsem_dredge <- dredge(wsem_mod_fd)
wsem_dredge #best model is full model


#----------------------

#analysis of wasp survival to emergence, with random effect of individual

#rescale load
ad_psub$resc_ld <- rescale(ad_psub$load, to=c(0,1))

wsem_re_mod <- glmer(cbind(num.em, num.unem) ~ temp.avg * resc_ld * expt + (1|bug.id),
                   family = binomial,
                   data = ad_psub,
                   na.action=na.omit,
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

anova(wsem_re_mod)
summary(wsem_re_mod)


#model selection using dredge
ad_wre <- select(ad_psub, bug.id, temp.avg, resc_ld, expt, num.em, num.unem)
ad_wre <- drop_na(ad_wre)


wsem_re_mod_fd <- glmer(cbind(num.em, num.unem) ~ temp.avg * resc_ld * expt + (1|bug.id),
                     family = binomial,
                     data = ad_wre,
                     na.action=na.fail,
                     control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

wsem_re_dredge <- dredge(wsem_re_mod_fd)
wsem_re_dredge  #best model excludes expt:load and expt:load:temp.avg interactions


#best model
wsem_best <- glmer(cbind(num.em, num.unem) ~ temp.avg + resc_ld + expt 
                   + temp.avg:resc_ld + temp.avg:expt + (1|bug.id),
                   family = binomial,
                   data = ad_psub,
                   na.action=na.omit,
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))



#testing reduced models (without term of interest [and associated interaction terms, if load or expt]) to best fit model

wsem_re_mod_ta <- glmer(cbind(num.em, num.unem) ~ resc_ld + expt + (1|bug.id),
                        family = binomial,
                        data = ad_psub,
                        na.action=na.omit,
                        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))



wsem_re_mod_ld <- glmer(cbind(num.em, num.unem) ~ temp.avg + expt + temp.avg:expt + (1|bug.id),
                        family = binomial,
                        data = ad_psub,
                        na.action=na.omit,
                        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))



wsem_re_mod_expt <- glmer(cbind(num.em, num.unem) ~ temp.avg + resc_ld 
                          + temp.avg:resc_ld + (1|bug.id),
                          family = binomial,
                          data = ad_psub,
                          na.action=na.omit,
                          control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))



wsem_re_mod_exptta <- glmer(cbind(num.em, num.unem) ~ temp.avg + resc_ld + expt 
                            + temp.avg:resc_ld + (1|bug.id),
                            family = binomial,
                            data = ad_psub,
                            na.action=na.omit,
                            control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))


wsem_re_mod_tald <- glmer(cbind(num.em, num.unem) ~ temp.avg + resc_ld + expt 
                          + temp.avg:expt + (1|bug.id),
                          family = binomial,
                          data = ad_psub,
                          na.action=na.omit,
                          control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))



anova(wsem_best, wsem_re_mod_ta, test="Chisq")
anova(wsem_best, wsem_re_mod_ld, test="Chisq")
anova(wsem_best, wsem_re_mod_expt, test="Chisq")
anova(wsem_best, wsem_re_mod_exptta, test="Chisq")
anova(wsem_best, wsem_re_mod_tald, test="Chisq")




#--------------------

#Joel noted two points in cpt data that have much lower num.em at 30 than the rest of the data, suggested 
#removing and seeing how that affected slope/difference between tv data

which(ad_psub$load > 100 & ad_psub$num.em < 25)

ad_psub2 <- ad_psub[-c(10, 43), ]


#plot of num em by load, colored by expt, faceted by temp--two outliers removed
com_wsnmem_plot2 <- ggplot(ad_psub2, aes(x=load, y=num.em, group=expt, color=expt))
num_em <- com_wsnmem_plot2 + geom_point(size=6
) + geom_smooth(method = "lm"
) + scale_color_manual(values=c("black", "#E69F00"),
                       breaks = c("cpt", "tvor"),
                       labels = c("Moore 2019", "Current Study"),
                       name = "Data"
) + labs(x="Total Load", y="Number Emerged"
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
          legend.position = "none")




#plot of time to emergence by load, colored by expt, faceted by temp--two outliers removed
com_wsttem_plot2 <- ggplot(ad_psub2, aes(x=load, y=ttem.w, group=expt, color=expt))
dev_time <- com_wsttem_plot2 + geom_point(size=6
) + geom_smooth(method = "lm"
) + scale_color_manual(values=c("black", "#E69F00"),
                       breaks = c("cpt", "tvor"),
                       labels = c("Moore 2019", "Current Study"),
                       name = "Data"
) + labs(x="Total Load", y="Dev. Time to Emergence"
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
          legend.position = c(.875, .8))



surv_dev <- plot_grid(num_em, dev_time, labels=c("A", "B"), align = "v", ncol=1)
surv_dev




#NUMBER EMERGED

ws_num_sum <- summarySE(ad_psub2, measurevar = "num.em",
                        groupvars = c("temp.avg", "expt"),
                        na.rm = TRUE)
ws_num_sum



#PROP EMERGED

ws_psem_sum <- summarySE(ad_psub2, measurevar = "ps.em",
                         groupvars = c("temp.avg", "expt"),
                         na.rm = TRUE)

ws_psem_sum




#quick and dirty analysis of wasp survival to emergence between the two experiments--two outliers removed

wsem_mod2 <- glm(cbind(num.em, num.unem) ~ temp.avg * load * expt,
                 family = quasibinomial,
                 data = ad_psub2,
                 na.action=na.omit)
anova(wsem_mod2)
summary(wsem_mod2)


#overdispersed, so adding random effect

wsem_re_mod2 <- glmer(cbind(num.em, num.unem) ~ temp.avg * resc_ld * expt + (1|bug.id),
                      family = binomial,
                      data = ad_psub2,
                      na.action=na.omit,
                      control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

anova(wsem_re_mod2)
summary(wsem_re_mod2)



#model selection using dredge
ad_wre2 <- select(ad_psub2, bug.id, temp.avg, resc_ld, expt, num.em, num.unem)
ad_wre2 <- drop_na(ad_wre2)


wsem_re_mod_fd2 <- glmer(cbind(num.em, num.unem) ~ temp.avg * resc_ld * expt + (1|bug.id),
                        family = binomial,
                        data = ad_wre2,
                        na.action=na.fail,
                        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

wsem_re_dredge2 <- dredge(wsem_re_mod_fd2)
wsem_re_dredge2  #best model excludes expt:load and expt:load:temp.avg interactions


#best model
wsem_best2 <- glmer(cbind(num.em, num.unem) ~ temp.avg + resc_ld + expt 
                   + temp.avg:resc_ld + temp.avg:expt + (1|bug.id),
                   family = binomial,
                   data = ad_psub2,
                   na.action=na.omit,
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))



#testing reduced models (without term of interest [and associated interaction terms, if load or expt]) to best fit model

wsem_re_mod_ta2 <- glmer(cbind(num.em, num.unem) ~ resc_ld + expt + (1|bug.id),
                        family = binomial,
                        data = ad_psub2,
                        na.action=na.omit,
                        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))



wsem_re_mod_ld2 <- glmer(cbind(num.em, num.unem) ~ temp.avg + expt + temp.avg:expt + (1|bug.id),
                        family = binomial,
                        data = ad_psub2,
                        na.action=na.omit,
                        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))



wsem_re_mod_expt2 <- glmer(cbind(num.em, num.unem) ~ temp.avg + resc_ld 
                          + temp.avg:resc_ld + (1|bug.id),
                          family = binomial,
                          data = ad_psub2,
                          na.action=na.omit,
                          control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))



wsem_re_mod_exptta2 <- glmer(cbind(num.em, num.unem) ~ temp.avg + resc_ld + expt 
                            + temp.avg:resc_ld + (1|bug.id),
                            family = binomial,
                            data = ad_psub2,
                            na.action=na.omit,
                            control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))


wsem_re_mod_tald2 <- glmer(cbind(num.em, num.unem) ~ temp.avg + resc_ld + expt 
                          + temp.avg:expt + (1|bug.id),
                          family = binomial,
                          data = ad_psub2,
                          na.action=na.omit,
                          control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))



anova(wsem_best2, wsem_re_mod_ta2, test="Chisq")
anova(wsem_best2, wsem_re_mod_ld2, test="Chisq")
anova(wsem_best2, wsem_re_mod_expt2, test="Chisq")
anova(wsem_best2, wsem_re_mod_exptta2, test="Chisq")
anova(wsem_best2, wsem_re_mod_tald2, test="Chisq")





#----------------------

#effects on wasp development time

#plot of time to emergence by load, colored by expt, faceted by temp--two outliers removed
com_wsttem_plot2 <- ggplot(ad_psub2, aes(x=load, y=ttem.w, group=expt, color=expt))
com_wsttem_plot2 + geom_point(size=6
) + geom_smooth(method = "lm"
) + scale_color_manual(values=c("black", "#E69F00"),
                       breaks = c("cpt", "tvor"),
                       labels = c("Moore 2019", "Current Study"),
                       name = "Data"
) + labs(x="Total Load", y="Number Emerged"
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
          legend.position = c(.9, .8))



#box plot of time to emergence, colored by expt, faceted by temp--two outliers removed
com_wsttem_boxplot <- ggplot(ad_psub2, aes(x=expt, y=ttem.w, group=expt, fill=expt))
com_wsttem_boxplot + geom_boxplot(
) + scale_color_manual(values=c("black", "#E69F00")
) + facet_wrap(~temp.avg)


#mean time to emergence
ws_ttem_sum <- summarySE(ad_psub2, measurevar = "ttem.w",
                        groupvars = c("temp.avg", "expt"),
                        na.rm = TRUE)
ws_ttem_sum


#---------------------

#analysis of wasp development time

wdev_mod <- lm(ttem.w ~ temp.avg * load *expt,
               data = ad_psub2,
               method = "ML",
               na.action=na.omit)

anova(wdev_mod)
summary(wdev_mod)


#model selection with dredge:
ad_dev <- select(ad_psub2, bug.id, temp.avg, load, expt, ttem.w)
ad_dev <- drop_na(ad_dev)


wdev_mod_fd <- lm(ttem.w ~ temp.avg * load *expt,
               data = ad_dev,
               method = "ML",
               na.action=na.fail)

wdev_dredge <- dredge(wdev_mod_fd)
wdev_dredge #best model removes all interactions


wdev_best <- lm(ttem.w ~ temp.avg + load + expt,
                data=ad_psub2,
                method = "ML",
                na.action = na.omit)

anova(wdev_best)



#--------------------------------

#combine wandering and emergence into end columm for preparation to convert data frame to long format

#remove para treatment that wandered
all_dat$date.wand[is.na(all_dat$date.wand)] <- 0
all_dat$keep <- ifelse(all_dat$treatment=="para" & all_dat$date.wand > 0, 0, 1)

all_dat <- subset(all_dat, keep==1)

all_dat$date.wand[all_dat$date.wand==0] <- NA
all_dat$keep <- NULL

#make combined data set into long format to compare host growth and development time

all_dat_sub <- subset(all_dat, temp.avg!=20)

ad_lng <- gather(all_dat_sub, instar, mass, mass.3, mass.4, mass.5, mass.end)
ad_lng$instar <- gsub("mass.", "", ad_lng$instar)
ad_lng$instar <- gsub("48", "", ad_lng$instar)


ad_age_lng <- gather(all_dat_sub, instar, age, tt3, tt4, tt5, ttend)
ad_age_lng$instar <- gsub("tt", "", ad_age_lng$instar)
ad_age_lng$instar <- gsub(".h", "", ad_age_lng$instar)


ad_age_lng <- select(ad_age_lng, bug.id, temp.avg, temp.var, instar, age)

ad_lng <- merge(ad_lng, ad_age_lng, by=c("bug.id", "temp.avg", "temp.var", "instar"))


#remove the orig 30 data
ad_lng$keep <- ifelse(ad_lng$temp.avg==30 & ad_lng$expt=="orig", 0, 1)

ad_lng <- subset(ad_lng, keep==1)


#rename orig and repl experiments to just be tvor
ad_lng$expt <- ifelse(ad_lng$expt=="orig" | ad_lng$expt=="repl", "tvor", ad_lng$expt)


#make temp.avg a factor
ad_lng$temp.avg <- as.factor(ad_lng$temp.avg)

#---------------------

#plot of growth and development for hosts between temp avgs and expts

#take log of mass
ad_lng$log_mss <- log(ad_lng$mass)

#means of mass and age
lmss_sum <- summarySE(ad_lng, measurevar = "log_mss",
                      groupvars = c("temp.avg", "treatment", "expt", "instar"),
                      na.rm = TRUE)
lmss_sum


age_sum <- summarySE(ad_lng, measurevar = "age",
                     groupvars = c("temp.avg", "treatment", "expt", "instar"),
                     na.rm = TRUE)
age_sum


lmss_sum$age <- age_sum[ , 6]
lmss_sum$age_se <- age_sum[ , 8]


facet_labels <- c("control" = "Unparasitized", "para" = "Parasitized", "25" = "25", "30" = "30")


#plot mean mass and mean age for each temp.avg, treatment and expt
mn_ma_plot <- ggplot(lmss_sum, aes(x=age, y=log_mss, group=expt, color=expt))
mn_ma_plot + geom_point(size=6
) + geom_line(size=2
) + geom_errorbar(aes(ymin=log_mss - se, ymax=log_mss + se),
                  width=1, size=1.2
) + geom_errorbarh(aes(xmin=age - age_se, xmax=age + age_se),
                   height=1, size=1.2
) + scale_color_manual(values=c("black", "#E69F00"),
                       breaks = c("cpt", "tvor"),
                       labels = c("Moore 2019", "Current Study"),
                       name = "Data"
) + labs(x="Age", y="Log(Mass)"
) + facet_wrap(treatment~temp.avg, labeller = as_labeller(facet_labels)
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
          legend.position = c(.9, .1))




#---------------------

#GAMM analysis of host growth and mass

#subset to only columns in model, remove rows with NAs
ad_mass <- select(ad_lng, bug.id, temp.avg, expt, treatment, log_mss, age)
ad_mass <- drop_na(ad_mass)

#make bug.id a factor so it will work as a random effect in the GAMM model
ad_mass$bug.id <- as.factor(ad_mass$bug.id)


#Full GAMM model 
#Response = log of mass. Smooth by age, with an interaction of parasitization treatment, mean temperature and
#temperature fluctuation. Smooth of individual added to act as a random effect. Parasitization treatment, mean
#temperature and temperature fluctuation used as fixed effects. 
gam_mass_mod<-gam(log_mss ~ s(age, by = interaction(treatment, temp.avg, expt, k=400, bs="ts")) 
                  + s(bug.id, bs ="re") + treatment * temp.avg * expt,
                  method="ML", data=ad_mass, na.action = na.omit)


anova(gam_mass_mod)



gam.check(gam_mass_mod)


#-------------------------

#calculate means and range of sample sizes for each treatment

temp.avg <- c(25, 25, 25, 25, 
              28, 28, 28, 28, 
              30, 30, 30, 30)

temp.var <- c(0, 10, 0, 10, 
              0, 10, 0, 10,  
              0, 10, 0, 10)

treatment <- c("control", "control", "para", "para",
               "control", "control", "para", "para",
               "control", "control", "para", "para")

n <- c(35, 32, 37, 38, 31, 26, 30, 24, 43, 42, 47, 31)


ncalc <- data.frame(temp.avg, temp.var, treatment, n)

View(ncalc)


nsum_ta <- summarySE(ncalc, measurevar = "n",
                  groupvars = "temp.avg")
nsum_ta


nsum_tv <- summarySE(ncalc, measurevar = "n",
                     groupvars = "temp.var")
nsum_tv


nsum_trt <- summarySE(ncalc, measurevar = "n",
                      groupvars = "treatment")
nsum_trt


range(ncalc$n)
mean(ncalc$n)


#------------------------

#analyze 25 and 30 const data from const mss expt along with other tv mss expt data, see if results differ in a way that
#would affect interpretation

#create column that combines mean temp, temp fluc and expt to make subsetting easier
all_dat <- unite(all_dat, temp_trt_expt, temp.avg, temp.var, expt, sep = "_", remove=FALSE)


all_dat$keep <- ifelse(all_dat$temp_trt_expt=="30_0_orig" | all_dat$temp_trt_expt=="30_0_repl" 
                       | all_dat$temp_trt_expt=="30_10_orig" | all_dat$temp_trt_expt=="25_0_orig" 
                       | all_dat$temp_trt_expt=="20_0_cpt", 0, 1)

ad_comp <- subset(all_dat, keep==1)


#make long data frame
adc_lng <- gather(ad_comp, instar, mass, mass.3, mass.4, mass.5, mass.end)
adc_lng$instar <- gsub("mass.", "", adc_lng$instar)

adc_age <- gather(ad_comp, instar, age, tt3, tt4, tt5, ttend)
adc_age$instar <- gsub("tt", "", adc_age$instar)

adc_age <- select(adc_age, bug.id, treatment, temp.avg, temp.var, instar, age)

adc_lng <- merge(adc_lng, adc_age, by=c("bug.id", "treatment", "temp.avg", "temp.var", "instar"))


#make sure columns are correct class
adc_lng$temp.avg <- as.factor(adc_lng$temp.avg)
adc_lng$temp.var <- as.factor(adc_lng$temp.var)



#----------------------------

#ANALYSIS OF HOST MASS USING GAMM MODELS--comparing the results using 25 const and 30 const from cpt data to full tv analysis

#calculate log mass
adc_lng$log_mss <- log(adc_lng$mass)

#subset to only columns in model, remove rows with NAs (so that predicted and fitted values can be added
#to the dataframe easily)
comp_mass <- select(adc_lng, bug.id, temp.avg, temp.var, treatment, log_mss, age)
comp_mass <- na.omit(comp_mass)

#make bug.id a factor so it will work as a random effect in the GAMM model
comp_mass$bug.id <- as.factor(comp_mass$bug.id)


#Full GAMM model 
#Response = log of mass. Smooth by age, with an interaction of parasitization treatment, mean temperature and
#temperature fluctuation. Smooth of individual added to act as a random effect. Parasitization treatment, mean
#temperature and temperature fluctuation used as fixed effects. 
gam_mass_mod_cpt<-gam(log_mss ~ s(age, by = interaction(treatment, temp.avg, temp.var, k=1000, bs="ts")) 
                  + s(bug.id, bs ="re") + treatment * temp.avg * temp.var,
                  method="ML", data=comp_mass, na.action = na.omit)


anova(gam_mass_mod_cpt)
summary(gam_mass_mod_cpt)

gam.check(gam_mass_mod, type = "deviance")


#See if tv gam model can have the fixed effect of fluc temp removed, and if that significantly worsens the model
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



#try removing the fluc temp fixed effect
gam_mass_mod_test<-gam(log_mss ~ s(age, by = interaction(treatment, temp.avg, temp.var, k=20, bs="ts")) 
                  + s(bug.id, bs ="re") + treatment * temp.avg + temp.avg:temp.var + temp.var:treatment 
                  + temp.avg:temp.var:treatment,
                  method="ML", data=tvor_mass, na.action = na.omit)


anova(gam_mass_mod_test)

anova(gam_mass_mod, gam_mass_mod_test, test="Chisq")
AIC(gam_mass_mod, gam_mass_mod_test)


#plot gamm model fit and residuals
#add predicted and residual values to model data set, plot results
comp_mass$pred <- predict(gam_mass_mod, level=0)
comp_mass$resid <- residuals(gam_mass_mod, level=0)



#make a label for the facet wrap panels
fnames <- c("control" = "NP", "para" = "P", "0" = "0", "10" = "10")



#plot model residuals against age, color by mean temperature, facet by parasitization treatment and fluctuation
md_gammod_ra<-ggplot(comp_mass, aes(x=age, y=resid, color=temp.avg))
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




md_gammod_fit<-ggplot(comp_mass, aes(x=age, y=log_mss, group=interaction(bug.id, temp.avg), color=temp.avg))
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





#MEAN MASS BY MEAN AGE FOR MANDUCA SEXTA

#Calcualte means and variance of mass and age data
lmass_sum <- summarySE(adc_lng, measurevar = "log_mss",
                       groupvars = c("temp.avg", "temp.var", "treatment", "instar"),
                       na.rm = TRUE)

lmass_sum


age_sum <- summarySE(adc_lng, measurevar = "age",
                     groupvars = c("temp.avg", "temp.var", "treatment", "instar"),
                     na.rm = TRUE)
age_sum


#add age and age standard error values to mass summary object
lmass_sum$age <- age_sum[, 6]
lmass_sum$age_se <- age_sum[, 8]



#Plot mean mass by mean age, facetted by mean reaing temperature, grouped by rearing temperature variation
##and parasitization treatment. Color by temp var, line type and point shape by parasitization treatment
##error bars are standard error for mass (y axis) and age (x axis)

mn_lma_plot <- ggplot(lmass_sum, aes(x=age, y=log_mss, group=interaction(temp.var, treatment), color=temp.var))
mn_lma_plot + geom_point(aes(shape=treatment),
                         size=6, stroke=2
) + geom_line(aes(linetype=treatment),
              size=2
) + geom_errorbar(aes(ymin=log_mss - se, ymax=log_mss + se),
                  width=2, size=1.2
) + geom_errorbarh(aes(xmin=age - age_se, xmax=age + age_se),
                   height=.4, size=1.2
) + scale_color_manual(values=c("#56B4E9","#D55E00"),name=c("Fluctuation [C]"),
                       breaks=c("0","10"),labels=c("0","10"),
                       guide=guide_legend(keywidth = 6, keyheight = 1.5)
) + scale_linetype_manual(values=c("solid","dashed"),name="Treatment",
                          breaks=c("control","para"),labels=c("NP","P"),
                          guide=guide_legend(keywidth = 6, keyheight = 1.5)
) + scale_shape_manual(values = c(16,2),name="Treatment",
                       breaks=c("control","para"),labels=c("NP","P"),
                       guide=guide_legend(keywidth = 6, keyheight = 1.5)
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


#---------------------

#gamm model with 25 const and 30 const from cpt experiment 
#Looking at parasitized treatment with effects of load

adcp_lng <- subset(adc_lng, treatment=="para")


#subset to only columns in model, and remove rows with NAs (so that predicted and fitted values can be
#added to the dataframe easily)
p_mass<-select(adcp_lng, bug.id, temp.avg, temp.var, load, log_mss, age)
p_mass<-na.omit(p_mass)

#convert bug.id to factor so it functions properly as a random effect
p_mass$bug.id<-as.factor(p_mass$bug.id)


#run a GAMM with age and load as separate smooths--does not have an interaction with temp this way (2 2D surfaces, instead
## of one 3D surface)
gam_pml_nointmod<-gam(log_mss ~ s(age, by=interaction(temp.avg, temp.var, bs="ts")) 
                      + s(load, by=interaction(temp.avg, temp.var, k=100, bs="ts")) + s(bug.id, bs="re") 
                      + temp.avg * temp.var, method="ML", data=p_mass, na.action = na.omit)

anova(gam_pml_nointmod)
summary(gam_pml_nointmod)

gam.check(gam_pml_nointmod)



#------------------------

#ANALYSES OF WASP SURVIVAL TO EMERGENCE AND ECLOSION

adcp <- subset(ad_comp, treatment=="para")

#make temp.avg a factor
adcp$temp.avg <- factor(adcp$temp.avg)

#Remove the 30+/-10 treatment, as no wasps survived, the mechanism is not the same as the other treatments

#create combo column with temp.avg and temp.var
adcp <- unite(adcp, "tatv", temp.avg, temp.var, remove = FALSE)

adcp_nw <- subset(adcp, tatv!="30_10")

#make temp var a character instead of factor
adcp_nw$temp.var <- as.character(adcp_nw$temp.var)


#rescale load 
adcp_nw$resc_ld <- rescale(adcp_nw$load, to=c(0,1))



#Full GLMER model, binomial distribution. Number survived to emergence (success) vs number that failed to emerge (failure) 
#as response variable
#Mean temperature, temperature fluctuation and rescaled load as fixed effects. Random intercept of individual
wsem_nw_re_mod1 <- glmer(cbind(num.em, num.unem) ~ temp.avg * temp.var * resc_ld + (1|bug.id),
                           family=binomial,
                           data=adcp_nw,
                           na.action=na.omit,
                           control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

anova(wsem_nw_re_mod1)
summary(wsem_nw_re_mod1)



#model selection using dredge() 

#dredge requires dataframe with no NAs--subsetting to only columns in the model
adcp_wsmod <- adcp_nw %>% select(bug.id, temp.avg, temp.var, resc_ld, num.em, num.unem)
adcp_wsmod <- drop_na(adcp_wsmod)

#model with data frame with no NAs, na.action set to na.fail
wsem_nw_re_mod_d <- glmer(cbind(num.em, num.unem) ~ temp.avg * temp.var * resc_ld + (1|bug.id),
                         family=binomial,
                         data=adcp_wsmod,
                         na.action=na.fail,
                         control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000)))


wsmod_dredge <- dredge(wsem_nw_re_mod_d)
wsmod_dredge


#reduced model from dredge--the next best model has almost the same terms as with the original data, and the difference in AIC
#is very small. Ask Joel whether it could be appropriate to use the model with comparable terms
wsem_nw_mod_rd <- glmer(cbind(num.em, num.unem) ~ temp.avg + temp.var + resc_ld + temp.avg:temp.var 
                        + temp.avg:resc_ld + temp.var:resc_ld + (1|bug.id),
                        family=binomial,
                        data=adcp_nw,
                        na.action=na.omit,
                        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

anova(wsem_nw_mod_rd)
summary(wsem_nw_mod_rd)



#compare models without terms of interest to model from dredge. 

wsem_nw_mod_ta <- glmer(cbind(num.em, num.unem) ~ temp.var + resc_ld + temp.var:resc_ld + (1|bug.id),
                        family=binomial,
                        data=adcp_nw,
                        na.action=na.omit,
                        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

wsem_nw_mod_tv <- glmer(cbind(num.em, num.unem) ~ temp.avg + resc_ld + temp.avg:resc_ld + (1|bug.id),
                        family=binomial,
                        data=adcp_nw,
                        na.action=na.omit,
                        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

wsem_nw_mod_ld <- glmer(cbind(num.em, num.unem) ~ temp.avg + temp.var + temp.avg:temp.var + (1|bug.id),
                        family=binomial,
                        data=adcp_nw,
                        na.action=na.omit,
                        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

wsem_nw_mod_tatv <- glmer(cbind(num.em, num.unem) ~ temp.avg + temp.var + resc_ld
                        + temp.avg:resc_ld + temp.var:resc_ld + (1|bug.id),
                        family=binomial,
                        data=adcp_nw,
                        na.action=na.omit,
                        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

wsem_nw_mod_tald <- glmer(cbind(num.em, num.unem) ~ temp.avg + temp.var + resc_ld + temp.avg:temp.var 
                        + temp.var:resc_ld + (1|bug.id),
                        family=binomial,
                        data=adcp_nw,
                        na.action=na.omit,
                        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))


wsem_nw_mod_tvld <- glmer(cbind(num.em, num.unem) ~ temp.avg + temp.var + resc_ld + temp.avg:temp.var 
                        + temp.avg:resc_ld + (1|bug.id),
                        family=binomial,
                        data=adcp_nw,
                        na.action=na.omit,
                        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))


anova(wsem_nw_mod_rd, wsem_nw_mod_ta, test="Chisq")
anova(wsem_nw_mod_rd, wsem_nw_mod_tv, test="Chisq")
anova(wsem_nw_mod_rd, wsem_nw_mod_ld, test="Chisq")
anova(wsem_nw_mod_rd, wsem_nw_mod_tatv, test="Chisq")
anova(wsem_nw_mod_rd, wsem_nw_mod_tald, test="Chisq")
anova(wsem_nw_mod_rd, wsem_nw_mod_tvld, test="Chisq")



#Plot mean wasp survival to emergence
#mean proportion surviving to emergence
psem_sum <- summarySE(adcp, measurevar = "ps.em",
                      groupvars = c("temp.avg", "temp.var"),
                      na.rm = TRUE)
psem_sum


#making temp.avg numeric instead of a factor
psem_sum$temp.avg <- as.numeric(psem_sum$temp.avg)
psem_sum$temp.avg <- ifelse(psem_sum$temp.avg==1, 25,
                            ifelse(psem_sum$temp.avg==2, 28, 30))



#plot mean proportion surviving to emergence, with mean temperature on the x axis, proportion emerged on y axis
#color by fluctuation. Error bars = SE
#saved as individual object for combining with prop ecl for full figure
mn_psem_plot <- ggplot(psem_sum, aes(x=temp.avg, y=ps.em, color=as.character(temp.var)))
mn_psem_plot <- mn_psem_plot + geom_point(size=6, shape=2, stroke=2
) + geom_line(size = 2, linetype="dashed"
) + geom_errorbar(aes(ymin = ps.em-se, ymax = ps.em+se),
                  width=.5, size=1.2
) + scale_color_manual(values=c("#56B4E9","#D55E00"),name=c("Fluctuation [C]"),
                       breaks=c("0","10"),labels=c("0","10"),
                       guide=guide_legend(keywidth = 6, keyheight = 1.5)
) + scale_x_continuous(limits=c(24.5,30.5),
                       breaks = c(25, 28, 30)
) + scale_y_continuous(limits = c(0, .95),
                       breaks = c(0, 0.2, 0.4, 0.6, 0.8)
) + labs(x="Mean Temperature [C]", y="Prop. Emergence"
) + theme(text = element_text(family=("Cambria")),
          axis.line.x=element_line(colour = 'black', size = 1),
          axis.line.y=element_line(colour = 'black', size = 1),
          axis.ticks = element_line(colour = 'black', size = 1),
          axis.ticks.length = unit(2, "mm"),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          legend.background = element_rect(color="black",linetype="solid"),
          legend.position = "none")

mn_psem_plot



#---------------------
#analysis of wasp dev time

wdev_em_mod1 <- lm(ttem.w ~ temp.avg * temp.var * load,
                   data = adcp_nw,
                   na.action = na.omit)

anova(wdev_em_mod1)


#model selection using dredge() 

#dredge requires dataframe with no NAs--subsetting to only columns in the model
adcp_wdevem_mod <- adcp_nw %>% select(bug.id, temp.avg, temp.var, load, ttem.w)
adcp_wdevem_mod <- drop_na(adcp_wdevem_mod)

#model with data frame with no NAs, na.action set to na.fail
wdev_em_mod <- lm(ttem.w ~ temp.avg * temp.var * load,
                  data = adcp_wdevem_mod,
                  na.action = na.fail)

#dredged models
wdev_em_dredge <- dredge(wdev_em_mod)
wdev_em_dredge


#model from dredge
wdev_em_mod_rd <- lm(ttem.w ~ temp.avg + temp.var + load +temp.avg:temp.var
                     + temp.avg:load + temp.var:load,
                   data = adcp_nw,
                   na.action = na.omit)

anova(wdev_em_mod_rd)




#---------------------------------

#switching best fit models for each data set -- WASP SURVIVAL TO EMERGENCE

#TV data, wasp survival to em

#create combo column with temp.avg and temp.var
tvor_p <- unite(tvor_p, "tatv", temp.avg, temp.var, remove = FALSE)

tvor_nw <- subset(tvor_p, tatv!="30_10")

#make temp var a character instead of factor
tvor_nw$temp.var <- as.character(tvor_nw$temp.var)


#rescale load 
tvor_nw$resc_ld <- rescale(tvor_nw$load, to=c(0,1))



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
wsem_nw_mod_tv <- glmer(cbind(num.em, num.unem) ~ temp.avg + temp.var + resc_ld + temp.avg:temp.var 
                        + (1|bug.id),
                        family=binomial,
                        data=tvor_nw,
                        na.action=na.omit,
                        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))





#Best fit model from CPT data, fit to TV only data
wsem_nw_comp_mod_tv <- glmer(cbind(num.em, num.unem) ~ temp.avg + temp.var + resc_ld + temp.avg:temp.var 
                          + temp.avg:resc_ld + temp.var:resc_ld + (1|bug.id),
                          family=binomial,
                          data=tvor_nw,
                          na.action=na.omit,
                          control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))


#compare models
anova(wsem_nw_mod_tv, wsem_nw_comp_mod_tv, test="Chisq")




#CPT data
#reduced model from dredge
wsem_nw_mod_cpt <- glmer(cbind(num.em, num.unem) ~ temp.avg + temp.var + resc_ld + temp.avg:temp.var 
                        + temp.avg:resc_ld + temp.var:resc_ld + (1|bug.id),
                        family=binomial,
                        data=adcp_nw,
                        na.action=na.omit,
                        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))



#best fit model of TV data
wsem_nw_mod_cpt_comp <- glmer(cbind(num.em, num.unem) ~ temp.avg + temp.var + resc_ld + temp.avg:temp.var 
                              + (1|bug.id),
                              family=binomial,
                              data=adcp_nw,
                              na.action=na.omit,
                              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

#compare models
anova(wsem_nw_mod_cpt, wsem_nw_mod_cpt_comp, test="Chisq")



#---------------------------------

#switching best fit models for each data set -- WASP DEVELOPMENT TO EMERGENCE


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
wdev_em_mod_tv <- lm(ttem.w ~ temp.avg + temp.var + load + temp.avg:temp.var + temp.avg:load,
                   data = tvor_nw,
                   na.action = na.omit)



#best fit model from CPT data
wdev_em_mod_tv_comp <- lm(ttem.w ~ temp.avg + temp.var + load +temp.avg:temp.var
                     + temp.avg:load + temp.var:load,
                     data = tvor_nw,
                     na.action = na.omit)


anova(wdev_em_mod_tv, wdev_em_mod_tv_comp, test="Chisq")
AIC(wdev_em_mod_tv, wdev_em_mod_tv_comp)

logLik(wdev_em_mod_tv)
logLik(wdev_em_mod_tv_comp)




#analysis of wasp dev time--CPT DATA

wdev_em_mod_cpt <- lm(ttem.w ~ temp.avg * temp.var * load,
                   data = adcp_nw,
                   na.action = na.omit)

anova(wdev_em_mod_cpt)


#model selection using dredge() 

#dredge requires dataframe with no NAs--subsetting to only columns in the model
adcp_wdevem_mod <- adcp_nw %>% select(bug.id, temp.avg, temp.var, load, ttem.w)
adcp_wdevem_mod <- drop_na(adcp_wdevem_mod)

#model with data frame with no NAs, na.action set to na.fail
wdev_em_mod <- lm(ttem.w ~ temp.avg * temp.var * load,
                  data = adcp_wdevem_mod,
                  na.action = na.fail)

#dredged models
wdev_em_dredge <- dredge(wdev_em_mod)
wdev_em_dredge


#model from dredge
wdev_em_mod_cpt <- lm(ttem.w ~ temp.avg + temp.var + load +temp.avg:temp.var
                     + temp.avg:load + temp.var:load,
                     data = adcp_nw,
                     na.action = na.omit)


#best fit model from TV data
wdev_em_mod_cpt_comp <- lm(ttem.w ~ temp.avg + temp.var + load + temp.avg:temp.var + temp.avg:load,
                           data = adcp_nw,
                           na.action = na.omit)


#compare models
anova(wdev_em_mod_cpt, wdev_em_mod_cpt_comp, test="Chisq")
AIC(wdev_em_mod_cpt, wdev_em_mod_cpt_comp)

logLik(wdev_em_mod_cpt)
logLik(wdev_em_mod_cpt_comp)
