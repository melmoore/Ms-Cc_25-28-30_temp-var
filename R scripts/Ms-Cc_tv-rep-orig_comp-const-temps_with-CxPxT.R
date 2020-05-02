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


#remove all temp treatments except 25C const and 30C const
tvor <- subset(tvor, temp.avg!=28 & temp.var!=10)
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
ad_psub <- subset(ad_p, temp.avg!=20)

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

wsem_mod1 <- glm(cbind(num.em, num.unem) ~ temp.avg * load * expt,
                 family = binomial,
                 data = ad_psub,
                 na.action=na.omit)
anova(wsem_mod1)
summary(wsem_mod1)


#testing significance of fixed effects

wsem_mod_ta <- glm(cbind(num.em, num.unem) ~ temp.avg,
                   family = binomial,
                   data = ad_psub,
                   na.action=na.omit)

wsem_mod_ld <- glm(cbind(num.em, num.unem) ~ load,
                   family = binomial,
                   data = ad_psub,
                   na.action=na.omit)


wsem_mod_expt <- glm(cbind(num.em, num.unem) ~ expt,
                     family = binomial,
                     data = ad_psub,
                     na.action=na.omit)


wsem_mod_tald <- glm(cbind(num.em, num.unem) ~ temp.avg:load,
                     family = binomial,
                     data = ad_psub,
                     na.action=na.omit)


wsem_mod_taex <- glm(cbind(num.em, num.unem) ~ temp.avg:expt,
                     family = binomial,
                     data = ad_psub,
                     na.action=na.omit)



wsem_mod_ldex <- glm(cbind(num.em, num.unem) ~ load:expt,
                     family = binomial,
                     data = ad_psub,
                     na.action=na.omit)


wsem_mod_taldex <- glm(cbind(num.em, num.unem) ~ temp.avg:load:expt,
                       family = binomial,
                       data = ad_psub,
                       na.action=na.omit)


wsem_mod_null <- glm(cbind(num.em, num.unem) ~ 1,
                     family = binomial,
                     data = ad_psub,
                     na.action=na.omit)

anova(wsem_mod_null, wsem_re_mod1, wsem_mod_ta, wsem_mod_ld, wsem_mod_expt, 
      wsem_mod_tald, wsem_mod_taex, wsem_mod_ldex, wsem_re_mod_taldex, test="Chisq")




#------------------------

#overdispersed, so adding a random effect of individual

#rescaling load
ad_psub$resc_ld <- rescale(ad_psub$load, to=c(0,1))

#model has a warning about large eigenvalues, something is off here
wsem_re_mod1 <- glmer(cbind(num.em, num.unem) ~ temp.avg * resc_ld * expt + (1|bug.id),
                 family = binomial,
                 data = ad_psub,
                 na.action=na.omit,
                 control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

anova(wsem_re_mod1)
summary(wsem_re_mod1)


#check singularity (following this guide :https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html)
tt <- getME(wsem_re_mod1,"theta")
ll <- getME(wsem_re_mod1,"lower")
min(tt[ll==0])  #not close to 0, so probably not a problem here

#checking gradients (from same guide)
#scaled gradient (??)
derivs1 <- wsem_re_mod1@optinfo$derivs
sc_grad1 <- with(derivs1,solve(Hessian,gradient))
max(abs(sc_grad1))

#absolute gradient (??)
max(pmin(abs(sc_grad1),abs(derivs1$gradient)))


#obviously this model isn't super great, but since it's just a quick check I'm going to leave it well 
#enough alone for now


#test significance of fixed effects
wsem_re_mod_ta <- glmer(cbind(num.em, num.unem) ~ temp.avg + (1|bug.id),
                        family = binomial,
                        data = ad_psub,
                        na.action=na.omit,
                        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))


wsem_re_mod_ld <- glmer(cbind(num.em, num.unem) ~ resc_ld + (1|bug.id),
                        family = binomial,
                        data = ad_psub,
                        na.action=na.omit,
                        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))


wsem_re_mod_expt <- glmer(cbind(num.em, num.unem) ~  expt + (1|bug.id),
                          family = binomial,
                          data = ad_psub,
                          na.action=na.omit,
                          control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))


wsem_re_mod_tald <- glmer(cbind(num.em, num.unem) ~ temp.avg:resc_ld + (1|bug.id),
                          family = binomial,
                          data = ad_psub,
                          na.action=na.omit,
                          control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))


wsem_re_mod_taex <- glmer(cbind(num.em, num.unem) ~ temp.avg:expt + (1|bug.id),
                          family = binomial,
                          data = ad_psub,
                          na.action=na.omit,
                          control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))


wsem_re_mod_ldex <- glmer(cbind(num.em, num.unem) ~ resc_ld:expt + (1|bug.id),
                          family = binomial,
                          data = ad_psub,
                          na.action=na.omit,
                          control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))


wsem_re_mod_taldex <- glmer(cbind(num.em, num.unem) ~ temp.avg:resc_ld:expt + (1|bug.id),
                            family = binomial,
                            data = ad_psub,
                            na.action=na.omit,
                            control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))


wsem_re_mod_null <- glmer(cbind(num.em, num.unem) ~ 1 + (1|bug.id),
                          family = binomial,
                          data = ad_psub,
                          na.action=na.omit,
                          control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))


anova(wsem_re_mod_null, wsem_re_mod1, wsem_re_mod_ta, wsem_re_mod_ld, wsem_re_mod_expt, 
      wsem_re_mod_tald, wsem_re_mod_taex, wsem_re_mod_ldex, wsem_re_mod_taldex, test="Chisq")



#--------------------

#Joel noted two points in cpt data that have much lower num.em at 30 than the rest of the data, suggested 
#removing and seeing how that affected slope/difference between tv data

which(ad_psub$load > 100 & ad_psub$num.em < 25)

ad_psub2 <- ad_psub[-c(10, 43), ]


#plot of num em by load, colored by expt, faceted by temp--two outliers removed
com_wsnmem_plot2 <- ggplot(ad_psub2, aes(x=load, y=num.em, group=expt, color=expt))
com_wsnmem_plot2 + geom_point(size=6
) + geom_smooth(method = "lm"
) + scale_color_manual(values=c("black", "#E69F00")
) + facet_wrap(~temp.avg)


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
                 family = binomial,
                 data = ad_psub2,
                 na.action=na.omit)
anova(wsem_mod1)
summary(wsem_mod1)



#----------------------

#effects on wasp development time

#plot of time to emergence by load, colored by expt, faceted by temp--two outliers removed
com_wsttem_plot2 <- ggplot(ad_psub2, aes(x=load, y=ttem.w, group=expt, color=expt))
com_wsttem_plot2 + geom_point(size=6
) + geom_smooth(method = "lm"
) + scale_color_manual(values=c("black", "#E69F00")
) + facet_wrap(~temp.avg)



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

#plot mean mass and mean age for each temp.avg, treatment and expt
mn_ma_plot <- ggplot(lmss_sum, aes(x=age, y=log_mss, group=expt, color=expt))
mn_ma_plot + geom_point(size=6
) + geom_line(size=2
) + geom_errorbar(aes(ymin=log_mss - se, ymax=log_mss + se),
                  width=.7, size=1.2
) + geom_errorbarh(aes(xmin=age - age_se, xmax=age + age_se),
                   height=.5, size=1.2
) + facet_wrap(temp.avg~treatment)



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
