#Ms Cc temp var replication (30C temp avg)--WOWE prelim figures and analysis

#load libraries
library(Rmisc)
library(readr)
library(nlme)
library(lme4)
library(ggplot2)
library(tidyr)
library(mgcv)
library(dplyr)
library(viridis)
library(cowplot)



#load data--REPLICATION ONLY
wower <- read_csv("data files/Ms-Cc_tv-rep-30_WOWES.csv")
View(wower)

wower_lng <- read_csv("data files/Ms-Cc_tv-rep-30_WOWES_lng.csv")
View(wower_lng)


#load data--orig and repl data combined, will need to subset to WOWE
wowero <- read_csv("data files/Ms-Cc_tv-orig-rep_comb_cl.csv")
View(wowero)

wowero_lng <- read_csv("data files/Ms-Cc_tv-orig-rep_comb_lng.csv")
View(wowero_lng)


#----------------------

#subset orig and repl comb data sets to only WOWEs (including only those at 30.10 for simplicity--N at
#other treatments extremely small (like 2 individuals))
wowero <- subset(wowero, temp.avg==30 & temp.var==10)
wowero <- subset(wowero, treatment=="para")

wowero_lng <- subset(wowero_lng, temp.avg==30 & temp.var==10)
wowero_lng <- subset(wowero_lng, treatment=="para")



#-----------------------

#plot mass end by ttend (for those that have both)
ma_end_plot <- ggplot(wowero, aes(x=ttend, y=mass.end, color=end.class))
ma_end_plot + geom_point(
) + geom_smooth(method="lm"
) + facet_wrap(~expt)



#boxplot of mass end by end.class
mass_end_boxplot <- ggplot(wowero, aes(x=end.class, y=mass.end, group=interaction(expt, end.class),
                                       fill=expt))
mass_end_boxplot + geom_boxplot(position = "dodge")



#distribution of ttend--subset to only cull or dead cull
ttend_dist_plot <- ggplot(subset(wowero, end.class=="cull" | end.class=="dead_cull"), 
                          aes(x=ttend, fill=expt))
ttend_dist_plot + geom_density(adjust=2/3, alpha=.5)















