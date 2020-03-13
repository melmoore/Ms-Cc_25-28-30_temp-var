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



#plot mean mass by mean age
mn_lma_plot <- ggplot(mn_lmss_sum, aes(x=age, y=log_mss, group=expt, color=expt))
mn_lma_plot + geom_point(
)+geom_line(
)+facet_wrap(temp.var~treatment)








