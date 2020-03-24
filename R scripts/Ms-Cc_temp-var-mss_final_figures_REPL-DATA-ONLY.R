#PLOTS FOR Ms Cc TEMPERATURE VARIATION MANUSCRIPT
##FINAL VERSIONS--USING REPL 30 DATA ONLY, NO ORIG 30 DATA



#load libraries

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

loadfonts(device="win")
font_import()


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




#--------------------------------

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


#set plot theme
theme_set(theme_classic())

#----------------------------------

#MEAN MASS BY MEAN AGE FOR MANDUCA SEXTA

#Calcualte means and variance of mass and age data
lmass_sum <- summarySE(tvor_lng, measurevar = "log_mss",
                       groupvars = c("temp.avg", "temp.var", "treatment", "instar"),
                       na.rm = TRUE)

lmass_sum


age_sum <- summarySE(tvor_lng, measurevar = "age",
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
                         size=5
) + geom_line(aes(linetype=treatment),
              size=1.7
) + geom_errorbar(aes(ymin=log_mss - se, ymax=log_mss + se),
                  width=1.2, size=1
) + geom_errorbarh(aes(xmin=age - age_se, xmax=age + age_se),
                   height=.4, size=1
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
        legend.position = c(.85, .25))




#-----------------------


#WASP SURVIVAL TO EMERGENCE AND ECLOSION:

#create subset with only parasitized individuals
tvor_p <- subset(tvor, treatment=="para")


#mean proportion surviving to emergence
psem_sum <- summarySE(tvor_p, measurevar = "ps.em",
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
mn_psem_plot <- ggplot(psem_sum, aes(x=temp.avg, y=ps.em, color=temp.var))
mn_psem_plot <- mn_psem_plot + geom_point(size=6
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


#Mean proportion survivng to eclosion

psecl_sum <- summarySE(tvor_p, measurevar = "ps.ecl",
                       groupvars = c("temp.avg", "temp.var"),
                       na.rm = TRUE)
psecl_sum


#making temp.avg numeric instead of a factor
psecl_sum$temp.avg <- as.numeric(psecl_sum$temp.avg)
psecl_sum$temp.avg <- ifelse(psecl_sum$temp.avg==1, 25,
                             ifelse(psecl_sum$temp.avg==2, 28, 30))



#Plot of mean proportion survivng to eclosion, with mean temperature on x axis, proportion eclosed on y axis
#color by fluctuation treatment
#saved as separate object for combining with prop em for full figure
mn_psecl_plot <- ggplot(psecl_sum, aes(x=temp.avg, y=ps.ecl, color=temp.var))
mn_psecl_plot <- mn_psecl_plot + geom_point(size=6
) + geom_line(size = 2
) + geom_errorbar(aes(ymin = ps.ecl-se, ymax = ps.ecl+se),
                  width=.5, size=1.2
) + scale_color_manual(values=c("#56B4E9","#D55E00"),name=c("Fluctuation [C]"),
                       breaks=c("0","10"),labels=c("0","10"),
                       guide=guide_legend(keywidth = 2.5)
) + scale_x_continuous(limits=c(24.5,30.5),
                       breaks = c(25, 28, 30)
) + scale_y_continuous(limits = c(0, 0.9),
                       breaks = c(0, 0.2, 0.4, 0.6, 0.8)) + labs(x="Mean Temperature [C]", y="Prop. Eclosion"
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
          legend.position = c(.8, .85))

mn_psecl_plot


#combine into one figure using cowplot
surv_fig <- plot_grid(mn_psem_plot, mn_psecl_plot, labels=c("A", "B"), align = "h")
surv_fig




#----------------------

#WASP DEVELOPMENT TIME TO EMERGENCE AND ECLOSION

#find mean development time to emergence
ttemw_sum <- summarySE(tvor_p, measurevar = "ttem.w",
                       groupvars = c("temp.avg", "temp.var"),
                       na.rm = TRUE)
ttemw_sum


#making temp.avg numeric instead of a factor
ttemw_sum$temp.avg <- as.numeric(ttemw_sum$temp.avg)
ttemw_sum$temp.avg <- ifelse(ttemw_sum$temp.avg==1, 25,
                             ifelse(ttemw_sum$temp.avg==2, 28, 30))

#Plotting development time to emergence with mean temperature on x axis, fluctuation as grouping variable
mn_ttem_plot <- ggplot(ttemw_sum, aes(x=temp.avg, y=ttem.w, color=temp.var))
mn_ttem_plot <- mn_ttem_plot + geom_point(size=5, shape=17
) + geom_line(size=2
) + geom_errorbar(aes(ymin=ttem.w - se, ymax=ttem.w + se),
                width=.5, size=1.2
) + scale_color_manual(values=c("#56B4E9","#D55E00"),name=c("Fluctuation [C]"),
                     breaks=c("0","10"),labels=c("0","10"),
                     guide=guide_legend(keywidth = 2.5)
) + scale_x_continuous(limits=c(24.5,30.5),
                     breaks = c(25, 28, 30)
) + scale_y_continuous(limits = c(10, 24),
                     breaks = c(10, 12, 14, 16, 18, 20, 22, 24)
) + labs(x="Mean Temperature [C]", y="Time to Emergence [days]"
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

mn_ttem_plot


#find mean development time to eclosion
ttecl_sum <- summarySE(tvor_p, measurevar = "ttecl",
                       groupvars = c("temp.avg", "temp.var"),
                       na.rm = TRUE)
ttecl_sum


#making temp.avg numeric instead of a factor
ttecl_sum$temp.avg <- as.numeric(ttecl_sum$temp.avg)
ttecl_sum$temp.avg <- ifelse(ttecl_sum$temp.avg==1, 25,
                             ifelse(ttecl_sum$temp.avg==2, 28, 30))


#Plotting development time to eclosion with mean temperature on x axis, fluctuation as grouping variable
mn_ttecl_plot <- ggplot(ttecl_sum, aes(x=temp.avg, y=ttecl, color=temp.var))
mn_ttecl_plot <- mn_ttecl_plot + geom_point(size=5, shape=17
) + geom_line(size=2
) + geom_errorbar(aes(ymin=ttecl - se, ymax=ttecl + se),
                  width=.5, size=1.2
) + scale_color_manual(values=c("#56B4E9","#D55E00"),name=c("Fluctuation [C]"),
                       breaks=c("0","10"),labels=c("0","10"),
                       guide=guide_legend(keywidth = 2.5)
) + scale_x_continuous(limits=c(24.5,30.5),
                       breaks = c(25, 28, 30)
) + scale_y_continuous(limits = c(10, 24),
                       breaks = c(10, 12, 14, 16, 18, 20, 22, 24)
) + labs(x="Mean Temperature [C]", y="Time to Eclosion [days]"
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
          legend.text = element_text(size=15),
          legend.title = element_text(size=18),
          legend.position = c(0.85, 0.15))

mn_ttecl_plot


#combine into one figure
dev_fig<-plot_grid(mn_ttem_plot, mn_ttecl_plot, labels=c("A", "B"), align="h")
dev_fig




#-----------------------------------

#DISTRIBUTION OF LIFE SPAN IN 30 MEAN TEMP TREATMENTS



#--------------------------------

#POSSIBLE SUPPLEMENTAL FIGURES
##EFFECTS OF LOAD ON SURVIVAL TO EMERGENCE AND ECLOSION
##EFFECTS OF LOAD ON HOST GROWTH AND DEVELOPMENT


#LOAD EFFECTS ON SURVIVAL TO EMERGENCE (NUMBER EMERGED)

#Color by mean temperature, facet by fluctuations
numem_ld_plot <- ggplot(tvor_p, aes(x=load, y=num.em, color=temp.avg))
numem_ld_plot + geom_point(size=6, alpha=.8
) + geom_smooth(method="lm", se=FALSE, size=2
) + scale_color_manual(values=c("#009E73","#E69F00","#000000"),name=c("Avg. Temp. [C]"),
                       breaks=c("25","28","30"),labels=c("25","28","30"),
                       guide=guide_legend(keywidth=3)       
) + facet_wrap(~temp.var
) + labs(x="Total Load", y="Number Emerged"
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
          legend.position = c(.9, .2))




#LOAD EFFECTS ON SURVIVAL TO ECLOSION (NUMBER ECLOSED)

#Color by mean temperature, facet by fluctuations
numecl_ld_plot <- ggplot(tvor_p, aes(x=load, y=num.ecl, color=temp.avg))
numecl_ld_plot + geom_point(size=6, alpha=.8
) + geom_smooth(method="lm", se=FALSE, size=2
) + scale_color_manual(values=c("#009E73","#E69F00","#000000"),name=c("Avg. Temp. [C]"),
                       breaks=c("25","28","30"),labels=c("25","28","30"),
                       guide=guide_legend(keywidth=3)       
) + facet_wrap(~temp.var
) + labs(x="Total Load", y="Number Eclosed"
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
          legend.position = c(.9, .2))



