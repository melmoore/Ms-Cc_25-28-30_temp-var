#PLOTS FOR Ms Cc TEMPERATURE VARIATION MANUSCRIPT
  ##FINAL VERSIONS


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


#-----------------------

#load data sets

tv <- read_csv("data files/25-28-30_tv-final_clean.csv", 
               col_types = cols(temp.avg = col_factor(levels = c("25", "28", "30")), 
                                temp.var = col_factor(levels = c("0", "5", "10")), 
                                treatment = col_factor(levels = c("control", "para"))))

View(tv)



tv.long <- read_csv("data files/25-28-30_tv-final_clean_LONG.csv", 
                    col_types = cols(temp.avg = col_factor(levels = c("25", "28", "30")), 
                                     temp.var = col_factor(levels = c("0", "5", "10")), 
                                     treatment = col_factor(levels = c("control", "para"))))

View(tv.long)


#-------------------------

#Configure data for use in creating plots

#creating a column for log.mass in tv.long

tv.long$log.mass<-log(tv.long$mass)


#There is one WOWE that has an erroneous date in the date.em.j column--I assume that it's an error, 
## and should be the date died, but am removing that individual for now
tv.long<-subset(tv.long, bug.id!="30.10_p_17")

tv<-subset(tv, bug.id!="30.10_p_17")

#Creating a few columns with treatment combinations to be able to sort out the 30+/-5 treatments

#Long dataset
tv.long<-unite(tv.long,tmp.trt,temp.avg,temp.var,sep=".",remove = FALSE)
tv.long<-unite(tv.long,all.trt,tmp.trt,treatment,sep=".",remove=FALSE)
tv.long$all.trt<-gsub(".para",".p",tv.long$all.trt)
tv.long$all.trt<-gsub(".control",".c",tv.long$all.trt)


tv.long.no5<-subset(tv.long,all.trt!="30.5.p")
tv.long.no5<-subset(tv.long.no5,all.trt!="30.5.c")

#Making a subset without field caterpillars
tv.long.wof<-subset(tv.long, pop=="lab")


#Wide dataset
tv<-unite(tv,tmp.trt,temp.avg,temp.var,sep=".",remove = FALSE)
tv<-unite(tv,all.trt,tmp.trt,treatment,sep=".",remove=FALSE)
tv$all.trt<-gsub(".para",".p",tv$all.trt)
tv$all.trt<-gsub(".control",".c",tv$all.trt)


tv.ab<-subset(tv,all.trt!="30.5.p")
tv.ab<-subset(tv.ab,all.trt!="30.5.c")

#Making a subset without field caterpillars
tv.wof<-subset(tv, pop=="lab")


#subset to only parasitized individuals for wasp figures
tv.para<-subset(tv.wof, treatment=="para")



#----------------------------------

#MEAN MASS BY MEAN AGE FOR MANDUCA SEXTA

#plotting caterpillar mass and age without field caterpillar data (for comparison)
##doesn't seem to make that much difference, decreases sample sizes in 30 treatments
##keep in mind for analyses and for manuscript--may just remove for simplicity

#PLOTTING HOST CATERPILLAR MASS X AGE (+/- 5 group excluded)

#remove individuals that were parasitized and wandered

tv.long.wof$keep.p<-ifelse(tv.long.wof$treatment=="para" & tv.long.wof$wander==1, 0, 1)

tv.long.wof<-subset(tv.long.wof, keep.p==1)


#Calcualte means and variance of mass and age data

#mean log(mass)
amass.wof.sum<-summarySE(tv.long.wof,measurevar = "log.mass",
                         groupvar=c("treatment","temp.avg","temp.var","instar"),
                         na.rm=TRUE)

amass.wof.sum


#removing end day.age for those without wasp emergence at 30.10.p--to calculate mean age only including 
  ##those with wasp emergence
tv.long.nomwof<-tv.long.wof
tv.long.nomwof$suc.ovp[is.na(tv.long.nomwof$suc.ovp)]<-0
tv.long.nomwof$end.use<-ifelse(tv.long.nomwof$instar=="end" & tv.long.nomwof$all.trt=="30.10.p" & tv.long.nomwof$suc.ovp=="0",0,1)
tv.long.nomwof$suc.ovp[(tv.long.nomwof$suc.ovp)==0]<-NA
tv.long.nomwof<-subset(tv.long.nomwof,end.use==1)

#calculate mean and variance of age 
dagem.wof.sum<-summarySE(tv.long.nomwof,measurevar = "day.age",
                         groupvar=c("treatment","temp.avg","temp.var","instar"),
                         na.rm=TRUE)


dagem.wof.sum

#add mean and standard error of age to mean mass data frame so they can be plotted together
amass.wof.sum$day.age<-dagem.wof.sum[,6]
amass.wof.sum$dage.se<-dagem.wof.sum[,8]

#remove the +/-5 treatment
amass.wofno5.sum<-subset(amass.wof.sum, temp.var!=5)
amass.wofno5.sum


#Plot mean mass by mean age, facetted by mean reaing temperature, grouped by rearing temperature variation
  ##and parasitization treatment. Color by temp var, line type and point shape by parasitization treatment
  ##error bars are standard error for mass (y axis) and age (x axis)

amass.wof.plot2<-ggplot(amass.wofno5.sum,aes(x=day.age,y=log.mass,
                                             group=interaction(temp.var,treatment),color=temp.var))
amass.wof.plot2+geom_point(aes(shape=treatment),size=5
)+geom_line(aes(linetype=treatment),size=1.7
)+geom_errorbar(aes(ymin=log.mass-se,ymax=log.mass+se),width=1.2,size=1
)+geom_errorbarh(aes(xmin=day.age-dage.se,xmax=day.age+dage.se), height=.4, size=1    
)+scale_color_manual(values=c("#56B4E9","#D55E00"),name=c("Fluctuation [C]"),
                     breaks=c("0","10"),labels=c("0","10"),
                     guide=guide_legend(keywidth = 2.5)
)+scale_linetype_manual(values=c("solid","dashed"),name="Treatment",
                        breaks=c("control","para"),labels=c("Cotrol","Parasitized"),
                        guide=guide_legend(keywidth = 2.5)
)+scale_shape_manual(values = c(16,17),name="Treatment",
                     breaks=c("control","para"),labels=c("Cotrol","Parasitized"),
                     guide=guide_legend(keywidth = 2.5)
)+labs(x="Age [days]",y="Log(Mass) [mg]"
)+facet_wrap(~temp.avg
)+theme(text = element_text(family=("Cambria")),
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

#----------------------------

#WASP SURVIVAL TO EMERGENCE AND ECLOSION

# % survival to emergence

percem.sum<-summarySE(tv.para,measurevar = "tot.elsurv",groupvars = c("temp.avg","temp.var"),na.rm=TRUE)
percem.sum

#plotting with mean temp on x axis, and color by temp var
#removing +/-5 temp treatment
percem.sum.no5<-subset(percem.sum, temp.var!=5)

#make mean temp numeric instead of factor
percem.sum.no5$temp.avg<-as.numeric(percem.sum.no5$temp.avg)
percem.sum.no5$temp.avg<-ifelse(percem.sum.no5$temp.avg==1, 25,
                                ifelse(percem.sum.no5$temp.avg==2, 28, 30))


emsurv.plot2<-ggplot(percem.sum.no5,aes(x=temp.avg,y=tot.elsurv,group=temp.var,color=temp.var))
percem_fig<-emsurv.plot2+geom_point(size=5,
                        shape=17
)+geom_line(size=2,
            linetype="dashed"
)+geom_errorbar(aes(ymin=tot.elsurv-se, ymax=tot.elsurv+se),
                width=.5, size=1.2
)+scale_color_manual(values=c("#56B4E9","#D55E00"),name=c("Fluctuation [C]"),
                     breaks=c("0","10"),labels=c("0","10"),
                     guide=guide_legend(keywidth = 2.5)
)+scale_x_continuous(limits=c(24.5,30.5),
                     breaks = c(25, 28, 30)
)+scale_y_continuous(limits = c(0, 0.9)
)+labs(x="Mean Temperature [C]", y="% Emergence"
)+theme(axis.line.x=element_line(colour = 'black', size = 1),
        axis.line.y=element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = 'black', size = 1),
        axis.ticks.length = unit(2, "mm"),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.background = element_rect(color="black",linetype="solid"),
        legend.position = "none")

percem_fig

# % survival to eclosion (divided by total load)

totsurv.sum<-summarySE(tv.para,measurevar = "tot.surv",
                       groupvars = c("temp.avg","temp.var"),
                       na.rm=TRUE)
totsurv.sum


#removing +/-5 treatment
totsurv.sum.no5<-subset(totsurv.sum, temp.var!=5)

#making temp avg be numeric
totsurv.sum.no5$temp.avg<-as.numeric(totsurv.sum.no5$temp.avg)
totsurv.sum.no5$temp.avg<-ifelse(totsurv.sum.no5$temp.avg==1, 25,
                                 ifelse(totsurv.sum.no5$temp.avg==2, 28, 30))

totsurv.plot2<-ggplot(totsurv.sum.no5,aes(x=temp.avg,y=tot.surv,group=temp.var,color=temp.var))
percecl_fig<-totsurv.plot2+geom_point(size=5,
                         shape=17
)+geom_line(size=2,
            linetype="dashed"
)+geom_errorbar(aes(ymin=tot.surv-se, ymax=tot.surv+se),
                width=.5, size=1.2
)+scale_color_manual(values=c("#56B4E9","#D55E00"),name=c("Fluctuation [C]"),
                     breaks=c("0","10"),labels=c("0","10"),
                     guide=guide_legend(keywidth = 2.5)
)+scale_x_continuous(limits=c(24.5,30.5),
                     breaks = c(25, 28, 30)
)+scale_y_continuous(limits = c(0, 0.9)
)+labs(x="Mean Temperature [C]", y="% Eclosion"
)+theme(axis.line.x=element_line(colour = 'black', size = 1),
        axis.line.y=element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = 'black', size = 1),
        axis.ticks.length = unit(2, "mm"),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.background = element_rect(color="black",linetype="solid"),
        legend.position = c(0.6, 0.85),
        legend.text = element_text(size=15),
        legend.title = element_text(size=18))

percecl_fig


#Combining perc_em and perc_ecl into one figure

surv_fig<-plot_grid(percem_fig, percecl_fig, labels=c("A", "B"), align="h")
surv_fig


#----------------------------

#WASP DEVELOPMENT TIME TO EMERGENCE AND ECLOSION

#Plotting dev time with temp avg on x axis, temp var as grouping variable

#Mean and variance of wasp development time to emergence
wdevint.sum<-summarySE(tv.para,measurevar = "waspdev.int",
                       groupvars = c("temp.avg","temp.var"),
                       na.rm=TRUE)
wdevint.sum

#removing +/-5 temp treatment
wdevint.sum.no5<-subset(wdevint.sum, temp.var!=5)

#make mean temp numeric instead of factor
wdevint.sum.no5$temp.avg<-as.numeric(wdevint.sum.no5$temp.avg)
wdevint.sum.no5$temp.avg<-ifelse(wdevint.sum.no5$temp.avg==1, 25,
                                 ifelse(wdevint.sum.no5$temp.avg==2, 28, 30))


wdevint.plot2<-ggplot(wdevint.sum.no5,aes(x=temp.avg,y=waspdev.int,group=temp.var,color=temp.var))
emdev_fig<-wdevint.plot2+geom_point(size=5, shape=17
)+geom_line(size=2, linetype="dashed"
)+geom_errorbar(aes(ymin=waspdev.int-se, ymax=waspdev.int+se),
                width=.5, size=1.2
)+scale_color_manual(values=c("#56B4E9","#D55E00"),name=c("Fluctuation [C]"),
                     breaks=c("0","10"),labels=c("0","10"),
                     guide=guide_legend(keywidth = 2.5)
)+scale_x_continuous(limits=c(24.5,30.5),
                     breaks = c(25, 28, 30)
)+scale_y_continuous(limits = c(10, 24),
                     breaks = c(10, 12, 14, 16, 18, 20, 22, 24)
)+labs(x="Mean Temperature [C]", y="Time to Emergence [days]"
)+theme(axis.line.x=element_line(colour = 'black', size = 1),
        axis.line.y=element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = 'black', size = 1),
        axis.ticks.length = unit(2, "mm"),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.background = element_rect(color="black",linetype="solid"),
        legend.position = "none")

emdev_fig


#Mean and variance of wasp development time to eclosion
wdevtot.sum<-summarySE(tv.para,measurevar = "waspdev.tot",
                       groupvars = c("temp.avg","temp.var"),
                       na.rm=TRUE)
wdevtot.sum

#removing +/-5 temp treatment
wdevtot.sum.no5<-subset(wdevtot.sum, temp.var!=5)

#make mean temp numeric instead of factor
wdevtot.sum.no5$temp.avg<-as.numeric(wdevtot.sum.no5$temp.avg)
wdevtot.sum.no5$temp.avg<-ifelse(wdevtot.sum.no5$temp.avg==1, 25,
                                 ifelse(wdevtot.sum.no5$temp.avg==2, 28, 30))


wdevtot.plot2<-ggplot(wdevtot.sum.no5,aes(x=temp.avg,y=waspdev.tot,group=temp.var,color=temp.var))
ecldev_fig<-wdevtot.plot2+geom_point(size=5,shape=17
)+geom_line(size=2,
            linetype="dashed"
)+geom_errorbar(aes(ymin=waspdev.tot-se, ymax=waspdev.tot+se),
                width=.5, size=1.2
)+scale_color_manual(values=c("#56B4E9","#D55E00"),name=c("Fluctuation [C]"),
                     breaks=c("0","10"),labels=c("0","10"),
                     guide=guide_legend(keywidth = 2.5)
)+scale_x_continuous(limits=c(24.5,30.5),
                     breaks = c(25, 28, 30)
)+scale_y_continuous(limits = c(10, 24),
                     breaks = c(10, 12, 14, 16, 18, 20, 22, 24)
)+labs(x="Mean Temperature [C]", y="Time to Eclosion [days]"
)+theme(axis.line.x=element_line(colour = 'black', size = 1),
        axis.line.y=element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = 'black', size = 1),
        axis.ticks.length = unit(2, "mm"),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.background = element_rect(color="black",linetype="solid"),
        legend.position = c(0.7, 0.15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=18))

ecldev_fig


#combine into one figure
dev_fig<-plot_grid(emdev_fig, ecldev_fig, labels=c("A", "B"), align="h")
dev_fig




#-----------------------------

#DISTRIBUTION OF LIFE SPAN IN 30 MEAN TEMP TREATMENTS

#subset to only 30 mean temp
high<-subset(tv.para, temp.avg==30)

#remove those that wandered
high<-subset(high, wander==0)

#subset to only costant or +/-10 treatments
high.no5<-subset(high, temp.var!=5)

#make a column for whether the individual had emergence or not
high.no5$date.em.j[is.na(high.no5$date.em.j)]<-0
high.no5$hadem<-ifelse(high.no5$date.em.j>0, 1, 0)

#need a column that combines temp.var and hadem
high.no5$emclass<-ifelse(high.no5$temp.var==0 & high.no5$hadem==1, "30_c_em",
                         ifelse(high.no5$temp.var==10 & high.no5$hadem==1, "30_f_em", "30_f_wowe"))


#histogram with density (no density curve)
mongo.hist.plot<-ggplot(high.no5, aes(x=ttend, fill=emclass))
mongo.hist.plot+geom_histogram(aes(y=..density..),
                                  binwidth = 1,
                                  position = "identity",
                                  col="black",
                                  alpha=.75
)+scale_fill_viridis(discrete = TRUE,
                     breaks=c("30_c_em", "30_f_em", "30_f_wowe"),
                     name="Treatment Outcome",
                     labels=c("30+/-0 emergence", "30+/10 emergence", "30+/-10 WOWE")
)+labs(x="Time [days]", y="Density"
)+theme(axis.line.x=element_line(colour = 'black', size = 1),
        axis.line.y=element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = 'black', size = 1),
        axis.ticks.length = unit(2, "mm"),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.background = element_rect(color="black",linetype="solid"),
        legend.position = c(0.65, 0.8),
        legend.text = element_text(size=15),
        legend.title = element_text(size=18))





#-------------------------

#POSSIBLE SUPPLEMENTAL FIGURES
  ##EFFECTS OF LOAD ON SURVIVAL TO EMERGENCE AND ECLOSION
  ##EFFECTS OF LOAD ON HOST GROWTH AND DEVELOPMENT














