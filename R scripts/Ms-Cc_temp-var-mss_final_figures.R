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

#load data

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

#Configure data for use in creatng plots

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

#x axis temp avg or temp var? eclosion / total load, or ecl / num em?


#----------------------------

#WASP DEVELOPMENT TIME TO EMERGENCE AND ECLOSION

#x axis temp avg or temp var?


#-----------------------------

#DISTRIBUTION OF LIFE SPAN IN 30 MEAN TEMP TREATMENTS

#frequency dist plot, or polygon histogram (lines)?



#-------------------------

#POSSIBLE SUPPLEMENTAL FIGURES
  ##EFFECTS OF LOAD ON SURVIVAL TO EMERGENCE AND ECLOSION
  ##EFFECTS OF LOAD ON HOST GROWTH AND DEVELOPMENT














