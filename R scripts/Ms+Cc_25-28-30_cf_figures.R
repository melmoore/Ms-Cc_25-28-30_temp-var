#Ms+Cc 25-28-30 const and fluc expt--FIGURES


#load libraries

library(readr)
library(plyr)
library(ggplot2)
library(Rmisc)
library(dplyr)
library(tidyr)
library(reshape2)
library(cowplot)


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



#creating a column for log.mass in tv.long

tv.long$log.mass<-log(tv.long$mass)




#Creating a few columns with treatment combinations to be able to sort out the 30+/-5 treatments

#Long dataset
tv.long<-unite(tv.long,tmp.trt,temp.avg,temp.var,sep=".",remove = FALSE)
tv.long<-unite(tv.long,all.trt,tmp.trt,treatment,sep=".",remove=FALSE)
tv.long$all.trt<-gsub(".para",".p",tv.long$all.trt)
tv.long$all.trt<-gsub(".control",".c",tv.long$all.trt)


tv.long.no5<-subset(tv.long,all.trt!="30.5.p")
tv.long.no5<-subset(tv.long.no5,all.trt!="30.5.c")



#Wide dataset
tv<-unite(tv,tmp.trt,temp.avg,temp.var,sep=".",remove = FALSE)
tv<-unite(tv,all.trt,tmp.trt,treatment,sep=".",remove=FALSE)
tv$all.trt<-gsub(".para",".p",tv$all.trt)
tv$all.trt<-gsub(".control",".c",tv$all.trt)


tv.ab<-subset(tv,all.trt!="30.5.p")
tv.ab<-subset(tv.ab,all.trt!="30.5.c")


#-----------------------------------

#PLOTTING HOST CATERPILLAR MASS X AGE (+/- 5 group included)


amass.sum<-summarySE(tv.long,measurevar = "log.mass",
                     groupvar=c("treatment","temp.avg","temp.var","instar"),
                     na.rm=TRUE)

amass.sum


#removing end day.age for those without wasp emergence at 30.10.p
tv.long.nom<-tv.long
tv.long.nom$suc.ovp[is.na(tv.long.nom$suc.ovp)]<-0
tv.long.nom$end.use<-ifelse(tv.long.nom$instar=="end" & tv.long.nom$all.trt=="30.10.p" & tv.long.nom$suc.ovp=="0",0,1)
tv.long.nom$suc.ovp[(tv.long.nom$suc.ovp)==0]<-NA
tv.long.nom<-subset(tv.long.nom,end.use==1)


dagem.sum<-summarySE(tv.long.nom,measurevar = "day.age",
                     groupvar=c("treatment","temp.avg","temp.var","instar"),
                     na.rm=TRUE)


dagem.sum


amass.sum$day.age<-dagem.sum[,6]
amass.sum$dage.se<-dagem.sum[,8]


amass.plot<-ggplot(amass.sum,aes(x=day.age,y=log.mass,
                                 group=interaction(temp.var,treatment),color=temp.var))
amass.plot+geom_point(aes(shape=treatment),size=5
)+geom_line(aes(linetype=treatment),size=1.7
)+geom_errorbar(aes(ymin=log.mass-se,ymax=log.mass+se),width=1.2,size=1
)+geom_errorbarh(aes(xmin=day.age-dage.se,xmax=day.age+dage.se), height=.4, size=1    
)+scale_color_manual(values=c("#56B4E9","#000000","#D55E00"),name=c("Fluctuation [C]"),
                     breaks=c("0","5","10"),labels=c("0","5","10"),
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





#---------------------------------------

#PLOTTING WASP SURVIVAL TO EMERGENCE AND ECLOSION

#Wasp total survival to adulthood plot

tv.para<-subset(tv,treatment=="para")

totsurv.sum<-summarySE(tv.para,measurevar = "tot.surv",groupvars = c("temp.avg","temp.var"),na.rm=TRUE)
totsurv.sum


totsurv.plot<-ggplot(totsurv.sum,aes(x=temp.var,y=tot.surv,group=temp.avg,color=temp.avg))
totsurv.plot<-totsurv.plot+geom_point(aes(shape=temp.avg),size=5
)+geom_line(aes(linetype=temp.avg),size=2
)+geom_errorbar(aes(ymin=tot.surv-se,ymax=tot.surv+se),width=.5, size=1.2
)+scale_color_manual(values=c("#009E73","#E69F00","#000000"),name=c("Avg. Temp. [C]"),
                     breaks=c("25","28","30"),labels=c("25","28","30"),
                     guide=guide_legend(keywidth=3)
)+scale_linetype_discrete(name=c("Avg. Temp. [C]"),
                          breaks=c("25","28","30"),labels=c("25","28","30"),
                          guide=guide_legend(keywidth=3)
)+scale_shape_manual(name="Avg. Temp. [C]",breaks=c("25","28","30"),values=c(15,16,17),
                     labels=c("25","28","30"),
                     guide=guide_legend(keywidth=3)
)+scale_y_continuous(limits=c(0,.9)
)+labs(x="Diurnal Fluctuation [C]", y="% Eclosion"
)+theme(text = element_text(family=("Cambria")),
        strip.background = element_rect(colour="black",linetype = "solid",fill="white",
                                        size = 1),
        strip.text = element_text(size=24),
        axis.line.x=element_line(colour = 'black', size = 1),
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

totsurv.plot


#Survival to emergence rxn norm

percem.sum<-summarySE(tv.para,measurevar = "tot.elsurv",groupvars = c("temp.avg","temp.var"),na.rm=TRUE)
percem.sum


emsurv.plot<-ggplot(percem.sum,aes(x=temp.var,y=tot.elsurv,group=temp.avg,color=temp.avg))
emsurv.plot<-emsurv.plot+geom_point(aes(shape=temp.avg),
                                    show.legend=FALSE,size=5
)+geom_line(aes(linetype=temp.avg),
            size=2, show.legend=FALSE
)+geom_errorbar(aes(ymin=tot.elsurv-se,ymax=tot.elsurv+se),width=.5, size=1.2
)+scale_color_manual(values=c("#009E73","#E69F00","#000000"),
                     breaks=c("25","28","30"),labels=c("25","28","30"),
                     guide=FALSE
)+scale_linetype_discrete(breaks=c("25","28","30"),labels=c("25","28","30")
)+scale_shape_manual(breaks=c("25","28","30"),values=c(15,16,17),
                     labels=c("25","28","30")
)+scale_y_continuous(limits=c(0,.90)
)+labs(x="Diurnal Fluctuation [C]", y="% Emergence"
)+theme(text = element_text(family=("Cambria")),
        strip.background = element_rect(colour="black",linetype = "solid",fill="white",
                                        size = 1),
        strip.text = element_text(size=24),
        axis.line.x=element_line(colour = 'black', size = 1),
        axis.line.y=element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = 'black', size = 1),
        axis.ticks.length = unit(2, "mm"),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.background = element_rect(color="black",linetype="solid"))

emsurv.plot


surv.fig<-plot_grid(emsurv.plot, totsurv.plot, labels=c("A", "B"), align="h")
surv.fig

#---------------------------------

#PLOTTING WASP DEVELOPMENT TIME TO EMERGENCE AND ECLOSION

#Wasp development time to eclosion

wdevtot.sum<-summarySE(tv.para,measurevar = "waspdev.tot",groupvars = c("temp.avg","temp.var"),na.rm=TRUE)
wdevtot.sum



wdevtot.plot<-ggplot(wdevtot.sum,aes(x=temp.var,y=waspdev.tot,group=temp.avg,color=temp.avg))
wdevtot.plot<-wdevtot.plot+geom_point(aes(shape=temp.avg),size=3)+geom_line(aes(linetype=temp.avg),size=1.2
)+geom_errorbar(aes(ymin=waspdev.tot-se,ymax=waspdev.tot+se),width=.5
)+scale_color_manual(values=c("#009E73","#E69F00","#000000"),name=c("Avg. Temperature [C]"),
                     breaks=c("25","28","30"),labels=c("25","28","30")
)+scale_linetype_discrete(name=c("Avg. Temperature [C]"),
                          breaks=c("25","28","30"),labels=c("25","28","30")
)+scale_shape_manual(name="Avg. Temperature [C]",breaks=c("25","28","30"),values=c(15,16,17),
                     labels=c("25","28","30")
)+scale_y_continuous(limits=c(8,25)
)+labs(x="Diurnal Fluctuation [C]", y="Days to Eclosion"
)+theme(legend.position = c(.6, .3))

wdevtot.plot



#Survival to emergence rxn norm

wdevint.sum<-summarySE(tv.para,measurevar = "waspdev.int",groupvars = c("temp.avg","temp.var"),na.rm=TRUE)
wdevint.sum


wdevint.plot<-ggplot(wdevint.sum,aes(x=temp.var,y=waspdev.int,group=temp.avg,color=temp.avg))
wdevint.plot<-wdevint.plot+geom_point(aes(shape=temp.avg),
                                      size=3, show.legend = FALSE
)+geom_line(aes(linetype=temp.avg),
            size=1.2, show.legend = FALSE
)+geom_errorbar(aes(ymin=waspdev.int-se,ymax=waspdev.int+se),width=.5
)+scale_color_manual(values=c("#009E73","#E69F00","#000000"),
                     breaks=c("25","28","30"),labels=c("25","28","30"),
                     guide=FALSE
)+scale_linetype_discrete(breaks=c("25","28","30"),labels=c("25","28","30")
)+scale_shape_manual(breaks=c("25","28","30"),values=c(15,16,17),
                     labels=c("25","28","30")
)+scale_y_continuous(limits=c(8,25)
)+labs(x="Diurnal Fluctuation [C]", y="Days to Emergence")

wdevint.plot


dev_fig<-plot_grid(wdevint.plot, wdevtot.plot, labels=c("A", "B"), align="h")
dev_fig




#------------------

#Creating a plot of conditional survival to eclosion--survival from emergence to eclosion
  ##Percent survival will be #ecl/#emerged

tv.para$emsecl<-tv.para$num.ecl/tv.para$num.em

emsecl.sum<-summarySE(tv.para, measurevar = "emsecl",
                      groupvars = c("temp.avg", "temp.var"),
                      na.rm = TRUE)


emsecl.plot<-ggplot(emsecl.sum,aes(x=temp.var,y=emsecl,group=temp.avg,color=temp.avg))
emsecl.plot+geom_point(aes(shape=temp.avg),size=5
)+geom_line(aes(linetype=temp.avg),size=2
)+geom_errorbar(aes(ymin=emsecl-se,ymax=emsecl+se),width=.5, size=1.2
)+scale_color_manual(values=c("#009E73","#E69F00","#000000"),name=c("Avg. Temp. [C]"),
                     breaks=c("25","28","30"),labels=c("25","28","30"),
                     guide=guide_legend(keywidth=3)
)+scale_linetype_discrete(name=c("Avg. Temp. [C]"),
                          breaks=c("25","28","30"),labels=c("25","28","30"),
                          guide=guide_legend(keywidth=3)
)+scale_shape_manual(name="Avg. Temp. [C]",breaks=c("25","28","30"),values=c(15,16,17),
                     labels=c("25","28","30"),
                     guide=guide_legend(keywidth=3)
)+scale_y_continuous(limits=c(0,.9)
)+labs(x="Diurnal Fluctuation [C]", y="% Eclosion"
)+theme(text = element_text(family=("Cambria")),
        strip.background = element_rect(colour="black",linetype = "solid",fill="white",
                                        size = 1),
        strip.text = element_text(size=24),
        axis.line.x=element_line(colour = 'black', size = 1),
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




