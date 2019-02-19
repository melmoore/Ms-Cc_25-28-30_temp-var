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
library(viridis)


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


#----------------------------------

#plotting caterpillar mass and age without field caterpillar data (for comparison)
  ##doesn't seem to make that much difference, decreases sample sizes in 30 treatments
  ##keep in mind for analyses and for manuscript--may just remove for simplicity

#PLOTTING HOST CATERPILLAR MASS X AGE (+/- 5 group included)


amass.wof.sum<-summarySE(tv.long.wof,measurevar = "log.mass",
                     groupvar=c("treatment","temp.avg","temp.var","instar"),
                     na.rm=TRUE)

amass.wof.sum


#removing end day.age for those without wasp emergence at 30.10.p
tv.long.nomwof<-tv.long.wof
tv.long.nomwof$suc.ovp[is.na(tv.long.nomwof$suc.ovp)]<-0
tv.long.nomwof$end.use<-ifelse(tv.long.nomwof$instar=="end" & tv.long.nomwof$all.trt=="30.10.p" & tv.long.nomwof$suc.ovp=="0",0,1)
tv.long.nomwof$suc.ovp[(tv.long.nomwof$suc.ovp)==0]<-NA
tv.long.nomwof<-subset(tv.long.nomwof,end.use==1)


dagem.wof.sum<-summarySE(tv.long.nomwof,measurevar = "day.age",
                     groupvar=c("treatment","temp.avg","temp.var","instar"),
                     na.rm=TRUE)


dagem.wof.sum


amass.wof.sum$day.age<-dagem.wof.sum[,6]
amass.wof.sum$dage.se<-dagem.wof.sum[,8]


amass.wof.plot<-ggplot(amass.wof.sum,aes(x=day.age,y=log.mass,
                                 group=interaction(temp.var,treatment),color=temp.var))
amass.wof.plot+geom_point(aes(shape=treatment),size=5
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



#-----------------------

#plotting wasp survival to emergence and eclosion by load

tv.para.wof<-subset(tv.para, pop=="lab")

tv.para.wof$percem<-tv.para.wof$num.em/tv.para.wof$load

tv.para.wof$percem[is.nan(tv.para.wof$percem)]<-NA

emload.plot<-ggplot(tv.para.wof, aes(x=load, y=percem, color=temp.var))
emload.plot+geom_point(
)+geom_smooth(method=lm, se=FALSE
)+facet_wrap(~temp.avg)

#survival to eclosion
eclload.plot<-ggplot(tv.para.wof, aes(x=load, y=tot.surv, color=temp.var))
eclload.plot+geom_point(
)+geom_smooth(method=lm, se=FALSE
)+facet_wrap(~temp.avg)


#---------------------

#plotting load effects on host mass by binning load into bins of 50 (?)

tv.lp.wof<-subset(tv.long, treatment=="para" & pop=="lab")

#subset out the unsuc ovp
tv.lp.wof<-subset(tv.lp.wof, suc.ovp==1)

range(tv.lp.wof$load)

#Making a column "bin" that puts hosts into bins determined by load, binned by 50

tv.lp.wof$bin<-ifelse(tv.lp.wof$load<=50 & tv.lp.wof$load!=0, 50,
                      ifelse(tv.lp.wof$load>50 & tv.lp.wof$load<=100, 100,
                             ifelse(tv.lp.wof$load>100 & tv.lp.wof$load<=150, 150,
                                    ifelse(tv.lp.wof$load>150 & tv.lp.wof$load<=200, 200,
                                           ifelse(tv.lp.wof$load>200 & tv.lp.wof$load<=250, 250,
                                                  ifelse(tv.lp.wof$load>250 & tv.lp.wof$load<=350, 300, 0))))))



#plotting raw data of mass by age, facet_warpped by temp.var, color by bin plotting each temp.avg separately

tv.lp.wof$bin<-as.factor(tv.lp.wof$bin)

tv.lp.wof25<-subset(tv.lp.wof, temp.avg==25)
tv.lp.wof28<-subset(tv.lp.wof, temp.avg==28)
tv.lp.wof30<-subset(tv.lp.wof, temp.avg==30)


rawma25.lbin.plot<-ggplot(tv.lp.wof25, aes(x=day.age, y=log.mass, group=interaction(bug.id, bin), 
                                           color=bin))
rawma25.lbin.plot+geom_point(
)+geom_line(
)+scale_color_viridis(discrete = TRUE
)+facet_wrap(~temp.var)


rawma28.lbin.plot<-ggplot(tv.lp.wof28, aes(x=day.age, y=log.mass, group=interaction(bug.id, bin),
                                           color=bin))
rawma28.lbin.plot+geom_point(
)+geom_line(
)+scale_color_viridis(discrete = TRUE
)+facet_wrap(~temp.var)


rawma30.lbin.plot<-ggplot(tv.lp.wof30, aes(x=day.age, y=log.mass, group=interaction(bug.id, bin),
                                           color=bin))
rawma30.lbin.plot+geom_point(
)+geom_line(
)+scale_color_viridis(discrete = TRUE
)+facet_wrap(~temp.var)




#plotting mean mass and age by bin

massbin.sum<-summarySE(tv.lp.wof, measurevar = "log.mass",
                       groupvars = c("temp.avg", "temp.var", "instar", "bin"),
                       na.rm = TRUE)
massbin.sum


agebin.sum<-summarySE(tv.lp.wof, measurevar = "day.age",
                      groupvars = c("temp.avg", "temp.var", "instar", "bin"),
                      na.rm = TRUE)
agebin.sum


massbin.sum$age<-agebin.sum[,6]
massbin.sum$age.se<-agebin.sum[,8]

massbin.sum25<-subset(massbin.sum, temp.avg==25)
massbin.sum28<-subset(massbin.sum, temp.avg==28)
massbin.sum30<-subset(massbin.sum, temp.avg==30)


mnbin25.plot<-ggplot(massbin.sum25, aes(x=age, y=log.mass, color=bin))
mnbin25.plot+geom_point(size=3
)+geom_line(size=1.2
)+geom_errorbar(aes(ymin=log.mass-se, ymax=log.mass+se),
                width=.5, size=1
)+geom_errorbarh(aes(xmin=age-age.se, xmax=age+age.se),
                 height=.5, size=1
)+facet_wrap(~temp.var)

mnbin28.plot<-ggplot(massbin.sum28, aes(x=age, y=log.mass, color=bin))
mnbin28.plot+geom_point(size=3
)+geom_line(size=1.2
)+geom_errorbar(aes(ymin=log.mass-se, ymax=log.mass+se),
                width=.5, size=1
)+geom_errorbarh(aes(xmin=age-age.se, xmax=age+age.se),
                 height=.5, size=1
)+facet_wrap(~temp.var)


mnbin30.plot<-ggplot(massbin.sum30, aes(x=age, y=log.mass, color=bin))
mnbin30.plot+geom_point(size=3
)+geom_line(size=1.2
)+geom_errorbar(aes(ymin=log.mass-se, ymax=log.mass+se),
                width=.5, size=1
)+geom_errorbarh(aes(xmin=age-age.se, xmax=age+age.se),
                 height=.5, size=1
)+facet_wrap(~temp.var)

