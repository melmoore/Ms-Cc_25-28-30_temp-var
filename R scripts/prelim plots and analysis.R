library(readr)
library(plyr)
library(ggplot2)
library(Rmisc)
library(dplyr)
library(tidyr)
library(reshape2)
library(cowplot)
library(nlme)


library(extrafont)

font_import()
loadfonts(device="win")

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

tv.long$log.mass<-log(tv.long$mass)





long.p<-subset(tv.long,treatment=="para")
long.c<-subset(tv.long,treatment=="control")


theme_set(theme_classic())


#Plotting age at each instar for each caterpillar--kind of messy, not sure if I will use this

age.plot<-ggplot(tv.long,aes(x=instar,y=day.age,group=interaction(bug.id, temp.avg),color=temp.avg))
age.plot+geom_point()+geom_line(aes(linetype=treatment))+geom_jitter()+facet_wrap(~temp.var)

age.plot<-ggplot(tv.long,aes(x=instar,y=day.age,group=interaction(bug.id, temp.avg),color=temp.avg))
age.plot+geom_line(aes(linetype=treatment))+facet_wrap(~temp.var)



#Making age at instar figures with average data

age.sum<-summarySE(tv.long,measurevar = "day.age",groupvars = c("temp.avg","temp.var","treatment","instar"),na.rm=TRUE)
age.sum

age.sum$temp.var<-as.factor(age.sum$temp.var)

age.sum.p<-subset(age.sum,treatment=="para")
age.sum.c<-subset(age.sum,treatment=="control")


ten.30<-subset(tv.long,treatment=="para" & temp.avg=="30" & temp.var=="10")
ten.30.end<-subset(ten.30,instar=="end")

old<-ten.30.end[,c("temp.var","instar","day.age","treatment","temp.avg")]

#Colors needs some work on this one--too many lines?
age.plot2<-ggplot(age.sum,aes(x=instar,y=day.age,group=interaction(temp.var,treatment), color=temp.var))
age.plot2+geom_point(aes(shape=treatment),size=4
        )+geom_line(aes(linetype=treatment),size=1.4
        )+geom_errorbar(aes(ymin=day.age-se,ymax=day.age+se),width=.3,size=1
        )+scale_color_manual(values=c("#56B4E9","#000000","#D55E00"),name=c("Fluctuation [C]"),
                             breaks=c("0","5","10"),labels=c("0","5","10")
        )+scale_linetype_manual(values=c("solid","dashed"),name="Treatment",
                                breaks=c("control","para"),labels=c("Cotrol","Parasitized")
        )+scale_shape_manual(values = c(16,17),name="Treatment",
                             breaks=c("control","para"),labels=c("Cotrol","Parasitized")
        )+labs(x="Instar",y="Age [days]"
        )+facet_wrap(~temp.avg
        )+theme(text = element_text(family=("Cambria")),
                strip.background = element_rect(colour="black",linetype = "solid",fill="white",
                                                size = 1),
                strip.text = element_text(size=16),
                axis.line.x=element_line(colour = 'black', size = 1),
                axis.line.y=element_line(colour = 'black', size = 1),
                axis.ticks = element_line(colour = 'black', size = 1),
                axis.ticks.length = unit(2, "mm"),
                axis.text.x = element_text(size = 14),
                axis.text.y = element_text(size = 14),
                axis.title.x = element_text(size = 16),
                axis.title.y = element_text(size = 16),
                legend.key.width=unit(15,"mm"),
                legend.background = element_rect(color="black",linetype="solid"))


#Day age by instar, red Xs for 30.10 group


age.plot3<-ggplot(age.sum,aes(x=instar,y=day.age,group=interaction(temp.var,treatment), color=temp.var))
age.plot3+geom_point(aes(shape=treatment),size=4
        )+geom_line(aes(linetype=treatment),size=1.5
        )+geom_errorbar(aes(ymin=day.age-se,ymax=day.age+se),width=.3,size=1
        )+scale_color_manual(values=c("#56B4E9","#000000","#D55E00"),name=c("Fluctuation [C]"),
                     breaks=c("0","5","10"),labels=c("0","5","10")
        )+scale_linetype_manual(values=c("solid","dashed"),name="Treatment",
                                breaks=c("control","para"),labels=c("Cotrol","Parasitized")
        )+scale_shape_manual(values = c(16,17),name="Treatment",
                             breaks=c("control","para"),labels=c("Cotrol","Parasitized")
        )+geom_jitter(data=old,aes(x=instar,y=day.age),color="red",shape=4,size=4             
        )+facet_wrap(~temp.avg
        )+labs(x="Instar",y="Age [days]"
        )+theme(text = element_text(family=("Cambria")),
                strip.background = element_rect(colour="black",linetype = "solid",fill="white",
                                                size = 1),
                strip.text = element_text(size=16),
                axis.line.x=element_line(colour = 'black', size = 1),
                axis.line.y=element_line(colour = 'black', size = 1),
                axis.ticks = element_line(colour = 'black', size = 1),
                axis.ticks.length = unit(2, "mm"),
                axis.text.x = element_text(size = 14),
                axis.text.y = element_text(size = 14),
                axis.title.x = element_text(size = 16),
                axis.title.y = element_text(size = 16),
                legend.key.width=unit(15,"mm"),
                legend.background = element_rect(color="black",linetype="solid"))


        
#Avg day age by instar, plotting raw points on top of average

age.plot3<-ggplot(age.sum,aes(x=instar,y=day.age,group=interaction(temp.var,treatment), color=temp.var))
age.plot3+geom_point(aes(shape=treatment),size=3)+geom_line(aes(linetype=treatment),size=1.2
        )+geom_errorbar(aes(ymin=day.age-se,ymax=day.age+se),width=.3
        )+scale_color_manual(values=c("#56B4E9","#000000","#D55E00"),name=c("Fluctuation [C]"),
                     breaks=c("0","5","10"),labels=c("0","5","10")
        )+scale_linetype_manual(values=c("solid","dashed"),name="Treatment",
                                breaks=c("control","para"),labels=c("Cotrol","Parasitized")
        )+scale_shape_manual(values = c(16,17),name="Treatment",
                             breaks=c("control","para"),labels=c("Cotrol","Parasitized")
        )+geom_jitter(data=tv.long,aes(x=instar,y=day.age,color=temp.var),shape=4,size=3,alpha=.4
        )+labs(x="Instar",y="Age [days]"
        )+facet_wrap(~temp.avg
        )+theme(text = element_text(family=("Cambria")),
                strip.background = element_rect(colour="black",linetype = "solid",fill="white",
                                        size = 1),
                strip.text = element_text(size=16),
                axis.line.x=element_line(colour = 'black', size = 1),
                axis.line.y=element_line(colour = 'black', size = 1),
                axis.ticks = element_line(colour = 'black', size = 1),
                axis.ticks.length = unit(2, "mm"),
                axis.text.x = element_text(size = 14),
                axis.text.y = element_text(size = 14),
                axis.title.x = element_text(size = 16),
                axis.title.y = element_text(size = 16),
                legend.key.width=unit(15,"mm"),
                legend.background = element_rect(color="black",linetype="solid"))





#Mass at instar plot with average data

#Logging mass

tv.long$log.mass<-log(tv.long$mass)

log.mass.sum<-summarySE(tv.long,measurevar = "log.mass",groupvars = c("temp.avg","temp.var","treatment","instar"),na.rm=TRUE)
log.mass.sum

#right now, have size of the point relative to sample size. Keep? Or have some other way to indicate low sample size in 30.10?
mass.plot<-ggplot(log.mass.sum,aes(x=instar,y=log.mass,group=interaction(temp.var,treatment), color=temp.var))
mass.plot+geom_point(aes(shape=treatment),size=4
)+geom_line(aes(linetype=treatment),size=1.4
)+geom_errorbar(aes(ymin=log.mass-se,ymax=log.mass+se),width=.3,size=1
)+scale_color_manual(values=c("#56B4E9","#000000","#D55E00"),name=c("Fluctuation [C]"),
                     breaks=c("0","5","10"),labels=c("0","5","10")
)+scale_linetype_manual(values=c("solid","twodash"),name="Treatment",
                        breaks=c("control","para"),labels=c("Cotrol","Parasitized")
)+scale_shape_manual(values = c(16,17),name="Treatment",
                     breaks=c("control","para"),labels=c("Cotrol","Parasitized")
)+labs(x="Instar",y="Log(Mass) [mg]"
)+facet_wrap(~temp.avg
)+theme(text = element_text(family=("Cambria")),
        strip.background = element_rect(colour="black",linetype = "solid",fill="white",
                                        size = 1),
        strip.text = element_text(size=16),
        axis.line.x=element_line(colour = 'black', size = 1),
        axis.line.y=element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = 'black', size = 1),
        axis.ticks.length = unit(2, "mm"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.key.width=unit(15,"mm"),
        legend.background = element_rect(color="black",linetype="solid"))



#Attempting mass by age plots--rather messy, haven't decided on best way to display. 


#Creating a few columns with treatment combinations to be able to sort out the 30.10 para group
tv.long<-unite(tv.long,tmp.trt,temp.avg,temp.var,sep=".",remove = FALSE)
tv.long<-unite(tv.long,all.trt,tmp.trt,treatment,sep=".",remove=FALSE)
tv.long$all.trt<-gsub(".para",".p",tv.long$all.trt)
tv.long$all.trt<-gsub(".control",".c",tv.long$all.trt)


tv.long.ab<-subset(tv.long,all.trt!="30.10.p")

tv.lab.p<-subset(tv.long.ab,treatment="para")
tv.lab.c<-subset(tv.long.ab,treatment=="control")


#Several plots with different colors, line types, and faceting options

am.plot<-ggplot(tv.long.ab,aes(x=day.age,y=log.mass,group=bug.id,color=treatment))
am.plot+geom_point()+geom_line(aes(linetype=temp.var,color=treatment)
      )+scale_color_manual(values=c("black","red"),breaks=c("control","para"),name=c("Treatment"),
                            labels=c("Control","Parasitized")
      )+facet_wrap(~temp.avg)


am.plot2<-ggplot(tv.long.ab,aes(x=day.age,y=log.mass,group=bug.id,color=temp.avg))
am.plot2+geom_point()+geom_line(aes(linetype=treatment,color=temp.avg)
       )+scale_color_manual(values=c("#009E73","#E69F00","#000000"),breaks=c("25","28","30"),name=c("Avg temp"),
                     labels=c("25","28","30")
       )+facet_wrap(~temp.var)



am.plot3<-ggplot(tv.long.ab,aes(x=day.age,y=log.mass,group=bug.id,color=temp.var))
am.plot3+geom_point()+geom_line(aes(linetype=treatment,color=temp.var)
)+scale_color_manual(values=c("#56B4E9","#000000","#D55E00"),breaks=c("0","5","10"),name=c("Fluctuation"),
                     labels=c("0","5","10")
)+facet_wrap(~temp.avg)


am.plot4<-ggplot(tv.long.ab,aes(x=day.age,y=log.mass,group=bug.id,color=treatment))
am.plot4+geom_point()+geom_line(aes(linetype=treatment,color=treatment)
)+scale_color_manual(values=c("black","red"),breaks=c("control","para"),name=c("Treatment"),
                     labels=c("Control","Parasitized")
)+facet_wrap(~temp.avg+temp.var)



am.plot5.c<-ggplot(tv.lab.c,aes(x=day.age,y=log.mass,group=bug.id,color=temp.var))
am.plot5.c+geom_point()+geom_line(aes(linetype=temp.var,color=temp.var)
       )+scale_color_manual(values=c("#56B4E9","#000000","#D55E00"),breaks=c("0","5","10"),
                            name=c("Fluctuation"),labels=c("0","5","10")
       )+facet_wrap(~temp.avg)


am.plot5.p<-ggplot(tv.lab.p,aes(x=day.age,y=log.mass,group=bug.id,color=temp.var))
am.plot5.p+geom_point()+geom_line(aes(linetype=temp.var,color=temp.var)
)+scale_color_manual(values=c("#56B4E9","#000000","#D55E00"),breaks=c("0","5","10"),
                     name=c("Fluctuation"),labels=c("0","5","10")
)+facet_wrap(~temp.avg)



#Trying to make average mass and age plot (minus 30.10.p)
#Find a way to include only 30.10.p that have final masses

amass.sum<-summarySE(tv.long.ab,measurevar = "log.mass",
                     groupvar=c("treatment","temp.avg","temp.var","instar"),
                     na.rm=TRUE)

amass.sum


dagem.sum<-summarySE(tv.long.ab,measurevar = "day.age",
                     groupvar=c("treatment","temp.avg","temp.var","instar"),
                     na.rm=TRUE)


dagem.sum


amass.sum$day.age<-dagem.sum[,6]
amass.sum$dage.se<-dagem.sum[,8]


colnames(amass.sum)


amass.plot<-ggplot(amass.sum,aes(x=day.age,y=log.mass,
                                 group=interaction(temp.var,treatment),color=temp.var))
amass.plot+geom_point(aes(shape=treatment),size=4
)+geom_line(aes(linetype=treatment),size=1.4
)+geom_errorbar(aes(ymin=log.mass-se,ymax=log.mass+se),width=.3,size=1
)+geom_errorbarh(aes(xmin=day.age-dage.se,xmax=day.age+dage.se)      
)+scale_color_manual(values=c("#56B4E9","#000000","#D55E00"),name=c("Fluctuation [C]"),
                     breaks=c("0","5","10"),labels=c("0","5","10")
)+scale_linetype_manual(values=c("solid","dashed"),name="Treatment",
                        breaks=c("control","para"),labels=c("Cotrol","Parasitized")
)+scale_shape_manual(values = c(16,17),name="Treatment",
                     breaks=c("control","para"),labels=c("Cotrol","Parasitized")
)+labs(x="Age [days]",y="Log(Mass) [mg]"
)+facet_wrap(~temp.avg)


#Wasp total survival to adulthood plot

tv.para<-subset(tv,treatment=="para")

totsurv.sum<-summarySE(tv.para,measurevar = "tot.surv",groupvars = c("temp.avg","temp.var"),na.rm=TRUE)
totsurv.sum


totsurv.plot<-ggplot(totsurv.sum,aes(x=temp.var,y=tot.surv,group=temp.avg,color=temp.avg))
totsurv.plot+geom_point(aes(shape=temp.avg),size=3)+geom_line(aes(linetype=temp.avg),size=1.2
                       )+geom_errorbar(aes(ymin=tot.surv-se,ymax=tot.surv+se),width=.5
                       )+scale_color_manual(values=c("#009E73","#E69F00","#000000"),name=c("Average Temperature [C]"),
                                            breaks=c("25","28","30"),labels=c("25","28","30")
                       )+scale_linetype_discrete(name=c("Average Temperature [C]"),
                                                 breaks=c("25","28","30"),labels=c("25","28","30")
                       )+scale_shape_manual(name="Average Temperature [C]",breaks=c("25","28","30"),values=c(15,16,17),
                                            labels=c("25","28","30"))
                       

#Survival to emergence rxn norm

percem.sum<-summarySE(tv.para,measurevar = "tot.elsurv",groupvars = c("temp.avg","temp.var"),na.rm=TRUE)
percem.sum


emsurv.plot<-ggplot(percem.sum,aes(x=temp.var,y=tot.elsurv,group=temp.avg,color=temp.avg))
emsurv.plot+geom_point(aes(shape=temp.avg),size=3)+geom_line(aes(linetype=temp.avg),size=1.2
)+geom_errorbar(aes(ymin=tot.elsurv-se,ymax=tot.elsurv+se),width=.5
)+scale_color_manual(values=c("#009E73","#E69F00","#000000"),name=c("Average Temperature [C]"),
                     breaks=c("25","28","30"),labels=c("25","28","30")
)+scale_linetype_discrete(name=c("Average Temperature [C]"),
                          breaks=c("25","28","30"),labels=c("25","28","30")
)+scale_shape_manual(name="Average Temperature [C]",breaks=c("25","28","30"),values=c(15,16,17),
                     labels=c("25","28","30"))




#Wasp development time graph

tv.para$ovp.day<-1

cc.long<-gather(tv.para,stage,devtime,waspdev.tot,waspdev.int,ovp.day)

cc.long$stage<-gsub("waspdev.tot", "eclosion",cc.long$stage)
cc.long$stage<-gsub("waspdev.int", "emergence",cc.long$stage)
cc.long$stage<-gsub("ovp.day", "oviposition",cc.long$stage)

cc.long$devtime[(cc.long$devtime=="1")]<-0

cc.long$stage<-factor(cc.long$stage,levels=c("oviposition","emergence","eclosion"))

dev.sum<-summarySE(cc.long, measurevar = "devtime", groupvars = c("temp.avg","temp.var","stage"),na.rm=TRUE)
dev.sum




waspdev.plot<-ggplot(dev.sum,aes(x=stage,y=devtime,group=temp.var,color=temp.var))
waspdev.plot+geom_point(size=4,shape=17
           )+geom_line(size=1.4
           )+geom_errorbar(aes(ymin=devtime-se,ymax=devtime+se),width=.5,size=1
           )+scale_color_manual(values=c("#56B4E9","#000000","#D55E00"),name=c("Fluctuation [C]"),
                     breaks=c("0","5","10"),labels=c("0","5","10")
           )+labs(x="Stage",y="Age [days]"
           )+facet_wrap(~temp.avg
           )+theme(text = element_text(family=("Cambria")),
                   strip.background = element_rect(colour="black",linetype = "solid",fill="white",
                                                   size = 1),
                   strip.text = element_text(size=16),
                   axis.line.x=element_line(colour = 'black', size = 1),
                   axis.line.y=element_line(colour = 'black', size = 1),
                   axis.ticks = element_line(colour = 'black', size = 1),
                   axis.ticks.length = unit(2, "mm"),
                   axis.text.x = element_text(size = 14),
                   axis.text.y = element_text(size = 14),
                   axis.title.x = element_text(size = 16),
                   axis.title.y = element_text(size = 16),
                   legend.key.width=unit(15,"mm"),
                   legend.background = element_rect(color="black",linetype="solid"))




waspdev.plot2<-ggplot(dev.sum,aes(x=temp.var,y=devtime,group=temp.avg,color=temp.avg))
waspdev.plot2+geom_point()+geom_line()+facet_wrap(~stage)



#Wasp dev time by load

wdev.plot<-ggplot(tv,aes(x=load,y=waspdev.tot,group=temp.var,color=temp.var))
wdev.plot+geom_point(
        )+geom_smooth(se=FALSE,method = "lm"
        )+facet_wrap(~temp.avg)



#Cat mass at em by load

massem.plot<-ggplot(tv,aes(x=load,y=mass.48em,group=temp.var,color=temp.var))
massem.plot+geom_point(
)+geom_smooth(se=FALSE,method = "lm"
)+facet_wrap(~temp.avg)





#Plotting "survivor" curve of wasp life stages.

#Making standardized wasp columns (standardized by dividing by load)

#removing individuals with unsuccessful ovp, while keeping 30.10 group
tv.para$suc.ovp[is.na(tv.para$suc.ovp)]<-0
tv.para$suc.ovp<-ifelse(tv.para$temp.avg=="30" & tv.para$temp.var=="10",.5,tv.para$suc.ovp)
tv.para<-subset(tv.para,suc.ovp!="0")
tv.para$suc.ovp[(tv.para$suc.ovp==".5")]<-NA


tv.para$stnd.load<-(tv.para$load/tv.para$load)
tv.para$stnd.load[is.nan(tv.para$stnd.load)]<-NA

tv.para$stnd.em<-(tv.para$num.em/tv.para$load)
tv.para$stnd.em[is.nan(tv.para$stnd.em)]<-NA

tv.para$stnd.coc<-(tv.para$num.coc/tv.para$load)
tv.para$stnd.coc[is.nan(tv.para$stnd.coc)]<-NA

tv.para$stnd.ecl<-(tv.para$num.ecl/tv.para$load)
tv.para$stnd.ecl[is.nan(tv.para$stnd.ecl)]<-NA


para.sub<-select(tv.para,bug.id,treatment,temp.avg,temp.var,stnd.load,stnd.em,stnd.coc,stnd.ecl)

para.sub$stnd.load[is.nan(para.sub$stnd.load)]<-0
para.sub$stnd.em[is.nan(para.sub$stnd.em)]<-0
para.sub$stnd.coc[is.nan(para.sub$stnd.coc)]<-0
para.sub$stnd.ecl[is.nan(para.sub$stnd.ecl)]<-0


para.sub<-gather(para.sub,stage,prop.surv,stnd.load,stnd.em,stnd.coc,stnd.ecl)

para.sub$stage<-factor(para.sub$stage,levels=c("stnd.load","stnd.em","stnd.coc","stnd.ecl"))

surv.plot<-ggplot(para.sub,aes(x=stage,y=prop.surv,group=interaction(bug.id,temp.var),color=temp.var))
surv.plot+geom_point()+geom_line()+facet_wrap(~temp.avg)




surv.sum<-summarySE(para.sub,measurevar = "prop.surv",groupvars = c("temp.avg","temp.var","stage"),na.rm=TRUE)
surv.sum

surv.plot2<-ggplot(surv.sum,aes(x=stage,y=prop.surv,group=temp.var,color=temp.var))
surv.plot2+geom_point(size=4,shape=17
         )+geom_line(size=1.4
         )+geom_errorbar(aes(ymin=prop.surv-se,ymax=prop.surv+se),width=.5,size=1
         )+scale_color_manual(values=c("#56B4E9","#000000","#D55E00"),name=c("Diurnal Fluctuation [C]"),
                     breaks=c("0","5","10"),labels=c("0","5","10")
         )+labs(x="Stage",y="Survival [%]"
         )+facet_wrap(~temp.avg
         )+theme(text = element_text(family=("Cambria")),
                 strip.background = element_rect(colour="black",linetype = "solid",fill="white",
                                                 size = 1),
                 strip.text = element_text(size=16),
                 axis.line.x=element_line(colour = 'black', size = 1),
                 axis.line.y=element_line(colour = 'black', size = 1),
                 axis.ticks = element_line(colour = 'black', size = 1),
                 axis.ticks.length = unit(2, "mm"),
                 axis.text.x = element_text(size = 14),
                 axis.text.y = element_text(size = 14),
                 axis.title.x = element_text(size = 16),
                 axis.title.y = element_text(size = 16),
                 legend.key.width=unit(15,"mm"),
                 legend.background = element_rect(color="black",linetype="solid"))








#Plotting a frequency distribution curve of 30.10 dev time to wandering and emergence

all.30<-subset(tv,temp.avg=="30")
all.30.c<-subset(all.30,treatment=="control")
all.30.p<-subset(all.30,treatment=="para")

fd.plot<-ggplot(all.30.c, aes(x=ttend,group=temp.var,color=temp.var)) 
fd.plot+geom_density()+scale_color_manual(values=c("#56B4E9","#000000","#D55E00"),name=c("Diurnal Fluctuation [C]"),
                                          breaks=c("0","5","10"),labels=c("0","5","10"))


fd.plot2<-ggplot(all.30.p, aes(x=ttend,group=temp.var,color=temp.var)) 
fd.plot2+geom_density()+scale_color_manual(values=c("#56B4E9","#000000","#D55E00"),name=c("Diurnal Fluctuation [C]"),
                                           breaks=c("0","5","10"),labels=c("0","5","10"))





#Plotting wasp survival by load

ldsurv.plot<-ggplot(tv.para,aes(x=load,y=tot.surv,group=temp.var,color=temp.var))
ldsurv.plot+geom_point()+geom_smooth(se=FALSE)+facet_wrap(~temp.avg)


#Plotting wasp dev by load

devload.plot<-ggplot(tv.para,aes(x=load,y=waspdev.int,group=temp.var,color=temp.var))
devload.plot+geom_point()+geom_jitter()+geom_smooth(se=FALSE)+facet_wrap(~temp.avg)


#try analysis with multinomial (like binomial, but with multiple possibilities)--but isn't contained within glmer






#Preliminary analysis


#log.mass (for Ms) model

tv.long$temp.var<-as.numeric(tv.long$temp.var)

#need to include mass at 3rd?  Does not seem to be running, get an error about singularity...maybe because not all temps
#have 5 fluc?--made temp.var continuous, and now it seems to run--is this ok?
mass.mod1<-lme(log.mass~treatment*temp.avg*temp.var,random=~1|bug.id,data=tv.long,method="ML",na.action=na.omit)
anova(mass.mod1)
summary(mass.mod1)


mass.mod1a<-lme(log.mass~treatment+temp.avg+temp.var+treatment:temp.var+temp.avg:temp.var+treatment:temp.avg:temp.var,
                random=~1|bug.id,data=tv.long,method="ML",na.action=na.omit)
anova(mass.mod1a)

anova(mass.mod1,mass.mod1a)


mass.mod1b<-lme(log.mass~treatment+temp.avg+temp.var+temp.avg:temp.var+treatment:temp.avg:temp.var,
                random=~1|bug.id,data=tv.long,method="ML",na.action=na.omit)
anova(mass.mod1b)

anova(mass.mod1,mass.mod1b)


mass.mod1c<-lme(log.mass~treatment+temp.avg+temp.var+temp.avg:temp.var,
                random=~1|bug.id,data=tv.long,method="ML",na.action=na.omit)
anova(mass.mod1c)

anova(mass.mod1,mass.mod1c)


mass.mod1d<-lme(log.mass~treatment+temp.avg+temp.var,
                random=~1|bug.id,data=tv.long,method="ML",na.action=na.omit)
anova(mass.mod1d)

anova(mass.mod1,mass.mod1d)


mass.mod1e<-lme(log.mass~treatment+temp.avg,
                random=~1|bug.id,data=tv.long,method="ML",na.action=na.omit)
anova(mass.mod1e)

anova(mass.mod1,mass.mod1e)




#prelim model for caterpillar dev time

age.mod1<-lme(day.age~treatment*temp.avg*temp.var,random=~1|bug.id,data=tv.long,method="ML",na.action=na.omit)
anova(age.mod1)
summary(age.mod1)




#Trying KU model:
#have to add temp.var as well
#Getting an error about singularity--it runs if you have temp.var as numeric--maybe because
#only 30 has a temp var of 5?

tv.long$temp.var.num<-as.numeric(tv.long$temp.var)

#Age as linear and quadratic in FE, only linear in RE

lm.mod1<-lme(log.mass~(day.age+I(day.age^2)):(temp.avg*temp.var.num*treatment)+temp.avg,
             random=~day.age|bug.id,data=tv.long,na.action=na.omit,method="ML")
anova(lm.mod1)
summary(lm.mod1)

VarCorr(lm.mod1)
          

#Age as only linear in FE and RE

lm.mod2<-lme(log.mass~day.age:(temp.avg*temp.var*treatment)+temp.avg,
             random=~day.age|bug.id,data=tv.long,na.action=na.omit,method="ML")
anova(lm.mod2)
summary(lm.mod2)

AIC(lm.mod1,lm.mod2) #Better with quadratic age term


#Age as linear and quadratic in both FE and RE--can't get this one to converge

lm.mod3<-lme(log.mass~(day.age+I(day.age^2)):(temp.avg*temp.var*treatment)+temp.avg,
             random=~day.age+I(day.age^2)|bug.id,data=tv.long,na.action=na.omit,method="ML",
             control=lmeControl(maxIter = 100000))





#Modelling age by temp.avg, temp.var and treatment

age.mod1<-lme(day.age~temp.avg*temp.var*treatment,random=~1|bug.id,data=tv.long,
              na.action=na.omit,method="ML")
anova(age.mod1)
summary(age.mod1)







ttend.sum<-summarySE(tv,measurevar = "ttend",groupvars = "all.trt",na.rm=TRUE)

ttend.sum


int5.sum<-summarySE(tv,measurevar = "int5",groupvars = "all.trt",na.rm=TRUE)

int5.sum


thrty<-subset(tv,temp.avg=="30")
thrty<-subset(thrty,temp.var!="5")


old.plot<-ggplot(thrty,aes(x=temp.var, y=ttend,fill=treatment))
old.plot+geom_boxplot(lwd=1.2
       )+scale_fill_manual(values = c("black", "red"),
                           breaks=c("control","para"),
                           labels=c("Control","Parasitized"),
                           name="Treatment"
       )+scale_x_discrete(labels=c("30+/-0","30+/-10")
       )+labs(x="Temperature",y="Age [days]"
       )+theme(text = element_text(family=("Cambria")),
               axis.line.x=element_line(colour = 'black', size = 1),
               axis.line.y=element_line(colour = 'black', size = 1),
               axis.ticks = element_line(colour = 'black', size = 1),
               axis.ticks.length = unit(2, "mm"),
               axis.text.x = element_text(size = 18),
               axis.text.y = element_text(size = 20),
               axis.title.x = element_text(size = 24),
               axis.title.y = element_text(size = 24),
               legend.key.width=unit(15,"mm"),
               legend.key.height = unit(10,"mm"),
               legend.text=element_text(size=12),
               legend.title=element_text(size=16),
               legend.position = c(.15,.8),
               legend.background = element_rect(color="black",linetype="solid",size=1))







#Attempting to plot age of 30.10.p on ageXmass graph

ten.30<-subset(tv.long,all.trt=="30.10.p")

old<-summarySE(ten.30,measurevar = "day.age",groupvar=c("treatment","temp.avg","temp.var","instar"),
               na.rm=TRUE)

old

old.mass<-summarySE(ten.30,measurevar = "log.mass",
                    groupvar=c("treatment","temp.avg","temp.var","instar"),
                    na.rm=TRUE)
old.mass


old$log.mass<-old.mass[,6]
old$lm.se<-old.mass[,8]


ten.30.end<-subset(ten.30,instar=="end")

ten.30.end<-ten.30.end[,c("temp.var","instar","day.age","treatment","temp.avg")]

ten.30.end$log.mass<-3.5



amass.plot2<-ggplot(amass.sum,aes(x=day.age,y=log.mass,
                                  group=interaction(temp.var,treatment),color=temp.var))
amass.plot2+geom_point(aes(shape=treatment),size=4
)+geom_line(aes(linetype=treatment),size=1.4
)+geom_errorbar(aes(ymin=log.mass-se,ymax=log.mass+se),width=.3,size=1
)+geom_errorbarh(aes(xmin=day.age-dage.se,xmax=day.age+dage.se)    
)+geom_point(data=old,aes(x=day.age,y=log.mass,color=temp.var),size=4, shape=17
)+geom_line(data=old,aes(x=day.age,y=log.mass),size=1.4, linetype="dashed"
)+geom_errorbar(data=old,aes(ymin=log.mass-lm.se,ymax=log.mass+lm.se),width=.3,size=1
)+geom_errorbarh(data=old,aes(xmin=day.age-se,xmax=day.age+se),height=.3,size=1
)+geom_jitter(data=ten.30.end,aes(x=day.age,y=log.mass),
              size=5,shape=18,color="red",width=.2,height=.1
)+scale_color_manual(values=c("#56B4E9","#D55E00"),name=c("Fluctuation [C]"),
                     breaks=c("0","10"),labels=c("0","10")
)+scale_linetype_manual(values=c("solid","dashed"),name="Treatment",
                        breaks=c("control","para"),labels=c("Cotrol","Parasitized")
)+scale_shape_manual(values = c(16,17),name="Treatment",
                     breaks=c("control","para"),labels=c("Cotrol","Parasitized")
)+labs(x="Age [days]",y="Log(Mass) [mg]"
)+facet_wrap(~temp.avg
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
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        legend.key.width=unit(15,"mm"),
        legend.background = element_rect(color="black",linetype="solid"))





#Plotting effects of load on survival to emergence

loadem.plot<-ggplot(tv.para, aes(x=load, y=num.em, group=temp.var, color=temp.var))
loadem.plot+geom_point(
)+geom_smooth(method="lm"
)+facet_wrap(~temp.avg)


loadecl.plot<-ggplot(tv.para, aes(x=load, y=num.ecl, group=temp.var, color=temp.var))
loadecl.plot+geom_point(
)+geom_smooth(method="lm"
)+facet_wrap(~temp.avg)


load.devem.plot<-ggplot(tv.para, aes(x=load, y=waspdev.int, group=temp.var, color=temp.var))
load.devem.plot+geom_point(
)+geom_smooth(method="lm"
)+facet_wrap(~temp.avg)


load.devecl.plot<-ggplot(tv.para, aes(x=load, y=waspdev.tot, group=temp.var, color=temp.var))
load.devecl.plot+geom_point(
)+geom_smooth(method="lm"
)+facet_wrap(~temp.avg)



load.percem.plot<-ggplot(tv.para, aes(x=load, y=tot.elsurv, group=temp.var, color=temp.var))
load.percem.plot+geom_point(
)+geom_smooth(method="lm"
)+facet_wrap(~temp.avg)


load.percecl.plot<-ggplot(tv.para, aes(x=load, y=tot.surv, group=temp.var, color=temp.var))
load.percecl.plot+geom_point(
)+geom_smooth(method="lm"
)+facet_wrap(~temp.avg)



#------------------------

#plotting mass at end of dev for all treats

theme_set(theme_classic())

lmend.sum<-summarySE(tv, measurevar = "lmend",
                     groupvars = c("temp.avg", "temp.var", "treatment"),
                     na.rm = TRUE)
lmend.sum


massend.plot<-ggplot(lmend.sum, aes(x=temp.var, y=lmend, group=interaction (temp.avg, treatment), color=temp.avg))
massend.plot+geom_point(size=4
           )+geom_line(aes(linetype=treatment),
                       size=1.2
           )+geom_errorbar(aes(ymin=lmend-se, ymax=lmend+se),
                           width=.4,
                           size=1)


lmend.sum.no5<-subset(lmend.sum, temp.var!=5)

massend.plot2<-ggplot(lmend.sum.no5, aes(x=temp.avg, y=lmend, group=interaction (temp.var, treatment), color=treatment))
massend.plot2+geom_point(size=4
)+geom_line(aes(linetype=temp.var),
            size=1.2
)+geom_errorbar(aes(ymin=lmend-ci, ymax=lmend+ci),
                width=.4,
                size=1)

#------------------------------
#comparing diff parts of my data to help make sense of model results

lmend.tatr.sum<-summarySE(tv, measurevar = "lmend",
                        groupvars = c("temp.avg", "treatment"),
                        na.rm = TRUE)
lmend.tatr.sum




lmend.tatv.sum<-summarySE(tv, measurevar = "lmend",
                          groupvars = c("temp.avg", "temp.var"),
                          na.rm = TRUE)
lmend.tatv.sum


lmend.tvtr.sum<-summarySE(tv, measurevar = "lmend",
                          groupvars = c("temp.var", "treatment"),
                          na.rm = TRUE)

tatr.plot<-ggplot(lmend.tatr.sum, aes(x=temp.avg, y=lmend, group=treatment, color=treatment))
tatr.plot+geom_point(
)+geom_line()


tatv.plot<-ggplot(lmend.tatv.sum, aes(x=temp.var, y=lmend, group=temp.avg, color=temp.avg))
tatv.plot+geom_point(
)+geom_line()


tvtr.plot<-ggplot(lmend.tvtr.sum, aes(x=temp.var, y=lmend, group=treatment, color=treatment))
tvtr.plot+geom_point(
)+geom_line()
