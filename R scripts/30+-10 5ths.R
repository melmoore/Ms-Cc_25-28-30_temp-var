library(readr)
library(plyr)
library(ggplot2)
library(Rmisc)
library(dplyr)
library(tidyr)
library(reshape2)


tv.ed <- read_csv("data files/incomplete data/25.28.30 tv 10.12 ed.csv", 
                  col_types = cols(temp.avg = col_factor(levels = c("25", "28", "30")), 
                                   temp.var = col_factor(levels = c("0", "5", "10")), 
                                   treatment = col_factor(levels = c("control", "para"))))


ten.30<-subset(tv.ed,treatment=="para" & temp.avg=="30" & temp.var=="10")


ten.30.raw <- read_csv("data files/incomplete data/30.10 5ths RAW.csv",
                       col_types = cols(temp.avg = col_factor(levels = c("25", "28", "30")), 
                                        temp.var = col_factor(levels = c("0", "5", "10")), 
                                        treat = col_factor(levels = c("control", "para"))))


View(ten.30.raw)

ten.30.raw<-rename(ten.30.raw,treatment=treat)
ten.30.raw<-select(ten.30.raw,-X20,-X21,-X22,-notes)



ten.30.raw$date.died.j<-strptime(ten.30.raw$date.died, "%m/%d")$yday+1
ten.30.raw$date.mass1.j<-strptime(ten.30.raw$date.mass.1, "%m/%d")$yday+1
ten.30.raw$date.5.j<-strptime(ten.30.raw$date.5, "%m/%d")$yday+1
ten.30.raw$date.6.j<-strptime(ten.30.raw$date.6, "%m/%d")$yday+1
ten.30.raw$date.mass2.j<-strptime(ten.30.raw$date.mass.2, "%m/%d")$yday+1



ten.30.fin<-merge(ten.30,ten.30.raw,by=c("bug.id"))



colnames(ten.30.fin)


ten.30.fin<-select(ten.30.fin,-treatment.y,-temp.avg.y,-temp.var.y,-date.died,-date.mass.1,-date.5,-date.6,
                   -date.mass.2,-date.died.j.y,-date.5.j.y)

ten.30.fin<-rename(ten.30.fin,temp.avg=temp.avg.x,temp.var=temp.var.x,treatment=treatment.x,date.died.j=date.died.j.x,date.5.j=date.5.j.x)

ten.30.fin$date.em.j[is.na(ten.30.fin$date.em.j)]<-0

ten.30.fin<-subset(ten.30.fin,date.em.j=="0")

ten.30.fin$date.em.j[(ten.30.fin$date.em.j=="0")]<-NA


ten.30.fin$daysince.5.0<-0


g1<-gather(ten.30.fin,measure,age.s5,int5,daysince.5.1,daysince.5.2,daysince.5.6,daysince.5.0)
g1$measure<-gsub("daysince.5.","",g1$measure)
g1$measure<-gsub("int5","end",g1$measure)


g2<-gather(ten.30.fin,measure,mass.s5,final.mass,mass.1,mass.2,mass.6,mass.5)
g2$measure<-gsub("mass.","",g2$measure)
g2$measure<-gsub("5","0",g2$measure)
g2$measure<-gsub("final.mass","end",g2$measure)

g1<-select(g1,bug.id,measure,age.s5)


ten.30.long<-merge(g1,g2,by=c("bug.id","measure"))
View(ten.30.long)

ten.30.long$after.10<-gsub("f","30.10",ten.30.long$after.10)

ten.30.long$mass.s5[is.na(ten.30.long$mass.s5)]<-0

no.mass<-subset(ten.30.long,measure=="end" & mass.s5=="0")

ten.30.long$mass.s5[(ten.30.long$mass.s5=="0")]<-NA


theme_set(theme_classic())


#Ask Joel about this plot--still not quite sure how best to present data

am.plot<-ggplot(ten.30.long,aes(x=age.s5,y=mass.s5,group=bug.id))
am.plot+geom_point(size=4
      )+geom_line(aes(linetype=after.10),size=1.4
      )+geom_point(data=no.mass,aes(x=age.s5,y=mass.s5),shape=17,color="red",size=5
      )+scale_linetype_manual(values=c("solid","dashed"),breaks=c("25.00","30.10"),
                              name="Test Temperature [C]",labels=c("25","30+/-10")
      )+labs(x="Age since molt to 5th [days]",y="Mass [mg]"
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




am.plot2<-ggplot(ten.30.long,aes(x=age.s5,y=mass.s5,group=bug.id))
am.plot2+geom_point(size=5)+geom_line(size=1.4
      )+geom_point(data=no.mass,aes(x=age.s5,y=mass.s5),shape=17,color="red",size=5
      )+facet_wrap(~after.10
      )+labs(x="Age since molt to 5th [days]",y="Mass [mg]"
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



                

















