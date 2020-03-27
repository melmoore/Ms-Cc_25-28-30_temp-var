#trying to visualize/test effects of diet and pop on Ms+Cc 25-28-30 expt


#load libraries

library(scales)
library(readr)
library(nlme)
library(lme4)
library(lmerTest)



#load data

tv <- read_csv("~/Manduca expts/Summer+Fall 2017/25.28.30 wasp temp var/data files/25-28-30_tv-final_clean.csv", 
               col_types = cols(temp.avg = col_factor(levels = c("25","28", "30")), 
                                temp.var = col_factor(levels = c("0", "5", "10")), 
                                treatment = col_factor(levels = c("control","para"))))

View(tv)



tv.long <- read_csv("~/Manduca expts/Summer+Fall 2017/25.28.30 wasp temp var/data files/25-28-30_tv-final_clean_LONG.csv", 
                    col_types = cols(temp.avg = col_factor(levels = c("25", "28", "30")), 
                                     temp.var = col_factor(levels = c("0","5", "10")), 
                                     treatment = col_factor(levels = c("control","para"))))

View(tv.long)



#creating a column for log.mass in tv.long

tv.long$log.mass<-log(tv.long$mass)






#quick test of effect of diet (colony or tobacco) on cat mass

diet.mod1<-lme(log.mass~pop+diet,
               random=~day.age|bug.id,
               data=tv.long,
               na.action=na.omit,
               method="ML")


anova(diet.mod1)


#subset for lab only

tv.lng.lb<-subset(tv.long,pop=="lab")

diet.mod2<-lme(log.mass~diet,
               random=~day.age|bug.id,
               data=tv.lng.lb,
               na.action=na.omit,
               method="ML")


anova(diet.mod2)
summary(diet.mod2)



diet.mod3<-lme(log.mass~diet*treatment,
               random=~day.age|bug.id,
               data=tv.lng.lb,
               na.action=na.omit,
               method="ML")

anova(diet.mod3)
summary(diet.mod3)






pop.mod1<-lme(log.mass~pop*treatment,
              random=~day.age|bug.id,
              data=tv.long,
              na.action=na.omit,
              method="ML")

anova(pop.mod1)





pd.plot<-ggplot(tv,aes(x=temp.avg,y=mass.5,
                       group=interaction(pop,temp.avg),
                       color=pop))
pd.plot+geom_boxplot()


pd.plot2<-ggplot(tv,aes(x=diet,y=mass.5,
                        group=interaction(temp.var,treatment),
                        color=temp.var,
                        fill=treatment))
pd.plot2+geom_boxplot(
       )+facet_wrap(~temp.avg)






#Size and diet are correlated, because only 30C is the only one fed tobacco diet, and they're smaller because of temperature
##ask Joel and James how to deal with this, as well as pop (also only in 30)



pd.sum<-summarySE(tv,measurevar = "ttend",
                  groupvars = c("temp.avg","temp.var","treatment","pop"),
                  na.rm=TRUE)


pd.sum


#looks like there are only a handful of field bugs, with the highest being 30.10.p grou--still only 14. Ask Joel if we should include, 
##or discard






