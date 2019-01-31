#Scrap notes from Ms+Cc_25-28-30_analysis script


ttend.plot<-ggplot(tv,aes(x=temp.avg,y=ttend,group=interaction(temp.avg,temp.var),color=temp.var))
ttend.plot+geom_jitter()



#attempting to see how the model differs if mongo caterpillars are removed

tv$date.em.j[is.na(tv$date.em.j)]<-0

tv$mongo<-ifelse(tv$treatment=="para" & tv$date.em.j==0, 1,0)

tv$mongo<-as.factor(tv$mongo)

tv.nom<-subset(tv,mongo==0)


ttend.plot2<-ggplot(tv,aes(x=temp.avg,y=ttend,group=interaction(temp.avg,temp.var),color=mongo))
ttend.plot2+geom_jitter(aes(shape=temp.var))


tv.long.no5$date.em.j[is.na(tv.long.no5$date.em.j)]<-0
tv.long.no5$mongo<-ifelse(tv.long.no5$treatment=="para" & tv.long.no5$date.em.j==0, 1,0)

tv.lng.no5nom<-subset(tv.long.no5,mongo==0)

#Model with no mongos:

lms.mod2<-lme(log.mass~(day.age+I(day.age^2)):(temp.var+treatment+temp.avg)^2+temp.avg,
              random=~day.age|bug.id,
              data=tv.lng.no5nom,
              na.action=na.omit,
              method="ML")

anova(lms.mod2)
summary(lms.mod2)




#getting rid of one unsuc.ovp that got missed in cleaning (fix this later)

tv.para<-subset(tv.para,bug.id!="30.5_p_75")




#Attempting to test the fixed effects of wasp total survival model


wtots.mod0<-glmer(cbind(tot.died,num.ecl)~1+(1|bug.id),
                  family=binomial,
                  data=tv.para,
                  na.action = na.omit,
                  control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))



wtots.mod2a<-glmer(cbind(tot.died,num.ecl)~temp.avg+(1|bug.id),
                   family=binomial,
                   data=tv.para,
                   na.action = na.omit,
                   control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))



wtots.mod2b<-glmer(cbind(tot.died,num.ecl)~temp.var.num+(1|bug.id),
                   family=binomial,
                   data=tv.para,
                   na.action = na.omit,
                   control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))



wtots.mod2c<-glmer(cbind(tot.died,num.ecl)~resc.ld+(1|bug.id),
                   family=binomial,
                   data=tv.para,
                   na.action = na.omit,
                   control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))



anova(wtots.mod0,wtots.mod2,wtots.mod2a,wtots.mod2b,wtots.mod2c)

#I think this means that both temp.var and temp.avg significantly affect the model, but resc.ld doesn't

#Attempting to use lmerTest package functions:
##Not working, look at later

lmerTest::anova(wtots.mod2,ddf="Kenward-Roger")




#testing wasp dev int model fixed effects

wdint.mod1a<-lme(waspdev.int~temp.avg,
                 random=~1|bug.id,
                 data=tv.para,
                 na.action = na.omit,
                 method="ML")


wdint.mod1b<-lme(waspdev.int~temp.var,
                 random=~1|bug.id,
                 data=tv.para,
                 na.action = na.omit,
                 method="ML")

wdint.mod1c<-lme(waspdev.int~resc.ld,
                 random=~1|bug.id,
                 data=tv.para,
                 na.action = na.omit,
                 method="ML")


anova(wdint.mod1,wdint.mod1a,wdint.mod1b,wdint.mod1c)  #best is full model




#Testing wasp dev tot model fixed effects

wdtot.mod1a<-lme(waspdev.tot~temp.avg,
                 random=~1|bug.id,
                 data=tv.para,
                 na.action=na.omit,
                 method="ML")

wdtot.mod1b<-lme(waspdev.tot~temp.var.num,
                 random=~1|bug.id,
                 data=tv.para,
                 na.action=na.omit,
                 method="ML")

wdtot.mod1c<-lme(waspdev.tot~resc.ld,
                 random=~1|bug.id,
                 data=tv.para,
                 na.action=na.omit,
                 method="ML")


anova(wdtot.mod1,wdtot.mod1a,wdtot.mod1b,wdtot.mod1c) #full model is best





