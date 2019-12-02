#Ms Cc temp var expt--analysis of wasp survival using brglm2


#load libraries
library(brglm2)
library(readr)


#load data
tv<-read_csv("data files/25-28-30_tv-final_clean.csv",
             col_types = cols(temp.avg = col_factor(levels = c("25","28", "30")), 
                              temp.var = col_factor(levels = c("0", "5", "10")), 
                              treatment = col_factor(levels = c("control","para"))))

View(tv)


tv.long <- read_csv("data files/25-28-30_tv-final_clean_LONG.csv",
                    col_types = cols(temp.avg = col_factor(levels = c("25", "28", "30")), 
                                     temp.var = col_factor(levels = c("0","5", "10")), 
                                     treatment = col_factor(levels = c("control","para"))))

View(tv.long)



#-------------------

#initial data wrangling

#Create datasets without the +/-5 treatment 
tv.long.no5<-subset(tv.long,temp.var!=5)
tv.no5<-subset(tv, temp.var!=5)


#create datasets with only parasitized hosts

#long dataframe
tv.long.para<-subset(tv.long.no5, treatment=="para")

#remove parasitized bugs that wandered
tv.long.para$date.wand.j[is.na(tv.long.para$date.wand.j)]<-0
tv.long.para<-subset(tv.long.para, date.wand.j==0)
tv.long.para$date.wand.j[tv.long.para$date.wand.j==0]<-NA


#wide dataframe
tv.para<-subset(tv.no5, treatment=="para")

#removing individuals that were parasitized but wandered
tv.para$date.wand.j[is.na(tv.para$date.wand.j)]<-0
tv.para<-subset(tv.para, date.wand.j==0)


#create a tot_died column for the "fail" category of the binomial model
tv.para$tot.died <- tv.para$load - tv.para$num.ecl


#------------------------------

#binomial glm model of wasp survival to eclosion

wsecl_mod1 <- glm(cbind(num.ecl, tot.died) ~ temp.avg*temp.var*load,
                  family = binomial,
                  data=tv.para,
                  na.action=na.omit)

summary(wsecl_mod1)


#checking for infinite estimates
inf_est_check <- check_infinite_estimates(wsecl_mod1)


#detect separation in in my model parameters
wsecl_sep <- glm(cbind(num.ecl, tot.died) ~ temp.avg*temp.var*load,
                 family = binomial,
                 method = "detect_separation",
                 data=tv.para,
                 na.action=na.omit)

wsecl_sep


#-----------------------

#attempting to use the brglm_fit function 

wsecl_br_mod1 <- glm(cbind(num.ecl, tot.died) ~ temp.avg*temp.var*load,
                     family = binomial,
                     method = "brglmFit",
                     type="AS_mixed",
                     data=tv.para,
                     na.action=na.omit)
