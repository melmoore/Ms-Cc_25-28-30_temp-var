#Ms Cc temp var expt--analysis of wasp survival using brglm2


#load libraries
library(brglm2)
library(readr)
library(logistf)


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

View(inf_est_check)

#detect separation in in my model parameters
wsecl_sep <- glm(cbind(num.ecl, tot.died) ~ temp.avg*temp.var*load,
                 family = binomial,
                 method = "detect_separation",
                 data=tv.para,
                 na.action=na.omit)

wsecl_sep


#-----------------------

#Using logistf to model the proportion of wasps that survived to eclosion
  ##Implements Firth's bias-Reduced penalized-likelihood logistic regression.
  ##Don't know why this one runs when the num.ecl model will not...both are numeric
  ##the only thing I can think of is tot.surv is a proportion (and therefor between 0 and 1),
  ##where as num.ecl is an integer?

wpsecl_lf_mod1 <- logistf(tot.surv ~ temp.avg*temp.var*load,
                         data=tv.para)

summary(wpsecl_lf_mod1)
logistftest(wpsecl_lf_mod1)



#using logistf to model the number of wasps that survived to eclosion
  ##this does not run, I think because my response variable must be a vector with 0,1 or TRUE,FALSE values
  ##is this appropriate for my analysis, since my response variable would be a binary response of 
  ##wasp survival, i.e. the host had any wasp emergence or none. That's not really the question I want to ask
  ## Instead, I want to see how wasp survival is affected by temp--are there fewer wasps at higher temps?
  ##unless I'm missing something (very possible), I don't think a logistic regression using logistf will 
  ##answer the question I'm interested in?

tv.para$num.ecl <- as.numeric(tv.para$num.ecl)

wsnecl_lf_mod1 <- logistf(num.ecl ~ temp.avg*temp.var*load,
                         data=tv.para)

summary(wsnecl_lf_mod1)


#creating a binary wasp survival column (0 if no wasps emerged, 1 if at least some did) to see what that analysis
  ##looks like using logistf
tv.para$bin.wsecl <- ifelse(tv.para$num.ecl>0, 1, 0)

#removed load from the model because it cannot inform survival in binary--an is strongly correlated?
  ##all 1s with have loads, almost all 0s will not
wsecl_bin_mod1 <- logistf(bin.wsecl ~ temp.avg*temp.var,
                          data=tv.para)
summary(wsecl_bin_mod1)
