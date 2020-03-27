#Ms Cc temp var expt--analysis of wasp survival using brglm2 and logistf


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

#initial data wrangling--ORIGINAL 30 DATA ONLY

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


#remove empty +/-5 factor level, I think it was slowing down the model a ton
tv.para$temp.var <- factor(tv.para$temp.var, levels=c("0", "10"))

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
  ##the only thing I can think of is ps.ecl is a proportion (and therefor between 0 and 1),
  ##where as num.ecl is an integer?

wpsecl_lf_mod1 <- logistf(tot.surv ~ temp.avg*temp.var*load,
                         data=tv.para)

summary(wpsecl_lf_mod1)
logistftest(wpsecl_lf_mod1)


#removing load 
wpsecl_lf_mod2 <- logistf(tot.surv ~ temp.avg*temp.var,
                          data=tv.para)

summary(wpsecl_lf_mod2)


#trying cbind of num.ecl and tot died as response variable
wsmatecl_lf_mod1 <- logistf(cbind(num.ecl, tot.died) ~ temp.avg*temp.var*load,
                            data=tv.para)


#not using firth methods--just returned NAN for everything

wpsecl_lf_mod1.5 <- logistf(tot.surv ~ temp.avg*temp.var*load,
                          data=tv.para,
                          firth=FALSE, pl=FALSE)

summary(wpsecl_lf_mod1.5)

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



#Is there a way to rearrange my data so that each wasp has it's own row, with a binary survival variable?
  ##that would be a huge data frame, but would maybe allow me to answer my research question with a 
  ##logisitic regression model?

#using function from https://sakai.unc.edu/access/content/group/3d1eb92e-7848-4f55-90c3-7c72a54e7e43/public/docs/lectures/lecture26.htm
  ##should generate a number of 1s = num.ecl, and a number of 0s = load - num.ecl (or tot.died) 
bnry.func <- function(x) rep(c(1,0),c(x[1],x[2]-x[1]))

# sample calculation
test <- apply(tv.para[1:4, c("num.ecl", "load")], 1, bnry.func)

#now apply to whole data frame: unlist the results, and add columns indicating the values of the 
  ##bug.id, temp.var, temp.avg, load (do I need this?) repeated the necessary number of times.
  ##wowes not included because load is 0--figure out how to fix

tvp.bnry <- data.frame(surv=unlist(apply(tv.para[,c("num.ecl", "load")], 1, bnry.func)),
                       bug.id=rep(tv.para$bug.id, tv.para$load), temp.avg=rep(tv.para$temp.avg, tv.para$load),
                       temp.var=rep(tv.para$temp.var, tv.para$load), load=rep(tv.para$load, tv.para$load))
dim(tvp.bnry)
View(tvp.bnry)


#now attempting to run a logistf model using my binary response data
  ##chisq are inf and p value is 0 for temp.avg 30....
wsbnry_lf_mod1 <- logistf(surv ~ temp.avg*temp.var*load,
                          data=tvp.bnry)

summary(wsbnry_lf_mod1)
exp(coef(wsbnry_lf_mod1))
logistftest(wsbnry_lf_mod1)


#trying model without penalized likelihood
wsbnry_lf_mod1.5 <- logistf(surv ~ temp.avg*temp.var*load,
                          data=tvp.bnry,
                          pl=FALSE)

summary(wsbnry_lf_mod1.5)
exp(coef(wsbnry_lf_mod1.5))
logistftest(wsbnry_lf_mod1.5)



#trying without load, since I don't quite know if it makes sense to include
##chisq are inf and p value is 0 for temp.avg 30 and 28:10......
wsbnry_lf_mod2 <- logistf(surv ~ temp.avg+temp.var,
                          data=tvp.bnry)

summary(wsbnry_lf_mod2)
logistftest(wsbnry_lf_mod2)


tvp.bnry.30 <- subset(tvp.bnry, temp.avg==30)



