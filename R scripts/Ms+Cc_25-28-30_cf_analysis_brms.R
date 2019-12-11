#Ms Cc temp var expt--analysis of wasp survival using brms

#load libraries
library(readr)
library(brms)



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



#----------------

#practice with building brms models, following https://www.rensvandeschoot.com/tutorials/brms-started/

#intercept only model, with num.ecl as response
  ##not enough iterations, small burn in 
wsnecl_int_mod1 <- brm(num.ecl ~ 1 + (1|bug.id),
                       data=tv.para,
                       warmup = 100, 
                       iter   = 200, 
                       chains = 2, 
                       inits  = "random",
                       cores  = 2)
  
  


#increase iterations and burn in 
  #still doesn't run, but am moving on to adding predictors
wsnecl_int_mod2 <- brm(num.ecl ~ 1 + (1|bug.id),
                       data=tv.para,
                       warmup = 2000, iter = 5000, 
                       cores = 2, chains = 2, 
                       seed = 123)  

#calling pairs() on model to diagnose sampling problem (from https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup)
#do not understand how to read this....
pairs(wsnecl_int_mod2)

#attempting to increase the adapt_delta() parameter in the model
#did not seem to help
wsnecl_int_mod2.5 <- brm(num.ecl ~ 1 + (1|bug.id),
                       data=tv.para,
                       warmup = 2000, iter = 5000, 
                       cores = 2, chains = 2, 
                       control = list(adapt_delta = 0.99),
                       seed = 123)  




#add predictors to model
wsnecl_pred_mod2 <- brm(num.ecl ~ 1 + temp.avg +temp.var +load + (1|bug.id),
                       data=tv.para,
                       warmup = 2000, iter = 5000, 
                       cores = 2, chains = 2, 
                       seed = 123) 




#-------------------------------

#moving on to try my actual model structure (binomial glm to see if I can get that to work)











