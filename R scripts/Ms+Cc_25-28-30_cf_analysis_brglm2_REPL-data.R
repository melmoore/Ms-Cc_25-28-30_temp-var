#Ms Cc temp var expt--analysis of wasp survival using brglm2 and logistf--REPL DATA


#load libraries
library(brglm2)
library(readr)
library(logistf)


#data with replication experiment

tvor <- read_csv("data files/Ms-Cc_tv-orig-rep_comb_cl.csv", 
                 col_types = cols(temp.avg = col_factor(levels = c("25", "28", "30")), 
                                  temp.var = col_factor(levels = c("0",  "5", "10"))))
View(tvor)



tvor_lng <- read_csv("data files/Ms-Cc_tv-orig-rep_comb_lng.csv", 
                     col_types = cols(temp.avg = col_factor(levels = c("25", "28", "30")), 
                                      temp.var = col_factor(levels = c("0", "5", "10"))))
View(tvor_lng)


#make sorting column to remove orig 30 data
tvor$keep <- ifelse(tvor$temp.avg==30 & tvor$expt=="orig", 0, 1)
tvor_lng$keep <- ifelse(tvor_lng$temp.avg==30 & tvor_lng$expt=="orig", 0, 1)

tvor <- subset(tvor, keep==1)
tvor_lng <- subset(tvor_lng, keep==1)


#----------------------

#necessary data cleaning and configuration for analyses

#log mass
tvor_lng$log_mss <- log(tvor_lng$mass)


#remove individuals with overly large loads (>300)
tvor$keep_ld <- ifelse(tvor$treatment=="para" & tvor$load > 300, 0, 1)
tvor_lng$keep_ld <- ifelse(tvor_lng$end.class=="em" & tvor_lng$load > 300, 0, 1)

tvor <- subset(tvor, keep_ld==1)
tvor_lng <- subset(tvor_lng, keep_ld==1)

#remove wanderers and WOWEs that wandered
tvor$keep_p <- ifelse(tvor$treatment=="para" & tvor$end.class=="wand", 0, 1)
tvor_lng$keep_p <- ifelse(tvor_lng$treatment=="para" & tvor_lng$end.class=="wand", 0, 1)

tvor <- subset(tvor, keep_p==1)
tvor_lng <- subset(tvor_lng, keep_p==1)


#remove WOWEs in 30C treatment
tvor$keep_p2 <- ifelse(tvor$temp.avg==30 & tvor$temp.var==0 & tvor$treatment=="para" & tvor$end.class!="em", 0, 1)
tvor_lng$keep_p2 <- ifelse(tvor_lng$temp.avg==30 & tvor_lng$temp.var==0 & tvor_lng$treatment=="para" & tvor_lng$end.class!="em", 0, 1)

tvor <- subset(tvor, keep_p2==1)
tvor_lng <- subset(tvor_lng, keep_p2==1)


#make temp.var a factor, removing the now non existant +/-5 treatment
tvor$temp.var <- factor(tvor$temp.var, levels = c(0, 10))
tvor_lng$temp.var <- factor(tvor_lng$temp.var, levels = c(0, 10))


#create tot.died column for analysing wasp survival 
tvor$tot.died <- tvor$load - tvor$num.ecl

#make subsets that have only parasitized treatments
tvor_p <- subset(tvor, treatment=="para")
tvor_lngp <- subset(tvor_lng, treatment=="para")



#-------------------

#binomial glm model of wasp survival to eclosion

wsecl_mod1 <- glm(cbind(num.ecl, tot.died) ~ temp.avg * temp.var * load,
                  family = binomial,
                  data=tvor_p,
                  na.action=na.omit)

summary(wsecl_mod1)


#checking for infinite estimates--not running for some reason, has error about not finding starting values
inf_est_check <- check_infinite_estimates(wsecl_mod1)

View(inf_est_check)



#detect separation in in my model parameters
wsecl_sep <- glm(cbind(num.ecl, tot.died) ~ temp.avg*temp.var*load,
                 family = binomial,
                 method = "detect_separation",
                 data=tvor_p,
                 na.action=na.omit)

wsecl_sep




#-----------

#attempting to use logistf model--will not run--left it running for 10 hours and it did not complete

wpsecl_lf_mod1 <- logistf(ps.ecl ~ temp.avg*temp.var*load,
                          data=tvor_p)



#----------------

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

tvorp.bnry <- data.frame(surv=unlist(apply(tvor_p[,c("num.ecl", "load")], 1, bnry.func)),
                         bug.id=rep(tvor_p$bug.id, tvor_p$load), temp.avg=rep(tvor_p$temp.avg, tvor_p$load),
                         temp.var=rep(tvor_p$temp.var, tvor_p$load), load=rep(tvor_p$load, tvor_p$load))
dim(tvorp.bnry)
View(tvorp.bnry)



#now attempting to run a logistf model using my binary response data
#this also seems to just run infinitely....not sure why this isn't working when the original data set did....
wsbnry_lf_mod1 <- logistf(surv ~ temp.avg*temp.var*load,
                          data=tvp.bnry)

summary(wsbnry_lf_mod1)
exp(coef(wsbnry_lf_mod1))
logistftest(wsbnry_lf_mod1)




