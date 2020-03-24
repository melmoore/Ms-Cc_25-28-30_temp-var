gam_mass_mod<-gam(log_mss ~ s(age, by = interaction(treatment, temp.avg, temp.var, k=20, bs="ts")) 
                  + s(bug.id, bs ="re") + treatment * temp.avg * temp.var,
                  method="ML", data=tvor_mass, na.action = na.omit)



#Make a null model for model testing, where the smooth of age is not affected by temp.var, temp.avg or treat
gam_mnull_mod<-gam(log_mss ~ s(age, k=10, bs="ts") 
                   + s(bug.id, bs="re") + treatment * temp.avg * temp.var,
                   method="ML", data=tvor_mass, na.action = na.omit)

anova(gam_mnull_mod)



#make reduced gam models to test and see which is best
gam_mass_mod_trta <- gam(log_mss ~ s(age, by = interaction(treatment, temp.avg, k=20, bs="ts")) 
                        + s(bug.id, bs ="re") + treatment * temp.avg * temp.var,
                        method="ML", data=tvor_mass, na.action = na.omit)



gam_mass_mod_trtv<-gam(log_mss ~ s(age, by = interaction(treatment, temp.var, k=20, bs="ts")) 
                  + s(bug.id, bs ="re") + treatment * temp.avg * temp.var,
                  method="ML", data=tvor_mass, na.action = na.omit)


gam_mass_mod_tatv<-gam(log_mss ~ s(age, by = interaction(temp.avg, temp.var, k=20, bs="ts")) 
                  + s(bug.id, bs ="re") + treatment * temp.avg * temp.var,
                  method="ML", data=tvor_mass, na.action = na.omit)



#make reduced models where fixed effects are reduced, but keep smooth interaction
gam_mass_mod_ftr<-gam(log_mss ~ s(age, by = interaction(treatment, temp.avg, temp.var, k=20, bs="ts")) 
                  + s(bug.id, bs ="re") + treatment,
                  method="ML", data=tvor_mass, na.action = na.omit)


gam_mass_mod_fta<-gam(log_mss ~ s(age, by = interaction(treatment, temp.avg, temp.var, k=20, bs="ts")) 
                     + s(bug.id, bs ="re") + temp.avg,
                     method="ML", data=tvor_mass, na.action = na.omit)


gam_mass_mod_ftv<-gam(log_mss ~ s(age, by = interaction(treatment, temp.avg, temp.var, k=20, bs="ts")) 
                      + s(bug.id, bs ="re") + temp.var,
                      method="ML", data=tvor_mass, na.action = na.omit)


gam_mass_mod_ftatv<-gam(log_mss ~ s(age, by = interaction(treatment, temp.avg, temp.var, k=20, bs="ts")) 
                      + s(bug.id, bs ="re") + temp.avg:temp.var,
                      method="ML", data=tvor_mass, na.action = na.omit)


gam_mass_mod_ftatr<-gam(log_mss ~ s(age, by = interaction(treatment, temp.avg, temp.var, k=20, bs="ts")) 
                        + s(bug.id, bs ="re") + temp.avg:treatment,
                        method="ML", data=tvor_mass, na.action = na.omit)


gam_mass_mod_ftrtv<-gam(log_mss ~ s(age, by = interaction(treatment, temp.avg, temp.var, k=20, bs="ts")) 
                        + s(bug.id, bs ="re") + treatment:temp.var,
                        method="ML", data=tvor_mass, na.action = na.omit)


gam_mass_mod_ftatvtr<-gam(log_mss ~ s(age, by = interaction(treatment, temp.avg, temp.var, k=20, bs="ts")) 
                        + s(bug.id, bs ="re") + temp.avg:temp.var:treatment,
                        method="ML", data=tvor_mass, na.action = na.omit)


#test the reduced models

anova(gam_mnull_mod, gam_mass_mod, gam_mass_mod_tatv, gam_mass_mod_trta, gam_mass_mod_trtv,
      gam_mass_mod_fta, gam_mass_mod_ftr, gam_mass_mod_ftv, gam_mass_mod_ftatr, gam_mass_mod_ftatv, 
      gam_mass_mod_trtv, gam_mass_mod_ftatvtr, test="Chisq")

AIC(gam_mnull_mod, gam_mass_mod, gam_mass_mod_tatv, gam_mass_mod_trta, gam_mass_mod_trtv)
