#### dCT summary stats ####
ref = "name_reference_gene"

dat2 <- dat1 %>% 
  group_by(Rep, Sample.Name, Target.Name) %>% 
  summarise(delta.CT = mean(CT)) %>% 
  mutate_at(vars(matches("delta")), funs(.- .[Target.Name == ref])) %>%  # calc delta.CT
  filter(Target.Name != ref) %>% 
  group_by(Sample.Name, Target.Name) %>% 
  summarise_each(funs(mean(., na.rm = T), 
                      n = sum(!is.na(.)), 
                      se = sd(., na.rm = T)/sqrt(sum(!is.na(.)))), 
                 delta.CT) %>% 
  mutate(TrA = 2^-mean, # transcript abundance (TrA)
         perc.TrA = TrA * 100, # percentage TrA
         lw.TrA = 2^ - (mean + se),  # lower error position
         up.TrA = 2^ - (mean - se)) %>%  # upper error position
  rename(dCT = mean) %>% 
  mutate_if(is.numeric, round, 5)


#### ddCT summary stats ####
dat3 <- dat1 %>% 
  group_by(Target.Name) %>%
  mutate_at(vars(matches("dCT")), funs(.- .[Sample.Name == "WT"])) %>% 
  select(1:5) %>% 
  rename(ddCT = dCT) %>% 
  mutate(fc = 2^-ddCT, # fold change ddCT
         fc.se = se^2) %>% # SE fold change
  mutate_at(vars(matches("fc.se")), funs(sqrt(.+. [Sample.Name == "WT"]))) %>% 
  mutate(lw.fc = 2^ - (ddCT + fc.se),  # lower error position
         up.fc = 2^ - (ddCT - fc.se)) %>% # upper error position
  mutate_if(is.numeric, round, 5)




