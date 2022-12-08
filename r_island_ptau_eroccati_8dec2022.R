#THIS IS THE SCRIPT USED IN THE ANALYSIS OF ISLAND P-TAU DATA
#ISLAND STUDY LINKING AGEING AND NEURODEGENERATIVE DISEASE
#WRITTEN BY DR EDDY ROCCATI, POST-DOC AT WICKING DEMENTIA CENTRE
#PLEASE CONTACT eddy.roccati@utas.edu.au for any queries

#### install packages ========================================================
library(tidyverse)
library(dplyr)
library(readr)
library(plyr)
library(easycsv)
library(table1)
library(tableone)
library(performance)
library(gam)
library(sjPlot)
library(plotly)
library(mgcv)
library(tidymv)
library(stargazer)
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)
library(ggdist)

#### extract ISLAND data =====================================================
island_df <- list(apoe, bkgd, cantab, drp, hads, 
                  participants, ptau, social)%>% reduce(full_join, by = "UID")

#### clean ISLAND data =======================================================
#let's first take out duplicate cases
island_df <- island_df %>% distinct(UID, .keep_all = TRUE)

#total school years
island_df$university_years[is.na(island_df$university_years)] <- 0
island_df$total_school_years = island_df$school_years + island_df$university_years

#remove outlier p-tau conc of over 13 pg/mL
island_df <- island_df %>% filter(!ptau181_conc > 13)

#change zeros so can fit log without negative infinity
island_df$PALTEA28_1 <- island_df$PALTEA28 + 1

#### assign labels for tables ================================================
label(island_df$age_at_bkgd) <- "Age in years"
label(island_df$marital_status) <- "What is your current marital status?"
label(island_df$university_qualification) <- "Have you completed university studies?"
label(island_df$university_years) <- "How many years of university education?"
label(island_df$PALTEA28) <- "CANTAB PAL Total Errors (Adjusted)"
label(island_df$PALMETS28) <- "CANTAB PAL Mean Errors to Success"
label(island_df$PALNPR28) <- "CANTAB PAL Number of Patterns Reached"
label(island_df$SWMPR) <- "CANTAB SWM Problem Reached"
label(island_df$SWMBE468) <- "CANTAB SWM Between Errors"
label(island_df$SWMSX) <- "CANTAB SWM Strategy"
label(island_df$school_years) <- "Years of education at the school level eg. up to year 12"
label(island_df$highest_education_level) <- "What is the highest level of education you have obtained?"
label(island_df$retired) <- "Have you retired?"
label(island_df$employed) <- "Are you currently employed?"
label(island_df$volunteer ) <- "Do you work as a volunteer?"
label(island_df$language_english_only) <- "Language - English only"
label(island_df$dementia_diagnosis) <- "Have you been told by a doctor that you have dementia?"
label(island_df$b12_deficiency_diagnosis) <- "Have you been diagnosed with a vitamin B12 deficiency?"
label(island_df$cns_diagnosis) <- "Have you been diagnosed with a central nervous system degenerative disease eg. Parkinson's, Huntington's, Multiple Sclerosis?"
label(island_df$cancer_diagnosis) <- "Have you been diagnosed with cancer?"
label(island_df$delerium_diagnosis) <- "Have you been diagnosed with delirium?"
label(island_df$dementia_family_history) <- "Is there a history of conditions such as dementia in your direct family  for example  siblings, parents, grandparents, aunties and uncles ?"
label(island_df$epilepsy_diagnosis) <- "Have you been diagnosed with epilepsy?"
label(island_df$head_injury) <- "Have you ever had a serious head injury?"
label(island_df$hearing_impairment) <- "Do you have a hearing impairment?"
label(island_df$heart_disease_diagnosis) <- "Have you been diagnosed with heart disease?"
label(island_df$kidney_disease_diagnosis) <- "Have you been diagnosed with kidney disease?"
label(island_df$liver_disease_diagnosis) <- "Have you been diagnosed with liver disease?"
label(island_df$memory_change) <- "Have you noticed a substantial change in your memory and mental function in recent years?"
label(island_df$memory_impairment_diagnosis) <- "Have you been told by your doctor that you have a memory impairment but they were uncertain if you have dementia?"
label(island_df$pysch_diagnosis) <- "Have you been diagnosed with a psychiatric disorder   eg. depression, psychosis, bipolar disorder, anxiety disorder?"
label(island_df$stroke_tia_attack) <- "Have you had a stroke or experienced transient ischaemic attacks?"
label(island_df$visual_legally_blind) <- "Are you legally blind?"
label(island_df$visual_corrective_glasses) <- "Do you wear corrective glasses?"
label(island_df$medications) <- "Are you currently taking any prescription medications or hormonal supplements?"
label(island_df$visual_colour_blind) <- "Are you colour blind?"
label(island_df$ptau181_conc) <- "Fitted concentration of ptau181 (pg/ml)"
label(island_df$final_genotype) <- "APOE"
label(island_df$total_school_years) <- "Years of education"
label(island_df$total_lubben) <- "Lubben Total Score"
label(island_df$liv_diabetes) <- "Diagnosed with diabetes"
label(island_df$liv_alcohol) <- "Alcohol more than 21 units per week"
label(island_df$liv_smok) <- "Current smoker"
label(island_df$liv_bp) <- "Diagnosed with diabetes"


#### descriptives===============================================================

#table1: demographic table
vars <- c("age_at_bkgd", "total_school_years", "highest_education_level", "gender","marital_status","employed","retired" ,"final_genotype", "PALTEA28","SWMSX", "ptau181_conc", "dementia_family_history")
cat_vars <- c("highest_education_level","marital_status", "employed", "retired"  ,"final_genotype","dementia_family_history" )
table_1 <- CreateTableOne(vars = vars, strata = "gender", data = island_df, factorVars = cat_vars, addOverall = TRUE)
table_1_export <- print(table_1, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(table_1_export, file = "table_1_export.csv")

#table2: modifiable risk factors
vars <- c("Alcohol_score " , " BloodPressure_score " , " BMI_score " , " Cholesterol_score " , " CognitiveActivity_score " , " Diabetes_score " , " Diet_score " , " PhysicalActivity_score" , " Smoking_score")
cat_vars <- c("liv_alcohol", "hearing_impairment", "head_injury", "liv_bp", "smoking", "liv_smok")
table_2 <- CreateTableOne(vars = vars, strata = "gender", data = island_df, factorVars = cat_vars, addOverall = TRUE)
table_2_export <- print(table_2, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(table_2_export, file = "table_2_export.csv")

#table2.1: drp
vars <- c("Alcohol_risk","BloodPressure_risk","BMI_risk","Cholesterol_risk","CognitiveActivity_risk","Diabetes_risk","Diet_risk","PhysicalActivity_risk","Smoking_risk")
cat_vars <- c("Alcohol_risk","BloodPressure_risk","BMI_risk","Cholesterol_risk","CognitiveActivity_risk","Diabetes_risk","Diet_risk","PhysicalActivity_risk","Smoking_risk")
table_2.1 <- CreateTableOne(vars = vars, strata = "gender", data = island_df, factorVars = cat_vars, addOverall = TRUE)
table_2.1_export <- print(table_2.1, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(table_2.1_export, file = "table_2.1_export.csv")


#### regression=================================================================
#dependent: ptau181_conc
#modifiable: total_school_years  +  hearing_impairment  +  head_injury  +  liv_bp + Alcohol_score + BMI_score + liv_smok + depression_total + total_lubben + PhysicalActivity_score + liv_diabetes
#unmodifiable: age_at_bkgd + gender + final_genotype
#cognition: PALTEA28 + SWMSX
#data=island_df
#DRP score: Alcohol_score + BloodPressure_score + BMI_score + Cholesterol_score + CognitiveActivity_score + Diabetes_score + Diet_score + PhysicalActivity_score+ Smoking_score#DRP score: Alcohol_score + BloodPressure_score + BMI_score + Cholesterol_score + CognitiveActivity_score + Diabetes_score + Diet_score + PhysicalActivity_score+ Smoking_score
#DRP risk: Alcohol_risk + BloodPressure_risk + BMI_risk + Cholesterol_risk + CognitiveActivity_risk + Diabetes_risk + Diet_risk + PhysicalActivity_risk+ Smoking_risk

###### COGNITION
#run linear models on PAL
lcm1 <- lm(log(PALTEA28_1) ~ ptau181_conc, data = island_df)
summary(lcm1)
lcm2 <- lm(log(PALTEA28_1) ~ ptau181_conc + age_at_bkgd + gender + total_school_years + e4, data = island_df)
summary(lcm2)
#check_model(lcm2)
#AIC(lcm1, lcm2)
stargazer(lcm1, lcm2, title="Results Table: PAL", align=TRUE, type= "html", out = "lcm.htm")

#run linear models on SWM
lcn1 <- lm(SWMSX ~ ptau181_conc, data = island_df)
summary(lcn1)
lcn2 <- lm(SWMSX ~ ptau181_conc + age_at_bkgd + gender + total_school_years + e4, data = island_df)
summary(lcn2)
#check_model(lcn2)
#AIC(lcn1, lcn2)
stargazer(lcn1, lcn2, title="Results Table: SWM", align=TRUE, type= "html", out = "lcn.htm")


###### P-TAU
gp1 <- lm(log(ptau181_conc) ~ Alcohol_score + BloodPressure_score + BMI_score + Cholesterol_score + CognitiveActivity_score + Diabetes_score + Diet_score + PhysicalActivity_score + Smoking_score, data = island_df)
summary(gp1)
#check_model(gp1)
gp2 <- lm(log(ptau181_conc) ~ Alcohol_score + BloodPressure_score + BMI_score + Cholesterol_score + CognitiveActivity_score + Diabetes_score + Diet_score + PhysicalActivity_score + Smoking_score + age_at_bkgd + gender + total_school_years + e4, data = island_df)
summary(gp2)
#check_model(gp2)
#AIC(gp1, gp2)
stargazer(gp1, gp2, title="Results Table: P-tau", align=TRUE, type= "html", out = "gcp.htm")


#### figures ==================================================================
###### DRP & P-TAU
p1 <- island_df %>% filter(Alcohol_risk=="high" | Alcohol_risk=="medium" | Alcohol_risk=="low") %>% 
  ggplot(aes(x= Alcohol_risk, y=ptau181_conc, group=Alcohol_risk, colour= Alcohol_risk))+
  ggdist::stat_halfeye(adjust = .5, width = .7, .width = 0, justification = -.2, point_colour = NA) + 
  geom_boxplot(width = .2, outlier.shape = NA) + 
  geom_jitter(width = .05, alpha = .2) +scale_x_discrete(limits = c("low", "medium", "high"))+
  theme(legend.position="none", axis.title.y=element_blank())+labs(x = "Alcohol")
p2 <- island_df %>% filter(BloodPressure_risk=="high" | BloodPressure_risk=="medium" | BloodPressure_risk=="low") %>% 
  ggplot(aes(x= BloodPressure_risk, y=ptau181_conc, group=BloodPressure_risk, colour= BloodPressure_risk))+
  ggdist::stat_halfeye(adjust = .5, width = .7, .width = 0, justification = -.2, point_colour = NA) + 
  geom_boxplot(width = .2, outlier.shape = NA) + 
  geom_jitter(width = .05, alpha = .2) +scale_x_discrete(limits = c("low", "medium", "high"))+
  theme(legend.position="none", axis.title.y=element_blank())+labs(x = "Blood Pressure")
p3 <- island_df %>% filter(BMI_risk=="high" | BMI_risk=="medium" | BMI_risk=="low") %>% 
  ggplot(aes(x= BMI_risk, y=ptau181_conc, group=BMI_risk, colour= BMI_risk))+
  ggdist::stat_halfeye(adjust = .5, width = .7, .width = 0, justification = -.2, point_colour = NA) + 
  geom_boxplot(width = .2, outlier.shape = NA) + 
  geom_jitter(width = .05, alpha = .2) +scale_x_discrete(limits = c("low", "medium", "high"))+
  theme(legend.position="none", axis.title.y=element_blank())+labs(x = "BMI")
p4 <- island_df %>% filter(Cholesterol_risk=="high" | Cholesterol_risk=="medium" | Cholesterol_risk=="low") %>% 
  ggplot(aes(x= Cholesterol_risk, y=ptau181_conc, group=Cholesterol_risk, colour= Cholesterol_risk))+
  ggdist::stat_halfeye(adjust = .5, width = .7, .width = 0, justification = -.2, point_colour = NA) + 
  geom_boxplot(width = .2, outlier.shape = NA) + 
  geom_jitter(width = .05, alpha = .2) +scale_x_discrete(limits = c("low", "medium", "high"))+
  theme(legend.position="none", axis.title.y=element_blank())+labs(x = "Cholesterol")
p5 <- island_df %>% filter(CognitiveActivity_risk=="high" | CognitiveActivity_risk=="medium" | CognitiveActivity_risk=="low") %>% 
  ggplot(aes(x= CognitiveActivity_risk, y=ptau181_conc, group=CognitiveActivity_risk, colour= CognitiveActivity_risk))+
  ggdist::stat_halfeye(adjust = .5, width = .7, .width = 0, justification = -.2, point_colour = NA) + 
  geom_boxplot(width = .2, outlier.shape = NA) + 
  geom_jitter(width = .05, alpha = .2) +scale_x_discrete(limits = c("low", "high"))+
  theme(legend.position="none", axis.title.y=element_blank())+labs(x = "Cognitive Activity")
p6 <- island_df %>% filter(Diabetes_risk=="high" | Diabetes_risk=="medium" | Diabetes_risk=="low") %>% 
  ggplot(aes(x= Diabetes_risk, y=ptau181_conc, group=Diabetes_risk, colour= Diabetes_risk))+
  ggdist::stat_halfeye(adjust = .5, width = .7, .width = 0, justification = -.2, point_colour = NA) + 
  geom_boxplot(width = .2, outlier.shape = NA) + 
  geom_jitter(width = .05, alpha = .2) +scale_x_discrete(limits = c("low", "medium", "high"))+
  theme(legend.position="none", axis.title.y=element_blank())+labs(x = "Diabetes")
p7 <- island_df %>% filter(Diet_risk=="high" | Diet_risk=="medium" | Diet_risk=="low") %>% 
  ggplot(aes(x= Diet_risk, y=ptau181_conc, group=Diet_risk, colour= Diet_risk)) +
  ggdist::stat_halfeye(adjust = .5, width = .7, .width = 0, justification = -.2, point_colour = NA) + 
  geom_boxplot(width = .2, outlier.shape = NA) + 
  geom_jitter(width = .05, alpha = .2) +scale_x_discrete(limits = c("low","medium", "high"))+
  theme(legend.position="none", axis.title.y=element_blank())+labs(x = "Diet")
p8 <- island_df %>% filter(PhysicalActivity_risk=="high" | PhysicalActivity_risk=="medium" | PhysicalActivity_risk=="low") %>% 
  ggplot(aes(x= PhysicalActivity_risk, y=ptau181_conc, group=PhysicalActivity_risk, colour= PhysicalActivity_risk))+
  ggdist::stat_halfeye(adjust = .5, width = .7, .width = 0, justification = -.2, point_colour = NA) + 
  geom_boxplot(width = .2, outlier.shape = NA) + 
  geom_jitter(width = .05, alpha = .2) +scale_x_discrete(limits = c("low", "high"))+
  theme(legend.position="none", axis.title.y=element_blank())+labs(x = "Physical Activity")
p9 <- island_df %>% filter(Smoking_risk=="high" | Smoking_risk=="medium" | Smoking_risk=="low") %>% 
  ggplot(aes(x= Smoking_risk, y=ptau181_conc, group=Smoking_risk, colour= Smoking_risk))+
  ggdist::stat_halfeye(adjust = .5, width = .7, .width = 0, justification = -.2, point_colour = NA) + 
  geom_boxplot(width = .2, outlier.shape = NA) + 
  geom_jitter(width = .05, alpha = .2) +scale_x_discrete(limits = c("low", "medium", "high"))+
  theme(legend.position="none", axis.title.y=element_blank())+labs(x = "Smoking")


p_all <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9,
                   ncol = 3, nrow = 3)
require(grid)
annotate_figure(p_all, left = textGrob("Phosphorylatd plasma p-tau 181 (pg/mL)", rot = 90, vjust = 0.5, gp = gpar(cex = 1.3)))

p_all

###### FOREST PLOT LM
fit1 <- lm(ptau181_conc ~ Alcohol_score + BloodPressure_score + BMI_score + Cholesterol_score + CognitiveActivity_score + Diabetes_score + Diet_score + PhysicalActivity_score + Smoking_score , data = island_df)
fit2 <- update(fit1, . ~ . + age_at_bkgd + gender + total_school_years + e4)

forest <- plot_models(fit1, fit2, std.est = "std2", axis.labels = c("APOE: e4", "Education", "Gender: Male", "Age", "Smoking", "Physical Activity", "Diet", "Diabetes", "Cognitive Activity" , "Cholesterol", "BMI","Blood Pressure", "Alcohol"),
                      m.labels = c("Unadjusted", "Adjusted"), colors = c("darkgreen", "darkorange"),
                      legend.title = "Models") + theme(legend.position="bottom") + ylim(-.5, .5)



