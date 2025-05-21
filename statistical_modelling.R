# Code used for Buron et al. (2025) - Ubiquitous and unexpected neonicotinoid contaminations 
# in agricultural soils: investigating the role of cropping history and spatial transport
# See attached text file for variable names detail

# Statistical modelling

library(readxl)
library(lmerTest)
library(dplyr)
library(glmmTMB)

# Load data
data = read_xlsx("data_Buron.xlsx")

# Statistical modelling
# Note: Bonferronni corrections are not applied in the output of the models shown here 

# Time since last treatment / Number of treatments --> Bonferronni correction of *4

## Time since last treatment 
# Clothianidin
mod1 = data %>% filter(!beet_year %in% c("No\noccurrence")) %>% mutate(beet_year = as.numeric(beet_year))
mod1 = lmer(Clothianidin ~ log(beet_year+1) + (1|Site), data = mod1)

summary(mod1)

# Imidacloprid
mod2 = data %>% filter(!beetveg_year %in% c("No\noccurrence")) %>% mutate(beetveg_year = as.numeric(beetveg_year))
mod2 = glmmTMB::glmmTMB(Imidacloprid ~ log(beetveg_year+1) + (1|Site), data = mod2)

summary(mod2)

## Number of treatment 
data_graph_occ = data %>% 
  mutate(beet_occurrence = case_when(beet_occurrence == 0 ~ "No\noccurrence", TRUE ~ as.character(beet_occurrence)),
         beetveg_occurrence = case_when(beetveg_occurrence == 0 ~ "No\noccurrence", TRUE ~ as.character(beetveg_occurrence)))

mod3 = lmer(Clothianidin ~ beet_occurrence + (1|Site), data = data_graph_occ %>% 
                         filter(beet_occurrence != "No\noccurrence") %>% 
                         mutate(beet_occurrence = as.numeric(beet_occurrence))) 

summary(mod3)

mod4 = lmer(Imidacloprid ~ beetveg_occurrence + (1|Site), data = data_graph_occ %>% 
                         filter(beetveg_occurrence != "No\noccurrence") %>% 
                         mutate(beetveg_occurrence = as.numeric(beetveg_occurrence)))

summary(mod4)

## Spatial transport - Effect of runoff and dust dispersion
# Bonferronni correction of *6

#Clothianidin
#1 years period
m11 = glmmTMB::glmmTMB(any_clo ~ scale(dust_clo1) + scale(runoff_clo1) + (1 | Site), data = data %>% 
                         mutate(any_clo = if_else(Clothianidin > 0, 1, 0)) %>% 
                         filter(beet_occurrence == 0), family = "binomial")
DHARMa::testDispersion(m11)
summary(m11)

#2 years period
m12 = glmmTMB::glmmTMB(any_clo ~ scale(dust_clo2) + scale(runoff_clo2) + (1 | Site), data = data %>% 
                         mutate(any_clo = if_else(Clothianidin > 0, 1, 0)) %>% 
                         filter(beet_occurrence == 0), family = "binomial")
DHARMa::testDispersion(m12)
summary(m12)

#3 years period
m13 = glmmTMB::glmmTMB(any_clo ~ scale(dust_clo3) + scale(runoff_clo3) + (1 | Site), data = data %>% 
                         mutate(any_clo = if_else(Clothianidin > 0, 1, 0)) %>% 
                         filter(beet_occurrence == 0), family = "binomial")
DHARMa::testDispersion(m13)
summary(m13)

#Imidacloprid 
#1 year period
m21 = glmmTMB::glmmTMB(any_imid ~ scale(dust_imid1) + scale(runoff_imid1) + (1 | Site), data = data %>% 
                         mutate(any_imid = if_else(Imidacloprid > 0, 1, 0)) %>% 
                         filter(beetveg_occurrence == 0), family = "binomial") 
DHARMa::testDispersion(m21)
summary(m21)

#2 years period
m22 = glmmTMB::glmmTMB(any_imid ~ scale(dust_imid2) + scale(runoff_imid2) + (1 | Site), data = data %>% 
                         mutate(any_imid = if_else(Imidacloprid > 0, 1, 0)) %>% 
                         filter(beetveg_occurrence == 0), family = "binomial") 
DHARMa::testDispersion(m22)
summary(m22)

#3 years period
m23 = glmmTMB::glmmTMB(any_imid ~ scale(dust_imid3) + scale(runoff_imid3) + (1 | Site), data = data %>% 
                         mutate(any_imid = if_else(Imidacloprid > 0, 1, 0)) %>% 
                         filter(beetveg_occurrence == 0), family = "binomial") 
DHARMa::testDispersion(m23)
summary(m23)