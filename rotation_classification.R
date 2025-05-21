# Code used for Buron et al. (2025) - Ubiquitous and unexpected neonicotinoid contaminations 
# in agricultural soils: investigating the role of cropping history and spatial transport

# Hierarchical clustering

library(sf)
library(dplyr)
library(readxl)

data = read_xlsx("data_Buron.xlsx")
classif_data = read_xlsx("classif_data.xlsx")

# Create sequence list 
data = data %>% 
  rowwise() %>% 
  mutate(list_cult = list(c(crop_2018, crop_2019, crop_2020, crop_2021, crop_2022, crop_2023))) %>% 
  ungroup()

## Levhenstein distance ---------------
# Function to compare sequence non-identities
calculate_similarity <- function(row1, row2) {
  sum(row1 != row2, na.rm = TRUE) 
}

data_lehven = data %>%
  as.data.frame() %>%
  select(crop_2018, crop_2019, crop_2020, crop_2021, crop_2022, crop_2023)

# Dissimilarity matrix
dissimilarity_levenshtein = outer(1:nrow(data_lehven), 1:nrow(data_lehven), 
                                  Vectorize(function(i, j) calculate_similarity(data_lehven[i, ], data_lehven[j, ])))

dissimilarity_levenshtein_matrix = as.dist(dissimilarity_levenshtein)
dissimilarity_levenshtein_matrix = as.matrix(dissimilarity_levenshtein_matrix)

## Number of crop type occurrences ---------------
dissimilarity_occurence = data
for (j in unique(classif_data$Agri_group)) {
  dissimilarity_occurence = dissimilarity_occurence %>% 
    rowwise() %>% 
    mutate({{j}} := sum(list_cult == {{j}}, na.rm = TRUE)) %>% 
    ungroup()
}

# Variable scaling
for (j in unique(classif_data$Agri_group)) {
  dissimilarity_occurence = dissimilarity_occurence %>%
    mutate(!!j := scale(.data[[j]]))
}

# Empty columns removal
dissimilarity_occurence = dissimilarity_occurence %>% 
  select(all_of(unique(cod_cult_sens4$Agri_group))) %>% 
  select(where(~ any(!is.na(.))))

# Dissimilarity matrix 
dissimilarity_occurence_matrix = dist(dissimilarity_occurence, method = "euclidean")
dissimilarity_occurence_matrix = as.matrix(dissimilarity_occurence_matrix)

## Longest common sequence ---------------
data_Lcommon = data %>%
  select(crop_2018, crop_2019, crop_2020, crop_2021, crop_2022, crop_2023)

Lcommon_matrix = matrix(0, nrow = nrow(data_Lcommon)+1, ncol = nrow(data_Lcommon)+1)
for(i in 1:nrow(data_Lcommon)) {
  for(j in (i):nrow(data_Lcommon)){
    Lcommon_matrix[i,j] = qualV::LCS(unlist(data_Lcommon[i,]), unlist(data_Lcommon[j,]))$LLCS
  }
}
Lcommon_matrix[lower.tri(Lcommon_matrix)] = t(Lcommon_matrix)[lower.tri(Lcommon_matrix)]

# Dissimilarity matrix 
Lcommon_matrix = as.dist(Lcommon_matrix)
Lcommon_matrix = as.matrix(Lcommon_matrix)

## Combine the three matrices ---------------
combined_dissimilarity <- matrix(0, nrow = nrow(data), ncol = nrow(data))

for (i in 1:nrow(data)) {
  for (j in 1:nrow(data)) {
    combined_dissimilarity[i, j] = mean(c(dissimilarity_levenshtein_matrix[i, j],
                                          dissimilarity_occurence_matrix[i, j],
                                          Lcommon_matrix[i, j]), na.rm = TRUE)
  }
}

## Clustering ---------------
hc = hclust(as.dist(combined_dissimilarity), method = "ward.D2")

# Plot inertie
inertie = sort(hc$height, decreasing = TRUE) 
plot(inertie[1:30], type ="s")

# Choose group number
x_clust = 12 # Change if needed
clusters = cutree(hc, k = x_clust)
plot(hc, main = "Hierarchical classification dendogram")
rect.hclust(hc, k = x_clust, border = 2:4)

# Add classes to data
data$cluster = as.factor(clusters)

####### Plot classes average composition ---------------
meann = data
for (j in unique(classif_data$Agri_group)) {
  meann = meann %>% 
    rowwise() %>% 
    mutate({{j}} := sum(list_cult == {{j}}, na.rm = TRUE)) %>% 
    ungroup()
}
meann = meann %>% 
  mutate("No data" = map_int(list_cult, ~ sum(is.na(.))))

meann = meann %>%
  group_by(cluster) %>% 
  dplyr::summarise(across(c(unique(classif_data$Agri_group), "No data"), ~mean(.x, na.rm = T))) %>% 
  pivot_longer(cols = c(unique(classif_data$Agri_group), "No data"), names_to = "culture", values_to = "mean") %>% 
  filter(mean != 0)

ggplot(meann)+
  geom_bar(data = meann, aes(x = cluster, y = mean, fill = culture), stat = "identity")+
  theme_classic()+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_text(face = "bold"))+
  ylab("Mean occurence over 6 years")+
  scale_x_discrete(drop = FALSE)