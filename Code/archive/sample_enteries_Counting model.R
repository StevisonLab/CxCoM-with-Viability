library(readr)
library(tidyverse)
library(magrittr)
library(dplyr)
library(emmeans)
library(ggplot2)
library(tidyr)
library(dplyr)
library(tidyr)
library(readr)

# Load merged data
merged <- read_csv("merged_data_output.csv", col_types = cols()) %>%
  mutate(
    VialID_Brood = paste(VialID, Brood, sep = "_"),
    F1.developmental = as.character(`F1 developmental`),
    F1.Genotype = as.character(`F1 Genotype`)
  )

# Load survival data
survival <- read_csv("mean_survivalmarkerfree.csv", col_types = cols()) %>%
  mutate(
    F1.developmental = as.character(F1developmental),
    F1.Genotype = as.character(F1Genotype)
  )

# Merge survival info
merged <- merged %>%
  left_join(survival, by = c("F1.Genotype", "F1.developmental"))

all_classes <- as.character(1:132)

# Sum counts per vial+class
entries <- merged %>%
  group_by(VialID_Brood, Class) %>%
  summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    id_cols = VialID_Brood,
    names_from = Class,
    values_from = Count,
    values_fill = 0
  )

# Ensure all class columns exist
missing <- setdiff(all_classes, names(entries))
if(length(missing)) entries[missing] <- 0

# Attach metadata and survival columns (just first value per vial)
entries_final <- entries %>%
  left_join(
    merged %>%
      group_by(VialID_Brood) %>%
      summarise(
        total_egg_count  = ceiling(first(Eggsumcount)),   # take first and round up
        Survival_Total   = first(Survival_Total),
        Mean_Female_Survival = first(Mean_Female_Survival),
        Mean_Male_Survival   = first(Mean_Male_Survival),
        F1.developmental = first(F1.developmental),
        F1.Genotype      = first(F1.Genotype),
        .groups = "drop"
      ),
    by = "VialID_Brood"
  ) %>%
  mutate(
    `129` = total_egg_count,
    `130` = if_else(!is.na(Survival_Total), pmin(Survival_Total / 100, 1), 0),
    `131` = if_else(!is.na(Mean_Female_Survival), pmin(Mean_Female_Survival / 100, 1), 0),
    `132` = if_else(!is.na(Mean_Male_Survival), pmin(Mean_Male_Survival / 100, 1), 0),
    treatment = paste(F1.developmental, F1.Genotype, sep = "_"),
    vial = VialID_Brood
  ) %>%
  #dplyr::select(VialID_Brood, dplyr::all_of(all_classes), `129`, `130`, `131`, `132`, treatment, vial)
  dplyr::select(VialID_Brood, dplyr::all_of(all_classes), `129`, `130`, `131`, `132`)
# save
write_csv(entries_final, "enteries_countingmodel.csv")
