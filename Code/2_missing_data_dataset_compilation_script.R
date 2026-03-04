#setwd("Raw Data/") # Please use Raw Data as working directory

# Read in csv files of raw count data from the marker free cross
EH_marker_free_egg_count <- read.csv("A_EH_raw_marker_free_egg_counts.csv", header=TRUE)
MS_marker_free_egg_count <- read.csv("B_MS_raw_marker_free_egg_counts.csv", header=TRUE)
SK_marker_free_adult_count <- read.csv("C_SK_raw_marker_free_adult_counts.csv", header=TRUE)

# Compile the raw data from different observers and life history stages
EH_raw_marker_free_total_egg_count <- EH_marker_free_egg_count$total_egg_count
MS_raw_marker_free_total_egg_count <- MS_marker_free_egg_count$total_egg_count
compiled_raw_marker_free_dataset <- cbind(SK_marker_free_adult_count,EH_raw_marker_free_total_egg_count,MS_raw_marker_free_total_egg_count)
write.csv(compiled_raw_marker_free_dataset,"../Summaries and Outputs/D_compiled_raw_marker_free_dataset.csv",row.names = TRUE)

# Transform individual raw egg counts to address undercounting the average and round these two counts to the nearest whole number
EH_transformed_marker_free_total_egg_count <-(EH_raw_marker_free_total_egg_count-2.383755)/0.9459906 # Slope and Intercept from Type II Regression on Schneider et al. 2025 dataset
MS_transformed_marker_free_total_egg_count <-(MS_raw_marker_free_total_egg_count-11.07247)/0.7537183 # Slope and Intercept from Type II Regression on Schneider et al. 2025 dataset
transformed_mean_marker_free_total_egg_count <-round((EH_transformed_marker_free_total_egg_count+MS_transformed_marker_free_total_egg_count)/2)
derived_marker_free_dataset <-cbind(compiled_raw_marker_free_dataset,transformed_mean_marker_free_total_egg_count)
write.csv(derived_marker_free_dataset,"../Summaries and Outputs/E_derived_marker_free_dataset.csv",row.names = TRUE)

# Calculate the per vial viability of females, males, and total under the 1:1 primary sex ratio assumption
derived_marker_free_female_viability <- derived_marker_free_dataset$female_adult_count/(derived_marker_free_dataset$transformed_mean_marker_free_total_egg_count/2)
derived_marker_free_male_viability <- derived_marker_free_dataset$male_adult_count/(derived_marker_free_dataset$transformed_mean_marker_free_total_egg_count/2)
derived_marker_free_total_viability <- derived_marker_free_dataset$total_adult_count/derived_marker_free_dataset$transformed_mean_marker_free_total_egg_count
marker_free_viability_dataset <-cbind(derived_marker_free_dataset,derived_marker_free_female_viability,derived_marker_free_male_viability,derived_marker_free_total_viability)
write.csv(marker_free_viability_dataset,"../Summaries and Outputs/F_marker_free_viability_dataset.csv",row.names = TRUE)

# Compute the overall viabilities for reporting in main text
overall_marker_free_viability <- mean(derived_marker_free_total_viability)
female_marker_free_viability <- mean(derived_marker_free_female_viability)
male_marker_free_viability <- mean(derived_marker_free_male_viability)

# Generate a data frame for sex-specific viability effects in the marker free crosses
female_vector <- rep("Female", 100)
male_vector <- rep("Male", 100)
marker_free_female_viabilities <- data.frame(female_vector,marker_free_viability_dataset$cross_id,marker_free_viability_dataset$brood_id,marker_free_viability_dataset$derived_marker_free_female_viability)
names(marker_free_female_viabilities) <- c("Sex","Cross","Brood","Viability")
marker_free_male_viabilities <- data.frame(male_vector,marker_free_viability_dataset$cross_id,marker_free_viability_dataset$brood_id,marker_free_viability_dataset$derived_marker_free_male_viability)
names(marker_free_male_viabilities) <- c("Sex","Cross","Brood","Viability")
marker_free_viability_anova_input <- rbind(marker_free_female_viabilities, marker_free_male_viabilities)

# Convert varible "Cross" to a categorical factor and "Brood" to a numeric value
marker_free_viability_anova_input$Cross <- factor(marker_free_viability_anova_input$Cross)
marker_free_viability_anova_input$Brood <- factor(marker_free_viability_anova_input$Brood)
marker_free_viability_anova_input$Brood <- as.numeric(marker_free_viability_anova_input$Brood)

# Read in csv files of raw count data from the multiply marked cross
EH_multiply_marked_egg_count <- read.csv("1_EH_raw_multiply_marked_egg_counts.csv", header=TRUE)
MS_multiply_marked_egg_count <- read.csv("2_MS_raw_multiply_marked_egg_counts.csv", header=TRUE)
SK_multiply_marked_adult_phenotypic_class_counts <- read.csv("3_SK_raw_multiply_marked_adult_phenotypic_class_counts.csv", header=TRUE)

# Compile the raw data from different observers and life history stages
EH_raw_multiply_marked_total_egg_count <- EH_multiply_marked_egg_count$total_egg_count
MS_raw_multiply_marked_total_egg_count <- MS_multiply_marked_egg_count$total_egg_count
compiled_raw_multiply_marked_dataset <- cbind(SK_multiply_marked_adult_phenotypic_class_counts,EH_raw_multiply_marked_total_egg_count,MS_raw_multiply_marked_total_egg_count)
write.csv(compiled_raw_multiply_marked_dataset,"../Summaries and Outputs/4_compiled_raw_multiply_marked_dataset.csv",row.names = TRUE)

# Transform individual raw egg counts to address undercounting the average and round these two counts to the nearest whole number
EH_transformed_multiply_marked_total_egg_count <-(EH_raw_multiply_marked_total_egg_count-2.383755)/0.9459906 # Slope and Intercept from Type II anova on Schneider et al. 2025 dataset
MS_transformed_multiply_marked_total_egg_count <-(MS_raw_multiply_marked_total_egg_count-11.07247)/0.7537183 # Slope and Intercept from Type II anova on Schneider et al. 2025 dataset
transformed_mean_multiply_marked_total_egg_count <-round((EH_transformed_multiply_marked_total_egg_count+MS_transformed_multiply_marked_total_egg_count)/2)
derived_multiply_marked_dataset <-cbind(compiled_raw_multiply_marked_dataset,transformed_mean_multiply_marked_total_egg_count)
write.csv(derived_multiply_marked_dataset,"../Summaries and Outputs/5_derived_multiply_marked_dataset.csv",row.names = TRUE)

# Calculate the per vial viability of females, males, and total under the 1:1 primary sex ratio assumption
derived_multiply_marked_female_viability <- derived_multiply_marked_dataset$female_adult_count/(derived_multiply_marked_dataset$transformed_mean_multiply_marked_total_egg_count/2)
derived_multiply_marked_male_viability <- derived_multiply_marked_dataset$male_adult_count/(derived_multiply_marked_dataset$transformed_mean_multiply_marked_total_egg_count/2)
derived_multiply_marked_total_viability <- derived_multiply_marked_dataset$total_adult_count/derived_multiply_marked_dataset$transformed_mean_multiply_marked_total_egg_count
multiply_marked_viability_dataset <-cbind(derived_multiply_marked_dataset,derived_multiply_marked_female_viability,derived_multiply_marked_male_viability,derived_multiply_marked_total_viability)
write.csv(multiply_marked_viability_dataset,"../Summaries and Outputs/6_multiply_marked_viability_dataset.csv",row.names = TRUE)

# Compute the overall mean viabilities for reporting in main text (excluding Broods G, H, and I where adult flies could not be scored)
overall_multiply_marked_viability <- mean(derived_multiply_marked_total_viability[-c(61:90)])
female_multiply_marked_viability <- mean(derived_multiply_marked_female_viability[-c(61:90)])
male_multiply_marked_viability <- mean(derived_multiply_marked_male_viability[-c(61:90)])

# Generate a data frame for sex-specific viability effects in the marker free crosses
female_vector <- rep("Female", 100)
male_vector <- rep("Male", 100)
multiply_marked_female_viabilities <- data.frame(female_vector,multiply_marked_viability_dataset$cross_id,multiply_marked_viability_dataset$brood_id,multiply_marked_viability_dataset$derived_multiply_marked_female_viability)
names(multiply_marked_female_viabilities) <- c("Sex","Cross","Brood","Viability")
multiply_marked_male_viabilities <- data.frame(male_vector,multiply_marked_viability_dataset$cross_id,multiply_marked_viability_dataset$brood_id,multiply_marked_viability_dataset$derived_multiply_marked_male_viability)
names(multiply_marked_male_viabilities) <- c("Sex","Cross","Brood","Viability")
multiply_marked_viability_anova_input <- rbind(multiply_marked_female_viabilities, multiply_marked_male_viabilities)

# Convert variable "Cross" to a categorical factor and "Brood" to a numeric value
multiply_marked_viability_anova_input$Cross <- factor(multiply_marked_viability_anova_input$Cross)
multiply_marked_viability_anova_input$Brood <- factor(multiply_marked_viability_anova_input$Brood)
multiply_marked_viability_anova_input$Brood <- as.numeric(multiply_marked_viability_anova_input$Brood)

# Remove Broods G, H, and I where adult flies could not be scored
cleaned_multiply_marked_viability_anova_input <- multiply_marked_viability_anova_input[-c(61:90,161:190), ]

# Create 6 single locus datasets as vectors
phenotypic_class_dataset <- multiply_marked_viability_dataset[,c(7:134)]
phenotypic_class_sums_vector <- colSums(phenotypic_class_dataset)
phenotypic_class_codes <- read.csv("phenotypic_class_codes.csv", header = TRUE)

counter <- 0
marker_vect <- c()
marker_matrix <- replicate(8,replicate(nrow(phenotypic_class_codes),0))

# Create dataframe of marker codes 
# Phenotypes : 0 - wildtype and 1 - mutant 
# Sex : 0 - female and 1 - male 
for (marker in phenotypic_class_codes$Marker.Binary) {
  counter <- counter + 1
  marker_code <- unlist(strsplit(marker, split = ""))
  if (marker_code[8] == "f") {
    sex <- 0
  } else {
    sex <- 1
  }
  marker_vect <- as.numeric(c(counter, marker_code[6:1], sex))
  marker_matrix[counter,] <- marker_vect
}
marker_df <- data.frame(marker_matrix)
colnames(marker_df) <- c("obsv", "scute", "crossveinless", "vermilion", "forked", "carnation", "yellow_plus", "sex")
phenotypes <- names(marker_df)
phenotypes <- phenotypes[-c(1,length(phenotypes))]

# Function to return a vector of the sums of each phenotype's counts at a single locus 
Single_Locus_Function <- function(x) {
  
  wt_female <- marker_df[marker_df[,x] == 0 & marker_df[,8] == 0, ]$obsv
  wt_male <- marker_df[marker_df[,x] == 0 & marker_df[,8] == 1, ]$obsv
  mutant_female <- marker_df[marker_df[,x] == 1 & marker_df[,8] == 0, ]$obsv
  mutant_male <- marker_df[marker_df[,x] == 1 & marker_df[,8] == 1, ]$obsv
  
  wt_female_sum <- sum(phenotypic_class_sums_vector[wt_female])
  wt_male_sum <- sum(phenotypic_class_sums_vector[wt_male])
  mutant_female_sum <- sum(phenotypic_class_sums_vector[mutant_female])
  mutant_male_sum <- sum(phenotypic_class_sums_vector[mutant_male])
  
  phenotype_single_locus_vector <- c(wt_female_sum, wt_male_sum, mutant_female_sum, mutant_male_sum, 4795)
  
  return(phenotype_single_locus_vector)
  
}

file_count <- 7
for (i in seq_len(length(phenotypes))) {
  phenotype_dataset <- Single_Locus_Function(i + 1)
  names(phenotype_dataset) <- c("female_wildtype", "male_wildtype", "female_mutant", "male_mutant", "lethal_zygote")
  file_name <- paste0("../Summaries and Outputs/", file_count, "_", phenotypes[i], "_single_locus_dataset.csv")
  write.csv(phenotype_dataset, file_name, row.names = TRUE)
  file_count <- file_count + 1
}

# Create input data file for fitting Cx(Co)M model to flies pooled by individual vials
# Remove Broods G, H, and I where adult flies could not be scored, remove extraneous meta-data columns
# Calculate the number of lethal zygotes by subtracting total adults from mean transformed egg count
# Add marker-free viabilities as classes 130, 131, 132, rename columns to conform to standard format
cleaned_derived_multiply_marked_dataset <- derived_multiply_marked_dataset[-c(61:90,161:190), ]
count_lethal_zygotes <- cleaned_derived_multiply_marked_dataset$transformed_mean_multiply_marked_total_egg_count - cleaned_derived_multiply_marked_dataset$total_adult_count
cleaned_derived_multiply_marked_dataset <- cbind(cleaned_derived_multiply_marked_dataset,count_lethal_zygotes)
pruned_derived_multiply_marked_dataset <- cleaned_derived_multiply_marked_dataset[,-c(1,2,3,5,6,135,136,137,138,139,140)]
total_viability_vector <- rep(0.921, 70)
female_viability_vector <- rep(0.921, 70)
male_viability_vector <- rep(0.921, 70)
multi_locus_individual_vials_dataset <- cbind(pruned_derived_multiply_marked_dataset,total_viability_vector,female_viability_vector,male_viability_vector)
names(multi_locus_individual_vials_dataset) <- c("ID", paste0("class_", 1:132))
write.csv(multi_locus_individual_vials_dataset,"../Summaries and Outputs/13_multi_locus_individual_vials_dataset.csv",row.names = TRUE)

# Create input data file for fitting Cx(Co)M model to flies pooled by replicate cross using individual_vials_dataset
unique_cross_id <- c("SK_14_1", "SK_14_3", "SK_14_5", "SK_14_6", "SK_14_8", "SK_14_12", "SK_14_14", "SK_14_17", "SK_14_18", "SK_14_19")
cross_pooling_matrix <- multi_locus_individual_vials_dataset[,-c(1,131,132,133)]

# Function to return a vector of multi-locus cross pooled data 
Cross_Pooling_Function <- function(x) {
  cross_pooling_rows <- seq( from = x, by = 10, length.out = 7)
  cross_pooling_vector <- colSums(cross_pooling_matrix[cross_pooling_rows, ])
  return(cross_pooling_vector)
} 

for (i in seq(1:10)) {
  assign(unique_cross_id[i], Cross_Pooling_Function(i))
}
multi_locus_cross_pooled_dataset <- rbind(SK_14_1, SK_14_3, SK_14_5, SK_14_6, SK_14_8, SK_14_12, SK_14_14, SK_14_17, SK_14_18, SK_14_19)
cross_viability_vector <- replicate(10, 0.921)
multi_locus_cross_pooled_dataset <- as.data.frame(cbind(unique_cross_id, multi_locus_cross_pooled_dataset, cross_viability_vector, cross_viability_vector, cross_viability_vector))
names(multi_locus_cross_pooled_dataset) <- c("ID", paste0("class_", 1:132))
write.csv(multi_locus_cross_pooled_dataset,"../Summaries and Outputs/14_multi_locus_cross_pooled_dataset.csv",row.names = TRUE)

# Create input data file for fitting Cx(Co)M model to flies pooled by brooding period using individual_vials_dataset
unique_brood_id <- c("SK_14_A", "SK_14_B", "SK_14_C", "SK_14_D", "SK_14_E", "SK_14_F", "SK_14_J")
brood_pooling_matrix <- multi_locus_individual_vials_dataset[,-c(1,131,132,133)]

# Function to return a vector of multi-locus brood pooled data 
Brood_Pooling_Function <- function(x) {
  brood_pooling_rows <- seq(from = 10*(x - 1) + 1, to = 10*x)
  brood_pooling_vector <- colSums(brood_pooling_matrix[brood_pooling_rows, ])
  return(brood_pooling_vector)
}

for (i in seq(1:7)) {
  assign(unique_brood_id[i], Brood_Pooling_Function(i))
}
multi_locus_brood_pooled_dataset <- rbind(SK_14_A, SK_14_B, SK_14_C, SK_14_D, SK_14_E, SK_14_F, SK_14_J)
brood_viability_vector <- replicate(7, 0.921)
multi_locus_brood_pooled_dataset <- as.data.frame(cbind(unique_brood_id, multi_locus_brood_pooled_dataset, brood_viability_vector, brood_viability_vector, brood_viability_vector))
names(multi_locus_brood_pooled_dataset) <- c("ID", paste0("class_", 1:132))
write.csv(multi_locus_brood_pooled_dataset,"15_multi_locus_brood_pooled_dataset.csv",row.names = TRUE)

# Create input data file for fitting Cx(Co)M model to flies pooled for the full experiment using individual_vials_dataset
full_experiment_matrix <- multi_locus_individual_vials_dataset[,-c(1,131,132,133)]
SK_14_rows <- c(1:70)
SK_14 <- colSums(full_experiment_matrix[SK_14_rows, ])
multi_locus_full_experiment_dataset <- rbind(SK_14, SK_14)
full_experiment_id <- c("full_experiment", "full_experiment")
full_experiment_vector <- c(0.921, 0.921)
multi_locus_full_experiment_dataset <- as.data.frame(cbind(full_experiment_id, multi_locus_full_experiment_dataset, full_experiment_vector, full_experiment_vector, full_experiment_vector))
multi_locus_full_experiment_dataset <- multi_locus_full_experiment_dataset[-2,]
names(multi_locus_full_experiment_dataset) <- c("ID", paste0("class_", 1:132))
write.csv(multi_locus_full_experiment_dataset,"../Summaries and Outputs/16_multi_locus_full_experiment_dataset.csv",row.names = TRUE)

# Create a vector of observed values for fitting CxCoM model
single_experimental_unit <- as.numeric(multi_locus_full_experiment_dataset)
observed_count <- c(single_experimental_unit[130],single_experimental_unit[2:129])

