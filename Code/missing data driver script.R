setwd("Raw Data/")

# Read in csv files of raw count data from the marker free cross
EH_marker_free_egg_count <- read.csv("Raw_Egg_Counts/A_EH_raw_marker_free_egg_counts.csv", header=TRUE)
MS_marker_free_egg_count <- read.csv("Raw_Egg_Counts/B_MS_raw_marker_free_egg_counts.csv", header=TRUE)
SK_marker_free_adult_count <- read.csv("Raw_Adult_Counts/C_SK_raw_marker_free_adult_counts.csv", header=TRUE)

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

# Perform Regression and write ANOVA table for supplemental material
marker_free_viability_anova.lm <- lm(Viability ~ Brood + Sex + Cross, data=marker_free_viability_anova_input)
marker_free_viability_anova.table <- anova(marker_free_viability_anova.lm)
marker_free_viability_regression_coefficients <- summary(marker_free_viability_anova.lm)
write.csv(marker_free_viability_anova.table,"../Summaries and Outputs/table_S1_marker_free_viability_regression.csv",row.names = TRUE)

# Read in csv files of raw count data from the multiply marked cross
EH_multiply_marked_egg_count <- read.csv("Raw_Egg_Counts/1_EH_raw_multiply_marked_egg_counts.csv", header=TRUE)
MS_multiply_marked_egg_count <- read.csv("Raw_Egg_Counts/2_MS_raw_multiply_marked_egg_counts.csv", header=TRUE)
SK_multiply_marked_adult_phenotypic_class_counts <- read.csv("Raw_Adult_Counts/3_SK_raw_multiply_marked_adult_phenotypic_class_counts.csv", header=TRUE)

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
write.csv(derived_multiply_marked_dataset,"/Users/spencerkoury/Desktop/local_missing_data/5_derived_multiply_marked_dataset.csv",row.names = TRUE)

# Calculate the per vial viability of females, males, and total under the 1:1 primary sex ratio assumption
derived_multiply_marked_female_viability <- derived_multiply_marked_dataset$female_adult_count/(derived_multiply_marked_dataset$transformed_mean_multiply_marked_total_egg_count/2)
derived_multiply_marked_male_viability <- derived_multiply_marked_dataset$male_adult_count/(derived_multiply_marked_dataset$transformed_mean_multiply_marked_total_egg_count/2)
derived_multiply_marked_total_viability <- derived_multiply_marked_dataset$total_adult_count/derived_multiply_marked_dataset$transformed_mean_multiply_marked_total_egg_count
multiply_marked_viability_dataset <-cbind(derived_multiply_marked_dataset,derived_multiply_marked_female_viability,derived_multiply_marked_male_viability,derived_multiply_marked_total_viability)
write.csv(multiply_marked_viability_dataset,"/Users/spencerkoury/Desktop/local_missing_data/6_multiply_marked_viability_dataset.csv",row.names = TRUE)

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

# Perform Regression and write ANOVA table for supplemental material
multiply_marked_viability_anova.lm <- lm(Viability ~ Brood + Sex + Cross, data=cleaned_multiply_marked_viability_anova_input)
multiply_marked_viability_anova.table <- anova(multiply_marked_viability_anova.lm)
multiply_marked_viability_regression_coefficients <- summary(multiply_marked_viability_anova.lm)
write.csv(multiply_marked_viability_anova.table,"/Users/spencerkoury/Desktop/local_missing_data/table_S2_multiply_marked_viability_regression.csv",row.names = TRUE)

#Create 6 single locus datasets as vectors
phenotypic_class_dataset <- multiply_marked_viability_dataset[,c(7:134)]
phenotypic_class_sums_vector <- colSums(phenotypic_class_dataset)

scute_single_locus_dataset <- c(sum(phenotypic_class_sums_vector[1],
                                    phenotypic_class_sums_vector[3],
                                    phenotypic_class_sums_vector[5],
                                    phenotypic_class_sums_vector[7],
                                    phenotypic_class_sums_vector[9],
                                    phenotypic_class_sums_vector[11],
                                    phenotypic_class_sums_vector[13],
                                    phenotypic_class_sums_vector[15],
                                    phenotypic_class_sums_vector[17],
                                    phenotypic_class_sums_vector[19],
                                    phenotypic_class_sums_vector[21],
                                    phenotypic_class_sums_vector[23],
                                    phenotypic_class_sums_vector[25],
                                    phenotypic_class_sums_vector[27],
                                    phenotypic_class_sums_vector[29],
                                    phenotypic_class_sums_vector[31],
                                    phenotypic_class_sums_vector[33],
                                    phenotypic_class_sums_vector[35],
                                    phenotypic_class_sums_vector[37],
                                    phenotypic_class_sums_vector[39],
                                    phenotypic_class_sums_vector[41],
                                    phenotypic_class_sums_vector[43],
                                    phenotypic_class_sums_vector[45],
                                    phenotypic_class_sums_vector[47],
                                    phenotypic_class_sums_vector[49],
                                    phenotypic_class_sums_vector[51],
                                    phenotypic_class_sums_vector[53],
                                    phenotypic_class_sums_vector[55],
                                    phenotypic_class_sums_vector[57],
                                    phenotypic_class_sums_vector[59],
                                    phenotypic_class_sums_vector[61],
                                    phenotypic_class_sums_vector[63]),
                                sum(phenotypic_class_sums_vector[2],
                                    phenotypic_class_sums_vector[4],
                                    phenotypic_class_sums_vector[6],
                                    phenotypic_class_sums_vector[8],
                                    phenotypic_class_sums_vector[10],
                                    phenotypic_class_sums_vector[12],
                                    phenotypic_class_sums_vector[14],
                                    phenotypic_class_sums_vector[16],
                                    phenotypic_class_sums_vector[18],
                                    phenotypic_class_sums_vector[20],
                                    phenotypic_class_sums_vector[22],
                                    phenotypic_class_sums_vector[24],
                                    phenotypic_class_sums_vector[26],
                                    phenotypic_class_sums_vector[28],
                                    phenotypic_class_sums_vector[30],
                                    phenotypic_class_sums_vector[32],
                                    phenotypic_class_sums_vector[34],
                                    phenotypic_class_sums_vector[36],
                                    phenotypic_class_sums_vector[38],
                                    phenotypic_class_sums_vector[40],
                                    phenotypic_class_sums_vector[42],
                                    phenotypic_class_sums_vector[44],
                                    phenotypic_class_sums_vector[46],
                                    phenotypic_class_sums_vector[48],
                                    phenotypic_class_sums_vector[50],
                                    phenotypic_class_sums_vector[52],
                                    phenotypic_class_sums_vector[54],
                                    phenotypic_class_sums_vector[56],
                                    phenotypic_class_sums_vector[58],
                                    phenotypic_class_sums_vector[60],
                                    phenotypic_class_sums_vector[62],
                                    phenotypic_class_sums_vector[64]),
                                sum(phenotypic_class_sums_vector[65],
                                    phenotypic_class_sums_vector[67],
                                    phenotypic_class_sums_vector[69],
                                    phenotypic_class_sums_vector[71],
                                    phenotypic_class_sums_vector[73],
                                    phenotypic_class_sums_vector[75],
                                    phenotypic_class_sums_vector[77],
                                    phenotypic_class_sums_vector[79],
                                    phenotypic_class_sums_vector[81],
                                    phenotypic_class_sums_vector[83],
                                    phenotypic_class_sums_vector[85],
                                    phenotypic_class_sums_vector[87],
                                    phenotypic_class_sums_vector[89],
                                    phenotypic_class_sums_vector[91],
                                    phenotypic_class_sums_vector[93],
                                    phenotypic_class_sums_vector[95],
                                    phenotypic_class_sums_vector[97],
                                    phenotypic_class_sums_vector[99],
                                    phenotypic_class_sums_vector[101],
                                    phenotypic_class_sums_vector[103],
                                    phenotypic_class_sums_vector[105],
                                    phenotypic_class_sums_vector[107],
                                    phenotypic_class_sums_vector[109],
                                    phenotypic_class_sums_vector[111],
                                    phenotypic_class_sums_vector[113],
                                    phenotypic_class_sums_vector[115],
                                    phenotypic_class_sums_vector[117],
                                    phenotypic_class_sums_vector[119],
                                    phenotypic_class_sums_vector[121],
                                    phenotypic_class_sums_vector[123],
                                    phenotypic_class_sums_vector[125],
                                    phenotypic_class_sums_vector[127]),
                                sum(phenotypic_class_sums_vector[66],
                                    phenotypic_class_sums_vector[68],
                                    phenotypic_class_sums_vector[70],
                                    phenotypic_class_sums_vector[72],
                                    phenotypic_class_sums_vector[74],
                                    phenotypic_class_sums_vector[76],
                                    phenotypic_class_sums_vector[78],
                                    phenotypic_class_sums_vector[80],
                                    phenotypic_class_sums_vector[82],
                                    phenotypic_class_sums_vector[84],
                                    phenotypic_class_sums_vector[86],
                                    phenotypic_class_sums_vector[88],
                                    phenotypic_class_sums_vector[90],
                                    phenotypic_class_sums_vector[92],
                                    phenotypic_class_sums_vector[94],
                                    phenotypic_class_sums_vector[96],
                                    phenotypic_class_sums_vector[98],
                                    phenotypic_class_sums_vector[100],
                                    phenotypic_class_sums_vector[102],
                                    phenotypic_class_sums_vector[104],
                                    phenotypic_class_sums_vector[106],
                                    phenotypic_class_sums_vector[108],
                                    phenotypic_class_sums_vector[110],
                                    phenotypic_class_sums_vector[112],
                                    phenotypic_class_sums_vector[114],
                                    phenotypic_class_sums_vector[116],
                                    phenotypic_class_sums_vector[118],
                                    phenotypic_class_sums_vector[120],
                                    phenotypic_class_sums_vector[122],
                                    phenotypic_class_sums_vector[124],
                                    phenotypic_class_sums_vector[126],
                                    phenotypic_class_sums_vector[128]),
                                4795)
names(scute_single_locus_dataset) <- c("female_wildtype","female_mutant","male_wildtype","male_mutant", "lethal_zygote")
write.csv(scute_single_locus_dataset,"/Users/spencerkoury/Desktop/local_missing_data/7_scute_single_locus_dataset.csv",row.names = TRUE)

# custom script for goodness-of-fit G-test for marker viability effects under six models of viability
# H0: no excess experimental mortality in multiply-marked crosses above and beyond the marker-free cross,
# H1: excess experimental mortality in multiply-marked crosses is random with respect to marker alleles,
# H2: excess experimental mortality in multiply-marked crosses is solely due to mutant alleles at marker loci, and
# H3: excess experimental mortality in multiply-marked crosses is due to viability effects associated with both mutant and wildtype alleles at marker loci.

# experimental mortality not associated with markers (v_1) is assumed to be 0.921 from external marker free experimental evidence
# with this assumption the marker associated effects can be estimated with method of least squares from system of linear equations
# without this assumption the marker associated viability effects must be solved for by method that permit system of nonlinear equations

# enter egg count sample size (n)
n <- 10952

# enter observed and unobserved counts
observed_female_wildtype_adult_count <- scute_single_locus_dataset[1]
observed_male_wildtype_adult_count   <- scute_single_locus_dataset[3]
observed_female_mutant_adult_count   <- scute_single_locus_dataset[2]
observed_male_mutant_adult_count     <- scute_single_locus_dataset[4]
observed_failed_development_count    <- scute_single_locus_dataset[5]

# solve linear system of equations to estimate viability effects
# H0: no experimental mortality (no coefficients to be estimated)
# H1: experimental mortality is random with respect to marker alleles (no coefficients to be estimated)
# H2: experimental mortality has both a random component and marker associated effects
H2_a <- c(0, 0, 2348.825, 2348.825, -4697.650)
H2_y <- c(observed_female_wildtype_adult_count-2348.825,
          observed_male_wildtype_adult_count-2348.825,
          observed_female_mutant_adult_count,
          observed_male_mutant_adult_count,
          observed_failed_development_count-5297.350)
H2_matrix <- as.data.frame(cbind(H2_y, H2_a))
H2_model <- lm(H2_y ~ H2_a, data=H2_matrix)
H2_vector <- as.vector(coef(H2_model))

# H3: experimental mortality has both a random component and sex-specific marker associated effects
H3_a <- c(0, 0, 2348.825, 0, -2348.825)
H3_b <- c(0, 0, 0, 2348.825, -2348.825)
H3_y <- c(observed_female_wildtype_adult_count-2348.825,
          observed_male_wildtype_adult_count-2348.825,
          observed_female_mutant_adult_count,
          observed_male_mutant_adult_count,
          observed_failed_development_count-5297.350)
H3_matrix <- as.data.frame(cbind(H3_y, H3_a, H3_b))
H3_model <- lm(H3_y ~ H3_a + H3_b, data=H3_matrix)
H3_vector <- as.vector(coef(H3_model))

# H4: experimental mortality has both a random compoment as well as sex and sex-specific marker associated effects
H4_a <- c(0, 2348.825, 0, 0, -2348.825)
H4_b <- c(0, 0, 2348.825, 0, -2348.825)
H4_c <- c(0, 0, 0, 2348.825, -2348.825)
H4_y <- c(observed_female_wildtype_adult_count-2348.825,
          observed_male_wildtype_adult_count,
          observed_female_mutant_adult_count,
          observed_male_mutant_adult_count,
          observed_failed_development_count-7646.175)
H4_matrix <- as.data.frame(cbind(H4_y, H4_a, H4_b, H4_c))
H4_model <- lm(H4_y ~ H4_a + H4_b + H4_c, data=H4_matrix)
H4_vector <- as.vector(coef(H4_model))

# H5: experimental mortality has both a random component and unique deterministic component for each adult phenotypic class
H5_a <- c(2348.825, 0, 0, 0, -2348.825)
H5_b <- c(0, 2348.825, 0, 0, -2348.825)
H5_c <- c(0, 0, 2348.825, 0, -2348.825)
H5_d <- c(0, 0, 0, 2348.825, -2348.825)
H5_y <- c(observed_female_wildtype_adult_count,
          observed_male_wildtype_adult_count,
          observed_female_mutant_adult_count,
          observed_male_mutant_adult_count,
          observed_failed_development_count-9995)
H5_matrix <- as.data.frame(cbind(H5_y, H5_a, H5_b, H5_c, H5_d))
H5_model <- lm(H5_y ~ H5_a + H5_b + H5_c + H5_d, data=H5_matrix)
H5_vector <- as.vector(coef(H5_model))

# save nls() outputs as viability effects for goodness-of-fit tests
H1_v_1 <- 0.94
H2_v_1 <- 0.94
H3_v_1 <- 0.94
H4_v_1 <- 0.94
H5_v_1 <- 0.94
H2_v_2 <- H2_vector[2]
H3_v_3 <- H3_vector[2]
H3_v_4 <- H3_vector[3]
H4_v_5 <- H4_vector[2]
H4_v_6 <- H4_vector[3]
H4_v_7 <- H4_vector[4]
H5_v_8 <- H5_vector[2]
H5_v_9 <- H5_vector[3]
H5_v_10<- H5_vector[4]
H5_v_11<- H5_vector[5]

# expected counts under specified hypotheses/viability models
# H0: no experimental mortality
H0_expected_female_wildtype_count <- n/4
H0_expected_male_wildtype_count   <- n/4
H0_expected_female_mutant_count   <- n/4
H0_expected_male_mutant_count     <- n/4
H0_expected_failed_development    <- 0

# H1: experimental mortality is random with respect to marker alleles
H1_expected_female_wildtype_count <- (n*H1_v_1)/4
H1_expected_male_wildtype_count   <- (n*H1_v_1)/4
H1_expected_female_mutant_count   <- (n*H1_v_1)/4
H1_expected_male_mutant_count     <- (n*H1_v_1)/4
H1_expected_failed_development    <- n*(1-H1_v_1)

# H2: experimental mortality has both a random component and marker associated effects
H2_expected_female_wildtype_count <- (n*H2_v_1)/4
H2_expected_male_wildtype_count   <- (n*H2_v_1)/4
H2_expected_female_mutant_count   <- (n*H2_v_1*H2_v_2)/4
H2_expected_male_mutant_count     <- (n*H2_v_1*H2_v_2)/4
H2_expected_failed_development    <- (n*(2-(H2_v_1*(1+H2_v_2))))/2

# H3: experimental mortality has both a random component and sex-specific marker associated effects
H3_expected_female_wildtype_count <- (n*H3_v_1)/4
H3_expected_male_wildtype_count   <- (n*H3_v_1)/4
H3_expected_female_mutant_count   <- (n*H3_v_1*H3_v_3)/4
H3_expected_male_mutant_count     <- (n*H3_v_1*H3_v_4)/4
H3_expected_failed_development    <- (n*(4-(H3_v_1*(2+H3_v_3+H3_v_4))))/4

# H4: experimental mortality has both a random compoment as well as sex and sex-specific marker associated effects
H4_expected_female_wildtype_count <- (n*H4_v_1)/4
H4_expected_male_wildtype_count   <- (n*H4_v_1*H4_v_5)/4
H4_expected_female_mutant_count   <- (n*H4_v_1*H4_v_6)/4
H4_expected_male_mutant_count     <- (n*H4_v_1*H4_v_7)/4
H4_expected_failed_development    <- (n*(4-(H4_v_1*(1+H4_v_5+H4_v_6+H4_v_7))))/4

# H5: experimental mortality has both a random component and unique deterministic component for each adult phenotypic class
H5_expected_female_wildtype_count <- (n*H5_v_1*H5_v_8)/4
H5_expected_male_wildtype_count   <- (n*H5_v_1*H5_v_9)/4
H5_expected_female_mutant_count   <- (n*H5_v_1*H5_v_10)/4
H5_expected_male_mutant_count     <- (n*H5_v_1*H5_v_11)/4
H5_expected_failed_development    <- (n*(4-(H5_v_1*(H5_v_8+H5_v_9+H5_v_10+H5_v_11))))/4

# Goodness-of-fit G-test for each hypothesis/viability model using parameters fit by method of least squares with nls()
G_statistic_H0 <- 2*sum(observed_female_wildtype_adult_count*(log((observed_female_wildtype_adult_count/H0_expected_female_wildtype_count))),
                        observed_male_wildtype_adult_count*(log((observed_male_wildtype_adult_count/H0_expected_male_wildtype_count))),
                        observed_female_mutant_adult_count*(log((observed_female_mutant_adult_count/H0_expected_female_mutant_count))),
                        observed_male_mutant_adult_count*(log((observed_male_mutant_adult_count/H0_expected_male_mutant_count))),
                        3838)

G_statistic_H1 <- 2*sum(observed_female_wildtype_adult_count*(log((observed_female_wildtype_adult_count/H1_expected_female_wildtype_count))),
                        observed_male_wildtype_adult_count*(log((observed_male_wildtype_adult_count/H1_expected_male_wildtype_count))),
                        observed_female_mutant_adult_count*(log((observed_female_mutant_adult_count/H1_expected_female_mutant_count))),
                        observed_male_mutant_adult_count*(log((observed_male_mutant_adult_count/H1_expected_male_mutant_count))),
                        observed_failed_development_count*(log((observed_failed_development_count/H1_expected_failed_development))))

G_statistic_H2 <- 2*sum(observed_female_wildtype_adult_count*(log((observed_female_wildtype_adult_count/H2_expected_female_wildtype_count))),
                        observed_male_wildtype_adult_count*(log((observed_male_wildtype_adult_count/H2_expected_male_wildtype_count))),
                        observed_female_mutant_adult_count*(log((observed_female_mutant_adult_count/H2_expected_female_mutant_count))),
                        observed_male_mutant_adult_count*(log((observed_male_mutant_adult_count/H2_expected_male_mutant_count))),
                        observed_failed_development_count*(log((observed_failed_development_count/H2_expected_failed_development))))

G_statistic_H3 <- 2*sum(observed_female_wildtype_adult_count*(log((observed_female_wildtype_adult_count/H3_expected_female_wildtype_count))),
                        observed_male_wildtype_adult_count*(log((observed_male_wildtype_adult_count/H3_expected_male_wildtype_count))),
                        observed_female_mutant_adult_count*(log((observed_female_mutant_adult_count/H3_expected_female_mutant_count))),
                        observed_male_mutant_adult_count*(log((observed_male_mutant_adult_count/H3_expected_male_mutant_count))),
                        observed_failed_development_count*(log((observed_failed_development_count/H3_expected_failed_development))))

G_statistic_H4 <- 2*sum(observed_female_wildtype_adult_count*(log((observed_female_wildtype_adult_count/H4_expected_female_wildtype_count))),
                        observed_male_wildtype_adult_count*(log((observed_male_wildtype_adult_count/H4_expected_male_wildtype_count))),
                        observed_female_mutant_adult_count*(log((observed_female_mutant_adult_count/H4_expected_female_mutant_count))),
                        observed_male_mutant_adult_count*(log((observed_male_mutant_adult_count/H4_expected_male_mutant_count))),
                        observed_failed_development_count*(log((observed_failed_development_count/H4_expected_failed_development))))

G_statistic_H5 <- 2*sum(observed_female_wildtype_adult_count*(log((observed_female_wildtype_adult_count/H5_expected_female_wildtype_count))),
                        observed_male_wildtype_adult_count*(log((observed_male_wildtype_adult_count/H5_expected_male_wildtype_count))),
                        observed_female_mutant_adult_count*(log((observed_female_mutant_adult_count/H5_expected_female_mutant_count))),
                        observed_male_mutant_adult_count*(log((observed_male_mutant_adult_count/H5_expected_male_mutant_count))),
                        observed_failed_development_count*(log((observed_failed_development_count/H5_expected_failed_development)))) 

H0_p_value <- pchisq(G_statistic_H0, 4, lower.tail=FALSE)
H1_p_value <- pchisq(G_statistic_H1, 4, lower.tail=FALSE)
H2_p_value <- pchisq(G_statistic_H2, 3, lower.tail=FALSE)
H3_p_value <- pchisq(G_statistic_H3, 2, lower.tail=FALSE)
H4_p_value <- pchisq(G_statistic_H4, 1, lower.tail=FALSE)
H5_p_value <- pchisq(G_statistic_H5, 0, lower.tail=FALSE)

G_statistic_H0
G_statistic_H1
G_statistic_H2
G_statistic_H3
G_statistic_H4
G_statistic_H5

H0_p_value
H1_p_value
H2_p_value
H3_p_value
H4_p_value
H5_p_value

H2_v_2
H3_v_3
H3_v_4
H4_v_5
H4_v_6
H4_v_7
H5_v_8
H5_v_9
H5_v_10
H5_v_11

crossveinless_single_locus_dataset <- c(sum(phenotypic_class_sums_vector[1],
                                           phenotypic_class_sums_vector[3],
                                           phenotypic_class_sums_vector[5],
                                           phenotypic_class_sums_vector[7],
                                           phenotypic_class_sums_vector[9],
                                           phenotypic_class_sums_vector[12],
                                           phenotypic_class_sums_vector[13],
                                           phenotypic_class_sums_vector[15],
                                           phenotypic_class_sums_vector[17],
                                           phenotypic_class_sums_vector[20],
                                           phenotypic_class_sums_vector[21],
                                           phenotypic_class_sums_vector[23],
                                           phenotypic_class_sums_vector[26],
                                           phenotypic_class_sums_vector[27],
                                           phenotypic_class_sums_vector[30],
                                           phenotypic_class_sums_vector[32],
                                           phenotypic_class_sums_vector[33],
                                           phenotypic_class_sums_vector[35],
                                           phenotypic_class_sums_vector[38],
                                           phenotypic_class_sums_vector[39],
                                           phenotypic_class_sums_vector[42],
                                           phenotypic_class_sums_vector[44],
                                           phenotypic_class_sums_vector[45],
                                           phenotypic_class_sums_vector[48],
                                           phenotypic_class_sums_vector[50],
                                           phenotypic_class_sums_vector[52],
                                           phenotypic_class_sums_vector[53],
                                           phenotypic_class_sums_vector[56],
                                           phenotypic_class_sums_vector[58],
                                           phenotypic_class_sums_vector[60],
                                           phenotypic_class_sums_vector[62],
                                           phenotypic_class_sums_vector[64]),
                                       sum(phenotypic_class_sums_vector[2],
                                           phenotypic_class_sums_vector[4],
                                           phenotypic_class_sums_vector[6],
                                           phenotypic_class_sums_vector[8],
                                           phenotypic_class_sums_vector[10],
                                           phenotypic_class_sums_vector[11],
                                           phenotypic_class_sums_vector[14],
                                           phenotypic_class_sums_vector[16],
                                           phenotypic_class_sums_vector[18],
                                           phenotypic_class_sums_vector[19],
                                           phenotypic_class_sums_vector[22],
                                           phenotypic_class_sums_vector[24],
                                           phenotypic_class_sums_vector[25],
                                           phenotypic_class_sums_vector[28],
                                           phenotypic_class_sums_vector[29],
                                           phenotypic_class_sums_vector[31],
                                           phenotypic_class_sums_vector[34],
                                           phenotypic_class_sums_vector[36],
                                           phenotypic_class_sums_vector[37],
                                           phenotypic_class_sums_vector[40],
                                           phenotypic_class_sums_vector[41],
                                           phenotypic_class_sums_vector[43],
                                           phenotypic_class_sums_vector[46],
                                           phenotypic_class_sums_vector[47],
                                           phenotypic_class_sums_vector[49],
                                           phenotypic_class_sums_vector[51],
                                           phenotypic_class_sums_vector[54],
                                           phenotypic_class_sums_vector[55],
                                           phenotypic_class_sums_vector[57],
                                           phenotypic_class_sums_vector[59],
                                           phenotypic_class_sums_vector[61],
                                           phenotypic_class_sums_vector[63]),
                                       sum(phenotypic_class_sums_vector[65],
                                           phenotypic_class_sums_vector[67],
                                           phenotypic_class_sums_vector[69],
                                           phenotypic_class_sums_vector[71],
                                           phenotypic_class_sums_vector[73],
                                           phenotypic_class_sums_vector[76],
                                           phenotypic_class_sums_vector[77],
                                           phenotypic_class_sums_vector[79],
                                           phenotypic_class_sums_vector[81],
                                           phenotypic_class_sums_vector[84],
                                           phenotypic_class_sums_vector[85],
                                           phenotypic_class_sums_vector[87],
                                           phenotypic_class_sums_vector[90],
                                           phenotypic_class_sums_vector[91],
                                           phenotypic_class_sums_vector[94],
                                           phenotypic_class_sums_vector[96],
                                           phenotypic_class_sums_vector[97],
                                           phenotypic_class_sums_vector[99],
                                           phenotypic_class_sums_vector[102],
                                           phenotypic_class_sums_vector[103],
                                           phenotypic_class_sums_vector[106],
                                           phenotypic_class_sums_vector[108],
                                           phenotypic_class_sums_vector[109],
                                           phenotypic_class_sums_vector[112],
                                           phenotypic_class_sums_vector[114],
                                           phenotypic_class_sums_vector[116],
                                           phenotypic_class_sums_vector[117],
                                           phenotypic_class_sums_vector[120],
                                           phenotypic_class_sums_vector[122],
                                           phenotypic_class_sums_vector[124],
                                           phenotypic_class_sums_vector[126],
                                           phenotypic_class_sums_vector[128]),
                                       sum(phenotypic_class_sums_vector[66],
                                           phenotypic_class_sums_vector[68],
                                           phenotypic_class_sums_vector[70],
                                           phenotypic_class_sums_vector[72],
                                           phenotypic_class_sums_vector[74],
                                           phenotypic_class_sums_vector[75],
                                           phenotypic_class_sums_vector[78],
                                           phenotypic_class_sums_vector[80],
                                           phenotypic_class_sums_vector[82],
                                           phenotypic_class_sums_vector[83],
                                           phenotypic_class_sums_vector[86],
                                           phenotypic_class_sums_vector[88],
                                           phenotypic_class_sums_vector[89],
                                           phenotypic_class_sums_vector[92],
                                           phenotypic_class_sums_vector[93],
                                           phenotypic_class_sums_vector[95],
                                           phenotypic_class_sums_vector[98],
                                           phenotypic_class_sums_vector[100],
                                           phenotypic_class_sums_vector[101],
                                           phenotypic_class_sums_vector[104],
                                           phenotypic_class_sums_vector[105],
                                           phenotypic_class_sums_vector[107],
                                           phenotypic_class_sums_vector[110],
                                           phenotypic_class_sums_vector[111],
                                           phenotypic_class_sums_vector[113],
                                           phenotypic_class_sums_vector[115],
                                           phenotypic_class_sums_vector[118],
                                           phenotypic_class_sums_vector[119],
                                           phenotypic_class_sums_vector[121],
                                           phenotypic_class_sums_vector[123],
                                           phenotypic_class_sums_vector[125],
                                           phenotypic_class_sums_vector[127]),
                                       4795)
names(crossveinless_single_locus_dataset) <- c("female_wildtype","female_mutant","male_wildtype","male_mutant", "lethal_zygote")
write.csv(crossveinless_single_locus_dataset,"/Users/spencerkoury/Desktop/local_missing_data/8_crossveinless_single_locus_dataset.csv",row.names = TRUE)

vermilion_single_locus_dataset <- c(sum(phenotypic_class_sums_vector[1],
                                         phenotypic_class_sums_vector[3],
                                         phenotypic_class_sums_vector[5],
                                         phenotypic_class_sums_vector[7],
                                         phenotypic_class_sums_vector[10],
                                         phenotypic_class_sums_vector[12],
                                         phenotypic_class_sums_vector[13],
                                         phenotypic_class_sums_vector[15],
                                         phenotypic_class_sums_vector[18],
                                         phenotypic_class_sums_vector[20],
                                         phenotypic_class_sums_vector[21],
                                         phenotypic_class_sums_vector[24],
                                         phenotypic_class_sums_vector[26],
                                         phenotypic_class_sums_vector[28],
                                         phenotypic_class_sums_vector[30],
                                         phenotypic_class_sums_vector[31],
                                         phenotypic_class_sums_vector[33],
                                         phenotypic_class_sums_vector[36],
                                         phenotypic_class_sums_vector[38],
                                         phenotypic_class_sums_vector[40],
                                         phenotypic_class_sums_vector[42],
                                         phenotypic_class_sums_vector[43],
                                         phenotypic_class_sums_vector[46],
                                         phenotypic_class_sums_vector[48],
                                         phenotypic_class_sums_vector[49],
                                         phenotypic_class_sums_vector[51],
                                         phenotypic_class_sums_vector[54],
                                         phenotypic_class_sums_vector[56],
                                         phenotypic_class_sums_vector[57],
                                         phenotypic_class_sums_vector[59],
                                         phenotypic_class_sums_vector[61],
                                         phenotypic_class_sums_vector[63]),
                                     sum(phenotypic_class_sums_vector[2],
                                         phenotypic_class_sums_vector[4],
                                         phenotypic_class_sums_vector[6],
                                         phenotypic_class_sums_vector[8],
                                         phenotypic_class_sums_vector[9],
                                         phenotypic_class_sums_vector[11],
                                         phenotypic_class_sums_vector[14],
                                         phenotypic_class_sums_vector[16],
                                         phenotypic_class_sums_vector[17],
                                         phenotypic_class_sums_vector[19],
                                         phenotypic_class_sums_vector[22],
                                         phenotypic_class_sums_vector[23],
                                         phenotypic_class_sums_vector[25],
                                         phenotypic_class_sums_vector[27],
                                         phenotypic_class_sums_vector[29],
                                         phenotypic_class_sums_vector[32],
                                         phenotypic_class_sums_vector[34],
                                         phenotypic_class_sums_vector[35],
                                         phenotypic_class_sums_vector[37],
                                         phenotypic_class_sums_vector[39],
                                         phenotypic_class_sums_vector[41],
                                         phenotypic_class_sums_vector[44],
                                         phenotypic_class_sums_vector[45],
                                         phenotypic_class_sums_vector[47],
                                         phenotypic_class_sums_vector[50],
                                         phenotypic_class_sums_vector[52],
                                         phenotypic_class_sums_vector[53],
                                         phenotypic_class_sums_vector[55],
                                         phenotypic_class_sums_vector[58],
                                         phenotypic_class_sums_vector[60],
                                         phenotypic_class_sums_vector[62],
                                         phenotypic_class_sums_vector[64]),
                                     sum(phenotypic_class_sums_vector[65],
                                         phenotypic_class_sums_vector[67],
                                         phenotypic_class_sums_vector[69],
                                         phenotypic_class_sums_vector[71],
                                         phenotypic_class_sums_vector[74],
                                         phenotypic_class_sums_vector[76],
                                         phenotypic_class_sums_vector[77],
                                         phenotypic_class_sums_vector[79],
                                         phenotypic_class_sums_vector[82],
                                         phenotypic_class_sums_vector[84],
                                         phenotypic_class_sums_vector[85],
                                         phenotypic_class_sums_vector[88],
                                         phenotypic_class_sums_vector[90],
                                         phenotypic_class_sums_vector[92],
                                         phenotypic_class_sums_vector[94],
                                         phenotypic_class_sums_vector[95],
                                         phenotypic_class_sums_vector[97],
                                         phenotypic_class_sums_vector[100],
                                         phenotypic_class_sums_vector[102],
                                         phenotypic_class_sums_vector[104],
                                         phenotypic_class_sums_vector[106],
                                         phenotypic_class_sums_vector[107],
                                         phenotypic_class_sums_vector[110],
                                         phenotypic_class_sums_vector[112],
                                         phenotypic_class_sums_vector[113],
                                         phenotypic_class_sums_vector[115],
                                         phenotypic_class_sums_vector[118],
                                         phenotypic_class_sums_vector[120],
                                         phenotypic_class_sums_vector[121],
                                         phenotypic_class_sums_vector[123],
                                         phenotypic_class_sums_vector[125],
                                         phenotypic_class_sums_vector[127]),
                                     sum(phenotypic_class_sums_vector[66],
                                         phenotypic_class_sums_vector[68],
                                         phenotypic_class_sums_vector[70],
                                         phenotypic_class_sums_vector[72],
                                         phenotypic_class_sums_vector[73],
                                         phenotypic_class_sums_vector[75],
                                         phenotypic_class_sums_vector[78],
                                         phenotypic_class_sums_vector[80],
                                         phenotypic_class_sums_vector[81],
                                         phenotypic_class_sums_vector[83],
                                         phenotypic_class_sums_vector[86],
                                         phenotypic_class_sums_vector[87],
                                         phenotypic_class_sums_vector[89],
                                         phenotypic_class_sums_vector[91],
                                         phenotypic_class_sums_vector[93],
                                         phenotypic_class_sums_vector[96],
                                         phenotypic_class_sums_vector[98],
                                         phenotypic_class_sums_vector[99],
                                         phenotypic_class_sums_vector[101],
                                         phenotypic_class_sums_vector[103],
                                         phenotypic_class_sums_vector[105],
                                         phenotypic_class_sums_vector[108],
                                         phenotypic_class_sums_vector[109],
                                         phenotypic_class_sums_vector[111],
                                         phenotypic_class_sums_vector[114],
                                         phenotypic_class_sums_vector[116],
                                         phenotypic_class_sums_vector[117],
                                         phenotypic_class_sums_vector[119],
                                         phenotypic_class_sums_vector[122],
                                         phenotypic_class_sums_vector[124],
                                         phenotypic_class_sums_vector[126],
                                         phenotypic_class_sums_vector[128]),
                                     4795)
names(vermilion_single_locus_dataset) <- c("female_wildtype","female_mutant","male_wildtype","male_mutant", "lethal_zygote")
write.csv(vermilion_single_locus_dataset,"/Users/spencerkoury/Desktop/local_missing_data/9_vermilion_single_locus_dataset.csv",row.names = TRUE)

forked_single_locus_dataset <-c(sum(phenotypic_class_sums_vector[1],
                                    phenotypic_class_sums_vector[3],
                                    phenotypic_class_sums_vector[5],
                                    phenotypic_class_sums_vector[8],
                                    phenotypic_class_sums_vector[10],
                                    phenotypic_class_sums_vector[12],
                                    phenotypic_class_sums_vector[13],
                                    phenotypic_class_sums_vector[16],
                                    phenotypic_class_sums_vector[18],
                                    phenotypic_class_sums_vector[20],
                                    phenotypic_class_sums_vector[22],
                                    phenotypic_class_sums_vector[24],
                                    phenotypic_class_sums_vector[26],
                                    phenotypic_class_sums_vector[27],
                                    phenotypic_class_sums_vector[29],
                                    phenotypic_class_sums_vector[31],
                                    phenotypic_class_sums_vector[34],
                                    phenotypic_class_sums_vector[36],
                                    phenotypic_class_sums_vector[38],
                                    phenotypic_class_sums_vector[39],
                                    phenotypic_class_sums_vector[41],
                                    phenotypic_class_sums_vector[43],
                                    phenotypic_class_sums_vector[45],
                                    phenotypic_class_sums_vector[47],
                                    phenotypic_class_sums_vector[49],
                                    phenotypic_class_sums_vector[52],
                                    phenotypic_class_sums_vector[53],
                                    phenotypic_class_sums_vector[55],
                                    phenotypic_class_sums_vector[57],
                                    phenotypic_class_sums_vector[60],
                                    phenotypic_class_sums_vector[62],
                                    phenotypic_class_sums_vector[64]),
                                sum(phenotypic_class_sums_vector[2],
                                    phenotypic_class_sums_vector[4],
                                    phenotypic_class_sums_vector[6],
                                    phenotypic_class_sums_vector[7],
                                    phenotypic_class_sums_vector[9],
                                    phenotypic_class_sums_vector[11],
                                    phenotypic_class_sums_vector[14],
                                    phenotypic_class_sums_vector[15],
                                    phenotypic_class_sums_vector[17],
                                    phenotypic_class_sums_vector[19],
                                    phenotypic_class_sums_vector[21],
                                    phenotypic_class_sums_vector[23],
                                    phenotypic_class_sums_vector[25],
                                    phenotypic_class_sums_vector[28],
                                    phenotypic_class_sums_vector[30],
                                    phenotypic_class_sums_vector[32],
                                    phenotypic_class_sums_vector[33],
                                    phenotypic_class_sums_vector[35],
                                    phenotypic_class_sums_vector[37],
                                    phenotypic_class_sums_vector[40],
                                    phenotypic_class_sums_vector[42],
                                    phenotypic_class_sums_vector[44],
                                    phenotypic_class_sums_vector[46],
                                    phenotypic_class_sums_vector[48],
                                    phenotypic_class_sums_vector[50],
                                    phenotypic_class_sums_vector[51],
                                    phenotypic_class_sums_vector[54],
                                    phenotypic_class_sums_vector[56],
                                    phenotypic_class_sums_vector[58],
                                    phenotypic_class_sums_vector[59],
                                    phenotypic_class_sums_vector[61],
                                    phenotypic_class_sums_vector[63]),
                                sum(phenotypic_class_sums_vector[65],
                                    phenotypic_class_sums_vector[67],
                                    phenotypic_class_sums_vector[69],
                                    phenotypic_class_sums_vector[72],
                                    phenotypic_class_sums_vector[74],
                                    phenotypic_class_sums_vector[76],
                                    phenotypic_class_sums_vector[77],
                                    phenotypic_class_sums_vector[80],
                                    phenotypic_class_sums_vector[82],
                                    phenotypic_class_sums_vector[84],
                                    phenotypic_class_sums_vector[86],
                                    phenotypic_class_sums_vector[88],
                                    phenotypic_class_sums_vector[90],
                                    phenotypic_class_sums_vector[91],
                                    phenotypic_class_sums_vector[93],
                                    phenotypic_class_sums_vector[95],
                                    phenotypic_class_sums_vector[98],
                                    phenotypic_class_sums_vector[100],
                                    phenotypic_class_sums_vector[102],
                                    phenotypic_class_sums_vector[103],
                                    phenotypic_class_sums_vector[105],
                                    phenotypic_class_sums_vector[107],
                                    phenotypic_class_sums_vector[109],
                                    phenotypic_class_sums_vector[111],
                                    phenotypic_class_sums_vector[113],
                                    phenotypic_class_sums_vector[116],
                                    phenotypic_class_sums_vector[117],
                                    phenotypic_class_sums_vector[119],
                                    phenotypic_class_sums_vector[121],
                                    phenotypic_class_sums_vector[124],
                                    phenotypic_class_sums_vector[126],
                                    phenotypic_class_sums_vector[128]),
                                sum(phenotypic_class_sums_vector[66],
                                    phenotypic_class_sums_vector[68],
                                    phenotypic_class_sums_vector[70],
                                    phenotypic_class_sums_vector[71],
                                    phenotypic_class_sums_vector[73],
                                    phenotypic_class_sums_vector[75],
                                    phenotypic_class_sums_vector[78],
                                    phenotypic_class_sums_vector[79],
                                    phenotypic_class_sums_vector[81],
                                    phenotypic_class_sums_vector[83],
                                    phenotypic_class_sums_vector[85],
                                    phenotypic_class_sums_vector[87],
                                    phenotypic_class_sums_vector[89],
                                    phenotypic_class_sums_vector[92],
                                    phenotypic_class_sums_vector[94],
                                    phenotypic_class_sums_vector[96],
                                    phenotypic_class_sums_vector[97],
                                    phenotypic_class_sums_vector[99],
                                    phenotypic_class_sums_vector[101],
                                    phenotypic_class_sums_vector[104],
                                    phenotypic_class_sums_vector[106],
                                    phenotypic_class_sums_vector[108],
                                    phenotypic_class_sums_vector[110],
                                    phenotypic_class_sums_vector[112],
                                    phenotypic_class_sums_vector[114],
                                    phenotypic_class_sums_vector[115],
                                    phenotypic_class_sums_vector[118],
                                    phenotypic_class_sums_vector[120],
                                    phenotypic_class_sums_vector[122],
                                    phenotypic_class_sums_vector[123],
                                    phenotypic_class_sums_vector[125],
                                    phenotypic_class_sums_vector[127]),
                                4795)
names(forked_single_locus_dataset) <- c("female_wildtype","female_mutant","male_wildtype","male_mutant", "lethal_zygote")
write.csv(forked_single_locus_dataset,"/Users/spencerkoury/Desktop/local_missing_data/10_forked_single_locus_dataset.csv",row.names = TRUE)

carnation_single_locus_dataset <- c(sum(phenotypic_class_sums_vector[1],
                                       phenotypic_class_sums_vector[3],
                                       phenotypic_class_sums_vector[6],
                                       phenotypic_class_sums_vector[8],
                                       phenotypic_class_sums_vector[10],
                                       phenotypic_class_sums_vector[12],
                                       phenotypic_class_sums_vector[14],
                                       phenotypic_class_sums_vector[16],
                                       phenotypic_class_sums_vector[18],
                                       phenotypic_class_sums_vector[20],
                                       phenotypic_class_sums_vector[21],
                                       phenotypic_class_sums_vector[23],
                                       phenotypic_class_sums_vector[25],
                                       phenotypic_class_sums_vector[27],
                                       phenotypic_class_sums_vector[29],
                                       phenotypic_class_sums_vector[31],
                                       phenotypic_class_sums_vector[33],
                                       phenotypic_class_sums_vector[35],
                                       phenotypic_class_sums_vector[37],
                                       phenotypic_class_sums_vector[39],
                                       phenotypic_class_sums_vector[41],
                                       phenotypic_class_sums_vector[43],
                                       phenotypic_class_sums_vector[46],
                                       phenotypic_class_sums_vector[48],
                                       phenotypic_class_sums_vector[50],
                                       phenotypic_class_sums_vector[52],
                                       phenotypic_class_sums_vector[54],
                                       phenotypic_class_sums_vector[56],
                                       phenotypic_class_sums_vector[58],
                                       phenotypic_class_sums_vector[60],
                                       phenotypic_class_sums_vector[61],
                                       phenotypic_class_sums_vector[63]),
                                   sum(phenotypic_class_sums_vector[2],
                                       phenotypic_class_sums_vector[4],
                                       phenotypic_class_sums_vector[5],
                                       phenotypic_class_sums_vector[7],
                                       phenotypic_class_sums_vector[9],
                                       phenotypic_class_sums_vector[11],
                                       phenotypic_class_sums_vector[13],
                                       phenotypic_class_sums_vector[15],
                                       phenotypic_class_sums_vector[17],
                                       phenotypic_class_sums_vector[19],
                                       phenotypic_class_sums_vector[22],
                                       phenotypic_class_sums_vector[24],
                                       phenotypic_class_sums_vector[26],
                                       phenotypic_class_sums_vector[28],
                                       phenotypic_class_sums_vector[30],
                                       phenotypic_class_sums_vector[32],
                                       phenotypic_class_sums_vector[34],
                                       phenotypic_class_sums_vector[36],
                                       phenotypic_class_sums_vector[38],
                                       phenotypic_class_sums_vector[40],
                                       phenotypic_class_sums_vector[42],
                                       phenotypic_class_sums_vector[44],
                                       phenotypic_class_sums_vector[45],
                                       phenotypic_class_sums_vector[47],
                                       phenotypic_class_sums_vector[49],
                                       phenotypic_class_sums_vector[51],
                                       phenotypic_class_sums_vector[53],
                                       phenotypic_class_sums_vector[55],
                                       phenotypic_class_sums_vector[57],
                                       phenotypic_class_sums_vector[59],
                                       phenotypic_class_sums_vector[62],
                                       phenotypic_class_sums_vector[64]),
                                   sum(phenotypic_class_sums_vector[65],
                                       phenotypic_class_sums_vector[67],
                                       phenotypic_class_sums_vector[70],
                                       phenotypic_class_sums_vector[72],
                                       phenotypic_class_sums_vector[74],
                                       phenotypic_class_sums_vector[76],
                                       phenotypic_class_sums_vector[78],
                                       phenotypic_class_sums_vector[80],
                                       phenotypic_class_sums_vector[82],
                                       phenotypic_class_sums_vector[84],
                                       phenotypic_class_sums_vector[85],
                                       phenotypic_class_sums_vector[87],
                                       phenotypic_class_sums_vector[89],
                                       phenotypic_class_sums_vector[91],
                                       phenotypic_class_sums_vector[93],
                                       phenotypic_class_sums_vector[95],
                                       phenotypic_class_sums_vector[97],
                                       phenotypic_class_sums_vector[99],
                                       phenotypic_class_sums_vector[101],
                                       phenotypic_class_sums_vector[103],
                                       phenotypic_class_sums_vector[105],
                                       phenotypic_class_sums_vector[107],
                                       phenotypic_class_sums_vector[110],
                                       phenotypic_class_sums_vector[112],
                                       phenotypic_class_sums_vector[114],
                                       phenotypic_class_sums_vector[116],
                                       phenotypic_class_sums_vector[118],
                                       phenotypic_class_sums_vector[120],
                                       phenotypic_class_sums_vector[122],
                                       phenotypic_class_sums_vector[124],
                                       phenotypic_class_sums_vector[125],
                                       phenotypic_class_sums_vector[127]),
                                   sum(phenotypic_class_sums_vector[66],
                                       phenotypic_class_sums_vector[68],
                                       phenotypic_class_sums_vector[69],
                                       phenotypic_class_sums_vector[71],
                                       phenotypic_class_sums_vector[73],
                                       phenotypic_class_sums_vector[75],
                                       phenotypic_class_sums_vector[77],
                                       phenotypic_class_sums_vector[79],
                                       phenotypic_class_sums_vector[81],
                                       phenotypic_class_sums_vector[83],
                                       phenotypic_class_sums_vector[86],
                                       phenotypic_class_sums_vector[88],
                                       phenotypic_class_sums_vector[90],
                                       phenotypic_class_sums_vector[92],
                                       phenotypic_class_sums_vector[94],
                                       phenotypic_class_sums_vector[96],
                                       phenotypic_class_sums_vector[98],
                                       phenotypic_class_sums_vector[100],
                                       phenotypic_class_sums_vector[102],
                                       phenotypic_class_sums_vector[104],
                                       phenotypic_class_sums_vector[106],
                                       phenotypic_class_sums_vector[108],
                                       phenotypic_class_sums_vector[109],
                                       phenotypic_class_sums_vector[111],
                                       phenotypic_class_sums_vector[113],
                                       phenotypic_class_sums_vector[115],
                                       phenotypic_class_sums_vector[117],
                                       phenotypic_class_sums_vector[119],
                                       phenotypic_class_sums_vector[121],
                                       phenotypic_class_sums_vector[123],
                                       phenotypic_class_sums_vector[126],
                                       phenotypic_class_sums_vector[128]),
                                   4795)
names(carnation_single_locus_dataset) <- c("female_wildtype","female_mutant","male_wildtype","male_mutant", "lethal_zygote")
write.csv(carnation_single_locus_dataset,"/Users/spencerkoury/Desktop/local_missing_data/11_carnation_single_locus_dataset.csv",row.names = TRUE)

yellow_plus_single_locus_dataset <- c(sum(phenotypic_class_sums_vector[1],
                                          phenotypic_class_sums_vector[4],
                                          phenotypic_class_sums_vector[6],
                                          phenotypic_class_sums_vector[8],
                                          phenotypic_class_sums_vector[10],
                                          phenotypic_class_sums_vector[12],
                                          phenotypic_class_sums_vector[13],
                                          phenotypic_class_sums_vector[15],
                                          phenotypic_class_sums_vector[17],
                                          phenotypic_class_sums_vector[19],
                                          phenotypic_class_sums_vector[21],
                                          phenotypic_class_sums_vector[23],
                                          phenotypic_class_sums_vector[25],
                                          phenotypic_class_sums_vector[27],
                                          phenotypic_class_sums_vector[29],
                                          phenotypic_class_sums_vector[31],
                                          phenotypic_class_sums_vector[34],
                                          phenotypic_class_sums_vector[36],
                                          phenotypic_class_sums_vector[38],
                                          phenotypic_class_sums_vector[40],
                                          phenotypic_class_sums_vector[42],
                                          phenotypic_class_sums_vector[44],
                                          phenotypic_class_sums_vector[46],
                                          phenotypic_class_sums_vector[48],
                                          phenotypic_class_sums_vector[50],
                                          phenotypic_class_sums_vector[52],
                                          phenotypic_class_sums_vector[53],
                                          phenotypic_class_sums_vector[55],
                                          phenotypic_class_sums_vector[57],
                                          phenotypic_class_sums_vector[59],
                                          phenotypic_class_sums_vector[61],
                                          phenotypic_class_sums_vector[64]),
                                      sum(phenotypic_class_sums_vector[2],
                                          phenotypic_class_sums_vector[3],
                                          phenotypic_class_sums_vector[5],
                                          phenotypic_class_sums_vector[7],
                                          phenotypic_class_sums_vector[9],
                                          phenotypic_class_sums_vector[11],
                                          phenotypic_class_sums_vector[14],
                                          phenotypic_class_sums_vector[16],
                                          phenotypic_class_sums_vector[18],
                                          phenotypic_class_sums_vector[20],
                                          phenotypic_class_sums_vector[22],
                                          phenotypic_class_sums_vector[24],
                                          phenotypic_class_sums_vector[26],
                                          phenotypic_class_sums_vector[28],
                                          phenotypic_class_sums_vector[30],
                                          phenotypic_class_sums_vector[32],
                                          phenotypic_class_sums_vector[33],
                                          phenotypic_class_sums_vector[35],
                                          phenotypic_class_sums_vector[37],
                                          phenotypic_class_sums_vector[39],
                                          phenotypic_class_sums_vector[41],
                                          phenotypic_class_sums_vector[43],
                                          phenotypic_class_sums_vector[45],
                                          phenotypic_class_sums_vector[47],
                                          phenotypic_class_sums_vector[49],
                                          phenotypic_class_sums_vector[51],
                                          phenotypic_class_sums_vector[54],
                                          phenotypic_class_sums_vector[56],
                                          phenotypic_class_sums_vector[58],
                                          phenotypic_class_sums_vector[60],
                                          phenotypic_class_sums_vector[62],
                                          phenotypic_class_sums_vector[63]),
                                      sum(phenotypic_class_sums_vector[65],
                                          phenotypic_class_sums_vector[67],
                                          phenotypic_class_sums_vector[70],
                                          phenotypic_class_sums_vector[72],
                                          phenotypic_class_sums_vector[74],
                                          phenotypic_class_sums_vector[76],
                                          phenotypic_class_sums_vector[78],
                                          phenotypic_class_sums_vector[80],
                                          phenotypic_class_sums_vector[82],
                                          phenotypic_class_sums_vector[84],
                                          phenotypic_class_sums_vector[85],
                                          phenotypic_class_sums_vector[87],
                                          phenotypic_class_sums_vector[89],
                                          phenotypic_class_sums_vector[91],
                                          phenotypic_class_sums_vector[93],
                                          phenotypic_class_sums_vector[95],
                                          phenotypic_class_sums_vector[97],
                                          phenotypic_class_sums_vector[99],
                                          phenotypic_class_sums_vector[101],
                                          phenotypic_class_sums_vector[103],
                                          phenotypic_class_sums_vector[105],
                                          phenotypic_class_sums_vector[107],
                                          phenotypic_class_sums_vector[110],
                                          phenotypic_class_sums_vector[112],
                                          phenotypic_class_sums_vector[114],
                                          phenotypic_class_sums_vector[116],
                                          phenotypic_class_sums_vector[118],
                                          phenotypic_class_sums_vector[120],
                                          phenotypic_class_sums_vector[122],
                                          phenotypic_class_sums_vector[124],
                                          phenotypic_class_sums_vector[125],
                                          phenotypic_class_sums_vector[127]),
                                      sum(phenotypic_class_sums_vector[66],
                                          phenotypic_class_sums_vector[67],
                                          phenotypic_class_sums_vector[69],
                                          phenotypic_class_sums_vector[71],
                                          phenotypic_class_sums_vector[73],
                                          phenotypic_class_sums_vector[75],
                                          phenotypic_class_sums_vector[78],
                                          phenotypic_class_sums_vector[80],
                                          phenotypic_class_sums_vector[82],
                                          phenotypic_class_sums_vector[84],
                                          phenotypic_class_sums_vector[86],
                                          phenotypic_class_sums_vector[88],
                                          phenotypic_class_sums_vector[90],
                                          phenotypic_class_sums_vector[92],
                                          phenotypic_class_sums_vector[94],
                                          phenotypic_class_sums_vector[96],
                                          phenotypic_class_sums_vector[97],
                                          phenotypic_class_sums_vector[99],
                                          phenotypic_class_sums_vector[101],
                                          phenotypic_class_sums_vector[103],
                                          phenotypic_class_sums_vector[105],
                                          phenotypic_class_sums_vector[107],
                                          phenotypic_class_sums_vector[109],
                                          phenotypic_class_sums_vector[111],
                                          phenotypic_class_sums_vector[113],
                                          phenotypic_class_sums_vector[115],
                                          phenotypic_class_sums_vector[118],
                                          phenotypic_class_sums_vector[120],
                                          phenotypic_class_sums_vector[122],
                                          phenotypic_class_sums_vector[124],
                                          phenotypic_class_sums_vector[126],
                                          phenotypic_class_sums_vector[127]),
                                      4795)
names(yellow_plus_single_locus_dataset) <- c("female_wildtype","female_mutant","male_wildtype","male_mutant", "lethal_zygote")
write.csv(yellow_plus_single_locus_dataset,"/Users/spencerkoury/Desktop/local_missing_data/12_yellow_plus_single_locus_dataset.csv",row.names = TRUE)

