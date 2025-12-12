setwd("...") # Please use Raw Data as working directory

# Load package dfoptim
library(dfoptim)

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

# Perform Regression and write ANOVA table for supplemental material
marker_free_viability_anova.lm <- lm(Viability ~ Brood + Sex + Cross, data=marker_free_viability_anova_input)
marker_free_viability_anova.table <- anova(marker_free_viability_anova.lm)
marker_free_viability_regression_coefficients <- summary(marker_free_viability_anova.lm)
write.csv(marker_free_viability_anova.table,"table_S1_marker_free_viability_regression.csv",row.names = TRUE)

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

# Perform Regression and write ANOVA table for supplemental material
multiply_marked_viability_anova.lm <- lm(Viability ~ Brood + Sex + Cross, data=cleaned_multiply_marked_viability_anova_input)
multiply_marked_viability_anova.table <- anova(multiply_marked_viability_anova.lm)
multiply_marked_viability_regression_coefficients <- summary(multiply_marked_viability_anova.lm)
write.csv(multiply_marked_viability_anova.table,"../Summaries and Outputs/table_S2_multiply_marked_viability_regression.csv",row.names = TRUE)

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

# Write csv for the single locus dataset
names(scute_single_locus_dataset) <- c("female_wildtype", "male_wildtype", "female_mutant", "male_mutant", "lethal_zygote")
write.csv(scute_single_locus_dataset,"../Summaries and Outputs/7_scute_single_locus_dataset.csv",row.names = TRUE)

# Perform single-locus goodness-of-fit G-tests for scute calling custom function from "single_locus_analysis_custom_function.R"
single_locus_input <- scute_single_locus_dataset
scute_G_test_table <- single_locus_analysis()
write.csv(scute_G_test_table,"../Summaries and Outputs/table_S3_scute_single_locus_G_tests.csv",row.names = TRUE)

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
names(crossveinless_single_locus_dataset) <- c("female_wildtype", "male_wildtype", "female_mutant", "male_mutant", "lethal_zygote")
write.csv(crossveinless_single_locus_dataset,"../Summaries and Outputs/8_crossveinless_single_locus_dataset.csv",row.names = TRUE)

# Perform single-locus goodness-of-fit G-tests for crossveinless calling custom function from "single_locus_analysis_custom_function.R"
single_locus_input <- crossveinless_single_locus_dataset
crossveinless_G_test_table <- single_locus_analysis()
write.csv(crossveinless_G_test_table,"../Summaries and Outputs/table_S4_crossveinless_single_locus_G_tests.csv",row.names = TRUE)

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
names(vermilion_single_locus_dataset) <- c("female_wildtype", "male_wildtype", "female_mutant", "male_mutant", "lethal_zygote")
write.csv(vermilion_single_locus_dataset,"../Summaries and Outputs/9_vermilion_single_locus_dataset.csv",row.names = TRUE)

# Perform single-locus goodness-of-fit G-tests for vermilion calling custom function from "single_locus_analysis_custom_function.R"
single_locus_input <- vermilion_single_locus_dataset
vermilion_G_test_table <- single_locus_analysis()
write.csv(vermilion_G_test_table,"../Summaries and Outputs/table_S5_vermilion_single_locus_G_tests.csv",row.names = TRUE)

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
names(forked_single_locus_dataset) <- c("female_wildtype", "male_wildtype", "female_mutant", "male_mutant", "lethal_zygote")
write.csv(forked_single_locus_dataset,"../Summaries and Outputs/10_forked_single_locus_dataset.csv",row.names = TRUE)

# Perform single-locus goodness-of-fit G-tests for forked calling custom function from "single_locus_analysis_custom_function.R"
single_locus_input <- forked_single_locus_dataset
forked_G_test_table <- single_locus_analysis()
write.csv(forked_G_test_table,"../Summaries and Outputs/table_S6_forked_single_locus_G_tests.csv",row.names = TRUE)

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
names(carnation_single_locus_dataset) <- c("female_wildtype", "male_wildtype", "female_mutant", "male_mutant", "lethal_zygote")
write.csv(carnation_single_locus_dataset,"../Summaries and Outputs/11_carnation_single_locus_dataset.csv",row.names = TRUE)

# Perform single-locus goodness-of-fit G-tests for carnation calling custom function from "single_locus_analysis_custom_function.R"
single_locus_input <- carnation_single_locus_dataset
carnation_G_test_table <- single_locus_analysis()
write.csv(carnation_G_test_table,"../Summaries and Outputs/table_S7_carnation_single_locus_G_tests.csv",row.names = TRUE)

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
names(yellow_plus_single_locus_dataset) <- c("female_wildtype", "male_wildtype", "female_mutant", "male_mutant", "lethal_zygote")
write.csv(yellow_plus_single_locus_dataset,"../Summaries and Outputs/12_yellow_plus_single_locus_dataset.csv",row.names = TRUE)

# Perform single-locus goodness-of-fit G-tests for yellow plus calling custom function from "single_locus_analysis_custom_function.R"
single_locus_input <- yellow_plus_single_locus_dataset
yellow_plus_G_test_table <- single_locus_analysis()
write.csv(yellow_plus_G_test_table,"../Summaries and Outputs/table_S8_yellow_plus_single_locus_G_tests.csv",row.names = TRUE)

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
cross_pooling_matrix <- multi_locus_individual_vials_dataset[,-c(1,131,132,133)]
SK_14_1_rows <- c(1, 11, 21, 31, 41, 51, 61)
SK_14_1 <- colSums(cross_pooling_matrix[SK_14_1_rows, ])
SK_14_3_rows <- c(2, 12, 22, 32, 42, 52, 62)
SK_14_3 <- colSums(cross_pooling_matrix[SK_14_3_rows, ])
SK_14_5_rows <- c(3, 13, 23, 33, 43, 53, 63)
SK_14_5 <- colSums(cross_pooling_matrix[SK_14_5_rows, ])
SK_14_6_rows <- c(4, 14, 24, 34, 44, 54, 64)
SK_14_6 <- colSums(cross_pooling_matrix[SK_14_6_rows, ])
SK_14_8_rows <- c(5, 15, 25, 35, 45, 55, 65)
SK_14_8 <- colSums(cross_pooling_matrix[SK_14_8_rows, ])
SK_14_12_rows <- c(6, 16, 26, 36, 46, 56, 66)
SK_14_12 <- colSums(cross_pooling_matrix[SK_14_12_rows, ])
SK_14_14_rows <- c(7, 17, 27, 37, 47, 57, 67)
SK_14_14 <- colSums(cross_pooling_matrix[SK_14_14_rows, ])
SK_14_17_rows <- c(8, 18, 28, 38, 48, 58, 68)
SK_14_17 <- colSums(cross_pooling_matrix[SK_14_17_rows, ])
SK_14_18_rows <- c(9, 19, 29, 39, 49, 59, 69)
SK_14_18 <- colSums(cross_pooling_matrix[SK_14_18_rows, ])
SK_14_19_rows <- c(10, 20, 30, 40, 50, 60, 70)
SK_14_19 <- colSums(cross_pooling_matrix[SK_14_19_rows, ])
multi_locus_cross_pooled_dataset <- rbind(SK_14_1, SK_14_3, SK_14_5, SK_14_6, SK_14_8, SK_14_12, SK_14_14, SK_14_17, SK_14_18, SK_14_19)
unique_cross_id <- c("SK_14_1", "SK_14_3", "SK_14_5", "SK_14_6", "SK_14_8", "SK_14_12", "SK_14_14", "SK_14_17", "SK_14_18", "SK_14_19")
cross_viability_vector <- c(0.921, 0.921, 0.921, 0.921, 0.921, 0.921, 0.921, 0.921, 0.921, 0.921)
multi_locus_cross_pooled_dataset <- as.data.frame(cbind(unique_cross_id, multi_locus_cross_pooled_dataset, cross_viability_vector, cross_viability_vector, cross_viability_vector))
names(multi_locus_cross_pooled_dataset) <- c("ID", paste0("class_", 1:132))
write.csv(multi_locus_cross_pooled_dataset,"../Summaries and Outputs/14_multi_locus_cross_pooled_dataset.csv",row.names = TRUE)

# Create input data file for fitting Cx(Co)M model to flies pooled by brooding period using individual_vials_dataset
brood_pooling_matrix <- multi_locus_individual_vials_dataset[,-c(1,131,132,133)]
SK_14_A_rows <- c(1:10)
SK_14_A <- colSums(brood_pooling_matrix[SK_14_A_rows, ])
SK_14_B_rows <- c(11:20)
SK_14_B <- colSums(brood_pooling_matrix[SK_14_B_rows, ])
SK_14_C_rows <- c(21:30)
SK_14_C <- colSums(brood_pooling_matrix[SK_14_C_rows, ])
SK_14_D_rows <- c(31:40)
SK_14_D <- colSums(brood_pooling_matrix[SK_14_D_rows, ])
SK_14_E_rows <- c(41:50)
SK_14_E <- colSums(brood_pooling_matrix[SK_14_E_rows, ])
SK_14_F_rows <- c(51:60)
SK_14_F <- colSums(brood_pooling_matrix[SK_14_F_rows, ])
SK_14_J_rows <- c(61:70)
SK_14_J <- colSums(brood_pooling_matrix[SK_14_J_rows, ])
multi_locus_brood_pooled_dataset <- rbind(SK_14_A, SK_14_B, SK_14_C, SK_14_D, SK_14_E, SK_14_F, SK_14_J)
unique_brood_id <- c("SK_14_A", "SK_14_B", "SK_14_C", "SK_14_D", "SK_14_E", "SK_14_F", "SK_14_J")
brood_viability_vector <- c(0.921, 0.921, 0.921, 0.921, 0.921, 0.921, 0.921)
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

# Begin algorithm for fitting CxCoM H0 to full experiment dataset
# Create matrix to store all optimized outputs of CxCoM.Likelihood.H0
# Log likelihoods calculated for each combination of y & p model parameters
MLE_H0_full_outputs <- matrix(nrow=0, ncol=39)
colnames(MLE_H0_full_outputs)<-c("ID",
                                "MLE_X_length",
                                "MLE_x1",
                                "MLE_x2",
                                "MLE_x3",
                                "MLE_x4",
                                "MLE_x5",
                                "MLE_p",
                                "MLE_y1",
                                "MLE_y2",
                                "MLE_y3",
                                "MLE_y4",
                                "MLE_y5",
                                "MLE_v1f+",
                                "MLE_v2f+",
                                "MLE_v3f+",
                                "MLE_v4f+",
                                "MLE_v5f+",
                                "MLE_v6f+",
                                "MLE_v1f-",
                                "MLE_v2f-",
                                "MLE_v3f-",
                                "MLE_v4f-",
                                "MLE_v5f-",
                                "MLE_v6f-",
                                "MLE_v1m+",
                                "MLE_v2m+",
                                "MLE_v3m+",
                                "MLE_v4m+",
                                "MLE_v5m+",
                                "MLE_v6m+",
                                "MLE_v1m-",
                                "MLE_v2m-",
                                "MLE_v3m-",
                                "MLE_v4m-",
                                "MLE_v5m-",
                                "MLE_v6m-",
                                "lnL",
                                "Convergence")


# For loop Uniform.CxCoM.Likelihood for all possible combinations of p
# Output model parameters and likelihoods in matrix MLE_H0_full_outputs
for (v in 1:max_p_value) {
  
  # Define parameters starting points for y, v, and p in single run of Nelder-Mead algorithm
  # When bounding Nelder-Mead using package dfoptim, starting points cannot be at the boundaries
  parameters <- as.numeric(c(y_uniform[v,1],
                             y_uniform[v,2],
                             y_uniform[v,3],
                             y_uniform[v,4],
                             y_uniform[v,5]))
  p <- p_uniform[v,]
  
  # Optimize CxCoM model parameters y and v with bounded Nelder-Mead method
  # Upper limit on y's is arbitrary, whereas upper and lower limit on v's is
  # given by definition of egg-to-adult viability (cannot exceed 0-1 range)
  MLE_H0_Full_CxCoM_Table <- nmkb(par = parameters,
                              fn = CxCoM.Likelihood.H0,
                              control = list(maxfeval = 50000),
                              lower = c(0,0,0,0,0),
                              upper = c(9,9,9,9,9))
  
  # Create output vector of y, p, likelihood, convergence code
  MLE_H0_outputs <- c(0,
                  (((MLE_H0_Full_CxCoM_Table$par[1]+MLE_H0_Full_CxCoM_Table$par[2]+MLE_H0_Full_CxCoM_Table$par[3]+MLE_H0_Full_CxCoM_Table$par[4]+MLE_H0_Full_CxCoM_Table$par[5])/(2*p[1]))*100),
                   MLE_H0_Full_CxCoM_Table$par[1]/(2*p[1]),
                   MLE_H0_Full_CxCoM_Table$par[2]/(2*p[1]),
                   MLE_H0_Full_CxCoM_Table$par[3]/(2*p[1]),
                   MLE_H0_Full_CxCoM_Table$par[4]/(2*p[1]),
                   MLE_H0_Full_CxCoM_Table$par[5]/(2*p[1]),
                   p[1],
                   MLE_H0_Full_CxCoM_Table$par[1],
                   MLE_H0_Full_CxCoM_Table$par[2],
                   MLE_H0_Full_CxCoM_Table$par[3],
                   MLE_H0_Full_CxCoM_Table$par[4],
                   MLE_H0_Full_CxCoM_Table$par[5],
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   MLE_H0_Full_CxCoM_Table$value,
                   MLE_H0_Full_CxCoM_Table$convergence)
  
  # Add vector of y, p, likelihood, code to MLE_Uniform_CxCoM
  MLE_H0_full_outputs <- rbind(MLE_H0_full_outputs, MLE_H0_outputs)
}

# Take output matrix MLE_Uniform_CxCoM and sort by likelihood
# Print the sorted table for records, make sure to RENAME FILE
MLE_H0_full<-as.data.frame(MLE_H0_full_outputs)
df0 <- MLE_H0_full[order(MLE_H0_full$lnL),]
write.csv(df0,"../Summaries and Outputs/17_multi_locus_full_experiment_H0_mle_output.csv",row.names = TRUE)


# Begin algorithm for fitting CxCoM H1 to full experiment dataset
# Create matrix to store all optimized outputs of CxCoM.Likelihood.H1
# Log likelihoods calculated for each combination of y & p model parameters
MLE_H1_full_outputs <- matrix(nrow=0, ncol=39)
colnames(MLE_H1_full_outputs)<-c("ID",
                                 "MLE_X_length",
                                 "MLE_x1",
                                 "MLE_x2",
                                 "MLE_x3",
                                 "MLE_x4",
                                 "MLE_x5",
                                 "MLE_p",
                                 "MLE_y1",
                                 "MLE_y2",
                                 "MLE_y3",
                                 "MLE_y4",
                                 "MLE_y5",
                                 "MLE_v1f+",
                                 "MLE_v2f+",
                                 "MLE_v3f+",
                                 "MLE_v4f+",
                                 "MLE_v5f+",
                                 "MLE_v6f+",
                                 "MLE_v1f-",
                                 "MLE_v2f-",
                                 "MLE_v3f-",
                                 "MLE_v4f-",
                                 "MLE_v5f-",
                                 "MLE_v6f-",
                                 "MLE_v1m+",
                                 "MLE_v2m+",
                                 "MLE_v3m+",
                                 "MLE_v4m+",
                                 "MLE_v5m+",
                                 "MLE_v6m+",
                                 "MLE_v1m-",
                                 "MLE_v2m-",
                                 "MLE_v3m-",
                                 "MLE_v4m-",
                                 "MLE_v5m-",
                                 "MLE_v6m-",
                                 "lnL",
                                 "Convergence")


# For loop Uniform.CxCoM.Likelihood for all possible combinations of p
# Output model parameters and likelihoods in matrix MLE_H1_full_outputs
for (v in 1:max_p_value) {
  
  # Define parameters starting points for y, v, and p in single run of Nelder-Mead algorithm
  # When bounding Nelder-Mead using package dfoptim, starting points cannot be at the boundaries
  parameters <- as.numeric(c(y_uniform[v,1],
                             y_uniform[v,2],
                             y_uniform[v,3],
                             y_uniform[v,4],
                             y_uniform[v,5],
                             0.5,
                             0.5))
                             
  p <- p_uniform[v,]
  
  # Optimize CxCoM model parameters y and v with bounded Nelder-Mead method
  # Upper limit on y's is arbitrary, whereas upper and lower limit on v's is
  # given by definition of egg-to-adult viability (cannot exceed 0-1 range)
  MLE_H1_Full_CxCoM_Table <- nmkb(par = parameters,
                                  fn = CxCoM.Likelihood.H1,
                                  control = list(maxfeval = 50000),
                                  lower = c(0,0,0,0,0,0,0),
                                  upper = c(9,9,9,9,9,1,1))
  
  # Create output vector of y, p, likelihood, convergence code
  MLE_H1_outputs <- c(1,
                   (((MLE_H1_Full_CxCoM_Table$par[1]+MLE_H1_Full_CxCoM_Table$par[2]+MLE_H1_Full_CxCoM_Table$par[3]+MLE_H1_Full_CxCoM_Table$par[4]+MLE_H1_Full_CxCoM_Table$par[5])/(2*p[1]))*100),
                   MLE_H1_Full_CxCoM_Table$par[1]/(2*p[1]),
                   MLE_H1_Full_CxCoM_Table$par[2]/(2*p[1]),
                   MLE_H1_Full_CxCoM_Table$par[3]/(2*p[1]),
                   MLE_H1_Full_CxCoM_Table$par[4]/(2*p[1]),
                   MLE_H1_Full_CxCoM_Table$par[5]/(2*p[1]),
                   p[1],
                   MLE_H1_Full_CxCoM_Table$par[1],
                   MLE_H1_Full_CxCoM_Table$par[2],
                   MLE_H1_Full_CxCoM_Table$par[3],
                   MLE_H1_Full_CxCoM_Table$par[4],
                   MLE_H1_Full_CxCoM_Table$par[5],
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   MLE_H1_Full_CxCoM_Table$value,
                   MLE_H1_Full_CxCoM_Table$convergence)
  
  # Add vector of y, p, likelihood, code to MLE_Uniform_CxCoM
  MLE_H1_full_outputs <- rbind(MLE_H1_full_outputs, MLE_H1_outputs)
}

# Take output matrix MLE_Uniform_CxCoM and sort by likelihood
# Print the sorted table for records, make sure to RENAME FILE
MLE_H1_full<-as.data.frame(MLE_H1_full_outputs)
df1 <- MLE_H1_full[order(MLE_H1_full$lnL),]
write.csv(df1,"../Summaries and Outputs/18_multi_locus_full_experiment_H1_mle_output.csv",row.names = TRUE)

# Begin algorithm for fitting CxCoM H2 to full experiment dataset
# Create matrix to store all optimized outputs of CxCoM.Likelihood.H2
# Log likelihoods calculated for each combination of y & p model parameters
MLE_H2_full_outputs <- matrix(nrow=0, ncol=39)
colnames(MLE_H2_full_outputs)<-c("ID",
                                 "MLE_X_length",
                                 "MLE_x1",
                                 "MLE_x2",
                                 "MLE_x3",
                                 "MLE_x4",
                                 "MLE_x5",
                                 "MLE_p",
                                 "MLE_y1",
                                 "MLE_y2",
                                 "MLE_y3",
                                 "MLE_y4",
                                 "MLE_y5",
                                 "MLE_v1f+",
                                 "MLE_v2f+",
                                 "MLE_v3f+",
                                 "MLE_v4f+",
                                 "MLE_v5f+",
                                 "MLE_v6f+",
                                 "MLE_v1f-",
                                 "MLE_v2f-",
                                 "MLE_v3f-",
                                 "MLE_v4f-",
                                 "MLE_v5f-",
                                 "MLE_v6f-",
                                 "MLE_v1m+",
                                 "MLE_v2m+",
                                 "MLE_v3m+",
                                 "MLE_v4m+",
                                 "MLE_v5m+",
                                 "MLE_v6m+",
                                 "MLE_v1m-",
                                 "MLE_v2m-",
                                 "MLE_v3m-",
                                 "MLE_v4m-",
                                 "MLE_v5m-",
                                 "MLE_v6m-",
                                 "lnL",
                                 "Convergence")


# For loop Uniform.CxCoM.Likelihood for all possible combinations of p
# Output model parameters and likelihoods in matrix MLE_H2_full_outputs
for (v in 1:max_p_value) {
  
  # Define parameters starting points for y, v, and p in single run of Nelder-Mead algorithm
  # When bounding Nelder-Mead using package dfoptim, starting points cannot be at the boundaries
  parameters <- as.numeric(c(y_uniform[v,1],
                             y_uniform[v,2],
                             y_uniform[v,3],
                             y_uniform[v,4],
                             y_uniform[v,5],
                             0.619,
                             0.615,
                             0.617,
                             0.584,
                             0.540,
                             0.591,
                             0.479,
                             0.479,
                             0.466,
                             0.412,
                             0.337,
                             0.428))
  
  p <- p_uniform[v,]
  
  # Optimize CxCoM model parameters y and v with bounded Nelder-Mead method
  # Upper limit on y's is arbitrary, whereas upper and lower limit on v's is
  # given by definition of egg-to-adult viability (cannot exceed 0-1 range)
  MLE_H2_Full_CxCoM_Table <- nmkb(par = parameters,
                                  fn = CxCoM.Likelihood.H2,
                                  control = list(maxfeval = 50000),
                                  lower = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                  upper = c(9,9,9,9,9,1,1,1,1,1,1,1,1,1,1,1,1))
  
  # Create output vector of y, p, likelihood, convergence code
  MLE_H2_outputs <- c(2,
                      (((MLE_H2_Full_CxCoM_Table$par[1]+MLE_H2_Full_CxCoM_Table$par[2]+MLE_H2_Full_CxCoM_Table$par[3]+MLE_H2_Full_CxCoM_Table$par[4]+MLE_H2_Full_CxCoM_Table$par[5])/(2*p[1]))*100),
                      MLE_H2_Full_CxCoM_Table$par[1]/(2*p[1]),
                      MLE_H2_Full_CxCoM_Table$par[2]/(2*p[1]),
                      MLE_H2_Full_CxCoM_Table$par[3]/(2*p[1]),
                      MLE_H2_Full_CxCoM_Table$par[4]/(2*p[1]),
                      MLE_H2_Full_CxCoM_Table$par[5]/(2*p[1]),
                      p[1],
                      MLE_H2_Full_CxCoM_Table$par[1],
                      MLE_H2_Full_CxCoM_Table$par[2],
                      MLE_H2_Full_CxCoM_Table$par[3],
                      MLE_H2_Full_CxCoM_Table$par[4],
                      MLE_H2_Full_CxCoM_Table$par[5],
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      MLE_H2_Full_CxCoM_Table$par[6],
                      MLE_H2_Full_CxCoM_Table$par[7],
                      MLE_H2_Full_CxCoM_Table$par[8],
                      MLE_H2_Full_CxCoM_Table$par[9],
                      MLE_H2_Full_CxCoM_Table$par[10],
                      MLE_H2_Full_CxCoM_Table$par[11],
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      MLE_H2_Full_CxCoM_Table$par[12],
                      MLE_H2_Full_CxCoM_Table$par[13],
                      MLE_H2_Full_CxCoM_Table$par[14],
                      MLE_H2_Full_CxCoM_Table$par[15],
                      MLE_H2_Full_CxCoM_Table$par[16],
                      MLE_H2_Full_CxCoM_Table$par[17],
                      MLE_H2_Full_CxCoM_Table$value,
                      MLE_H2_Full_CxCoM_Table$convergence)
  
  # Add vector of y, p, likelihood, code to MLE_Uniform_CxCoM
  MLE_H2_full_outputs <- rbind(MLE_H2_full_outputs, MLE_H2_outputs)
}

# Take output matrix MLE_Uniform_CxCoM and sort by likelihood
# Print the sorted table for records, make sure to RENAME FILE
MLE_H2_full<-as.data.frame(MLE_H2_full_outputs)
df2 <- MLE_H2_full[order(MLE_H2_full$lnL),]
write.csv(df2,"../Summaries and Outputs/19_multi_locus_full_experiment_H2_mle_output.csv",row.names = TRUE)

# Begin algorithm for fitting CxCoM H3 to full experiment dataset
# Create matrix to store all optimized outputs of CxCoM.Likelihood.H3
# Log likelihoods calculated for each combination of y & p model parameters
MLE_H3_full_outputs <- matrix(nrow=0, ncol=39)
colnames(MLE_H3_full_outputs)<-c("ID",
                                 "MLE_X_length",
                                 "MLE_x1",
                                 "MLE_x2",
                                 "MLE_x3",
                                 "MLE_x4",
                                 "MLE_x5",
                                 "MLE_p",
                                 "MLE_y1",
                                 "MLE_y2",
                                 "MLE_y3",
                                 "MLE_y4",
                                 "MLE_y5",
                                 "MLE_v1f+",
                                 "MLE_v2f+",
                                 "MLE_v3f+",
                                 "MLE_v4f+",
                                 "MLE_v5f+",
                                 "MLE_v6f+",
                                 "MLE_v1f-",
                                 "MLE_v2f-",
                                 "MLE_v3f-",
                                 "MLE_v4f-",
                                 "MLE_v5f-",
                                 "MLE_v6f-",
                                 "MLE_v1m+",
                                 "MLE_v2m+",
                                 "MLE_v3m+",
                                 "MLE_v4m+",
                                 "MLE_v5m+",
                                 "MLE_v6m+",
                                 "MLE_v1m-",
                                 "MLE_v2m-",
                                 "MLE_v3m-",
                                 "MLE_v4m-",
                                 "MLE_v5m-",
                                 "MLE_v6m-",
                                 "lnL",
                                 "Convergence")


# For loop Uniform.CxCoM.Likelihood for all possible combinations of p
# Output model parameters and likelihoods in matrix MLE_H3_full_outputs
for (v in 1:max_p_value) {
  
  # Define parameters starting points for y, v, and p in single run of Nelder-Mead algorithm
  # When bounding Nelder-Mead using package dfoptim, starting points cannot be at the boundaries
  parameters <- as.numeric(c(y_uniform[v,1],
                             y_uniform[v,2],
                             y_uniform[v,3],
                             y_uniform[v,4],
                             y_uniform[v,5],
                             0.738,
                             0.741,
                             0.739,
                             0.772,
                             0.817,
                             0.766,
                             0.619,
                             0.615,
                             0.617,
                             0.584,
                             0.540,
                             0.591,
                             0.606,
                             0.606,
                             0.619,
                             0.673,
                             0.748,
                             0.657,
                             0.479,
                             0.479,
                             0.466,
                             0.412,
                             0.337,
                             0.428))
  
  p <- p_uniform[v,]
  
  # Optimize CxCoM model parameters y and v with bounded Nelder-Mead method
  # Upper limit on y's is arbitrary, whereas upper and lower limit on v's is
  # given by definition of egg-to-adult viability (cannot exceed 0-1 range)
  MLE_H3_Full_CxCoM_Table <- nmkb(par = parameters,
                                  fn = CxCoM.Likelihood.H3,
                                  control = list(maxfeval = 50000),
                                  lower = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                  upper = c(9,9,9,9,9,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))
  
  # Create output vector of y, p, likelihood, convergence code
  MLE_H3_outputs <- c(3,
                      (((MLE_H3_Full_CxCoM_Table$par[1]+MLE_H3_Full_CxCoM_Table$par[2]+MLE_H3_Full_CxCoM_Table$par[3]+MLE_H3_Full_CxCoM_Table$par[4]+MLE_H3_Full_CxCoM_Table$par[5])/(2*p[1]))*100),
                      MLE_H3_Full_CxCoM_Table$par[1]/(2*p[1]),
                      MLE_H3_Full_CxCoM_Table$par[2]/(2*p[1]),
                      MLE_H3_Full_CxCoM_Table$par[3]/(2*p[1]),
                      MLE_H3_Full_CxCoM_Table$par[4]/(2*p[1]),
                      MLE_H3_Full_CxCoM_Table$par[5]/(2*p[1]),
                      p[1],
                      MLE_H3_Full_CxCoM_Table$par[1],
                      MLE_H3_Full_CxCoM_Table$par[2],
                      MLE_H3_Full_CxCoM_Table$par[3],
                      MLE_H3_Full_CxCoM_Table$par[4],
                      MLE_H3_Full_CxCoM_Table$par[5],
                      MLE_H3_Full_CxCoM_Table$par[6],
                      MLE_H3_Full_CxCoM_Table$par[7],
                      MLE_H3_Full_CxCoM_Table$par[8],
                      MLE_H3_Full_CxCoM_Table$par[9],
                      MLE_H3_Full_CxCoM_Table$par[10],
                      MLE_H3_Full_CxCoM_Table$par[11],
                      MLE_H3_Full_CxCoM_Table$par[12],
                      MLE_H3_Full_CxCoM_Table$par[13],
                      MLE_H3_Full_CxCoM_Table$par[14],
                      MLE_H3_Full_CxCoM_Table$par[15],
                      MLE_H3_Full_CxCoM_Table$par[16],
                      MLE_H3_Full_CxCoM_Table$par[17],
                      MLE_H3_Full_CxCoM_Table$par[18],
                      MLE_H3_Full_CxCoM_Table$par[19],
                      MLE_H3_Full_CxCoM_Table$par[20],
                      MLE_H3_Full_CxCoM_Table$par[21],
                      MLE_H3_Full_CxCoM_Table$par[22],
                      MLE_H3_Full_CxCoM_Table$par[23],
                      MLE_H3_Full_CxCoM_Table$par[24],
                      MLE_H3_Full_CxCoM_Table$par[25],
                      MLE_H3_Full_CxCoM_Table$par[26],
                      MLE_H3_Full_CxCoM_Table$par[27],
                      MLE_H3_Full_CxCoM_Table$par[28],
                      MLE_H3_Full_CxCoM_Table$par[29],
                      MLE_H3_Full_CxCoM_Table$value,
                      MLE_H3_Full_CxCoM_Table$convergence)
  
  # Add vector of y, p, likelihood, code to MLE_Uniform_CxCoM
  MLE_H3_full_outputs <- rbind(MLE_H3_full_outputs, MLE_H3_outputs)
}

# Take output matrix MLE_Uniform_CxCoM and sort by likelihood
# Print the sorted table for records, make sure to RENAME FILE
MLE_H3_full<-as.data.frame(MLE_H3_full_outputs)
df3 <- MLE_H3_full[order(MLE_H3_full$lnL),]
write.csv(df3,"../Summaries and Outputs/20_multi_locus_full_experiment_H3_mle_output.csv",row.names = TRUE)

# RUN H3 ON VIAL, PERFORM REGRESSION

# Begin algorithm for fitting CxCoM H3 to individual vial dataset
# Create matrix to store all optimized outputs of looped CxCoM.Likelihood.H3
# Log likelihoods calculated for each combination of y & p model parameters
MLE_H3_looped_vials_outputs <- matrix(nrow=0, ncol=39)
colnames(MLE_H3_looped_vials_outputs)<-c("ID",
                                 "MLE_X_length",
                                 "MLE_x1",
                                 "MLE_x2",
                                 "MLE_x3",
                                 "MLE_x4",
                                 "MLE_x5",
                                 "MLE_p",
                                 "MLE_y1",
                                 "MLE_y2",
                                 "MLE_y3",
                                 "MLE_y4",
                                 "MLE_y5",
                                 "MLE_v1f+",
                                 "MLE_v2f+",
                                 "MLE_v3f+",
                                 "MLE_v4f+",
                                 "MLE_v5f+",
                                 "MLE_v6f+",
                                 "MLE_v1f-",
                                 "MLE_v2f-",
                                 "MLE_v3f-",
                                 "MLE_v4f-",
                                 "MLE_v5f-",
                                 "MLE_v6f-",
                                 "MLE_v1m+",
                                 "MLE_v2m+",
                                 "MLE_v3m+",
                                 "MLE_v4m+",
                                 "MLE_v5m+",
                                 "MLE_v6m+",
                                 "MLE_v1m-",
                                 "MLE_v2m-",
                                 "MLE_v3m-",
                                 "MLE_v4m-",
                                 "MLE_v5m-",
                                 "MLE_v6m-",
                                 "lnL",
                                 "Convergence")

for (g in 1:70) {
  single_experimental_unit <- as.numeric(multi_locus_individual_vials_dataset[g,])
  observed_count <- c(single_experimental_unit[130],single_experimental_unit[2:129])
  
  

  # Create matrix to store all optimized outputs of Uniform.CxCoM.Likelihood
  # Log likelihoods calculated for each combination of y & p model parameters
  MLE_H3_looped_vials_internal <- matrix(nrow=0, ncol=39)
  colnames(MLE_H3_looped_vials_internal)<-c("ID",
                                 "MLE_X_length",
                                 "MLE_x1",
                                 "MLE_x2",
                                 "MLE_x3",
                                 "MLE_x4",
                                 "MLE_x5",
                                 "MLE_p",
                                 "MLE_y1",
                                 "MLE_y2",
                                 "MLE_y3",
                                 "MLE_y4",
                                 "MLE_y5",
                                 "MLE_v1f+",
                                 "MLE_v2f+",
                                 "MLE_v3f+",
                                 "MLE_v4f+",
                                 "MLE_v5f+",
                                 "MLE_v6f+",
                                 "MLE_v1f-",
                                 "MLE_v2f-",
                                 "MLE_v3f-",
                                 "MLE_v4f-",
                                 "MLE_v5f-",
                                 "MLE_v6f-",
                                 "MLE_v1m+",
                                 "MLE_v2m+",
                                 "MLE_v3m+",
                                 "MLE_v4m+",
                                 "MLE_v5m+",
                                 "MLE_v6m+",
                                 "MLE_v1m-",
                                 "MLE_v2m-",
                                 "MLE_v3m-",
                                 "MLE_v4m-",
                                 "MLE_v5m-",
                                 "MLE_v6m-",
                                 "lnL",
                                 "Convergence")
  
  # For loop Uniform.CxCoM.Likelihood for all possible combinations of p
  # Output model parameters and likelihoods in matrix MLE_H3_full_outputs
  for (v in 1:max_p_value) {
    
    # Define parameters starting points for y, v, and p in single run of Nelder-Mead algorithm
    # When bounding Nelder-Mead using package dfoptim, starting points cannot be at the boundaries
    parameters <- as.numeric(c(y_uniform[v,1],
                               y_uniform[v,2],
                               y_uniform[v,3],
                               y_uniform[v,4],
                               y_uniform[v,5],
                               0.738,
                               0.741,
                               0.739,
                               0.772,
                               0.817,
                               0.766,
                               0.619,
                               0.615,
                               0.617,
                               0.584,
                               0.540,
                               0.591,
                               0.606,
                               0.606,
                               0.619,
                               0.673,
                               0.748,
                               0.657,
                               0.479,
                               0.479,
                               0.466,
                               0.412,
                               0.337,
                               0.428))
    
    p <- p_uniform[v,]
    
    # Optimize CxCoM model parameters y and v with bounded Nelder-Mead method
    # Upper limit on y's is arbitrary, whereas upper and lower limit on v's is
    # given by definition of egg-to-adult viability (cannot exceed 0-1 range)
    MLE_H3_Vial_CxCoM_Table <- nmkb(par = parameters,
                                    fn = CxCoM.Likelihood.H3,
                                    control = list(maxfeval = 50000),
                                    lower = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                    upper = c(9,9,9,9,9,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))
    
    # Create output vector of y, p, likelihood, convergence code
    MLE_Vial_outputs <- c(3,
                        (((MLE_H3_Vial_CxCoM_Table$par[1]+MLE_H3_Vial_CxCoM_Table$par[2]+MLE_H3_Vial_CxCoM_Table$par[3]+MLE_H3_Vial_CxCoM_Table$par[4]+MLE_H3_Vial_CxCoM_Table$par[5])/(2*p[1]))*100),
                        MLE_H3_Vial_CxCoM_Table$par[1]/(2*p[1]),
                        MLE_H3_Vial_CxCoM_Table$par[2]/(2*p[1]),
                        MLE_H3_Vial_CxCoM_Table$par[3]/(2*p[1]),
                        MLE_H3_Vial_CxCoM_Table$par[4]/(2*p[1]),
                        MLE_H3_Vial_CxCoM_Table$par[5]/(2*p[1]),
                        p[1],
                        MLE_H3_Vial_CxCoM_Table$par[1],
                        MLE_H3_Vial_CxCoM_Table$par[2],
                        MLE_H3_Vial_CxCoM_Table$par[3],
                        MLE_H3_Vial_CxCoM_Table$par[4],
                        MLE_H3_Vial_CxCoM_Table$par[5],
                        MLE_H3_Vial_CxCoM_Table$par[6],
                        MLE_H3_Vial_CxCoM_Table$par[7],
                        MLE_H3_Vial_CxCoM_Table$par[8],
                        MLE_H3_Vial_CxCoM_Table$par[9],
                        MLE_H3_Vial_CxCoM_Table$par[10],
                        MLE_H3_Vial_CxCoM_Table$par[11],
                        MLE_H3_Vial_CxCoM_Table$par[12],
                        MLE_H3_Vial_CxCoM_Table$par[13],
                        MLE_H3_Vial_CxCoM_Table$par[14],
                        MLE_H3_Vial_CxCoM_Table$par[15],
                        MLE_H3_Vial_CxCoM_Table$par[16],
                        MLE_H3_Vial_CxCoM_Table$par[17],
                        MLE_H3_Vial_CxCoM_Table$par[18],
                        MLE_H3_Vial_CxCoM_Table$par[19],
                        MLE_H3_Vial_CxCoM_Table$par[20],
                        MLE_H3_Vial_CxCoM_Table$par[21],
                        MLE_H3_Vial_CxCoM_Table$par[22],
                        MLE_H3_Vial_CxCoM_Table$par[23],
                        MLE_H3_Vial_CxCoM_Table$par[24],
                        MLE_H3_Vial_CxCoM_Table$par[25],
                        MLE_H3_Vial_CxCoM_Table$par[26],
                        MLE_H3_Vial_CxCoM_Table$par[27],
                        MLE_H3_Vial_CxCoM_Table$par[28],
                        MLE_H3_Vial_CxCoM_Table$par[29],
                        MLE_H3_Vial_CxCoM_Table$value,
                        MLE_H3_Vial_CxCoM_Table$convergence)
    
    # Add vector of y, p, likelihood, code to MLE_Uniform_CxCoM
    MLE_H3_looped_vials_internal <- rbind(MLE_H3_looped_vials_internal, MLE_Vial_outputs)
  }
  
  # Take output matrix MLE_Uniform_CxCoM and sort by likelihood
  # Print the sorted table for records, make sure to RENAME FILE
  MLE_H3_vials <- as.data.frame(MLE_H3_looped_vials_internal)
  df4 <- MLE_H3_vials[order(MLE_H3_vials$lnL),]
  MLE_H3_looped_vials_outputs<- rbind(MLE_H3_looped_vials_outputs, df4[1,])
  
}

MLE_H3_looped_vials_outputs

write.csv(MLE_H3_looped_vials_outputs,"../Summaries and Outputs/21_multi_locus_individual_vials_H3_mle_output.csv",row.names = TRUE)

# Perform ANOVA on individual vial estimates fitting an intercept only model
vials_pooled_anova_input <- as.data.frame(MLE_H3_looped_vials_outputs)
vials_pooled_anova.lm <- lm(vials_pooled_anova_input$MLE_X_length ~ 1, data=vials_pooled_anova_input)
vials_pooled_anova_anova.table <- anova(vials_pooled_anova.lm)
write.csv(vials_pooled_anova_anova.table,"../Summaries and Outputs/table_S10_vials_pooled_anova.csv",row.names = TRUE)



# Begin algorithm for fitting CxCoM H3 to individual vial dataset
# Create matrix to store all optimized outputs of looped CxCoM.Likelihood.H3
# Log likelihoods calculated for each combination of y & p model parameters
MLE_H3_looped_cross_outputs <- matrix(nrow=0, ncol=39)
colnames(MLE_H3_looped_cross_outputs)<-c("ID",
                                         "MLE_X_length",
                                         "MLE_x1",
                                         "MLE_x2",
                                         "MLE_x3",
                                         "MLE_x4",
                                         "MLE_x5",
                                         "MLE_p",
                                         "MLE_y1",
                                         "MLE_y2",
                                         "MLE_y3",
                                         "MLE_y4",
                                         "MLE_y5",
                                         "MLE_v1f+",
                                         "MLE_v2f+",
                                         "MLE_v3f+",
                                         "MLE_v4f+",
                                         "MLE_v5f+",
                                         "MLE_v6f+",
                                         "MLE_v1f-",
                                         "MLE_v2f-",
                                         "MLE_v3f-",
                                         "MLE_v4f-",
                                         "MLE_v5f-",
                                         "MLE_v6f-",
                                         "MLE_v1m+",
                                         "MLE_v2m+",
                                         "MLE_v3m+",
                                         "MLE_v4m+",
                                         "MLE_v5m+",
                                         "MLE_v6m+",
                                         "MLE_v1m-",
                                         "MLE_v2m-",
                                         "MLE_v3m-",
                                         "MLE_v4m-",
                                         "MLE_v5m-",
                                         "MLE_v6m-",
                                         "lnL",
                                         "Convergence")

for (g in 1:10) {
  single_experimental_unit <- as.numeric(multi_locus_cross_pooled_dataset[g,])
  observed_count <- c(single_experimental_unit[130],single_experimental_unit[2:129])
  
  
  
  # Create matrix to store all optimized outputs of Uniform.CxCoM.Likelihood
  # Log likelihoods calculated for each combination of y & p model parameters
  MLE_H3_looped_cross_internal <- matrix(nrow=0, ncol=39)
  colnames(MLE_H3_looped_cross_internal)<-c("ID",
                                            "MLE_X_length",
                                            "MLE_x1",
                                            "MLE_x2",
                                            "MLE_x3",
                                            "MLE_x4",
                                            "MLE_x5",
                                            "MLE_p",
                                            "MLE_y1",
                                            "MLE_y2",
                                            "MLE_y3",
                                            "MLE_y4",
                                            "MLE_y5",
                                            "MLE_v1f+",
                                            "MLE_v2f+",
                                            "MLE_v3f+",
                                            "MLE_v4f+",
                                            "MLE_v5f+",
                                            "MLE_v6f+",
                                            "MLE_v1f-",
                                            "MLE_v2f-",
                                            "MLE_v3f-",
                                            "MLE_v4f-",
                                            "MLE_v5f-",
                                            "MLE_v6f-",
                                            "MLE_v1m+",
                                            "MLE_v2m+",
                                            "MLE_v3m+",
                                            "MLE_v4m+",
                                            "MLE_v5m+",
                                            "MLE_v6m+",
                                            "MLE_v1m-",
                                            "MLE_v2m-",
                                            "MLE_v3m-",
                                            "MLE_v4m-",
                                            "MLE_v5m-",
                                            "MLE_v6m-",
                                            "lnL",
                                            "Convergence")
  
  # For loop Uniform.CxCoM.Likelihood for all possible combinations of p
  # Output model parameters and likelihoods in matrix MLE_H3_full_outputs
  for (v in 1:max_p_value) {
    
    # Define parameters starting points for y, v, and p in single run of Nelder-Mead algorithm
    # When bounding Nelder-Mead using package dfoptim, starting points cannot be at the boundaries
    parameters <- as.numeric(c(y_uniform[v,1],
                               y_uniform[v,2],
                               y_uniform[v,3],
                               y_uniform[v,4],
                               y_uniform[v,5],
                               0.738,
                               0.741,
                               0.739,
                               0.772,
                               0.817,
                               0.766,
                               0.619,
                               0.615,
                               0.617,
                               0.584,
                               0.540,
                               0.591,
                               0.606,
                               0.606,
                               0.619,
                               0.673,
                               0.748,
                               0.657,
                               0.479,
                               0.479,
                               0.466,
                               0.412,
                               0.337,
                               0.428))
    
    p <- p_uniform[v,]
    
    # Optimize CxCoM model parameters y and v with bounded Nelder-Mead method
    # Upper limit on y's is arbitrary, whereas upper and lower limit on v's is
    # given by definition of egg-to-adult viability (cannot exceed 0-1 range)
    MLE_H3_Cross_CxCoM_Table <- nmkb(par = parameters,
                                    fn = CxCoM.Likelihood.H3,
                                    control = list(maxfeval = 50000),
                                    lower = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                    upper = c(9,9,9,9,9,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))
    
    # Create output vector of y, p, likelihood, convergence code
    MLE_Cross_outputs <- c(3,
                          (((MLE_H3_Cross_CxCoM_Table$par[1]+MLE_H3_Cross_CxCoM_Table$par[2]+MLE_H3_Cross_CxCoM_Table$par[3]+MLE_H3_Cross_CxCoM_Table$par[4]+MLE_H3_Cross_CxCoM_Table$par[5])/(2*p[1]))*100),
                          MLE_H3_Cross_CxCoM_Table$par[1]/(2*p[1]),
                          MLE_H3_Cross_CxCoM_Table$par[2]/(2*p[1]),
                          MLE_H3_Cross_CxCoM_Table$par[3]/(2*p[1]),
                          MLE_H3_Cross_CxCoM_Table$par[4]/(2*p[1]),
                          MLE_H3_Cross_CxCoM_Table$par[5]/(2*p[1]),
                          p[1],
                          MLE_H3_Cross_CxCoM_Table$par[1],
                          MLE_H3_Cross_CxCoM_Table$par[2],
                          MLE_H3_Cross_CxCoM_Table$par[3],
                          MLE_H3_Cross_CxCoM_Table$par[4],
                          MLE_H3_Cross_CxCoM_Table$par[5],
                          MLE_H3_Cross_CxCoM_Table$par[6],
                          MLE_H3_Cross_CxCoM_Table$par[7],
                          MLE_H3_Cross_CxCoM_Table$par[8],
                          MLE_H3_Cross_CxCoM_Table$par[9],
                          MLE_H3_Cross_CxCoM_Table$par[10],
                          MLE_H3_Cross_CxCoM_Table$par[11],
                          MLE_H3_Cross_CxCoM_Table$par[12],
                          MLE_H3_Cross_CxCoM_Table$par[13],
                          MLE_H3_Cross_CxCoM_Table$par[14],
                          MLE_H3_Cross_CxCoM_Table$par[15],
                          MLE_H3_Cross_CxCoM_Table$par[16],
                          MLE_H3_Cross_CxCoM_Table$par[17],
                          MLE_H3_Cross_CxCoM_Table$par[18],
                          MLE_H3_Cross_CxCoM_Table$par[19],
                          MLE_H3_Cross_CxCoM_Table$par[20],
                          MLE_H3_Cross_CxCoM_Table$par[21],
                          MLE_H3_Cross_CxCoM_Table$par[22],
                          MLE_H3_Cross_CxCoM_Table$par[23],
                          MLE_H3_Cross_CxCoM_Table$par[24],
                          MLE_H3_Cross_CxCoM_Table$par[25],
                          MLE_H3_Cross_CxCoM_Table$par[26],
                          MLE_H3_Cross_CxCoM_Table$par[27],
                          MLE_H3_Cross_CxCoM_Table$par[28],
                          MLE_H3_Cross_CxCoM_Table$par[29],
                          MLE_H3_Cross_CxCoM_Table$value,
                          MLE_H3_Cross_CxCoM_Table$convergence)
    
    # Add vector of y, p, likelihood, code to MLE_Uniform_CxCoM
    MLE_H3_looped_cross_internal <- rbind(MLE_H3_looped_cross_internal, MLE_Cross_outputs)
  }
  
  # Take output matrix MLE_Uniform_CxCoM and sort by likelihood
  # Print the sorted table for records, make sure to RENAME FILE
  MLE_H3_cross <- as.data.frame(MLE_H3_looped_cross_internal)
  df5 <- MLE_H3_cross[order(MLE_H3_cross$lnL),]
  MLE_H3_looped_cross_outputs<- rbind(MLE_H3_looped_cross_outputs, df5[1,])
  
}

MLE_H3_looped_cross_outputs
write.csv(MLE_H3_looped_cross_outputs,"../Summaries and Outputs/22_multi_locus_cross_pooled_H3_mle_output.csv",row.names = TRUE)

# Perform ANOVA on cross pooled estimates fitting an intercept only model
cross_pooled_anova_input <- as.data.frame(MLE_H3_looped_cross_outputs)
cross_pooled_anova.lm <- lm(cross_pooled_anova_input$MLE_X_length ~ 1, data=cross_pooled_anova_input)
cross_pooled_anova_anova.table <- anova(cross_pooled_anova.lm)
write.csv(cross_pooled_anova_anova.table,"../Summaries and Outputs/table_S11_cross_pooled_anova.csv",row.names = TRUE)


# RUN H3 on Brood
# Begin algorithm for fitting CxCoM H3 to individual vial dataset
# Create matrix to store all optimized outputs of looped CxCoM.Likelihood.H3
# Log likelihoods calculated for each combination of y & p model parameters
MLE_H3_looped_brood_outputs <- matrix(nrow=0, ncol=39)
colnames(MLE_H3_looped_brood_outputs)<-c("ID",
                                         "MLE_X_length",
                                         "MLE_x1",
                                         "MLE_x2",
                                         "MLE_x3",
                                         "MLE_x4",
                                         "MLE_x5",
                                         "MLE_p",
                                         "MLE_y1",
                                         "MLE_y2",
                                         "MLE_y3",
                                         "MLE_y4",
                                         "MLE_y5",
                                         "MLE_v1f+",
                                         "MLE_v2f+",
                                         "MLE_v3f+",
                                         "MLE_v4f+",
                                         "MLE_v5f+",
                                         "MLE_v6f+",
                                         "MLE_v1f-",
                                         "MLE_v2f-",
                                         "MLE_v3f-",
                                         "MLE_v4f-",
                                         "MLE_v5f-",
                                         "MLE_v6f-",
                                         "MLE_v1m+",
                                         "MLE_v2m+",
                                         "MLE_v3m+",
                                         "MLE_v4m+",
                                         "MLE_v5m+",
                                         "MLE_v6m+",
                                         "MLE_v1m-",
                                         "MLE_v2m-",
                                         "MLE_v3m-",
                                         "MLE_v4m-",
                                         "MLE_v5m-",
                                         "MLE_v6m-",
                                         "lnL",
                                         "Convergence")

for (g in 1:7) {
  single_experimental_unit <- as.numeric(multi_locus_brood_pooled_dataset[g,])
  observed_count <- c(single_experimental_unit[130],single_experimental_unit[2:129])
  
  
  
  # Create matrix to store all optimized outputs of Uniform.CxCoM.Likelihood
  # Log likelihoods calculated for each combination of y & p model parameters
  MLE_H3_looped_brood_internal <- matrix(nrow=0, ncol=39)
  colnames(MLE_H3_looped_brood_internal)<-c("ID",
                                            "MLE_X_length",
                                            "MLE_x1",
                                            "MLE_x2",
                                            "MLE_x3",
                                            "MLE_x4",
                                            "MLE_x5",
                                            "MLE_p",
                                            "MLE_y1",
                                            "MLE_y2",
                                            "MLE_y3",
                                            "MLE_y4",
                                            "MLE_y5",
                                            "MLE_v1f+",
                                            "MLE_v2f+",
                                            "MLE_v3f+",
                                            "MLE_v4f+",
                                            "MLE_v5f+",
                                            "MLE_v6f+",
                                            "MLE_v1f-",
                                            "MLE_v2f-",
                                            "MLE_v3f-",
                                            "MLE_v4f-",
                                            "MLE_v5f-",
                                            "MLE_v6f-",
                                            "MLE_v1m+",
                                            "MLE_v2m+",
                                            "MLE_v3m+",
                                            "MLE_v4m+",
                                            "MLE_v5m+",
                                            "MLE_v6m+",
                                            "MLE_v1m-",
                                            "MLE_v2m-",
                                            "MLE_v3m-",
                                            "MLE_v4m-",
                                            "MLE_v5m-",
                                            "MLE_v6m-",
                                            "lnL",
                                            "Convergence")
  
  # For loop Uniform.CxCoM.Likelihood for all possible combinations of p
  # Output model parameters and likelihoods in matrix MLE_H3_full_outputs
  for (v in 1:max_p_value) {
    
    # Define parameters starting points for y, v, and p in single run of Nelder-Mead algorithm
    # When bounding Nelder-Mead using package dfoptim, starting points cannot be at the boundaries
    parameters <- as.numeric(c(y_uniform[v,1],
                               y_uniform[v,2],
                               y_uniform[v,3],
                               y_uniform[v,4],
                               y_uniform[v,5],
                               0.738,
                               0.741,
                               0.739,
                               0.772,
                               0.817,
                               0.766,
                               0.619,
                               0.615,
                               0.617,
                               0.584,
                               0.540,
                               0.591,
                               0.606,
                               0.606,
                               0.619,
                               0.673,
                               0.748,
                               0.657,
                               0.479,
                               0.479,
                               0.466,
                               0.412,
                               0.337,
                               0.428))
    
    p <- p_uniform[v,]
    
    # Optimize CxCoM model parameters y and v with bounded Nelder-Mead method
    # Upper limit on y's is arbitrary, whereas upper and lower limit on v's is
    # given by definition of egg-to-adult viability (cannot exceed 0-1 range)
    MLE_H3_Brood_CxCoM_Table <- nmkb(par = parameters,
                                    fn = CxCoM.Likelihood.H3,
                                    control = list(maxfeval = 50000),
                                    lower = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                    upper = c(9,9,9,9,9,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))
    
    # Create output vector of y, p, likelihood, convergence code
    MLE_Brood_outputs <- c(3,
                          (((MLE_H3_Brood_CxCoM_Table$par[1]+MLE_H3_Brood_CxCoM_Table$par[2]+MLE_H3_Brood_CxCoM_Table$par[3]+MLE_H3_Brood_CxCoM_Table$par[4]+MLE_H3_Brood_CxCoM_Table$par[5])/(2*p[1]))*100),
                          MLE_H3_Brood_CxCoM_Table$par[1]/(2*p[1]),
                          MLE_H3_Brood_CxCoM_Table$par[2]/(2*p[1]),
                          MLE_H3_Brood_CxCoM_Table$par[3]/(2*p[1]),
                          MLE_H3_Brood_CxCoM_Table$par[4]/(2*p[1]),
                          MLE_H3_Brood_CxCoM_Table$par[5]/(2*p[1]),
                          p[1],
                          MLE_H3_Brood_CxCoM_Table$par[1],
                          MLE_H3_Brood_CxCoM_Table$par[2],
                          MLE_H3_Brood_CxCoM_Table$par[3],
                          MLE_H3_Brood_CxCoM_Table$par[4],
                          MLE_H3_Brood_CxCoM_Table$par[5],
                          MLE_H3_Brood_CxCoM_Table$par[6],
                          MLE_H3_Brood_CxCoM_Table$par[7],
                          MLE_H3_Brood_CxCoM_Table$par[8],
                          MLE_H3_Brood_CxCoM_Table$par[9],
                          MLE_H3_Brood_CxCoM_Table$par[10],
                          MLE_H3_Brood_CxCoM_Table$par[11],
                          MLE_H3_Brood_CxCoM_Table$par[12],
                          MLE_H3_Brood_CxCoM_Table$par[13],
                          MLE_H3_Brood_CxCoM_Table$par[14],
                          MLE_H3_Brood_CxCoM_Table$par[15],
                          MLE_H3_Brood_CxCoM_Table$par[16],
                          MLE_H3_Brood_CxCoM_Table$par[17],
                          MLE_H3_Brood_CxCoM_Table$par[18],
                          MLE_H3_Brood_CxCoM_Table$par[19],
                          MLE_H3_Brood_CxCoM_Table$par[20],
                          MLE_H3_Brood_CxCoM_Table$par[21],
                          MLE_H3_Brood_CxCoM_Table$par[22],
                          MLE_H3_Brood_CxCoM_Table$par[23],
                          MLE_H3_Brood_CxCoM_Table$par[24],
                          MLE_H3_Brood_CxCoM_Table$par[25],
                          MLE_H3_Brood_CxCoM_Table$par[26],
                          MLE_H3_Brood_CxCoM_Table$par[27],
                          MLE_H3_Brood_CxCoM_Table$par[28],
                          MLE_H3_Brood_CxCoM_Table$par[29],
                          MLE_H3_Brood_CxCoM_Table$value,
                          MLE_H3_Brood_CxCoM_Table$convergence)
    
    # Add vector of y, p, likelihood, code to MLE_Uniform_CxCoM
    MLE_H3_looped_brood_internal <- rbind(MLE_H3_looped_brood_internal, MLE_Brood_outputs)
  }
  
  # Take output matrix MLE_Uniform_CxCoM and sort by likelihood
  # Print the sorted table for records, make sure to RENAME FILE
  MLE_H3_brood <- as.data.frame(MLE_H3_looped_brood_internal)
  df6 <- MLE_H3_brood[order(MLE_H3_brood$lnL),]
  MLE_H3_looped_brood_outputs<- rbind(MLE_H3_looped_brood_outputs, df6[1,])
  
}

MLE_H3_looped_brood_outputs
write.csv(MLE_H3_looped_brood_outputs,"../Summaries and Outputs/23_multi_locus_brood_pooled_H3_mle_output.csv",row.names = TRUE)

# Perform ANOVA on brood pooled estimates fitting an intercept only model
brood_pooled_anova_input <- as.data.frame(MLE_H3_looped_brood_outputs)
brood_pooled_anova.lm <- lm(brood_pooled_anova_input$MLE_X_length ~ 1, data=brood_pooled_anova_input)
brood_pooled_anova_anova.table <- anova(brood_pooled_anova.lm)
write.csv(brood_pooled_anova_anova.table,"../Summaries and Outputs/table_S12_brood_pooled_anova.csv",row.names = TRUE)
