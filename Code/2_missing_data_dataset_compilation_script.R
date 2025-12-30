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

