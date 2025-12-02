# Set working directory
setwd("/Users/spencerkoury/Desktop/local_protocols/missing_data_analyses")

# Load package dfoptim
library(dfoptim)

# Read in raw counts from csv file for a single recombination experiment. 
# Data structure of row x column is experimental unit x phenotypic class.
# Each cell represents an observed, sex-specific, phenotypic class count.
# Phenotypic classes are coded as 1,0 by presence of six visible markers.
# The order of phenotypic classes is pre-defined (see example dataset).
# One column added for experimental unit-specific, marker-present egg counts.
# Three additional columns added for the total/female/male treatment-specific,
# independently-estimated, marker-free, egg-to-adult development failure rates.
full_experiment_dataset <- read.csv("SK14_full_experiment.csv")

# Extract raw count data for single experimental unit
single_experimental_unit <- as.numeric(full_experiment_dataset) []#pull out a single row, and this is the place to iterate

# Calculate total and sex-specific sums of adult counts
adult_total_raw_count <- sum(single_experimental_unit[1:128])
adult_female_raw_count <- sum(single_experimental_unit[1:64])
adult_male_raw_count <- sum(single_experimental_unit[65:128])

# Calculate sex-specific egg derived counts assuming X and Y sperm have 
# equal chance of fertilizing eggs, sex-specific counts must be integers
egg_total_raw_count <- single_experimental_unit[129]
egg_female_derived_count <- round(egg_total_raw_count/2, 0)
egg_male_derived_count <- round(egg_total_raw_count/2, 0)

# Calculate the total and sex-specific derived counts of observed eggs 
# failing to develop into adults scorable for evidence of recombination
failed_total_derived_count <- egg_total_raw_count - adult_total_raw_count
failed_female_derived_count <- egg_female_derived_count - adult_female_raw_count
failed_male_derived_count <- egg_male_derived_count - adult_male_raw_count

# For Hypothesis H1 only, frequency of failed development when markers are present
failed_total_frequency <- 1-(adult_total_raw_count/egg_total_raw_count)
failed_female_frequency <- 1-(adult_female_raw_count/egg_female_derived_count)
failed_male_frequency <- 1-(adult_male_raw_count/egg_male_derived_count)

# Define the observed count vector for multinomial likelihood calculations
# Structure is: count of female marker-present failed egg-to-adult development
#.              count of female phenotypic classes 1-64 (in predefined order)
#               count of male marker-present failed egg-to-adult development
#               count of male phenotypic classes 65-128 (in predefined order)
# optional pooling of high order exchange classes (see rec in Zhao et al 1995)
# obs_high_order_exchange_f <- sum(single_experimental_unit[33:64])
# obs_high_order_exchange_m <- sum(single_experimental_unit[97:128])
observed_count <- c(failed_female_derived_count,
                    single_experimental_unit[1:64],
                    failed_male_derived_count,
                    single_experimental_unit[65:128])

# Define number of marker based intervals and max p&k value in an interval
# Parameter p gives strength of renewal process (p=m+1) in the CxCoM model
# Parameter k gives the number of crossovers in an interval in CxCoM model
no_interval <- 5
max_p_value <- 10
max_k_value <- 3

# Enter estimate of genetic length (x) in each interval (index) in Morgans
# Use to calculate initial value of model parameter y optimized by optim()
# Can use published genetic maps or recombination fractions from raw data
x1 <- 0.03
x2 <- 0.06
x3 <- 0.24
x4 <- 0.21
x5 <- 0.12
genetic_lengths <- c(x1,x2,x3,x4,x5)

# Generate combinations of the model parameter p in marker-based intervals 
# Parameter p gives strength of renewal process (p=m+1) in the CxCoM model
# Evaluate values of p from "no interference" (p=1) up to (p=max_p_value)
for (l in 1:no_interval) {
  assign(paste0("p", l), c(1:max_p_value))
}
p_uniform <- cbind(p1,p2,p3,p4,p5)
colnames(p_uniform) <- c("p1","p2","p3","p4","p5")

# Generate starting point for parameter y in each interval to be optimized
# Calculate as 2*p*x using the p_uniform matrix and genetic_length vector
y_uniform <- as.matrix(cbind(2*p_uniform[,1]*x1,
                             2*p_uniform[,2]*x2,
                             2*p_uniform[,3]*x3,
                             2*p_uniform[,4]*x4,
                             2*p_uniform[,5]*x5))
colnames(y_uniform) <- c("y1","y2","y3","y4","y5")

# Define the probability of an egg being fertilized by an X or Y bearing sperm
Pr_X_bearing_sperm <- 0.5
Pr_Y_bearing_sperm <- 0.5

# Define marker viabilities external to the custom function (i.e. constrained model)
# These default values are superceded by fitting viabilities inside custom function

v1_wildtype_f <- 1.00
v2_wildtype_f <- 1.00
v3_wildtype_f <- 1.00
v4_wildtype_f <- 1.00
v5_wildtype_f <- 1.00
v6_wildtype_f <- 1.00

v1__mutant__f <- 1.00
v2__mutant__f <- 1.00
v3__mutant__f <- 1.00
v4__mutant__f <- 1.00
v5__mutant__f <- 1.00
v6__mutant__f <- 1.00

v1_wildtype_m <- 1.00
v2_wildtype_m <- 1.00
v3_wildtype_m <- 1.00
v4_wildtype_m <- 1.00
v5_wildtype_m <- 1.00
v6_wildtype_m <- 1.00

v1__mutant__m <- 1.00
v2__mutant__m <- 1.00
v3__mutant__m <- 1.00
v4__mutant__m <- 1.00
v5__mutant__m <- 1.00
v6__mutant__m <- 1.00

# Custom function to calculate the likelihood of specific Cx(Co)M model
# Corresponds to the analysis of chi-square model from Zhao et al. 1995
# Required inputs are "parameters" a vector of y values to be optimized
# Probabilities for classes will need to be tailored to each experiment
Uniform.CxCoM.Likelihood <- function(parameters) {
  
  # Vectors for internal use in the custom function
  # Dimensions auto-adjusted for combinations of p
  y <- parameters[1:5]
  v1__mutant__f <- parameters[6]
  v2__mutant__f <- parameters[7]
  v3__mutant__f <- parameters[8]
  v4__mutant__f <- parameters[9]
  v5__mutant__f <- parameters[10]
  v6__mutant__f <- parameters[11]
  v1__mutant__m <- parameters[12]
  v2__mutant__m <- parameters[13]
  v3__mutant__m <- parameters[14]
  v4__mutant__m <- parameters[15]
  v5__mutant__m <- parameters[16]
  v6__mutant__m <- parameters[17]
  stationary_vector <- matrix(data=1, nrow=1, ncol=p[1])
  summation_vector <- matrix(data=1, nrow=p[no_interval], ncol=1)
  
  # Matrices for internal use in custom function, stored as lists
  # Dimensions auto-adjusted to allow valid matrix multiplication
  kclass <- list(matrix(data=1, nrow=ncol(stationary_vector), ncol=p[1]))
  for (l in 2:no_interval) {
    kclass[[l]] <- matrix(data=1, nrow=p[l-1], ncol=p[l])
  }
  
  k0class <- list(matrix(data=1, nrow=ncol(stationary_vector), ncol=p[1]))
  for (l in 2:no_interval) {
    k0class[[l]] <- matrix(data=1, nrow=p[l-1], ncol=p[l])
  }
  
  # CASE OF k > 0
  #-------------------------------------------------------------------------
  # For loops giving probability of "k" number of crossovers in interval "l"
  #   1) assuming "s" is the number of "C" precursor events in interval "l" 
  #      and events are produced by poisson process with rate parameter "y"
  #   2) assuming "C" precursor events mature to "k" number of "Cx" events
  #      with "s-k" number of "Co" events in the counting process "Cx(Co)^m"
  #   3) assuming the "Cx(Co)^m" counting process is stationary in that the
  #      first "C" equally likely to mature into "Cx" or any one of the "Co"
  #   4) defining "i" as the identity of the first "C" in interval "l" that
  #      numerically is position in series CxCoCoCo plus 1 (i for Cx = 1)
  #   5) defining "j" as the identity of the last "C" in interval "l" that
  #      numerically is p minus position in series CoCoCoCx (j for Cx = 0)
  #
  # R script uses pk-p+i+j definition from Zhao et al. 1995 appendix proof
  # Dk_yl is p(l-1) by p(l) matrix of crossover probabilities in interval l
  # Summing over all k>0 (max_k_value=3) and assuming s>0 in the interval l
  # Probability of k=0 treated separately
  #
  # Matrix multiplication using interval and parameter notation requires 
  # First interval treated separately for valid matrix and multiplication
  # First interval is square matrix by virtue of stationary assumption with
  # The implied definition of p[l-1] as the number of columns in 0th interval
  
  # Nested for loops for the first interval
  interval_matrix <- as.matrix(kclass[[1]])
  for (k in 1:max_k_value) {
    for (i in 1:p[1]) {
      for (j in (p[1]-1):0) {
        interval_matrix[i,(p[1]-j)] <- ((exp(-y[1]))*(y[1]^((p[1]*k)-p[1]+i+j)))/(factorial((p[1]*k)-p[1]+i+j))
      } 
    }
    assign(paste0("Dk", k), interval_matrix)
  }
  assign(paste0("Dk_y", 1), Dk1+Dk2+Dk3)
  
  # Nested for loops for all other intervals
  
  for (l in 2:no_interval) {
    interval_matrix <- kclass[[l]]
    for (k in 1:max_k_value) {
      for (i in 1:p[l-1]) {
        for (j in (p[l]-1):0) {
          interval_matrix[i,(p[l]-j)] <- ((exp(-y[l]))*(y[l]^((p[l]*k)-p[l]+i+j)))/(factorial((p[l]*k)-p[l]+i+j))
        } 
      }
      assign(paste0("Dk", k), interval_matrix)
    }
    assign(paste0("Dk_y", l), Dk1+Dk2+Dk3)
  }
  
  # R_interval is a p(l-1) by p(l) matrix describing the probability of
  # observing a crossover in interval l with model parameters yl and pl
  R_1 <- (1/2)*Dk_y1
  R_2 <- (1/2)*Dk_y2
  R_3 <- (1/2)*Dk_y3
  R_4 <- (1/2)*Dk_y4
  R_5 <- (1/2)*Dk_y5
  
  # CASE OF k = 0
  #-------------------------------------------------------------------------
  # For loops giving probability of "k" number of crossovers in interval "l"
  #   1) assuming "s" is the number of "C" precursor events in interval "l" 
  #      and events are produced by poisson process with rate parameter "y"
  #   2) assuming "C" precursor events mature to "k" number of "Cx" events
  #      with "s-k" number of "Co" events in the counting process "Cx(Co)^m"
  #   3) assuming the "Cx(Co)^m" counting process is stationary in that the
  #      first "C" equally likely to mature into "Cx" or any one of the "Co"
  #   4) defining "i" as the identity of the first "C" in interval "l" that
  #      numerically is position in series CxCoCoCo plus 1 (i for Cx = 1)
  #   5) defining "j" as the identity of the last "C" in interval "l" that
  #      numerically is p minus position in series CoCoCoCx (j for Cx = 0)
  #
  # R script uses pk-p+i+j definition from Zhao et al. 1995 appendix proof
  # Lower triangular matrix D0 in Zhao et al. 1995 was not explicitly defined
  # R script uses implied definition in Zhao et al. 1995 appendix theorem 1
  # D0_yl is p(l-1) by p(l) matrix of probability no crossover in interval l
  #
  # Matrix multiplication using interval and parameter notation requires 
  # First interval treated separately for valid matrix and multiplication
  # First interval is square matrix by virtue of stationary assumption with
  # The implied definition of p[l-1] as the number of columns in 0th interval
  
  # Nested for loops for the first interval
  interval_matrix <- as.matrix(k0class[[1]])
  for (i in 1:p[1]) {
    for (j in (p[1]-1):0) {
      if ((0-p[1]+i+j)>-1) {
        interval_matrix[i,(p[1]-j)] <- ((exp(-y[1]))*(y[1]^(0-p[1]+i+j)))/(factorial(0-p[1]+i+j))
      } 
      else {
        interval_matrix[i,(p[1]-j)]<-0 
      } 
    } 
  }
  assign(paste0("D0_y", 1), interval_matrix)
  
  # Nested for loops for all other intervals
  for (l in 2:no_interval) {
    interval_matrix <- as.matrix(k0class[[l]])
    for (i in 1:p[l-1]) {
      for (j in (p[l]-1):0) {
        if ((0-p[l]+i+j)>-1) {
          interval_matrix[i,p[l]-j] <- ((exp(-y[l]))*(y[l]^(0-p[l]+i+j)))/(factorial(0-p[l]+i+j))
        } 
        else {
          interval_matrix[i,(p[l]-j)]<-0 
        } 
      } 
    }
    assign(paste0("D0_y", l), interval_matrix)
  }
  
  # N_interval is a p(l-1) by p(l) matrix describing the probability of
  # Not observing a crossover in interval l with the parameters yl and pl
  N_1 <- D0_y1+((1/2)*Dk_y1)
  N_2 <- D0_y2+((1/2)*Dk_y2)
  N_3 <- D0_y3+((1/2)*Dk_y3)
  N_4 <- D0_y4+((1/2)*Dk_y4)
  N_5 <- D0_y5+((1/2)*Dk_y5)
  
  # System of equations defining the probabilities of marker patterns generated by the specified crossovers
  # Probability of oogenesis generating an haploid gamete with the marker pattern from a noncrossover event
  
  Pr_000000_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%N_1%*%N_2%*%N_3%*%N_4%*%N_5%*%summation_vector)
  Pr_111111_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%N_1%*%N_2%*%N_3%*%N_4%*%N_5%*%summation_vector)
  
  # Probability of oogenesis generating an haploid gamete with marker pattern from a single crossover event
  
  Pr_100000_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%R_1%*%N_2%*%N_3%*%N_4%*%N_5%*%summation_vector)
  Pr_011111_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%R_1%*%N_2%*%N_3%*%N_4%*%N_5%*%summation_vector)
  
  Pr_110000_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%N_1%*%R_2%*%N_3%*%N_4%*%N_5%*%summation_vector)
  Pr_001111_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%N_1%*%R_2%*%N_3%*%N_4%*%N_5%*%summation_vector)
  
  Pr_111000_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%N_1%*%N_2%*%R_3%*%N_4%*%N_5%*%summation_vector)
  Pr_000111_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%N_1%*%N_2%*%R_3%*%N_4%*%N_5%*%summation_vector)
  
  Pr_111100_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%N_1%*%N_2%*%N_3%*%R_4%*%N_5%*%summation_vector)
  Pr_000011_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%N_1%*%N_2%*%N_3%*%R_4%*%N_5%*%summation_vector)
  
  Pr_111110_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%N_1%*%N_2%*%N_3%*%N_4%*%R_5%*%summation_vector)
  Pr_000001_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%N_1%*%N_2%*%N_3%*%N_4%*%R_5%*%summation_vector)
  
  # Probability of oogenesis generating an haploid gamete with marker pattern from a double crossover event
  
  Pr_010000_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%R_1%*%R_2%*%N_3%*%N_4%*%N_5%*%summation_vector)
  Pr_101111_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%R_1%*%R_2%*%N_3%*%N_4%*%N_5%*%summation_vector)
  
  Pr_011000_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%R_1%*%N_2%*%R_3%*%N_4%*%N_5%*%summation_vector)
  Pr_100111_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%R_1%*%N_2%*%R_3%*%N_4%*%N_5%*%summation_vector)
  
  Pr_011100_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%R_1%*%N_2%*%N_3%*%R_4%*%N_5%*%summation_vector)
  Pr_100011_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%R_1%*%N_2%*%N_3%*%R_4%*%N_5%*%summation_vector)
  
  Pr_011110_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%R_1%*%N_2%*%N_3%*%N_4%*%R_5%*%summation_vector)
  Pr_100001_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%R_1%*%N_2%*%N_3%*%N_4%*%R_5%*%summation_vector)
  
  Pr_001000_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%N_1%*%R_2%*%R_3%*%N_4%*%N_5%*%summation_vector)
  Pr_110111_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%N_1%*%R_2%*%R_3%*%N_4%*%N_5%*%summation_vector)
  
  Pr_001100_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%N_1%*%R_2%*%N_3%*%R_4%*%N_5%*%summation_vector)
  Pr_110011_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%N_1%*%R_2%*%N_3%*%R_4%*%N_5%*%summation_vector)
  
  Pr_001110_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%N_1%*%R_2%*%N_3%*%N_4%*%R_5%*%summation_vector)
  Pr_110001_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%N_1%*%R_2%*%N_3%*%N_4%*%R_5%*%summation_vector)
  
  Pr_000100_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%N_1%*%N_2%*%R_3%*%R_4%*%N_5%*%summation_vector)
  Pr_111011_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%N_1%*%N_2%*%R_3%*%R_4%*%N_5%*%summation_vector)
  
  Pr_000110_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%N_1%*%N_2%*%R_3%*%N_4%*%R_5%*%summation_vector)
  Pr_111001_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%N_1%*%N_2%*%R_3%*%N_4%*%R_5%*%summation_vector)
  
  Pr_000010_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%N_1%*%N_2%*%N_3%*%R_4%*%R_5%*%summation_vector)
  Pr_111101_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%N_1%*%N_2%*%N_3%*%R_4%*%R_5%*%summation_vector)
  
  # Probability of oogenesis generating an haploid gamete with marker pattern from a triple crossover event
  
  Pr_101000_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%R_1%*%R_2%*%R_3%*%N_4%*%N_5%*%summation_vector)
  Pr_010111_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%R_1%*%R_2%*%R_3%*%N_4%*%N_5%*%summation_vector)
  
  Pr_101100_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%R_1%*%R_2%*%N_3%*%R_4%*%N_5%*%summation_vector)
  Pr_010011_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%R_1%*%R_2%*%N_3%*%R_4%*%N_5%*%summation_vector)
  
  Pr_101110_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%R_1%*%R_2%*%N_3%*%N_4%*%R_5%*%summation_vector)
  Pr_010001_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%R_1%*%R_2%*%N_3%*%N_4%*%R_5%*%summation_vector)
  
  Pr_100100_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%R_1%*%N_2%*%R_3%*%R_4%*%N_5%*%summation_vector)
  Pr_011011_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%R_1%*%N_2%*%R_3%*%R_4%*%N_5%*%summation_vector)
  
  Pr_100110_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%R_1%*%N_2%*%R_3%*%N_4%*%R_5%*%summation_vector)
  Pr_011001_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%R_1%*%N_2%*%R_3%*%N_4%*%R_5%*%summation_vector)
  
  Pr_100010_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%R_1%*%N_2%*%N_3%*%R_4%*%R_5%*%summation_vector)
  Pr_011101_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%R_1%*%N_2%*%N_3%*%R_4%*%R_5%*%summation_vector)
  
  Pr_110100_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%N_1%*%R_2%*%R_3%*%R_4%*%N_5%*%summation_vector)
  Pr_001011_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%N_1%*%R_2%*%R_3%*%R_4%*%N_5%*%summation_vector)
  
  Pr_110110_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%N_1%*%R_2%*%R_3%*%N_4%*%R_5%*%summation_vector)
  Pr_001001_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%N_1%*%R_2%*%R_3%*%N_4%*%R_5%*%summation_vector)
  
  Pr_110010_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%N_1%*%R_2%*%N_3%*%R_4%*%R_5%*%summation_vector)
  Pr_001101_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%N_1%*%R_2%*%N_3%*%R_4%*%R_5%*%summation_vector)
  
  Pr_111010_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%N_1%*%N_2%*%R_3%*%R_4%*%R_5%*%summation_vector)
  Pr_000101_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%N_1%*%N_2%*%R_3%*%R_4%*%R_5%*%summation_vector)
  
  # Probability of oogenesis generating an haploid gamete with marker pattern from quadruple crossover event
  
  Pr_010100_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%R_1%*%R_2%*%R_3%*%R_4%*%N_5%*%summation_vector)
  Pr_101011_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%R_1%*%R_2%*%R_3%*%R_4%*%N_5%*%summation_vector)
  
  Pr_010110_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%R_1%*%R_2%*%R_3%*%N_4%*%R_5%*%summation_vector)
  Pr_101001_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%R_1%*%R_2%*%R_3%*%N_4%*%R_5%*%summation_vector)
  
  Pr_010010_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%R_1%*%R_2%*%N_3%*%R_4%*%R_5%*%summation_vector)
  Pr_101101_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%R_1%*%R_2%*%N_3%*%R_4%*%R_5%*%summation_vector)
  
  Pr_011010_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%R_1%*%N_2%*%R_3%*%R_4%*%R_5%*%summation_vector)
  Pr_100101_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%R_1%*%N_2%*%R_3%*%R_4%*%R_5%*%summation_vector)
  
  Pr_001010_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%N_1%*%R_2%*%R_3%*%R_4%*%R_5%*%summation_vector)
  Pr_110101_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%N_1%*%R_2%*%R_3%*%R_4%*%R_5%*%summation_vector)
  
  # Probability of oogenesis generating an haploid gamete with marker pattern from quintuple crossover event
  
  Pr_101010_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%R_1%*%R_2%*%R_3%*%R_4%*%R_5%*%summation_vector)
  Pr_010101_egg <- (1/2) * ((1/p[1])%*%stationary_vector%*%R_1%*%R_2%*%R_3%*%R_4%*%R_5%*%summation_vector)
  
  # System of equations defining the egg-to-adult survival probabilities of zygotes with given marker patterns
  # A multiplicative fitness function is assumed with no second order or higher interactions included in model
  # Probability that a no exchange gamete is fertilized by X bearing sperm and survives to adulthood
  
  Viability_000000_f <- v1_wildtype_f * v2_wildtype_f * v3_wildtype_f * v4_wildtype_f * v5_wildtype_f * v6_wildtype_f * single_experimental_unit[131]
  Viability_111111_f <- v1__mutant__f * v2__mutant__f * v3__mutant__f * v4__mutant__f * v5__mutant__f * v6__mutant__f * single_experimental_unit[131]
  # Probability that a single exchange gamete is fertilized by X bearing sperm and survives to adulthood
  
  Viability_100000_f <- v1__mutant__f * v2_wildtype_f * v3_wildtype_f * v4_wildtype_f * v5_wildtype_f * v6_wildtype_f * single_experimental_unit[131]
  Viability_011111_f <- v1_wildtype_f * v2__mutant__f * v3__mutant__f * v4__mutant__f * v5__mutant__f * v6__mutant__f * single_experimental_unit[131]
  
  Viability_110000_f <- v1__mutant__f * v2__mutant__f * v3_wildtype_f * v4_wildtype_f * v5_wildtype_f * v6_wildtype_f * single_experimental_unit[131]
  Viability_001111_f <- v1_wildtype_f * v2_wildtype_f * v3__mutant__f * v4__mutant__f * v5__mutant__f * v6__mutant__f * single_experimental_unit[131]
  
  Viability_111000_f <- v1__mutant__f * v2__mutant__f * v3__mutant__f * v4_wildtype_f * v5_wildtype_f * v6_wildtype_f * single_experimental_unit[131]
  Viability_000111_f <- v1_wildtype_f * v2_wildtype_f * v3_wildtype_f * v4__mutant__f * v5__mutant__f * v6__mutant__f * single_experimental_unit[131]
  
  Viability_111100_f <- v1__mutant__f * v2__mutant__f * v3__mutant__f * v4__mutant__f * v5_wildtype_f * v6_wildtype_f * single_experimental_unit[131]
  Viability_000011_f <- v1_wildtype_f * v2_wildtype_f * v3_wildtype_f * v4_wildtype_f * v5__mutant__f * v6__mutant__f * single_experimental_unit[131]
  
  Viability_111110_f <- v1__mutant__f * v2__mutant__f * v3__mutant__f * v4__mutant__f * v5__mutant__f * v6_wildtype_f * single_experimental_unit[131]
  Viability_000001_f <- v1_wildtype_f * v2_wildtype_f * v3_wildtype_f * v4_wildtype_f * v5_wildtype_f * v6__mutant__f * single_experimental_unit[131]
  
  # Probability that a double exchange gamete is fertilized by X bearing sperm and survives to adulthood
  
  Viability_010000_f <- v1_wildtype_f * v2__mutant__f * v3_wildtype_f * v4_wildtype_f * v5_wildtype_f * v6_wildtype_f * single_experimental_unit[131]
  Viability_101111_f <- v1__mutant__f * v2_wildtype_f * v3__mutant__f * v4__mutant__f * v5__mutant__f * v6__mutant__f * single_experimental_unit[131]
  
  Viability_011000_f <- v1_wildtype_f * v2__mutant__f * v3__mutant__f * v4_wildtype_f * v5_wildtype_f * v6_wildtype_f * single_experimental_unit[131]
  Viability_100111_f <- v1__mutant__f * v2_wildtype_f * v3_wildtype_f * v4__mutant__f * v5__mutant__f * v6__mutant__f * single_experimental_unit[131]
  
  Viability_011100_f <- v1_wildtype_f * v2__mutant__f * v3__mutant__f * v4__mutant__f * v5_wildtype_f * v6_wildtype_f * single_experimental_unit[131]
  Viability_100011_f <- v1__mutant__f * v2_wildtype_f * v3_wildtype_f * v4_wildtype_f * v5__mutant__f * v6__mutant__f * single_experimental_unit[131]
  
  Viability_011110_f <- v1_wildtype_f * v2__mutant__f * v3__mutant__f * v4__mutant__f * v5__mutant__f * v6_wildtype_f * single_experimental_unit[131]
  Viability_100001_f <- v1__mutant__f * v2_wildtype_f * v3_wildtype_f * v4_wildtype_f * v5_wildtype_f * v6__mutant__f * single_experimental_unit[131]
  
  Viability_001000_f <- v1_wildtype_f * v2_wildtype_f * v3__mutant__f * v4_wildtype_f * v5_wildtype_f * v6_wildtype_f * single_experimental_unit[131]
  Viability_110111_f <- v1__mutant__f * v2__mutant__f * v3_wildtype_f * v4__mutant__f * v5__mutant__f * v6__mutant__f * single_experimental_unit[131]
  
  Viability_001100_f <- v1_wildtype_f * v2_wildtype_f * v3__mutant__f * v4__mutant__f * v5_wildtype_f * v6_wildtype_f * single_experimental_unit[131]
  Viability_110011_f <- v1__mutant__f * v2__mutant__f * v3_wildtype_f * v4_wildtype_f * v5__mutant__f * v6__mutant__f * single_experimental_unit[131]
  
  Viability_001110_f <- v1_wildtype_f * v2_wildtype_f * v3__mutant__f * v4__mutant__f * v5__mutant__f * v6_wildtype_f * single_experimental_unit[131]
  Viability_110001_f <- v1__mutant__f * v2__mutant__f * v3_wildtype_f * v4_wildtype_f * v5_wildtype_f * v6__mutant__f * single_experimental_unit[131]
  
  Viability_000100_f <- v1_wildtype_f * v2_wildtype_f * v3_wildtype_f * v4__mutant__f * v5_wildtype_f * v6_wildtype_f * single_experimental_unit[131]
  Viability_111011_f <- v1__mutant__f * v2__mutant__f * v3__mutant__f * v4_wildtype_f * v5__mutant__f * v6__mutant__f * single_experimental_unit[131]
  
  Viability_000110_f <- v1_wildtype_f * v2_wildtype_f * v3_wildtype_f * v4__mutant__f * v5__mutant__f * v6_wildtype_f * single_experimental_unit[131]
  Viability_111001_f <- v1__mutant__f * v2__mutant__f * v3__mutant__f * v4_wildtype_f * v5_wildtype_f * v6__mutant__f * single_experimental_unit[131]
  
  Viability_000010_f <- v1_wildtype_f * v2_wildtype_f * v3_wildtype_f * v4_wildtype_f * v5__mutant__f * v6_wildtype_f * single_experimental_unit[131]
  Viability_111101_f <- v1__mutant__f * v2__mutant__f * v3__mutant__f * v4__mutant__f * v5_wildtype_f * v6__mutant__f * single_experimental_unit[131]
  
  # Probability that a triple exchange gamete is fertilized by X bearing sperm and survives to adulthood
  
  Viability_101000_f <- v1__mutant__f * v2_wildtype_f * v3__mutant__f * v4_wildtype_f * v5_wildtype_f * v6_wildtype_f * single_experimental_unit[131]
  Viability_010111_f <- v1_wildtype_f * v2__mutant__f * v3_wildtype_f * v4__mutant__f * v5__mutant__f * v6__mutant__f * single_experimental_unit[131]
  
  Viability_101100_f <- v1__mutant__f * v2_wildtype_f * v3__mutant__f * v4__mutant__f * v5_wildtype_f * v6_wildtype_f * single_experimental_unit[131]
  Viability_010011_f <- v1_wildtype_f * v2__mutant__f * v3_wildtype_f * v4_wildtype_f * v5__mutant__f * v6__mutant__f * single_experimental_unit[131]
  
  Viability_101110_f <- v1__mutant__f * v2_wildtype_f * v3__mutant__f * v4__mutant__f * v5__mutant__f * v6_wildtype_f * single_experimental_unit[131]
  Viability_010001_f <- v1_wildtype_f * v2__mutant__f * v3_wildtype_f * v4_wildtype_f * v5_wildtype_f * v6__mutant__f * single_experimental_unit[131]
  
  Viability_100100_f <- v1__mutant__f * v2_wildtype_f * v3_wildtype_f * v4__mutant__f * v5_wildtype_f * v6_wildtype_f * single_experimental_unit[131]
  Viability_011011_f <- v1_wildtype_f * v2__mutant__f * v3__mutant__f * v4_wildtype_f * v5__mutant__f * v6__mutant__f * single_experimental_unit[131]
  
  Viability_100110_f <- v1__mutant__f * v2_wildtype_f * v3_wildtype_f * v4__mutant__f * v5__mutant__f * v6_wildtype_f * single_experimental_unit[131]
  Viability_011001_f <- v1_wildtype_f * v2__mutant__f * v3__mutant__f * v4_wildtype_f * v5_wildtype_f * v6__mutant__f * single_experimental_unit[131]
  
  Viability_100010_f <- v1__mutant__f * v2_wildtype_f * v3_wildtype_f * v4_wildtype_f * v5__mutant__f * v6_wildtype_f * single_experimental_unit[131]
  Viability_011101_f <- v1_wildtype_f * v2__mutant__f * v3__mutant__f * v4__mutant__f * v5_wildtype_f * v6__mutant__f * single_experimental_unit[131]
  
  Viability_110100_f <- v1__mutant__f * v2__mutant__f * v3_wildtype_f * v4__mutant__f * v5_wildtype_f * v6_wildtype_f * single_experimental_unit[131]
  Viability_001011_f <- v1_wildtype_f * v2_wildtype_f * v3__mutant__f * v4_wildtype_f * v5__mutant__f * v6__mutant__f * single_experimental_unit[131]
  
  Viability_110110_f <- v1__mutant__f * v2__mutant__f * v3_wildtype_f * v4__mutant__f * v5__mutant__f * v6_wildtype_f * single_experimental_unit[131]
  Viability_001001_f <- v1_wildtype_f * v2_wildtype_f * v3__mutant__f * v4_wildtype_f * v5_wildtype_f * v6__mutant__f * single_experimental_unit[131]
  
  Viability_110010_f <- v1__mutant__f * v2__mutant__f * v3_wildtype_f * v4_wildtype_f * v5__mutant__f * v6_wildtype_f * single_experimental_unit[131]
  Viability_001101_f <- v1_wildtype_f * v2_wildtype_f * v3__mutant__f * v4__mutant__f * v5_wildtype_f * v6__mutant__f * single_experimental_unit[131]
  
  Viability_111010_f <- v1__mutant__f * v2__mutant__f * v3__mutant__f * v4_wildtype_f * v5__mutant__f * v6_wildtype_f * single_experimental_unit[131]
  Viability_000101_f <- v1_wildtype_f * v2_wildtype_f * v3_wildtype_f * v4__mutant__f * v5_wildtype_f * v6__mutant__f * single_experimental_unit[131]
  
  # Probability that a quadruple exchange gamete is fertilized by X bearing sperm and survives to adulthood
  
  Viability_010100_f <- v1_wildtype_f * v2__mutant__f * v3_wildtype_f * v4__mutant__f * v5_wildtype_f * v6_wildtype_f * single_experimental_unit[131]
  Viability_101011_f <- v1__mutant__f * v2_wildtype_f * v3__mutant__f * v4_wildtype_f * v5__mutant__f * v6__mutant__f * single_experimental_unit[131]
  
  Viability_010110_f <- v1_wildtype_f * v2__mutant__f * v3_wildtype_f * v4__mutant__f * v5__mutant__f * v6_wildtype_f * single_experimental_unit[131]
  Viability_101001_f <- v1__mutant__f * v2_wildtype_f * v3__mutant__f * v4_wildtype_f * v5_wildtype_f * v6__mutant__f * single_experimental_unit[131]
  
  Viability_010010_f <- v1_wildtype_f * v2__mutant__f * v3_wildtype_f * v4_wildtype_f * v5__mutant__f * v6_wildtype_f * single_experimental_unit[131]
  Viability_101101_f <- v1__mutant__f * v2_wildtype_f * v3__mutant__f * v4__mutant__f * v5_wildtype_f * v6__mutant__f * single_experimental_unit[131]
  
  Viability_011010_f <- v1_wildtype_f * v2__mutant__f * v3__mutant__f * v4_wildtype_f * v5__mutant__f * v6_wildtype_f * single_experimental_unit[131]
  Viability_100101_f <- v1__mutant__f * v2_wildtype_f * v3_wildtype_f * v4__mutant__f * v5_wildtype_f * v6__mutant__f * single_experimental_unit[131]
  
  Viability_001010_f <- v1_wildtype_f * v2_wildtype_f * v3__mutant__f * v4_wildtype_f * v5__mutant__f * v6_wildtype_f * single_experimental_unit[131]
  Viability_110101_f <- v1__mutant__f * v2__mutant__f * v3_wildtype_f * v4__mutant__f * v5_wildtype_f * v6__mutant__f * single_experimental_unit[131]
  
  # Probability that a quintuple exchange gamete is fertilized by X bearing sperm and survives to adulthood
  
  Viability_101010_f <- v1__mutant__f * v2_wildtype_f * v3__mutant__f * v4_wildtype_f * v5__mutant__f * v6_wildtype_f * single_experimental_unit[131]
  Viability_010101_f <- v1_wildtype_f * v2__mutant__f * v3_wildtype_f * v4__mutant__f * v5_wildtype_f * v6__mutant__f * single_experimental_unit[131]
  
  # Probability that a no exchange gamete is fertilized by Y bearing sperm and survives to adulthood
  
  Viability_000000_m <- v1_wildtype_m * v2_wildtype_m * v3_wildtype_m * v4_wildtype_m * v5_wildtype_m * v6_wildtype_m * single_experimental_unit[132]
  Viability_111111_m <- v1__mutant__m * v2__mutant__m * v3__mutant__m * v4__mutant__m * v5__mutant__m * v6__mutant__m * single_experimental_unit[132]
  
  # Probability that a single exchange gamete is fertilized by Y bearing sperm and survives to adulthood
  
  Viability_100000_m <- v1__mutant__m * v2_wildtype_m * v3_wildtype_m * v4_wildtype_m * v5_wildtype_m * v6_wildtype_m * single_experimental_unit[132]
  Viability_011111_m <- v1_wildtype_m * v2__mutant__m * v3__mutant__m * v4__mutant__m * v5__mutant__m * v6__mutant__m * single_experimental_unit[132]
  
  Viability_110000_m <- v1__mutant__m * v2__mutant__m * v3_wildtype_m * v4_wildtype_m * v5_wildtype_m * v6_wildtype_m * single_experimental_unit[132]
  Viability_001111_m <- v1_wildtype_m * v2_wildtype_m * v3__mutant__m * v4__mutant__m * v5__mutant__m * v6__mutant__m * single_experimental_unit[132]
  
  Viability_111000_m <- v1__mutant__m * v2__mutant__m * v3__mutant__m * v4_wildtype_m * v5_wildtype_m * v6_wildtype_m * single_experimental_unit[132]
  Viability_000111_m <- v1_wildtype_m * v2_wildtype_m * v3_wildtype_m * v4__mutant__m * v5__mutant__m * v6__mutant__m * single_experimental_unit[132]
  
  Viability_111100_m <- v1__mutant__m * v2__mutant__m * v3__mutant__m * v4__mutant__m * v5_wildtype_m * v6_wildtype_m * single_experimental_unit[132]
  Viability_000011_m <- v1_wildtype_m * v2_wildtype_m * v3_wildtype_m * v4_wildtype_m * v5__mutant__m * v6__mutant__m * single_experimental_unit[132]
  
  Viability_111110_m <- v1__mutant__m * v2__mutant__m * v3__mutant__m * v4__mutant__m * v5__mutant__m * v6_wildtype_m * single_experimental_unit[132]
  Viability_000001_m <- v1_wildtype_m * v2_wildtype_m * v3_wildtype_m * v4_wildtype_m * v5_wildtype_m * v6__mutant__m * single_experimental_unit[132]
  
  # Probability that a double exchange gamete is fertilized by Y bearing sperm and survives to adulthood
  
  Viability_010000_m <- v1_wildtype_m * v2__mutant__m * v3_wildtype_m * v4_wildtype_m * v5_wildtype_m * v6_wildtype_m * single_experimental_unit[132]
  Viability_101111_m <- v1__mutant__m * v2_wildtype_m * v3__mutant__m * v4__mutant__m * v5__mutant__m * v6__mutant__m * single_experimental_unit[132]
  
  Viability_011000_m <- v1_wildtype_m * v2__mutant__m * v3__mutant__m * v4_wildtype_m * v5_wildtype_m * v6_wildtype_m * single_experimental_unit[132]
  Viability_100111_m <- v1__mutant__m * v2_wildtype_m * v3_wildtype_m * v4__mutant__m * v5__mutant__m * v6__mutant__m * single_experimental_unit[132]
  
  Viability_011100_m <- v1_wildtype_m * v2__mutant__m * v3__mutant__m * v4__mutant__m * v5_wildtype_m * v6_wildtype_m * single_experimental_unit[132]
  Viability_100011_m <- v1__mutant__m * v2_wildtype_m * v3_wildtype_m * v4_wildtype_m * v5__mutant__m * v6__mutant__m * single_experimental_unit[132]
  
  Viability_011110_m <- v1_wildtype_m * v2__mutant__m * v3__mutant__m * v4__mutant__m * v5__mutant__m * v6_wildtype_m * single_experimental_unit[132]
  Viability_100001_m <- v1__mutant__m * v2_wildtype_m * v3_wildtype_m * v4_wildtype_m * v5_wildtype_m * v6__mutant__m * single_experimental_unit[132]
  
  Viability_001000_m <- v1_wildtype_m * v2_wildtype_m * v3__mutant__m * v4_wildtype_m * v5_wildtype_m * v6_wildtype_m * single_experimental_unit[132]
  Viability_110111_m <- v1__mutant__m * v2__mutant__m * v3_wildtype_m * v4__mutant__m * v5__mutant__m * v6__mutant__m * single_experimental_unit[132]
  
  Viability_001100_m <- v1_wildtype_m * v2_wildtype_m * v3__mutant__m * v4__mutant__m * v5_wildtype_m * v6_wildtype_m * single_experimental_unit[132]
  Viability_110011_m <- v1__mutant__m * v2__mutant__m * v3_wildtype_m * v4_wildtype_m * v5__mutant__m * v6__mutant__m * single_experimental_unit[132]
  
  Viability_001110_m <- v1_wildtype_m * v2_wildtype_m * v3__mutant__m * v4__mutant__m * v5__mutant__m * v6_wildtype_m * single_experimental_unit[132]
  Viability_110001_m <- v1__mutant__m * v2__mutant__m * v3_wildtype_m * v4_wildtype_m * v5_wildtype_m * v6__mutant__m * single_experimental_unit[132]
  
  Viability_000100_m <- v1_wildtype_m * v2_wildtype_m * v3_wildtype_m * v4__mutant__m * v5_wildtype_m * v6_wildtype_m * single_experimental_unit[132]
  Viability_111011_m <- v1__mutant__m * v2__mutant__m * v3__mutant__m * v4_wildtype_m * v5__mutant__m * v6__mutant__m * single_experimental_unit[132]
  
  Viability_000110_m <- v1_wildtype_m * v2_wildtype_m * v3_wildtype_m * v4__mutant__m * v5__mutant__m * v6_wildtype_m * single_experimental_unit[132]
  Viability_111001_m <- v1__mutant__m * v2__mutant__m * v3__mutant__m * v4_wildtype_m * v5_wildtype_m * v6__mutant__m * single_experimental_unit[132]
  
  Viability_000010_m <- v1_wildtype_m * v2_wildtype_m * v3_wildtype_m * v4_wildtype_m * v5__mutant__m * v6_wildtype_m * single_experimental_unit[132]
  Viability_111101_m <- v1__mutant__m * v2__mutant__m * v3__mutant__m * v4__mutant__m * v5_wildtype_m * v6__mutant__m * single_experimental_unit[132]
  
  # Probability that a triple exchange gamete is fertilized by Y bearing sperm and survives to adulthood
  
  Viability_101000_m <- v1__mutant__m * v2_wildtype_m * v3__mutant__m * v4_wildtype_m * v5_wildtype_m * v6_wildtype_m * single_experimental_unit[132]
  Viability_010111_m <- v1_wildtype_m * v2__mutant__m * v3_wildtype_m * v4__mutant__m * v5__mutant__m * v6__mutant__m * single_experimental_unit[132]
  
  Viability_101100_m <- v1__mutant__m * v2_wildtype_m * v3__mutant__m * v4__mutant__m * v5_wildtype_m * v6_wildtype_m * single_experimental_unit[132]
  Viability_010011_m <- v1_wildtype_m * v2__mutant__m * v3_wildtype_m * v4_wildtype_m * v5__mutant__m * v6__mutant__m * single_experimental_unit[132]
  
  Viability_101110_m <- v1__mutant__m * v2_wildtype_m * v3__mutant__m * v4__mutant__m * v5__mutant__m * v6_wildtype_m * single_experimental_unit[132]
  Viability_010001_m <- v1_wildtype_m * v2__mutant__m * v3_wildtype_m * v4_wildtype_m * v5_wildtype_m * v6__mutant__m * single_experimental_unit[132]
  
  Viability_100100_m <- v1__mutant__m * v2_wildtype_m * v3_wildtype_m * v4__mutant__m * v5_wildtype_m * v6_wildtype_m * single_experimental_unit[132]
  Viability_011011_m <- v1_wildtype_m * v2__mutant__m * v3__mutant__m * v4_wildtype_m * v5__mutant__m * v6__mutant__m * single_experimental_unit[132]
  
  Viability_100110_m <- v1__mutant__m * v2_wildtype_m * v3_wildtype_m * v4__mutant__m * v5__mutant__m * v6_wildtype_m * single_experimental_unit[132]
  Viability_011001_m <- v1_wildtype_m * v2__mutant__m * v3__mutant__m * v4_wildtype_m * v5_wildtype_m * v6__mutant__m * single_experimental_unit[132]
  
  Viability_100010_m <- v1__mutant__m * v2_wildtype_m * v3_wildtype_m * v4_wildtype_m * v5__mutant__m * v6_wildtype_m * single_experimental_unit[132]
  Viability_011101_m <- v1_wildtype_m * v2__mutant__m * v3__mutant__m * v4__mutant__m * v5_wildtype_m * v6__mutant__m * single_experimental_unit[132]
  
  Viability_110100_m <- v1__mutant__m * v2__mutant__m * v3_wildtype_m * v4__mutant__m * v5_wildtype_m * v6_wildtype_m * single_experimental_unit[132]
  Viability_001011_m <- v1_wildtype_m * v2_wildtype_m * v3__mutant__m * v4_wildtype_m * v5__mutant__m * v6__mutant__m * single_experimental_unit[132]
  
  Viability_110110_m <- v1__mutant__m * v2__mutant__m * v3_wildtype_m * v4__mutant__m * v5__mutant__m * v6_wildtype_m * single_experimental_unit[132]
  Viability_001001_m <- v1_wildtype_m * v2_wildtype_m * v3__mutant__m * v4_wildtype_m * v5_wildtype_m * v6__mutant__m * single_experimental_unit[132]
  
  Viability_110010_m <- v1__mutant__m * v2__mutant__m * v3_wildtype_m * v4_wildtype_m * v5__mutant__m * v6_wildtype_m * single_experimental_unit[132]
  Viability_001101_m <- v1_wildtype_m * v2_wildtype_m * v3__mutant__m * v4__mutant__m * v5_wildtype_m * v6__mutant__m * single_experimental_unit[132]
  
  Viability_111010_m <- v1__mutant__m * v2__mutant__m * v3__mutant__m * v4_wildtype_m * v5__mutant__m * v6_wildtype_m * single_experimental_unit[132]
  Viability_000101_m <- v1_wildtype_m * v2_wildtype_m * v3_wildtype_m * v4__mutant__m * v5_wildtype_m * v6__mutant__m * single_experimental_unit[132]
  
  # Probability that a quadruple exchange gamete is fertilized by Y bearing sperm and survives to adulthood
  
  Viability_010100_m <- v1_wildtype_m * v2__mutant__m * v3_wildtype_m * v4__mutant__m * v5_wildtype_m * v6_wildtype_m * single_experimental_unit[132]
  Viability_101011_m <- v1__mutant__m * v2_wildtype_m * v3__mutant__m * v4_wildtype_m * v5__mutant__m * v6__mutant__m * single_experimental_unit[132]
  
  Viability_010110_m <- v1_wildtype_m * v2__mutant__m * v3_wildtype_m * v4__mutant__m * v5__mutant__m * v6_wildtype_m * single_experimental_unit[132]
  Viability_101001_m <- v1__mutant__m * v2_wildtype_m * v3__mutant__m * v4_wildtype_m * v5_wildtype_m * v6__mutant__m * single_experimental_unit[132]
  
  Viability_010010_m <- v1_wildtype_m * v2__mutant__m * v3_wildtype_m * v4_wildtype_m * v5__mutant__m * v6_wildtype_m * single_experimental_unit[132]
  Viability_101101_m <- v1__mutant__m * v2_wildtype_m * v3__mutant__m * v4__mutant__m * v5_wildtype_m * v6__mutant__m * single_experimental_unit[132]
  
  Viability_011010_m <- v1_wildtype_m * v2__mutant__m * v3__mutant__m * v4_wildtype_m * v5__mutant__m * v6_wildtype_m * single_experimental_unit[132]
  Viability_100101_m <- v1__mutant__m * v2_wildtype_m * v3_wildtype_m * v4__mutant__m * v5_wildtype_m * v6__mutant__m * single_experimental_unit[132]
  
  Viability_001010_m <- v1_wildtype_m * v2_wildtype_m * v3__mutant__m * v4_wildtype_m * v5__mutant__m * v6_wildtype_m * single_experimental_unit[132]
  Viability_110101_m <- v1__mutant__m * v2__mutant__m * v3_wildtype_m * v4__mutant__m * v5_wildtype_m * v6__mutant__m * single_experimental_unit[132]
  
  # Probability that a quintuple exchange gamete is fertilized by Y bearing sperm and survives to adulthood
  
  Viability_101010_m <- v1__mutant__m * v2_wildtype_m * v3__mutant__m * v4_wildtype_m * v5__mutant__m * v6_wildtype_m * single_experimental_unit[132]
  Viability_010101_m <- v1_wildtype_m * v2__mutant__m * v3_wildtype_m * v4__mutant__m * v5_wildtype_m * v6__mutant__m * single_experimental_unit[132]
  
  # System of equations defining expected frequency of the sex-specific observable adult phenotypic classes
  # Expected frequency of each class is the product of probability X or Y bearing sperm fertilizing the egg,
  # The probability of each specific marker complement is formed and segregates to the functional egg pole,
  # And the probability that each type of sex-specific zygote will survive to be scored as an F2 adult fly.
  # Probability of observing a female with recombination pattern from no exchange
  
  exp_000000_f <- (Pr_X_bearing_sperm) * (Pr_000000_egg) * (Viability_000000_f)
  exp_111111_f <- (Pr_X_bearing_sperm) * (Pr_111111_egg) * (Viability_111111_f)
  
  # Probability of observing a female with single exchange recombination pattern
  
  exp_100000_f <- (Pr_X_bearing_sperm) * (Pr_100000_egg) * (Viability_100000_f)
  exp_011111_f <- (Pr_X_bearing_sperm) * (Pr_011111_egg) * (Viability_011111_f)
  
  exp_110000_f <- (Pr_X_bearing_sperm) * (Pr_110000_egg) * (Viability_110000_f)
  exp_001111_f <- (Pr_X_bearing_sperm) * (Pr_001111_egg) * (Viability_001111_f)
  
  exp_111000_f <- (Pr_X_bearing_sperm) * (Pr_111000_egg) * (Viability_111000_f)
  exp_000111_f <- (Pr_X_bearing_sperm) * (Pr_000111_egg) * (Viability_000111_f)
  
  exp_111100_f <- (Pr_X_bearing_sperm) * (Pr_111100_egg) * (Viability_111100_f)
  exp_000011_f <- (Pr_X_bearing_sperm) * (Pr_000011_egg) * (Viability_000011_f)
  
  exp_111110_f <- (Pr_X_bearing_sperm) * (Pr_111110_egg) * (Viability_111110_f)
  exp_000001_f <- (Pr_X_bearing_sperm) * (Pr_000001_egg) * (Viability_000001_f)
  
  # Probability of observing a female with double exchange recombination pattern
  
  exp_010000_f <- (Pr_X_bearing_sperm) * (Pr_010000_egg) * (Viability_010000_f)
  exp_101111_f <- (Pr_X_bearing_sperm) * (Pr_101111_egg) * (Viability_101111_f)
  
  exp_011000_f <- (Pr_X_bearing_sperm) * (Pr_011000_egg) * (Viability_011000_f)
  exp_100111_f <- (Pr_X_bearing_sperm) * (Pr_100111_egg) * (Viability_100111_f)
  
  exp_011100_f <- (Pr_X_bearing_sperm) * (Pr_011100_egg) * (Viability_011100_f)
  exp_100011_f <- (Pr_X_bearing_sperm) * (Pr_100011_egg) * (Viability_100011_f)
  
  exp_011110_f <- (Pr_X_bearing_sperm) * (Pr_011110_egg) * (Viability_011110_f)
  exp_100001_f <- (Pr_X_bearing_sperm) * (Pr_100001_egg) * (Viability_100001_f)
  
  exp_001000_f <- (Pr_X_bearing_sperm) * (Pr_001000_egg) * (Viability_001000_f)
  exp_110111_f <- (Pr_X_bearing_sperm) * (Pr_110111_egg) * (Viability_110111_f)
  
  exp_001100_f <- (Pr_X_bearing_sperm) * (Pr_001100_egg) * (Viability_001100_f)
  exp_110011_f <- (Pr_X_bearing_sperm) * (Pr_110011_egg) * (Viability_110011_f)
  
  exp_001110_f <- (Pr_X_bearing_sperm) * (Pr_001110_egg) * (Viability_001110_f)
  exp_110001_f <- (Pr_X_bearing_sperm) * (Pr_110001_egg) * (Viability_110001_f)
  
  exp_000100_f <- (Pr_X_bearing_sperm) * (Pr_000100_egg) * (Viability_000100_f)
  exp_111011_f <- (Pr_X_bearing_sperm) * (Pr_111011_egg) * (Viability_111011_f)
  
  exp_000110_f <- (Pr_X_bearing_sperm) * (Pr_000110_egg) * (Viability_000110_f)
  exp_111001_f <- (Pr_X_bearing_sperm) * (Pr_111001_egg) * (Viability_111001_f)
  
  exp_000010_f <- (Pr_X_bearing_sperm) * (Pr_000010_egg) * (Viability_000010_f)
  exp_111101_f <- (Pr_X_bearing_sperm) * (Pr_111101_egg) * (Viability_111101_f)
  
  # Probability of observing a female with triple exchange recombination pattern
  
  exp_101000_f <- (Pr_X_bearing_sperm) * (Pr_101000_egg) * (Viability_101000_f)
  exp_010111_f <- (Pr_X_bearing_sperm) * (Pr_010111_egg) * (Viability_010111_f)
  
  exp_101100_f <- (Pr_X_bearing_sperm) * (Pr_101100_egg) * (Viability_101100_f)
  exp_010011_f <- (Pr_X_bearing_sperm) * (Pr_010011_egg) * (Viability_010011_f)
  
  exp_101110_f <- (Pr_X_bearing_sperm) * (Pr_101110_egg) * (Viability_101110_f)
  exp_010001_f <- (Pr_X_bearing_sperm) * (Pr_010001_egg) * (Viability_010001_f)
  
  exp_100100_f <- (Pr_X_bearing_sperm) * (Pr_100100_egg) * (Viability_100100_f)
  exp_011011_f <- (Pr_X_bearing_sperm) * (Pr_011011_egg) * (Viability_011011_f)
  
  exp_100110_f <- (Pr_X_bearing_sperm) * (Pr_100110_egg) * (Viability_100110_f)
  exp_011001_f <- (Pr_X_bearing_sperm) * (Pr_011001_egg) * (Viability_011001_f)
  
  exp_100010_f <- (Pr_X_bearing_sperm) * (Pr_100010_egg) * (Viability_100010_f)
  exp_011101_f <- (Pr_X_bearing_sperm) * (Pr_011101_egg) * (Viability_011101_f)
  
  exp_110100_f <- (Pr_X_bearing_sperm) * (Pr_110100_egg) * (Viability_110100_f)
  exp_001011_f <- (Pr_X_bearing_sperm) * (Pr_001011_egg) * (Viability_001011_f)
  
  exp_110110_f <- (Pr_X_bearing_sperm) * (Pr_110110_egg) * (Viability_110110_f)
  exp_001001_f <- (Pr_X_bearing_sperm) * (Pr_001001_egg) * (Viability_001001_f)
  
  exp_110010_f <- (Pr_X_bearing_sperm) * (Pr_110010_egg) * (Viability_110010_f)
  exp_001101_f <- (Pr_X_bearing_sperm) * (Pr_001101_egg) * (Viability_001101_f)
  
  exp_111010_f <- (Pr_X_bearing_sperm) * (Pr_111010_egg) * (Viability_111010_f)
  exp_000101_f <- (Pr_X_bearing_sperm) * (Pr_000101_egg) * (Viability_000101_f)
  
  # Probability of observing female with quadruple exchange recombination pattern
  
  exp_010100_f <- (Pr_X_bearing_sperm) * (Pr_010100_egg) * (Viability_010100_f)
  exp_101011_f <- (Pr_X_bearing_sperm) * (Pr_101011_egg) * (Viability_101011_f)
  
  exp_010110_f <- (Pr_X_bearing_sperm) * (Pr_010110_egg) * (Viability_010110_f)
  exp_101001_f <- (Pr_X_bearing_sperm) * (Pr_101001_egg) * (Viability_101001_f)
  
  exp_010010_f <- (Pr_X_bearing_sperm) * (Pr_010010_egg) * (Viability_010010_f)
  exp_101101_f <- (Pr_X_bearing_sperm) * (Pr_101101_egg) * (Viability_101101_f)
  
  exp_011010_f <- (Pr_X_bearing_sperm) * (Pr_011010_egg) * (Viability_011010_f)
  exp_100101_f <- (Pr_X_bearing_sperm) * (Pr_100101_egg) * (Viability_100101_f)
  
  exp_001010_f <- (Pr_X_bearing_sperm) * (Pr_001010_egg) * (Viability_001010_f)
  exp_110101_f <- (Pr_X_bearing_sperm) * (Pr_110101_egg) * (Viability_110101_f)
  
  # Probability of observing female with quintuple exchange recombination pattern
  
  exp_101010_f <- (Pr_X_bearing_sperm) * (Pr_101010_egg) * (Viability_101010_f)
  exp_010101_f <- (Pr_X_bearing_sperm) * (Pr_010101_egg) * (Viability_010101_f)
  
  # Cumulative probability of observing females with higher order exchanges (>2)
  
  exp_high_order_exchange_f <- sum(exp_101000_f,
                                   exp_010111_f,
                                   exp_101100_f,
                                   exp_010011_f,
                                   exp_101110_f,
                                   exp_010001_f,
                                   exp_100100_f,
                                   exp_011011_f,
                                   exp_100110_f,
                                   exp_011001_f,
                                   exp_100010_f,
                                   exp_011101_f,
                                   exp_110100_f,
                                   exp_001011_f,
                                   exp_110110_f,
                                   exp_001001_f,
                                   exp_110010_f,
                                   exp_001101_f,
                                   exp_111010_f,
                                   exp_000101_f,
                                   exp_010100_f,
                                   exp_101011_f,
                                   exp_010110_f,
                                   exp_101001_f,
                                   exp_010010_f,
                                   exp_101101_f,
                                   exp_011010_f,
                                   exp_100101_f,
                                   exp_001010_f,
                                   exp_110101_f,
                                   exp_101010_f,
                                   exp_010101_f)
  
  # Probability of a female zygote failing to develop to an observable adult fly
  
  exp_failed_development_f <- sum((Pr_X_bearing_sperm) * (Pr_000000_egg) * (1-Viability_000000_f),
                                  (Pr_X_bearing_sperm) * (Pr_111111_egg) * (1-Viability_111111_f),
                                  (Pr_X_bearing_sperm) * (Pr_100000_egg) * (1-Viability_100000_f),
                                  (Pr_X_bearing_sperm) * (Pr_011111_egg) * (1-Viability_011111_f),
                                  (Pr_X_bearing_sperm) * (Pr_110000_egg) * (1-Viability_110000_f),
                                  (Pr_X_bearing_sperm) * (Pr_001111_egg) * (1-Viability_001111_f),
                                  (Pr_X_bearing_sperm) * (Pr_111000_egg) * (1-Viability_111000_f),
                                  (Pr_X_bearing_sperm) * (Pr_000111_egg) * (1-Viability_000111_f),
                                  (Pr_X_bearing_sperm) * (Pr_111100_egg) * (1-Viability_111100_f),
                                  (Pr_X_bearing_sperm) * (Pr_000011_egg) * (1-Viability_000011_f),
                                  (Pr_X_bearing_sperm) * (Pr_111110_egg) * (1-Viability_111110_f),
                                  (Pr_X_bearing_sperm) * (Pr_000001_egg) * (1-Viability_000001_f),
                                  (Pr_X_bearing_sperm) * (Pr_010000_egg) * (1-Viability_010000_f),
                                  (Pr_X_bearing_sperm) * (Pr_101111_egg) * (1-Viability_101111_f),
                                  (Pr_X_bearing_sperm) * (Pr_011000_egg) * (1-Viability_011000_f),
                                  (Pr_X_bearing_sperm) * (Pr_100111_egg) * (1-Viability_100111_f),
                                  (Pr_X_bearing_sperm) * (Pr_011100_egg) * (1-Viability_011100_f),
                                  (Pr_X_bearing_sperm) * (Pr_100011_egg) * (1-Viability_100011_f),
                                  (Pr_X_bearing_sperm) * (Pr_011110_egg) * (1-Viability_011110_f),
                                  (Pr_X_bearing_sperm) * (Pr_100001_egg) * (1-Viability_100001_f),
                                  (Pr_X_bearing_sperm) * (Pr_001000_egg) * (1-Viability_001000_f),
                                  (Pr_X_bearing_sperm) * (Pr_110111_egg) * (1-Viability_110111_f),
                                  (Pr_X_bearing_sperm) * (Pr_001100_egg) * (1-Viability_001100_f),
                                  (Pr_X_bearing_sperm) * (Pr_110011_egg) * (1-Viability_110011_f),
                                  (Pr_X_bearing_sperm) * (Pr_001110_egg) * (1-Viability_001110_f),
                                  (Pr_X_bearing_sperm) * (Pr_110001_egg) * (1-Viability_110001_f),
                                  (Pr_X_bearing_sperm) * (Pr_000100_egg) * (1-Viability_000100_f),
                                  (Pr_X_bearing_sperm) * (Pr_111011_egg) * (1-Viability_111011_f),
                                  (Pr_X_bearing_sperm) * (Pr_000110_egg) * (1-Viability_000110_f),
                                  (Pr_X_bearing_sperm) * (Pr_111001_egg) * (1-Viability_111001_f),
                                  (Pr_X_bearing_sperm) * (Pr_000010_egg) * (1-Viability_000010_f),
                                  (Pr_X_bearing_sperm) * (Pr_111101_egg) * (1-Viability_111101_f),
                                  (Pr_X_bearing_sperm) * (Pr_101000_egg) * (1-Viability_101000_f),
                                  (Pr_X_bearing_sperm) * (Pr_010111_egg) * (1-Viability_010111_f),
                                  (Pr_X_bearing_sperm) * (Pr_101100_egg) * (1-Viability_101100_f),
                                  (Pr_X_bearing_sperm) * (Pr_010011_egg) * (1-Viability_010011_f),
                                  (Pr_X_bearing_sperm) * (Pr_101110_egg) * (1-Viability_101110_f),
                                  (Pr_X_bearing_sperm) * (Pr_010001_egg) * (1-Viability_010001_f),
                                  (Pr_X_bearing_sperm) * (Pr_100100_egg) * (1-Viability_100100_f),
                                  (Pr_X_bearing_sperm) * (Pr_011011_egg) * (1-Viability_011011_f),
                                  (Pr_X_bearing_sperm) * (Pr_100110_egg) * (1-Viability_100110_f),
                                  (Pr_X_bearing_sperm) * (Pr_011001_egg) * (1-Viability_011001_f),
                                  (Pr_X_bearing_sperm) * (Pr_100010_egg) * (1-Viability_100010_f),
                                  (Pr_X_bearing_sperm) * (Pr_011101_egg) * (1-Viability_011101_f),
                                  (Pr_X_bearing_sperm) * (Pr_110100_egg) * (1-Viability_110100_f),
                                  (Pr_X_bearing_sperm) * (Pr_001011_egg) * (1-Viability_001011_f),
                                  (Pr_X_bearing_sperm) * (Pr_110110_egg) * (1-Viability_110110_f),
                                  (Pr_X_bearing_sperm) * (Pr_001001_egg) * (1-Viability_001001_f),
                                  (Pr_X_bearing_sperm) * (Pr_110010_egg) * (1-Viability_110010_f),
                                  (Pr_X_bearing_sperm) * (Pr_001101_egg) * (1-Viability_001101_f),
                                  (Pr_X_bearing_sperm) * (Pr_111010_egg) * (1-Viability_111010_f),
                                  (Pr_X_bearing_sperm) * (Pr_000101_egg) * (1-Viability_000101_f),
                                  (Pr_X_bearing_sperm) * (Pr_010100_egg) * (1-Viability_010100_f),
                                  (Pr_X_bearing_sperm) * (Pr_101011_egg) * (1-Viability_101011_f),
                                  (Pr_X_bearing_sperm) * (Pr_010110_egg) * (1-Viability_010110_f),
                                  (Pr_X_bearing_sperm) * (Pr_101001_egg) * (1-Viability_101001_f),
                                  (Pr_X_bearing_sperm) * (Pr_010010_egg) * (1-Viability_010010_f),
                                  (Pr_X_bearing_sperm) * (Pr_101101_egg) * (1-Viability_101101_f),
                                  (Pr_X_bearing_sperm) * (Pr_011010_egg) * (1-Viability_011010_f),
                                  (Pr_X_bearing_sperm) * (Pr_100101_egg) * (1-Viability_100101_f),
                                  (Pr_X_bearing_sperm) * (Pr_001010_egg) * (1-Viability_001010_f),
                                  (Pr_X_bearing_sperm) * (Pr_110101_egg) * (1-Viability_110101_f),
                                  (Pr_X_bearing_sperm) * (Pr_101010_egg) * (1-Viability_101010_f),
                                  (Pr_X_bearing_sperm) * (Pr_010101_egg) * (1-Viability_010101_f))
  
  # Probability of observing a male with recombination pattern from no exchange
  
  exp_000000_m <- (Pr_Y_bearing_sperm) * (Pr_000000_egg) * (Viability_000000_m)
  exp_111111_m <- (Pr_Y_bearing_sperm) * (Pr_111111_egg) * (Viability_111111_m)
  
  # Probability of observing a male with single exchange recombination pattern
  
  exp_100000_m <- (Pr_Y_bearing_sperm) * (Pr_100000_egg) * (Viability_100000_m)
  exp_011111_m <- (Pr_Y_bearing_sperm) * (Pr_011111_egg) * (Viability_011111_m)
  
  exp_110000_m <- (Pr_Y_bearing_sperm) * (Pr_110000_egg) * (Viability_110000_m)
  exp_001111_m <- (Pr_Y_bearing_sperm) * (Pr_001111_egg) * (Viability_001111_m)
  
  exp_111000_m <- (Pr_Y_bearing_sperm) * (Pr_111000_egg) * (Viability_111000_m)
  exp_000111_m <- (Pr_Y_bearing_sperm) * (Pr_000111_egg) * (Viability_000111_m)
  
  exp_111100_m <- (Pr_Y_bearing_sperm) * (Pr_111100_egg) * (Viability_111100_m)
  exp_000011_m <- (Pr_Y_bearing_sperm) * (Pr_000011_egg) * (Viability_000011_m)
  
  exp_111110_m <- (Pr_Y_bearing_sperm) * (Pr_111110_egg) * (Viability_111110_m)
  exp_000001_m <- (Pr_Y_bearing_sperm) * (Pr_000001_egg) * (Viability_000001_m)
  
  # Probability of observing a male with double exchange recombination pattern
  
  exp_010000_m <- (Pr_Y_bearing_sperm) * (Pr_010000_egg) * (Viability_010000_m)
  exp_101111_m <- (Pr_Y_bearing_sperm) * (Pr_101111_egg) * (Viability_101111_m)
  
  exp_011000_m <- (Pr_Y_bearing_sperm) * (Pr_011000_egg) * (Viability_011000_m)
  exp_100111_m <- (Pr_Y_bearing_sperm) * (Pr_100111_egg) * (Viability_100111_m)
  
  exp_011100_m <- (Pr_Y_bearing_sperm) * (Pr_011100_egg) * (Viability_011100_m)
  exp_100011_m <- (Pr_Y_bearing_sperm) * (Pr_100011_egg) * (Viability_100011_m)
  
  exp_011110_m <- (Pr_Y_bearing_sperm) * (Pr_011110_egg) * (Viability_011110_m)
  exp_100001_m <- (Pr_Y_bearing_sperm) * (Pr_100001_egg) * (Viability_100001_m)
  
  exp_001000_m <- (Pr_Y_bearing_sperm) * (Pr_001000_egg) * (Viability_001000_m)
  exp_110111_m <- (Pr_Y_bearing_sperm) * (Pr_110111_egg) * (Viability_110111_m)
  
  exp_001100_m <- (Pr_Y_bearing_sperm) * (Pr_001100_egg) * (Viability_001100_m)
  exp_110011_m <- (Pr_Y_bearing_sperm) * (Pr_110011_egg) * (Viability_110011_m)
  
  exp_001110_m <- (Pr_Y_bearing_sperm) * (Pr_001110_egg) * (Viability_001110_m)
  exp_110001_m <- (Pr_Y_bearing_sperm) * (Pr_110001_egg) * (Viability_110001_m)
  
  exp_000100_m <- (Pr_Y_bearing_sperm) * (Pr_000100_egg) * (Viability_000100_m)
  exp_111011_m <- (Pr_Y_bearing_sperm) * (Pr_111011_egg) * (Viability_111011_m)
  
  exp_000110_m <- (Pr_Y_bearing_sperm) * (Pr_000110_egg) * (Viability_000110_m)
  exp_111001_m <- (Pr_Y_bearing_sperm) * (Pr_111001_egg) * (Viability_111001_m)
  
  exp_000010_m <- (Pr_Y_bearing_sperm) * (Pr_000010_egg) * (Viability_000010_m)
  exp_111101_m <- (Pr_Y_bearing_sperm) * (Pr_111101_egg) * (Viability_111101_m)
  
  # Probability of observing a male with triple exchange recombination pattern
  
  exp_101000_m <- (Pr_Y_bearing_sperm) * (Pr_101000_egg) * (Viability_101000_m)
  exp_010111_m <- (Pr_Y_bearing_sperm) * (Pr_010111_egg) * (Viability_010111_m)
  
  exp_101100_m <- (Pr_Y_bearing_sperm) * (Pr_101100_egg) * (Viability_101100_m)
  exp_010011_m <- (Pr_Y_bearing_sperm) * (Pr_010011_egg) * (Viability_010011_m)
  
  exp_101110_m <- (Pr_Y_bearing_sperm) * (Pr_101110_egg) * (Viability_101110_m)
  exp_010001_m <- (Pr_Y_bearing_sperm) * (Pr_010001_egg) * (Viability_010001_m)
  
  exp_100100_m <- (Pr_Y_bearing_sperm) * (Pr_100100_egg) * (Viability_100100_m)
  exp_011011_m <- (Pr_Y_bearing_sperm) * (Pr_011011_egg) * (Viability_011011_m)
  
  exp_100110_m <- (Pr_Y_bearing_sperm) * (Pr_100110_egg) * (Viability_100110_m)
  exp_011001_m <- (Pr_Y_bearing_sperm) * (Pr_011001_egg) * (Viability_011001_m)
  
  exp_100010_m <- (Pr_Y_bearing_sperm) * (Pr_100010_egg) * (Viability_100010_m)
  exp_011101_m <- (Pr_Y_bearing_sperm) * (Pr_011101_egg) * (Viability_011101_m)
  
  exp_110100_m <- (Pr_Y_bearing_sperm) * (Pr_110100_egg) * (Viability_110100_m)
  exp_001011_m <- (Pr_Y_bearing_sperm) * (Pr_001011_egg) * (Viability_001011_m)
  
  exp_110110_m <- (Pr_Y_bearing_sperm) * (Pr_110110_egg) * (Viability_110110_m)
  exp_001001_m <- (Pr_Y_bearing_sperm) * (Pr_001001_egg) * (Viability_001001_m)
  
  exp_110010_m <- (Pr_Y_bearing_sperm) * (Pr_110010_egg) * (Viability_110010_m)
  exp_001101_m <- (Pr_Y_bearing_sperm) * (Pr_001101_egg) * (Viability_001101_m)
  
  exp_111010_m <- (Pr_Y_bearing_sperm) * (Pr_111010_egg) * (Viability_111010_m)
  exp_000101_m <- (Pr_Y_bearing_sperm) * (Pr_000101_egg) * (Viability_000101_m)
  
  # Probability of observing male with quadruple exchange recombination pattern
  
  exp_010100_m <- (Pr_Y_bearing_sperm) * (Pr_010100_egg) * (Viability_010100_m)
  exp_101011_m <- (Pr_Y_bearing_sperm) * (Pr_101011_egg) * (Viability_101011_m)
  
  exp_010110_m <- (Pr_Y_bearing_sperm) * (Pr_010110_egg) * (Viability_010110_m)
  exp_101001_m <- (Pr_Y_bearing_sperm) * (Pr_101001_egg) * (Viability_101001_m)
  
  exp_010010_m <- (Pr_Y_bearing_sperm) * (Pr_010010_egg) * (Viability_010010_m)
  exp_101101_m <- (Pr_Y_bearing_sperm) * (Pr_101101_egg) * (Viability_101101_m)
  
  exp_011010_m <- (Pr_Y_bearing_sperm) * (Pr_011010_egg) * (Viability_011010_m)
  exp_100101_m <- (Pr_Y_bearing_sperm) * (Pr_100101_egg) * (Viability_100101_m)
  
  exp_001010_m <- (Pr_Y_bearing_sperm) * (Pr_001010_egg) * (Viability_001010_m)
  exp_110101_m <- (Pr_Y_bearing_sperm) * (Pr_110101_egg) * (Viability_110101_m)
  
  # Probability of observing male with quintuple exchange recombination pattern
  
  exp_101010_m <- (Pr_Y_bearing_sperm) * (Pr_101010_egg) * (Viability_101010_m)
  exp_010101_m <- (Pr_Y_bearing_sperm) * (Pr_010101_egg) * (Viability_010101_m)
  
  # Cumulative probability of observing males from higher order exchanges (>2)
  # Optional pooling of high order exchange classes (see rec in Zhao et al 1995)
  
  exp_high_order_exchange_m <- sum(exp_101000_m,
                                   exp_010111_m,
                                   exp_101100_m,
                                   exp_010011_m,
                                   exp_101110_m,
                                   exp_010001_m,
                                   exp_100100_m,
                                   exp_011011_m,
                                   exp_100110_m,
                                   exp_011001_m,
                                   exp_100010_m,
                                   exp_011101_m,
                                   exp_110100_m,
                                   exp_001011_m,
                                   exp_110110_m,
                                   exp_001001_m,
                                   exp_110010_m,
                                   exp_001101_m,
                                   exp_111010_m,
                                   exp_000101_m,
                                   exp_010100_m,
                                   exp_101011_m,
                                   exp_010110_m,
                                   exp_101001_m,
                                   exp_010010_m,
                                   exp_101101_m,
                                   exp_011010_m,
                                   exp_100101_m,
                                   exp_001010_m,
                                   exp_110101_m,
                                   exp_101010_m,
                                   exp_010101_m)
  
  # Probability of a female zygote failing to develop to an observable adult fly
  
  exp_failed_development_m <- sum((Pr_Y_bearing_sperm) * (Pr_000000_egg) * (1-Viability_000000_m),
                                  (Pr_Y_bearing_sperm) * (Pr_111111_egg) * (1-Viability_111111_m),
                                  (Pr_Y_bearing_sperm) * (Pr_100000_egg) * (1-Viability_100000_m),
                                  (Pr_Y_bearing_sperm) * (Pr_011111_egg) * (1-Viability_011111_m),
                                  (Pr_Y_bearing_sperm) * (Pr_110000_egg) * (1-Viability_110000_m),
                                  (Pr_Y_bearing_sperm) * (Pr_001111_egg) * (1-Viability_001111_m),
                                  (Pr_Y_bearing_sperm) * (Pr_111000_egg) * (1-Viability_111000_m),
                                  (Pr_Y_bearing_sperm) * (Pr_000111_egg) * (1-Viability_000111_m),
                                  (Pr_Y_bearing_sperm) * (Pr_111100_egg) * (1-Viability_111100_m),
                                  (Pr_Y_bearing_sperm) * (Pr_000011_egg) * (1-Viability_000011_m),
                                  (Pr_Y_bearing_sperm) * (Pr_111110_egg) * (1-Viability_111110_m),
                                  (Pr_Y_bearing_sperm) * (Pr_000001_egg) * (1-Viability_000001_m),
                                  (Pr_Y_bearing_sperm) * (Pr_010000_egg) * (1-Viability_010000_m),
                                  (Pr_Y_bearing_sperm) * (Pr_101111_egg) * (1-Viability_101111_m),
                                  (Pr_Y_bearing_sperm) * (Pr_011000_egg) * (1-Viability_011000_m),
                                  (Pr_Y_bearing_sperm) * (Pr_100111_egg) * (1-Viability_100111_m),
                                  (Pr_Y_bearing_sperm) * (Pr_011100_egg) * (1-Viability_011100_m),
                                  (Pr_Y_bearing_sperm) * (Pr_100011_egg) * (1-Viability_100011_m),
                                  (Pr_Y_bearing_sperm) * (Pr_011110_egg) * (1-Viability_011110_m),
                                  (Pr_Y_bearing_sperm) * (Pr_100001_egg) * (1-Viability_100001_m),
                                  (Pr_Y_bearing_sperm) * (Pr_001000_egg) * (1-Viability_001000_m),
                                  (Pr_Y_bearing_sperm) * (Pr_110111_egg) * (1-Viability_110111_m),
                                  (Pr_Y_bearing_sperm) * (Pr_001100_egg) * (1-Viability_001100_m),
                                  (Pr_Y_bearing_sperm) * (Pr_110011_egg) * (1-Viability_110011_m),
                                  (Pr_Y_bearing_sperm) * (Pr_001110_egg) * (1-Viability_001110_m),
                                  (Pr_Y_bearing_sperm) * (Pr_110001_egg) * (1-Viability_110001_m),
                                  (Pr_Y_bearing_sperm) * (Pr_000100_egg) * (1-Viability_000100_m),
                                  (Pr_Y_bearing_sperm) * (Pr_111011_egg) * (1-Viability_111011_m),
                                  (Pr_Y_bearing_sperm) * (Pr_000110_egg) * (1-Viability_000110_m),
                                  (Pr_Y_bearing_sperm) * (Pr_111001_egg) * (1-Viability_111001_m),
                                  (Pr_Y_bearing_sperm) * (Pr_000010_egg) * (1-Viability_000010_m),
                                  (Pr_Y_bearing_sperm) * (Pr_111101_egg) * (1-Viability_111101_m),
                                  (Pr_Y_bearing_sperm) * (Pr_101000_egg) * (1-Viability_101000_m),
                                  (Pr_Y_bearing_sperm) * (Pr_010111_egg) * (1-Viability_010111_m),
                                  (Pr_Y_bearing_sperm) * (Pr_101100_egg) * (1-Viability_101100_m),
                                  (Pr_Y_bearing_sperm) * (Pr_010011_egg) * (1-Viability_010011_m),
                                  (Pr_Y_bearing_sperm) * (Pr_101110_egg) * (1-Viability_101110_m),
                                  (Pr_Y_bearing_sperm) * (Pr_010001_egg) * (1-Viability_010001_m),
                                  (Pr_Y_bearing_sperm) * (Pr_100100_egg) * (1-Viability_100100_m),
                                  (Pr_Y_bearing_sperm) * (Pr_011011_egg) * (1-Viability_011011_m),
                                  (Pr_Y_bearing_sperm) * (Pr_100110_egg) * (1-Viability_100110_m),
                                  (Pr_Y_bearing_sperm) * (Pr_011001_egg) * (1-Viability_011001_m),
                                  (Pr_Y_bearing_sperm) * (Pr_100010_egg) * (1-Viability_100010_m),
                                  (Pr_Y_bearing_sperm) * (Pr_011101_egg) * (1-Viability_011101_m),
                                  (Pr_Y_bearing_sperm) * (Pr_110100_egg) * (1-Viability_110100_m),
                                  (Pr_Y_bearing_sperm) * (Pr_001011_egg) * (1-Viability_001011_m),
                                  (Pr_Y_bearing_sperm) * (Pr_110110_egg) * (1-Viability_110110_m),
                                  (Pr_Y_bearing_sperm) * (Pr_001001_egg) * (1-Viability_001001_m),
                                  (Pr_Y_bearing_sperm) * (Pr_110010_egg) * (1-Viability_110010_m),
                                  (Pr_Y_bearing_sperm) * (Pr_001101_egg) * (1-Viability_001101_m),
                                  (Pr_Y_bearing_sperm) * (Pr_111010_egg) * (1-Viability_111010_m),
                                  (Pr_Y_bearing_sperm) * (Pr_000101_egg) * (1-Viability_000101_m),
                                  (Pr_Y_bearing_sperm) * (Pr_010100_egg) * (1-Viability_010100_m),
                                  (Pr_Y_bearing_sperm) * (Pr_101011_egg) * (1-Viability_101011_m),
                                  (Pr_Y_bearing_sperm) * (Pr_010110_egg) * (1-Viability_010110_m),
                                  (Pr_Y_bearing_sperm) * (Pr_101001_egg) * (1-Viability_101001_m),
                                  (Pr_Y_bearing_sperm) * (Pr_010010_egg) * (1-Viability_010010_m),
                                  (Pr_Y_bearing_sperm) * (Pr_101101_egg) * (1-Viability_101101_m),
                                  (Pr_Y_bearing_sperm) * (Pr_011010_egg) * (1-Viability_011010_m),
                                  (Pr_Y_bearing_sperm) * (Pr_100101_egg) * (1-Viability_100101_m),
                                  (Pr_Y_bearing_sperm) * (Pr_001010_egg) * (1-Viability_001010_m),
                                  (Pr_Y_bearing_sperm) * (Pr_110101_egg) * (1-Viability_110101_m),
                                  (Pr_Y_bearing_sperm) * (Pr_101010_egg) * (1-Viability_101010_m),
                                  (Pr_Y_bearing_sperm) * (Pr_010101_egg) * (1-Viability_010101_m))
  
  # Vector of expected frequencies of each all sex-specific phenotypic classes
  
  expected_freqs <- c(exp_failed_development_f,
                      exp_000000_f,
                      exp_111111_f,
                      exp_100000_f,
                      exp_011111_f,
                      exp_110000_f,
                      exp_001111_f,
                      exp_111000_f,
                      exp_000111_f,
                      exp_111100_f,
                      exp_000011_f,
                      exp_111110_f,
                      exp_000001_f,
                      exp_010000_f,
                      exp_101111_f,
                      exp_011000_f,
                      exp_100111_f,
                      exp_011100_f,
                      exp_100011_f,
                      exp_011110_f,
                      exp_100001_f,
                      exp_001000_f,
                      exp_110111_f,
                      exp_001100_f,
                      exp_110011_f,
                      exp_001110_f,
                      exp_110001_f,
                      exp_000100_f,
                      exp_111011_f,
                      exp_000110_f,
                      exp_111001_f,
                      exp_000010_f,
                      exp_111101_f,
                      exp_101000_f,
                      exp_010111_f,
                      exp_101100_f,
                      exp_010011_f,
                      exp_101110_f,
                      exp_010001_f,
                      exp_100100_f,
                      exp_011011_f,
                      exp_100110_f,
                      exp_011001_f,
                      exp_100010_f,
                      exp_011101_f,
                      exp_110100_f,
                      exp_001011_f,
                      exp_110110_f, 
                      exp_001001_f,
                      exp_110010_f,
                      exp_001101_f,
                      exp_111010_f,
                      exp_000101_f,
                      exp_010100_f,
                      exp_101011_f,
                      exp_010110_f,
                      exp_101001_f,
                      exp_010010_f,
                      exp_101101_f,
                      exp_011010_f, 
                      exp_100101_f,
                      exp_001010_f, 
                      exp_110101_f,
                      exp_101010_f, 
                      exp_010101_f,
                      exp_failed_development_m,
                      exp_000000_m,
                      exp_111111_m,
                      exp_100000_m,
                      exp_011111_m,
                      exp_110000_m,
                      exp_001111_m,
                      exp_111000_m,
                      exp_000111_m,
                      exp_111100_m,
                      exp_000011_m,
                      exp_111110_m,
                      exp_000001_m,
                      exp_010000_m,
                      exp_101111_m,
                      exp_011000_m,
                      exp_100111_m,
                      exp_011100_m,
                      exp_100011_m,
                      exp_011110_m,
                      exp_100001_m,
                      exp_001000_m,
                      exp_110111_m,
                      exp_001100_m,
                      exp_110011_m,
                      exp_001110_m,
                      exp_110001_m,
                      exp_000100_m,
                      exp_111011_m,
                      exp_000110_m,
                      exp_111001_m,
                      exp_000010_m,
                      exp_111101_m,
                      exp_101000_m,
                      exp_010111_m,
                      exp_101100_m,
                      exp_010011_m,
                      exp_101110_m,
                      exp_010001_m,
                      exp_100100_m,
                      exp_011011_m,
                      exp_100110_m,
                      exp_011001_m,
                      exp_100010_m,
                      exp_011101_m,
                      exp_110100_m,
                      exp_001011_m,
                      exp_110110_m,
                      exp_001001_m,
                      exp_110010_m,
                      exp_001101_m,
                      exp_111010_m,
                      exp_000101_m,
                      exp_010100_m,
                      exp_101011_m,
                      exp_010110_m,
                      exp_101001_m,
                      exp_010010_m,
                      exp_101101_m,
                      exp_011010_m,
                      exp_100101_m,
                      exp_001010_m,
                      exp_110101_m,
                      exp_101010_m,
                      exp_010101_m)
  
  # raw_log_likelihoods for each of the individual recombination classes
  # log_likelihood is the sum of the individual recombination classes
  # Function returns negative log-likelihood value because optimization
  # with the Nelder-Mead method in optim() is, by default, minimization 
  raw_log_likelihoods <- observed_count*log(expected_freqs)
  log_likelihood <- sum(raw_log_likelihoods)
  return(-log_likelihood)
}

# Create matrix to store all optimized outputs of Uniform.CxCoM.Likelihood
# Log likelihoods calculated for each combination of y & p model parameters
MLE_Uniform_CxCoM <- matrix(nrow=0, ncol=26)
colnames(MLE_Uniform_CxCoM)<-c("MLE_X_length",
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
                               "MLE_v1f-",
                               "MLE_v2f-",
                               "MLE_v3f-",
                               "MLE_v4f-",
                               "MLE_v5f-",
                               "MLE_v6f-",
                               "MLE_v1m-",
                               "MLE_v2m-",
                               "MLE_v3m-",
                               "MLE_v4m-",
                               "MLE_v5m-",
                               "MLE_v6m-",
                               "lnL",
                               "Convergence")

# For loop Uniform.CxCoM.Likelihood for all possible combinations of p
# Output model parameters and likelihoods in matrix MLE_Uniform_CxCoM
for (v in 1:max_p_value) {
  
  # Define parameters starting points for y, v, and p in single run of Nelder-Mead algorithm
  # When bounding Nelder-Mead using package dfoptim, starting points cannot be at the boundaries
  parameters <- as.numeric(c(y_uniform[v,1],
                             y_uniform[v,2],
                             y_uniform[v,3],
                             y_uniform[v,4],
                             y_uniform[v,5],
                             0.5,
                             0.5,
                             0.5,
                             0.5,
                             0.5,
                             0.5,
                             0.5,
                             0.5,
                             0.5,
                             0.5,
                             0.5,
                             0.5))
  p <- p_uniform[v,]
  
  # Optimize CxCoM model parameters y and v with bounded Nelder-Mead method
  # Upper limit on y's is arbitary, whereas upper and lower limit on v's is
  # given by definition of egg-to-adult viablity (cannot excede 0-1 range)
  Uniform_CxCoM_Table <- nmkb(par = parameters,
                              fn = Uniform.CxCoM.Likelihood,
                              lower = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                              upper = c(9,9,9,9,9,1,1,1,1,1,1,1,1,1,1,1,1))
  
  # Create output vector of y, p, likelihood, convergence code
  MLE_outputs <- c((((Uniform_CxCoM_Table$par[1]+Uniform_CxCoM_Table$par[2]+Uniform_CxCoM_Table$par[3]+Uniform_CxCoM_Table$par[4]+Uniform_CxCoM_Table$par[5])/(2*p[1]))*100),
                   Uniform_CxCoM_Table$par[1]/(2*p[1]),
                   Uniform_CxCoM_Table$par[2]/(2*p[1]),
                   Uniform_CxCoM_Table$par[3]/(2*p[1]),
                   Uniform_CxCoM_Table$par[4]/(2*p[1]),
                   Uniform_CxCoM_Table$par[5]/(2*p[1]),
                   p[1],
                   Uniform_CxCoM_Table$par[1],
                   Uniform_CxCoM_Table$par[2],
                   Uniform_CxCoM_Table$par[3],
                   Uniform_CxCoM_Table$par[4],
                   Uniform_CxCoM_Table$par[5],
                   Uniform_CxCoM_Table$par[6],
                   Uniform_CxCoM_Table$par[7],
                   Uniform_CxCoM_Table$par[8],
                   Uniform_CxCoM_Table$par[9],
                   Uniform_CxCoM_Table$par[10],
                   Uniform_CxCoM_Table$par[11],
                   Uniform_CxCoM_Table$par[12],
                   Uniform_CxCoM_Table$par[13],
                   Uniform_CxCoM_Table$par[14],
                   Uniform_CxCoM_Table$par[15],
                   Uniform_CxCoM_Table$par[16],
                   Uniform_CxCoM_Table$par[17],
                   Uniform_CxCoM_Table$value,
                   Uniform_CxCoM_Table$convergence)
  
  # Add vector of y, p, likelihood, code to MLE_Uniform_CxCoM
  MLE_Uniform_CxCoM <- rbind(MLE_Uniform_CxCoM, MLE_outputs)
}

# Take output matrix MLE_Uniform_CxCoM and sort by likelihood
# Print the sorted table for records, make sure to RENAME FILE
output<-as.data.frame(MLE_Uniform_CxCoM)
df2 <- output[order(output$lnL),]
write.csv(df2,"/Users/spencerkoury/Desktop/SK14_H2_full_output.csv",row.names = TRUE)


 