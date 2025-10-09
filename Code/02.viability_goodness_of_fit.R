
# enter egg count sample size (n)
n <- 9995
  
# enter observed and unobserved counts
observed_female_wildtype_adult_count <-
observed_male_wildtype_adult_count   <-
observed_female_mutant_adult_count   <-
observed_male_mutant_adult_count     <-
observed_failed_development_count    <-

# solve nonlinear system of equations to estimate viability effects
nls(observed_counts_vector~H0)
  
  

# save nls() outputs as viability effects for goodness-of-fit tests
H1_v_1 <- 
H2_v_1 <- 
H3_v_1 <- 
H4_v_1 <- 
H5_v_1 <- 
H2_v_2 <- 
H3_v_3 <- 
H3_v_4 <- 
H4_v_5 <-
H4_v_6 <-
H4_v_7 <-
H5_v_8 <-
H5_v_9 <-
H5_v_10<-
H5_v_11<-

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
H2_expected_failed_development    <- (n*(1-H2_v_1)*(2-H2_v_2))/2

# H3: experimental mortality has both a random component and sex-specific marker associated effects
H3_expected_female_wildtype_count <- (n*H3_v_1)/4
H3_expected_male_wildtype_count   <- (n*H3_v_1)/4
H3_expected_female_mutant_count   <- (n*H3_v_1*H3_v_3)/4
H3_expected_male_mutant_count     <- (n*H3_v_1*H3_v_4)/4
H3_expected_failed_development    <- (n*(1-H3_v_1)*(4-H3_v_3-H3_v_4))/4

# H4: experimental mortality has both a random compoment as well as sex and sex-specific marker associated effects
H4_expected_female_wildtype_count <- (n*H4_v_1)/4
H4_expected_male_wildtype_count   <- (n*H4_v_1*H4_v_5)/4
H4_expected_female_mutant_count   <- (n*H4_v_1*H4_v_6)/4
H4_expected_male_mutant_count     <- (n*H4_v_1*H4_v_7)/4
H4_expected_failed_development    <- (n*(1-H4_v_1)*(4-H4_v_5-H4_v_6-H4_v_7))/4

# H5: experimental mortality has both a random component and unique deterministic component for each adult phenotypic class
H5_expected_female_wildtype_count <- (n*H5_v_1*H5_v_8)/4
H5_expected_male_wildtype_count   <- (n*H5_v_1*H5_v_9)/4
H5_expected_female_mutant_count   <- (n*H5_v_1*H5_v_10)/4
H5_expected_male_mutant_count     <- (n*H5_v_1*H5_v_11)/4
H5_expected_failed_development    <- (n*(1-H5_v_1)*(4-H5_v_8-H5_v_9-H5_v_10-H5_v_11))/4

# Goodness-of-fit G-test for each hypothesis/viability model using parameters fit by method of least squares with nls()
G_statistic_H0 <- 2*sum(observed_female_wildtype_adult_count*(ln((observed_female_wildtype_adult_count/H0_expected_female_wildtype_count))),
                   observed_male_wildtype_adult_count*(ln((observed_male_wildtype_adult_count/H0_expected_male_wildtype_count))),
                   observed_female_mutant_adult_count*(ln((observed_female_mutant_adult_count/H0_expected_female_mutant_count))),
                   observed_male_mutant_adult_count*(ln((observed_male_mutant_adult_count/H0_expected_male_mutant_count))))

G_statistic_H1 <- 2*sum(observed_female_wildtype_adult_count*(ln((observed_female_wildtype_adult_count/H0_expected_female_wildtype_count))),
                   observed_male_wildtype_adult_count*(ln((observed_male_wildtype_adult_count/H0_expected_male_wildtype_count))),
                   observed_female_mutant_adult_count*(ln((observed_female_mutant_adult_count/H0_expected_female_mutant_count))),
                   observed_male_mutant_adult_count*(ln((observed_male_mutant_adult_count/H0_expected_male_mutant_count))),
                   observed_failed_development_count*(ln((observed_failed_development_count/H0_expected_failed_development))))

G_statistic_H2 <- 2*sum(observed_female_wildtype_adult_count*(ln((observed_female_wildtype_adult_count/H0_expected_female_wildtype_count))),
                   observed_male_wildtype_adult_count*(ln((observed_male_wildtype_adult_count/H0_expected_male_wildtype_count))),
                   observed_female_mutant_adult_count*(ln((observed_female_mutant_adult_count/H0_expected_female_mutant_count))),
                   observed_male_mutant_adult_count*(ln((observed_male_mutant_adult_count/H0_expected_male_mutant_count))),
                   observed_failed_development_count*(ln((observed_failed_development_count/H0_expected_failed_development))))

G_statistic_H3 <- 2*sum(observed_female_wildtype_adult_count*(ln((observed_female_wildtype_adult_count/H0_expected_female_wildtype_count))),
                   observed_male_wildtype_adult_count*(ln((observed_male_wildtype_adult_count/H0_expected_male_wildtype_count))),
                   observed_female_mutant_adult_count*(ln((observed_female_mutant_adult_count/H0_expected_female_mutant_count))),
                   observed_male_mutant_adult_count*(ln((observed_male_mutant_adult_count/H0_expected_male_mutant_count))),
                   observed_failed_development_count*(ln((observed_failed_development_count/H0_expected_failed_development))))

G_statistic_H4 <- 2*sum(observed_female_wildtype_adult_count*(ln((observed_female_wildtype_adult_count/H0_expected_female_wildtype_count))),
                   observed_male_wildtype_adult_count*(ln((observed_male_wildtype_adult_count/H0_expected_male_wildtype_count))),
                   observed_female_mutant_adult_count*(ln((observed_female_mutant_adult_count/H0_expected_female_mutant_count))),
                   observed_male_mutant_adult_count*(ln((observed_male_mutant_adult_count/H0_expected_male_mutant_count))),
                   observed_failed_development_count*(ln((observed_failed_development_count/H0_expected_failed_development))))

G_statistic_H5 <- 2*sum(observed_female_wildtype_adult_count*(ln((observed_female_wildtype_adult_count/H0_expected_female_wildtype_count))),
                   observed_male_wildtype_adult_count*(ln((observed_male_wildtype_adult_count/H0_expected_male_wildtype_count))),
                   observed_female_mutant_adult_count*(ln((observed_female_mutant_adult_count/H0_expected_female_mutant_count))),
                   observed_male_mutant_adult_count*(ln((observed_male_mutant_adult_count/H0_expected_male_mutant_count))),
                   observed_failed_development_count*(ln((observed_failed_development_count/H0_expected_failed_development)))) 
  
  
  
  
  
  
  

