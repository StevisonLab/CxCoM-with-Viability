# Statistical analysis of meiotic events with viability effects

Classic recombination experiments designed to test genetic and environmental treatments do not directly measure crossing-over, instead rates and distribution of meiotic events in F<sub>1</sub> oocytes is inferred from genetic markers in F<sub>2</sub> adults. In _Drosophila melanogaster_ this procedure introduces substantial "missing data problems" because 75% of meiotic chromatids segregate to polar body nuclei and another 11% are transmitted to inviable F<sub>2</sub> zygotes which cannot be scored for recombination. To address these sources of uncertainty and bias we extend the _Cx(Co)<sup>m</sup>_ model of the data-generating process by assuming: 1) double strand breaks occur as a Poisson point process, 2) crossover maturation is a stationary renewal process, 3) chromosome segregation is random one-half thinning of this process, 4) fertilization by X- versus Y-bearing sperm is mendelian, and 5) egg-to-adult survival is binomially distributed with a rate parameter determined by F<sub>2</sub> marker alleles. To quantify experimental mortality, we performed egg counts in 6-point X chromosome testcrosses and marker-free controls on identical genetic backgrounds under standard laboratory conditions. Likelihood ratio tests using the 19,927 fly dataset reveals 44% F<sub>2</sub> experimental mortality, and support a model where 36 of the 44% is due to sex-specific viability effects. Variability in X chromosome genetic lengths with experimental mortality can be simulated and we provide case-control 80% power curves to guide experimental design. We propose that differential mortality should be the _de facto_ null hypothesis when comparing F<sub>2</sub> recombinant fractions and provide a probabilistic model of the data-generating process to improve characterization of patterns in F<sub>1</sub> meiotic events.

# Code Directory

To perform all R analyses in "An extension of the Cx(Co)<sup>m</sup> model of crossover patterning to account for experimental mortality in _Drosophila melanogaster_" there are six **raw data** input files in "Raw Data" folder. Please set this as your working directory:

There are eight **scripts** to run the analyses in "Code" folder but you only need to run ONE.

## Primary Script Files:

1. **1_missing\_data\_driver\_script.R** == THIS IS THE PRIMARY DRIVER SCRIPT!!! Others will be sourced in the main script and do not need to be run separately. This script also sets the working directory and loads the package dfoptim().
2. 2_missing\_data\_dataset\_compilation\_script.R
3. 3_missing\_data\_single\_locus\_custom\_function.R
4. 4_missing\_data\_multi\_locus\_custom\_function.R
5. 5_missing\_data\_organismal\_analysis.R
6. 6_missing\_data\_single\_locus\_analysis.R
7. 7_missing\_data\_multi\_locus\_analysis.R
8. 8_missing\_data\_pooling\_analysis.R

## Raw Data Files:

1. 1\_EH\_raw\_multiply\_marked\_egg\_count.csv
2. 2\_MS\_raw\_multiply\_marked\_egg\_count.csv
3. 3\_SK\_raw\_multiply\_marked\_adult\_phenotypic\_class\_counts.csv
4. A\_EH\_raw\_marker\_free\_egg\_counts.csv
5. B\_MS\_raw\_marker\_free\_egg\_counts.csv
6. C\_SK\_raw\_marker\_free\_adult\_counts.csv

## Output files:

The following output files will be written into the **"Summaries and Outputs"** folder:

- 4\_compiled\_raw\_multiply\_marked\_dataset.csv
- 5\_derived\_multiply\_marked\_dataset.csv
- 6\_multiply\_marked\_viability\_dataset.csv
- 7\_scute\_single\_locus\_dataset.csv
- 8\_crossveinless\_single\_locus\_dataset.csv
- 9\_vermilion\_single\_locus\_dataset.csv
- 10\_forked\_single\_locus\_dataset.csv
- 11\_carnation\_single\_locus\_dataset.csv
- 12\_yellow\_plus\_single\_locus\_dataset.csv
- 13\_multi\_locus\_individual\_vials\_dataset.csv
- 14\_multi\_locus\_cross\_pooled\_dataset.csv
- 15\_multi\_locus\_brood\_pooled\_dataset.csv
- 16\_multi\_locus\_full\_experiment\_dataset.csv
- 17\_multi\_locus\_full\_experiment\_H0\_mle\_output.csv
- 18\_multi\_locus\_full\_experiment\_H1\_mle\_output.csv
- 19\_multi\_locus\_full\_experiment\_H2\_mle\_output.csv
- 20\_multi\_locus\_full\_experiment\_H3\_mle\_output.csv
- 21\_multi\_locus\_individual\_vials\_H3\_mle\_output.csv
- 22\_multi\_locus\_cross\_pooled\_H3\_mle\_output.csv
- 23\_multi\_locus\_brood\_pooled\_H3\_mle\_output.csv
- D\_compiled\_raw\_marker\_free\_dataset.csv
- E\_derived\_marker\_free\_dataset.csv
- F\_marker\_free\_viability\_dataset.csv
- table\_S1\_marker\_free\_viability\_regression.csv
- table\_S2\_multiply\_marked\_viability\_regression.csv
- table\_S3\_scute\_single\_locus\_G\_tests.csv.csv
- table\_S4\_crossveinless\_single\_locus\_G\_tests.csv
- table\_S5\_vermilion\_single\_locus\_G\_tests.csv
- table\_S6\_forked\_single\_locus\_G\_tests.csv
- table\_S7\_carnation\_single\_locus\_G\_tests.csv
- table\_S8\_yellow\_plus\_single\_locus\_G\_tests.csv
- table\_S10\_vials\_pooled\_anova.csv
- table\_S11\_cross\_pooled\_anova.csv
- table\_S12\_brood\_pooled\_anova.csv
  
