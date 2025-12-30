# Code Directory

To perform all R analyses in "An extension of the Cx(Co)^m model of crossover patterning to account for experimental mortality in \_Drosophila melanogaster\_" there are six **raw data** input files in "Raw Data" folder. Please set this as your working directory:

## Raw Data Files:

1. 1\_EH\_raw\_multiply\_marked\_egg\_count.csv
2. 2\_MS\_raw\_multiply\_marked\_egg\_count.csv
3. 3\_SK\_raw\_multiply\_marked\_adult\_phenotypic\_class\_counts.csv
4. A\_EH\_raw\_marker\_free\_egg\_counts.csv
5. B\_MS\_raw\_marker\_free\_egg\_counts.csv
6. C\_SK\_raw\_marker\_free\_adult\_counts.csv


There are eight **scripts** to run the analyses in "Code" folder

## Primary Script files:

1. 1_missing\_data\_driver\_script.R
2. 2_missing\_data\_dataset\_compilation\_script.R
3. 3_missing\_data\_single\_locus\_custom\_function.R
4. 4_missing\_data\_multi\_locus\_custom\_function.R
5. 5_missing\_data\_organismal\_analysis.R
6. 6_missing\_data\_single\_locus\_analysis.R
7. 7_missing\_data\_multi\_locus\_analysis.R
8. 8_missing\_data\_pooling\_analysis.R


1_missing\_data\_driver\_script.R sources each of these scripts, as well as loading the package dfoptim()

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














