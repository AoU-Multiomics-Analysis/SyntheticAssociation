# template
# Overview: 
This pipeline takes finemapping results at different thresholds (typically should be MAF > 1% and MAF > 0.1% or lower) and first gets the set of variants that are uniquely identified from each analysis. 
Then using the sets of variants it performs a conditional regression analysis. The goal of this is to find variants that have forms of synthetic associations. 


## Inputs: 
`GenotypeMatrix`: Genotype dosages for the cis window generated from BCFtools for the cis window. This is generated in the susieR workflow 
`GenotypeMatrixIndex`: Index file for the above file 
`QTLCovariates`: QTL covariates used in tensorQTL and susie 
`ThresholdFinemapping`: parquet file for fine mapping of a single gene that is at a higher MAF threshold then the full analysis 
`FullFinemapping`: TSV file that has annotated allele frequencies for all variants of all genes. This is generated from the Aggregate and Annotate susie pipelines. This is the standard fine-mapping analysis file 
`PhenotypeBed`: Bed file for molecular data, same input as tensorQTL 

## Outputs:
