library(tidyverse)
library(data.table)
library(magrittr)

###### PARSE COMMAND LINE ARGUMENTS ########## 
option_list <- list(
  # TODO look around if there is a package recognizing delimiter in dataset
  optparse::make_option(c("--genotype_matrix"), type="character", default=NULL,
    help="Genotype dosage matrix extracted from VCF", metavar="type"),
  optparse::make_option(c("--covariates"), type="character", default=NULL,
    help="Path to covariates file in QTLtools format", metavar="type"),
  optparse::make_option(c("--out_prefix"), type="character", default="./finemapping_output",
    help="Prefix of output files", metavar="type"),
  optparse::make_option(c("--thresholded_finemapping"), type="character", default="./finemapping_output",
    help="Path to parquet file for finemapping data that has been thresholded on allele frequency", metavar="type"),
  optparse::make_option(c("--full_finemapping"), type="character", default="./finemapping_output",
    help="Path to parquet file for finemapping data that has been run across all variants", metavar="type"),
  optparse::make_option(c("--phenotype_bed"), type="character", default="./finemapping_output",
    help="Path to bed file that contains the phenotype data", metavar="type")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

########## FUNCTIONS ##############
parse_genotype_data <- function(genotype_path) {
genotype_dosages_df <- fread(basename(genotype_path)) %>% 
    mutate(variant = paste(CHROM,POS,REF,ALT,sep = '_')) %>% 
    select(-CHROM,-POS,-REF,-ALT) %>% 
    select(variant,everything()) %>%
    column_to_rownames('variant')
  
genotype_dosages_df 
}


parse_expression_data <- function (expression_df, input_gene_id) {
    gene_vector <- expression_df %>% 
        filter(gene_id == input_gene_id) %>% 
        select(-1, -2, -3) %>% 
        column_to_rownames("gene_id") %>% 
        t() %>% 
        data.frame()
    gene_vector
}

parse_expression_covars <- function(covars_path) {
expression_covars_df <- fread(covars_path) %>%  janitor::row_to_names(row_number = 1)
rownames(expression_covars_df) <- NULL
expression_covars_df %<>% column_to_rownames('ID')
expression_covars_df %<>% t()    
expression_covars_df
}

extract_variants <- function(genotype_df,variant_list){
genotypes <- genotype_df[variant_list,] %>% t() %>% data.frame()
genotypes
}

merge_genotype_covariate_data <- function(covar_df,genotype_df) {
merged_data <- left_join(
                covar_df %>% data.frame() %>% rownames_to_column('ID'),
                genotype_df %>% data.frame() %>% rownames_to_column('ID'),
                by = 'ID'
                ) %>% 
                column_to_rownames('ID')   
merged_data  
}

create_variant_list <- function(variant_list1,variant_list2) {
out <- Map(list, rep(variant_list1, each = length(variant_list2)), 
           rep(variant_list2, times = length(variant_list1)))
out   
}

create_pairwise_variant_list <- function(x,y) {
#out <- map2(x, y, ~ list(x, y))
out <- map2(x, y, c)
out

}

run_conditional_regression <- function(genotype_dosages,
                                      expression_covars,
                                      gene_vector,
                                      variant_list) {
message('Extracting genotype data')
#variant_genotype_data <- extract_variants(genotype_dosages,variant_list)
variant_genotype_data <- genotype_dosages %>% select(all_of(variant_list))
merged_data <- merge_genotype_covariate_data(expression_covars,variant_genotype_data)

message('Computing correlation between variants')

if (length(variant_list) == 2) {
variant_cor <- cor(variant_genotype_data) %>% 
    data.frame() %>% 
    rownames_to_column('variant') %>% 
    pivot_longer(!variant) %>% 
    filter(variant != name) %>% 
    distinct(value) %>% 
    pull(value)    
} else {
variant_cor <- NA 
    
}
    
message('Running conditional regression')
conditional_regression_result <-  broom::tidy(lm(gene_vector[,1] ~ .,data = merged_data )) %>% 
            filter(str_detect(term,'chr')) %>% 
            mutate(type = 'conditional',
                   variants = paste0(colnames(variant_genotype_data),collapse=';'),
                   LD = variant_cor)  %>% 
            mutate(variants = str_remove(variants,term)) %>% 
            mutate(variants = str_remove(variants,'^;|;$')) %>% 
            dplyr::rename('conditional_variant' = 'variants')

    

message('Running marginal regressions ')
marginal_regression_results <- data.frame()
for (x in (colnames(variant_genotype_data))) {
    subset_genotype_data <- variant_genotype_data %>% select(x)
    subset_merged_data <- merge_genotype_covariate_data(expression_covars_df,subset_genotype_data)
    marginal_res <-  broom::tidy(lm(gene_vector[,1] ~ .,data = subset_merged_data )) %>% 
            filter(str_detect(term,'chr')) %>% 
            mutate(type = 'marginal')  
    marginal_regression_results <- bind_rows(marginal_res,marginal_regression_results)

} 
    
output <- bind_rows(conditional_regression_result,marginal_regression_results)
output
}

create_variant_df <- function(variant_sets) {
variant_df <- unlist(variant_sets) %>% 
    unique() %>% 
    data.frame() %>% 
    dplyr::rename('variant' = 1) %>% 
    separate(variant,into = c('chrom','pos','ref','alt'),remove= FALSE) %>% 
    mutate(pos = as.numeric(pos))   
variant_df   
}


extract_genotype_data_tabix <- function(variant_sets,genotype_dosages_path) {
variant_df <- create_variant_df(variant_sets)
start_pos <- min(variant_df$pos)
end_pos <- max(variant_df$pos)
chrom <- variant_df %>% distinct(chrom) %>% pull(chrom)
genotype_dat <- eQTLUtils::extractGenotypeMatrixFromDosage(chrom, 
                                           start_pos, 
                                           end_pos,
                                           genotype_dosages
                                           ) %>% 
                    data.frame() %>% 
                    rownames_to_column('variant') %>% 
                    mutate(variant = str_replace(variant,'chrchr','chr')) %>%
                    column_to_rownames('variant') %>% 
                    t() %>% 
                    data.frame() %>% 
                    select(all_of(variant_df$variant))
rownames(genotype_dat) <- str_remove(rownames(genotype_dat),'X')
genotype_dat   
    
}


create_common_with_all_rare <- function(rare_variants, common_variants) {
  Map(
    function(common) list( common, rare_variants),
    common_variants
  )
}

create_variant_list <- function(variant_list1,variant_list2) {
out <- Map(list, rep(variant_list1, each = length(variant_list2)), 
           rep(variant_list2, times = length(variant_list1)))
out   
}
compare_variant_sets <- function(all_variant_finemapping,thresholded_finemapping,input_gene_id) {
rare_variants <- all_variant_finemapping %>% 
                    filter(molecular_trait_id == input_gene_id)  %>%
                    filter(!variant %in% thresholded_finemapping$variant) %>% 
                    filter(MAF < .01 & pip > 0.9)%>% 
                    pull(variant)
common_lead_cs_variants <- thresholded_finemapping %>% 
                    filter(molecular_trait_id == input_gene_id) %>% 
                    group_by(cs_id) %>% 
                    filter(pip == max(pip)) %>% 
                    pull(variant)
variant_list<- create_variant_list(rare_variants,common_lead_cs_variants)
common_with_all_rare_list <- create_common_with_all_rare(rare_variants,common_lead_cs_variants)
conditional_variant_sets <- c(variant_list,common_with_all_rare_list)
#outlist <- list(common_variants = common_lead_cs_variants,conditional_sets = conditional_variant_sets)
outlist <- c(conditional_variant_sets,common_with_all_rare_list)
outlist
}



#run_conditional_analysis_pipeline <- function(gene_id,
                                              #full_finemapped_data,
                                              #thresholded_fm_data,
                                              #genotype_dosages_path,
                                              #expression_df,
                                              #covars_df,
                                              #out_dir ='conditonal_analysis'
                                              #) {
    #system(paste0('mkdir -p ',out_dir))
    #outfile_name <- paste0(gene_id,'_conditional_analysis.tsv')
    #full_outfile <- paste0(out_dir,'/',outfile_name)
    
    #message('Extracting variant sets')
    #variant_sets <- compare_variant_sets(full_finemapped_data,thresholded_fm_data,gene_id)
    
    #message('Extracting genotype dosages')
    #variant_list_unique <- variant_sets %>% flatten() %>% unlist() %>% unique()
    #genotype_dosages_df <- extract_genotype_data_tabix(variant_list_unique,genotype_dosages)
    ##genotype_dosages_df <- extract_genotype_data_tabix(variant_sets,genotype_dosages)
    
    #message('Extracting molecular trait data')
    #gene_vector <- parse_expression_data(expression_df,gene_id)  
    
    #message('Running conditional analysis')
    #conditional_analysis <- variant_sets %>% 
        #map_dfr(~run_conditional_regression(genotype_dosages_df,
                                    #covars_df,
                                    #gene_vector,
                                    #unlist(.))) %>% 
        #distinct() %>% 
        #mutate(gene_id = gene_id)
    
    #conditional_analysis %>% write_tsv(full_outfile)
    #conditional_analysis
    #}

run_conditional_analysis_pipeline <- function(gene_id,
                                              full_finemapped_data,
                                              thresholded_fm_data,
                                              genotype_dosages_path,
                                              expression_df,
                                              covars_df,
                                              out_dir ='conditonal_analysis',
                                              write_data = TRUE
                                              ) {

    
    message('Extracting variant sets')
    variant_sets <- compare_variant_sets(full_finemapped_data,thresholded_fm_data,gene_id)
    
    message('Extracting genotype dosages')
    genotype_dosages_df <- extract_genotype_data_tabix(variant_sets,genotype_dosages)
    
    message('Extracting molecular trait data')
    gene_vector <- parse_expression_data(expression_df,gene_id)  
    
    message('Running conditional analysis')
    conditional_analysis <- variant_sets %>% 
        map_dfr(~run_conditional_regression(genotype_dosages_df,
                                    covars_df,
                                    gene_vector,
                                    unlist(.))) %>% 
        distinct() %>% 
        mutate(gene_id = gene_id)
    if (write_data == TRUE) {
    conditional_analysis %>% write_tsv(full_outfile)
    system(paste0('mkdir -p ',out_dir))
    outfile_name <- paste0(gene_id,'_conditional_analysis.tsv')
    full_outfile <- paste0(out_dir,'/',outfile_name)
        }
    conditional_analysis
        
    }

########## PARSE ARGUMENTS  ###########

GenotypeDosages <- opt$genotype_matrix 
CovarsPath <- opt$covariates
PhenotypeBed <- opt$phenotype_bed
FullFinemapping <- opt$full_finemapping
OutputPrefix <- opt$out_prefix
ThresholdFinemapping <- opt$thresholded_finemapping


OutFile <- paste0(OutputPrefix,'_conditional_analysis.tsv')
########### LOAD DATA #########

message('Loading threshold data')
ThresholdData <- arrow::read_parquet(ThresholdFinemapping) %>%
    mutate(variant = str_replace(variant,'chrchr','chr'))
GeneList <- ThresholdData %>% distinct(molecular_trait_id) %>% dplyr::rename('gene_id'  = 1) 

message('Loading covariates')
CovarsDf <- parse_expression_covars(CovarsPath)

message('Loading molecular data')
PhenotypeBedDf <- fread(PhenotypeBed)

message('Loading full finemapping')
FullFinemapping <- fread(FullFinemapping)

######### RUN ANALYSIS #############
ConditionalAnalysisRes <- GeneList %>% 
    pull(gene_id) %>% 
    map_dfr(~run_conditional_analysis_pipeline(
                    .,
                    FullFinemapping,
                    ThresholdData,
                    GenotypeDosages,
                    PhenotypeBedDf,
                    CovarsDf,
                    write_data = FALSE
            ))

message('Writing results to output')
ConditionalAnalysisRes %>% write_tsv(OutFile)
