version 1.0 


task ConditionalRegression {
    input {
        File GenotypeMatrix
        File GenotypeMatrixIndex
        File QTLCovariates
        String OutputPrefix
        File ThresholdFinemapping
        File FullFinemapping
        File PhenotypeBed
        Int Memory 
        Int NumPrempt
    }
    
    command <<<
    Rscript /SyntheticAssociations.R \
        --genotype_matrix ~{GenotypeMatrix} \
        --covariates ~{QTLCovariates} \
        --out_prefix ~{OutputPrefix} \
        --thresholded_finemapping ~{ThresholdFinemapping} \
        --full_finemapping  ~{FullFinemapping} \
        --phenotype_bed ~{PhenotypeBed}
    >>>
    
    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/SyntheticAssociation:main"
        memory: "~{Memory}GB"
        disk: "local-disk 500 SSD"
        bootDiskSizeGb: 25
        preemptible: "${NumPrempt}"
        cpu: "1"
    }

    output {
        File ConditionalRes = "~{OutputPrefix}_conditional_analysis.tsv"
    }
}

workflow SyntheticAssociations {
    input {
        File GenotypeMatrix
        File GenotypeMatrixIndex
        File QTLCovariates
        String OutputPrefix
        File ThresholdFinemapping
        File FullFinemapping
        File PhenotypeBed
        Int Memory 
        Int NumPrempt
    }

    call ConditionalRegression {
        input:
            GenotypeMatrix = GenotypeMatrix,
            GenotypeMatrixIndex = GenotypeMatrixIndex,
            QTLCovariates = QTLCovariates,
            OutputPrefix = OutputPrefix,
            ThresholdFinemapping = ThresholdFinemapping,
            FullFinemapping = FullFinemapping,
            PhenotypeBed = PhenotypeBed,
            Memory = Memory,
            NumPrempt = NumPrempt
    }

    output {
        File ConditionalRes = ConditionalRegression.ConditionalRes
    }
}
