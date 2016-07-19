#!/bin/bash

FLS="/Users/simondi/PHD/data/data/sysvirdrug/integrated_data_files/rnai.screen.offtc.tsv"

OUT="/Users/simondi/PHD/data/data/sysvirdrug/genes"
CHIKV_OUT="${OUT}/chikv_kinase_genes.tsv"
DENV_DUG_OUT="${OUT}/denv_druggable_genes.tsv"
DENV_KINA_OUT="${OUT}/denv_kinase_genes.tsv"
HCV_DUG_OUT="${OUT}/hcv_druggable_genes.tsv"
HCV_KINA_OUT="${OUT}/hcv_kinase_genes.tsv"
SARS_KINA_OUT="${OUT}/sars_kinase_genes.tsv"


cut -f1,10,12,16 ${FLS} | egrep "^DENV.+Kinome.+" | cut -f4 | egrep "\d+"  | sort -u > ${DENV_KINA_OUT}
cut -f1,10,12,16 ${FLS} | egrep "^DENV.+DruggableGenome.+" | cut -f4 | egrep "\d+" | sort -u > ${DENV_DUG_OUT}

cut -f1,10,12,16 ${FLS} | egrep "^HCV.+Kinome.+" | cut -f4 | egrep "\d+" | sort -u > ${HCV_KINA_OUT}
cut -f1,10,12,16 ${FLS} | egrep "^HCV.+DruggableGenome.+" | cut -f4 | egrep "\d+" | sort -u > ${HCV_DUG_OUT}

cut -f1,10,12,16 ${FLS} | grep "CHIKV" | cut -f4 | egrep "\d+" | sort -u > ${CHIKV_OUT}
cut -f1,10,12,16 ${FLS} | grep "SARS" | cut -f4 | egrep "\d+" | sort -u > ${SARS_KINA_OUT}