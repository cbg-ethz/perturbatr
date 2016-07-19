#!/bin/bash

FLS="/Users/simondi/PHD/data/data/sysvirdrug/integrated_data_files"

CHIKV="${FLS}/chikv_kinase_screen.tsv"
DENV_DRUG="${FLS}/denv_druggable_genome_screen.tsv"
DENV_KINA="${FLS}/denv_kinase_screen.tsv"
HCV_DRUG="${FLS}/hcv_druggable_genome_screen.tsv"
HCV_KINA="${FLS}/hcv_kinase_screen.tsv"
SARS_KINA="${FLS}/sars_kinase_screen.tsv"


OUT="/Users/simondi/PHD/data/data/sysvirdrug/siRNAs"
CHIKV_OUT="${OUT}/chikv_kinase_sirnaID2Seq.tsv"
DENV_DUG_OUT="${OUT}/denv_druggable_sirnaID2Seq.tsv"
DENV_KINA_OUT="${OUT}/denv_kinase_sirnaID2Seq.tsv"
HCV_DUG_OUT="${OUT}/hcv_druggable_sirnaID2Seq.tsv"
HCV_KINA_OUT="${OUT}/hcv_kinase_sirnaID2Seq.tsv"
SARS_KINA_OUT="${OUT}/sars_kinase_sirnaID2Seq.tsv"


cut -f 22,21 ${SARS_KINA}  | grep "$[A|C|G|U]"  | awk 'BEGIN { OFS="\t" } { split($1, seqs, ","); split($2, ids, ","); for (i=1; i <= length(ids);i++) { print ids[i], seqs[i]} }' | sort -u > ${SARS_KINA_OUT};

cut -f8,11 ${HCV_KINA} | egrep "^\\d+"  | sort -u > ${HCV_KINA_OUT}

cut -f 14,12 ${HCV_DRUG} | 
egrep "^\\d+"  | 
awk 'BEGIN { OFS="\t" } { print $1, $2 }' | sort -u > ${HCV_DUG_OUT}

cut -f8,11 ${DENV_KINA} | egrep "^\\d+"  | sort -u > ${DENV_KINA_OUT}

cut -f 14,12 ${DENV_DRUG} | 
egrep "^\\d+"  | 
awk 'BEGIN { OFS="\t" } { print $1, $2 }' | sort -u > ${DENV_DUG_OUT}

cut -f 22,21 $CHIKV  | grep "$[A|C|G|U]"  | awk 'BEGIN { OFS="\t" } { split($1, seqs, ","); split($2, ids, ","); for (i=1; i <= length(ids);i++) { print ids[i], seqs[i]} }' | sort -u > ${CHIKV_OUT};

