#!/bin/bash
PTH="/Users/simondi/PHD/data/data/sysvirdrug/"
MAP="${PTH}maps/"
FL="${MAP}entrez2hugo.tsv"

SCR=$(greadlink -f $0)
DIR=`dirname $SCR`

echo "Doing DENV ambion"
cut -f2,3 "${PTH}dengue_kinase_screen/kinases_ambion_siRNA_info.txt" | awk '{if (NR != 1) print ;}' | sort -u | tee -a ${FL} "${MAP}denv_kinome_entrez2hugo.tsv" > /dev/null

echo "Doing HCV ambion"
cut -f2,3 "${PTH}hcv_kinase_screen/kinases_ambion_siRNA_info.txt" | awk '{if (NR != 1) print ;}' | sort -u | tee -a ${FL}  "${MAP}hcv_kinome_entrez2hugo.tsv" > /dev/null

echo "Doing CHIKV dharmacon"
cut -f5-6 "${PTH}chikv_kinase_screen/kinases_dharmacon_siRNA_info.txt" | awk 'BEGIN{ OFS="\t" }{if (NR != 1) print $2, $1;}' | sort -u | tee -a ${FL}  "${MAP}chikv_kinome_entrez2hugo.tsv" > /dev/null

echo "Doing CVB ambion"
cut -f3 "${PTH}cvb_screen/cvb_replicon_screen.txt" | gawk ' BEGIN{OFS="\t"} { if(NR != 1)  print "NA", $0; }' | sort -u | tee -a ${FL} "${MAP}cvb_entrez2hugo.tsv" > /dev/null

echo "Doing DENV/HCV druggable ambion"
python "${DIR}/entrez2Hugo.py" "${PTH}dengue_druggable_genome_screen/ambion_library_layout" | sort -u | tee -a ${FL} "${MAP}denv_druggable_entrez2hugo.tsv" > /dev/null
python "${DIR}/entrez2Hugo.py" "${PTH}hcv_druggable_genome_screen/ambion_library_layout"    | sort -u | tee -a ${FL} "${MAP}hcv_druggable_entrez2hugo.tsv"  > /dev/null

echo "Doing SARS dharmacon"
cut -f5-6 "${PTH}sarscov_kinase_screen/kinases_dharmacon_siRNA_info.txt" | awk 'BEGIN{ OFS="\t" }{ if (NR != 1) print $2, $1; }' | sort -u | tee -a ${FL}  "${MAP}sarscov_kinome_entrez2hugo.tsv" > /dev/null

sort -u ${FL} > /tmp/dummy
mv /tmp/dummy ${FL}
wc -l ${FL}
