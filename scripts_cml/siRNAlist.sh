#!/bin/bash

FLS="/Users/simondi/PHD/data/data/sysvirdrug/siRNAs"

for p in `find $FLS -type f -name "*ID2Seq*"`
do
	base=$(basename $p)
	pref=($(IFS="_"; echo $base))
	flout="${FLS}/${pref[0]}_${pref[1]}_sirnas.tsv"
	cut -f1 $p > ${flout}
done	