#!/bin/bash

awk '{print $1}' MAGstravaganza/metabolic_gene_info.txt > markers.txt

for file in *blast; do
	touch results.txt;
	while read gene; do
		grep $gene $file > subtable.out;
		while read line;
			do echo $line > row.txt;
			PID=$(awk '{print $3}' row.txt);
			compare=$(echo "${PID} > 30" | bc);
			echo $compare;
			done < subtable.out > hits.txt;
		counts=$(awk '{ sum += $1 } END { print sum }' hits.txt);
		echo $gene $counts >> results.txt;
	done < markers.txt;
	awk '{if (!$2) {print $1, "0"} else {print $1, $2}}' results.txt > $file.results.txt;
	rm results.txt;
done

awk '{print $1}' PUBU.len150.blast.results.txt > rownames.txt
results=$(echo *.results.txt)
awk '{print $2)' $results > allcounts.txt
paste rownames.txt allcounts.txt > marker_gene_counts.txt
