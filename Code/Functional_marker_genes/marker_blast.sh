#!/bin/bash

#Copy input file from gluster and unzip BLAST and database

tar xvzf ncbi-blast-2.7.1+-x64-linux.tar.gz
tar xvzf metabolic_genes.tar.gz

#Build db - this is super fast, so doing it each time
./ncbi-blast-2.7.1+/bin/makeblastdb -in metabolic_genes/Karthik_metabolic_genes.faa -dbtype prot  
#Run blastp
./ncbi-blast-2.7.1+/bin/blastp -query $1 -db metabolic_genes/Karthik_metabolic_genes.faa -outfmt 6 -out $1.out -evalue 0.0001

# remove anything below 30PIDthreshold=30
awk -v threshold="$threshold" $3 > threshold' $1.out > temp.txt
mv temp.txt $1.out




rm *txt
rm *tar.gz
rm -r ncbi-blast-2.7.1+/
rm -r metabolic_genes
rm *faa
