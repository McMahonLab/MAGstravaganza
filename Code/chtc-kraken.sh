#!/bin/bash
#CHTC kraken run on lake metagenomes

tar -zxvf chtc_testing.tgz
cd chtc_testing
tar -zxvf  minikraken.tgz
tar -zxvf kraken-0.10.5-beta.tgz
cd kraken-0.10.5-beta
./install_kraken.sh kraken_scripts
cd kraken_scripts
./kraken --threads 3 --preload --db ../../minikraken_20141208/ ../../genome.fasta > chtc_testkraken.output
./kraken-report --db ../../minikraken_20141208/ chtc_testkraken.output > chtc_testkraken.report

