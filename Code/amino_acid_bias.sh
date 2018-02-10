#!/bin/bash
genome=$1
grep -v ">" $genome > aminos.txt
num1=$(grep -o "R" aminos.txt | wc -c | cut  -f1 -d' ')
num2=$(grep -o "H" aminos.txt | wc -c | cut  -f1 -d' ')
num3=$(grep -o -E "N|Q|K|W" aminos.txt | wc -c | cut  -f1 -d' ')
total=$(wc -c aminos.txt | cut  -f1 -d' ')
num4=$(($total - ($num1 + $num2 + $num3)))

sum=$(($num1 * 4 + $num2 * 3 + $num3 * 2 + $num4))
index=$(echo "scale=3; $sum/$total" | bc -l)

echo $1 $index
