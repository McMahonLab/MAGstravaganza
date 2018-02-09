# Unzip IMG files and move the fasta file to a different folder

for NUM in *; do
	tar -xzf $NUM
done

rm *.tar.gz

for DIR in *; do
	NEWDIR="$(basename $DIR)"
	mv $NEWDIR/*.fna ../$1
done

cd ../$1
rm *gen*.fna
