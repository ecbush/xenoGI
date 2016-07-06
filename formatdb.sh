# copy protein fasta files

cp seq/broad/*simp.faa blast/databases/
cp seq/ncbi/*simp.faa blast/databases/

# use the same file list as before with the protein names

cd blast/databases/

ls *simp.faa > protsToFormat.txt

while read query
do
    formatdb -i $query
done < protsToFormat.txt


