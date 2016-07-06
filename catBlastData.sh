
strain_file=$1

while read strain1
do
    while read strain2
    do
	if [ $strain1 != $strain2 ]
	then
            temp='blast/out/'$strain1'-'$strain2'.out'
            cat $temp
	fi
    done < $strain_file 
done < $strain_file
