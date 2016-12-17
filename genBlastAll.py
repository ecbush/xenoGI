import sys

strainProtFNL = [line.rstrip() for line in open(sys.argv[1],'r')]

suffixToRemove = sys.argv[2]

for query in strainProtFNL:
    for db in strainProtFNL:
        if query != db:
            qstem=query.split(suffixToRemove)[0]
            qstem=qstem.split("fasta/")[1]
            
            dbstem=db.split(suffixToRemove)[0]
            dbstem=dbstem.split("fasta/")[1]
            
            file_out = "blast/" + qstem + '-' + dbstem + ".out"
            file_err = "blast/" + qstem + '-' + dbstem + ".err"


            print('qsub -cwd -V -e ' + file_err +\
                          ' -o ' + file_out +\
                          ' -b y /usr/bin/blastall -p blastp -M BLOSUM62 -G 11 -E 1 -e 0.01 -v 10000 -b 10000 -F \\"m S\\" -m 8 -z 300000000' +\
                          ' -d ' + db + ' -i ' + query)


# got these blast parms from Duret et colleagues. http://bioinformatics.oxfordjournals.org/content/28/8/1078.full
# -M BLOSUM62 -G 11 -E 1 -e 1e-04 -v 10000 -b 10000 -F "m S" -m 8 -z 300000000
# I've just changed the threshold, to make it more permissive.

# OLD parms I used                          " -b y /usr/bin/blastall -p blastp -e 0.01 -m 8" +\
