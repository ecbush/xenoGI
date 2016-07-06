import os

strainProtFNL = [line.rstrip() for line in open("blast/databases/simpProtFiles.txt",'r')]

outfile = open("code/blastAll.sh", 'w')

for query in strainProtFNL:
    for db in strainProtFNL:
        if query != db:
            qstem=query.split('.simp.faa')[0]
            qstem=qstem.split("blast/databases/")[1]
            
            dbstem=db.split('.simp.faa')[0]
            dbstem=dbstem.split("blast/databases/")[1]
            
            file_out = "blast/out/" + qstem + '-' + dbstem + ".out"
            file_err = "blast/out/" + qstem + '-' + dbstem + ".err"


            outfile.write("qsub -cwd -V -e " + file_err +\
                          " -o " + file_out +\
                          " -b y /usr/bin/blastall -p blastp -e 0.01 -m 8" +\
                          " -d " + db + " -i " + query + '\n')


outfile.close()
