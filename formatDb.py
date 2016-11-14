import sys,os


for dbFileName in sys.argv[1:]:
    os.system('formatdb -i '+dbFileName)
