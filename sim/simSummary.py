import sys

simLogFN=sys.argv[1]

f=open(simLogFN,'r')

eventLenD={}

while True:
    s=f.readline()
    if s =='':
        break
    s=s.rstrip()
    L=s.split('\t')
    eventType = L[1]
    genesL = L[2].split(' ')
    
    if eventType in eventLenD:
        eventLenD[eventType].append(len(genesL))
    else:
        eventLenD[eventType] = [len(genesL)]
    
f.close()


print("Summaries")

for key in eventLenD:
    eventL=eventLenD[key]

    print(key)
    print("  number:",len(eventL))
    print("  min size:",min(eventL))
    print("  max size:",max(eventL))
