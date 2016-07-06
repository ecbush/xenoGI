def load(filename):
    """Load fasta or multifasta, return list of tuples (header,seq)."""
    f=open(filename,"r")
    header=None
    outL=[]
    tempSeqL=[]
    while True:
        Str=f.readline()
        if Str=="":
            outSeq="".join(tempSeqL)
            tempSeqL=[]
            outSeq="".join(outSeq.split()) # remove all whitespace
            outL.append((header,outSeq))
            break
        if Str[0]==">":
            # this is new header, put together previous header,seq, and move one
            if header!=None:
                outSeq="".join(tempSeqL)
                tempSeqL=[]
                outSeq="".join(outSeq.split()) # remove all whitespace
                outL.append((header,outSeq))
            header=Str[:-1]
        else:
            # this is a seq line,
            tempSeqL.append(Str)
    f.close()
    return(outL)    
