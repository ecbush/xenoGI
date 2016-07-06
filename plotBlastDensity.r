
coreSc=scan("coreSimScores.txt")
allSc=scan("allSimScores.txt")

coreDens=density(coreSc)
allDens=density(allSc)

# scale the core density so its comparable to full
sc=length(coreSc)/length(allSc)
coreDensScaled=coreDens
coreDensScaled$y=sc*coreDens$y

pdf("coreSim.pdf",width = 7, height = 5)
plot(coreDensScaled,col='red')
points(allDens,type="l")
dev.off()
