# Reads quality
library(ShortRead)
fqSample <- FastqSampler("SRR3020039.fastq.gz", n = 10^6)  
fastq <- yield(fqSample)
fastq             
readQuality <- quality(fastq)
readQualities<-alphabetScore(readQuality)
library(ggplot2)
toPlot <- data.frame(ReadQ = readQualities)
ggplot(toPlot, aes(x = ReadQ)) + geom_histogram() + theme_minimal()
readSequences <- sread(fastq)
readSequences_AlpFreq <- alphabetFrequency(readSequences)
summed__AlpFreq <- colSums(readSequences_AlpFreq)
summed__AlpFreq[c("A", "C", "G", "T", "N")]
readSequences_AlpbyCycle <- alphabetByCycle(readSequences)
AFreq <- readSequences_AlpbyCycle["A", ]
CFreq <- readSequences_AlpbyCycle["C", ]
GFreq <- readSequences_AlpbyCycle["G", ]
TFreq <- readSequences_AlpbyCycle["T", ]
toPlot <- data.frame(Count = c(AFreq, CFreq, GFreq, TFreq), Cycle = rep(1:51, 4), Base = rep(c("A", "C", "G", "T"), each = 51))
ggplot(toPlot, aes(y = Count, x = Cycle, colour = Base)) + geom_line() + theme_bw()
qualAsMatrix <- as(readQuality, "matrix")
boxplot(qualAsMatrix[1:100000, ])

# data filter
fqStreamer <- FastqStreamer("SRR3020076.fastq.gz")
TotalReads <- 0
TotalReadsFilt <- 0
while (length(fq <- yield(fqStreamer)) > 0) {
    TotalReads <- TotalReads + length(fq)
    filt1 <- fq[alphabetScore(fq) > 300]
    filt2 <- filt1[alphabetFrequency(sread(filt1))[, "N"] < 5]
    TotalReadsFilt <- TotalReadsFilt + length(filt2)
    writeFastq(filt2, "filtered_SRR3020076.fastq.gz", mode = "a")
}
