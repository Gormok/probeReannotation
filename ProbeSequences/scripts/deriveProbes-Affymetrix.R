library(stringi) # stri_dup
library(magrittr)
source("common/probeSequences.R") # wl.file

# Get output files
allProbes.file = snakemake@log[["allProbes"]]
aSet.file = snakemake@input[["aSet"]]
bSet.file = snakemake@input[["bSet"]]

fastq.file = snakemake@output[["fastq"]]
print(fastq.file)
# data -> fastq file by lines
affToFastq = function(data)
{
   probe.names = data$Probe.Set.Name
   fastq.lines1 = paste("@", probe.names, sep="")
   fastq.lines3 = paste("+", probe.names, sep="")
   fastq.lines2 = data[, "sequence"]
   fastq.lines4 = sapply(fastq.lines2, nchar) %>% stri_dup(str="~")

   fastq.lines.mat = mapply(c, fastq.lines1, fastq.lines2, fastq.lines3, fastq.lines4)
   fastq.lines = fastq.lines.mat %>% as.character
   return(fastq.lines)
}

# retrieve Affymetrix data
hgu133aprobe = read.table(aSet.file, header = T, stringsAsFactors = F)
hgu133bprobe = read.table(bSet.file, header = T, stringsAsFactors = F)

# Write Affymetrix data to file
write.table(rbind(hgu133aprobe, hgu133bprobe), file = allProbes.file)

# Convert Affymetrix data to fastq file
a.fastq = affToFastq(hgu133aprobe)
b.fastq = affToFastq(hgu133bprobe)
fastq=c(a.fastq, b.fastq)
wl.file(fastq, fastq.file)
