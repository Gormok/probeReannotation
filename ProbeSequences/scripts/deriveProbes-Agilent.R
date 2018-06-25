library(stringi) # stri_dup
library(magrittr)
source("common/probeSequences.R") # wl.file

# Retrieve files
input.file = snakemake@input[["txt"]]

probesFile = snakemake@log[["probes"]]
fastqFull.file = snakemake@log[["fastqFull"]]

fastq.file = snakemake@output[["fastq"]]



# Count initial comment lines
lines = readLines(input.file)
skip = 0
for (l in lines) if (substr(l,1,1)=="#") skip=skip+1 else break

# Read the file from start
a=read.table(input.file, header=T, stringsAsFactors=F, quote="", comment.char="",sep="\t", skip=skip)

# Keep only certain columns
desiredCols = c("ID", "COL", "ROW", "SPOT_ID", "CONTROL_TYPE", "GENE_SYMBOL", "SEQUENCE")
cleaned = a[,desiredCols]

# Save the kept table
write.table(cleaned, file=probesFile)

# Prepare fastq file
fastq.lines1 = paste("@", cleaned[,"ID"], sep="")
fastq.lines3 = paste("+", cleaned[,"ID"], sep="")
fastq.lines2 = cleaned[, "SEQUENCE"]
fastq.lines4 = sapply(fastq.lines2, nchar) %>% stri_dup(str="~")

fastq.lines.mat = mapply(c, fastq.lines1, fastq.lines2, fastq.lines3, fastq.lines4)
fastq.lines = fastq.lines.mat %>% as.character

# Save full fastq for Agilent
wl.file(fastq.lines, fastqFull.file)

# Filter out lines without sequence (certain tools do not expect this)
fastq.lines.mat.clean = fastq.lines.mat[,sapply(fastq.lines2, nchar)>0]
fastq.lines.clean = fastq.lines.mat.clean %>% as.character

# Save cleaned fastq for Agilent
wl.file(fastq.lines.clean, fastq.file)
