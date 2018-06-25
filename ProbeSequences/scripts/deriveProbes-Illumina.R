library(stringi) # stri_dup
library(magrittr)
source("common/probeSequences.R") # wl.file

##
dataFile = snakemake@input[["txt"]]

probesAll.file = snakemake@log[["probesAll"]]
probesGenuine.file = snakemake@log[["probesGenuine"]]
probesControl.file = snakemake@log[["probesControl"]] 
rest.file = snakemake@log[["rest"]]

fastq.file = snakemake@output[["fastq"]]

### GET FILE STRUCTURE ###

lines = readLines(dataFile)

# starts of file sections
sectionStarts = sapply(lines, function(x) (substr(x,1,1)=="[" && substr(x, nchar(x), nchar(x)) =="]")) %>% which
names(sectionStarts) = names(sectionStarts) %>% sapply(FUN = function(x) substr(x, 2, nchar(x)-1))

# add proxy 'Start' and 'End'
sectionStarts.ext = sectionStarts%>% c(Start = 1, . , End = length(lines)+1)

lines.factor = cut(seq_along(lines) %>% {.[sectionStarts]=-1;.}, breaks = sectionStarts.ext)
lines.split = split(lines, lines.factor) %>% setNames(nm = sectionStarts.ext %>% names %>% head(., -1))
sections.len = sapply(lines.split, length)

#borders = sectionStarts.ext %>% {mapply(c, head(.,-1)+1, tail(.,-1)-1)} %>% {lapply(colnames(.) %>% setNames(nm=.), function(i) .[,i])}


### DERIVE PROBE MATRIX ###
df.headers = c(Probes = "Probes", Controls = "Controls")

matrices = lapply(df.headers, 
      function(x) read.table(dataFile, stringsAsFactors = F, header = T, skip = sectionStarts[x], nrows = sections.len[x]-1, quote="", comment.char="", sep="\t")
      )

# Use only selected columns
desiredCols = c("Probe_Id", "Array_Address_Id", "Probe_Sequence")
matrices.trimmed = lapply(matrices, function(x) x[,desiredCols])
full.matrix = matrices.trimmed %>% do.call(what=rbind)

# Write the tables read
write.table(matrices[["Probes"]], file = probesGenuine.file)
write.table(matrices[["Controls"]], file = probesControl.file)
write.table(full.matrix, file = probesAll.file)


# Save non-probe data (metadata)
lines.factor = cut(seq_along(lines), breaks = sectionStarts.ext, right=F) 
lines.split = split(lines, lines.factor) %>% setNames(nm = sectionStarts.ext %>% names %>% head(., -1))
lines.remainder = names(sectionStarts.ext) %>% head(-1) %>% {lines.split[.[!(. %in% df.headers)]]} %>% unlist
wl.file(lines.remainder, rest.file)

#### SANITY CHECK
heading =  read.table(dataFile, stringsAsFactors = F, header = F, skip = sectionStarts["Heading"], nrows = sections.len["Heading"], quote="", comment.char="", sep="\t")
rownames(heading) = heading[,1]
b1 = heading["Number of Probes",2] == nrow(matrices[["Probes"]])
b2 = heading["Number of Controls", 2] == nrow(matrices[["Controls"]])

if (!b1) warning("Loading may not have been successfull. Number of Probes in header mismatches data frame.")
if (!b2) warning("Loading may not have been successfull. Number of Controls in header mismatches data frame.")
if (b1 && b2) print("Loading successfull.")

### FASTQ FILE
# Create fastq representation
fastq.lines1 = paste("@", full.matrix[,"Array_Address_Id"], sep="")
fastq.lines3 = paste("+", full.matrix[,"Array_Address_Id"], sep="")
fastq.lines2 = full.matrix[, "Probe_Sequence"]
fastq.lines4 = sapply(fastq.lines2, nchar) %>% stri_dup(str="~") 

fastq.lines.mat = mapply(c, fastq.lines1, fastq.lines2, fastq.lines3, fastq.lines4)
fastq.lines = fastq.lines.mat %>% as.character

# Save fastq file
wl.file(fastq.lines, fastq.file)
