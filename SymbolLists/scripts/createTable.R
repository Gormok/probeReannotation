library("optparse")
library(tools)
library(magrittr)
library(dplyr)

option_list = list(
   make_option(c("-f", "--file"), type="character", default=NULL,  
   help="Input file. Required; no default.", metavar="character"),
   
   make_option(c("-o", "--output"), type="character", default=NULL,
   help="Output file. Default: <input>_result.txt"),

   make_option(c("-m", "--multimaps"), type="character", default=NULL, 
   help="Multimaps output file. Default: NULL", metavar="character")
);  

getOpt = function()
{
   opt_parser = OptionParser(option_list=option_list);
   opt = parse_args(opt_parser);
   print(opt) 

   if (is.null(opt$file)) stop("No file to analyse")
   if (is.null(opt$output)) opt$output = tools::file_path_sans_ext(opt$file) %>% paste("_result.txt", sep="")

   return(opt)
}

long.file = function(file) file %>% {paste(file_path_sans_ext(.), "_long.", file_ext(.), sep = "")} 

executeScript = function(opt)
{ 
   breakdown = getBreakdown(opt$file)
   result = processBreakdown(breakdown, skip_multimaps = opt$multimaps %>% is.null)
   uniq.short = result$uniq.short
   uniq.short = uniq.short[order(uniq.short$probeID),]
   write.table(uniq.short, file = opt$output)
   
   uniq.long = result$uniq.long
   uniq.long = uniq.long[order(uniq.long$probeID),]
   write.table(uniq.long, file = long.file(opt$output))
   
   paste("Result table saved to", opt$output) %>% print

   if (is.null(opt$multimaps)) return
   
   multi.short = result$multi.short
   multi.short = multi.short[order(multi.short$probeID),]
   write.table(multi.short, file = opt$multimaps)

   multi.long = result$multi.long
   multi.long = multi.long[order(multi.long$probeID),]
   write.table(multi.long, file = long.file(opt$multimaps))

   paste("Multi table saved to", opt$multimaps) %>% print
}

processBreakdown = function(breakdown, skip_multimaps = F)
{
   short = function(df) df[c("probeID", "Ensembl")]
# unique autoincluded
   uniq = breakdown$unique %>% do.call(what=rbind) 
# unique.mix is also unique
   uniq.mix = lapply(breakdown$unique.mix, filter, nchar(Ensembl) == 15) %>% do.call(what = rbind)
   uniq.long = rbind(uniq, uniq.mix)
  
   uniq.short = uniq.long %>% short %>% unique 
   if (skip_multimaps) return(list(uniq = uniq.long, multi = NULL));
# unmatched probes 
   unmatched = c(breakdown$NF, breakdown$LQ, breakdown$NA_, breakdown$residue) %>% lapply(FUN = function(x) x[1, ]) %>% do.call(what=rbind)

# multimaps.s
   multimaps.s.mix.cleaned = lapply(breakdown$multimaps.s.mix, function(x) x[nchar(x[,"Ensembl"])>=15,])
   multimaps.s.all = append(breakdown$multimaps.s, multimaps.s.mix.cleaned)
# multimaps.g
   multimaps.g.mix.cleaned = lapply(breakdown$multimaps.g.mix, function(x) x[nchar(x[,"Ensembl"])>=15,])
   multimaps.g.all = append(breakdown$multimaps.g, multimaps.g.mix.cleaned)

   multimaps.s.all.short = lapply(multimaps.s.all, short) 
   multimaps.s.all.uniq = lapply(multimaps.s.all.short, unique)
   multimaps.g.all.short = lapply(multimaps.g.all, short)
   multimaps.g.all.uniq = lapply(multimaps.g.all.short, function(x)
           {
              ens = strsplit(x$Ensembl, split="\\|") %>% unlist %>% unique
              return(data.frame(probeID = rep(x$probeID[1], length(ens)), Ensembl = ens))
           })
   multimaps.df = c(multimaps.s.all.uniq, multimaps.g.all.uniq) %>% do.call(what = rbind)
   uncovered = setdiff(multimaps.df$Ensembl, uniq.short$Ensembl)
   
   extra.short = multimaps.df[multimaps.df$Ensembl %in% uncovered,]

   extra.long.s = multimaps.s.all %>% .[names(.) %in% extra.short$probeID] %>%
                do.call(what = rbind) %>% merge(extra.short) 
   extra.long.g = multimaps.g.all %>% .[names(.) %in% extra.short$probeID] %>%
                  do.call(what = rbind) %>% 
                  (function(df)
                  {
                      ens = strsplit(df$Ensembl, split = "\\|")
                      l = lapply(ens, length)

                      lines = mapply(rep, seq_along(df$Ensembl), l, SIMPLIFY = F) %>% do.call(what = c)
                      df2 = df[lines,]
                      df2$Ensembl = ens %>% unlist
                      return(df2)
                  }) %>% merge(extra.short)

   extra.long = rbind(extra.long.s, extra.long.g)

   extra.geneCount = extra.short$Ensembl %>% unique %>% length
   unused = setdiff(multimaps.df$probeID, extra.short$probeID)

   cat(uniq.short %>% nrow, "probes were matched uniquely.\n")
   cat(extra.short %>% nrow, "extra symbolList entries were added for uncovered genes.\n")
   cat(extra.geneCount, "genes were covered by the entries.\n")
   cat(unused %>% length, "probes were left unused.\n")
   cat("Other entries:")
   for (n in breakdown %>% names %>% grep(pattern="XF:Z", value=T) %>% c("residue.mix"))
   {
       cat(n,":", length(breakdown[[n]]), "\n")
   }

   multi.long = rbind(uniq.long, extra.long)
   multi.short = rbind(uniq.short, extra.short)
   #foo = function(x)
   #{
   #   max = summary(x[,"Ensembl"] %>% as.factor) %>% sort(decreasing=T)
   #   if (max[1]>=3*max[2]) return(names(max)[1])
   #   else return(unique(x[,"Ensembl"]) %>% paste(collapse="|"))
   #}

   result = list(uniq.long = uniq.long, uniq.short = uniq.short, multi.long = multi.long, multi.short = multi.short)

   return(result)
}



getBreakdown = function(file)
{
   orig.df = read.table(file, header=F, comment.char="", quote="", stringsAsFactors = F)
   colnames(orig.df) = c("probeID", "start", "CIGAR", "Ensembl")

   cleaned.df = orig.df
   cleaned.df$Ensembl = cleaned.df$Ensembl %>% gsub(pattern = ".*?(ENSG[[:digit:]]*)", replacement="\\1\\|") %>% gsub(pattern = "(.*)\\|.*", replacement = "\\1")
   residue.tags = c("XF:Z:__no_feature" = "NFE",
           "XF:Z:__too_low_aQual" = "TKQ",
           "XF:Z:__not_aligned" = "NAL",
           "XF:Z:__alignment_not_unique" = "ANU")
   for (n in names(residue.tags))
   {
       cleaned.df$Ensembl[cleaned.df$Ensembl == n] = residue.tags[n]
   }
   split.df = split(cleaned.df, f=cleaned.df$probeID %>% as.factor)

   residue.pure = lapply(residue.tags, function(tag) sapply(split.df, function(probe) all(probe[,"Ensembl"]==tag)))
   residue.mix = sapply(split.df, function(x) all(x[,"Ensembl"] %in% residue.tags) && !all(x[1,"Ensembl"]==x[,"Ensembl"]))

   unique = sapply(split.df, function(x) all(x[1,"Ensembl"]==x[,"Ensembl"]) && all(nchar(x[,"Ensembl"])==15)) 
   multimaps.s = sapply(split.df, function(x) any(x[1,"Ensembl"]!=x[,"Ensembl"]) && all(nchar(x[,"Ensembl"])==15)) 
   multimaps.g = sapply(split.df, function(x) all(nchar(x[,"Ensembl"])>=15) && any(nchar(x[,"Ensembl"])>15)) 

   unique.mix = sapply(split.df, function(x) {y=x[nchar(x[,"Ensembl"])==15,]; all(y[1,"Ensembl"]==y[, "Ensembl"]) && any(nchar(x[,"Ensembl"])<15) && any(nchar(x[,"Ensembl"])==15) && all(nchar(x[,"Ensembl"])<=15)})
   multimaps.s.mix = sapply(split.df, function(x) {y=x[nchar(x[,"Ensembl"])==15,]; any(y[1,"Ensembl"]!=y[, "Ensembl"]) && all(nchar(x[,"Ensembl"])<=15) && any(nchar(x[,"Ensembl"])!=15)})
   multimaps.g.mix = sapply(split.df, function(x) {y=x[nchar(x[,"Ensembl"])==15,]; any(nchar(x[,"Ensembl"])<15) && any(nchar(x[,"Ensembl"])>15)})


   flatten = function(l)
   {
       if (class(l) != "list") return(list(l))
       else res = lapply(l, flatten) %>% do.call(what=c)
       return(res)
   }
   union.all = function(...)
   {
      args = list(...) %>% flatten 
      if (length(args) == 1L) return(which(args[[1]]))
      else return(union(which(args[[1]]), do.call(union.all, args[-1] )));
   }
   indices = union.all(residue.pure, residue.mix, unique, multimaps.s, multimaps.g, unique.mix, multimaps.s.mix, multimaps.g.mix)
   diff = setdiff(seq_along(split.df), indices)

   breakdown = residue.pure %>% c(., list(residue.mix = residue.mix, unique = unique, multimaps.s = multimaps.s,
       multimaps.g = multimaps.g, unique.mix = unique.mix, multimaps.s.mix = multimaps.s.mix, multimaps.g.mix = multimaps.g.mix))
   breakdown %<>% lapply(FUN=which)
   breakdown.df = lapply(breakdown, function(x) split.df[x])
   unique.df.mult = breakdown.df$unique[sapply(breakdown.df$unique, nrow)>1]
      
   all.indices = do.call(what=c, breakdown)
   uniqueness = anyDuplicated(all.indices)
   if (uniqueness==0 && length(diff)==0)
   {
      print("SUCCESS: Breakdown is a disjoint union and covers all indices.")
   } else if (uniqueness!=0 && length(diff)==0)
   {
      warning("Warning: Breakdown reuses certain indices. Check for new cases.")
      print("Warning: Breakdown reuses certain indices. Check for new cases.")
      print("Breakdown covers all indices.")
   } else if (uniqueness==0 && length(diff)!=0)
   {
      warning("Warning: Not all indices have been covered. Check for new cases.")
      print("Warning: Not all indices have been covered. Check for new cases.")
      print("Breakdown is a disjoint union.")
   } else
   {
      stop("Both uniqueness and disjoint union checks failed. Check implementation thoroughly.")
   }
   return(breakdown.df)
}


####
#### code to run

opt = getOpt()
executeScript(opt)
