library("optparse")
library(venn)
library(ggplot2)
library(magrittr)

option_list = list(
   make_option(c("-c", "--compare"), type="character", default=NULL,  
   help="Probe symbolList files. Required. Use ':' as delimiter.", metavar="character"),
   
   make_option(c("-n", "--names"), type="character", default=NULL,
   help="Probe symbolList file labels.", metavar="character"),

   make_option(c("-o", "--odir"), type="character", default=NULL,
   help="Output directory. Required."),

   make_option(c("-r", "--reference"), type="character", default=NULL,
   help="Reference probe symbolList. Probe -> Ensembl, or Probe -> Entrez. Default: none."),

   make_option(c("-m", "rm", "--map"), type="character", default=NULL, 
   help="Table relating Entrez with Ensembl. Use if reference is Probe -> Entrez.")

);  

getOpt = function()
{
   opt_parser = OptionParser(option_list=option_list);
   opt = parse_args(opt_parser);
   print(opt) 

   if (is.null(opt$compare)) stop("No files to analyse.")
   if (is.null(opt$odir)) stop("Output directory must be specified")
   disassemble = function(string, delimiter = ":") strsplit(string, split=delimiter) %>% unlist
   opt$compare = disassemble(opt$compare)
   if (!is.null(opt$names))
   {
      nm=disassemble(opt$names)
      if (length(nm)!=length(opt$compare)) stop("'compare' and 'names' do not match in length")
      opt$compare = setNames(opt$compare, nm = nm)
   } else opt$compare = setNames(opt$compare, nm = letters[seq_along(opt$compare)])

   if (length(opt$compare)==1) stop("Single comparison file.")
   
   if (length(opt$compare)>3) warning("More than three comparison files can result in illegible plots.")
   print(opt)
   return(opt)
}

getOpt2 = function()
{
   commandArgs(trailingOnly=F) %>% print
   args = commandArgs(trailingOnly=TRUE)
   args.raw = sapply(args, strsplit, split="=")
   args.names = sapply(args.raw, function(x) x[[1]][1])
   args.vals = lapply(args.raw, function(x) 
            if (length(x)>1) paste(x %>% tail(-1), collapse="=")
            else x
         )
   args.result = setNames(args.vals, nm=args.names)
   return(args.result)
}

# assumes same tool being used
executeScript = function(args)
{
   col.tools = c("red", "blue", "green")
   # symbolList tables
   lists.df = lapply(args$compare, read.table, stringsAsFactors=F, header=T)
   # retrieve probeIDs present in individual symbolLists
   probes = lapply(lists.df, subset, select = "probeID", drop=T)

   # probes shared across all datasets
   int = Reduce(intersect, probes)
   # Ensembl IDs assigned to each probe by each dataset, put in a data frame
   match = lapply(lists.df, function(x) x[int, "Ensembl"]) %>% do.call(what=data.frame) 
   rownames(match) = int

   # list of form probe#match for each comparison file
   match.str = lapply(setNames(nm=colnames(match)),
                   function(x) paste(rownames(match), match[,x], sep="#"))
   match.str.all = lapply(lists.df, function(x) paste(x[,"probeID"], x[,"Ensembl"], sep="#"))

   # comparison of assigned probes 
   pdf(paste(args$odir, "assigned.pdf", sep="/"))
   venn(probes, snames = names(probes), zcolor = col.tools)
   dev.off()

   # mismatch comparison
   pdf(paste(args$odir, "matches.pdf", sep="/"))
   venn(match.str.all, snames = names(probes), zcolor = col.tools)
   venn(match.str, snames = names(probes), zcolor = col.tools)
   dev.off()

   if (!is.null(args$reference))
   {
      symbolList = read.table(args$reference, stringsAsFactors = F, header = T) 
      if (!is.null(args$map))
      {
         # remap symbolList
         map = read.table(args$map, stringsAsFactors = F, header = T)
         symbolList = merge(map, symbolList)
      } 
   }

   match2 = match
   match2$probeID = rownames(match2)
   match.ref = merge(match2, symbolList)
   for (name in colnames(match))
   {
      boolCol = match.ref[[name]]==match.ref[["Ensembl"]]
      match.ref[[name]] = boolCol
   }
   match.ref = t(match.ref)


   match.str.ref = paste(symbolList[,"probeID"], symbolList[,"Ensembl"], sep="#")
   match.vals = lapply(match.str.all, intersect, match.str.ref)
   match.count = sapply(match.vals, length)
   match.total = sapply(match.str.all, length)
   match.probes = sapply(lists.df, . %$% probeID %>% unique %>% length)
   match.genes = sapply(lists.df, . %$% Ensembl %>% unique %>% length)
   match.df = data.frame(Tool = names(match.str.all), Matches = match.count, Total = match.total, Probes = match.probes, Genes = match.genes)

   match.probeID = match.vals %>% lapply(FUN = . %>% sub(pattern="#.*", replacement=""))
   mismatches = sapply(seq_along(lists.df), function(i) lists.df[[i]] %$% probeID %>% unique %>% {! . %in% match.probeID[[i]]} %>% sum)
   match.df$Mismatches = mismatches

   print(match.df)
   lims = function(data, cols, dif = 5)
   {
      ydiff = max(data[cols]) - min(data[cols])
      ymin = min(data[cols]) - ydiff %/% dif
      ymax = max(data[cols]) + ydiff %/% dif
      return (c(ymin, ymax))
   }
   g.m = ggplot(data=match.df, aes(x = Tool, y = Matches, fill = Tool)) +
         geom_bar(stat="identity") +
         coord_cartesian(ylim=lims(match.df, "Matches"))
   ggsave(filename = paste(args$odir, "refmatch.pdf", sep="/"), plot = g.m)

   g.t = ggplot(data=match.df, aes(x = Tool, y = Total, fill = Tool)) +
         geom_bar(stat="identity") +
         coord_cartesian(ylim=lims(match.df, "Total"))
   ggsave(filename = paste(args$odir, "total.pdf", sep="/"), plot = g.t)


   g.p = ggplot(data=match.df, aes(x = Tool, y = Probes, fill = Tool)) +
         geom_bar(stat = "identity") +
         coord_cartesian(ylim = lims(match.df, "Probes"))
   ggsave(filename = paste(args$odir, "probes.pdf", sep="/"), plot = g.p)

   g.g = ggplot(data=match.df, aes(x = Tool, y = Genes, fill = Tool)) +
         geom_bar(stat = "identity") +
         coord_cartesian(ylim = lims(match.df, "Genes"))
   ggsave(filename = paste(args$odir, "genes.pdf", sep="/"), plot = g.g)

   g.mm = ggplot(data=match.df, aes(x = Tool, y = Mismatches, fill = Tool)) +
         geom_bar(stat = "identity") +
         coord_cartesian(ylim = lims(match.df, "Mismatches"))
   ggsave(filename = paste(args$odir, "mismatches.pdf", sep="/"), plot = g.mm)


   # match.df was wide format
   # we need long format for this (we could have used melt from reshape2)
print(match.df)
   g.mt = ggplot(data=match.df, aes(x = Total, y = Matches, colour = Tool)) +
          geom_point() +
          coord_cartesian(xlim = lims(match.df, "Total", 8), ylim=lims(match.df, "Matches", 8))
   ggsave(filename = paste(args$odir, "total-refmatch.pdf", sep="/"), plot = g.mt)


}

####
#### code to run

args = getOpt()
executeScript(args)
