strandedBamImport <- function (file, selection)
{
  if (!file.exists(paste(file, "bai", sep = ".")))
    stop("Unable to find index for BAM file '", file, "'. You can
         build an index using the following command:\n\t",
         "library(Rsamtools)\n\tindexBam(\"", file, "\")")
  sinfo <- scanBamHeader(file)[[1]]
  res <- if (!as.character(seqnames(selection)[1]) %in%
             names(sinfo$targets))
  {
    mcols(selection) <- DataFrame(score = 0)
    selection
  }else
  {
    param <- ScanBamParam(what = c("pos", "qwidth", "strand"),
                          which = selection, flag = scanBamFlag(isUnmappedQuery = FALSE, isDuplicate=FALSE))
    x <- scanBam(file, param = param)[[1]]
    if(length(x$pos)==0)
    {
      mcols(selection) <- DataFrame(plus=0, minus=0)
      selection
    }else
    {
      gr <- GRanges(strand=x[["strand"]],ranges=IRanges(x[["pos"]],width = x[["qwidth"]]),seqnames=seqnames(selection)[1])
      grs <- split(gr, strand(gr))
      cov <- lapply(grs[c("+", "-")],function(y){coverage(ranges(y),width=end(selection))})
      pos <- sort(unique(unlist(lapply(cov, function(y) c(start(y),end(y))))))
      
      GRanges(seqnames = seqnames(selection)[1],
              ranges=IRanges(start=head(pos, -1), end=tail(pos, -1)),
              plus=as.numeric(cov[["+"]][head(pos, -1)]),
              minus=-as.numeric(cov[["-"]][head(pos, -1)]),
              both=as.numeric(cov[["+"]][head(pos, -1)])+as.numeric(cov[["-"]][head(pos, -1)]))
    }
  }
  return(res)
}