
pwm_db <- snakemake@params[["db"]]
snp_bed <- snakemake@input[["snp_bed"]]
outfile <- gzfile(snakemake@output[["out"]], open="w")
logfile <- file(snakemake@log[["logfile"]], open="w")
sink(logfile, type="message")
sink(logfile, type="output") 


library("motifbreakR")
library("MotifDb")
library("BSgenome.Hsapiens.UCSC.hg19")

snps <- snps.from.file(file = snp_bed, search.genome = BSgenome.Hsapiens.UCSC.hg19,format = "bed")

db <- NULL
if (pwm_db == "all") {
    pwm_human <- query(MotifDb, "hsapiens")
    db <- c(pwm_human, motifbreakR_motif)
} else {
    db <-  eval(as.symbol(pwm_db))
}

results <- motifbreakR(snpList = snps, threshold=0.9, pwmList = db, method = "default", bkg = c(A=0.25, C=0.25, G=0.25, T=0.25))
#results <- motifbreakR(snpList = snps, threshold=0.85, pwmList = db, method = "default", bkg = c(A=0.25, C=0.25, G=0.25, T=0.25))
df <- as.data.frame(results, row.names=NULL)
#df$snps <- names(results) # checked in the source code with  getMethod(as.data.frame,  "GenomicRanges"), the order is maintained
#> colnames(df)
# [1] "seqnames"     "start"        "end"          "width"        "strand"       "REF"          "ALT"          "snpPos"       "motifPos"    
#[10] "geneSymbol"   "dataSource"   "providerName" "providerId"   "seqMatch"     "pctRef"       "pctAlt"       "alleleRef"    "alleleAlt"   
#[19] "effect"       "sns"  

dfres <- data.frame(snpid=names(results), gene=df$geneSymbol, db=df$dataSource, name=df$providerName, id=df$providerId, scoreRef=df$pctRef, scoreAlt=df$pctAlt, effect=df$effect, seqMatch=df$seqMatch, motifPos= df$motifPos)
write.table(dfres, file=outfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
close(outfile)
