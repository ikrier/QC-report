setwd("/data/generate_QC_reports/temp_temp_dir/")
costable=read.csv("cosmtable.csv",header=F)

fastqdir="/data/fastq_gwendal/"
libname="49"
cosmname=costable$V2[as.numeric(libname)]
if(cosmname=="NONE"){cosmname="Aucun"}

regionsname="../3rdDesignRegions.bed"

reportname=paste("Report",libname,"_draft",sep="")

Sweave(paste(reportname,".Rnw",sep=""))
tools::texi2pdf(paste(reportname,".tex",sep=""))
