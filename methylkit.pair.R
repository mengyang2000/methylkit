library("methylKit")

#data_path     = "/home/meng_yang/project/201704/CR-methylkit/meth_count"
#compareGroup  =  "S10-DPM-CR003-002T-FR1-IRM_S10_R1_001,S12-DPM-CR003-003T2-FR1-IRM_S12_R1_001,S02-DPM-CR003-006T-FR1-IRM_S02_R1_001--S02-DPM-CR003-002N-FR1-IRM_S02_R1_001,S11-DPM-CR003-003T1-FR1-IRM_S11_R1_001,S03-DPM-CR003-003N-FR1-IRM_S03_R1_001,S05-DPM-CR003-006N-FR1-IRM_S05_R1_001"
#compareGroupname       =  "tumor--normal"
#### usage: Rscript script.R  count_path comparegroup outname out_path
args <- commandArgs(TRUE)
data_path     = args[1] #
compareGroup  = args[2] # 
compareGroupname       = args[3] # 
patient_ID = args[4]
out=args[5]

groupname = unlist(strsplit(compareGroupname,"--"))
sampleName  = unlist(strsplit(compareGroup,"--|,"))
compareGroupNames = unlist(strsplit(compareGroup,"--"))
sampleName1 = unlist(strsplit(compareGroupNames[1],","))
sampleName2 = unlist(strsplit(compareGroupNames[2],","))
treatment =  c(rep(1,length(sampleName1)),rep(0,length(sampleName2)))
print (treatment)

patientID = unlist(strsplit(patient_ID,"--|,"))
covariates = data.frame(patientID=patientID)
print(covariates)
outdir = paste(out,compareGroupname,sep="/")
dir.create(outdir,recursive = TRUE)
outname = paste(outdir,compareGroupname,sep="/")
file.list <- lapply(sampleName, function(f){list.files(path=data_path, pattern = f, full.names = TRUE)})
myobj=methRead(file.list,sample.id=as.list(sampleName),assembly="hg19",treatment=treatment,context="CpG",mincov=10)

#=================QC======================#
getMethylationStats(myobj[[2]],plot=FALSE,both.strands=FALSE)
file.pdf=paste(outname,".plot.pdf",sep="")
pdf(file.pdf,width=10,height=80)
par(mfrow=c(9,1))
getMethylationStats(myobj[[2]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj[[2]],plot=TRUE,both.strands=FALSE)
#filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)
meth=unite(myobj, destrand=FALSE,mc.cores=10)
#n=as.integer(round(0.8*min(length(sampleName1),length(sampleName2))))
#meth = unite(myobj, destrand=FALSE,min.per.group=n,mc.cores=10)
write.table(data.frame(meth), file=paste(outname,".meth.data.xls",sep=""),sep="\t",row.names=F,quote=F)
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
#hc = clusterSamples(meth, dist="correlation", method="ward", plot=FALSE)
PCASamples(meth, screeplot=TRUE)
PCASamples(meth)
dev.off()
#============calculate diff===============#
myDiff=calculateDiffMeth(meth,covariates=covariates,overdispersion="MN",mc.cores=10)
write.table(myDiff, file=paste(outname,".alldiff.meth.xls",sep=""),sep="\t",row.names=F,quote=F)
for (i in c(0.25,0.5,0.6,0.7,0.75,0.8,0.85,0.9,0.95)){
	for (j in c(0.5,0.3,0.2,0.1,0.05,0.01)){
	cutoff_i=quantile(data.frame(myDiff)[,"meth.diff"],i)
	select_diff=getMethylDiff(myDiff,difference=cutoff_i,qvalue=j)
	tmp=paste(paste(paste(i,cutoff_i,sep=":"),j,sep="-"),".dml.xls",sep="")
	write.table(data.frame(select_diff),file=paste(outname,tmp,sep="."),sep="\t",row.names=F,quote=F)
	}
}

#=====================plot=======================#
myDiff=data.frame(myDiff)
file.png=paste(outname,".diff.p.png",sep="")
png(file.png,width = 480, height = 480)
par(mfrow=c(2,2))
mygrey <- rgb(89/255,87/255,87/255)
par(family="sans",col.main=mygrey,col=mygrey,col.lab=mygrey,col.axis=mygrey,cex.main=2,cex.lab=1.5,cex.axis=1.5,mgp=c(2.5,1,0))
myDiff$logp = -log10(myDiff$pvalue)
myDiff$logq = -log10(myDiff$qvalue)
plot(myDiff[,"meth.diff"],myDiff[,"logp"],,xlab="meth.diff",ylab="-log10(pvalue)",cex.lab =1.7,pch=20,cex=0.5,main="Volcano plot for Two Groups(p)",cex.main=2,cex.main=2,xlim=range(myDiff[,"meth.diff"]),ylim=range(myDiff$logp))
plot(myDiff[,"meth.diff"],myDiff[,"logq"],,xlab="meth.diff",ylab="-log10(qvalue)",cex.lab =1.7,pch=20,cex=0.5,main="Volcano plot for Two Groups(q)",cex.main=2,cex.main=2,xlim=range(myDiff[,"meth.diff"]),ylim=range(myDiff$logp))
hist(myDiff$meth.diff,col=mygrey,main="diff",xlab="diff",ylab="Frequency")
hist(myDiff$pvalue,col=mygrey,main="p-value",xlab="p-value",ylab="Frequency")

dev.off()
