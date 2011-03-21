# ===============================================================
# Pipeline for correlation coefficient analysis between histone modification and expression
# expression data: CAGE/DiTAGS/RNAseq
# memory-efficient version
# Author: Xianjun.Dong @ umassmed.edu
# Date: Thu Mar  3 15:41:35 EST 2011
# ===============================================================


# load("correlation2.RData")
# save(list = ls(all=TRUE), file = "correlation3.RData")

source("mylib.R")

# ===============================================================
# ================ setting for input data ======================
# ===============================================================

# --------for histone modification bins

JOBID="4325288"; # for 40+40 bins PClongestPC
#JOBID="4325408"; # for 40+1 bins PClongestPC


# ------- transcript quantity (TSS-based)
tab5rows=read.table("../data/Gencodev3c_hg19_TSS_norm_withheader_cage_pet_small_long_02092011.gff", nrow=150)
classes <- sapply(tab5rows, class)
data01 <- read.table("../data/Gencodev3c_hg19_TSS_norm_withheader_cage_pet_small_long_02092011.gff", colClasses = classes)
# defat the table
RNAseq=data01[,c(1,2,4,7, c(5:round(ncol(data01)/2))*2)]  # TO EDIT: mean is c(77:108)*3-1, max is c(77:108)*3
names(RNAseq)=c("chr","source","TSS","strand",as.vector(t(data01[1,c(5:round(ncol(data01)/2))*2-1])))
RNAseq=RNAseq[,c(colnames(RNAseq)[1:9], sort(colnames(RNAseq)[10:round(ncol(data01)/2)]))]

# RNAseq only
tab5rows=read.table("../data/gencode_v3c_hg19_tr_with115_cshl_long_quantif.gff2", nrow=150)
classes <- sapply(tab5rows, class)
data02 <- read.table("../data/gencode_v3c_hg19_tr_with115_cshl_long_quantif.gff2", colClasses = classes)
# defat the table
RNAseq=data02[,c(1,2,4,5,7, c(5:round(ncol(data02)/2))*2)]  # TO EDIT: mean is c(77:108)*3-1, max is c(77:108)*3
names(RNAseq)=c("chr","source","start","end","strand",as.vector(t(data02[1,c(5:round(ncol(data02)/2))*2-1])))
RNAseq=RNAseq[,c(colnames(RNAseq)[1:8], sort(colnames(RNAseq)[9:round(ncol(data02)/2)]))]


# ------- gene annotation [Optional]
genes=read.table("../data/gencode_v3c_hg19_transcripts.gtf.cpg.PCPC.Level12.L4k.longestTransID.tab", header = F)  # V14 = normalizedCpG
# ------- gene annotation [Optional]
genes=read.table("../data/gencode_v3c_hg19_transcripts.gtf.cpg.tab", header = T)  # V14 = normalizedCpG
genes0=genes[genes$level<3 & abs(genes$end-genes$start)>4000,]  # only level12 genes (>4k bp)
genes0=unique(data.frame(chr=genes0$chr, TSS=ifelse(genes0$strand=='+',genes0$start,genes0$end), gene_id=genes0$gene_id, normalizedCpG=genes0$normalizedCpG))
rownames(genes0) = paste(genes0$gene_id, genes0$TSS, sep='.')  # normalizedCpG is defined based on TSS, so ...

# for Bill
genes1=genes[genes$level<3 & abs(genes$end-genes$start)>4000,]  # only level12 genes (>4k bp)
RNAseq0= cbind(genes1[match(intersect(RNAseq$transcript_id, genes1$trans_id), genes1$trans_id),], RNAseq[match(intersect(RNAseq$transcript_id, genes1$trans_id), RNAseq$transcript_id),c(2,9:ncol(RNAseq))])
write.table(RNAseq0, file="gencode_v3c_hg19_transcripts.gtf.cpg.level12.L4k.tr_with115_cshl_long_quantif.tab", quote=F, sep="\t", row.names=F, col.names=T)

# ===============================================================
# ======== construct data frame for regression model ============
# ===============================================================

# -------- mean histone density for the most predictive bins

hist.expr=c()
celllines=c()  # to record which celllines have been checked
for (i in c(0:148))
{
             # read each histone modification signals
             filepath = paste("../data/jobid_", "/tab.output_", paste("_", i+1, sep=""), sep=JOBID)
             tab5rows=read.table(filepath, nrow=100)
             classes <- sapply(tab5rows, class)
             histdata=read.table(filepath, colClasses = classes, comment.char ="")
             rownames(histdata)=histdata[,1] # gene ID (or transID if geneid is not unique)

             expname = histdata[1,3]  # e.g. Gm12878.Control
             expnames= unlist(strsplit(as.character(expname), ".", fixed=T))  # e.g. c("Gm12878", "Control")
             expname2= paste(expnames[2], expnames[1], sep=".")  # e.g. Control.Gm12878

             histdata = histdata[,c(-1:-3)]  # only numeric data for histone singla. e.g. Nx80

             # intersection: subset of RNAseq (unique geneid)
             RNAseq0 = RNAseq[match(intersect(RNAseq$gene_id, rownames(histdata)), RNAseq$gene_id),]
             histdata = histdata[as.character(RNAseq0$gene_id),]

             #RNAseq0.cellline = RNAseq0[ ,grep(paste("RNAseq.PolyA\\+", expnames[1], sep="\\\n"), label.tr(colnames(RNAseq0)))]
             index=grep(paste(expnames[1], ".*1", sep=""), ignore.case =T, label.tr(colnames(RNAseq0)))  # only for replicate 1
             #index=grep(paste(expnames[1], ".*", sep=""), ignore.case =T, label.tr(colnames(RNAseq0)))
             if(length(index)<2) next;
             if(!(expnames[1] %in% celllines)) {celllines=c(celllines, expnames[1])}
             RNAseq0.cellline = RNAseq0[ ,index]  # RNAseq data for specific cellline

             histcors = apply(histdata, 2, function(x) histcor(hists=x, RNA=RNAseq0.cellline))
             rownames(histcors) = colnames(RNAseq0.cellline)

             # 1. aggregation plot
             draw_aggregation_plot(expname = expname, path= paste("../result", JOBID, "Agg.byBIN", sep="/"), histcors = histcors, histdata=histdata);

             # 2. heatmap
             draw_correlation_heatmap(expname = expname, path=paste("../result", JOBID, "Heatmap.byBIN", sep="/"), histcors = histcors, histdata=histdata);

             ## mean density for bins with |r|>0.2
             HISTCORS = abs(histcors)>0.2  # CxB where C: number of celllines; B: number of bins
             NUM = apply(HISTCORS, 1, sum)  # Cx1
             NUM[NUM==0]=1
             mean.density = (as.matrix(histdata) %*% t(HISTCORS)) / NUM  # NxB * BxC = NxC
             colnames(mean.density) = paste(expname, colnames(RNAseq0.cellline), sep="_")

             hist.expr=cbind(hist.expr, mean.density);
}

# ----- add expression data to the hist.expr matrix
for(ex in celllines){
             index=grep(paste(ex, ".*1", sep="\\."), ignore.case =T, label.tr(colnames(RNAseq0)))
             if(length(index)>=2) hist.expr=cbind(hist.expr, RNAseq0[ ,index])
}

# ===============================================================
# ========= correlation analysis: regression model  =============
# ===============================================================

# --- calculate the correlation R between expression and histone;

lrm.r.result=c();
for(ex in colnames(hist.expr)[grep("_", colnames(hist.expr), invert=T)]){
             if(ex == 'rp2w1') next;
             if(ex == 'dtkl1') next;

             dataset=hist.expr[, grep(ex, colnames(hist.expr))]
             dataset=dataset[apply(dataset,2,sum)>0]

             # take out H3K36me3 and H3K79me2
             #dataset=dataset[grep('H3K36me3|H3k79me2', colnames(dataset), ignore.case =T, invert=T)]

             if(ncol(dataset)==1) next;  # if only one column left

             lrm.r=c(unlist(strsplit(sub("\\n","\\.", label.tr(ex)),"\\.")))

             # add normalized.CpG column
             dataset = cbind(dataset[intersect(genes$V5, rownames(dataset)),], CpG=genes$V14[genes$V5 %in% intersect(genes$V5, rownames(dataset))])

             subsets = c('all', 'HCP','LCP');
             for(ss in subsets)
             {
                          dataset0=dataset;
                          if(ss == 'HCP') dataset0=dataset[dataset$CpG >0.4,]
                          if(ss == 'LCP') dataset0=dataset[dataset$CpG<=0.4,]

                          re=lrm(expr=ex, data=dataset0, path=paste("../result", JOBID, "lrm", sep="/"), filename=paste(gsub("\n", ".",label.tr(ex)), "lrmbestbin",ss,sep="."), filetype="png")

                          lmg = sort(re[[2]], decreasing=T)
                          top75 = lmg[1:min(which(sapply(c(1:length(lmg)), function(x) sum(lmg[1:x]))>=0.75))]
                          #top75.string = sapply(c(1:length(top75)), function(x) paste(names(top75)[x], " (", round(as.numeric(top75)[x]*100, 1), ")", sep=""))
                          #lrm.r=c(lrm.r,re[[1]], list(names(top75)), list(as.numeric(top75)*100))
                          lrm.r=c(lrm.r,re[[1]], paste(names(top75),sep="",collapse=";"), paste(round(as.numeric(top75)*100, 1),sep="",collapse=";"))

                          ## -- plot the relative importance for different expression subdivisions
                          #relimp.bin.by.expression(dataset0=dataset0, ex=ex, N=10, path=paste("../result", JOBID, "lrm", sep="/"), filename=paste("relimp",ex, ss,"png", sep="."))
             }
             lrm.r.result=rbind(lrm.r.result, lrm.r)
}
colnames(lrm.r.result)= c("Techique", "RNA_Extract", "Cell_Line", "Compartment", "SCC.all", "SCC.cv.all", "Top75_predictors.all", "Top75_predictors_relimp.all", "SCC.HCP", "SCC.cv.HCP", "Top75_predictors.HCP", "Top75_predictors_relimp.HCP", "SCC.LCP", "SCC.cv.LCP", "Top75_predictors.LCP", "Top75_predictors_relimp.LCP")
write.table(lrm.r.result,  file=paste("../result", JOBID, "lrm/correlation.result.txt", sep="/"), quote=F, sep="\t", row.names=F, col.names=T)
lrm.r.result = read.delim(paste("../result", JOBID, "lrm/correlation.result.txt", sep="/"))
lrm.r.result = lrm.r.result[order(-lrm.r.result$SCC.all),]

# plot the bestbin.spearman
subsets = c('all', 'HCP','LCP');
for(ss in subsets)
{
             lrm.r.sorted = lrm.r.result[,paste("SCC", ss,sep=".")]
             names = apply(t(lrm.r.result[,1:4]), 2, function(x) paste(x,collapse="."))

             png(paste("../result", JOBID, "lrm", paste("lrm.bestbin.spearman", ss, "png", sep="."), sep="/"), width=1000, height=1000)
             par(mar=c(5,14,2,2))
             barplot(lrm.r.sorted, names.arg=names, horiz=T, las=1, cex.names =1, main=paste("Correlation Coeffecient for Different Experiments\nlongest PC transcripts of", ss, "genes", sep=" "), xlab="Spearman's rho from LRM with best bins")
             abline(v=0.5, col='red')
             dev.off()
}

# ------------------------------ END

## do lrm...

data.0=data0[,c(1,2,4,7, c(5:146)*2)]
names(data.0)=c("chr","source","TSS","strand",as.vector(t(data0[1,c(5:146)*2-1])))
data.0=data.0[,c(colnames(data.0)[1:8], sort(colnames(data.0)[9:114]), sort(colnames(data.0)[115:ncol(data.0)]))]
colnames(data.0)[72]="dpgc2"


# with mean/max
tab5rows=read.table("../data/Gencodev3c_hg19_TSS_norm_withtrlist_with28cshlrnaseq.BroadTier1.jiebw.mean.max.gff", nrow=150)
classes <- sapply(tab5rows, class)
data01 <- read.table("../data/Gencodev3c_hg19_TSS_norm_withtrlist_with28cshlrnaseq.BroadTier1.jiebw.mean.max.gff", colClasses = classes)
# defat the table
data.1=data01[,c(1,2,4,7, c(5:114)*2, c(77:108)*3)]  # TO EDIT: mean is c(77:108)*3-1, max is c(77:108)*3
names(data.1)=c("chr","source","TSS","strand",as.vector(t(data01[1,c(5:114)*2-1])), as.vector(t(data01[1,c(77:108)*3-2])))
data.1=data.1[,c(colnames(data.1)[1:8], sort(colnames(data.1)[9:114]), sort(colnames(data.1)[115:ncol(data.1)]))]
colnames(data.1)[72]="dpgc2"

# new data with 35RNAseq and all Borad histmod data(tier1+tier2)
tab5rows=read.table("../data/Gencodev3c_hg19_TSS_norm_withheader_withtrlist_andbiotypes_with35cshlrnaseq.hismod.BroadAll.ucscBigWig.4kTSS.mean.max.gff", nrow=150)
classes <- sapply(tab5rows, class)
data02 <- read.table("../data/Gencodev3c_hg19_TSS_norm_withheader_withtrlist_andbiotypes_with35cshlrnaseq.hismod.BroadAll.ucscBigWig.4kTSS.mean.max.gff", colClasses = classes)
data.2=data02[,c(1,2,4,7, c(5:122)*2, c(83:208)*3)]  # TO EDIT: mean is c(83:208)*3, max is c(83:208)*3?1
names(data.2)=c("chr","source","TSS","strand",as.vector(t(data02[1,c(5:122)*2-1])), as.vector(t(data02[1,c(83:208)*3-1])))
data.2=data.2[,c(colnames(data.2)[1:9], sort(colnames(data.2)[10:122]), sort(colnames(data.2)[123:ncol(data.2)]))]

# erro in data.2

data = data.0

## ------------  correlation vs. expression levels ---------------
# for cpkc1
#expts = colnames(data)[grep("^r.?g.?2", colnames(data))]
expts = colnames(data)[grep("^[a-z].?[1|k|g].?[1-2]", colnames(data))]
R1=c()
#R2=c()
for(ex in expts)
{
             ex='rpkn2'
             # discard bottom 40% of CAGE data
             d1=data[do.call(order, list(data[ex], decreasing=T)),]
             d1=d1[d1[ex]>0,] # filter out those genes not expressed

             ###### bin by log2(expression)
             d1.min=log2(min(d1[ex]));
             d1.max=log2(max(d1[ex]));
             interval = (d1.max-d1.min)/10

             r1=c();r11=c()
             r2=c();r22=c()
             for(i in c(1:10)){
                          n=round(i*nrow(d1)/10) # each bin has same number of genes
                          d1i=d1[1:n,]
                          # bin by log2(expression)
                          #d1i=d1[log2(d1[ex]) >= (d1.max-i*interval) & log2(d1[ex]) <= (d1.max-(i-1)*interval), ]
                          re=lrm(expr=ex, data=d1i, onlyR=TRUE)
                          r1=c(r1, re[[1]])
                          re=re[[2]]
                          r11=rbind(r11, re[sort(names(re))])
                          #r2=c(r2, lrm(expr=ex, data=d1[i:nrow(d1),], onlyR=TRUE))
             }
             R1=cbind(R1, r1)

             ## -- plot the relative importance
             png(paste("relimp", ex,"png",sep="."), width=600, height=600)
             r11=t(r11)
             colnames(r11)=paste(c(0:9)*10, "-", c(1:10)*10, "%", sep="")
             # Expand right side of clipping rect to make room for the legend
             par(xpd=T, mar=par()$mar+c(0,0,0,7))
             cols=topo.colors(nrow(r11))
             pm = barplot(r11, col=cols, space=0.1, cex.axis=0.8, cex.names=0.6, las=1, cex.main=0.8,
                     main=paste("Relative importance of predictors for ",sub("\\n","\\.", label.tr(ex)),sep=""),
                     xlab="intervel of decreasing-sorted expression levels",
                     ylab="percentages of relative importance")
                     #legend.text=T,
                     #args.legend = list(x=11, y=1, cex=0.8, fill=sort(rainbow(ncol(r11)), decreasing=F)))
             lines(pm, r1, type='b', pch='*', col='red')
             # Place the legend at (6,30) using heat colors
             legend(11, 1, c(rev(rownames(r11)), "R(measured, predicted)"), cex=0.8, pt.cex=c(rep(2,nrow(r11)),0.6), col=c(rev(cols),'red'), lty=c(rep(-1,nrow(r11)),1), pch=c(rep(15,nrow(r11)),8));
             dev.off()

             # check top 10%, bottom 10% and middle 10%
             #d1.1=d1[1:round(nrow(d1)/10), ]
             #d1.2=d1[round(nrow(d1)*0.9):nrow(d1), ]
             #lrm(expr=ex, data=d1.1, onlyR=TRUE)
             #plot(r1, type='b',pch='*', xlab="Top % of lowly expressed genes", ylab="Correlation of linear regression model", main="Prediction accuracy Vs. expression levels")
}
colnames(R1)=expts
#colnames(R2)=expts
png("prediction.vs.expression.png", width=1000, height=800)
# plot
par(xpd=T, mar=par()$mar+c(4,0,0,12))
cols=sort(rainbow(ncol(R1)), decreasing=F)
plot(R1[,1], type='b',pch='*', ylim=c(min(R1),max(R1)), col=cols[1],
     xlab="Top x*10% of highly-expressed genes",
     ylab="Correlation of linear regression model",
     main="Prediction accuracy Vs. expression levels",
    )
for(i in c(2:ncol(R1))){ lines(c(1:nrow(R1)), R1[,i], type='b',pch='*', col=cols[i])}
#lines(c(1:100), r2, type='b',pch='*', col='gray')
#axis(3, labels=seq(100,0,-20), at=seq(0,100,20), line=1, tick=F, col.axis='gray')
legend(10.5, max(R1), sub("\\n","\\.", label.tr(colnames(R1))), cex=0.8, col=cols, lty=2, pch='*', bty='n')
dev.off()



# merge tran_biotype from data.2 to data.0
data = cbind(data, trbiotype=data.2$trbiotype)

# ---------- get normalized.CpG.content for each 3k region centered on TSS
cpg=read.table("../data/gencode_v3c_hg19_transcripts.cpg.bed", header = T)  # NB: this file is 0-based, while the Gencodev3c is 1-based.
x=cbind(data, TSS_id=paste(data$gene_id, data$TSS, sep="_"))
y=unique(as.data.frame(cbind(TSS_id=paste(cpg$genename, ifelse(cpg$strand=='-', cpg$end, cpg$start+1), sep="_"), normalized.CpG.content=cpg$normalized.CpG.content)))
# merge CpG to the big table
data=cbind(x[match(intersect(x$TSS_id, y$TSS_id),x$TSS_id),], CpG=y[match(intersect(x$TSS_id, y$TSS_id),y$TSS_id),2])
data$CpG=as.numeric(as.character(data$CpG))  # de-factor

png(filename="normalized.CpG.content.hg19.2.png", width=2400, height=800)
split.screen(c(1,3))
screen(1)
hist(data$CpG, breaks=190, xlim=c(0,1),xlab="normalized CpG content", main="all genes(hg19)")  # biomode
split.screen(c(2,1), screen=2)
screen(4)
hist(data$CpG[data$gene_biotype == "protein_coding"],breaks=190, xlim=c(0,1),xlab="normalized CpG content", main="protein_coding genes(hg19)")  # biomode
screen(5)
hist(data$CpG[data$gene_biotype != "protein_coding"],breaks=190, xlim=c(0,1),xlab="normalized CpG content", main="non_coding genes(hg19)")  # biomode
split.screen(c(2,1), screen=3)
screen(6)
hist(data$CpG[data$gene_biotype == "protein_coding" & grepl("protein_coding", data$trbiotype)],breaks=190, xlim=c(0,1),xlab="normalized CpG content", main="protein_coding transcripts \nof protein_coding genes(hg19)")  # biomode
screen(7)
hist(data$CpG[data$gene_biotype == "protein_coding" & !grepl("protein_coding", data$trbiotype)],breaks=190, xlim=c(0,1),xlab="normalized CpG content", main="non_coding transcripts \nof protein_coding genes(hg19)")  # biomode
close.screen(all = TRUE)
dev.off()

# ---------- filter1: only protein-coding genes
data.nc = data[data$gene_biotype != "protein_coding",]
data.pc = data[data$gene_biotype == "protein_coding", ]
data.pc.pc = data[data$gene_biotype == "protein_coding" & grepl("protein_coding", data$trbiotype),]
data.pc.nc = data[data$gene_biotype == "protein_coding" & !grepl("protein_coding", data$trbiotype),]

data=data.pc.nc

# save the table
write.table(data, file="../data/Gencodev3c_hg19_TSS_norm_withheader_withtrlist_andbiotypes_with35cshlrnaseq.hismod.Broadtier1.jiebw.4kTSS.mean.tab", quote=F, row.names=F, sep="\t")

# subset of data

data2=data[,c(grep("^cpk", colnames(data)), grep("^dpk", colnames(data)), grep("^rpk", colnames(data)))]  # CAGE ployA+ K562, vs. DiTAGS PolyA+ K562, vs. RNAseq polyA+ K562
corplot(data2, grep_expre="", grep_histmod="", filename="CAGE.vs.DiTAGS.vs.RNAseq__PolyA+.K562.jiebw.4kTSS.mean",filetype='png')

data3=data[,c(grep("^cp", colnames(data)))]  # CAGE Expression data btw different cell lines
corplot(data3, grep_expre="1", grep_histmod="", filename="Different.Celllines__CAGE.PolyA+.jiebw.4kTSS.mean", filetype="png")

data1=data[,c(grep("^cpk", colnames(data)), grep("K562", colnames(data)))]  # CAGE ployA+ K562, vs. K562 hist
corplot(data1, grep_expre="2", grep_histmod="H", filename="CAGE2.vs.hist__PolyA+.K562.jiebw.4kTSS.mean", filetype="png")

data7=data[,c(grep("^rpk", colnames(data)), grep("K562", colnames(data)))]  # RNAseq ployA+ vs. hist, K562
corplot(data7, grep_expre="", grep_histmod="H", filename="RNAseq.vs.hist__PolyA+.K562.jiebw.4kTSS.mean", filetype="png")
data8=data[,c(grep("^rpg", colnames(data)), grep("Gm", colnames(data)))]  # RNAseq ployA+ vs. hist, Gm
corplot(data8, grep_expre="2", grep_histmod="H", filename="RNAseq.vs.hist__PolyA+.Gm21878.jiebw.4kTSS.mean", filetype="png")
data9=data[,c(grep("^rp1", colnames(data)), grep("H1hesc", colnames(data)))]  # RNAseq ployA+ vs. hist, H1hesc
corplot(data8, grep_expre="", grep_histmod="H", filename="RNAseq.vs.hist__PolyA+.H1hesc.jiebw.4kTSS.mean", filetype="png")


data5=data[,c(grep("^cpg", colnames(data)), grep("Gm", colnames(data)))]  # CAGE polyA+ Gm12878, vs. Gm12878 hist
corplot(data5, grep_expre="2", grep_histmod="H", filename="CAGE.vs.hist__PolyA+.Gm12878.jiebw.4kTSS.mean", filetype="png")

data6=data[,c(grep("^cp1", colnames(data)), grep("H1hesc", colnames(data)))]  # CAGE polyA+ H1hesc, vs. H1hesc hist
corplot(data6, grep_expre="", grep_histmod="H", filename="CAGE.vs.hist__PolyA+.H1hesc.jiebw.4kTSS.mean", filetype="png")

# ---------- estimation of pseudocount
# 1. split the dataset into two subsets randomly, one for estimation of pseudocount, the other for prediction
n = nrow(data)
n1 = sample(c(1:n), round(n/3), replace=F)
n2 = setdiff(c(1:n),n1)
# 2. estimate the pseudocount for each modification
data.n1=data[n1,c(grep("cpkc1", colnames(data)),grep("K562", colnames(data)))]  # test for cpkc1
data.n2=data[n2, c(grep("cpkc1", colnames(data)),grep("K562", colnames(data)))]

# filter out data points with 0 expression value
data.n1=data.n1[data.n1[,1]>0,]  # around 50% are not expressed...
#data.n1[data.n1[,1]==0,1]=1  # add pseudocount of 1 for expression
#data.n2[data.n2[,1]==0,1]=1  # add pseudocount of 1 for expression

for(i in c(2:ncol(data.n1))){
             m = c(1:max(data.n1[,i]))
             a = m[which.max(sapply(m, function(x) cor(log2(data.n1[,1]), log2(data.n1[,i]+x))))]

             # add the pseudocount to correpsonding hist.mod in dataset2
             data.n2[,i] = data.n2[,i]+a
}

# do regression for dataset2
lrm(expr="cpkc1", data=data.n2)

## conclusion: R is 0.6, same as the one without adding pseudocount...
## but the relative important modification are bit different..


# -----------------linear regression model
tech = unlist(strsplit("cdr", ""))  # cage, ditag, rna-seq
#cellline = unlist(strsplit("kg1uehmnp", ""))  # k = "K562",g = "Gm12878",'1' = "H1hesc",u = "HUVEC",e = "Hela-S3",h = "HepG2",m = "MCF7",n = "NHEK",p = "Pros",
cellline = unlist(strsplit("kg1", ""))  # k = "K562",g = "Gm12878",'1' = "H1hesc",u = "HUVEC",e = "Hela-S3",h = "HepG2",m = "MCF7",n = "NHEK",p = "Pros",
results=c()
for (i in tech){
             for (j in cellline){
                          reg.exp.code = paste("^", i, ".?", j, sep="");
                          exp.tech.cellline = colnames(data)[intersect(grep(reg.exp.code, colnames(data)), grep("\\.", colnames(data), invert=T))]
                          res = c(tech.tr(i), cellline.tr(j))
                          for (ex in exp.tech.cellline)
                          {
                                       #ex='cpkc1'
                                       # individual R for top predictors

                                       # R for lm without intercross

                                       # R for lm with intercross

                                       r = lrm(expr=ex, data=data, filename=paste("lrm", sub("\\n","\\.", label.tr(ex)), sep="."), filetype="png")
                                       res1=c(res, unlist(strsplit(sub("\\n","\\.", label.tr(ex)),"\\."))[c(2,4)], r)
                                       r = lrm(expr=ex, data=data[data$CpG > 0.4,], filename=paste("lrm", sub("\\n","\\.", label.tr(ex)), "HCP", sep="."), filetype="png")
                                       res1=c(res1, r)
                                       r = lrm(expr=ex, data=data[data$CpG <= 0.4,], filename=paste("lrm", sub("\\n","\\.", label.tr(ex)), "LCP", sep="."), filetype="png")
                                       res1=c(res1, r)

                                       print(res1, quote=F)

                                       results=rbind(results, res1)
                          }
             }
}
# write results to table
write.table(results, file="correlation.result.max.txt", quote=F, sep="\t", row.names=F, col.names=F)

## -------- statistics for the correlation --------
results = read.table(file="pc.pc/correlation.result.max.txt", quote="", sep="\t", row.names=NULL)
results.ordered = results[with(results,order(V5, V7, V9)), ]
results.ordered = results.ordered[grep("^Nuc[0-9]|Wcell|Cy", results.ordered$V4),] # only Cy, Nuc, Wcell
results.ordered = results.ordered[grep("RNAseq", results.ordered$V1, invert=T),] # only CAGE, DiTag
plot(x=c(1:nrow(results.ordered)), results.ordered$V5, ylim=c(min(results.ordered[, c(5, 7, 9)]), max(results.ordered[, c(5, 7, 9)])), col='blue', type='b', pch=19, xlab="Experiment ID [ordered by R(all genes)]", ylab="R(measured, predicted)")
lines(x=c(1:nrow(results.ordered)), results.ordered$V7, col='red', type='b', pch=sapply(results.ordered$V4, function(x) ifelse(grepl("Wcell", x), 8, ifelse(grepl("Nuc", x), 1,19))))
lines(x=c(1:nrow(results.ordered)), results.ordered$V9, col='pink', type='b', pch=sapply(results.ordered$V4, function(x) ifelse(grepl("Wcell", x), 8, ifelse(grepl("Nuc", x), 1,19))))
legend("bottomright",
       c("all genes", "HCP genes","LCP genes", "Whole cell", "Nucleus", "Cytosol"),
       col=c("blue","red", "pink", "black","black","black"),
       lty=c(1,1,1,0,0,0),
       pch=c(-1,-1,-1,8,1,19))

# ---------- apply coefficients from cellline A to cellline B
# cpkc2/n2/w2, cpgc2/n2/w2, cp1c2/n2/w2
png(filename="k562.to.others.png",width=1000, height=500)
TAB=c()
r=lrm.AtoB(expA="cpkc2", expB=c("cpgc2", "cp1c2"), data=data)
TAB=cbind(Cytosol=r, Nucleus=lrm.AtoB(expA="cpkn2", expB=c("cpgn2", "cp1n2"), data=data))
TAB=cbind(TAB, Whole.Cell=lrm.AtoB(expA="cpkw2", expB=c("cpgw2", "cp1w2"), data=data))
rownames(TAB)=c("K562", "K562->Gm12878", "K562->H1hesc")
mp <- barplot(TAB, beside = TRUE, axisnames = FALSE, ylab="Cor(measurement, prediction)", main="Cell line specificity of K562")
# Add the individual bar labels
mtext(1, at = mp, text = rownames(TAB),line = 0, cex = 0.8)
mtext(1, at = colMeans(mp), text = colnames(TAB), line = 2)
dev.off()

# cpkc2/n2/w2, cpgc2/n2/w2, cp1c2/n2/w2
png(filename="k562.to.others.lcp.png", width=1000, height=500)
data.lcp=data[data$CpG <= 0.4, ]
TAB=c()
r=lrm.AtoB(expA="cpkc2", expB=c("cpgc2", "cp1c2"), data=data.lcp)
TAB=cbind(Cytosol=r, Nucleus=lrm.AtoB(expA="cpkn2", expB=c("cpgn2", "cp1n2"), data=data.lcp))
TAB=cbind(TAB, Whole.Cell=lrm.AtoB(expA="cpkw2", expB=c("cpgw2", "cp1w2"), data=data.lcp))
rownames(TAB)=c("K562", "K562->Gm12878", "K562->H1hesc")
mp <- barplot(TAB, beside = TRUE, axisnames = FALSE, ylab="Cor(measurement, prediction)", main="Cell line specificity of K562 (LCP genes)")
# Add the individual bar labels
mtext(1, at = mp, text = rownames(TAB),line = 0, cex = 0.8)
mtext(1, at = colMeans(mp), text = colnames(TAB), line = 2)
dev.off()

# cpkc2/n2/w2, cpgc2/n2/w2, cp1c2/n2/w2
png(filename="k562.to.others.hcp.png", width=1000, height=500)
data.hcp=data[data$CpG > 0.4, ]
TAB=c()
r=lrm.AtoB(expA="cpkc2", expB=c("cpgc2", "cp1c2"), data=data.hcp)
TAB=cbind(Cytosol=r, Nucleus=lrm.AtoB(expA="cpkn2", expB=c("cpgn2", "cp1n2"), data=data.hcp))
TAB=cbind(TAB, Whole.Cell=lrm.AtoB(expA="cpkw2", expB=c("cpgw2", "cp1w2"), data=data.hcp))
rownames(TAB)=c("K562", "K562->Gm12878", "K562->H1hesc")
mp <- barplot(TAB, beside = TRUE, axisnames = FALSE, ylab="Cor(measurement, prediction)", main="Cell line specificity of K562 (HCP genes)")
# Add the individual bar labels
mtext(1, at = mp, text = rownames(TAB),line = 0, cex = 0.8)
mtext(1, at = colMeans(mp), text = colnames(TAB), line = 2)
dev.off()
# -----------
# cpkc2/n2/w2, cpgc2/n2/w2, cp1c2/n2/w2
png(filename="GM.to.others.png",width=1000, height=500)
TAB=c()
r=lrm.AtoB(expA="cpgc2", expB=c("cpkc2", "cp1c2"), data=data)
TAB=cbind(Cytosol=r, Nucleus=lrm.AtoB(expA="cpgn2", expB=c("cpkn2", "cp1n2"), data=data))
TAB=cbind(TAB, Whole.Cell=lrm.AtoB(expA="cpgw2", expB=c("cpkw2", "cp1w2"), data=data))
rownames(TAB)=c("Gm12878", "Gm->K562", "Gm->H1hesc")
mp <- barplot(TAB, beside = TRUE, axisnames = FALSE, ylab="Cor(measurement, prediction)", main="Cell line specificity of Gm12878")
# Add the individual bar labels
mtext(1, at = mp, text = rownames(TAB),line = 0, cex = 0.8)
mtext(1, at = colMeans(mp), text = colnames(TAB), line = 2)
dev.off()

png(filename="GM.to.others.lcp.png",width=1000, height=500)
data.lcp=data[data$CpG<0.4, ]
TAB=c()
r=lrm.AtoB(expA="cpgc2", expB=c("cpkc2", "cp1c2"), data=data.lcp)
TAB=cbind(Cytosol=r, Nucleus=lrm.AtoB(expA="cpgn2", expB=c("cpkn2", "cp1n2"), data=data.lcp))
TAB=cbind(TAB, Whole.Cell=lrm.AtoB(expA="cpgw2", expB=c("cpkw2", "cp1w2"), data=data.lcp))
rownames(TAB)=c("Gm12878", "Gm->K562", "Gm->H1hesc")
mp <- barplot(TAB, beside = TRUE, axisnames = FALSE, ylab="Cor(measurement, prediction)", main="Cell line specificity of Gm12878 (LCP genes)")
# Add the individual bar labels
mtext(1, at = mp, text = rownames(TAB),line = 0, cex = 0.8)
mtext(1, at = colMeans(mp), text = colnames(TAB), line = 2)
dev.off()

png(filename="GM.to.others.hcp.png",width=1000, height=500)
data.hcp=data[data$CpG>=0.4, ]
TAB=c()
r=lrm.AtoB(expA="cpgc2", expB=c("cpkc2", "cp1c2"), data=data.hcp)
TAB=cbind(Cytosol=r, Nucleus=lrm.AtoB(expA="cpgn2", expB=c("cpkn2", "cp1n2"), data=data.hcp))
TAB=cbind(TAB, Whole.Cell=lrm.AtoB(expA="cpgw2", expB=c("cpkw2", "cp1w2"), data=data.hcp))
rownames(TAB)=c("Gm12878", "Gm->K562", "Gm->H1hesc")
mp <- barplot(TAB, beside = TRUE, axisnames = FALSE, ylab="Cor(measurement, prediction)", main="Cell line specificity of Gm12878 (HCP genes)")
# Add the individual bar labels
mtext(1, at = mp, text = rownames(TAB),line = 0, cex = 0.8)
mtext(1, at = colMeans(mp), text = colnames(TAB), line = 2)
dev.off()


## ---- ncRNA vs. pcRNA



k562 = colnames(data)[intersect(grep("^c.?k", colnames(data)), grep("\\.", colnames(data), invert=T))]
for (ex in k562) {
             lrm(expr=ex, data=data, filename=paste("lrm", sub("\\n","\\.", label.tr(ex)), sep="."), filetype="png")
             lrm(expr=ex, data=data[data$CpG > 0.4,], filename=paste("lrm", sub("\\n","\\.", label.tr(ex)), "HCP", sep="."), filetype="png")
             lrm(expr=ex, data=data[data$CpG <=0.4,], filename=paste("lrm", sub("\\n","\\.", label.tr(ex)), "LCP", sep="."), filetype="png")
}
# for Gm CAGE
Gm = colnames(data)[intersect(grep("^c.?g", colnames(data)), grep("\\.", colnames(data), invert=T))]
for (ex in Gm) {
             lrm(expr=ex, data=data, filename=paste("lrm", sub("\\n","\\.", label.tr(ex)), sep="."), filetype="png")
             lrm(expr=ex, data=data[data$CpG > 0.4,], filename=paste("lrm", sub("\\n","\\.", label.tr(ex)), "HCP", sep="."), filetype="png")
             lrm(expr=ex, data=data[data$CpG <=0.4,], filename=paste("lrm", sub("\\n","\\.", label.tr(ex)), "LCP", sep="."), filetype="png")
}
# for H1Esc CAGE
H2esc = colnames(data)[intersect(grep("^c.?1", colnames(data)), grep("\\.", colnames(data), invert=T))]
for (ex in H2esc) {
             lrm(expr=ex, data=data, filename=paste("lrm", sub("\\n","\\.", label.tr(ex)), sep="."))
             lrm(expr=ex, data=data[data$CpG > 0.4,], filename=paste("lrm", sub("\\n","\\.", label.tr(ex)), "HCP", sep="."), filetype="png")
             lrm(expr=ex, data=data[data$CpG <=0.4,], filename=paste("lrm", sub("\\n","\\.", label.tr(ex)), "LCP", sep="."), filetype="png")
}

pairs(cbind(data1$cpkc1, fitted(fit1)),
      labels=c(label.tr("cpkc1"), "Predicted values (log2)"),
      lower.panel=function(...) {par(new=TRUE);smoothScatter(...)},
      upper.panel=panel.cor)
# diagnostic plots
layout(matrix(c(1,2,3,4),2,2)) # optional 4 graphs/page
plot(fit1)
# Calculate Relative Importance for Each Predictor
library(relaimpo)
re=calc.relimp(fit1,type=c("lmg"),rela=TRUE)
# create extra margin room on the right for an axis
par(mar=c(5, 8, 4, 4) + 0.1)
barplot(sort(re$lmg, decreasing=F), horiz=T, las=1, cex.names =0.7, xlab="relative importance (% to R^2)")

# Bootstrap Measures of Relative Importance (1000 samples)
boot <- boot.relimp(fit, b = 1000, type = c("lmg",
  "last", "first", "pratt"), rank = TRUE,
  diff = TRUE, rela = TRUE)
booteval.relimp(boot) # print result
plot(booteval.relimp(boot,sort=TRUE)) # plot result


# ---------- heatmap of RNA-seq
library(gplots)
#pdf(file="RNAseq.K562.RPKM.heatmap.2.pdf", height =10, width = 8)
heatmap.2(as.matrix(rnaseqlog2_hist[3000:5000,grep("Bio", colnames(rnaseqlog2_hist))]),
          #Rowv = as.vector(sort(cl$cluster)), Colv=T, #dendrogram = "both",
          col=greenred, #colorpanel(8,"green","white","red"),
          trace='none', margins=c(8, 8),
          key=T, keysize = 1, scale='none', cexRow=0.5, density.info= 'none',
          main = "RNAseq.K562.RPKM.heatmap",
           xlab = "Experiments (Cell line + Compartments)",
           ylab = "Transcript ID")
#dev.off()
