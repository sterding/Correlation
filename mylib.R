library(relaimpo)  # if first time
library(e1071) # interface for SVMLIB

get_histcors <- function (histdata=histdata, RNAseq0.cellline=RNAseq0.cellline)
{
             histcors=c();
             for(j in seq(ncol(histdata))){
                          r1=c()
                          for(k in seq(ncol(RNAseq0.cellline))){
                                       y=RNAseq0.cellline[,k]; x=histdata[,j];
                                       r1=c(r1,round(cor(log2(x[x>0 & y>0]), log2(y[x>0 & y>0]), method="spearman"), 3))
                          }
                          histcors=cbind(histcors, r1)
             }
             rownames(histcors) = colnames(RNAseq0.cellline)
             colnames(histcors) = paste("bin", seq(ncol(histdata)), sep="")
             return(histcors)
}

histcor <- function (hists=x, RNA=data0.cellline)
{
             r.array=c();
             for(i in c(1:ncol(RNA)))
             {
                          dat=data.frame(expr=RNA[,i], hist=hists)

                          #dat0=dat[dat$expr>0,] # filter
                          #dat0$expr=log2(dat0$expr) # log2
                          # both should be log(x)
                          dat=dat[apply(dat>0, 1, sum)==ncol(dat),]
                          dat=log2(dat)

                          # # linear regression model
                          #fit <- lm(expr ~ hist, data=dat0)  # full model
                          #pred = fitted(fit)
                          #r=round(cor(pred, dat0$expr, method="pearson"), 3)  # godness of fit
                          r2=round(cor(dat$expr, dat$hist, method="spearman"), 3)  # Pearson correalation coefficient (linear dependance), different from above R
                          r.array=c(r.array, r2)

                          ## svm regression model
                          #obj = tune.svm(expr ~ hist, data = iris, gamma = 2^(-1:1), cost = 2^(2:4))  # looking for best model
                          #pred = fitted(obj$best.model)

             }
             return(r.array)
}

draw_aggregation_plot <- function (expname = expname, path=".", histcors = histcors, histdata=histdata)
{
             # transfer range of histdata to range of histcors (e.g. [-1, 1])
             X=range(histcors)
             Y=range(apply(histdata, 2, mean))
             a=(Y[2]-Y[1])/(X[2]-X[1])
             b=(X[2]*Y[1]-X[1]*Y[2])/(X[2]-X[1])
             histcors2=histcors*a+b

             if (file.exists(path) == FALSE) dir.create(path, showWarnings=F, recursive=T)

             ## 1. aggregation plot
             ##  path= "../result/Agg.byBIN"
             png(filename=paste(path, paste(expname,"png", sep="."), sep="/"))
             # create extra margin room on the right for an axis
             par(mar=c(5, 4, 4, 5) + 0.1)

             x=seq(ncol(histcors))

             plot(x, apply(histdata, 2, mean),
                  type='n', pch=20, main=expname,
                  xaxt="n", ylab="mean density", xlab="")


             abline(v=20, col="gray", lty=3)
             abline(v=60, col="gray", lty=3)


             # http://www.statmethods.net/advgraphs/axes.html

             ## add x vs. 1/x
             cols=rainbow(nrow(histcors2))
             for(l in c(1:nrow(histcors2))){
                          lines(x, histcors2[l,], type="b", pch=20, col=cols[l], lty=3)
             }
             #apply(histcors2, 1, function(x) lines(c(1:80), x, type="b", pch=20, col="lightblue", lty=3))
             #lines(c(1:80), histcors2[1,], type="b", pch=20, col="blue", lty=2)

             lines(x, apply(histdata, 2, mean),type='o', pch=20)

             ## draw an axis on the bottom
             if(ncol(histcors)==80) axis(1, at=20*c(0:4),labels=c("-2kb", "TSS", "2kb...-2kb", "TTS", "2kb"))
             if(ncol(histcors)==41) axis(1, at=c(0,10,20,30,40,41),labels=c("-2k", "-1k", "TSS", "1k", "2kb", "TTS"))
             #
             ## draw an axis on the right, with smaller text and ticks
             axis(4, at=seq(range(histcors2)[1], range(histcors2)[2], 1), labels=round(seq(range(histcors)[1], range(histcors)[2], 1/a), digits=2), col.axis="blue", las=2, cex.axis=0.7, tck=-.01)
             #
             ## add a title for the right axis
             mtext("correlation coefficient (Spearman's rho)", side=4, line=3, cex.lab=1,las=0, col="blue")
             legend("topright", gsub("\n", ".", label.tr(rownames(histcors))), col=cols, pch=20, lty=3, cex=0.7)
             #
             ## add a main title and bottom and left axis labels
             #title(expname, ylab="mean density")

             dev.off()
}

draw_correlation_heatmap <- function (expname = expname, path=".", histcors = histcors, histdata=histdata)
{
             if (file.exists(path) == FALSE) dir.create(path, showWarnings=F, recursive=T)

             png(paste(path, paste(expname,"png", sep="."), sep="/"), height=850, width=950)

             #png(paste("../result/Heatmap.byBIN",expname,"png", sep="."), height=850, width=950)

             #par(oma=c(2,2,2,2))
             layout(matrix(seq(4), nrow=2, ncol=2, byrow=T), widths=c(1,10), heights=c(10,1), TRUE)

             t = as.matrix(histcors)
             collist <- c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#F7F7F7","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F")
             ColorRamp<-colorRampPalette(collist, bias=1)(10000)
             #Color <- ColorRamp  # use all colors
             #Color <- ColorRamp[round((0.5-max(abs(t))/2)*10000) : round((0.5+max(abs(t))/2)*10000)]  # only use the corresponding part [-max, +max]
             Color <- ColorRamp[round(1+(0.5+min(t)/2)*10000) : round((0.5+max(t)/2)*10000)]  # only use the corresponding part [-min, +max]
             #ColorLevels <- seq(min(t), max(t), length=length(Color))
             ColorLevels <- seq(-1, 1, length=length(ColorRamp))

              # Color Scale
             par(mar = c(1,4,3,1))
             image(1, ColorLevels,
                  matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
                  col=ColorRamp,
                  xlab="",ylab="Correlation Coefficient (Spearman's rho)",
                  xaxt="n")
             axis(LEFT<-2, at=range(t), labels=c(paste("min", min(t), sep="\n"), paste("max", max(t), sep="\n")), cex.axis=1, col.axis="blue", col.ticks="blue", lwd=0, lwd.ticks=4)

             # correlation heatmap
             par(mar=c(1,1,3,9))
             image(1:ncol(t), 1:nrow(t), t(t) , axes=FALSE, xlab="", ylab="", col=Color, main=expname)
             #axis(BELOW<-1, at=1:ncol(t), labels=colnames(t), las=3, cex.axis=1)
             axis(RIGHT<-4, at=1:nrow(t), labels=gsub("\n", ".", label.tr(rownames(histcors))), las= HORIZONTAL<-1, cex.axis=0.8)

             # Color Scale
             t=matrix(apply(histdata, 2, mean), byrow=T, nrow=1)
             ColorRamp<-rev(gray(1:1000 /1000))
             Color <- ColorRamp#[round(1+(0.5+min(t)/2)*10000) : round((0.5+max(t)/2)*10000)]  # only use the corresponding part [-min, +max]
             ColorLevels <- seq(min(t), max(t), length=length(Color))
             par(mar = c(3,4,1,1))
             image(1, ColorLevels,
                  matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
                  xlab="",ylab="Density",
                  xaxt="n")

             # correlation heatmap
             par(mar=c(3,1,1,9))
             image(1:ncol(t), 1:nrow(t), t(t) , axes=FALSE, xlab="", ylab="", col=Color, main="")
             if(ncol(t)==41) axis(1, at=c(1,10,20,30,40,41),labels=c("-2k", "-1k", "TSS", "1k", "2kb", "TTS"), las=HORIZONTAL<-1, cex.axis=1)
             if(ncol(t)==80) axis(BELOW<-1, at=seq(1, ncol(t), length.out=5), labels=c("-2k","TSS","2k..-2k","TTS","2k"), las=HORIZONTAL<-1, cex.axis=1)
             axis(RIGHT<-4, at=1:nrow(t), labels=expname, las= HORIZONTAL<-1, cex.axis=0.8)
             dev.off()
}

# ---------- linear regression model
lrm <- function(expr="cpkc1", dataset=dataset0, path="./", filename="", filetype="", onlyR=FALSE, histoneOnly=FALSE)
{
             # for test
             #expr=ex; dataset=data0; path=paste("."); filename=paste(gsub("\n", ".",label.tr(ex)), "lrmbestbin",sep="."); filetype="png"; onlyR=F;

             # -----------------------------------------------------
             # ---------- data preparation (cleanup, log etc.)
             # -----------------------------------------------------

             cellline = cellline.tr(substr(expr,3,3))  # k-->K562
             dataset=dataset[,c(which(expr==colnames(dataset)), grep(paste("^", cellline, sep=""), colnames(dataset)))]  # include all marks for expr
             # if only include histone marks
             if(histoneOnly==TRUE)  dataset=dataset[,c(which(expr==colnames(dataset)), grep(paste("^", cellline, "\\.H", sep=""), colnames(dataset)))]  # only include hist mod

             dataset=dataset[apply(dataset>0, 1, sum)==ncol(dataset),]  # only rows with all non-zero values. TODO: offset with max correlation
             #dataset=dataset[apply(dataset==0, 1, sum)<ncol(dataset),]  # allow some (not all) columns to be zero
             #dataset[dataset[,expr]<0.0001,expr]=0.0001

             dataset=log2(dataset)  # both histone and RNAseq signal are normal distribution after log(x)

             # -----------------------------------------------------
             #----------- linear regression model
             # -----------------------------------------------------

             ## Create a formula for a model with a large number of variables:
             xnam <- colnames(dataset)[grep("\\.", colnames(dataset))]
             # without intercross
             fmla <- as.formula(paste(expr, "~", paste(xnam, collapse= "+")))

             # with all first-order effects and interactions up to the nth order, where n is given by ( )^n:
             # y= bo+ b1*A + b2*B + b3*C + b4*AB + b5*AC + b6*BC
             # [REF: http://www.montefiore.ulg.ac.be/~kvansteen/GBIO0009-1/ac20092010/Class8/Using%20R%20for%20linear%20regression.pdf]
             # fmla <- as.formula(paste(expr, "~", "(", paste(xnam, collapse= "+"), ")^2"))

             fit <- lm(fmla, data=dataset)  # full model

             # # svm
             # fit1 = svm(fmla, data=dataset, kernel='linear', cross=10)
            # Squared correlation coefficient (of the predicted and the true values of the dependent variable)
            # fit1$scorrcoef

             mear = dataset[,expr]
             pred = fitted(fit)
             #Spearman's rho statistic is used to estimate a rank-based measure of association,
             # which is more robust and have been recommended if the data do not necessarily come from a bivariate normal distribution.
             #r0=round(cor(pred, mear, method='spearman'), 3)
             r0=round(cor(pred, mear, method='pearson'), 3)

             # -----------------------------------------------------
             # ----------- cross-validation
             # -----------------------------------------------------

            library(bootstrap)
            # define functions
            theta.fit <- function(x,y){lsfit(x,y)}
            theta.predict <- function(fit,x){cbind(1,x)%*%fit$coef}

            # matrix of predictors
            X <- as.matrix(dataset[xnam])
            # vector of predicted values
            y <- as.matrix(dataset[expr])

            results <- crossval(X,y,theta.fit,theta.predict,ngroup=3)  # 3-fold cross-validataion
            #r1=round(cor(mear, results$cv.fit, method='spearman'), 3)
            r1=round(cor(mear, results$cv.fit, method='pearson'), 3)

             # -----------------------------------------------------
             # ----------- relative importance
             # -----------------------------------------------------

            if(length(fit$coefficients)==2) {
                lmg=fit$coefficients[2]
                lmg[1]=sign(fit$coefficients[2])  # -1 or 1, 100%
                cof=fit$coefficients[2]
            }
            else if(nrow(dataset)< (ncol(dataset)+4))
            {
                lmg=c(too.few.obs=1)
                cof=c(too.few.obs=0)
            }
            else{
                re=calc.relimp(fit,type=c("lmg"),rela=TRUE)

                # add sign
                cof=(fit$coefficients)[-1]
                signs = sign(cof)
                lmg=sort(re$lmg)
                cof = cof[names(lmg)]
                lmg=lmg*signs[names(lmg)]

            }
             names(lmg) = as.vector(sapply(names(lmg), function(x) sub("^.*\\.(.*)_.*","\\1",x)))

             #lmg = sort(re$lmg, decreasing=T)
             #lmg.top75 = lmg[1:min(which(sapply(c(1:length(lmg)), function(x) sum(lmg[1:x]))>=0.75))]

             r = list(c(r0,r1), lmg)

             if(onlyR==TRUE) return(r)

             # -----------------------------------------------------
             # ---------- Otherwise, plot the figure
             # -----------------------------------------------------
             if (file.exists(path) == FALSE) dir.create(path, showWarnings=F, recursive=T)

             if(filename==""){
                          filename = sub("\\n","\\.", label.tr(expr))
             }
             if(filetype=="pdf"){
                          pdf(file=paste(path, paste(filename,filetype,sep="."), sep="/"), title=filename, height =8, width = 16)
             }
             if(filetype=="png"){
                          png(filename=paste(path, paste(filename,filetype,sep="."), sep="/"), height=800, width = 1600)
             }
             # plot 1 : prediction vs. measurement
             split.screen(c(1,2))
             screen(1)
             par(mar=c(5, 5, 4, 4) + 0.1)
             smoothScatter(pred, mear, xlab="predicted log2(expression)", ylab="measured log2(expression)", main=filename, cex.lab=1.5)
             abline(lm(mear ~ pred), col="red")
             legend("topleft", paste("Pearson's R = ", r0, "\nR for 3-fold cross-validation = ", r1), bty="n", cex=1.5)  #no border

             #plot 2: relative important predictors
             screen(2)
             # library(relaimpo)  # if first time
             # create extra margin room on the right for an axis
             par(mar=c(5, 12, 4, 4) + 0.1)
             barplot(abs(lmg*100), horiz=T, las=1, border =3-sign(lmg), cex.names =1.8, cex.lab =1.5, xlab=expression(paste("Relative importance (% to ",R^2,")")))
             text(1, seq(length(cof))*1.2-0.5, round(cof, 2))
             legend("bottomright", c("positive coefficient (+)", "negative coefficient (-)"), cex=2, col=c(2, 4), bty="n", pch=c(22,22), pt.bg = "gray")
             close.screen(all = TRUE)

             if(filetype!=""){
                          dev.off()
             }

             return(r)
}

relimp.bin.by.expression <- function(dataset0=dataset0, ex=ex, N=10, path=paste("../result", JOBID, "lrm", sep="/"), filename=paste("relimp",ex, ss,"png", sep="."))
{
            #N=10 # cut intervals
            #path=
            #filename=
            ind = intersect(grep(ex,colnames(dataset0)),grep("\\.",colnames(dataset0),invert=T))
            dataset0=dataset0[order(dataset0[,ind], decreasing=T),]  # order dataset by expression value
            expr.range = range(dataset0[,ind])
            if(expr.range[1]==0) expr.range[1]=0.00001  # pseudo to avoid log(0)
            expr.range = log2(expr.range)

            bin.lrm=c();
            for(k in c(1:N)){
                         #d0=dataset0[round((k-1)*nrow(dataset0)/N):round(k*nrow(dataset0)/N),]   # each bin has same number of genes
                         expr.range.min = expr.range[1]+(k-1)*(expr.range[2]-expr.range[1])/N
                         expr.range.max = expr.range[1]+k*(expr.range[2]-expr.range[1])/N
                         d0=dataset0[log2(dataset0[,ind])>=expr.range.min & log2(dataset0[,ind])<=expr.range.max, ]  # each bin has equal range of log2(expression values)

                         re0=lrm(expr=ex, data=d0, onlyR=TRUE)
                         HCP=sum(d0$CpG>0.4)/nrow(d0)
                         LCP=1-HCP
                         names(re0[[2]])=paste(cellline.tr(substr(ex,3,3)), names(re0[[2]]), sep=".")  # change H3k4me1 to NHEK.H3k4me1
                         r1=c(R=re0[[1]][1], Rcv=re0[[1]][2], re0[[2]], hCpG=HCP, lCpG=LCP)
                         bin.lrm=rbind(bin.lrm, r1)
            }
            bin.lrm=t(bin.lrm)
            colnames(bin.lrm)=paste(c(0:9)*10, "-", c(1:10)*10, "%", sep="")

            png(paste("relimp", ex,"png",sep="."), width=750, height=750)

            layout(matrix(seq(2), nrow=2, ncol=1, byrow=T), widths=5, heights=c(4,1), TRUE)

            # Expand right side of clipping rect to make room for the legend
            par(xpd=T, mar=c(1,4,4,8))
            bin.lrm.relimp = bin.lrm[grep("\\.", rownames(bin.lrm)),]
            cols=topo.colors(nrow(bin.lrm.relimp))
            pm = barplot(bin.lrm.relimp, width=1, col=cols, space=0.1, cex.axis=0.8, cex.names=0.6, las=1, cex.main=0.8,
                    main=paste("Relative importance of predictors for ",sub("\\n","\\.", label.tr(ex)),sep=""),
                    xlab="", xaxt="n",
                    ylab="percentages of relative importance")
                    #legend.text=T,
                    #args.legend = list(x=11, y=1, cex=0.8, fill=sort(rainbow(ncol(r11)), decreasing=F)))
            lines(pm, bin.lrm['R',], type='b', pch='*', col='red')
            lines(pm, bin.lrm['Rcv',], type='b', pch=1, col='blue')
            # Place the legend at (6,30) using heat colors
            legend(11,1, c(rev(rownames(bin.lrm.relimp)), "Spearman's rho", "Rho for cross-validation", "HCP", "LCP"),
                   cex=0.7,
                   pt.cex=c(rep(2,nrow(bin.lrm.relimp)),0.6,0.6,2,2),
                   col=c(rev(cols),'red','blue','black','gray'),
                   lty=c(rep(-1,nrow(bin.lrm.relimp)),1,1,-1,-1),
                   pch=c(rep(15,nrow(bin.lrm.relimp)),8,1,15,15));

            # add boxplot for CpG (either % of HCP/LCP or normalizedCpG value)
            par(xpd=T, mar=c(5,4,1,8))
            bin.lrm.cpg=bin.lrm[grep("CpG", rownames(bin.lrm)),]
            barplot(bin.lrm.cpg,  width=1, col=c("black","gray"), space=0.1, cex.axis=0.8, cex.names=0.8, las=1, cex.main=0.8,
                    xlab="Bins with equal range of log2(gene expression values) (in increasing order)",
                    ylab="CpG(%)")
            dev.off()
}

lrm.AtoB <- function(expA="cpkc2", expB=c("cpgc2", "cp1c2"), data=data)
{
             histones = c("Control", "Ctcf", "H3k27ac", "H3k27me3", "H3k36me3", "H3k4me1", "H3k4me2", "H3k4me3", "H3k9ac", "H4k20me1")
             cellline = cellline.tr(substr(expA,3,3))
             dat1=data[,c(expA, paste(cellline, ".", histones, sep=""))]  # CAGE ployA+ K562, vs. K562 hist

             dat1=dat1[apply(dat1>0, 1, sum)==ncol(dat1),]  # TODO: offset with max correlation
             dat1=log2(dat1)  # to be confirmed?

             ## Create a formula for a model with a large number of variables:
             xnam <- colnames(dat1)[grep("\\.", colnames(dat1))]
             # without intercross
             fmla <- as.formula(paste(expA, "~", paste(xnam, collapse= "+")))

             # # linear regression model
             fit1 <- lm(fmla, data=dat1)  # full model
             mear = dat1[,grep(expA,colnames(dat1))]
             pred = fitted(fit1)

             #r = sprintf("%.2f", cor(pred, mear))
             r = cor(pred, mear, method="spearman")

             for(exB in expB){
                          celllineB = cellline.tr(substr(exB,3,3))
                          dat2=data[, c(exB, paste(celllineB, ".", histones, sep=""))]  # CAGE ployA+ K562, vs. K562 hist
                          dat2=dat2[apply(dat2>0, 1, sum)==ncol(dat2),]  # TODO: offset with max correlation
                          dat2=log2(dat2)  # to be confirmed?
                          newdata = dat2[, grep("\\.", colnames(dat2))]
                          colnames(newdata)=sub(celllineB,cellline,colnames(newdata))
                          #r=c(r, sprintf("%.3f", cor(dat2[, grep(exB, colnames(dat2))], predict(fit1, newdata))))
                          r=c(r, cor(dat2[, grep(exB, colnames(dat2))], predict(fit1, newdata), method='spearman'))
             }
             return(r)
}


corplot <- function(dat, grep_expre="", grep_histmod="", filename="", filetype="")
{
             #debug
             #dat=data6
             #grep_expre=""
             #grep_histmod="H3k4"
             #pdf_filename=""
             #png_filename=""
             # filter all-zero lines
             rnaseq_hist =  dat[apply(dat>0, 1, sum)==ncol(dat),]
             # ---------- log2 transformation for RNA-seq data
             #rnaseqlog2_hist = cbind(log2(rnaseq_hist[,grep("\\.", colnames(rnaseq_hist), invert=T)]), rnaseq_hist[,grep("\\.", colnames(rnaseq_hist))])
             rnaseqlog2_hist = log2(rnaseq_hist)

             if(filename==""){
                          mytitle = "Smoothscatter of correlation analysis"
             }
             else{
                          mytitle = paste("Smoothscatter of correlation analysis", filename, sep=":")
             }

             if(filetype=="pdf"){
                          pdf(file=paste(filename,filetype,sep="."), title=mytitle, height =8, width = 8)
             }
             if(filetype=="png"){
                          png(filename=paste(filename,filetype,sep="."), height=1400, width = 1400)
             }

             panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
             {
                 usr <- par("usr"); on.exit(par(usr))
                 par(usr = c(0, 1, 0, 1))
                 r <- abs(cor(x, y))
                 txt <- format(c(r, 0.123456789), digits=digits)[1]
                 txt <- paste(prefix, txt, sep="")
                 if(missing(cex.cor)) cex <- 0.8/strwidth(txt)

                 test <- cor.test(x,y)
                 # borrowed from printCoefmat
                 Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
                               cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                               symbols = c("***", "**", "*", ".", " "))

                 text(0.5, 0.5, txt, cex = cex * r)
                 text(.8, .8, Signif, cex=cex, col=2)
             }
             ind_expr = grep("\\.", colnames(rnaseqlog2_hist), invert=T)
             ind_hist = grep("\\.", colnames(rnaseqlog2_hist))
             if(grep_expre!=""){
                          ind_expr=intersect(grep(grep_expre, colnames(rnaseqlog2_hist)),ind_expr)
             }
             if(grep_histmod!=""){
                          ind_hist=intersect(grep(grep_histmod, colnames(rnaseqlog2_hist)),ind_hist)
             }

             pairs(#rnaseqlog2_hist,
                   rnaseqlog2_hist[,c(ind_expr, ind_hist)],
                   labels=label.tr(colnames(rnaseqlog2_hist[,c(ind_expr, ind_hist)])),
                   lower.panel=function(...) {par(new=TRUE);smoothScatter(...)},
                   #lower.panel= panel.smooth,
                   upper.panel=panel.cor,
                   main = mytitle,
             )

             if(filetype!=""){
                          dev.off()
             }
}


tech.tr <- function(x)
{
             y=c();
             for(l in x){
                       y=c(y,switch(l,
                                       c = "CAGE",
                                       d = "DiTag",
                                       r = "RNAseq",
                                       l
                                       )
                          )
                       }
             return(y);
}

cellline.tr <- function(x)
{
             y=c();
             for(l in x){
                       y=c(y,switch(l,
                                       g = "Gm12878",
                                       '1' = "H1hesc",
                                       u = "HUVEC",
                                       e = "HelaS3",
                                       h = "HepG2",
                                       k = "K562",
                                       m = "MCF7",
                                       n = "NHEK",
                                       p = "Pros",
                                       '2'= "NHEK26",
                                       s = "HSMM",
                                       l
                                       )
                          )
                       }
             return(y);
}


# trick to get it in batch
# head -n1 tab.file | cut -f10-122 | sed 's/\t/\n/g' > experiments.id
# sed experiments.id 's/^r/RNAseq./g;s/^d/DiTag./g;s/^c/CAGE./g;s/\.p/.PolyA+\\n/g;s/\.t/.totalRNA\\n/g;s/\.m/.PolyA-\\n/g;s/nk/nK562./g;s/ng/nGm12878./g;s/\\nn/\\nNHEK./g;s/nu/nHUVEC./g;s/nh/nHepG2./g;s/ne/nHelaS3./g;s/\\n1/\\nH1hESC./g;s/np/nPros./g;s/nm/nMCF7./g;s/\.c/.Cy/g;s/\.n/.Nuc/g;s/\.w/.Wcell/g;s/\.h/.Chromatin/g;s/\.l/.Nucleoplasm/g;s/\.u/.Nucleolus/g;s/\.p/.Polysome/;' > experiments.full
# paste experiments.id experiments.full | sed 's/\t/ = "/g;s/$/",/g'

label.tr <- function(x)
{
             y=c();
             for(l in x){
                       y=c(y,switch(l,
cmec1 = "CAGE.PolyA-\nHelaS3.Cy1",
cmgc1 = "CAGE.PolyA-\nGm12878.Cy1",
cmgn1 = "CAGE.PolyA-\nGm12878.Nuc1",
cmhc1 = "CAGE.PolyA-\nHepG2.Cy1",
cmhn1 = "CAGE.PolyA-\nHepG2.Nuc1",
cmkp1 = "CAGE.PolyA-\nK562.Polysome1",
cmnc1 = "CAGE.PolyA-\nNHEK.Cy1",
cmnn1 = "CAGE.PolyA-\nNHEK.Nuc1",
cmuc1 = "CAGE.PolyA-\nHUVEC.Cy1",
cp1c2 = "CAGE.PolyA+\nH1hESC.Cy2",
cp1n2 = "CAGE.PolyA+\nH1hESC.Nuc2",
cp1w1 = "CAGE.PolyA+\nH1hESC.Wcell1",
cp1w2 = "CAGE.PolyA+\nH1hESC.Wcell2",
cpec1 = "CAGE.PolyA+\nHelaS3.Cy1",
cpec2 = "CAGE.PolyA+\nHelaS3.Cy2",
cpen1 = "CAGE.PolyA+\nHelaS3.Nuc1",
cpen2 = "CAGE.PolyA+\nHelaS3.Nuc2",
cpew1 = "CAGE.PolyA+\nHelaS3.Wcell1",
cpew2 = "CAGE.PolyA+\nHelaS3.Wcell2",
cpgc2 = "CAGE.PolyA+\nGm12878.Cy2",
cpgn1 = "CAGE.PolyA+\nGm12878.Nuc1",
cpgn2 = "CAGE.PolyA+\nGm12878.Nuc2",
cpgw1 = "CAGE.PolyA+\nGm12878.Wcell1",
cpgw2 = "CAGE.PolyA+\nGm12878.Wcell2",
cphc1 = "CAGE.PolyA+\nHepG2.Cy1",
cphc2 = "CAGE.PolyA+\nHepG2.Cy2",
cphn1 = "CAGE.PolyA+\nHepG2.Nuc1",
cphn2 = "CAGE.PolyA+\nHepG2.Nuc2",
cphw1 = "CAGE.PolyA+\nHepG2.Wcell1",
cphw2 = "CAGE.PolyA+\nHepG2.Wcell2",
cpkc1 = "CAGE.PolyA+\nK562.Cy1",
cpkc2 = "CAGE.PolyA+\nK562.Cy2",
cpkn1 = "CAGE.PolyA+\nK562.Nuc1",
cpkn2 = "CAGE.PolyA+\nK562.Nuc2",
cpkw1 = "CAGE.PolyA+\nK562.Wcell1",
cpkw2 = "CAGE.PolyA+\nK562.Wcell2",
cpmw1 = "CAGE.PolyA+\nMCF7.Wcell1",
cpnc3 = "CAGE.PolyA+\nNHEK.Cy3",
cpnn3 = "CAGE.PolyA+\nNHEK.Nuc3",
cpnw1 = "CAGE.PolyA+\nNHEK.Wcell1",
cpnw2 = "CAGE.PolyA+\nNHEK.Wcell2",
cpuc3 = "CAGE.PolyA+\nHUVEC.Cy3",
cpuc4 = "CAGE.PolyA+\nHUVEC.Cy4",
cpun3 = "CAGE.PolyA+\nHUVEC.Nuc3",
cpun4 = "CAGE.PolyA+\nHUVEC.Nuc4",
cpuw1 = "CAGE.PolyA+\nHUVEC.Wcell1",
cpuw2 = "CAGE.PolyA+\nHUVEC.Wcell2",
ct1w1 = "CAGE.totalRNA\nH1hESC.Wcell1",
cteu1 = "CAGE.totalRNA\nHelaS3.Nucleolus1",
ctgn1 = "CAGE.totalRNA\nGm12878.Nuc1",
ctgu1 = "CAGE.totalRNA\nGm12878.Nucleolus1",
cthn1 = "CAGE.totalRNA\nHepG2.Nuc1",
cthu1 = "CAGE.totalRNA\nHepG2.Nucleolus1",
ctkc1 = "CAGE.totalRNA\nK562.Cy1",
ctkh1 = "CAGE.totalRNA\nK562.Chromatin1",
ctkl1 = "CAGE.totalRNA\nK562.Nucleoplasm1",
ctkn1 = "CAGE.totalRNA\nK562.Nuc1",
ctku1 = "CAGE.totalRNA\nK562.Nucleolus1",
ctpw1 = "CAGE.totalRNA\nPros.Wcell1",
dp1w1 = "DiTag.PolyA+\nH1hESC.Wcell1",
dpec1 = "DiTag.PolyA+\nHelaS3.Cy1",
dpen1 = "DiTag.PolyA+\nHelaS3.Nuc1",
dpgc1 = "DiTag.PolyA+\nGm12878.Cy1",
dpgc2 = "DiTag.PolyA+\nGm12878.Cy2",
dpgn1 = "DiTag.PolyA+\nGm12878.Nuc1",
dphc1 = "DiTag.PolyA+\nHepG2.Cy1",
dphn1 = "DiTag.PolyA+\nHepG2.Nuc1",
dpkc1 = "DiTag.PolyA+\nK562.Cy1",
dpkn1 = "DiTag.PolyA+\nK562.Nuc1",
dpkp1 = "DiTag.PolyA+\nK562.Polysome1",
dpnc1 = "DiTag.PolyA+\nNHEK.Cy1",
dpnn1 = "DiTag.PolyA+\nNHEK.Nuc1",
dppw1 = "DiTag.PolyA+\nPros.Wcell1",
dpuc1 = "DiTag.PolyA+\nHUVEC.Cy1",
dpun1 = "DiTag.PolyA+\nHUVEC.Nuc1",
dtkh1 = "DiTag.totalRNA\nK562.Chromatin1",
dtkl1 = "DiTag.totalRNA\nK562.Nucleoplasm1",
dtku1 = "DiTag.totalRNA\nK562.Nucleolus1",
rp1c2 = "RNAseq.PolyA+\nH1hESC.Cy2",
rp1n2 = "RNAseq.PolyA+\nH1hESC.Nuc2",
rp1w1 = "RNAseq.PolyA+\nH1hESC.Wcell1",
rp1w2 = "RNAseq.PolyA+\nH1hESC.Wcell2",
rpec1 = "RNAseq.PolyA+\nHelaS3.Cy1",
rpec2 = "RNAseq.PolyA+\nHelaS3.Cy2",
rpen1 = "RNAseq.PolyA+\nHelaS3.Nuc1",
rpen2 = "RNAseq.PolyA+\nHelaS3.Nuc2",
rpew1 = "RNAseq.PolyA+\nHelaS3.Wcell1",
rpew2 = "RNAseq.PolyA+\nHelaS3.Wcell2",
rpgc1 = "RNAseq.PolyA+\nGm12878.Cy1",
rpgc2 = "RNAseq.PolyA+\nGm12878.Cy2",
rpgn1 = "RNAseq.PolyA+\nGm12878.Nuc1",
rpgn2 = "RNAseq.PolyA+\nGm12878.Nuc2",
rpgw1 = "RNAseq.PolyA+\nGm12878.Wcell1",
rpgw2 = "RNAseq.PolyA+\nGm12878.Wcell2",
rphc1 = "RNAseq.PolyA+\nHepG2.Cy1",
rphc2 = "RNAseq.PolyA+\nHepG2.Cy2",
rphn1 = "RNAseq.PolyA+\nHepG2.Nuc1",
rphn2 = "RNAseq.PolyA+\nHepG2.Nuc2",
rphw1 = "RNAseq.PolyA+\nHepG2.Wcell1",
rphw2 = "RNAseq.PolyA+\nHepG2.Wcell2",
rpkc1 = "RNAseq.PolyA+\nK562.Cy1",
rpkc2 = "RNAseq.PolyA+\nK562.Cy2",
rpkn1 = "RNAseq.PolyA+\nK562.Nuc1",
rpkn2 = "RNAseq.PolyA+\nK562.Nuc2",
rpkw1 = "RNAseq.PolyA+\nK562.Wcell1",
rpkw2 = "RNAseq.PolyA+\nK562.Wcell2",
rpnw1 = "RNAseq.PolyA+\nNHEK.Wcell1",
rpnw2 = "RNAseq.PolyA+\nNHEK.Wcell2",
rpuc4 = "RNAseq.PolyA+\nHUVEC.Cy4",
rpun3 = "RNAseq.PolyA+\nHUVEC.Nuc3",
rpun4 = "RNAseq.PolyA+\nHUVEC.Nuc4",
rpuw1 = "RNAseq.PolyA+\nHUVEC.Wcell1",
rpuw2 = "RNAseq.PolyA+\nHUVEC.Wcell2",
rm1w1 = "RNAseq.PolyA-\nH1hESC.Wcell1",
rm1w2 = "RNAseq.PolyA-\nH1hESC.Wcell2",
rmgw2 = "RNAseq.PolyA-\nGm12878.Wcell2",
rmhw1 = "RNAseq.PolyA-\nHepG2.Wcell1",
rmhw2 = "RNAseq.PolyA-\nHepG2.Wcell2",
rmkw1 = "RNAseq.PolyA-\nK562.Wcell1",
rmkw2 = "RNAseq.PolyA-\nK562.Wcell2",
rmnw1 = "RNAseq.PolyA-\nNHEK.Wcell1",
rmnw2 = "RNAseq.PolyA-\nNHEK.Wcell2",
rp2w1 = "RNAseq.PolyA+\nNHEK26.Wcell1",  # original from Sarah,where they distinguish NHEK and NHEK26
rplw1 = "RNAseq.PolyA+\nNHLF.Wcell1",
rplw2 = "RNAseq.PolyA+\nNHLF.Wcell2",
rpmw1 = "RNAseq.PolyA+\nMCF7.Wcell1",
rpmw2 = "RNAseq.PolyA+\nMCF7.Wcell2",
rpnn4 = "RNAseq.PolyA+\nNHEK.Nuc4",
rpsw1 = "RNAseq.PolyA+\nHSMM.Wcell1",
rpsw2 = "RNAseq.PolyA+\nHSMM.Wcell2",
l
))
                       }
             return(y);
}
