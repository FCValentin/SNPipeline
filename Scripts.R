############
##Packages##
############

source("http://www.bioconductor.org/biocLite.R")
biocLite("illuminaio")
biocLite("crlmm")
biocLite("oligoClasses")
biocLite("SNPchip")
biocLite("stringr")
install.packages("png")
library(png)
library(illuminaio)
library(stringr)
library(oligoClasses)
library(matrixStats)
library(crlmm)
library(ff)
library(SNPchip)
library(tools)
library(grid)


###################################################################################  
##path to idatfiles:Here your directory with manifest file,idat and HC documents.##  
###################################################################################


setwd("XXX")
datadir<-getwd()


###########################################################
##Name of the Samplesheet link (read .idat files further)##
###########################################################

Sample<-"XXX"

###############################################################################
##For each cell,patient (family) and typeCell (Parental/hiPSC) in tsv format ##
###############################################################################

Family<-"XXX.txt"


###############################################################################################################
##For annotation, put the Manifest link(csv) with ";" as field separator & release control section at the end##
###############################################################################################################

manifest<-"XXX.csv"
genome<-'hg19'


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
########################
##Pre-process Analysis##
########################

#CreateffFiles
outdir<-paste(datadir,"ff_files",sep="/")
outdir<-paste(outdir,Sample,sep="_")
#Create Directory so as to save ff files
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
ldPath(outdir)
#ReadSamples
name<-read.csv(file.path(datadir,Sample),header=FALSE,as.is=TRUE)
name<-name[9,]
samplesheet=read.csv(file.path(datadir,Sample),header=FALSE,col.names=name,as.is=TRUE)
samplesheet<-samplesheet[-c(1:9),-c(4:ncol(samplesheet))]
#We create 2 vectors for Genotype.Illumina so as to read .idat files
arrayNames<-file.path(datadir,paste((samplesheet[,"SentrixBarcode_A"]),(samplesheet[,"SentrixPosition_A"]),sep="_"))
arrayInfo<-list(barcode="SentrixBarcode_A",position="SentrixPosition_A")
#ReadManifest
manifest<-read.csv2(manifest,header=F,fill=T,stringsAsFactors=FALSE)
manifest<-manifest[-c(1:8),]
names(manifest)<-c("IlmnID","Name","IlmnStrand","SNP","AddressA_ID","AlleleA_ProbeSeq","AddressB_ID","AlleleB_ProbeSeq","GenomeBuild","chromosome","MapInfo","Ploidy","Species","Source","SourceVersion","SourceStrand","SourceSeq","TopGenomicSeq","BeadSetID","Exp_Clusters","RefStrand"
)
#Convert sexual chromosome into numeric
manifest$chromosome[manifest$chromosome=="X"]=23
manifest$chromosome[manifest$chromosome=="Y"]=24
manifest$chromosome[manifest$chromosome=="XY"]=25
manifest$chromosome[manifest$chromosome=="MT"]=26
manifest$chromosome<-as.integer(manifest$chromosome)
#Create 3 vectors for Genotype.Illumina isSNP,position and chromosome
manifest$isSnp<-TRUE
Usefull<-manifest[,c(2,10,11)]
names(Usefull)<-c("featureNames","chr","position")
Usefull$position<-as.integer(Usefull$position)
manifest<-cbind(manifest,Usefull)
manifest$position<-as.integer(manifest$position)
anno<-manifest
chr<-as.character(manifest$chromosome)

############
##Analysis##
############

cnReSet<-genotype.Illumina(sampleSheet=samplesheet,arrayNames=arrayNames,ids=NULL,path="",arrayInfoColNames=arrayInfo,
highDensity=FALSE,sep="_",fileExt=list(green="Grn.idat",red="Red.idat"),XY=NULL,anno=manifest,genome=genome,
call.method="krlmm",trueCalls=NULL,cdfName='nopackage',copynumber=TRUE,batch=NULL,saveDate=TRUE,stripNorm=TRUE,
useTarget=TRUE,quantile.method="between",nopackage.norm="loess",mixtureSampleSize=100,fitMixture=TRUE,
eps=0.1,verbose=TRUE,seed=10,probs=rep(1/3,3),DF=6,SNRMin=5,recallMin=2,recallRegMin=100,gender=NULL,returnParams=TRUE,badSNP=0.7)
open(cnReSet)
#A SNP with a conf score less than the Threshold will be exclude
#Change Threshold value, decrease Threshold so as to be less stringent
thr<-0.90
#Calculate BAF and CNV
crlmmCopynumber(cnReSet, MIN.SAMPLES=10, SNRMin = 5, MIN.OBS = 1,
                DF.PRIOR = 50, bias.adj = FALSE,
                prior.prob = rep(1/4, 4), seed = 1, verbose = TRUE,
                GT.CONF.THR = thr, MIN.NU = 8, MIN.PHI = 8,
                THR.NU.PHI = TRUE, type=c("SNP", "NP", "X.SNP", "X.NP"),
                fit.linearModel=TRUE)
label<-FALSE

#######################################
##Plot and CNV/BAF Analysis Function ##
#######################################

exp<-function(cnReSet,mychr,mysample){ marker.index <- which(chromosome(cnReSet) == mychr)
  #BAF
  A <- CA(cnReSet, i = marker.index, j = mysample)
  B <- CB(cnReSet, i = marker.index, j = mysample)
  BAF<-B/(A+B)
  #CNV
  tot<<-tot+1
  cn <- totalCopynumber(cnReSet, i = marker.index, j = mysample)
  #Position chromosome
  x <- position(cnReSet)[marker.index]
  #Title of png and his index
  mychro<-paste("/ chr.",mychr)
  title<-paste(paste(Sample,mychr,sep="_Chr_"),"png",sep=".")
  pngfile<-title
  pngF[tot]<<-pngfile
  png(title)
  #Divide output graphic
  par(mfcol=c(2,1),mex=0.5,cex=0.5,mar=c(4.1,0, 6.6, 8.0))
  #Chromosomes 13/14/15/21/22 too long, we adapt the plot
    if((mychr==13)|(mychr==14)|(mychr==15)|(mychr==21)|(mychr==22)){
      #Plot BAF
      mychro1<-paste(mychro,"BAF",sep="   :  ")
      plot(x, BAF,cex.main=4, pch = 16,las=1,yaxt="n",xlab='', cex = 1.4, xaxt = "n", col = "chartreuse4",ylim = c(0, 1),xlim=c(0,max(x)), ylab = "BAF", main = paste(samplesheet$Sample_ID[mysample], mychro1,sep=" /"))
      ytick<-seq(0, 1, by=0.5)
      axis(side=4, at=ytick, labels = T,las=1,cex.axis=3)
      #Plot CNV
      mychro2<-paste(mychro,"CNV",sep="   :  ")
      plot(x, cn,cex.main=4, pch = 16,las=1,yaxt="n",xlab='', cex = 1.4, xaxt = "n", col = "blue",ylim = c(-2, 4),xlim=c(0,max(x)), ylab = "CNV",main = paste(samplesheet$Sample_ID[mysample], mychro2,sep=" /"))
      ytick<-seq(0, 4, by=2)
      axis(side=4, at=ytick, labels = T,las=1,cex.axis=3)
      #Plot Chromosome
      plotIdiogram(mychr,"hg19",unit =c("bp","Mb"), new = F, label.cytoband = label,cytoband.ycoords=c(-2,-1),label.y=c(-3),verbose=TRUE)
      }else{
      #Plot BAF
      mychro1<-paste(mychro,"BAF",sep="   :  ")
      plot(x, BAF,cex.main=4, pch = 16,yaxt="n", cex = 1.4,xlab='', xaxt = "n", col = "chartreuse4",ylim = c(0, 1), ylab = "BAF",main = paste(samplesheet$Sample_ID[mysample], mychro1,sep=" /"))
      ytick<-seq(0, 1, by=0.5)
      axis(side=4, at=ytick, labels = T,las=1,cex.axis=3)
      #Plot CNV
      mychro2<-paste(mychro,"CNV",sep="   :  ")
      plot(x, cn,cex.main=4, pch = 16,yaxt="n",las=1, cex = 1.4,xlab='', xaxt = "n", col = "blue",ylim = c(-2, 4), ylab = "CNV",main = paste(samplesheet$Sample_ID[mysample], mychro2,sep=" /"))
      ytick<-seq(0, 4, by=2)
      axis(side=4, at=ytick, labels = T,las=1,cex.axis=3)
      #Plot Chromosome
      plotIdiogram(mychr,"hg19",unit =c("bp","Mb"), new = F, label.cytoband = label,cytoband.ycoords=c(-2,-1),label.y=c(-3),verbose=TRUE)
      }
  dev.off()
  #We release forbidden character for Windows(file's title)
  if(mychr==23){
  Export<-as.character(samplesheet$Sample_ID[mysample])
  Export<-str_replace(Export,"/","_")
  Export<-str_replace(Export,":","_")
  Export<-str_replace(Export,"<","_")
  Export<-str_replace(Export,">","_")
  Export<-str_replace(Export,"|","_")
  Export<-str_replace(Export,'"',"_")
  Result<-paste(datadir,"Resultats",sep="/")
  #Create Directory "Resultats" so as to save datas
  dir.create(Result,recursive=TRUE,showWarnings=FALSE)     
  setwd(Result)
  #Export BAF/CNV TSV format
  chr<-2:23
  marker<- which(chromosome(cnReSet) == 1)
  for(y in chr){
  marker.index <- which(chromosome(cnReSet)== y)
  marker<-c(marker,marker.index)
  }
  Ae <- CA(cnReSet, i = marker, j = mysample)
  Be <- CB(cnReSet, i = marker, j = mysample)
  BAFe<-Be/(Ae+Be)
  cne <- totalCopynumber(cnReSet, i = marker, j = mysample)
  table<-cbind(BAFe,cne)
  colnames(table)<-c("BAF","CNV")
  write.table(table,file=paste(Export,"tsv",sep="."),sep="\t",col.names=NA,dec=".")
  read.table(file=paste(Export,"tsv",sep="."),header=TRUE,row.names=1)
  setwd(datadir)
  }
  #Create a new vector with BAF,CNV and position for this chromosome and sample
  res<-cbind(BAF,cn,x)
  res<-data.frame(res)
  names(res)<-c("BAF","CNV","x")
  #Remove lines with NA
  res<-res[-which(is.na(res)),]
  res<-res[order(res[,"x"],decreasing=F), ]
  res$x<-res$x/1000000
  value<-NULL
  warmi<-0
  warma<-0
  memo<-0
  #We check if we have some BAF/CNV datas (Threshold limit)
  if(length(res$x)>2000){
    mini<-res[which(res$BAF<0.405&res$BAF>0.19),]
    mini<-mini[order(mini[,"x"],decreasing=F), ]
    maxi<-res[which(res$BAF<0.81&res$BAF>0.595),]
    maxi<-maxi[order(maxi[,"x"],decreasing=F), ]
    del<-res[which(res$BAF<0.8&res$BAF>0.20),]
    del<-del[order(del[,"x"],decreasing=F), ]
    biais<-res[which(res$BAF<=0.50&res$BAF>=0),]
    bar<-length(res$x)/700
    if(bar<7){
      bar<-7
    }
    if(((abs(length(biais$x)/(length(res$x))))<0.55)&(abs(length(biais$x)/(length(res$x))))>0.45){
      #We check if we have many points in a small interval (duplication)
      if(length(mini$x)>5&length(maxi$x)>5){
        pos<--500000000000
        posb<--500000000000
        for(tai in 1:length(mini$x)){
          posmin<-pos-0.08
          posmax<-pos+0.08
          pos<-mini$x[tai]
          if(pos<posmax&pos>posmin){
            warmi<-warmi+1
          }else{
            warmi<-0
          }
          if(warmi>bar&memo==0){
            #We save the Chromosome with duplication
            coucou<-paste(paste("Attention, il y a une diminution de la BAF sur le génome et chromosome : ",(samplesheet$Sample_ID[mysample]),mychro,sep="/"))
            warning(paste(coucou))
            memo<-1
            value<-mychr
          }
        }
        for(taib in 1:length(maxi$x)){
          posminb<-posb-0.08
          posmaxb<-posb+0.08
          posb<-maxi$x[taib]
          if(posb<posmaxb&posb>posminb){
            warma<-warma+1
          }else{
            warma<-0
          }
          if(warma>bar&memo==0){
            #We save the Chromosome with duplication
            coucou<-paste(paste("Attention, il y a une augmentation de la BAF sur le génome et chromosome : ",(samplesheet$Sample_ID[mysample]),mychro,sep=" / "))
            warning(paste(coucou))
            memo<-1
            value<-mychr
          }
        }
      }
      #Check if we have deletion of chromsome and save it
      if(((length(del$x))/length(res$x))<0.17){
        coucou<-paste(paste("Attention, il y a une délétion sur le génome et chromosome : ",(samplesheet$Sample_ID[mysample]),mychro,sep=" / "))
        warning(paste(coucou))
        value<-mychr
      }
    }
  }  
  #Return chromosome if he had a problem.
  if(!is.null(value)){
    return(value)
  }
}

##############################
##Create PDF for each family##
##############################

setwd(datadir)
#Select different family and number of family
Group<-read.table(Family,header=T,sep="\t")
images<-paste(datadir,"images",sep="/")
#Create Directory so as to save images
dir.create(images,recursive=TRUE,showWarnings=FALSE)
setwd(images)
compt<-levels(Group$Patient)
borne<-nlevels(Group$Patient)
  for(n in 1:borne){
    #Read each family of sample, save Fibroblaste (Parental) and iPSC (hiPSC)
    setwd(images)
    Famille<-Group[which(Group$Patient==compt[n]),]
    Fibro<-Famille[which(Group$CellType=="Parental"),]
    Fibro<-subset(Fibro,!is.na(Fibro$Sample_ID))
    iPSC<-Famille[which(Group$CellType=="hiPSC"),]
    iPSC<-subset(iPSC,!is.na(iPSC$Sample_ID))
    a<-as.numeric(row.names(iPSC))
    b<-as.numeric(row.names(Fibro))
    lines<-c(b,a)
    #Create PDF name
    remove(Echant,name,lieu)
    Echant<-as.character(unique(Famille$Patient))
    Echant<-str_replace(Echant,"/","_")
    Echant<-str_replace(Echant,":","_")
    Echant<-str_replace(Echant,"<","_")
    Echant<-str_replace(Echant,">","_")
    Echant<-str_replace(Echant,"|","_")
    Echant<-str_replace(Echant,'"',"_")
    name<-paste(Echant,"pdf",sep=".")
    lieu<-paste(images,Echant,sep="/")
    #Create Directory so as to save images for each family
    dir.create(lieu,recursive=TRUE,showWarnings=FALSE)
    setwd(lieu)
    #Plot
    #Create pdf for this family
    pdf(name,width=10,height=10)
      for(indice in lines){
        tot<<-0
        val<-NULL 
        pngF<<-NULL
        #For each sample, plot BAF and CNV for 23 chromosomes in bitmap format and save chromosome with problem
            for(j in 1:23){
              setwd(lieu)
              valu<-exp(cnReSet,j,indice)
              if(is.null(val)){
                val<-valu
              }else{
                val<-c(val,valu)
              }
            }
        #Zoom on chromosome with a problem, 1/page
      par(mfcol=c(4,6),mex=0.5,cex=0.5,mar=c(4.1,0, 3.1, 2.1))
        for(j in 1:23){
          #Read png and save 23 bitmap plot ( 1 page in pdf in bitmap format )
          setwd(lieu)
          pngFile <- pngF[j]
          pn <- readPNG(pngFile)
          #Erase png
          unlink(pngFile)
          setwd(images)
          plot.new()
          rasterImage(pn,0,0,1,1)
        }
      #Zoom (plot) for each chromosome with an anomaly ( 1 page/anomaly in vectorial format)
        for(monchro in val){
        par(mfrow=c(2,1),mar=c(5.1, 4.1, 4.1, 2.1))
          marker.index <- which(chromosome(cnReSet) == monchro)
          A <- CA(cnReSet, i = marker.index, j = indice)
          B <- CB(cnReSet, i = marker.index, j = indice)
          BAF<-B/(A+B)
          cn <- totalCopynumber(cnReSet, i = marker.index, j = indice)
          x <- position(cnReSet)[marker.index]
          mychro<-paste("/ chr.",monchro)
          #Chromosomes 13/14/15/21/22 too long, we adapt the plot
          if((monchro==13)|(monchro==14)|(monchro==15)|(monchro==21)|(monchro==22)){
            #Plot BAF
            mychro1<-paste(mychro,"BAF",sep="   :  ")
            plot(x, BAF,cex.main=0.9, pch = 16,las=1,yaxt="n",xlab='', cex = 1, xaxt = "n", col = "chartreuse4",ylim = c(0, 1),xlim=c(0,max(x)), ylab = "BAF", main = paste(samplesheet$Sample_ID[indice], mychro1,sep=" /"))
            ytick<-seq(0, 1, by=0.25)
            axis(1, at = pretty(x), labels = pretty(x)/1e+06)
            axis(side=2, at=ytick, labels = T,las=1,cex.axis=1)
            #Plot CNV
            mychro2<-paste(mychro,"CNV",sep="   :  ")
            plot(x, cn,cex.main=0.9, pch = 16,las=1,yaxt="n",xlab='', cex = 1, xaxt = "n", col = "blue",ylim = c(-2, 4),xlim=c(0,max(x)), ylab = "CNV",main = paste(samplesheet$Sample_ID[indice], mychro2,sep=" /"))
            ytick<-seq(0, 4, by=2)
            axis(1, at = pretty(x), labels = pretty(x)/1e+06)
            axis(side=2, at=ytick, labels = T,las=1,cex.axis=1)
            #Plot Chromosome
            plotIdiogram(monchro,"hg19",unit =c("bp","Mb"), new = F, label.cytoband = label,cytoband.ycoords=c(-2,-1),label.y=c(-3),verbose=TRUE)
          }else{
            #Plot BAF
            mychro1<-paste(mychro,"BAF",sep="   :  ")
            plot(x, BAF,cex.main=0.9, pch = 16,yaxt="n", cex = 1,xlab='', xaxt = "n", col = "chartreuse4",ylim = c(0, 1), ylab = "BAF",main = paste(samplesheet$Sample_ID[indice], mychro1,sep=" /"))
            ytick<-seq(0, 1, by=0.25)
            axis(side=2, at=ytick, labels = T,las=1,cex.axis=1)
            axis(1, at = pretty(x), labels = pretty(x)/1e+06)
            #Plot CNV
            mychro2<-paste(mychro,"CNV",sep="   :  ")
            plot(x, cn,cex.main=0.9, pch = 16,yaxt="n",las=1, cex = 1,xlab='', xaxt = "n", col = "blue",ylim = c(-2, 4), ylab = "CNV",main = paste(samplesheet$Sample_ID[indice], mychro2,sep=" /"))
            ytick<-seq(0, 4, by=2)
            axis(1, at = pretty(x), labels = pretty(x)/1e+06)
            axis(side=2, at=ytick, labels = T,las=1,cex.axis=1)
            #Plot Chromosome
            plotIdiogram(monchro,"hg19",unit =c("bp","Mb"), new = F, label.cytoband = label,cytoband.ycoords=c(-2,-1),label.y=c(-3),verbose=TRUE)
          }
        }
      }
    #Close pdf for this family of sample
    dev.off()
    setwd(datadir)
  }
