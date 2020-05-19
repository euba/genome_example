#package installations
install.packages("RJSONIO")
install.packages("seqinr")

install.packages("RColorBrewer")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("OmicCircos")
BiocManager::install("KEGGREST")

#load the packages
library(RJSONIO)
library(seqinr)
library(RColorBrewer)
library(OmicCircos)
library(KEGGREST)

setwd("/Users/euba/Downloads/Dvul_symbiont_annotate.JSON/") #set your working directory

#######################################################################################
########################### create multi fasta file of amino acids
#######################################################################################

#read in the json annotation file from kbase:
json_file <- fromJSON("26.json")

annot <- list()
for(i in 1:length(json_file$features)){
  prot <- json_file$features[[i]]$protein_translation
  if(is.null(prot)){prot <- ""}
  annot[[json_file$features[[i]]$id]] = prot
}

#write a multi-fasta file with amino acid sequences for all genes
write.fasta(annot,names(annot), file="aasequences_genome.fasta")

#######################################################################################
########################### get categories for KO ids
#######################################################################################

kotab <- read.table("user_ko.txt",sep="\t",fill=T,row.names=NULL) #read in the KO ids assigned to each gene by BLASTKOALA
dupli <- which(duplicated(kotab[,1])) #remove duplicate ids for the same genes
if(length(dupli)!=0){kotab <- kotab[-dupli,]}
rownames(kotab) = kotab[,1]

for(i in 1:nrow(kotab)){ #look up the pathway categories for each KO id - this could take a while
  print(i)
  if(as.character(kotab[i,2])!=""){
    kresult <- keggGet(paste("ko",kotab[i,2],sep=":"))
    kotab[i,"lvl1"] <- kresult[[1]]$BRITE[2]
    kotab[i,"lvl2"] <- kresult[[1]]$BRITE[3]
    kotab[i,"lvl3"] <- kresult[[1]]$BRITE[4]
    kotab[i,"lvl4"] <- kresult[[1]]$BRITE[5]
  }
}
#this is just to format the categories and remove the first part of the string
kotab[,3] = unlist(lapply(strsplit(kotab[,3]," "),function(x){paste(x[3:length(x)],collapse=" ")}))
kotab[,4] = unlist(lapply(strsplit(kotab[,4]," "),function(x){paste(x[4:length(x)],collapse=" ")}))
kotab[,5] = unlist(lapply(strsplit(kotab[,5]," "),function(x){paste(x[5:length(x)],collapse=" ")}))
kotab[,6] = unlist(lapply(strsplit(kotab[,6]," "),function(x){paste(x[7:length(x)],collapse=" ")}))

write.table(kotab,file="kotab.tsv",row.names=TRUE,sep="\t")
kotab <- read.table("kotab.tsv",header=T,row.names=1,sep="\t")

#######################################################################################
########################### assemble all genome info - use this part if you have the annotation as json from kbase
#######################################################################################

lvl1cat = c("Genetic Information Processing","Metabolism",
             "Environmental Information Processing") #select 3 categories from level 1 of the annotation
lvl2cat = c("Metabolism of cofactors and vitamins","Amino acid metabolism",
            "Carbohydrate metabolism") #select 3 categories from level 2 of the annotation

catcol = c(brewer.pal(7,"Accent"),"grey90") #define the colors to be used for the genes
names(catcol) = c(lvl1cat,"RNA",lvl2cat,"Other") #define which colors correspond to which category

json_file <- fromJSON("26.json") #import gene annotation

annot <- data.frame(id=NA,chr=NA,annot=NA,start=NA,end=NA,length=NA,mid=NA,
                    ko=NA,lvl1=NA,lvl2=NA,lvl3=NA) #create a table with all neccessary gene information
allfeat = c(json_file$features,json_file$non_coding_features[grep("RNA",unlist(lapply(json_file$non_coding_features,function(x){x$type})),ignore.case=T)])
for(i in 1:length(allfeat)){ #iterate through all features and fill the table with information
  annot[i,"gsize"] = json_file$dna_size
  annot[i,"id"] = allfeat[[i]]$id
  if(length(allfeat[[i]]$functions)!=0){
    annot[i,"annot"] = allfeat[[i]]$functions[1]
  }
  annot[i,"chr"] = allfeat[[i]]$location[[1]][[1]]
  if(allfeat[[i]]$location[[1]][[3]]=="+"){
    annot[i,"start"] = allfeat[[i]]$location[[1]][[2]]
    annot[i,"end"] = allfeat[[i]]$location[[1]][[2]]+allfeat[[i]]$location[[1]][[4]]
  }else{
    annot[i,"start"] = allfeat[[i]]$location[[1]][[2]]-allfeat[[i]]$location[[1]][[4]]
    annot[i,"end"] = allfeat[[i]]$location[[1]][[2]]
  }
  annot[i,"orientation"] = allfeat[[i]]$location[[1]][[3]]
  annot[i,"length"] = allfeat[[i]]$location[[1]][[4]]
  annot[i,"mid"] = annot[i,"start"]+(annot[i,"length"]/2)
  kos <- kotab[annot[i,"id"],]
  if(!is.na(kos[,1])){
    annot[i,"ko"] = as.character(kos[,2])
    annot[i,"lvl1"] = as.character(kos$lvl1)
    annot[i,"lvl2"] = as.character(kos$lvl2)
    annot[i,"lvl3"] = as.character(kos$lvl3)
  }
  if(is.null(allfeat[[i]]$type)){allfeat[[i]]$type="other"}
  if(allfeat[[i]]$type %in% c("tRNA","RNA","rRNA","rna","rrna","trna")){
    annot[i,"lvl1"] = "RNA"
    annot[i,"lvl2"] = "RNA"
    annot[i,"lvl3"] = "RNA"
  }
}
#Assign color of right category for each gene
annot$class = NA
for(i in lvl2cat){annot$class[which(annot$lvl2==i)]=i}
for(i in lvl1cat){annot$class[intersect(which(annot$lvl1==i),which(is.na(annot$class)))]=i}
annot$class[intersect(which(annot$lvl1=="RNA"),which(is.na(annot$class)))] = "RNA"
annot$class[which(is.na(annot$class))] = "Other"
annot$color = catcol[as.character(annot$class)]

#######################################################################################
########################### get GC skew information
#######################################################################################

gseq = read.fasta("Dvul_gc_selected_contigs.fa") #import the genome fasta file
slide = 1000 #define the sliding window for GC calculation
gcdat = data.frame(chr=NA,seg.po=NA,name1=NA)
for(i in 1:length(gseq)){
  numelement <- ceiling(length(gseq[[i]])/slide)
  pos1 = 0
  pos2 = 0
  start = 0
  for(j in 1:numelement){
    pos1[j] = start
    pos2[j] = pos1[j] + slide
    start = start + slide
  }
  pos2[length(pos2)] = length(gseq[[i]])
  gcs = vector()
  for(j in 1:length(pos1)){
    gcs[j] = GC(gseq[[i]][pos1[j]:pos2[j]])
  }
  gccont = data.frame(chr=names(gseq)[i],seg.po=1:length(gcs),name1=gcs)
  gcdat = rbind(gcdat,gccont)
}
gcdat = gcdat[-1,]

#######################################################################################
########################### make circos plot
#######################################################################################

par(mar=c(2, 2, 2, 2));# create a canvas
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");

#add the different circles of the genome plot
clengths = unlist(lapply(gseq,length))
#set a list of very short contigs you want to remove --> use at your own risk!
tabu <- names(which(clengths<1000))
if(length(tabu) != 0){
  clengths = clengths[-which(names(clengths) %in% tabu)]
  annotsel = annot[-which(annot$chr %in% tabu),]
}else{annotsel = annot}

degree2 = 270+cumsum(360*(clengths/sum(clengths)))
degree1 = c(270,degree2)
degree1 = degree1[-length(degree1)]
#set a short separator degree between contigs --> uset at your own risk!
degree2 = degree2 - 0.5

chrdat = cbind(names(clengths),
               as.character(degree1),
               as.character(degree2),
               "0",clengths/slide,
               "0",clengths/slide)
colnames(chrdat) = c("seg.name","angle.start","angle.end","seg.sum.start","seg.sum.end","seg.start","seg.end")
genomedat <- data.frame("chr"=annotsel$chr,"start"=annotsel$start/1000,
                      "end"=annotsel$end/1000,"value"=factor(annotsel$class),
                      "orientation"=annotsel$orientation,"color"=annotsel$color)
circos(R=355,W=20,cir=chrdat,type="chr",col="black",scale=F,print.chr.lab=F,cex=4,lwd=10)
circos(R=320,cir=chrdat,mapping=genomedat[which(genomedat$orientation=="+"),],W=0,type="arc2",scale=T,print.chr.lab=F,cex=4,
       col=as.character(genomedat[which(genomedat$orientation=="+"),"color"]),lwd=50,cutoff=0)#8x8
circos(R=270,cir=chrdat,mapping=genomedat[which(genomedat$orientation=="-"),],W=0,type="arc2",scale=T,print.chr.lab=F,cex=4,
       col=as.character(genomedat[which(genomedat$orientation=="-"),"color"]),lwd=50,cutoff=0)#8x8
circos(R=175,cir=chrdat,W=65,mapping=gcdat,col.v=3,type="ls",B=F,col="grey",lwd=1.5,scale=F);
#save this plot as 8x8 PDF file

#plot the legend
par(mar=c(0,0,0,0))
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");
legend("top",legend=names(catcol),col=catcol,lty=1,cex=1,box.lwd=0,lwd=14)

