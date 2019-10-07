library(OmicCircos)
library(RJSONIO)
library(seqinr)
library(KEGGREST)
library(RColorBrewer)

setwd("U:/genome_example") #set your working directory

#######################################################################################
########################### create multi fasta file of amino acids
#######################################################################################

#read in the json annotation file from kbase:
json_file <- fromJSON("Stammera_capleta_CP024013_annotation.json")

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

#######################################################################################
########################### assemble all genome info
#######################################################################################

lvl1cat = c("Genetic Information Processing","Metabolism",
             "Environmental Information Processing") #select 3 categories from level 1 of the annotation
lvl2cat = c("Metabolism of cofactors and vitamins","Amino acid metabolism",
            "Carbohydrate metabolism") #select 3 categories from level 2 of the annotation

catcol = c(brewer.pal(7,"Accent"),"grey90") #define the colors to be used for the genes
names(catcol) = c(lvl1cat,"RNA",lvl2cat,"Other") #define which colors correspond to which category

json_file <- fromJSON("Stammera_capleta_CP024013_annotation.json") #import gene annotation

annot <- data.frame(id=NA,annot=NA,start=NA,end=NA,length=NA,mid=NA,
                    ko=NA,lvl1=NA,lvl2=NA,lvl3=NA) #create a table with all neccessary gene information
allfeat = c(json_file$features,json_file$non_coding_features[grep("RNA",unlist(lapply(json_file$non_coding_features,function(x){x$type})),ignore.case=T)])
for(i in 1:length(allfeat)){ #iterate through all features and fill the table with information
  annot[i,"gsize"] = json_file$dna_size
  annot[i,"id"] = allfeat[[i]]$id
  if(length(allfeat[[i]]$functions)!=0){
    annot[i,"annot"] = allfeat[[i]]$functions[1]
  }
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

gseq = read.fasta("Stammera_capleta_CP024013.fasta") #import the genome fasta file
slide = 1000 #define the sliding mirror for GC calculation
numelement <- ceiling(length(gseq[[1]])/slide)
pos1 = 0
pos2 = 0
start = 0
for(i in 1:numelement){
  pos1[i] = start
  pos2[i] = pos1[i] + slide
  start = start + slide
}
pos2[length(pos2)] = length(gseq[[1]])
gcs = vector()
for(i in 1:length(pos1)){
  gcs[i] = GC(gseq[[1]][pos1[i]:pos2[i]])
}
gcdat = data.frame(chr="chr1",seg.po=1:length(gcs),name1=gcs)

#######################################################################################
########################### make circos plot
#######################################################################################

par(mar=c(2, 2, 2, 2));# create a canvas
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");

#add the different circles of the genome plot
chrdat = cbind("chr1","270","630","0",annot$gsize[1]/1000,"0",annot$gsize[1]/1000)
colnames(chrdat) = c("seg.name","angle.start","angle.end","seg.sum.start","seg.sum.end","seg.start","seg.end")
genomedat <- data.frame("chr"="chr1","start"=annot$start/1000,
                      "end"=annot$end/1000,"value"=factor(annot$class),
                      "orientation"=annot$orientation,"color"=annot$color)
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

