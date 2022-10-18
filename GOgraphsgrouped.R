eggnog <- read.csv('eggnoggGOtermscoma', header = FALSE, sep = "\t")


#Instalation of Go.db
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GO.db")

#Only once 
library(GO.db)


##Function to get the terms of each category and then maybe sort by just the level 2 Goterm in each gene non repetible
TerminMF <- c()
TerminBP <- c()
TerminCC <- c()
for (j in 1:length(eggnog$V3)){
  GoID<- as.list(strsplit(eggnog$V3[j], ",")[[1]])
  for( i in 1:length(GoID)){
    if (is.na(Ontology(GoID[[i]])[[1]])){
      
    }
    else {if(Ontology(GoID[[i]])[[1]] == "MF"){
      TerminMF[length(TerminMF)+1]<- Term(GoID[[i]])[[1]]
    }
    if(Ontology(GoID[[i]])[[1]]== "CC"){
      TerminCC[length(TerminCC)+1]<- Term(GoID[[i]])[[1]]
    }
    if(Ontology(GoID[[i]])[[1]]== "BP"){
      TerminBP[length(TerminBP)+1]<- Term(GoID[[i]])[[1]]
    }
    }
}
}
GotermsBP <- data.frame(TerminBP, stringsAsFactors=FALSE)
GotermsCC <- data.frame(TerminCC, stringsAsFactors=FALSE)
GotermsMF <- data.frame(TerminMF, stringsAsFactors=FALSE)
save(GotermsBP,file="GotermsclasBP.Rda")
save(GotermsMF,file="GotermsclasMF.Rda")
save(GotermsCC,file="GotermsclasCC.Rda")

#Pie Graphs
total <- length(TerminBP)+length(TerminCC)+length(TerminMF)
BP <- length(TerminBP)/total
MF <- length(TerminMF)/total
CC <- length(TerminCC)/total
pie(c(BP,MF,CC),labels=c(paste(round(100*BP,2),"%"),paste(round(100*MF,1),"%"),paste(round(100*CC,1),"%")),explode=0.1, main="",col=c("blue", "lightblue", "cyan"))
legend("topright", c("Biological Process","Molecular Function","Cellular Component"), cex = 0.58,fill = c("blue", "lightblue", "cyan"),box.lty=0,inset = -0.010)
title("GO Ontology Eggnog")

TermBP<-25081
TermCC<-16284
TermMF<-30185
total <- TermBP+TermCC+TermMF
BP <- TermBP/total
MF <- TermMF/total
CC <- TermCC/total
pie(c(BP,MF,CC),labels=c(paste(round(100*BP,2),"%"),paste(round(100*MF,1),"%"),paste(round(100*CC,1),"%")),explode=0.1, main="",col=c("blue", "lightblue", "cyan"))
legend("topleft", c("Biological Process","Molecular Function","Cellular Component"), cex = 0.58,fill = c("blue", "lightblue", "cyan"),box.lty=0,inset = -0.010)
title("GO Ontology HMMER2GO")

#Histograms

#Dataframesofcounts
Topgenes <- 15
countsBP <- dplyr::count(GotermsBP, TerminBP, sort = TRUE)[1:Topgenes,]
countsMF <- dplyr::count(GotermsMF, TerminMF, sort = TRUE)[1:Topgenes,]
countsCC <- dplyr::count(GotermsCC, TerminCC, sort = TRUE)[1:Topgenes,]

colnames(countsBP) <- c("Goterms", "counts")
colnames(countsMF) <- c("Goterms", "counts")
colnames(countsCC) <- c("Goterms", "counts")

countsTot <-  rbind(countsBP,countsMF,countsCC)

#plotting the top 15 Goterms of each Ontology

library(dplyr)
library(forcats)
library(ggplot2)

countsTot$Goterms <- factor(countsTot$Goterms , levels = countsTot$Goterms )
Gene_Ontology <- c(rep("Biological Process",Topgenes),rep("Molecular Function",Topgenes),rep("Cellular Component",Topgenes))
  
countsTot %>%#reorder( x,y)
  ggplot( aes(x=Goterms, y=counts, color=Gene_Ontology)) +
  geom_bar(stat="identity", fill='white', alpha=.6, width=.4) +
 # coord_flip() +
  #scale_x_discrete(guide = guide_axis(angle = 90)) +
  xlab("Go Terms") +
  ggtitle('Top 15 Go terms per ontology') +
  coord_flip() +
  ylab("Number of genes")+
  theme(plot.title = element_text(family = "Helvetica", face = "bold", size = (15), hjust = 0.5))




#Function to save all in a data frame
na.pad <- function(x,len){
  x[1:len]
}

makePaddedDataFrame <- function(l,...){
  maxlen <- max(sapply(l,length))
  data.frame(lapply(l,na.pad,len=maxlen),...)
}
makePaddedDataFrame(list(x=x,y=y,z=z))


#refreneces for more things
  GOID(xx[[1]])
  Term(xx[[1]])
  Synonym(xx[[1]])
  Secondary(xx[[1]])
  Definition(xx[[1]])
  Ontology(xx[[1]])
  
  
  
  #####Graphing each Ontology and comparing it with others:
TotalgenesPfam<-38952# this is equal to a hundred percent or 48 pixel for 50 percent
MFdammitnames<-c('binding','catalytic activity','transporter activity','molecular function activity','structural molecule activity','transcription regulation activity','molecular transducer activty','antioxidant activity','molecular carrier activity','translation regulator activity','nutrient reservoir activity','toxin activity','cargo receptor activity','synapse','nucleoid')
MFdammitvals<-c(57,125,263,280,281,286,291,296,298,300,300,300,300,301,301)
BPdammitnames<-c('cellular process','metabolic process','biological regulation process','regulation of biological process','localization','response to stimulus','biogenesis','signaling','multi organism process','multicellular organismal process','developmental process','biological adhesion','locomotion','negative regulation of biological process','positive regulaton of biological process')
BPdammitvals<-c(85,94,223,227,234,251,264,273,276,287,292,294,295,295,296)
CCdammitnames<-c('membrane','cell part', 'cell','organelle','membrane part','protein containing complex','organelle part','extracellular region','membrane enclosed lumen','virion','virion part','membrane complex','cell junction','other organism','other organism part')
CCdammitvals<-c(183,192,196,219,226,239,250,277,283,289,290,297,300,300,300)
totgenesHermodice<-32500
#MF 182-535 == 25000
MFHermodicenames<-c('molecular_function','binding','protein binding','hydrolase activity','nuclear binding','nucleic acid binding','transferase activity','transporter activity','DNA binding','molecular transducer activity','signal transducer activity','transfering phosphorus','kinase activity','peptidase activity','enzyme regulation activity')
MFHermodicevalues<-c(564,485,361,259,248,243,239,215,210,205,205,202,202,198,196)
#BP 181-571 20000
BPHermodicenames<-c('biological_process','metabolic process','cellular process','primary metabolic process','cellular metabolic process','biological regulation','regulation of biological process','macromolecule metabollic process','cellular component molecule metabolic process','nucleic acid metabolic process','cellular nitrogen component process','biosynthetic process','localization','transport','protein metabolic process')
BPHermodicevalues<-c(614,459,458,416,362,252,351,332,307,299,298,285,282,281,186)
# CC 181-531 20000
CCHermodicenames<-c('cellular_component','cell','cell part','intracellular','intracellular part','organelle','intracellular organelle','cytoplasm','intracellular membrane-bounded organelle','membrane-bounded organelle','cytoplasm part','nucleus','macromolecular complex','protein complex','intracellular organelle part')
CCHermodicevalues<-c(540,532,470,446,407,381,370,347,346,345,292,282,268,258,250)
totgeneseggnog<-43862

# AMking a consensus of the axis of the histogram
MF <- unique(append(countsMF$Goterms, MFdammitnames))
MF <- unique(append(MF, MFHermodicenames))
BP <- unique(append(countsBP$Goterms, BPdammitnames))
BP <- unique(append(BP, BPHermodicenames))
CC <- unique(append(countsCC$Goterms, CCdammitnames))
CC <- unique(append(CC, CCHermodicenames))

#########CC
#Making the same data frame for all three objects
##Starting with dammit
scaledammit<-2*(302 - 48)
CCdf <- data.frame(CCdammitnames,(302-CCdammitvals)*(100/scaledammit))
colnames(CCdf)<-c('Goterms','% genes dammit') 
diffCC<-setdiff(CC, CCdammitnames)
CCdammit<-data.frame(diffCC,rep(0,length(diffCC)))
colnames(CCdammit)<-c('Goterms','% genes dammit') 
CCdammit<-rbind(CCdf,CCdammit)

##Next is Hermodice
scaleHermoCC<-(531 - 181)#61.15384%
CCdfH <- data.frame(CCHermodicenames,(CCHermodicevalues-181)*(61.15384/scaleHermoCC))
colnames(CCdfH)<-c('Goterms','% genes Hermodice curunculata') 
diffCCHermo<-setdiff(CC, CCHermodicenames)
CCHermo<-data.frame(diffCCHermo,rep(0,length(diffCCHermo)))
colnames(CCHermo)<-c('Goterms','% genes Hermodice curunculata') 
CCHermo<-rbind(CCdfH,CCHermo)

#Next is Eggnog
CCdfEggnog <- data.frame(countsCC$Goterms,100*(countsCC$counts / totgeneseggnog))
colnames(CCdfEggnog)<-c('Goterms','% genes Eggnog') 
diffCCEggnog<-setdiff(CC, countsCC$Goterms)
CCEggnog<-data.frame(diffCCEggnog,rep(0,length(diffCCEggnog)))
colnames(CCEggnog)<-c('Goterms','% genes Eggnog') 
CCdfEggnog<-rbind(CCdfEggnog,CCEggnog)

#Merging all 3 data frames
CCHermo$type <- '% genes Hermodice curunculata'
CCdammit$type <- '% genes dammit'
CCdfEggnog$type <-'% genes Eggnog'

colnames(CCHermo)<-c('Goterms','% genes','type') 
colnames(CCdammit)<-c('Goterms','% genes','type') 
colnames(CCdfEggnog)<-c('Goterms','% genes','type') 

CCotherfinal<-rbind(CCHermo,CCdammit,CCdfEggnog)

library(dplyr)
library(forcats)
library(ggplot2)

e<-  ggplot(data=CCotherfinal, aes(x=reorder(Goterms, -CCotherfinal$`% genes`), y=CCotherfinal$`% genes`,fill=type)) +
  geom_bar( stat="identity",inherit.aes=TRUE, width=.4) +#,alpha = 0.5
  labs(color="Source of GO Terms")+
  #scale_x_discrete(guide = guide_axis(angle = 90)) +
  xlab("Go Terms") +
  ggtitle('Top 15 Go terms Cellular component') +
  coord_flip() +
  ylab("% of genes")+
  theme(plot.title = element_text(family = "Helvetica", face = "bold", size = (17), hjust = 0.5))+
  theme(axis.text.y = element_text(size = 14, hjust = 1, family = " Helvetica", face = "bold",colour="black"),#, face = "bold"
        plot.margin = margin(rep(15, 4)))+
  theme(legend.position = c(0.7, 0.8))+
  theme(axis.title.x = element_text(size = 15, family = " Helvetica", face = "bold",colour="black"),
        axis.title.y = element_text(size = 15, family = " Helvetica", face = "bold",colour="black"))+
  theme(legend.key.size = unit(0.8, "cm"))+
  theme(legend.text=element_text(color="black",size=16))+
  theme(axis.text.x = element_text(size = 14, hjust = 1, family = " Helvetica", face = "bold",colour="black"),#, face = "bold"
        plot.margin = margin(rep(15, 4)))
e +           # Remove background elements manually
  theme(legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA))
#########MF
  ##Starting with dammit
  scaledammit<-2*(302 - 48)
  MMdf <- data.frame(MFdammitnames,(302-MFdammitvals)*(100/scaledammit))
  colnames(MMdf)<-c('Goterms','% genes dammit') 
  diffMF<-setdiff(MF, MFdammitnames)
  MFdammit<-data.frame(diffMF,rep(0,length(diffMF)))
  colnames(MFdammit)<-c('Goterms','% genes dammit') 
  MFdammit<-rbind(MMdf,MFdammit)
  
  ##Next is Hermodice
  scaleHermoMF<-(535 - 182)#61.15384% 25000
  MFdfH <- data.frame(MFHermodicenames,(MFHermodicevalues-181)*(76.923076/scaleHermoMF))
  colnames(MFdfH)<-c('Goterms','% genes Hermodice curunculata') 
  diffMFHermo<-setdiff(MF, MFHermodicenames)
  MFHermo<-data.frame(diffMFHermo,rep(0,length(diffMFHermo)))
  colnames(MFHermo)<-c('Goterms','% genes Hermodice curunculata') 
  MFHermo<-rbind(MFdfH,MFHermo)
  
  #Next is Eggnog
  MFdfEggnog <- data.frame(countsMF$Goterms,100*(countsMF$counts / totgeneseggnog))
  colnames(MFdfEggnog)<-c('Goterms','% genes Eggnog') 
  diffMFEggnog<-setdiff(MF, countsMF$Goterms)
  MFEggnog<-data.frame(diffMFEggnog,rep(0,length(diffMFEggnog)))
  colnames(MFEggnog)<-c('Goterms','% genes Eggnog') 
  MFdfEggnog<-rbind(MFdfEggnog,MFEggnog)
  
  #Merging all 3 data frames
  MFHermo$type <- '% genes Hermodice curunculata'
  MFdammit$type <- '% genes dammit'
  MFdfEggnog$type <-'% genes Eggnog'
  
  colnames(MFHermo)<-c('Goterms','% genes','type') 
  colnames(MFdammit)<-c('Goterms','% genes','type') 
  colnames(MFdfEggnog)<-c('Goterms','% genes','type') 
  
  MFotherfinal<-rbind(MFHermo,MFdammit,MFdfEggnog)
  
  library(dplyr)
  library(forcats)
  library(ggplot2)
  
d<-  ggplot(data=MFotherfinal, aes(x=reorder(Goterms, -MFotherfinal$`% genes`), y=MFotherfinal$`% genes`,fill=type)) +
    geom_bar( stat="identity",inherit.aes=TRUE, width=.4) +
    labs(color="Source of GO Terms")+
    #scale_x_discrete(guide = guide_axis(angle = 90)) +
    xlab("Go Terms") +
    ggtitle('Top 15 Go terms Molecular Function') +
  coord_flip() +
  ylab("% of genes")+
  theme(plot.title = element_text(family = "Helvetica", face = "bold", size = (17), hjust = 0.5))+
  theme(axis.text.y = element_text(size = 14, hjust = 1, family = " Helvetica", face = "bold",colour="black"),#, face = "bold"
        plot.margin = margin(rep(15, 4)))+
  theme(legend.position = c(0.7, 0.8))+
  theme(axis.title.x = element_text(size = 15, family = " Helvetica", face = "bold",colour="black"),
        axis.title.y = element_text(size = 15, family = " Helvetica", face = "bold",colour="black"))+
  theme(legend.key.size = unit(0.8, "cm"))+
  theme(legend.text=element_text(color="black",size=16))+
  theme(axis.text.x = element_text(size = 14, hjust = 1, family = " Helvetica", face = "bold",colour="black"),#, face = "bold"
        plot.margin = margin(rep(15, 4)))


d +           # Remove background elements manually
  theme(legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA))

#########BP
  scaledammit<-2*(302 - 48)
  BPdf <- data.frame(BPdammitnames,(302-BPdammitvals)*(100/scaledammit))
  colnames(BPdf)<-c('Goterms','% genes dammit') 
  diffBP<-setdiff(BP, BPdammitnames)
  BPdammit<-data.frame(diffBP,rep(0,length(diffBP)))
  colnames(BPdammit)<-c('Goterms','% genes dammit') 
  BPdammit<-rbind(BPdf,BPdammit)
  
  ##Next is Hermodice
  scaleHermoBP<-(571 - 181)#61.15384%
  BPdfH <- data.frame(BPHermodicenames,(BPHermodicevalues-181)*(61.15384/scaleHermoBP))
  colnames(BPdfH)<-c('Goterms','% genes Hermodice curunculata') 
  diffBPHermo<-setdiff(BP, BPHermodicenames)
  BPHermo<-data.frame(diffBPHermo,rep(0,length(diffBPHermo)))
  colnames(BPHermo)<-c('Goterms','% genes Hermodice curunculata') 
  BPHermo<-rbind(BPdfH,BPHermo)
  
  #Next is Eggnog
  BPdfEggnog <- data.frame(countsBP$Goterms,100*(countsBP$counts / totgeneseggnog))
  colnames(BPdfEggnog)<-c('Goterms','% genes Eggnog') 
  diffBPEggnog<-setdiff(BP, countsBP$Goterms)
  BPEggnog<-data.frame(diffBPEggnog,rep(0,length(diffBPEggnog)))
  colnames(BPEggnog)<-c('Goterms','% genes Eggnog') 
  BPdfEggnog<-rbind(BPdfEggnog,BPEggnog)
  
  #Merging all 3 data frames
  BPHermo$type <- '% genes Hermodice curunculata'
  BPdammit$type <- '% genes dammit'
  BPdfEggnog$type <-'% genes Eggnog'
  
  colnames(BPHermo)<-c('Goterms','% genes','type') 
  colnames(BPdammit)<-c('Goterms','% genes','type') 
  colnames(BPdfEggnog)<-c('Goterms','% genes','type') 
  
  BPotherfinal<-rbind(BPHermo,BPdammit,BPdfEggnog)
  
  library(dplyr)
  library(forcats)
  library(ggplot2)
  
  p <-  ggplot(data=BPotherfinal, aes(x=reorder(Goterms, -BPotherfinal$`% genes`), y=BPotherfinal$`% genes`,fill=type)) +
    geom_bar( stat="identity",inherit.aes=TRUE, width=.4) +
    labs(color="Source of GO Terms")+
    xlab("Go Terms") +
    ggtitle('Top 15 Go terms Biological process') +
    coord_flip() +
    ylab("% of genes")+
    theme(plot.title = element_text(family = "Helvetica", face = "bold", size = (17), hjust = 0.5))+
    theme(axis.text.y = element_text(size = 14, hjust = 1, family = " Helvetica", face = "bold",colour="black"),#, face = "bold"
        plot.margin = margin(rep(15, 4)))+
    theme(legend.position = c(0.7, 0.8))+
    theme(axis.title.x = element_text(size = 15, family = " Helvetica", face = "bold",colour="black"),
          axis.title.y = element_text(size = 15, family = " Helvetica", face = "bold",colour="black"))+
    theme(legend.key.size = unit(0.8, "cm"))+
    theme(legend.text=element_text(color="black",size=16))+
    theme(axis.text.x = element_text(size = 14, hjust = 1, family = " Helvetica", face = "bold",colour="black"),#, face = "bold"
          plot.margin = margin(rep(15, 4)))


p +           # Remove background elements manually
    theme(legend.background = element_rect(fill = "transparent"),
          legend.box.background = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA))


###############################################
#Pygospio moreii
zero<-347
oneh<-38
total<-zero-oneh

##Biological process
# Adjusting the column names



GoTermBPPygospio<-c("metabolic process", "developmental process", "cell differentiation", "cell morphogenesis", "protein metabolic process", "cell organization and biogenesis","biosynthetic process","nucleic acid metabolic process","transport","catabolic process","primary metabolic process","organic substance metabolic process","cellular metabolic process","cellular process")
ValueBPPygospio<-c(197,317,330,334,297,304, 305, 305, 312, 314, 315,314,314,279)

ValueBPPygospio<- 100*(total-(ValueBPPygospio-oneh))/total

BPPygo<-data.frame(`Goterms`=GoTermBPPygospio,`% genes`=ValueBPPygospio)
BPPygo[nrow(BPPygo) + 1,] <- c("biological_process",65.7)


diffBPPygo<-setdiff(BP,BPPygo$Goterms)
valBP<-rep(0,length(diffBPPygo))
missingBP<-data.frame(`Goterms`=diffBPPygo,`% genes`=valBP)

BPPygo<-rbind(BPPygo,missingBP)




###Cellular component

GotermCCPygospio<-c("cell","intracellular","extracellular region","cytoplasm","nucleus","extracellular matrix", "chromosome","ribosome", "cytoskeleton","plasma membrane")

ValueCCPygospio<-c(72, 141, 347, 263, 298, 347, 341, 314, 323, 324)

ValueCCPygospio<-100*(total-(ValueCCPygospio-oneh))/total

CCPygo<-data.frame(`Goterms`=GotermCCPygospio,`% genes`=ValueCCPygospio)
CCPygo[nrow(CCPygo) + 1,] <- c("cellular_component", 9.54)

diffCCPygo<-setdiff(CC,CCPygo$Goterms)
valCC<-rep(0,length(diffCCPygo))
missingCC<-data.frame(`Goterms`=diffCCPygo,`% genes`=valCC)

CCPygo<-rbind(CCPygo,missingCC)




### Molecular function

GoTermMFPygospio<-c("catalytic activity","hydrolase activity","binding","transferase activity","transporter activity","nucleic acid binding","peptidase activity","DNA binding","protein binding","RNA binding")

ValueMFPygospio<- c(161, 260, 281, 305, 287, 329, 333, 347, 332, 333)


ValueMFPygospio<- 100*(total-(ValueMFPygospio- oneh))/total

MFPygo<-data.frame(`Goterms`=GoTermMFPygospio,`% genes`=ValueMFPygospio)

MFPygo[nrow(MFPygo) + 1,] <- c("molecular_function", 24.635)
diffMFPygo<-setdiff(MF,MFPygo$Goterms)## Maybe some of them are missing
valMF<-rep(0,length(diffMFPygo))
missingMF<-data.frame(`Goterms`=diffMFPygo,`% genes`=valMF)

MFPygo<-rbind(MFPygo,missingMF)


library(ggplot2)
Hermo<-expression(italic("Hermodice curunculata"))
Moreii<-expression(italic("Microspio moreii"))


# Graphing them all in one figure
BPHermo$type<-c(rep("Hermodice curunculata",length(BPHermo$type)))
BPdfEggnog$type<-c(rep("Microspio moreii",length(BPdfEggnog$type)))
BPPygo$type<-c(rep("Pygospio elegans",length(BPPygo$Goterms)))

 
colnames(BPPygo)<-c('Goterms','% genes','type')
#Selecting only the top 10
BPtermsgraph<-BPdfEggnog[1:10,]
BPHermographHermo<-BPHermo[BPHermo$Goterms %in% BPtermsgraph$Goterms,]
BPtermsgraph<-rbind(BPtermsgraph,BPHermographHermo)
BPPygospiograph<-BPPygo[BPPygo$Goterms %in% BPtermsgraph$Goterms,]


BPtermsgraph<-rbind(BPtermsgraph,BPPygospiograph)
BPfinal <-BPtermsgraph[order(BPtermsgraph$`% genes`,decreasing = TRUE),]

# This is to make an space in the xticks
BPfinal[nrow(BPfinal) + 1,] <- c("",0,"Microspio moreii")

## Now with MF

MFHermo$type<-c(rep("Hermodice curunculata",length(MFHermo$type)))
MFdfEggnog$type<-c(rep("Microspio moreii",length(MFdfEggnog$type)))
MFPygo$type<-c(rep("Pygospio elegans",length(MFPygo$Goterms)))

colnames(MFPygo)<-c('Goterms','% genes','type')

#Selecting only the top 10
MFtermsgraph<-MFdfEggnog[1:10,]
MFtermsgraphHermo<-MFHermo[MFHermo$Goterms %in% MFtermsgraph$Goterms,]
MFtermsgraph<-rbind(MFtermsgraph,MFtermsgraphHermo)
MFPygospiograph<-MFPygo[MFPygo$Goterms %in% MFtermsgraph$Goterms,]
MFtermsgraph<-rbind(MFtermsgraph,MFPygospiograph)
MFfinal <-MFtermsgraph[order(MFtermsgraph$`% genes`,decreasing = TRUE),]

### NOW for CC

CCHermo$type<-c(rep("Hermodice curunculata",length(CCHermo$type)))
CCdfEggnog$type<-c(rep("Microspio moreii",length(CCdfEggnog$type)))
CCPygo$type<-c(rep("Pygospio elegans",length(CCPygo$Goterms)))

colnames(CCPygo)<-c('Goterms','% genes','type')
#Selecting only the top 10
CCtermsgraph<-CCdfEggnog[1:10,]
CCtermsgraphHermo<-CCHermo[CCHermo$Goterms %in% CCtermsgraph$Goterms,]
CCtermsgraph<-rbind(CCtermsgraph,CCtermsgraphHermo)
CCPygospiograph<-CCPygo[CCPygo$Goterms %in% CCtermsgraph$Goterms,]
CCtermsgraph<-rbind(CCtermsgraph,CCPygospiograph)

CCfinal <-CCtermsgraph[order(CCtermsgraph$`% genes`,decreasing = TRUE),]

CCfinal[nrow(CCfinal) + 1,] <- c("a",0,"Microspio moreii")

allGOterms<-rbind(BPfinal,CCfinal,MFfinal)
colnames(allGOterms)<-c("Go_terms","Percentage","Organism")

allGOterms$Percentage <- as.double(allGOterms$Percentage)
row.names(allGOterms) <- NULL

group.colors <- c(`Hermodice curunculata` ="#6c6c6c" ,`Microspio moreii`="#a1a1a1",`Pygospio elegans`="#e1e1e1")

BPorder<-BPdfEggnog[1:10,]$Goterms
CCorder<-CCdfEggnog[1:10,]$Goterms
MForder<-MFdfEggnog[1:10,]$Goterms

level_order <- factor(allGOterms$Go_terms, levels=c(BPorder,"",CCorder,"a",MForder))
lev_ord<-c(BPorder,"",CCorder,"a",MForder)
labels_order<-c(BPorder,"",CCorder,"",MForder)

p <-  ggplot(data=allGOterms, aes(x=level_order, y=Percentage ,fill=Organism)) +
  geom_bar(position = position_dodge2(width = 0.9, preserve = "single"), stat="identity",inherit.aes=TRUE) +
  labs(fill = "")+
  xlab("Go Terms") +
  ylab("% of GO terms")+
   theme( legend.background = element_rect(fill = "transparent",linetype="blank"),
    legend.box.background = element_rect(fill = "transparent",linetype="blank"),
    panel.background = element_rect(fill = "transparent",linetype="blank"),
    panel.grid.major.x = element_blank() ,
    # explicitly set the horizontal lines (or they will disappear too)
    panel.grid.major.y = element_line( size=.05, color="black", ),
    axis.text.x = element_text(angle = 45, hjust=1,size= 10,color='black'),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.text = element_text(face = "italic"))+
    #scale_x_discrete(breaks=level_order[nchar(as.character(level_order))!=0],)+
    #scale_x_discrete(breaks=lev_ord,labels=labels_order)+
    scale_fill_manual(values=group.colors)+
    theme(legend.position = c(0.9, 0.82))+
    scale_y_continuous(breaks=seq(0,100,10))+
    annotate(geom="text",x="organic cyclic compound binding", label='atop(bold("Molecular function"))' , y = 88, parse = TRUE,size=6)+
  annotate(geom="text",x="metabolic process", label='atop(bold("Biological process"))' , y = 88, parse = TRUE,size=6)+
  annotate(geom="text",x="intracellular organelle", label= 'atop(bold("Cellular components"))' , y = 88, parse = TRUE,size=6)+
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14))
 
  ###Add the texts manueally

p
#missing labell space and title

