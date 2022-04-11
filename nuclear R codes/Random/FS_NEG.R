library(Matrix)
library(kmed)
library(pROC)
library(class)
library(caret)
library(kernlab)

setwd("f:/finalresults/nuclear/randomnuclear")
#---------------------------------------
mynormalize <- function(x) 
{ 
  for(j in 1:length(x[1,]))
  {
    # print(j)
    min <- min(x[,j])
    max <- max(x[,j])
    for(i in 1:length(x[,j])){
      x[i,j] <-  (x[i,j] - min)/( max - min) 
    }
  }
  return(x)
}

#-----------------------------------------------------------
dclustno<-6
pclustno<-4
negno<-180
dpclustdim<-dclustno*pclustno
x=read.delim("nucleardescriptors.txt",header = T)
y<-x[,1]
dnum<-length(y)
dfnum<-length(x[1,])
#-------------------------------------------
x=read.delim("nuclearfeatures.txt",header = T)
y<-x[,1]
pfnum<-length(x[1,])
pnum<-length(y)
#------
dclusterfile<-read.delim(file="dcluster.txt",sep="\t",header = T)
dcluster<-dclusterfile[,1]
names(dcluster)<-rownames(dclusterfile)

pclusterfile<-read.delim(file="pcluster.txt",sep="\t",header = T)
pcluster<-pclusterfile[,1]
names(pcluster)<-rownames(pclusterfile)

dclusterfile2<-read.delim(file="dcluster2.txt",sep="\t",header = T)
dcluster2<-dclusterfile2[,1]
names(dcluster2)<-rownames(dclusterfile2)

pclusterfile2<-read.delim(file="pcluster.txt",sep="\t",header = T)
pcluster2<-pclusterfile2[,1]
names(pcluster2)<-rownames(pclusterfile2)

#------------------------------------------------------

dfeatures<-read.delim("nucleardescriptors.txt",header = T)
pfeatures<-read.delim("nuclearfeatures.txt",header = T)
dfnum<-length(dfeatures[1,])
dnum<-length(dfeatures[,1])
pnum<-length(pfeatures[,1])
pfnum<-length(pfeatures[1,])
allcol<-dfnum+pfnum+1  # drugname + drug features+ target name + taeget features + label column
allcol
dpfeatures<-merge(dfeatures,pfeatures,by=NULL)
dpfeatures<-dpfeatures[,c(1,(dfnum+1),(2:dfnum),(dfnum+2):(allcol-1))]
#write.table(dpfeatures,file="dpfeatures.txt",sep="\t",row.names = F,quote = F)
dpfeatures<-cbind(dpfeatures,dplabel=0)
poslabel<-read.delim("nuclear receptor.txt",header = F)
poslabel<-poslabel[c(2,1)]
posrows<-length(poslabel[,1])
posindex<-c()

#---------------------------------------------------------
for(i in 1:posrows)
{
  d<-poslabel[i,1]
  p<-poslabel[i,2]
  findrow<-which((dpfeatures[,1]==as.character(d)) & (dpfeatures[,2]==as.character(p)))
  dpfeatures[findrow[1],allcol]=1
  posindex<-c(posindex,findrow[1])
}
#---------------------------------
allrowno<-c(1:length(dpfeatures[,1]))
notpos<-setdiff(allrowno,posindex)
negindex<-sample(notpos,negno)
###### in 2 khat 2-2-98 ezafeh shod
neglist<-data.frame("drugname"=dpfeatures[negindex,1],"proteinname"=dpfeatures[negindex,2])
write.table(neglist,file="negname.txt",sep="\t",row.names = F,quote = F)
#####################################
dpfeatures[negindex,allcol]=-1
randomsubdpfeatures<-subset(dpfeatures,dplabel==1 | dplabel==-1)
dim(randomsubdpfeatures)
samplenum<-length(randomsubdpfeatures[,1])
samplenum

write.table(randomsubdpfeatures,file ="randomsubdpfeatures.txt",quote = F,sep = "\t",row.names = F)
#----------------        section 2 ----------------------
labeledsamples<-read.delim("randomsubdpfeatures.txt",sep="\t",header = T)
normallabeles<-mynormalize(as.matrix(labeledsamples[,-c(1:2,ncol(labeledsamples))]))
normallabeles<-cbind(labeledsamples[,1:2],normallabeles)
normallabeles<-cbind(normallabeles,labeledsamples[ncol(labeledsamples)])
normallabeles[,ncol(normallabeles)]<-as.factor(normallabeles[,ncol(normallabeles)])
samples<-normallabeles 
######### calculate S matrix ###############
samples[,1]<-as.character(samples[,1])
samples[,2]<-as.character(samples[,2])
S<-matrix(c(rep(0,samplenum*samplenum)),nr=samplenum)
for(i in 2:samplenum)
{
  #print("i=")
  #print(i)
  
  for(j in 1:(i-1))
  {
    #print("j=")
    #  print(j)
    
    
    if((samples[i,"dplabel"]==samples[j,"dplabel"]) && (dcluster[samples[i,1]]==dcluster[samples[j,1]] || dcluster2[samples[i,1]]==dcluster2[samples[j,1]]) && (pcluster[samples[i,2]]==pcluster[samples[j,2]] || pcluster2[samples[i,2]]==pcluster2[samples[j,2]]))
      S[i,j]=1
    
    #else if((samples[i,"dplabel"]!=samples[j,"dplabel"])&& (drugknn[samples[i,1],samples[j,1]]==-targetknn[samples[i,2],samples[j,2]] || drugknn[samples[i,1],samples[j,1]]==-targetknn[samples[j,2],samples[i,2]]))
    # S[i,j]=0
    
    else
      S[i,j]=0
    
    
  }
}
diag(S)<-0
S<-as.matrix(forceSymmetric(S,uplo = 'L'))
write.table(S,file ="Smatrix25nn.txt",quote = F,sep = "\t")

##### calculate laplacian matrix ##################
yekmat<-matrix(rep(1,samplenum),nr=samplenum)
D<-S%*%yekmat
D<-diag(as.vector(D),samplenum,samplenum)
L<-D-S
######  calculate Lr  ###############################
Lr<-c(rep(0,(allcol-3)))
names(Lr)<-colnames(samples)[3:(allcol-1)]
for(i in 3:(allcol-1))
{
  # print(i)
  fr<-samples[,i]
  frT<-t(fr)
  yekT<-t(yekmat)
  frtild<-fr-((frT%*%D)%*%yekmat)/((yekT%*%D)%*%yekmat)
  frtildT<-t(frtild)
  Lr[i-2]<-((frtildT%*%L)%*%frtild)/((frtildT%*%D)%*%frtild)
}


names(Lr)<-c(1:length(Lr))
sortfeaturesASC<-sort(Lr)

write.table(sortfeaturesASC,file ="sortfeaturesASC.txt",quote = F,sep = "\t")
