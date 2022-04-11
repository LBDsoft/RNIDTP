library(Matrix)
library(kmed)
library(pROC)
library(class)
library(caret)
library(kernlab)

setwd("f:/finalresults/nuclear/newnuclear")
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

pclusterfile2<-read.delim(file="pcluster2.txt",sep="\t",header = T)
pcluster2<-pclusterfile2[,1]
names(pcluster2)<-rownames(pclusterfile2)

#------------------------------------------------------
dpcluster<-c(rep(0,dpclustdim))
dim(dpcluster)<-c(dclustno,pclustno)
rownames(dpcluster)<-paste("dcluster",1:dclustno,sep="")
colnames(dpcluster)<-paste("pcluster",1:pclustno,sep="")

x=read.delim("nuclear receptor.txt",header = F)
d=read.delim("nucleardescriptors.txt",header = T)
p=read.delim("nuclearfeatures.txt",header = T)
dpinteract<-c(rep(0,dnum*pnum))
dim(dpinteract)<-c(dnum,pnum)
colnames(dpinteract)<-p[,1]
rownames(dpinteract)<-d[,1]
for( i in 1:length(x[,1]))
{
  dpinteract[as.character(x[i,2]),as.character(x[i,1])]=1
}
#----------------------------------------------------------------------
for(i in 1:dnum)
{
  for(j in 1:pnum)
  {
    if(dpinteract[i,j]==1)
    {
      dpcluster[dcluster[rownames(dpinteract)[i]],pcluster[colnames(dpinteract)[j]]]=dpcluster[dcluster[rownames(dpinteract)[i]],pcluster[colnames(dpinteract)[j]]]+1
    }
  }
}  
#------------------------
maxclustinteraction<-matrix(0, nrow = dclustno, ncol = pclustno)
for (i in 1:dclustno) {
  x1=length(which(dcluster==i))
  
  for(j in 1:pclustno)
  {
    x2=length(which(pcluster==j)) 
    maxclustinteraction[i,j]=x1*x2
  }
}

for (i in 1:dclustno) 
{
  for(j in 1:pclustno)
    dpcluster[i,j]=dpcluster[i,j]/maxclustinteraction[i,j]
}
#--------------------------------------------------------------------------------
maxclusinteract<-max(dpcluster)
for(i in 1:dnum)
{
  for(j in 1:pnum)
  {
    if(dpinteract[i,j]==0)
    {
      dpinteract[i,j]=dpcluster[dcluster[rownames(dpinteract)[i]],pcluster[colnames(dpinteract)[j]]]-maxclusinteract
    }
  }
}  
write.table(dpinteract,file ="dpinteractkmedoid.txt",quote = F,sep = "\t")

#-----------------------------------------------------------------------------
# repeat code for kmeans clustering
#------------------------------------------------------
dpcluster2<-c(rep(0,dpclustdim))
dim(dpcluster2)<-c(dclustno,pclustno)
rownames(dpcluster2)<-paste("dcluster",1:dclustno,sep="")
colnames(dpcluster2)<-paste("pcluster",1:pclustno,sep="")

#------------------------------------------
x=read.delim("nuclear receptor.txt",header = F)
d=read.delim("nucleardescriptors.txt",header = T)
p=read.delim("nuclearfeatures.txt",header = T)
dpinteract2<-c(rep(0,dnum*pnum))
dim(dpinteract2)<-c(dnum,pnum)
colnames(dpinteract2)<-p[,1]
rownames(dpinteract2)<-d[,1]
for( i in 1:length(x[,1]))
{
  dpinteract2[as.character(x[i,2]),as.character(x[i,1])]=1
}
#----------------------------------------------------------------------
for(i in 1:dnum)
{
  for(j in 1:pnum)
  {
    if(dpinteract2[i,j]==1)
    {
      dpcluster2[dcluster2[rownames(dpinteract2)[i]],pcluster2[colnames(dpinteract2)[j]]]=dpcluster2[dcluster2[rownames(dpinteract2)[i]],pcluster2[colnames(dpinteract2)[j]]]+1
    }
  }
}  
#------------------------
maxclustinteraction2<-matrix(0, nrow = dclustno, ncol = pclustno)
for (i in 1:dclustno) {
  x1=length(which(dcluster2==i))
  
  for(j in 1:pclustno)
  {
    x2=length(which(pcluster2==j)) 
    maxclustinteraction2[i,j]=x1*x2
  }
}

for (i in 1:dclustno) 
{
  for(j in 1:pclustno)
    dpcluster2[i,j]=dpcluster2[i,j]/maxclustinteraction2[i,j]
}

#--------------------------------------------------------------------------------
maxclusinteract2<-max(dpcluster2)
for(i in 1:dnum)
{
  for(j in 1:pnum)
  {
    if(dpinteract2[i,j]==0)
    {
      dpinteract2[i,j]=dpcluster2[dcluster2[rownames(dpinteract2)[i]],pcluster2[colnames(dpinteract2)[j]]]-maxclusinteract2
    }
  }
}  
write.table(dpinteract2,file ="dpinteractkmeans.txt",quote = F,sep = "\t")

#---------------------------------------------
#code for sum of 2 dpinteraction matrix
#---------------------------------------------
dpinteraction<-dpinteract+dpinteract2
write.table(dpinteraction,file ="dpinteractionkmeanskmedoid.txt",quote = F,sep = "\t")

#----------------------------------------------------------------------------
dpvector<-sort(as.vector(dpinteraction))
mindpinteract<-dpvector[3*negno]
for(i in 1:dnum)
{
  for(j in 1:pnum)
  {
    if(dpinteraction[i,j]>0)
      
      dpinteraction[i,j]=1
    else if(dpinteraction[i,j]<=mindpinteract) 
      dpinteraction[i,j]=-1
    else
      dpinteraction[i,j]=0
    
  }
}  
write.table(dpinteraction,file ="dpinteractionkmeanskmedoid.txt",quote = F,sep = "\t")
#------------------------------------------------------------------
#------------------------------------------------
newrow<-matrix(c(rep(0,pnum)),nr=1)
dpinteraction<-rbind(dpinteraction,newrow)
newcol<-matrix(c(rep(0,dnum+1)),nr=dnum+1 )
dpinteraction<-cbind(dpinteraction,newcol)
#-----------------------------------------------------------
for(i in 1:pnum)
{
  dpinteraction[dnum+1,i]<-length(which(dpinteraction[,i]==-1))
}

for(i in 1:dnum)
{
  dpinteraction[i,pnum+1]<-length(which(dpinteraction[i,]==-1))
}
dpinteraction[dnum+1,]<-LICORS::normalize(dpinteraction[dnum+1,],byrow = TRUE)
dpinteraction[,pnum+1]<-LICORS::normalize(dpinteraction[,pnum+1],byrow = F)
#-----------------------------------------------------------------------
for(i in 1:dnum)
{
  for(j in 1:pnum)
  {
    if(dpinteraction[i,j]==-1)
    {
      dpinteraction[i,j]=dpinteraction[dnum+1,j]*dpinteraction[i,pnum+1]
    }
  }
}
write.table(dpinteraction,file ="dpinteractionkmeanskmedoid.txt",quote = F,sep = "\t")
#-------------------------------------------------
dpinteraction<-dpinteraction[-(dnum+1),]
dpinteraction<-dpinteraction[,-(pnum+1)]
neglikely<-which((dpinteraction!=0) & (dpinteraction!=1),arr.ind = T)
probcol<-matrix(c(rep(0,length(neglikely[,1]))),byrow = T)
neglikely<-cbind(neglikely,probcol)
for(i in 1:length(neglikely[,1]))
{
  neglikely[i,3]<-dpinteraction[neglikely[i,1],neglikely[i,2]]
}
write.table(neglikely,file="neglikely.txt",sep="\t",row.names = F)


neglist<-neglikely[order(neglikely[,3],decreasing = TRUE),]
neglist<-as.data.frame(neglist)
neglist<-neglist[1:negno,1:2]
colnames(neglist)<-c("row","col")
neglist$drugname<-c(rep("",negno))
neglist$proteinname<-c(rep("",negno))
for(i in 1:length(neglist[,1]))
{
  neglist[i,3]<-rownames(dpinteraction)[neglist[i,1]]
  neglist[i,4]<-colnames(dpinteraction)[neglist[i,2]]
  
}
write.table(neglist[,3:4],file="negname.txt",sep="\t",row.names = F,quote = F)
######################################################################################
dfeatures<-read.delim("nucleardescriptors.txt",header = T)
pfeatures<-read.delim("nuclearfeatures.txt",header = T)
dpfeatures<-merge(dfeatures,pfeatures,by=NULL)
allcol<-dfnum+pfnum+1  # drugname + drug features+ target name + taeget features + label column
allcol
dpfeatures<-dpfeatures[,c(1,(dfnum+1),(2:dfnum),(dfnum+2):(allcol-1))]
#write.table(dpfeatures,file="dpfeatures.txt",sep="\t",row.names = F,quote = F)
dpfeatures<-cbind(dpfeatures,dplabel=0)
poslabel<-read.delim("nuclear receptor.txt",header = F)
neglabel<-read.delim("negname.txt",header = T)
poslabel<-poslabel[c(2,1)]
posrows<-length(poslabel[,1])
negrows<-length(neglabel[,1])
samplenum<-posrows+negrows
#---------------------------------------------------------
for(i in 1:posrows)
{
  d<-poslabel[i,1]
  p<-poslabel[i,2]
  findrow<-which((dpfeatures[,1]==as.character(d)) & (dpfeatures[,2]==as.character(p)))
  dpfeatures[findrow[1],allcol]=1
}
#---------------------------------
for(i in 1:negrows)
{
  
  d<-neglabel[i,1]
  p<-neglabel[i,2]
  findrow<-which((dpfeatures[,1]==as.character(d)) & (dpfeatures[,2]==as.character(p)))
  dpfeatures[findrow[1],allcol]=-1
}
subdpfeatures<-subset(dpfeatures,dplabel==1 | dplabel==-1)
write.table(subdpfeatures,file ="subdpfeatures.txt",quote = F,sep = "\t",row.names = F)
            #----------------        section 2 ----------------------
labeledsamples<-read.delim("subdpfeatures.txt",sep="\t",header = T)
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

write.table(S,file ="Smatrix.txt",quote = F,sep = "\t")

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

