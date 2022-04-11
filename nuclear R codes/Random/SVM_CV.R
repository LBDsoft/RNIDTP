library(randomForest)
library(cvTools)
library(pROC)
library(e1071)
library(PRROC)
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

labeledsamples=read.delim("k:/finalresults/nuclear/randomnuclear/randomsubdpfeatures.txt",header = T)
sortfeatures=read.delim("k:/finalresults/nuclear/randomnuclear/sortfeaturesASC.txt",header =F)
normallabeles<-mynormalize(as.matrix(labeledsamples[,-c(1:2,ncol(labeledsamples))]))
normallabeles<-cbind(labeledsamples[,1:2],normallabeles)
normallabeles<-cbind(normallabeles,labeledsamples[ncol(labeledsamples)])
dim(normallabeles)
trainset<-normallabeles[,-c(1:2)]
n=round(length(normallabeles[1,])*0.1)
#n<-scan()
fs<-sortfeatures[1:n,1]
trainset<-trainset[,c(fs,ncol(trainset))]
dim(trainset)
manfi<-which(trainset[,"dplabel"]==-1)
trainset[manfi,"dplabel"]=0
samplesnum<-nrow(trainset)
trainset<-cbind(normallabeles[,c(1:2)],trainset)
delrows<-samplesnum-floor(samplesnum/10)*10 # baraye round kardan hengame cross validation
delrows
if(delrows>0)
{
delrowsmanfi<-sample(1:length(manfi),delrows,replace = F)
delrowsindex<-manfi[delrowsmanfi]
trainset<-trainset[-delrowsindex,]
samplesnum<-samplesnum-delrows
}
dim(trainset)

samplesnum<-nrow(trainset)
posfoldmat<-read.delim("k:/finalresults/nuclear/newnuclear/AUC_AUPR/ten_fold_newprogram/posfoldmat.txt",header = T)
dim(posfoldmat)
#folds<-cvFolds(samplesnum, K = 10, type = "random")
#dfcross<-data.frame(fold= folds$which,index=folds$subsets)
foldsmat<-matrix(rep(0,samplesnum),nr=samplesnum/10)
for (j in 1:10)
{
  k=1
  fld=which(posfoldmat[,3]==j)
  for (i in 1:length(fld))
  {
    d=posfoldmat[fld[i],1]
    p=posfoldmat[fld[i],2]
    finddp<-which(trainset[,1]==d & trainset[,2]==p)
    foldsmat[k,j]=finddp[1]
    k=k+1
  }
}
foldsmat
samplespool=c(1:nrow(trainset))
selectedrows=which(trainset[,ncol(trainset)]==1)
for (j in 1:10) 
  {
  k=length(which(foldsmat[,j]==0))
  negfold=sample(samplespool[-selectedrows],k)
  firstnegcell<-nrow(foldsmat)-k+1
  foldsmat[firstnegcell:nrow(foldsmat),j]=negfold
  selectedrows<-c(selectedrows,negfold)
  
}
foldsmat
trainset<-trainset[,-c(1,2)]
trainset[,ncol(trainset)]<-as.factor(trainset[,ncol(trainset)])

F1_M<-c()
Recallmetrics<-c()
precision<-c()
accuracy=c()
#AUCresults<-c()
#AUPRresults<-c()
for (i in 1:10) {
  print(i)
  testset<-trainset[foldsmat[,i],]
  Nexttrainset<-trainset[-foldsmat[,i],]
  dim(testset)
  dim(Nexttrainset)
  svm.model <- svm(dplabel ~ ., data = Nexttrainset,probability = TRUE)
  svm.pred <- predict(svm.model, testset[,-ncol(testset)],decision.values = TRUE, probability = TRUE)
  confusionmat<-table(pred = svm.pred, true = as.factor(testset[,ncol(testset)]))
  acc<-(confusionmat[1,1]+confusionmat[2,2])/nrow(testset)
  accuracy<-c(accuracy,acc)
  
  pre<-confusionmat[1,1]/(confusionmat[1,1]+confusionmat[2,1])
  precision<-c(precision,pre)
  
  rec<-pre<-confusionmat[1,1]/(confusionmat[1,1]+confusionmat[1,2])
  Recallmetrics<-c(Recallmetrics,rec)
  
  F1<-(2*pre*rec)/(pre+rec)
  F1_M<-c(F1_M,F1)
 
}
resultsrownames<-c()
for(i in 1:10)
{
  resultsrownames<-c(resultsrownames,paste("testfold",toString(i),sep = ""))
}
resultsrownames<-c(resultsrownames,"means")

dfresults<-data.frame(Accuracy=accuracy,precision=precision,recall=Recallmetrics,F1_measure=F1_M)
dfresults<-rbind(dfresults,colMeans(dfresults))
row.names(dfresults)<-resultsrownames


write.table(dfresults,file ="k:/finalresults/nuclear/randomnuclear/AUC_AUPR/ten_fold_newprogram/dfresults_SVM.txt",quote = F,sep = "\t",row.names = T,col.names = T)



