#Data pre-processing
  #Reading data downloaded from UCSC
  #You can download the data you want from UCSC and process it in the following way, 
  #For those who want to understand the algorithmic process, a small dataset is provided in this paper damo
  #A<-read.table("TCGA-LUSC.txt",header = F)
  #B<-read.table("TCGA-LUAD.txt",header = F)
  #datanomal1<-A
  #datadisease1<-A
  #datanomal2<-B
  #datadisease2<-B
  #Separating paracancer samples from cancer samples
  #for(i in 2:length(datanomal1))
  #if(strsplit(datanomal1[1,i],"-")[[1]][4]=="01"||strsplit(as.character(datanomal1[1,i]),"-")[[1]][4]=="02"||strsplit(datanomal1[1,i],"-")[[1]][4]=="06")
  #  datanomal1<-datanomal1[-i]

  #for(i in 2:length(datadisease1))
  #if(strsplit(datadisease1[1,i],"-")[[1]][4]=="11")
  #  datadisease1<-datadisease1[-i] 
  #Delete genes with garbled gene names
  #genename<-datadisease1[,1]
  #genename<-genename1[-1]
  #for (i in 1:length(genename)) 
  #if(substr(genename[i],1,1)=="?")
  #{    genename<-genename[-i]
  #datadisease1<- datadisease1[-i,]
  #datanomal1<- datanomal1[-i,]
  #}
  #Save processed data
  #write.table(datadisease1,"data(disease)1.txt",row.names = F,col.names = F)
  #write.table(datanomal1,"data(nomal)1.txt",row.names = F,col.names = F)
  #for(i in 2:length(datanomal2))
  #if(strsplit(datanomal2[1,i],"-")[[1]][4]=="01"||strsplit(as.character(datanomal2[1,i]),"-")[[1]][4]=="02"||strsplit(datanomal2[1,i],"-")[[1]][4]=="06")
  #  datanomal2<-datanomal2[-i]

  #for(i in 2:length(datadisease2))
  #if(strsplit(datadisease2[1,i],"-")[[1]][4]=="11")
  #  datadisease2<-datadisease2[-i] 
  #Delete genes with garbled gene names
  #genename<-datadisease1[,1]
  #genename<-genename1[-1]
  #for (i in 1:length(genename)) 
  #if(substr(genename[i],1,1)=="?")
  #{    genename<-genename[-i]
  #datadisease2<- datadisease2[-i,]
  #datanomal2<- datanomal2[-i,]
  #}
  # #Save processed data
  # write.table(datadisease1,"data(disease)2.txt",row.names = F,col.names = F)
  # write.table(datanomal1,"data(nomal)2.txt",row.names = F,col.names = F)
  # 
  # ## Combine disease samples from two cancers with normal samples
  # nomal1<-read.table("data(nomal)1.txt",header = F)
  # nomal2<-read.table("data(nomal)2.txt",header = F)
  # disease1<-read.table("data(disease)1.txt",header = F)
  # disease2<-read.table("data(disease)2.txt",header = F)
  # nomal2<-nomal2[,-1]
  # nomal<-cbind(nomal1,nomal2)
  # #2:1 division of the training and test sets
  # disease<-cbind(disease1[,(2/3*length(disease1))],disease2[,(2/3*length(disease2))])
  # testdisease<-cbind(disease1[,(1/3*length(disease1))],disease2[,(1/3*length(disease2))])
  # write.table(disease,"data(disease).txt",row.names = F,col.names = F)
  # write.table(nomal,"data(nomal).txt",row.names = F,col.names = F)
  # write.table(disease,"testdata(disease).txt",row.names = F,col.names = F)
  # 
  #   #Matching cancer type
  # cancer_type<-c(rep(1,length(disease1[,(2/3*length(disease1))])),rep(2,length(disease2[,(2/3*length(disease1))])))#cancer type
  #   
  # write.table(cancer_type, file = "cancer_type.txt", col.names = F, quote = F, row.names = F)#Writing the cancer type to the file
  #If it is a mix of different cancers(three or more) you can use the following code
  #cancer_type<-c(rep(1,length(disease1)),rep(2,length(disease2)),rep(3,length(disease3)))#the cancer type


#The process of building a single gene network(The Methods starts here,there is a little damo data)
  #Read data
nomal<-read.table("data(nomal).txt",header = F)
disease<-read.table("data(disease).txt",header = F)
cancer_type<-read.table("cancer_type.txt",header = F)
cancer_type<-cancer_type[[1]]
genename<-disease[,1]
genename<-genename[-1]
nomal<-nomal[-1,-1]
disease<-disease[-1,-1]
  #Calculate the vector of differences between cancer samples and all normal samples at each gene
  #(the vectors of all cancer samples are clustered into a matrix under each gene)
for (j in 1:length(genename)) {
  m<-matrix(rep(0,length(nomal)*length(disease)),nrow = length(disease))
  for(i in 1:length(disease))
    m[i,]<-c(as.numeric(disease[j,i])-as.numeric(nomal[j,]))
  
  write.table(as.matrix(m),paste("./file/",genename[j],".txt",sep = ""),sep = "\t",col.names = F,row.names = F)
}
   #Calculating the sample network under a single gene
for(j in 1:length(genename))
{
  cd<-read.table(paste("./file/",genename[j],".txt",sep = ""),header = FALSE)
  mat <- matrix(as.matrix(cd),nrow=length(cd[[1]]))
  mat7<- as.matrix(dist(mat))
  write.table(mat7, file = paste("./similar/",genename[j],".txt"), sep = "\t", col.names = F, quote = F, row.names = F)
}


#Calculation of single gene heterogeneity scores based on the pseudo-F statistic
source('Fstat.R')
Flist<-list(name="score")
for(i in 1:length(genename))
{
  weight1<-read.table(paste("./similar/",genename[i],".txt"),header = FALSE)
  obs.combine.dist<-as.dist(as.matrix(weight1))
  obs.combine.res<-Fstat(obs.combine.dist~cancer_type, method="euclidean")
  Flist[[1]][i]<-obs.combine.res$F.Model[1]
}
genname<-genename
score<-Flist[[1]]
data<-data.frame(genname,score)
data<-data[order(data$score,decreasing=TRUE),]
write.csv(data, file = paste("./Fstat(order).csv"))


##Accuracy of adding different numbers of gene clusters according to the size of the pseudo-F statistic
Fdata<-read.csv("./Fstat(order).csv",head=TRUE)
fD<-data.frame(Fdata)
##Remove genes with null Fstat values (i.e. genes where gene expression is identical in normal and cancer samples)
for(i in 1:length(fD[[3]]))
  if(is.nan(fD[i,3])|is.na(fD[i,3]))
    fD<-fD[-i,]
disease<-read.table("data(disease).txt",header = T)

accList<-c()
numList<-c()
i<-1
while(i<=length(fD[[1]]))
{
    disease6<-disease[fD[1:i,1],]
    data<-t(as.matrix(disease6[,-1]))
    Hdata<-data.frame(data)
    names(Hdata)<-fD[1:i,2]
    Hdata$Species<-cancer_type
    #Hdata$Species<-cancer_type
    Hdata.kmeans<-kmeans(Hdata[,1:i],2)
    a<-table(Hdata$Species,Hdata.kmeans$cluster)
    acc<-(a[1,1]+a[2,2])/(length(disease)-1)#分为两类的方法
    if(acc<0.5) acc<-1-acc
    #p1<-a[1,1]+max((a[2,2]+a[3,3]),(a[2,3]+a[3,2]))#分为3类的方法
    #p2<-a[1,2]+max((a[2,1]+a[3,3]),(a[2,3]+a[3,1]))
    #p3<-a[1,3]+max((a[2,1]+a[3,2]),(a[2,2]+a[3,1]))
    #p[j]<-max(p1,p2,p3)
    #acc<-max(p[1:10])/(length(disease)-1)
  accList<-c(accList,acc)
  numList<-c(numList,i)
  if(i<100)
    i=i+10
  else if(i<1000)
    i<-i+100
  else if(i>=1000)
    i<-i+1000
}
data<-data.frame(numList,accList)
write.csv(data,"./F3acc.csv")


#Building a gene-gene network
library(readr)
Number_selected_build_gene_network<-50#This number is adjusted for accuracy
fDselect<-fD[1:Number_selected_build_gene_network,]
write.csv(fDselect,"./fDselect.csv")
gennet<-matrix(data=c(rep(0,length(fDselect[[1]])*length(fDselect[[1]]))), nrow = length(fDselect[[1]]), ncol = length(fDselect[[1]]), byrow = FALSE, dimnames = NULL)
Onelist<-list()
for(i in 1:length(fDselect[[1]]))
{
  weight1<-read_delim(paste("./similar/",genename[fDselect[[1]][i]],".txt"),,col_names = F,delim = "\t")
  Onelist[[i]]<-as.matrix(weight1)
}
for(i in 1:(length(Onelist)-1))
{ for(j in (i+1):length(Onelist))
{
  A<-as.matrix(Onelist[[i]])-as.matrix(Onelist[[j]])
  d<-norm(A,"F")
  gennet[i,j]=gennet[j,i]<-d
}
}

##Network Data Normalization
normalize <- function(x) {
  x<-as.matrix(x)
  return (x / (max(x) - min(x)))
}
NGennet=as.data.frame(normalize(gennet))
write_csv(NGennet, file = "./NGennet.csv")

##Set a threshold of 0.7 and remove the edges with distances greater than 0.7 to obtain the degree distribution of the graph
NGennet<-read_csv("./NGennet.csv")
degree<-c()##Used to calculate gene node degree
chosege<-c()## used to Record the order of select Genes for clustering
for(i in 1:length(NGennet))
{m<-0
for(j in 1:length(NGennet))
{
  if(i==j)
    next
  else if(NGennet[[i]][j]<0.7) 
    m<-m+1
}
degree<-c(degree,m)
if(m==(length(NGennet)-1))
  chosege<-c(chosege,i)
}

Cgen<-read.csv("./fDselect.csv")
Cgen<-Cgen[chosege,]##Selection of genes for clustering
write.csv(Cgen,"Cgen.csv",row.names = F)

##Clustering of samples
Cgen<-read.csv("Cgen.csv")
disease<-read.table("data(disease).txt",header = T)
disease7<-disease[Cgen[[2]],]
data<-t(as.matrix(disease7[,-1]))
Hdata<-data.frame(data)
names(Hdata)<-Cgen[[3]]
Hdata$Species<-cancer_type
#There are several categories according to how many different diseases there are, here are 2 categories
Hdata.kmeans<-kmeans(Hdata[,1:length(Cgen[[1]])],2)
a<-table(Hdata$Species,Hdata.kmeans$cluster)

NGennet<-read_csv("./NGennet.csv")
NGennet<-NGennet[1:500,1:500]
x<-matrix(rep(0,length(NGennet)*length(NGennet)),nrow=length(NGennet),ncol=length(NGennet))##Creating an adjacency matrix from a distance matrix
for(i in 1:(length(NGennet)-1))
  for(j in (i+1):length(NGennet))
  {if(NGennet[i,j]<0.7)
    x[i,j]=x[j,i]<-1
  }
library(igraph)
g=graph.adjacency(x,mode="undirected",weighted=T)

v<-NGennet[lower.tri(NGennet)]

for(j in 1:10)
{for(i in 1:length(v))
  {
  if(is.na(v[i])|v[i]>=0.7)
    v<-v[-i]
  }
}

E(g)$weight<-(1-v)
l<-cluster_walktrap(g,weights=E(g)$weight,7)#Set random wander step length to 7

#Ranking of the different associations according to their average scores on the pseudo-F statistic
Sf<-read.csv("fDselect.csv")
Sf<-Sf[,-1]
avegeF<-c()
for(i in 1:length(l))
{a<-0
for(j in 1:length(l[[i]]))
  a<-(a+Sf[l[[i]][j],3])

m<-(a/(length(l[[i]])))
avegeF<-c(avegeF,m)
}
p<-l[order(avegeF,decreasing = TRUE)]

##Selection of societies based on mcc or ARI factors
m<-c()
maxmcc<-0  #MCC factor
maxRI<-0  #RI factor
nua<-0
disease<-read.table("data(disease).txt",header = T)
for(i in 1:length(p))
{
  t<-m
  m<-c(m,p[[i]])
  disease9<-disease[Sf[m,1],]
  data<-t(as.matrix(disease9[,-1]))
  Hdata<-data.frame(data)
  names(Hdata)<-Sf[m,2]
  Hdata$Species<-cancer_type
  Hdata.kmeans<-kmeans(Hdata[,1:length(m)],2)
  a<-table(Hdata.kmeans$cluster,Hdata$Species)
  mcc<-((a[1,1]*a[2,2]-a[1,2]*a[2,1])/(sqrt(a[1,1]+a[1,2])*sqrt(a[1,1]+a[2,1])*sqrt(a[2,2]+a[2,1])*sqrt(a[2,2]+a[1,2])))
  if(mcc<=maxmcc) m<-t
  else
    maxmcc<-mcc
  #The following is based on the Rand factor(RI)
  # q<-0
  # for(i in 1:3)
  #   for(j in 1:3)
  #   {if(a[i,j]>1)
  #     q<-q+(a[i,j]*(a[i,j]-1)/2)
  #   }
  # h<-((a[1,1]+a[1,2]+a[1,3])*(a[1,1]+a[1,2]+a[1,3]-1)/2)+((a[2,1]+a[2,2]+a[2,3])*(a[2,1]+a[2,2]+a[2,3]-1)/2)+((a[3,1]+a[3,2]+a[3,3])*(a[3,1]+a[3,2]+a[3,3]-1)/2)
  # v<-((a[1,1]+a[2,1]+a[3,1])*(a[1,1]+a[2,1]+a[3,1]-1)/2)+((a[1,2]+a[2,2]+a[3,2])*(a[1,2]+a[2,2]+a[3,2]-1)/2)+((a[1,3]+a[2,3]+a[3,3])*(a[1,3]+a[2,3]+a[3,3]-1)/2)
  # s<-h-q
  # x<-v-q
  #all<-(length(cancer_type[[2]])*(length(cancer_type[[2]])-1)/2)
  # ar<-all-s-x-q
  # RI<-(ar+q)/all
  #if(RI<=maxRI) 
  # m<-t
  # else
  # {maxRI<-RI
  # nua<-(nua+1)
  #}
}
write.csv(Sf[m,],"./GEN(Societies discover genes).csv")


