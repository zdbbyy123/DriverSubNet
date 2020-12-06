load("tumor_exp.rda")
load("mutation.rda")
mutation<-as.matrix(mutation)
#To obtain the intersectsect samples of tumor expression samples and mutation samples.
common_sample<-intersect(colnames(tumor_exp),colnames(mutation))
BRCA_reseqtrans<-tumor_exp[,common_sample]
mutation<-mutation[,common_sample]

#Data process from UCSC download
BRCA_reseq_trans<-2^(BRCA_reseq)
BRCA_reseq_trans[1,1]
BRCA_reseq_trans<-round(BRCA_reseq_trans)-1
BRCA_reseq_trans[1,1]
load("normal_exp.rda")
BRCAnormal_reseq<-normal_exp
BRCAnormal_reseq_trans<-2^(BRCAnormal_reseq)
BRCAnormal_reseq_trans[1,1]
BRCAnormal_reseq_trans<-round(BRCAnormal_reseq_trans)-1
BRCAnormal_reseq_trans[1,1]

#obtain the list of differentially expressed genes
BRCA_RNAmerg<-cbind(BRCA_reseq_trans,BRCAnormal_reseq_trans)
length(colnames(BRCA_RNAmerg))
condition <- factor(c(rep("control",dim(BRCA_reseq_trans)[2]),rep("treat",dim(BRCAnormal_reseq_trans)[2])), levels = c("control","treat"))
BRCAcolData <- data.frame(row.names=colnames(BRCA_RNAmerg), condition)
library(DESeq2)
dds<-DESeqDataSetFromMatrix(BRCA_RNAmerg, BRCAcolData, design= ~ condition)
dds<-DESeq(dds)
res=results(dds, contrast=c("condition", "control", "treat"))
res=res[order(res$padj),]
res1<-as.matrix(res)
DEgene<-which(res$padj<0.05,arr.ind=T)
DE_mRNA<-rownames(res)[DEgene]


#enrichment analysis for each submetwork
load("HPRD_network.rda")
HPRD_net<-HPRD_network
HPRD_net<-as.matrix(HPRD_net)
FGs<-read.csv("FGs.csv",header=FALSE)
FGs1<-as.matrix(FGs)
result1_HPRD<-intersect(colnames(HPRD_net),DE_mRNA)
a<-length(result1_HPRD)
b<-length(intersect(result1_HPRD,FGs1))
c<-(a-b)
num<-rep(0,nrow(HPRD_net))
nei<-rep(0,nrow(HPRD_net))
DE_mRNAfoldsum<-rep(0,nrow(HPRD_net))
pvalue<-rep(0,nrow(HPRD_net))
for (i in 1:nrow(HPRD_net))
{
  a=HPRD_net[i,]
  name<-which(a==1,arr.ind=T)
  h<-colnames(HPRD_net)[name]
  k<-intersect(DE_mRNA,h)
  DE_mRNAfold<-res1[k,2]
  DE_mRNAfoldsum[i]<-sum(abs(DE_mRNAfold))
  f<-intersect(FGs1,k)
  num[i]<-length(f)
  nei[i]<-length(unique(k))
  cha<-nei[i]-num[i]
  x= c(num[i],cha,b,c)
alle<-matrix(x, nrow=2)
pvalue[i]<-fisher.test(alle)$p.value
}
BRCAnetscore<-cbind(colnames(HPRD_net),num,nei,pvalue,DE_mRNAfoldsum)

#to obtain the ESg score
BRCAnetscore<-cbind(num,nei,pvalue,DE_mRNAfoldsum)
BRCAnetscore<-apply(BRCAnetscore,2,as.numeric)
rownames(BRCAnetscore)<-colnames(HPRD_net)
BRCAnetscore[1,]
BRCAenrichscorenorm<-(BRCAnetscore[,4]-min(BRCAnetscore[,4]))/(max(BRCAnetscore[,4])-min(BRCAnetscore[,4]))
BRCAenrichscorenorm<-as.matrix(BRCAenrichscorenorm)

#normalizes mutatio data to 0-1, which acted as mutational score.
mutatedsum<-apply(mutation,1,sum)
mutatedsum<-as.matrix(mutatedsum)
mutatnorm<-(mutatedsum[,1]-min(mutatedsum[,1]))/(max(mutatedsum[,1])-min(mutatedsum[,1]))
mutatnorm<-as.matrix(mutatnorm)


#To calculate the driver gene score and priorize the driver genes.
sigGene<-which(BRCAnetscore[,3]<5e-6,arr.ind=T)
sigGene<-rownames(BRCAenrichscorenorm)[sigGene]
commongene<-intersect(sigGene,rownames(mutatnorm))
BRCAenrichscorenorm_order<-BRCAenrichscorenorm[commongene,]
mutatnorm_order<-mutatnorm[commongene,]
combscorematrix<-cbind(BRCAenrichscorenorm_order,mutatnorm_order)
combscore<-apply(combscorematrix,1,sum)
combscore<-as.matrix(combscore)
combscore<-combscore/2
ranking<-matrix(1:length(commongene),ncol=1)
generank<-cbind(ranking,commongene,combscore)
generank<-as.matrix(generank)
colnames(generank)<-c("rank","gene","gene score")
generank<-generank[order(generank[,3], decreasing= T),]
rownames(generank)<-NULL
generank[,1]<-matrix(1:length(commongene),ncol=1)
write.csv(generank, file="generank.csv")
