# NOTAS:

# data.raw -> dataset completo
# panel.raw -> (mini) dataset dos pares de interesse
# gene.panel -> lista em que cada elemento ? um array com as strings dos pares de interesse

#descomentar para instalar
#install.packages("writexl")

rm(list=ls(all=T))
library(qvalue) 

#colocar a diretoria certa
setwd("~/cancro")

#cria pasta output na diretoria 
if(!file.exists("output")) dir.create("output")

# alterar para ler a partir do disco
data.raw=readRDS('data-raw.RDS')
panel.raw=readRDS('top_pares.RDS')

# limitar apenas ?s primeiras N combina??es
# panel.raw = panel.raw[1:11,]


# 30 primeiros -> T
# 31:43 -> N
for (i in 1:30) {
  colnames(data.raw)[i+1] = paste(c("T.", as.character(i)), collapse="");
}

for (i in 1:13) {
  colnames(data.raw)[i+31] = paste(c("N.", as.character(i)), collapse="");
}



#include panel data
# inicializa gene.panel como lista vazia com espa?os
gene.panel=vector("list",dim(panel.raw)[1])

# inicializa measure.type.data como array numerico de tamanho 2
measure.type.data=numeric(dim(panel.raw)[2])

# inicializa a NULL, mas vai c? colocar os pares do "resultado"
name.string=NULL

# para cada linha i
#as.character transforma os dados em carateres  vai so buscar a vari?vel v1 #rbind junta os vetores
for(i in 1:dim(panel.raw)[1]) {
  gene.panel[[i]]=c(as.character(panel.raw$X1[i]),as.character(panel.raw$X2[i]))
  name.string=rbind(name.string,paste(as.character(panel.raw$X1[i])," & ",as.character(panel.raw$X2[i]),sep=""))
}


T.index=which(substr(names(data.raw),1,1)=="T") #vai buscar os q t?m T e tira do dataset pq vai cont?-los
N.index=which(substr(names(data.raw),1,1)=="N") 
gene.T.n=length(T.index) #v? aqui quantos s?o
gene.N.n=length(N.index)
gene.name=as.character(data.raw[,1]) #vai buscar apenas todas as linhas da primeira coluna
gene.data=data.raw[,(2:dim(data.raw)[2])] #vai aos dados e extrai as outras q sao numeros
fold.change=log2(2)   ##Setting initial value##  # valor inicial usado na razao =1 

#introduz uma fun??o que calcula a m?dia do n e o log r   normaliza?ao 
Ratio.fun=function(dataT,dataN){
  mean.N=mean(2^dataN)
  log.r=dataT-log2(mean.N)
  return(log.r)
}

# conta os q t?m gene t+n  e coloca os dados numa matriz come?a com uma nula toda a zeros e depois preenche
#n de colunas t+n mix e n de linhas ? quantas vezes isso acontece
# se o i for 1 preenche a primeira linha, conta sem reposi??o False
gene.perm.index=function(gene.T.n,gene.N.n,times){
  mix.gene.n=gene.T.n+gene.N.n
  index.m=matrix(0,ncol=mix.gene.n,nrow=times)
  for(i in 1:times) index.m[i,]=sample(1:mix.gene.n,mix.gene.n,replace=F)
  return(index.m)
}

### measure type up-up=1 up-down=2 down-up=3 down-down=4 # index.m ? a matriz o w ? o n
Permutation.s.test=function(gene1.t,gene1.w,gene2.t,gene2.w,index.m){
  gene.T.n=length(gene1.t)
  gene.N.n=length(gene1.w)
  mix.gene1=c(gene1.t,gene1.w)
  mix.gene2=c(gene2.t,gene2.w)
  mix.gene.n=length(mix.gene1)
  Ratio.1=Ratio.fun(gene1.t,gene1.w) #razao 
  Ratio.2=Ratio.fun(gene2.t,gene2.w)
  obs_prob.temp=numeric(4) # obs vai ter comprimento 4
  obs_prob.temp[1]=sum(Ratio.1>=fold.change & Ratio.2>=fold.change)/gene.T.n
  obs_prob.temp[2]=sum(Ratio.1>=fold.change & Ratio.2<=-fold.change)/gene.T.n
  obs_prob.temp[3]=sum(Ratio.1<=-fold.change & Ratio.2>=fold.change)/gene.T.n
  obs_prob.temp[4]=sum(Ratio.1<=-fold.change & Ratio.2<=-fold.change)/gene.T.n
  names(obs_prob.temp) <- c("up,up","up,down","down,up","down,down")
  obs_prob=max(obs_prob.temp)
  measure.type=which.max(obs_prob.temp)   #posi?ao indice do valor maximo
  times=dim(index.m)[1] #n de linhas
  perm_prob=numeric(times)
  for(i in 1:times){
    index=index.m[i,] #ver as linhas da matriz
    Ratio.1=Ratio.fun(mix.gene1[index[1:gene.T.n]],mix.gene1[index[(gene.T.n+1):mix.gene.n]])
    # print(Ratio.1) #vai ? linda da matriz e vai buscar os dados de 1 at? ao genes t e depois vai buscar o resto sepadado
    Ratio.2=Ratio.fun(mix.gene2[index[1:gene.T.n]],mix.gene2[index[(gene.T.n+1):mix.gene.n]])
    if (measure.type==1) perm_prob[i]=sum(Ratio.1>=fold.change & Ratio.2>=fold.change)/gene.T.n # se o valor max for na posi??o 1 vai fazer esta linha
    if (measure.type==2) perm_prob[i]=sum(Ratio.1>=fold.change & Ratio.2<=-fold.change)/gene.T.n
    if (measure.type==3) perm_prob[i]=sum(Ratio.1<=-fold.change & Ratio.2>=fold.change)/gene.T.n
    if (measure.type==4) perm_prob[i]=sum(Ratio.1<=-fold.change & Ratio.2<=-fold.change)/gene.T.n
  }
  p.value=(sum(perm_prob>=obs_prob)+1)/(times+1)
  return(list(p.hat=obs_prob.temp,p.value=p.value,perm_prob=perm_prob))
}
record.time=numeric(dim(panel.raw)[1])
B.times=9999   ##Setting initial value##
index.m=gene.perm.index(gene.T.n,gene.N.n,times=B.times) #vai buscar os daods da matriz
record.panel=vector("list",dim(panel.raw)[1])
for(i in 1:dim(panel.raw)[1]){
  gene1.index=which(gene.name==gene.panel[[i]][1])
  gene2.index=which(gene.name==gene.panel[[i]][2])
  gene1.T.data=as.numeric(data.raw[gene1.index,T.index])
  gene1.N.data=as.numeric(data.raw[gene1.index,N.index])
  gene2.T.data=as.numeric(data.raw[gene2.index,T.index])
  gene2.N.data=as.numeric(data.raw[gene2.index,N.index])
  record.panel[[i]]=Permutation.s.test(gene1.t=gene1.T.data,gene1.w=gene1.N.data,gene2.t=gene2.T.data,gene2.w=gene2.N.data,index.m=index.m)
} # teste de permuta??o binomial de propor?a?ao h0 ? serem iguais h1 n serem iguais alpha 5%


p.value=numeric(dim(panel.raw)[1])
p.hat=matrix(nrow=dim(panel.raw)[1],ncol=4) #propor??es matriz com elas
colnames(p.hat) <- c("up,up","up,down","down,up","down,down")
for(i in 1:dim(panel.raw)[1]) {
  p.value[i]=record.panel[[i]]$p.value
  p.hat[i,]=record.panel[[i]]$p.hat
} #preencher a matriz

#  FDR.value=p.adjust(p.value,method="BH")
# ORIGINAL: FDR.value=qvalue(p.value)$qvalues  
FDR.value=qvalue(p.value, pi0=1)$qvalues  
output.data=cbind(p.hat,p.vaule=p.value,FDR=FDR.value)
rownames(output.data)=name.string

#cria ficheiro com o resultado na pasta output
write.csv(output.data,file=paste("output/1-output of fractions of gene expression patterns.csv",sep=""))
require(writexl)
#descomentar para obter ficheiro xlsx
#output.data=as.data.frame(output.data)
#write_xlsx(output.data,"paresSL.xlsx")
