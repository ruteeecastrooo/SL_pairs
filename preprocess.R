#remove todas as variáveis do ambiente
rm(list=ls(all=T))

#packages correr se n?o tiver instalado
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.14")

#BiocManager::install("Biobase")
#BiocManager::install("BiocGenerics")
#BiocManager::install("GEOquery")
#BiocManager::install("qvalue")

library(GEOquery)

#colocar a diretoria certa
setwd("~/Downloads/cancro")

#fun??es auxiliares para, a partir da posicao na matriz (de tamanho n*n), 
#dizer qual ? o gene correto a que aquela linha/coluna corresponde

get_nome_coluna_por_indice <- function(corrs, indice, n) {
  posicao_coluna = indice %% n;
  if(posicao_coluna == 0) {
    posicao_coluna = 10;
  }
  colnames(corrs)[posicao_coluna]
} 

get_nome_linha_por_indice <- function(corrs, indice, n) {
  resto_divisao = indice %% n;
  if ( resto_divisao == 0) {
    posicao_linha = floor(indice / n);
  }else {
    posicao_linha = floor(indice / n) + 1;
  }
  colnames(corrs)[posicao_linha]
} 


#acesso ao nosso dataset
GSE38959<-getGEO("GSE38959", destdir=".")

info = GSE38959[["GSE38959_series_matrix.txt.gz"]]@featureData@data

x = data.frame(GSE38959[["GSE38959_series_matrix.txt.gz"]]@assayData[["exprs"]])

# eliminar ?ltimas 4 colunas porque sao amostras de cora??o, pulm?o, f?gado, rim
x = x[,-44:-47]

# nome dos genes
gene = info$GENE_SYMBOL

# caso existam genes sem nome no top, ? necess?rio ver pelo id:
# gene = GSE38959[["GSE38959_series_matrix.txt.gz"]]@featureData@data[["ID"]]

#juntar o nome dos gens ao dataframe
x = cbind(gene, x)

# apagar as linhas que n?o tenham nome do gene (que veio do GENE_SYMBOL)
x = x[x$gene != "",]

# apaga genes duplicados
x = x[!duplicated(x[,c('gene')]),]

# descomentar isto s? para correr os primeiros 1000 caso esteja muito lento!!!
# n = 1000
n = dim(x)[1]

# usar apenas n linhas:
x = x[1:n,]

print("A calcular correla??es...")
corrs = cor(t(x[-1]))

# truque para n?o devolver a diagonal!
for(i in 1:n) {
  corrs[i,i] = 0;
}

colnames(corrs) = x$gene[1:n]
rownames(corrs) = x$gene[1:n]


num_top_a_reter = 50

print("A ordenar resultados...")
posicoes_top = order(corrs, decreasing = TRUE)[1:(num_top_a_reter*2)]

# resultado = data.frame(matrix(ncol = 2, nrow = num_top_a_reter), V1 = character(), V2 = character(), stringsAsFactors = FALSE);
resultado = data.frame(matrix(ncol = 2, nrow = num_top_a_reter));

# vai buscar os nomes dos genes da coluna/linha correspondentes ?s posi??es top
for(i in 1:num_top_a_reter*2) {
  if (i %% 2 == 0) {
    k = floor(i / 2);
    posicao_atual = posicoes_top[i];
    nome_linha = get_nome_linha_por_indice(corrs, posicao_atual, n);
    nome_coluna = get_nome_coluna_por_indice(corrs, posicao_atual, n);
    resultado[k,] = c(nome_linha, nome_coluna);
  }
}

# guarda a lista para dos pares para n?o termos que correr isto de novo
print("Saving top_pares.RDS")
saveRDS(resultado, 'top_pares.RDS')
print("Saving data-raw.RDS")
saveRDS(x, 'data-raw.RDS')

