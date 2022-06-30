######
rm(list=ls(all=T))

#colocar a diretoria certa
#setwd("~/cancro")
setwd("~/Downloads/cancro")


#descomentar caso n?o esteja instalado
#BiocManager::install("GOSemSim")
#BiocManager::install("org.Hs.eg.db")
require(BiSEp)
require(writexl)

#ler data-raw (dataset completo e tratado) guardado
x = readRDS("data-raw.RDS")
rownames(x) = x$gene

# apaga a primeira coluna (gene)
x = x[,-1]

# mantem apenas os triple negative
x = x[,1:30]
x = log(x,2)
#BiSEp (Bimodality subsetting expression)
# # basta correr da primeira vez
# BISEP_data <- BISEP(x)
# saveRDS(BISEP_data, 'BISEP_data.RDS')

BISEP_data = readRDS('BISEP_data.RDS')

#identificacao de genes bimodais no dataset
#biIndex cont?m os scores bimodais
#bisepIndex cont?m o p-value para non-normality e o ponto medio entre os dois modos de expressao para um gene
biIndex <- BISEP_data$BI
bisepIndex <- BISEP_data$BISEP

tmp = bisepIndex
tmp["gene"] = rownames(bisepIndex)

#graficos para os tr?s genes com os valores mais altos
plot(density(as.numeric(x["LIPC",])), main="LIPC Density Distribution")
plot(density(as.numeric(x["TAT",])), main="TAT Density Distribution")
plot(density(as.numeric(x["NMNAT3",])), main="NMNAT3 Density Distribution")

#BIGEE: Bimodality in Gene Expression Exclusivity
BIGEE_out <- BIGEE(BISEP_data, sampleType="patient_low")
View(BIGEE_out)
mp=BIGEE_out[1:5,]

#visualizacao das 5 rela??es de SL candidatas usando o expressionPlot
expressionPlot(BISEP_data, gene1="ZNF516", gene2="LOC100128098")
expressionPlot(BISEP_data, gene1="SEDLP", gene2="HOXC10")
expressionPlot(BISEP_data, gene1="VEGFC", gene2="HOXC10")
expressionPlot(BISEP_data, gene1="MARCH8", gene2="HOXC10")
expressionPlot(BISEP_data, gene1="TCTN2", gene2="HOXC10")

#FURE: Functional redundancy between synthetic lethal genes
#fOut -> todos os pares candidatos 
#fOut1 -> os 5 com maior score
#fOut <- FURE(BIGEE_out, inputType="BIGEE")
#saveRDS(fOut, 'fOut.RDS')
fOut = readRDS('fOut.RDS')

#fOut1 <- FURE(mp, inputType="BIGEE")
frPairs <- fOut$funcRedundantPairs
allPairs <- fOut$allPairs

allPairs[1,]
write_xlsx(allPairs,"allPairs.xlsx")