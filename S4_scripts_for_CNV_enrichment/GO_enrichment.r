library(topGO)
library(fdrtool)

# Load the CNV_Contingencies environment 
load('/home/eric/Liverpool/CNV_v2/CNV_contingencies_alt_mapq002/CNV_contingencies_alt.Rdata')

# Load all of the gene annotations
annotations.table <- read.table('/home/eric/Liverpool/AR3/geneset/gene_annotation_fullgenetable.csv', header = T, sep = '\t', quote = '', row.names = 1, stringsAsFactors = F)
# Extract the GO term annotations
GO.terms <- tapply(annotations.table$GO_terms, rownames(annotations.table), function(x) strsplit(x, ';'))

# Load the table of duplicated genes just so we can get their names
cnv.genes <- unique(unlist(CNV.encompass.list.allchrom))

# We need to give a score to every gene to tell topgo which ones are in our list of interest
gene.scores <- rep(0, length(GO.terms))
names(gene.scores) <- names(GO.terms)
gene.scores[cnv.genes] <- 1

# We define the function to determine whether a gene is significant.
cnvGenes <- function(score){
	return(score==1)
}

GOdata.BP <- new("topGOdata", ontology='BP', allGenes=gene.scores, geneSel=cnvGenes, nodeSize=5, annot=annFUN.gene2GO, gene2GO=GO.terms)
fisher.test.BP <- runTest(GOdata.BP, algorithm='classic', statistic='fisher')
p.BP <- score(fisher.test.BP)
fdr.BP <- fdrtool(p.BP, statistic='pvalue')$qval
table.BP <- GenTable(GOdata.BP, classic=fisher.test.BP, orderBy='classic', topNodes=50)

GOdata.MF <- new("topGOdata", ontology='MF', allGenes=gene.scores, geneSel=cnvGenes, nodeSize=5, annot=annFUN.gene2GO, gene2GO=GO.terms)
fisher.test.MF <- runTest(GOdata.MF, algorithm='classic', statistic='fisher')
p.MF <- score(fisher.test.MF)
fdr.MF <- fdrtool(p.MF, statistic='pvalue')$qval
table.MF <- GenTable(GOdata.MF, classic=fisher.test.MF, orderBy='classic', topNodes=50)

GOdata.CC <- new("topGOdata", ontology='CC', allGenes=gene.scores, geneSel=cnvGenes, nodeSize=5, annot=annFUN.gene2GO, gene2GO=GO.terms)
fisher.test.CC <- runTest(GOdata.CC, algorithm='classic', statistic='fisher')
p.CC <- score(fisher.test.CC)
fdr.CC <- fdrtool(p.CC, statistic='pvalue')$qval
table.CC <- GenTable(GOdata.CC, classic=fisher.test.CC, orderBy='classic', topNodes=50)

# Only the MF ontology have any significant results, and this includes GO terms attributed to P450s.

# Save the workspace
save.image('GO_enrichment.Rdata')

