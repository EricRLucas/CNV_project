# This script reports the number of CNVs found, where they are found, how many genes are duplicated and how many are detox genes.

load('/home/eric/Liverpool/CNV_v2/CNV_contingencies_alt_mapq002/CNV_contingencies_alt.Rdata')

# How many CNVs did we find? 
cat('We found a total of ', nrow(good.CNV.ranges.allchrom), ' CNVs:\n',
    '\t', nrow(good.CNV.ranges.2L), ' on chromosome 2L\n',
    '\t', nrow(good.CNV.ranges.2R), ' on chromosome 2R\n',
    '\t', nrow(good.CNV.ranges.3L), ' on chromosome 3L\n',
    '\t', nrow(good.CNV.ranges.3R), ' on chromosome 3R\n',
    '\t', nrow(good.CNV.ranges.X), ' on chromosome X\n\n', sep = '')

# How were the CNVs distributed between the populations?
CNVs.bypop.table <- rbind(pop.cnv.totals.2L, pop.cnv.totals.2R, pop.cnv.totals.3L, pop.cnv.totals.3R, pop.cnv.totals.X, pop.cnv.totals) 
cat('Here are the number of CNVs found in each population:\n')
print(CNVs.bypop.table)
cat('\n\n')
write.table(CNVs.bypop.table, file = 'CNVs_bypop_table.csv', sep = '\t', quote = F, col.names = NA)

# Are the detox CNVs equally distributed among populations?
fisher.test.detox.by.pop.table <- rbind(pop.cnv.totals - pop.cnv.detox.totals, pop.cnv.detox.totals)
rownames(fisher.test.detox.by.pop.table) <- c('non-detox', 'detox')
cat('Running fisher test on the following table:\n')
print(fisher.test.detox.by.pop.table)
fisher.test.detox.by.pop <- fisher.test(fisher.test.detox.by.pop.table, simulate.p.value = T, B = 1000000)
print(fisher.test.detox.by.pop)
cat('The proportions of detox CNVs in each population are as follows:\n')
pop.cnv.detox.proportions <- pop.cnv.detox.totals / pop.cnv.totals
print(pop.cnv.detox.proportions)
colnames(fisher.test.detox.by.pop.table) <- paste(colnames(fisher.test.detox.by.pop.table), '(', pop.sizes.reduced, ')', sep = '')
write.table(rbind(fisher.test.detox.by.pop.table, pop.cnv.detox.proportions), file = 'CNVs_bypop_detox_proportions.csv', sep = '\t', quote = F, col.names = NA)

# What about the proportions of detox genes relative to genes?
fisher.test.detox.vs.gene.table <- rbind(pop.cnv.withgene.totals - pop.cnv.detox.totals, pop.cnv.detox.totals)
rownames(fisher.test.detox.vs.gene.table) <- c('non-detox gene', 'detox gene')
cat('Running fisher test on the following table:\n')
print(fisher.test.detox.vs.gene.table)
fisher.test.detox.vs.gene <- fisher.test(fisher.test.detox.vs.gene.table, simulate.p.value = T, B = 1000000)
print(fisher.test.detox.vs.gene)

# How many genes were amplified by at least one CNV (and how many were eligible), and how many of these were detox?
duplicated.genes <- unique(unlist(CNV.encompass.list.allchrom))
duplicated.detox.genes <- intersect(duplicated.genes, rownames(all.genes)[all.genes$is.detox])
cat('Out of ', length(eligible.gene.names), ' eligible genes, ', length(duplicated.genes), ' were found in at least one of the ', nrow(good.CNV.ranges.allchrom), ' CNVs.\n',
    length(duplicated.detox.genes), ' of these are detox genes.\n\n', sep = '')

# How many CNVs contained at least one gene (/ detox gene)
cat('Out of ', nrow(good.CNV.ranges.allchrom), ' CNVs, ', nrow(CNV.encompass.table.allchrom), ' contained at least one gene, \n',
    'of which ', sum(CNV.encompass.table.allchrom$has.detox), ' contained at least 1 detox gene.\n\n', sep = '')

# Load the gene annotations
annot.info <- read.table('/home/eric/Liverpool/data_tables/all_ID2name.csv', header = T, sep = '\t', row.names = 1, stringsAsFactors = F)
# Create a table recording the genes that were duplicated, there genomic position and their annotation
# If a gene name is not in the annotation table, then it gets given NA
all.eligible.genes$annotation <- annot.info[rownames(all.eligible.genes),'Gene.name']
output.gene.table <- all.eligible.genes[duplicated.genes, c('Chrom', 'start', 'end', 'is.detox', 'annotation')]
output.gene.table$annotation[output.gene.table$annotation == '' | is.na(output.gene.table$annotation)] <- 'None'
# Write that table to file
write.table(output.gene.table, file = 'duplicated_genes.csv', quote = F, sep = '\t', col.names = NA)

# We also want a table that lists each CNVs, giving the proportion of individuals in which it is found in each 
# population, and any genes that it contains. 
output.cnv.table <- good.CNV.ranges.allchrom
# Add the chromosome information
output.cnv.table$Chromosome <- sub('\\d*$', '', sub('CNV_', '', rownames(output.cnv.table)))
output.cnv.table <- output.cnv.table[,c('Chromosome', 'start', 'end')]
# Add the gene information
output.cnv.table$Genes <- CNV.encompass.table.allchrom[rownames(output.cnv.table), 'Features']
output.cnv.table$Genes[is.na(output.cnv.table$Genes)] <- ''
output.cnv.table$Genes <- sub('Genes: ', '', output.cnv.table$Genes)
output.cnv.table$has.detox <- CNV.encompass.table.allchrom[rownames(output.cnv.table), 'has.detox']
output.cnv.table$has.detox[is.na(output.cnv.table$has.detox)] <- F
# Get the population proportions for each CNV
pop.cnv.counts.allchrom <- rbind(pop.cnv.counts.2L, pop.cnv.counts.2R, pop.cnv.counts.3L, pop.cnv.counts.3R, pop.cnv.counts.X) 
pop.good.cnv.counts <- pop.cnv.counts.allchrom[rownames(good.CNV.ranges.allchrom),]
pop.good.cnv.props <- pop.good.cnv.counts / matrix(pop.sizes.reduced, nrow(pop.good.cnv.counts), ncol(pop.good.cnv.counts), byrow = T)
# Sanity check to make sure rownames match
if (any(rownames(pop.good.cnv.props) != rownames(output.cnv.table)))
	stop('Rownames do not match.')
# Add the cnv.proportions to the table
output.cnv.table <- cbind(output.cnv.table, pop.good.cnv.props)
# Write that to file
write.table(output.cnv.table, file = 'all_detected_CNVs.csv', quote = F, sep = '\t', col.names = NA)

