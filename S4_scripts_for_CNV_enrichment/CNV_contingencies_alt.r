# The point of this script is to get some objective measure of the number of genes that have a duplication
# at an appreciable frequency in each population, and then compare this with the same data using just 
# metabolic genes, and then just the four regions of interest. 

# We load the list of detox genes. 
detox.genes <- read.table('/home/eric/Manuscripts/GSTE/CNV_enrichment/detox_genes_output2.txt', as.is = T)[,1]

# Now list the AGAP numbers for the focal genes of interest
focal.genes <- list('GST' = character(), 'CYP6P' = character(), 'CYP9K1' = character(), 'CYP6M' = character())
focal.genes[['GST']]['GSTE1'] <- 'AGAP009195'
focal.genes[['GST']]['GSTE2'] <- 'AGAP009194'
focal.genes[['GST']]['GSTE3'] <- 'AGAP009197'
focal.genes[['GST']]['GSTE4'] <- 'AGAP009193'
focal.genes[['GST']]['GSTE5'] <- 'AGAP009192'
focal.genes[['GST']]['GSTE6'] <- 'AGAP009191'
focal.genes[['GST']]['GSTE7'] <- 'AGAP009196'
focal.genes[['GST']]['GSTU4'] <- 'AGAP009190'
focal.genes[['CYP6P']]['CYP6AA1'] <- 'AGAP002862'
focal.genes[['CYP6P']]['CYP6AA2'] <- 'AGAP013128'
focal.genes[['CYP6P']]['COEAE6O'] <- 'AGAP002863'
focal.genes[['CYP6P']]['CYP6P1'] <- 'AGAP002868'
focal.genes[['CYP6P']]['CYP6P2'] <- 'AGAP002869'
focal.genes[['CYP6P']]['CYP6P3'] <- 'AGAP002865'
focal.genes[['CYP6P']]['CYP6P4'] <- 'AGAP002867'
focal.genes[['CYP6P']]['CYP6P5'] <- 'AGAP002866'
focal.genes[['CYP6P']]['CYP6AD1'] <- 'AGAP002870'
focal.genes[['CYP9K1']]['CYP9K1'] <- 'AGAP000818'
focal.genes[['CYP6M']]['CYP6M2'] <- 'AGAP008212'
focal.genes[['CYP6M']]['CYP6M3'] <- 'AGAP008213'
focal.genes[['CYP6M']]['CYP6M4'] <- 'AGAP008214'
focal.genes[['CYP6M']]['CYP6Z1'] <- 'AGAP008219'
focal.genes[['CYP6M']]['CYP6Z2'] <- 'AGAP008218'
focal.genes[['CYP6M']]['CYP6Z3'] <- 'AGAP008217'

# Get the names and positions of all genes in the genome
all.genes <- read.table('/home/eric/Liverpool/CNV_/find_gene_regions_output.csv', header = T, row.names = 1)
all.genes$is.detox <- rownames(all.genes) %in% detox.genes
all.genes$covered.windows <- 0
all.genes$full.windows <- 0

# Get the annotation information
annot.info <- read.table('/home/eric/Liverpool/data_tables/all_ID2name.csv', header = T, sep = '\t', row.names = 1, stringsAsFactors = F)

# Write a function to plot the data for a given gene
pp <- function(gene.name, list.of.samples = names(counts.list)){
	spos <- all.genes[gene.name, 'start']
	epos <- all.genes[gene.name, 'end']
	sindex <- min(which(counts.list[[1]]$Position >= spos)) - 50
	if (sindex < 1)
		sindex <- 1
	eindex <- max(which(counts.list[[1]]$Position <= epos)) + 50
	if (eindex > max(counts.list[[1]]$Position))
		eindex <- max(counts.list[[1]]$Position)
	x.midpoint <- mean(spos, epos)
	i <- 1
	while(1){
		if (i < 1)
			i <- length(list.of.samples)
		this.sample <- list.of.samples[i]
		plot(counts.list[[this.sample]]$Position[sindex:eindex], counts.list[[this.sample]]$Normalised_coverage[sindex:eindex], main = this.sample, ylim = c(0,12)) 
		these.cnv.states <- counts.list[[this.sample]]$CNV[sindex:eindex]
		lines(counts.list[[this.sample]]$Position[sindex:eindex], these.cnv.states, col = 2)
		abline(v = spos)
		abline(v = epos)
		text(mean(c(spos, epos)), 9, gene.name, srt = 90, adj = 0)
		x <- locator(1)$x
		if (x <= x.midpoint)
			i <- i-1
		else
			i <- i+1
	}
}

# Load the results from 2L. 
load('/home/eric/Desktop/Liverpool/CNV_v2_bigfiles/counting_output_v4_2L/phase2/alt_fullgcnorm_nomapq_mapq002_varvar_trans000001/CNV_analysis_alt_2L.Rdata')
encompass.output.by.gene.2L <- encompass.output.by.gene
counts.list.positions.2L <- counts.list.positions
eligible.gene.position.indices.2L <- eligible.gene.position.indices
good.CNV.ranges.2L <- good.CNV.ranges
pop.cnv.2L <- pop.cnv
pop.cnv.counts.2L <- pop.cnv.counts
clustered.list.2L <- clustered.list
duplications.by.cluster.2L <- duplications.by.cluster
CNV.encompass.table.2L <- CNV.encompass.table
CNV.encompass.list.2L <- CNV.encompass.list

# Let's get the sum of the pop.cnv columns 
pop.cnv.totals.2L <- apply(pop.cnv.2L, 2, sum)

# As a sanity check, get the number of windows covering each gene and check this matches with the
# eligible.gene.positions.indices table
for (this.gene in rownames(all.genes)[all.genes$Chrom == '2L']){
	# Get the number of windows that are fully covered by this gene
	this.gene.covered.windows <- sum((counts.list.positions.2L >= all.genes[this.gene, 'start']) & (counts.list.positions.2L <= all.genes[this.gene, 'end'] - coverage.window))
	# Get the number of windows that would be fully covered by this gene without filtering
	this.gene.full.windows <- floor(all.genes[this.gene, 'end']/coverage.window) - ceiling(all.genes[this.gene, 'start']/coverage.window)
	# Check this is consistent with our previous calculation
	this.gene.full.windows[this.gene.full.windows < 0] <- 0
	if ((is.na(eligible.gene.position.indices.2L[this.gene, 'start']) & (this.gene.covered.windows >= (min.win.cov * this.gene.full.windows) & this.gene.full.windows != 0)) | (!(is.na(eligible.gene.position.indices.2L[this.gene, 'start'])) & (this.gene.covered.windows < (min.win.cov * this.gene.full.windows) | this.gene.full.windows == 0)))
		stop('Window coverage data inconsistent.')
	else{
		all.genes[this.gene, 'covered.windows'] <- this.gene.covered.windows
		all.genes[this.gene, 'full.windows'] <- this.gene.full.windows
	}
}

# Now let's do the same with the other chromosomes
#2R
load('/home/eric/Desktop/Liverpool/CNV_v2_bigfiles/counting_output_v4_2R/phase2/alt_fullgcnorm_nomapq_mapq002_varvar_trans000001/CNV_analysis_alt_2R.Rdata')
encompass.output.by.gene.2R <- encompass.output.by.gene
counts.list.positions.2R <- counts.list.positions
eligible.gene.position.indices.2R <- eligible.gene.position.indices
good.CNV.ranges.2R <- good.CNV.ranges
pop.cnv.2R <- pop.cnv
pop.cnv.counts.2R <- pop.cnv.counts
clustered.list.2R <- clustered.list
duplications.by.cluster.2R <- duplications.by.cluster
CNV.encompass.table.2R <- CNV.encompass.table
CNV.encompass.list.2R <- CNV.encompass.list
pop.cnv.totals.2R <- apply(pop.cnv.2R, 2, sum)
for (this.gene in rownames(all.genes)[all.genes$Chrom == '2R']){
	this.gene.covered.windows <- sum((counts.list.positions.2R >= all.genes[this.gene, 'start']) & (counts.list.positions.2R <= all.genes[this.gene, 'end'] - coverage.window))
	this.gene.full.windows <- floor(all.genes[this.gene, 'end']/coverage.window) - ceiling(all.genes[this.gene, 'start']/coverage.window)
	this.gene.full.windows[this.gene.full.windows < 0] <- 0
	if ((is.na(eligible.gene.position.indices.2R[this.gene, 'start']) & (this.gene.covered.windows >= (min.win.cov * this.gene.full.windows) & this.gene.full.windows != 0)) | (!(is.na(eligible.gene.position.indices.2R[this.gene, 'start'])) & (this.gene.covered.windows < (min.win.cov * this.gene.full.windows) | this.gene.full.windows == 0)))
		stop('Window coverage data inconsistent.')
	else{
		all.genes[this.gene, 'covered.windows'] <- this.gene.covered.windows
		all.genes[this.gene, 'full.windows'] <- this.gene.full.windows
	}
}

#3L
load('/home/eric/Desktop/Liverpool/CNV_v2_bigfiles/counting_output_v4_3L/phase2/alt_fullgcnorm_nomapq_mapq002_varvar_trans000001/CNV_analysis_alt_3L.Rdata')
encompass.output.by.gene.3L <- encompass.output.by.gene
counts.list.positions.3L <- counts.list.positions
eligible.gene.position.indices.3L <- eligible.gene.position.indices
good.CNV.ranges.3L <- good.CNV.ranges
pop.cnv.3L <- pop.cnv
pop.cnv.counts.3L <- pop.cnv.counts
clustered.list.3L <- clustered.list
duplications.by.cluster.3L <- duplications.by.cluster
CNV.encompass.table.3L <- CNV.encompass.table
CNV.encompass.list.3L <- CNV.encompass.list
pop.cnv.totals.3L <- apply(pop.cnv.3L, 2, sum)
for (this.gene in rownames(all.genes)[all.genes$Chrom == '3L']){
	this.gene.covered.windows <- sum((counts.list.positions.3L >= all.genes[this.gene, 'start']) & (counts.list.positions.3L <= all.genes[this.gene, 'end'] - coverage.window))
	this.gene.full.windows <- floor(all.genes[this.gene, 'end']/coverage.window) - ceiling(all.genes[this.gene, 'start']/coverage.window)
	this.gene.full.windows[this.gene.full.windows < 0] <- 0
	if ((is.na(eligible.gene.position.indices.3L[this.gene, 'start']) & (this.gene.covered.windows >= (min.win.cov * this.gene.full.windows) & this.gene.full.windows != 0)) | (!(is.na(eligible.gene.position.indices.3L[this.gene, 'start'])) & (this.gene.covered.windows < (min.win.cov * this.gene.full.windows) | this.gene.full.windows == 0)))
		stop('Window coverage data inconsistent.')
	else{
		all.genes[this.gene, 'covered.windows'] <- this.gene.covered.windows
		all.genes[this.gene, 'full.windows'] <- this.gene.full.windows
	}
}

#3R
load('/home/eric/Desktop/Liverpool/CNV_v2_bigfiles/counting_output_v4_3R/phase2/alt_fullgcnorm_nomapq_mapq002_varvar_trans000001/CNV_analysis_alt_3R.Rdata')
encompass.output.by.gene.3R <- encompass.output.by.gene
counts.list.positions.3R <- counts.list.positions
eligible.gene.position.indices.3R <- eligible.gene.position.indices
good.CNV.ranges.3R <- good.CNV.ranges
pop.cnv.3R <- pop.cnv
pop.cnv.counts.3R <- pop.cnv.counts
clustered.list.3R <- clustered.list
duplications.by.cluster.3R <- duplications.by.cluster
CNV.encompass.table.3R <- CNV.encompass.table
CNV.encompass.list.3R <- CNV.encompass.list
pop.cnv.totals.3R <- apply(pop.cnv.3R, 2, sum)
for (this.gene in rownames(all.genes)[all.genes$Chrom == '3R']){
	this.gene.covered.windows <- sum((counts.list.positions.3R >= all.genes[this.gene, 'start']) & (counts.list.positions.3R <= all.genes[this.gene, 'end'] - coverage.window))
	this.gene.full.windows <- floor(all.genes[this.gene, 'end']/coverage.window) - ceiling(all.genes[this.gene, 'start']/coverage.window)
	this.gene.full.windows[this.gene.full.windows < 0] <- 0
	if ((is.na(eligible.gene.position.indices.3R[this.gene, 'start']) & (this.gene.covered.windows >= (min.win.cov * this.gene.full.windows) & this.gene.full.windows != 0)) | (!(is.na(eligible.gene.position.indices.3R[this.gene, 'start'])) & (this.gene.covered.windows < (min.win.cov * this.gene.full.windows) | this.gene.full.windows == 0)))
		stop('Window coverage data inconsistent.')
	else{
		all.genes[this.gene, 'covered.windows'] <- this.gene.covered.windows
		all.genes[this.gene, 'full.windows'] <- this.gene.full.windows
	}
}

#X
load('/home/eric/Desktop/Liverpool/CNV_v2_bigfiles/counting_output_v4_X/phase2/alt_fullgcnorm_nomapq_mapq002_varvar_trans000001/CNV_analysis_alt_X.Rdata')
encompass.output.by.gene.X <- encompass.output.by.gene
counts.list.positions.X <- counts.list.positions
eligible.gene.position.indices.X <- eligible.gene.position.indices
good.CNV.ranges.X <- good.CNV.ranges
pop.cnv.X <- pop.cnv
pop.cnv.counts.X <- pop.cnv.counts
clustered.list.X <- clustered.list
duplications.by.cluster.X <- duplications.by.cluster
CNV.encompass.table.X <- CNV.encompass.table
CNV.encompass.list.X <- CNV.encompass.list
pop.cnv.totals.X <- apply(pop.cnv.X, 2, sum)
for (this.gene in rownames(all.genes)[all.genes$Chrom == 'X']){
	this.gene.covered.windows <- sum((counts.list.positions.X >= all.genes[this.gene, 'start']) & (counts.list.positions.X <= all.genes[this.gene, 'end'] - coverage.window))
	this.gene.full.windows <- floor(all.genes[this.gene, 'end']/coverage.window) - ceiling(all.genes[this.gene, 'start']/coverage.window)
	this.gene.full.windows[this.gene.full.windows < 0] <- 0
	if ((is.na(eligible.gene.position.indices.X[this.gene, 'start']) & (this.gene.covered.windows >= (min.win.cov * this.gene.full.windows) & this.gene.full.windows != 0)) | (!(is.na(eligible.gene.position.indices.X[this.gene, 'start'])) & (this.gene.covered.windows < (min.win.cov * this.gene.full.windows) | this.gene.full.windows == 0)))
		stop('Window coverage data inconsistent.')
	else{
		all.genes[this.gene, 'covered.windows'] <- this.gene.covered.windows
		all.genes[this.gene, 'full.windows'] <- this.gene.full.windows
	}
}

# Combine the CNV lists for all chromosomes
pop.cnv.allchrom <- rbind(pop.cnv.2L,
                          pop.cnv.2R,
                          pop.cnv.3L,
                          pop.cnv.3R,
                          pop.cnv.X)
good.CNV.ranges.allchrom <- rbind(good.CNV.ranges.2L,
                                  good.CNV.ranges.2R,
                                  good.CNV.ranges.3L,
                                  good.CNV.ranges.3R,
                                  good.CNV.ranges.X)
CNV.encompass.list.allchrom <- c(CNV.encompass.list.2L,
                                 CNV.encompass.list.2R,
                                 CNV.encompass.list.3L,
                                 CNV.encompass.list.3R,
                                 CNV.encompass.list.X)
CNV.encompass.table.allchrom <- rbind(CNV.encompass.table.2L,
                                      CNV.encompass.table.2R,
                                      CNV.encompass.table.3L,
                                      CNV.encompass.table.3R,
                                      CNV.encompass.table.X)

# Find out which CNVs contain detox genes. 
has.detox <- unlist(lapply(CNV.encompass.list.allchrom, function(x) any(x %in% detox.genes)))
# Sanity check. Verify that the names of is.detox are the same as the rownames of CNV.encompass.list.allchrom
if (!all(rownames(CNV.encompass.table.allchrom) == names(has.detox)))
	stop("Fail. The names of has.detox don't match those of CNV.encompass.list.allchrom.")
CNV.encompass.table.allchrom$has.detox <- has.detox
# Now do the same for each chrom separately
has.detox.2L <- unlist(lapply(CNV.encompass.list.2L, function(x) any(x %in% detox.genes)))
CNV.encompass.table.2L$has.detox <- has.detox.2L
has.detox.2R <- unlist(lapply(CNV.encompass.list.2R, function(x) any(x %in% detox.genes)))
CNV.encompass.table.2R$has.detox <- has.detox.2R
has.detox.3L <- unlist(lapply(CNV.encompass.list.3L, function(x) any(x %in% detox.genes)))
CNV.encompass.table.3L$has.detox <- has.detox.3L
has.detox.3R <- unlist(lapply(CNV.encompass.list.3R, function(x) any(x %in% detox.genes)))
CNV.encompass.table.3R$has.detox <- has.detox.3R
has.detox.X <- unlist(lapply(CNV.encompass.list.X, function(x) any(x %in% detox.genes)))
CNV.encompass.table.X$has.detox <- has.detox.X

# Sum cnv counts over all chromosomes:
pop.cnv.totals <- pop.cnv.totals.2L + 
                  pop.cnv.totals.2R + 
                  pop.cnv.totals.3L + 
                  pop.cnv.totals.3R + 
                  pop.cnv.totals.X 
pop.cnv.totals.table <- rbind(pop.cnv.totals.2L, pop.cnv.totals.2R, pop.cnv.totals.3L, pop.cnv.totals.3R, pop.cnv.totals.X, pop.cnv.totals)

# Now get a similar table for CNVs that contain at least one gene.
pop.cnv.withgene.2L <- pop.cnv.2L[rownames(CNV.encompass.table.2L),]
pop.cnv.withgene.totals.2L <- apply(pop.cnv.withgene.2L, 2, sum)
pop.cnv.withgene.2R <- pop.cnv.2R[rownames(CNV.encompass.table.2R),]
pop.cnv.withgene.totals.2R <- apply(pop.cnv.withgene.2R, 2, sum)
pop.cnv.withgene.3L <- pop.cnv.3L[rownames(CNV.encompass.table.3L),]
pop.cnv.withgene.totals.3L <- apply(pop.cnv.withgene.3L, 2, sum)
pop.cnv.withgene.3R <- pop.cnv.3R[rownames(CNV.encompass.table.3R),]
pop.cnv.withgene.totals.3R <- apply(pop.cnv.withgene.3R, 2, sum)
pop.cnv.withgene.X <- pop.cnv.X[rownames(CNV.encompass.table.X),]
pop.cnv.withgene.totals.X <- apply(pop.cnv.withgene.X, 2, sum)
pop.cnv.withgene.totals <- pop.cnv.withgene.totals.2L + 
                           pop.cnv.withgene.totals.2R + 
                           pop.cnv.withgene.totals.3L + 
                           pop.cnv.withgene.totals.3R + 
                           pop.cnv.withgene.totals.X 
pop.cnv.withgene.totals.table <- rbind(pop.cnv.withgene.totals.2L, pop.cnv.withgene.totals.2R, pop.cnv.withgene.totals.3L, pop.cnv.withgene.totals.3R, pop.cnv.withgene.totals.X, pop.cnv.withgene.totals)

# Now same again for detox genes
pop.cnv.detox.2L <- pop.cnv.2L[rownames(CNV.encompass.table.2L)[CNV.encompass.table.2L$has.detox], ,drop = F]
pop.cnv.detox.totals.2L <- apply(pop.cnv.detox.2L, 2, sum)
pop.cnv.detox.2R <- pop.cnv.2R[rownames(CNV.encompass.table.2R)[CNV.encompass.table.2R$has.detox], ,drop = F]
pop.cnv.detox.totals.2R <- apply(pop.cnv.detox.2R, 2, sum)
pop.cnv.detox.3L <- pop.cnv.3L[rownames(CNV.encompass.table.3L)[CNV.encompass.table.3L$has.detox], ,drop = F]
pop.cnv.detox.totals.3L <- apply(pop.cnv.detox.3L, 2, sum)
pop.cnv.detox.3R <- pop.cnv.3R[rownames(CNV.encompass.table.3R)[CNV.encompass.table.3R$has.detox], ,drop = F]
pop.cnv.detox.totals.3R <- apply(pop.cnv.detox.3R, 2, sum)
pop.cnv.detox.X <- pop.cnv.X[rownames(CNV.encompass.table.X)[CNV.encompass.table.X$has.detox], ,drop = F]
pop.cnv.detox.totals.X <- apply(pop.cnv.detox.X, 2, sum)
pop.cnv.detox.totals <- pop.cnv.detox.totals.2L + 
                        pop.cnv.detox.totals.2R + 
                        pop.cnv.detox.totals.3L + 
                        pop.cnv.detox.totals.3R + 
                        pop.cnv.detox.totals.X 
pop.cnv.detox.totals.table <- rbind(pop.cnv.detox.totals.2L, pop.cnv.detox.totals.2R, pop.cnv.detox.totals.3L, pop.cnv.detox.totals.3R, pop.cnv.detox.totals.X, pop.cnv.detox.totals)


###############################
# GENE CLUSTERING NORMAL FREQ #
###############################

# Let's now try to see if detox genes are enriched among the CNVs that we have found. 
# We could simulate all of the CNVs we found, and see how many of them contain genes / detox genes. Alternatively,
# we could look at the number of genes duplicated by each CNV and simulate that. Let's try the first option for 
# now. 

# Create a single list containing the useable positions on each chromosome
all.chrom.positions <- list('2L' = counts.list.positions.2L,
                            '2R' = counts.list.positions.2R,
                            '3L' = counts.list.positions.3L,
                            '3R' = counts.list.positions.3R,
                            'X' = counts.list.positions.X)
# Do the same with eligible gene positions 
eligible.gene.position.indices.allchrom <- list('2L' = eligible.gene.position.indices.2L,
                                                '2R' = eligible.gene.position.indices.2R,
                                                '3L' = eligible.gene.position.indices.3L,
                                                '3R' = eligible.gene.position.indices.3R,
                                                'X' = eligible.gene.position.indices.X)

# Write a function that will do the simulation
simulate.cnvs <- function(input.cnv.table, chrom.positions, gene.pos, detox.genes.list){
	# We will store our randomly generated CNVs in the following table
	random.cnvs <- data.frame(chrom = character(), start = numeric(), end = numeric(), Features = character(), stringsAsFactors = F)
	# Create an object listing the number of useable windows in each chromosome
	useable.windows <- unlist(lapply(chrom.positions, length))
	# For each CNV, give it a new location at random in the genome. 
	for (i in 1:nrow(input.cnv.table)){
		# First, make a weighted random choice of chromosome 
		chosen.chrom <- sample(names(useable.windows), 1, T, useable.windows)
		# Next, choose a start index, making sure it's not too large to accomodate the full size of the
		# duplication
		start.pos <- sample(1:(useable.windows[chosen.chrom] - input.cnv.table[i,'end'] + input.cnv.table[i,'start']), 1)
		# The end position is just the start position plus the size of the CNV
		end.pos <- start.pos + input.cnv.table[i,'end'] - input.cnv.table[i,'start']
		# Find the genes included in this interval
		included.genes <- find.features.encompassed(matrix(c(start.pos, end.pos),1,2), genepos.table = gene.pos[[chosen.chrom]], filtered.windows = chrom.positions[[chosen.chrom]], cov.win = coverage.window)
		if (nrow(included.genes))
			random.cnvs[i,] <- data.frame(chrom = chosen.chrom, start = start.pos, end = end.pos, Features = included.genes$Features, stringsAsFactors = F)[1,]
		else
			random.cnvs[i,] <- data.frame(chrom = chosen.chrom, start = start.pos, end = end.pos, Features = NA, stringsAsFactors = F)[1,]
	}
	# Find out which CNVs contain detox genes
	random.cnvs.list <- strsplit(sub('^Genes: ', '', random.cnvs$Features), ';')
	random.cnvs.has.detox <- unlist(lapply(random.cnvs.list, function(x) any(x %in% detox.genes.list)))
	random.cnvs$has.detox <- random.cnvs.has.detox
	random.cnvs
}

# Simulate K cnv sets and record the number of cnvs that have genes, and that have detox genes. 
simulation.one.output <- matrix(0,0,3, dimnames = list(c(), c('all.cnvs', 'gene.cnvs', 'detox.cnvs')))
K = 1000
cat('\nSimulation type 1:\n')
for (k in 1:K){
	this.simulation <- simulate.cnvs(good.CNV.ranges.allchrom, all.chrom.positions, eligible.gene.position.indices.allchrom, detox.genes)
	simulation.one.output <- rbind(simulation.one.output, c(nrow(this.simulation), sum(!is.na(this.simulation$Features)), sum(this.simulation$has.detox)))
	cat(k, ' ', sep = '')
}
cat('\n')
sim.above.gene.th <- simulation.one.output[,'gene.cnvs'] >= nrow(CNV.encompass.table.allchrom)
num.above.gene.th <- sum(sim.above.gene.th)
prop.above.gene.th <- 100 * num.above.gene.th / K
sim.above.detox.th <- (simulation.one.output[,'detox.cnvs'] / simulation.one.output[,'gene.cnvs']) > (sum(CNV.encompass.table.allchrom$has.detox) / nrow(CNV.encompass.table.allchrom))
num.above.detox.th <- sum(sim.above.detox.th)
prop.above.detox.th <- 100 * num.above.detox.th / K
cat('Out of ', K, ' simulations, ', num.above.gene.th, ' (', prop.above.gene.th, '%) had at least as many cnvs containing genes as we observed in the real data. ', sep = '')
cat('In ', num.above.detox.th, ' (', prop.above.detox.th, '%) of simulations, the proportion of gene-containing cnvs that contained detox genes was at least as high as observed in the real data.\n', sep = '')

# What if we simulate the cnvs by just shifting their position a bit (say by moving the start point to the end point
# or vice-versa). 
simulate.type.2 <- function(input.cnv.table, chrom.positions, gene.pos, detox.genes.list, direction = 'random'){
	# We will store our randomly generated CNVs in the following table
	random.cnvs <- data.frame(chrom = character(), start = numeric(), end = numeric(), Features = character(), stringsAsFactors = F)
	# Create an object listing the number of useable windows in each chromosome
	useable.windows <- unlist(lapply(chrom.positions, length))
	# For each CNV, give it a new location at random in the genome. 
	for (i in 1:nrow(input.cnv.table)){
		# The chromosome will be the same as the current one
		chosen.chrom <- sub('CNV_', '', rownames(input.cnv.table)[i])
		chosen.chrom <- sub('\\d+$', '', chosen.chrom)
		start.to.end <- NA
		if (direction == 'random'){
			# Now, we either move the start position to the end position or vice versa. If we are too close to the 
			# start or end for one of these options, we choose the other. If both are possible, we choose randomly. 
			if (input.cnv.table[i,'end'] > (2*input.cnv.table[i,'start']))
				start.to.end <- 1
			else if ((input.cnv.table[i,'end'] - input.cnv.table[i,'start']) > (useable.windows[[chosen.chrom]] - input.cnv.table[i,'end']))
				start.to.end <- -1
			else 
				start.to.end <- sample(c(1,-1), 1)
			# Set either the start of end point, then set the other accordingly
		}
		else if (direction == 'back'){
			# We move the end position to the start position or, if the start is too close to the start of the 
			# chromosome, we leave it where it is.
			if (input.cnv.table[i,'end'] > (2*input.cnv.table[i,'start']))
				start.to.end <- 0
			else 
				start.to.end <- -1
		}
		else if (direction == 'forward'){
			# We move the end position to the start position or, if the start is too close to the start of the 
			# chromosome, we leave it where it is.
			if ((input.cnv.table[i,'end'] - input.cnv.table[i,'start']) > (useable.windows[[chosen.chrom]] - input.cnv.table[i,'end']))
				start.to.end <- 0
			else 
				start.to.end <- 1
		}
		else {
			stop('"direction" should be one of "random", "back" or "forward".')
		}
		# Now set the new start and end positions based on the value of start.to.end
		if (start.to.end == 1){
			start.pos <- input.cnv.table[i,'end']
			end.pos <- start.pos + input.cnv.table[i,'end'] - input.cnv.table[i,'start']
		}
		else if (start.to.end == -1){
			end.pos <- input.cnv.table[i,'start']
			start.pos <- end.pos - input.cnv.table[i,'end'] + input.cnv.table[i,'start']
		}
		else if (start.to.end == 0){
			start.pos <- input.cnv.table[i,'start']
			end.pos <- input.cnv.table[i,'end']
		}
		# Find the genes included in this interval
		included.genes <- find.features.encompassed(matrix(c(start.pos, end.pos),1,2), genepos.table = gene.pos[[chosen.chrom]], filtered.windows = chrom.positions[[chosen.chrom]], cov.win = coverage.window)
		if (nrow(included.genes))
			random.cnvs[i,] <- data.frame(chrom = chosen.chrom, start = start.pos, end = end.pos, Features = included.genes$Features, stringsAsFactors = F)[1,]
		else
			random.cnvs[i,] <- data.frame(chrom = chosen.chrom, start = start.pos, end = end.pos, Features = NA, stringsAsFactors = F)[1,]
	}
	# Find out which CNVs contain detox genes
	random.cnvs.list <- strsplit(sub('^Genes: ', '', random.cnvs$Features), ';')
	random.cnvs.has.detox <- unlist(lapply(random.cnvs.list, function(x) any(x %in% detox.genes.list)))
	random.cnvs$has.detox <- random.cnvs.has.detox
	random.cnvs
}

# Simulate K.sim2 cnv sets and record the number of cnvs that have genes, and that have detox genes. 
simulation.two.output <- matrix(0,0,3, dimnames = list(c(), c('all.cnvs', 'gene.cnvs', 'detox.cnvs')))
K.sim2 = 100
cat('\nSimulation type 2:\n')
for (k in 1:K.sim2){
	this.simulation <- simulate.type.2(good.CNV.ranges.allchrom, all.chrom.positions, eligible.gene.position.indices.allchrom, detox.genes)
	simulation.two.output <- rbind(simulation.two.output, c(nrow(this.simulation), sum(!is.na(this.simulation$Features)), sum(this.simulation$has.detox)))
	cat(k, ' ', sep = '')
}
cat('\n')
sim.above.gene.th.2 <- simulation.two.output[,'gene.cnvs'] >= nrow(CNV.encompass.table.allchrom)
num.above.gene.th.2 <- sum(sim.above.gene.th.2)
prop.above.gene.th.2 <- 100 * num.above.gene.th.2 / K.sim2
sim.above.detox.th.2 <- (simulation.two.output[,'detox.cnvs'] / simulation.two.output[,'gene.cnvs']) > (sum(CNV.encompass.table.allchrom$has.detox) / nrow(CNV.encompass.table.allchrom))
num.above.detox.th.2 <- sum(sim.above.detox.th.2)
prop.above.detox.th.2 <- 100 * num.above.detox.th.2 / K.sim2
cat('Out of ', K.sim2, ' simulations, ', num.above.gene.th.2, ' (', prop.above.gene.th.2, '%) had at least as many cnvs containing genes as we observed in the real data. ', sep = '')
cat('In ', num.above.detox.th.2, ' (', prop.above.detox.th.2, '%) of simulations, the proportion of gene-containing cnvs that contained detox genes was at least as high as observed in the real data.\n', sep = '')

# The issue with the above simulations is that they are not independent (for a given CNV, there are only two possible 
# places it can move in a given simulation, so the simulations are only assortments of the same positions. So, instead,
# Let's just do two tests, one where we move everything back, and one where we move everything forward. CNVs that are
# too close to the border to move in a given direction are left where they are.
move.backwards <- simulate.type.2(good.CNV.ranges.allchrom, all.chrom.positions, eligible.gene.position.indices.allchrom, detox.genes, direction = 'back')
move.backwards.output <- c(nrow(move.backwards), sum(!is.na(move.backwards$Features)), sum(move.backwards$has.detox))
move.forwards <- simulate.type.2(good.CNV.ranges.allchrom, all.chrom.positions, eligible.gene.position.indices.allchrom, detox.genes, direction = 'forward')
move.forwards.output <- c(nrow(move.forwards), sum(!is.na(move.forwards$Features)), sum(move.forwards$has.detox))

# For each of the above, we can at least do a Fisher's exact test to find out if the new proportions of genes among 
# cnvs / detox genes among genes are significant. 
move.backwards.fisher.table.1 <- matrix(c(nrow(good.CNV.ranges.allchrom) - nrow(CNV.encompass.table.allchrom), 
                                          nrow(CNV.encompass.table.allchrom), 
                                          move.backwards.output[1] - move.backwards.output[2],
                                          move.backwards.output[2]), 2, 2, dimnames = list(c('no_gene', 'gene'), c('observed', 'simulated')))
cat('\nContingency table for CNVs containing or not containing genes, after moving all CNVs backwards:\n')
print(move.backwards.fisher.table.1)
print(fisher.test(move.backwards.fisher.table.1))
#
move.backwards.fisher.table.2 <- matrix(c(nrow(CNV.encompass.table.allchrom) - sum(CNV.encompass.table.allchrom$has.detox), 
                                          sum(CNV.encompass.table.allchrom$has.detox), 
                                          move.backwards.output[2] - move.backwards.output[3],
                                          move.backwards.output[3]), 2, 2, dimnames = list(c('non_detox', 'detox'), c('observed', 'simulated')))
cat('\nContingency table for CNVs containing genes including or not including detox genes, after moving all CNVs backwards:\n')
print(move.backwards.fisher.table.2)
print(fisher.test(move.backwards.fisher.table.2))
#
move.forwards.fisher.table.1 <- matrix(c(nrow(good.CNV.ranges.allchrom) - nrow(CNV.encompass.table.allchrom), 
                                          nrow(CNV.encompass.table.allchrom), 
                                          move.forwards.output[1] - move.forwards.output[2],
                                          move.forwards.output[2]), 2, 2, dimnames = list(c('no_gene', 'gene'), c('observed', 'simulated')))
cat('\nContingency table for CNVs containing or not containing genes, after moving all CNVs forwards:\n')
print(move.forwards.fisher.table.1)
print(fisher.test(move.forwards.fisher.table.1))
#
move.forwards.fisher.table.2 <- matrix(c(nrow(CNV.encompass.table.allchrom) - sum(CNV.encompass.table.allchrom$has.detox), 
                                          sum(CNV.encompass.table.allchrom$has.detox), 
                                          move.forwards.output[2] - move.forwards.output[3],
                                          move.forwards.output[3]), 2, 2, dimnames = list(c('non_detox', 'detox'), c('observed', 'simulated')))
cat('\nContingency table for CNVs containing genes including or not including detox genes, after moving all CNVs forwards:\n')
print(move.forwards.fisher.table.2)
print(fisher.test(move.forwards.fisher.table.2))

# Now do a third type of simulation where we simulate CNVs of duplicated genes and see whether they are enriched for
# detox genes. 
simulate.by.gene <- function(input.cnv.gene.list, gene.table){
	randomised.clusters <- detox.randomised.clusters <-  list()
	num.detox.genes <- 0
	# First, get a vector of the number of genes in each CNV
	gene.numbers <- unlist(lapply(input.cnv.gene.list, length))
	for (i in 1:length(gene.numbers)){
		# We create a while loop which we break out of once we have found a suitable cluster.
		while(1) {
			this.gene <- sample.int(nrow(gene.table), 1)
			# If the resulting cluster would exceed the size of the table, try again
			if ((this.gene + gene.numbers[i] - 1) > nrow(gene.table))
				next
			this.cluster <- gene.table[this.gene:(this.gene + gene.numbers[i] - 1), ]
			# Check that these genes aren't on different chromosomes
			if (this.cluster$Chrom[1] == this.cluster$Chrom[nrow(this.cluster)]){
				# As we go along, we check whether these clusters contain detox genes 
				if (sum(this.cluster$is.detox) > 0)
					detox.randomised.clusters <- c(detox.randomised.clusters, list(this.cluster))
				randomised.clusters[[i]] <- this.cluster
				# Also get the total number of detox genes in the clusters
				num.detox.genes <- num.detox.genes + sum(this.cluster$is.detox)
				break
			}
		}
	}
	list(all = randomised.clusters, detox = detox.randomised.clusters, all.num = length(randomised.clusters), detox.num = length(detox.randomised.clusters), detox.gene.num = num.detox.genes)
}

# Get a version of the all.genes object containing only the eligible genes
eligible.gene.names <- c(rownames(eligible.gene.position.indices.2L)[!is.na(eligible.gene.position.indices.2L[,'start'])],
                         rownames(eligible.gene.position.indices.2R)[!is.na(eligible.gene.position.indices.2R[,'start'])],
                         rownames(eligible.gene.position.indices.3L)[!is.na(eligible.gene.position.indices.3L[,'start'])],
                         rownames(eligible.gene.position.indices.3R)[!is.na(eligible.gene.position.indices.3R[,'start'])],
                         rownames(eligible.gene.position.indices.X)[!is.na(eligible.gene.position.indices.X[,'start'])])
all.eligible.genes <- all.genes[rownames(all.genes) %in% eligible.gene.names, ]

# Do that simulation K.simgene times 
simulation.gene.output <- matrix(0,0,3, dimnames = list(c(), c('all.clusters', 'detox.clusters', 'detox.genes')))
K.simgene = 1000
cat('\nSimulation by gene-containing CNVs:\n')
for (k in 1:K.simgene){
	this.simulation <- simulate.by.gene(CNV.encompass.list.allchrom, all.eligible.genes)
	simulation.gene.output <- rbind(simulation.gene.output, c(this.simulation$all.num, this.simulation$detox.num, this.simulation$detox.gene.num))
	cat(k, ' ', sep = '')
}
cat('\n')
sim.bygene.above.th <- simulation.gene.output[,'detox.clusters'] >= sum(CNV.encompass.table.allchrom$has.detox)
num.bygene.above.th <- sum(sim.bygene.above.th)
prop.bygene.above.th <- 100 * num.bygene.above.th / K.simgene
cat('Out of ', K.simgene, ' simulations, ', num.bygene.above.th, ' (', prop.bygene.above.th, '%) had at least as many detox clusters as we observed in the real data.\n', sep = '')

# We now want to create two tables showing the genes that are duplicated by at least one CNV that has achieved 
# reasonable frequency in that population. We create one table for all genes, and another for only detox genes.
empty.table <- all.genes[c(),]
empty.table$gene <- character()
empty.table$pop <- character()
empty.table$cnv.code <- character()
cnv.table.bypop <- cnv.table.bypop.detox <- empty.table
for (this.pop.name in pops){
	this.pop.table <- this.pop.detox.table <- empty.table
	# Get the CNVs that exist at high enough frequency in this population
	this.pop.cnvs <- rownames(pop.cnv.allchrom)[pop.cnv.allchrom[,this.pop.name]]
	# Find the ones of these that contain at least one gene
	this.pop.cnvs.withgenes <- intersect(this.pop.cnvs, names(CNV.encompass.list.allchrom))
	# Now go through each of these genes in turn and add the genes to the tables
	for (this.cnv in this.pop.cnvs.withgenes){
		genes.in.this.cnv <- CNV.encompass.list.allchrom[[this.cnv]]
		for (this.gene in genes.in.this.cnv){
			# Check if this gene has already been found in this population
			if (this.gene %in% rownames(this.pop.table)){
				# If so, just add this cnv.code to that row in the table
				this.pop.table[this.gene, 'cnv.code'] <- paste(this.pop.table[this.gene, 'cnv.code'], this.cnv, sep = ';') 
				# And, if the gene is detox, to the row in the detox table. 
				if (this.gene %in% detox.genes)
					this.pop.detox.table[this.gene, 'cnv.code'] <- paste(this.pop.detox.table[this.gene, 'cnv.code'], this.cnv, sep = ';') 
			}
			else{
				# Otherwise, add a row to the main table
				this.pop.table <- rbind(this.pop.table, cbind(all.genes[this.gene,], gene = this.gene, pop = this.pop.name, cnv.code = this.cnv, stringsAsFactors = F))
				# And, if the gene is detox, add a row to the detox table. 
				if (this.gene %in% detox.genes)
					this.pop.detox.table <- rbind(this.pop.detox.table, cbind(all.genes[this.gene,], gene = this.gene, pop = this.pop.name, cnv.code = this.cnv, stringsAsFactors = F))
			}
		}
	}
	# Add the resulting tables to their respective master tables, after adjusting the row names to include the 
	# population
	if (nrow(this.pop.table))
		rownames(this.pop.table) <- paste(rownames(this.pop.table), this.pop.name, sep = '_')
	if (nrow(this.pop.detox.table))
		rownames(this.pop.detox.table) <- paste(rownames(this.pop.detox.table), this.pop.name, sep = '_')
	cnv.table.bypop <- rbind(cnv.table.bypop, this.pop.table)
	cnv.table.bypop.detox <- rbind(cnv.table.bypop.detox, this.pop.detox.table)
}
# Add the annotations for all these genes
cnv.table.bypop$annotation <- annot.info[cnv.table.bypop$gene,]
cnv.table.bypop.detox$annotation <- annot.info[cnv.table.bypop.detox$gene,]

# Now output the tables and save the workspace.
write.table(cnv.table.bypop, file = 'CNV_table.csv', quote = F, sep = '\t', col.names = NA)
write.table(cnv.table.bypop.detox, file = 'CNV_detox_table.csv', quote = F, sep = '\t', col.names = NA)

# Remove the counts.list object to make the workspace smaller
rm(counts.list)

save.image('CNV_contingencies_alt.Rdata')

