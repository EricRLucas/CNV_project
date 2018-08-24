# Load a library that will allow padding 0s in front of a number
library(stringr)

# First, we load up some functions to help us
source('/home/eric/Liverpool/liv_scripts/mHMM_analysis_functions_v3.r')

# This is the minimum number of consecutive duplicated windows required to call a duplication
min.dup.length <- 5
# This is the minimum copy number state required to consider a duplication (2 = "normal" diploid state).
threshold.copy.number <- 3
haploid.threshold.copy.number <- 2
# This is the coverage window that was used to calculate the coverage for the samples we are loading. 
# This parameter should not be changed unless using input data calculated with a different window size.
coverage.window <- 300 
# Set the minimum proportion of windows in a gene that need to have some coverage data before making a duplication call.
min.win.cov <- 0.5
# Set the threshold number of samples we need in a given population
th.num <- 3
# Set the threshold proportion of samples we need in a given population
th <- 0.05
high.freq.th <- 0.5
# Set the likelihood ratio threshold at which we will accept a gene
lik.thresh <- 10000

# The chromsome can only be X
chrom <- 'X'

# Load the metadata
meta <- read.table('/home/eric/Liverpool/phase2.AR1/samples/samples.meta.txt', header = T, row.names = 1, sep = '\t', quote = "", comment.char = "")
rownames(meta) <- sub('-', '_', rownames(meta))
pops <- levels(meta$population)
pop.sizes <- table(meta$population)

# Load the coverage variance data
cov.var <- read.table('/home/eric/Liverpool/CNV_v2_bigfiles/coverage_variance_masked_09.csv', header = T, sep = '\t', row.names = 1)
# Identify samples with coverage variance > 0.2 for exclusion
high.var.samples <- rownames(cov.var)[cov.var$autosomes > 0.2]
# Get the meta data and population sizes after exclusion of high variance samples
meta.reduced <- meta[!(rownames(meta) %in% high.var.samples), ]
pop.sizes.reduced <- table(meta.reduced$population)
# Get the threshold number of samples required for each pop.
pop.th <- pop.sizes.reduced * th
pop.high.freq.th <- pop.sizes.reduced * high.freq.th

# Write a function to calculate the likelihood ratio of a duplication vs no duplication in a given region in
# a given sample
lik.ratio <- function(sample.name, region.index.range, variance.table = cov.var, base.variance = 0.01, trans.prob = 0.00001, ploidy = 2){
	# Narrow down the counts table for this sample to this region
	these.counts <- counts.list[[sample.name]][region.index.range, ]
	# Get the variance around autosomal coverage normalised to 2 in this sample
	this.var <- variance.table[sample.name, 'autosomes']
	# Get the probability of not changing state and the ratio of the two probabilities
	non.trans.prob <- 1 - 12*trans.prob
	# Get the number of transitions in the cnv states
	if (length(region.index.range) > 1)
		transitions <- sum(these.counts$CNV[2:length(region.index.range)] - these.counts$CNV[1:(length(region.index.range)-1)] != 0)
	else
		transitions <- 0
	# get the ratio of transition probabilities
	trans.ratio <- (trans.prob^transitions)/(non.trans.prob^transitions)
	# get the likelihood of the null model. 
	lik.null <- dnorm(these.counts$Normalised_coverage, ploidy, sqrt(base.variance + (this.var/2)*ploidy))
	# get the likelihood of the duplication model. 
	lik.dup <- dnorm(these.counts$Normalised_coverage, these.counts$CNV, sqrt(base.variance + (this.var/2)*these.counts$CNV))
	# If the likelihood of the duplicated model at any given window is 0, then coverage must have been so high
	# that even CNV state 12 could not capture it, so we report a very large likelihood ratio
	if (any(lik.dup == 0))
		return(Inf)
	# Otherwise, get the product of the ratios multiplied by the transition probability
	else
		return(trans.ratio * prod(lik.dup/lik.null))
}

# Write a function that will go through all of the duplications in a list and merge them if they have 
# similar start and end positions. Let's write it such that the function can be run on either indices or
# positions, making the permissiveness around the start and end points an argument, so that appropriate 
# values can be given depending on whether indices or positions are used. 
merge.duplications <- function(dup.list, leeway){
	# First, make a copy of the input list
	dup.list.out <- dup.list
	# Turn the matrices that make up this list into dataframes
	for (i in 1:length(dup.list.out)){
		dup.list.out[[i]] <- as.data.frame(dup.list.out[[i]])
		colnames(dup.list.out[[i]]) <- c('start', 'end')
		dup.list.out[[i]][,'CNV_code'] <- 0
	}
	# Here is the code we will assign the next CNV group we find
	current.CNV.code = 1
	# Now go through every dup, looking for matching ranges
	for (i in 1:length(dup.list.out)){
		cat('Analysing sample ', names(dup.list.out)[i], '.\n', sep = '')
		for (k in 1:nrow(dup.list.out[[i]])){
			this.start <- dup.list.out[[i]][k,1]
			this.end <- dup.list.out[[i]][k,2]
			# Go through every previous sample, finding matching duplications. We record matching duplications 
			# by the sample and row in which they are found, and their current CNV_code.
			matching.dups <- matrix(0,0,3)
			for (j in (1:i)[-i]){
				# find any dups that match the start and end points of this one, including the leeway
				match.start <- (dup.list.out[[j]][,1] >= this.start - leeway) & (dup.list.out[[j]][,1] <= this.start + leeway)
				match.end <- (dup.list.out[[j]][,2] >= this.end - leeway) & (dup.list.out[[j]][,2] <= this.end + leeway)
				matching.dup.index <- which(match.start & match.end)
				# We should never find more than one matching dup in a given sample. Check that this is true
				if (length(matching.dup.index) > 1) 
					stop(paste('Fail. Found more than one matching duplications for i = ', i, ', k = ', k, ', j = ', j, '.'))
				else if (length(matching.dup.index) == 1)
					matching.dups <- rbind(matching.dups, c(j, matching.dup.index, dup.list.out[[j]][matching.dup.index, 3]))
			}
			# Now we have identified all of the dups that match this one, we need to merge them all together
			# If there are no matching dups, then give the current CNV the next available code
			if (nrow(matching.dups) == 0){
				dup.list.out[[i]][k,3] <- current.CNV.code
				current.CNV.code <- current.CNV.code + 1
			}
			else {
				# If there is only one matching dup, just give the same code to the current dup
				matching.dup.codes <- sort(unique(matching.dups[,3]))
				if (length(matching.dup.codes) == 1)
					dup.list.out[[i]][k,3] <- matching.dups[1,3]
				# Otherwise, find the smallest of the matching codes and give all matching dups that code. We
				# also need to find all the dups that have the matching code (not just the ones that directly
				# matched this dup) and change the code for those dups. 
				else {
					smallest.code <- matching.dup.codes[1]
					cat('\tMerging CNV codes ', paste(matching.dup.codes, collapse = ', '), '.\n', sep = '')
					# Give this code to the current dup
					dup.list.out[[i]][k,3] <- smallest.code
					# For each other dup code, find all the samples that have that code and change it to the 
					# smallest code
					for (j in (1:i)[-i]){
						dup.list.out[[j]][dup.list.out[[j]][,3] %in% matching.dup.codes[-1], 3] <- smallest.code
					}
				}
			}
		}
	}
	# As a final check, we need to make sure that no sample has the same duplication code twice. 
	for (d in names(dup.list.out))
		if (length(unique(dup.list.out[[d]][,3])) < nrow(dup.list.out[[d]]))
			stop(paste('Fail. The same duplication appears more than once in sample ', d, '.', sep = ''))
	dup.list.out
}

# Load up the information on gene and exon positions. 
cat('Loading gene.positions table\n')
gene.positions <- read.table('/home/eric/Liverpool/CNV_/find_gene_regions_output.csv')
gene.positions <- gene.positions[gene.positions$Chrom == chrom, ]

# Identify all the files containing mHMM output. 
all.files <- list.files()
counts.files <- all.files[grepl('^counts_for_mHMM_python2_.*\\.csv$', all.files)]
# Load up each file and store the results in a list
cat('Loading counts files\n')
counts.list <- list()
for (this.file in counts.files){
	# Extract the sample name information
	sample.name <- regmatches(this.file, regexpr('(?<=_)A.\\d{4}_C.?(?=_)', this.file, perl=T))
	# Load the file and add the table to the counts.list
	cat('\t', this.file, '\n', sep = '')
	counts.list[[sample.name]] <- read.table(this.file, header = T, row.names = 1)
}
counts.list.positions <- counts.list[[1]]$Position

# For each of those tables, we find the duplications of a minimum length and store their positions.
cat('Finding regions of duplication at least', min.dup.length*coverage.window, 'bases long.\n')
duplications.list.indices <- list()
for (this.sample in names(counts.list)){
	this.thresh <- ifelse(meta[this.sample, 'sex'] == 'M', haploid.threshold.copy.number, threshold.copy.number)
	duplications.list.indices[[this.sample]] <- find.full.dup.sequences(counts.list[[this.sample]], n = min.dup.length, threshold = this.thresh)
}

# That gives us the index for the first and last position in the duplication, but we want to turn this into 
# positions on the genome. We do this by indexing the positions vector for this chromosome. We increase the 
# value of the end point by the size of the window, since we want the end position of the last window where
# we see an increase in coverage, rather than the start position of that window. 
cat('Calculating genomic positions for these regions\n')
duplications.list <- duplications.list.indices
for (i in 1:length(duplications.list.indices)){
	duplications.list[[i]][,1] <- counts.list.positions[duplications.list[[i]][,1]]
	duplications.list[[i]][,2] <- counts.list.positions[duplications.list[[i]][,2]] + coverage.window
}

# Now remove samples with high variance
goodvar.duplications.list.indices <- duplications.list.indices[!(names(duplications.list.indices) %in% high.var.samples)]

# Now go through all of these duplications and keep only the ones with a good likelihood ratio.
cat('Filtering raw CNV calls by likelihood ratio.\n')
goodvar.goodlr.duplications.list.indices <- list()
for (this.sample in names(goodvar.duplications.list.indices)){
	filtered.table <- matrix(0,0,2)
	this.cnv.table <- goodvar.duplications.list.indices[[this.sample]]
	for (k in 1:nrow(this.cnv.table)){
		this.raw.cnv <- this.cnv.table[k,1]:this.cnv.table[k,2]
		if (meta[this.sample, 'sex'] == 'M')
			this.ploidy <- 1
		else 
			this.ploidy <- 2
		this.lik.ratio <- lik.ratio(this.sample, this.raw.cnv, ploidy = this.ploidy)
		if (this.lik.ratio > lik.thresh)
			filtered.table <- rbind(filtered.table, this.cnv.table[k,])
	}
	goodvar.goodlr.duplications.list.indices[[this.sample]] <- filtered.table
}

cat('Merging CNV clusters.\n')
clustered.list <- merge.duplications(goodvar.goodlr.duplications.list.indices, 1)

# Order these data by CNV cluster, instead of by sample
combined.clustered.list <- do.call(rbind, clustered.list)
all.CNVs <- unique(combined.clustered.list[,3])
duplications.by.cluster <- list()
for (this.CNV in all.CNVs){
	this.CNV.table <- combined.clustered.list[combined.clustered.list[,3] == this.CNV,]
	rownames(this.CNV.table) <- sub('\\..+', '', rownames(this.CNV.table))
	duplications.by.cluster[[paste('CNV_', chrom, str_pad(this.CNV, 4, pad = '0'), sep = '')]] <- this.CNV.table[,1:2]
}

# Now, within a population, find CNVs that exist at a minimum freq 
# Create a dataframe with set columns and no rows. 
pop.cnv.counts <- data.frame()
for (pop in pops)
	pop.cnv.counts[[pop]] <- numeric()
# Now loop through the duplicated genes and add them to the table
for (cnv in names(duplications.by.cluster)){
	# Get the samples with this duplication
	duplicated.samples <- rownames(duplications.by.cluster[[cnv]])
	# Count the number of individuals with the duplication in each population
	pop.cnv.counts[cnv,] <- table(meta.reduced[duplicated.samples, 'population'])
}

# Turn that table into logical values based on whether they have met the required threshold frequencies. 
pop.cnv <- (pop.cnv.counts > matrix(pop.th, nrow(pop.cnv.counts), ncol(pop.cnv.counts), byrow = T)) & 
           (pop.cnv.counts > th.num)

# Let's get a list of the CNVs that are at the minimum frequency in at least one population
good.CNVs <- rownames(pop.cnv)[apply(pop.cnv, 1, any)]

# Create a table showing the start and end points of these CNVs. For each start point, we take the median of 
# the start points of all the samples that have that CNV. We then turn the indices into genomic positions.
good.CNV.ranges <- data.frame(start = numeric(), end = numeric())
for (this.cnv in good.CNVs){
	good.CNV.ranges[this.cnv, ] <- c(ceiling(median(duplications.by.cluster[[this.cnv]][, 'start'])), floor(median(duplications.by.cluster[[this.cnv]][, 'end'])))
}

# For each duplications table, identify any genes that are encompassed by any duplications, with a minimum
# number of windows providing useable coverage data. 
# First, convert the positions of the gene table into indices of windows that passed filtering. 
eligible.gene.position.indices <- gene.positions
for (this.gene in rownames(eligible.gene.position.indices)){
	this.gene.start <- eligible.gene.position.indices[this.gene, 'start']
	this.gene.end <- eligible.gene.position.indices[this.gene, 'end']
	# Get the windows included in this gene
	this.gene.covered.windows <- which((counts.list.positions >= this.gene.start) & (counts.list.positions <= this.gene.end - coverage.window))
	# Get the windows that would be included in this gene if none had been filtered
	this.gene.full.windows <- floor(this.gene.end/coverage.window) - ceiling(this.gene.start/coverage.window)
	# If there are no potential windows fully encompassed in the gene, declare the gene ineligible
	if (this.gene.full.windows <= 0){
		eligible.gene.position.indices[this.gene, 'start'] <- NA
		eligible.gene.position.indices[this.gene, 'end'] <- NA
	}
	else if (length(this.gene.covered.windows) >= (min.win.cov * this.gene.full.windows)){
		eligible.gene.position.indices[this.gene, 'start'] <- min(this.gene.covered.windows)
		eligible.gene.position.indices[this.gene, 'end'] <- max(this.gene.covered.windows)
	}
	else {
		eligible.gene.position.indices[this.gene, 'start'] <- NA
		eligible.gene.position.indices[this.gene, 'end'] <- NA
	}
}

# Now look for encompassed genes
cat('Finding genes encompassed within CNVs\n')
# Note that the eligible.gene.position.indices table records the indices of the windows that are fully contained
# inside each gene. Thus, implicitly, we are only requiring that these windows are in a cnv in order to call a
# gene as duplicated. If the windows either side (overlapping the ends of the gene) are not in the cnv, this will
# not prevent the duplication call. 
CNV.encompass.table <- find.features.encompassed(good.CNV.ranges, genepos.table = eligible.gene.position.indices, filtered.windows = counts.list.positions, cov.win = coverage.window)
# Create a list that captures the genes present in each CNV in a more effective way.
CNV.encompass.list <- strsplit(sub('^Genes: ', '', CNV.encompass.table$Features), ';')
names(CNV.encompass.list) <- rownames(CNV.encompass.table)

# Now we use the gene-centric approach to get the genes duplicated in each sample, irrespective of CNV code. 
cat('\nGene-centric approach:\n')
cat('Finding genes encompassed within duplications\n')
encompass.list <- lapply(duplications.list.indices, find.features.encompassed, genepos.table = eligible.gene.position.indices, filtered.windows = counts.list.positions, cov.win = coverage.window)

# From all of this, we can create a single table, where all the elements of encompass.list are rbinded, with an extra
# first column that records the sample in which the duplication was found. 
cat('Preparing encompass.output.by.sample table\n')
encompass.output.by.sample <- data.frame(Sample = factor(), start = numeric(), end = numeric(), Features = factor())
for (i in 1:length(encompass.list)){
	this.sample.name <- names(encompass.list)[i]
	this.table <- encompass.list[[i]]
	this.table <- cbind(data.frame(Sample = rep(this.sample.name, nrow(this.table))), this.table)
	encompass.output.by.sample <- rbind(encompass.output.by.sample, this.table)
}

# Alternatively, for each gene that we found that overlapped with a duplication, we can list the samples in which 
# that was the case.
cat('Preparing output.by.gene lists\n')
# We get a list of all the different genes covered by our duplications
encompass.gene.column <- as.character(encompass.output.by.sample$Features[grepl('^Genes:', encompass.output.by.sample$Features)])
# Split up all the genes within those lists (some entries are several features separated by semi-colons).
all.encompass.genes <- unique(unlist(strsplit(sub('^Genes:', '', encompass.gene.column), ';')))
all.encompass.genes <- sub(' ', '', all.encompass.genes)
# Create a list with an entry for each gene
encompass.output.by.gene <- list()
encompass.output.by.gene[all.encompass.genes] <- list(character())
# Go through each entry in the output.by.sample table and add that sample to the genes it is associated with.
for (i in 1:nrow(encompass.output.by.sample)){
	cat('Making encompass list for', this.sample.name, '\n')
	this.sample.name <- as.character(encompass.output.by.sample[i,'Sample'])
	these.features <- encompass.output.by.sample[i,'Features']
	if (grepl('^Genes:', these.features)){
		these.genes <- unlist(strsplit(sub('^Genes:', '', these.features), ';'))
		these.genes <- sub(' ', '', these.genes)
		for (this.gene in these.genes){
			encompass.output.by.gene[[this.gene]] <- c(encompass.output.by.gene[[this.gene]], this.sample.name)
		}
	}
}
# Finally, we go through those lists and remove any duplicates (Note: I don't think duplicates are possible
# given the method used above, but leaving this code here in case I missed something). 
for (i in 1:length(encompass.output.by.gene)){
	encompass.output.by.gene[[i]] <- unique(encompass.output.by.gene[[i]])
}

# Here are all the sample names
sample.names <- names(counts.list)

save.image(paste('CNV_analysis_alt_', chrom, '.Rdata', sep = ''))

