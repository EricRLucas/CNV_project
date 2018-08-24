# Function for getting the runs of duplicated sequences
find.full.dup.sequences <- function(cnv, n, threshold){
	# If cnv is a dataframe, extract the cnv vector
	if (is.data.frame(cnv))
		cnv.vector <- cnv$CNV
	else
		cnv.vector <- cnv
	current.run <- 0
	output <- matrix(0,0,2)
	for (i in 1:length(cnv.vector)){
		if (cnv.vector[i] >= threshold) {
			if (current.run == 0)
				start.pos <- i
			current.run <- current.run + 1
		}
		else{
			if (current.run >= n)
				output <- rbind(output, c(start.pos, i-1))
			current.run <- 0
		}
	}
	if (current.run >= n)
		output <- rbind(output, c(start.pos, i-1))
	return(output)
}

# This function takes two sets of regions (as matrices where the first column represents start positions
# and the second column represents corresponding end positions) and determines for each region defined in 
# the first matrix whether it overlaps with any regions in the second. The two matrices should be sorted 
# by their first column.
find.features.overlaps <- function(pos.table, genepos.table, exonpos.table){
	output.table <- data.frame(matrix(0,0,3, dimnames = list(character(),c('start', 'end', 'Features'))))
	if (nrow(pos.table) == 0){
		return(output.table)
	}
	for (i in 1:nrow(pos.table)){
		# Find the regions where the end position in the genepos.matrix is larger than the start position of the 
		# current region
		gene.req1 <- genepos.table[,'end'] > pos.table[i,1]
		# Find the regions where the start position in the genepos.matrix is smaller than the end position of the 
		# current region
		gene.req2 <- genepos.table[,'start'] < pos.table[i,2]
		# The regions where both of these requirements are met represent an overlap
		gene.overlaps <- genepos.table[gene.req1 & gene.req2,]
		# Do the same again for exons
		exon.req1 <- exonpos.table[,'end'] > pos.table[i,1]
		exon.req2 <- exonpos.table[,'start'] < pos.table[i,2]
		exon.overlaps <- exonpos.table[exon.req1 & exon.req2,]
		# If there is at least one overlap, add the pos.table region to the output table and record the name of the 
		# genepos.table regions that overlapped it
		if (nrow(gene.overlaps)){
			region.names <- paste('Genes:', paste(rownames(gene.overlaps), collapse = ";"))
			output.table <- rbind(output.table, data.frame(start = pos.table[i,1], end = pos.table[i,2], Features = region.names))
		}
		# And again for exons
		if (nrow(exon.overlaps)){
			region.names <- paste('Exons:', paste(rownames(exon.overlaps), collapse = ";"))
			output.table <- rbind(output.table, data.frame(start = pos.table[i,1], end = pos.table[i,2], Features = region.names))
		}
	}
	return(output.table)
}


# This function takes two sets of regions (as matrices where the first column represents start positions
# and the second column represents corresponding end positions) and determines for each region defined in 
# the first matrix whether it encompasses any region in the second. The two matrices should be sorted by
# their first column.
find.features.encompassed <- function(pos.table, genepos.table, filtered.windows, cov.win){
	output.table <- data.frame(matrix(0,0,3, dimnames = list(character(),c('start', 'end', 'Features'))))
	if (nrow(pos.table) == 0){
		return(output.table)
	}
	for (i in 1:nrow(pos.table)){
		if (!is.null(rownames(pos.table)))
			i <- rownames(pos.table)[i]
		# Find the genes whose post-filtering windows are all included in this duplication
		gene.req <- genepos.table[,'start'] >= pos.table[i,1] & genepos.table[,'end'] <= pos.table[i,2]
		# Gene that do not have enough covered windows to be considered will produce NAs in the above line. Turn 
		# these to FALSE
		gene.req[is.na(gene.req)] <- F
		# The regions where both of these requirements are met represent a gene encompassed in the duplication
		genes.encompassed <- genepos.table[gene.req,]
		# If there is at least one overlap, add the pos.table region to the output table and record the name of the 
		# genepos.table regions that overlapped it
		if (nrow(genes.encompassed)){
			region.names <- paste('Genes:', paste(rownames(genes.encompassed), collapse = ";"))
			output.table[i,] <- data.frame(start = pos.table[i,1], end = pos.table[i,2], Features = region.names, stringsAsFactors = F)[1,]
		}
	}
	# Now convert the indices to genomic positions in the output table
	output.table[,1] <- filtered.windows[output.table[,1]]
	output.table[,2] <- filtered.windows[output.table[,2]] + cov.win
	return(output.table)
}





