# This script reduces the output from the CNV analysis to reduce its memory footprint 

load('/home/eric/Liverpool/CNV_v2_bigfiles/counting_output_v4_X/phase2/alt_fullgcnorm_nomapq_mapq002_varvar_trans000001/CNV_analysis_alt_X.Rdata')

cyp9k1.region.coord <- c(15220000, 15255000) 

	
for (sample.name in names(counts.list)){
	cat('\t', sample.name, '\n', sep = '')
	counts.list[[sample.name]] <- subset(counts.list[[sample.name]], (Position >= cyp9k1.region.coord[1]) & (Position <= cyp9k1.region.coord[2]))
}
	
# Now save this to an .Rdata file
save.image('CNV_analysis_phase2_for_CYP9K1.Rdata')

