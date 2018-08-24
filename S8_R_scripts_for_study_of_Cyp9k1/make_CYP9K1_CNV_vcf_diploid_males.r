# This script outputs a vcf file containing the information regarding the different CNVs in the CYP9K1 region. 

library(Biostrings)

load('CYP9K1_analysis_shrunk_data.Rdata')

# Get the reference genome
chromX <- genome[6]
# We have manually determined that the following duplications can be considered to be single duplications.
single.dups <- c('Dup1', 'Dup2', 'Dup3', 'Dup4', 'Dup5', 'Dup6', 'Dup8', 'Dup9', 'Dup12', 'Dup16')
# Rename a couple of objects for the for loop later. 
Dup3.bp <- Dup3.end.bp
Dup16.bp <- Dup16.end.bp

# The following object is a manually created object of the samples in the Ag1000G vcf files. 
vcf.header.samples <- c('AA0040-C','AA0041-C','AA0042-C','AA0043-C','AA0044-C','AA0048-C','AA0049-C','AA0050-C','AA0051-C','AA0052-C','AA0053-C','AA0054-C','AA0055-C','AA0056-C','AA0060-C','AA0061-C','AA0063-C','AA0064-C','AA0066-C','AA0067-C','AA0068-C','AA0072-C','AA0073-C','AA0074-C','AA0075-C','AA0076-C','AA0077-C','AA0080-C','AA0084-C','AA0085-C','AA0086-C','AA0087-C','AA0088-C','AA0089-C','AA0090-C','AA0091-C','AA0096-C','AA0097-C','AA0098-C','AA0099-C','AA0100-C','AA0101-C','AA0102-C','AA0103-C','AA0104-C','AA0107-C','AA0108-C','AA0109-C','AA0110-C','AA0111-C','AA0113-C','AA0114-C','AA0115-C','AA0116-C','AA0122-C','AA0123-C','AA0124-C','AA0125-C','AA0127-C','AA0132-C','AA0133-C','AA0134-C','AA0135-C','AA0136-C','AA0139-C','AA0140-C','AA0141-C','AB0085-C','AB0087-C','AB0088-C','AB0089-C','AB0090-C','AB0091-C','AB0092-C','AB0094-C','AB0095-C','AB0097-C','AB0098-C','AB0099-C','AB0100-C','AB0101-C','AB0103-C','AB0104-C','AB0108-C','AB0109-C','AB0110-C','AB0111-C','AB0112-C','AB0113-C','AB0114-C','AB0115-C','AB0117-C','AB0118-C','AB0119-C','AB0122-C','AB0123-C','AB0124-C','AB0126-C','AB0127-C','AB0128-C','AB0129-C','AB0130-C','AB0133-C','AB0134-C','AB0135-C','AB0136-C','AB0137-C','AB0138-C','AB0139-C','AB0140-C','AB0142-C','AB0143-C','AB0145-C','AB0146-C','AB0147-C','AB0148-C','AB0150-C','AB0151-C','AB0153-C','AB0155-C','AB0157-C','AB0158-C','AB0159-C','AB0160-C','AB0161-C','AB0162-C','AB0164-C','AB0165-C','AB0166-C','AB0167-C','AB0169-C','AB0170-C','AB0171-C','AB0172-C','AB0173-C','AB0174-C','AB0175-C','AB0176-C','AB0177-C','AB0178-C','AB0179-C','AB0181-C','AB0182-C','AB0183-C','AB0184-C','AB0185-C','AB0186-C','AB0187-C','AB0188-C','AB0189-C','AB0190-C','AB0191-C','AB0192-C','AB0195-C','AB0196-C','AB0197-C','AB0198-C','AB0199-C','AB0200-C','AB0201-C','AB0202-C','AB0203-C','AB0204-C','AB0205-C','AB0206-C','AB0207-C','AB0208-C','AB0209-C','AB0210-C','AB0211-C','AB0212-C','AB0213-C','AB0215-C','AB0217-C','AB0218-C','AB0219-C','AB0221-C','AB0222-C','AB0223-C','AB0224-C','AB0226-C','AB0227-C','AB0228-C','AB0229-C','AB0231-C','AB0232-C','AB0233-C','AB0234-C','AB0235-C','AB0236-C','AB0237-C','AB0238-C','AB0239-C','AB0240-C','AB0241-C','AB0242-C','AB0243-C','AB0244-C','AB0246-C','AB0247-C','AB0248-C','AB0249-C','AB0250-C','AB0251-C','AB0252-C','AB0253-C','AB0255-C','AB0256-C','AB0257-C','AB0258-C','AB0260-C','AB0261-C','AB0262-C','AB0263-C','AB0264-C','AB0265-C','AB0266-C','AB0267-C','AB0268-C','AB0270-C','AB0271-C','AB0272-C','AB0273-C','AB0274-C','AB0275-C','AB0276-C','AB0277-C','AB0278-C','AB0279-C','AB0280-C','AB0281-C','AB0282-C','AB0283-C','AB0284-C','AC0089-C','AC0090-C','AC0091-C','AC0092-C','AC0093-C','AC0094-C','AC0095-C','AC0096-C','AC0097-C','AC0098-C','AC0099-C','AC0100-C','AC0101-C','AC0102-C','AC0103-C','AC0104-C','AC0105-C','AC0106-C','AC0107-C','AC0108-C','AC0109-C','AC0110-C','AC0111-C','AC0112-C','AC0113-C','AC0114-C','AC0115-C','AC0116-C','AC0117-C','AC0118-C','AC0119-C','AC0120-C','AC0121-C','AC0122-C','AC0123-C','AC0124-C','AC0125-C','AC0126-C','AC0127-C','AC0128-C','AC0129-C','AC0130-C','AC0131-C','AC0132-C','AC0133-C','AC0134-C','AC0135-C','AC0136-C','AC0137-C','AC0138-C','AC0139-C','AC0140-C','AC0141-C','AC0142-C','AC0143-C','AC0144-C','AC0145-C','AC0146-C','AC0147-C','AC0148-C','AC0149-C','AC0150-C','AC0151-C','AC0152-C','AC0153-C','AC0154-C','AC0155-C','AC0156-C','AC0157-C','AC0158-C','AC0159-C','AC0160-C','AC0161-C','AC0162-C','AC0163-C','AC0164-C','AC0166-C','AC0167-C','AC0168-C','AC0169-C','AC0170-C','AC0171-C','AC0172-C','AC0173-C','AC0174-C','AC0176-C','AC0177-C','AC0178-C','AC0179-C','AC0180-C','AC0181-C','AC0182-C','AC0183-C','AC0184-C','AC0185-C','AC0186-C','AC0187-C','AC0188-C','AC0189-C','AC0190-C','AC0191-C','AC0192-C','AC0193-C','AC0194-C','AC0195-C','AC0196-C','AC0197-C','AC0199-C','AC0200-C','AC0201-C','AC0202-C','AC0203-C','AG0082-C','AG0085-C','AG0089-C','AG0096-C','AG0097-C','AG0098-C','AG0100-C','AG0102-C','AG0104-C','AG0106-C','AG0108-C','AG0109-C','AG0111-C','AG0118-C','AG0120-C','AG0121-C','AG0123-C','AG0125-C','AG0126-C','AG0127-C','AG0128-C','AG0129-C','AG0133-C','AG0134-C','AG0136-C','AG0137-C','AG0138-C','AG0139-C','AG0141-C','AG0142-C','AG0143-C','AG0144-C','AG0145-C','AG0146-C','AG0147-C','AG0148-C','AG0152-C','AG0153-C','AG0156-C','AG0159-C','AG0162-C','AG0163-C','AG0169-C','AG0170-C','AG0172-C','AG0178-C','AG0179-C','AG0181-C','AG0183-C','AG0195-C','AG0197-C','AG0202-C','AG0203-C','AG0204-C','AG0206-C','AG0208-C','AG0214-C','AG0221-C','AG0223-C','AG0227-C','AG0229-C','AG0230-C','AG0231-C','AG0232-C','AG0263-C','AJ0023-C','AJ0024-C','AJ0028-C','AJ0032-C','AJ0035-C','AJ0036-C','AJ0037-C','AJ0038-C','AJ0039-C','AJ0043-C','AJ0044-C','AJ0045-C','AJ0047-C','AJ0051-C','AJ0052-C','AJ0056-C','AJ0059-C','AJ0060-C','AJ0061-C','AJ0063-C','AJ0064-C','AJ0066-C','AJ0068-C','AJ0070-C','AJ0071-C','AJ0072-C','AJ0074-C','AJ0075-C','AJ0076-C','AJ0077-C','AJ0078-C','AJ0080-C','AJ0081-C','AJ0084-C','AJ0085-C','AJ0086-C','AJ0087-C','AJ0088-C','AJ0090-C','AJ0092-C','AJ0093-C','AJ0095-C','AJ0096-C','AJ0097-C','AJ0098-C','AJ0099-C','AJ0100-C','AJ0101-C','AJ0102-C','AJ0103-C','AJ0105-C','AJ0107-C','AJ0109-C','AJ0113-C','AJ0115-C','AJ0116-C','AJ0117-C','AJ0119-C','AJ0128-C','AJ0129-C','AJ0130-C','AJ0131-C','AJ0132-C','AJ0133-C','AJ0134-C','AJ0135-C','AJ0136-C','AJ0137-C','AJ0138-C','AJ0139-C','AJ0140-C','AJ0141-C','AJ0142-C','AJ0143-C','AJ0144-C','AJ0145-C','AJ0146-C','AJ0147-C','AJ0148-C','AJ0149-C','AJ0150-C','AJ0151-C','AJ0152-C','AJ0153-C','AJ0154-C','AJ0155-C','AJ0156-C','AJ0157-C','AJ0158-C','AJ0159-C','AJ0161-C','AK0060-C','AK0062-C','AK0065-C','AK0066-C','AK0067-C','AK0068-C','AK0069-C','AK0070-C','AK0072-C','AK0073-C','AK0074-C','AK0075-C','AK0076-C','AK0077-C','AK0078-C','AK0079-C','AK0080-C','AK0081-C','AK0082-C','AK0085-C','AK0086-C','AK0087-C','AK0088-C','AK0089-C','AK0090-C','AK0091-C','AK0092-C','AK0093-C','AK0094-C','AK0095-C','AK0096-C','AK0098-C','AK0099-C','AK0100-C','AK0101-C','AK0102-C','AK0103-C','AK0104-C','AK0105-C','AK0106-C','AK0107-C','AK0108-C','AK0109-C','AK0110-C','AK0116-C','AK0117-C','AK0119-C','AK0127-C','AN0007-C','AN0008-C','AN0009-C','AN0010-C','AN0011-C','AN0012-C','AN0013-C','AN0014-C','AN0015-C','AN0016-C','AN0017-C','AN0018-C','AN0019-C','AN0020-C','AN0022-C','AN0023-C','AN0024-C','AN0025-C','AN0026-C','AN0027-C','AN0028-C','AN0029-C','AN0030-C','AN0031-C','AN0032-C','AN0033-C','AN0034-C','AN0035-C','AN0036-C','AN0037-C','AN0038-C','AN0039-C','AN0040-C','AN0041-C','AN0042-C','AN0043-C','AN0044-C','AN0045-C','AN0046-C','AN0047-C','AN0048-C','AN0049-C','AN0050-C','AN0051-C','AN0053-C','AN0054-C','AN0055-C','AN0056-C','AN0057-C','AN0058-C','AN0059-C','AN0060-C','AN0061-C','AN0062-C','AN0063-C','AN0064-C','AN0065-C','AN0066-C','AN0067-C','AN0068-C','AN0069-C','AN0070-C','AN0071-C','AN0072-C','AN0073-C','AN0074-C','AN0075-C','AN0076-C','AN0077-C','AN0079-C','AN0080-C','AN0081-C','AN0082-C','AN0083-C','AN0084-C','AN0085-C','AN0086-C','AN0087-C','AN0088-C','AN0089-C','AN0090-C','AN0091-C','AN0092-C','AN0093-C','AN0094-C','AN0095-C','AN0096-C','AN0097-C','AN0098-C','AN0099-C','AN0100-C','AN0101-C','AN0102-C','AN0103-C','AN0104-C','AN0105-C','AN0106-C','AN0107-C','AN0108-C','AN0109-C','AN0111-C','AN0112-C','AN0113-C','AN0114-C','AN0115-C','AN0117-C','AN0118-C','AN0119-C','AN0120-C','AN0121-C','AN0122-C','AN0123-C','AN0124-C','AN0125-C','AN0126-C','AN0127-C','AN0128-C','AN0129-C','AN0130-C','AN0131-C','AN0132-C','AN0134-C','AN0135-C','AN0136-C','AN0137-C','AN0138-C','AN0139-C','AN0140-C','AN0141-C','AN0142-C','AN0143-C','AN0144-C','AN0147-C','AN0149-C','AN0151-C','AN0152-C','AN0153-C','AN0154-C','AN0155-C','AN0156-C','AN0157-C','AN0158-C','AN0159-C','AN0160-C','AN0162-C','AN0163-C','AN0164-C','AN0165-C','AN0166-C','AN0167-C','AN0168-C','AN0169-C','AN0170-C','AN0171-C','AN0172-C','AN0173-C','AN0174-C','AN0175-C','AN0176-C','AN0177-C','AN0178-C','AN0179-C','AN0180-C','AN0181-C','AN0182-C','AN0183-C','AN0184-C','AN0185-C','AN0186-C','AN0187-C','AN0188-C','AN0189-C','AN0190-C','AN0191-C','AN0192-C','AN0193-C','AN0194-C','AN0196-C','AN0197-C','AN0198-C','AN0199-C','AN0200-C','AN0201-C','AN0202-C','AN0203-C','AN0204-C','AN0205-C','AN0206-C','AN0207-C','AN0208-C','AN0209-C','AN0210-C','AN0212-C','AN0213-C','AN0214-C','AN0215-C','AN0217-C','AN0218-C','AN0219-C','AN0220-C','AN0221-C','AN0222-C','AN0223-C','AN0224-C','AN0225-C','AN0226-C','AN0227-C','AN0228-C','AN0229-C','AN0230-C','AN0231-C','AN0232-C','AN0233-C','AN0234-C','AN0235-C','AN0236-C','AN0237-C','AN0238-C','AN0239-C','AN0240-C','AN0241-C','AN0242-C','AN0243-C','AN0244-C','AN0245-C','AN0246-C','AN0247-C','AN0248-C','AN0250-C','AN0251-C','AN0252-C','AN0253-C','AN0254-C','AN0255-C','AN0256-C','AN0258-C','AN0259-C','AN0260-C','AN0261-C','AN0262-C','AN0263-C','AN0264-C','AN0265-C','AN0266-C','AN0267-C','AN0268-C','AN0269-C','AN0270-C','AN0271-C','AN0272-C','AN0275-C','AN0276-C','AN0277-C','AN0278-C','AN0279-C','AN0280-C','AN0281-C','AN0282-C','AN0283-C','AN0284-C','AN0285-C','AN0286-C','AN0287-C','AN0288-C','AN0289-C','AN0290-C','AN0291-C','AN0292-C','AN0293-C','AN0294-C','AN0295-C','AN0296-C','AN0297-C','AN0298-C','AN0299-C','AN0300-C','AN0301-C','AN0302-C','AN0303-C','AN0304-C','AN0305-C','AN0306-C','AN0307-C','AN0308-C','AN0309-C','AN0310-C','AN0311-C','AN0312-C','AN0313-C','AN0314-C','AN0315-C','AN0316-C','AN0317-C','AN0318-C','AN0319-C','AN0320-C','AN0321-C','AP0002-C','AP0005-C','AP0006-C','AP0007-C','AP0008-C','AP0009-C','AP0010-C','AP0011-C','AP0014-C','AP0017-C','AP0018-C','AP0019-C','AP0020-C','AP0021-C','AP0022-C','AP0023-C','AP0024-C','AP0025-C','AP0030-C','AP0031-C','AP0032-C','AP0033-C','AP0034-C','AP0035-C','AQ0001-C','AQ0002-C','AQ0004-C','AQ0005-C','AQ0011-C','AQ0012-C','AQ0013-C','AQ0014-C','AQ0015-C','AR0001-C','AR0007-C','AR0008-C','AR0009-C','AR0010-C','AR0011-C','AR0012-C','AR0013-C','AR0014-C','AR0015-C','AR0016-C','AR0017-C','AR0018-C','AR0019-C','AR0020-C','AR0021-C','AR0022-C','AR0023-C','AR0024-C','AR0026-C','AR0027-C','AR0034-C','AR0035-C','AR0036-C','AR0038-C','AR0040-C','AR0042-C','AR0043-C','AR0044-C','AR0045-C','AR0046-C','AR0047-C','AR0048-C','AR0049-C','AR0050-C','AR0051-C','AR0052-C','AR0053-C','AR0054-C','AR0057-C','AR0059-C','AR0060-C','AR0061-C','AR0062-C','AR0063-C','AR0064-C','AR0065-C','AR0066-C','AR0069-C','AR0070-C','AR0071-C','AR0072-C','AR0073-C','AR0074-C','AR0075-C','AR0076-C','AR0077-C','AR0078-C','AR0079-C','AR0080-C','AR0081-C','AR0082-C','AR0083-C','AR0084-C','AR0085-C','AR0086-C','AR0087-C','AR0088-C','AR0089-C','AR0090-C','AR0092-C','AR0093-C','AR0095-C','AR0096-C','AR0097-C','AR0098-C','AR0099-C','AR0100-C','AS0001-C','AS0002-C','AS0003-C','AS0004-C','AS0005-C','AS0006-C','AS0007-C','AS0008-C','AS0009-C','AS0010-C','AS0011-C','AS0012-C','AS0013-C','AS0014-C','AS0015-C','AS0016-C','AS0017-C','AS0018-C','AS0019-C','AS0020-C','AS0021-C','AS0022-C','AS0024-C','AS0025-C','AS0026-C','AS0027-C','AS0028-C','AS0029-C','AS0030-C','AS0032-C','AS0033-C','AS0034-C','AS0035-C','AS0036-C','AS0037-C','AS0039-C','AS0040-C','AS0041-C','AS0042-C','AS0044-C','AS0045-C','AS0046-C','AS0047-C','AS0048-C','AS0049-C','AS0051-C','AS0052-C','AS0053-C','AS0054-C','AS0055-C','AS0056-C','AS0057-C','AS0058-C','AS0059-C','AS0060-C','AS0062-C','AS0064-C','AS0065-C','AS0066-C','AS0068-C','AS0069-C','AS0070-C','AS0071-C','AS0072-C','AS0073-C','AS0074-C','AS0076-C','AS0077-C','AS0078-C','AV0001-C','AV0002-C','AV0003-C','AV0004-C','AV0005-C','AV0006-C','AV0007-C','AV0008-C','AV0009-C','AV0010-C','AV0011-C','AV0012-C','AV0013-C','AV0014-C','AV0015-C','AV0017-C','AV0018-C','AV0020-C','AV0021-C','AV0022-C','AV0024-C','AV0025-C','AV0026-C','AV0027-C','AV0028-C','AV0029-C','AV0030-C','AV0031-C','AV0032-C','AV0033-C','AV0034-C','AV0035-C','AV0036-C','AV0037-C','AV0038-C','AV0039-C','AV0040-C','AV0041-C','AV0042-C','AV0043-C','AV0044-C','AV0045-C','AV0046-C','AV0047-C','AY0006-C','AY0007-C','AY0010-C','AY0011-C','AY0012-C','AY0013-C','AY0015-C','AY0016-C','AY0017-C','AY0018-C','AY0019-C','AY0020-C','AY0021-C','AY0023-C','AY0024-C','AY0025-C','AY0026-C','AY0027-C','AY0029-C','AY0031-C','AY0032-C','AY0033-C','AY0034-C','AY0035-C','AY0036-C','AY0038-C','AY0039-C','AY0040-C','AY0041-C','AY0042-C','AY0043-C','AY0045-C','AY0046-C','AY0047-C','AY0048-C','AY0049-C','AY0050-C','AY0052-C','AY0053-C','AY0054-C','AY0055-C','AY0056-C','AY0057-C','AY0058-C','AY0059-C','AY0060-C','AY0061-C','AY0062-C','AY0063-C','AY0064-C','AY0065-C','AY0066-C','AY0067-C','AY0068-C','AY0069-C','AY0070-C','AY0072-C','AY0074-C','AY0076-C','AY0077-C','AY0078-C','AY0079-C','AY0080-C','AY0082-C','AY0083-C','AY0085-C','AY0087-C','AY0088-C','AY0089-C','AY0090-C','AY0091-C')
# Sanity check. Make sure the order of samples in the genotype table is the same as in the Ag1000G vcf headers
if (any(vcf.header.samples != sub('_', '-', rownames(genotype.calls))))
	stop('Order of samples in genotype call table does not match Ag1000G vcf files.')

vcf.header.samples.males <- vcf.header.samples[meta$sex == 'M']
vcf.header.samples.females <- vcf.header.samples[meta$sex == 'F']
genotype.calls.males <- genotype.calls[meta$sex == 'M', ]
genotype.calls.females <- genotype.calls[meta$sex == 'F', ]

# Make the first few lines of the vcf
vcf.header <- '##fileformat=VCFv4.0
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=.,Type=Integer,Description="Phred-scaled Genotype Likelihoods">'
# We first put all of the females, then all of the males, since the shapeit file only contains the females, and 
# thus will not have the samples in the normal order
vcf.column.heads <- paste('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT', paste(vcf.header.samples.females, collapse = '\t'), paste(vcf.header.samples.males, collapse = '\t'), sep = '\t')
vcf.txt <- paste(vcf.header, vcf.column.heads, sep = '\n')

fem.empty.allele.table <- matrix('0', length(vcf.header.samples.females), 2, dimnames = list(rownames(genotype.calls.females), c('gen1', 'gen2')))
male.empty.allele.table <- matrix('0', length(vcf.header.samples.males), 2, dimnames = list(rownames(genotype.calls.males), c('gen1', 'gen2')))
call.list <- list()
for (this.dup in single.dups){
	# Get the position of this duplication (taken as its start point).
	this.dup.pos <- get(paste(this.dup, 'bp', sep = '.'))[1]
	call.list[[this.dup]] <- list()
	call.list[[this.dup]]$females <- fem.empty.allele.table
	call.list[[this.dup]]$males <- male.empty.allele.table
	# Where a sample is classed as carrying a dup, assume it has at least one copy of the CNV. For the males,
	# This means coding it as a homozygous diploid
	this.dup.samples <- get(paste(this.dup, 'samples', sep = '.'))
	this.dup.fem.samples <- intersect(this.dup.samples, rownames(genotype.calls.females))
	this.dup.male.samples <- intersect(this.dup.samples, rownames(genotype.calls.males))
	call.list[[this.dup]]$females[this.dup.fem.samples, 'gen2'] <- '1'
	call.list[[this.dup]]$males[this.dup.male.samples, ] <- '1'
	# Where the genotype calling has claimed that coverage was higher than 3 in females, call a homozygote
	fem.hom.mut <- rownames(genotype.calls.females)[genotype.calls.females[,this.dup] > 3]
	fem.hom.mut <- fem.hom.mut[!is.na(fem.hom.mut)]
	call.list[[this.dup]]$females[fem.hom.mut, 'gen1'] <- '1'
	# Where the coverage call is NA, record the allele as .
	call.list[[this.dup]]$females[is.na(genotype.calls.females[,this.dup])] <- '.'
	call.list[[this.dup]]$males[is.na(genotype.calls.males[,this.dup])] <- '.'
	# Sanity check. No female row should have column1 = 1 and column2 = 0
	if (any(call.list[[this.dup]]$females[,'gen1'] == '1' & call.list[[this.dup]]$females[,'gen2'] == '0'))
		stop('Fail. The first genotype column cannot be 1 if the second is 0.')
	# Turn the genotypes into character vectors
	this.dup.calls.females <- apply(call.list[[this.dup]]$females, 1, paste, collapse = '/')
	this.dup.calls.males <- apply(call.list[[this.dup]]$males, 1, paste, collapse = '/')
	# For the MVNcall software, we need to assign phred-scaled likelihoods to each of the three possible genotypes
	# (wt/wt, wt/mut, mut/mut) for each sample. We simply assign 0 for the genotype we believe is true, and 50 for
	# the other genotypes.
	for (i in 1:length(this.dup.calls.females)){
		if (this.dup.calls.females[i] == '0/0')
			this.dup.calls.females[i] <- paste(this.dup.calls.females[i], ':0,50,50', sep = '')
		else if (this.dup.calls.females[i] == '0/1')
			this.dup.calls.females[i] <- paste(this.dup.calls.females[i], ':50,0,50', sep = '')
		else if (this.dup.calls.females[i] == '1/1')
			this.dup.calls.females[i] <- paste(this.dup.calls.females[i], ':50,50,0', sep = '')
		else if (this.dup.calls.females[i] == './.')
			this.dup.calls.females[i] <- paste(this.dup.calls.females[i], ':.', sep = '')
	}
	# And again for males, but this time there is no heterzygote option
	for (i in 1:length(this.dup.calls.males)){
		if (this.dup.calls.males[i] == '0/0')
			this.dup.calls.males[i] <- paste(this.dup.calls.males[i], ':0,50,50', sep = '')
		else if (this.dup.calls.males[i] == '1/1')
			this.dup.calls.males[i] <- paste(this.dup.calls.males[i], ':50,50,0', sep = '')
		else if (this.dup.calls.males[i] == './.')
			this.dup.calls.males[i] <- paste(this.dup.calls.males[i], ':.', sep = '')
	}
	# Get the reference base at this position
	ref <- substr(chromX, this.dup.pos, this.dup.pos)
	# Create a vcf line for this CNV. For the MVNcall software, we need to give values to REF and ALT. For REF, 
	# we take the value from the reference genome. For Alt, we just put that same character twice. 
	# For the chromosome name, the haplotype scaffold file codes X as 0, so we use this here as well to match.
	vcf.line <- paste('0', this.dup.pos, this.dup, ref, paste(rep(ref,2), collapse = ''), '.', '.', '.', 'GT:PL', paste(this.dup.calls.females, collapse = '\t'), paste(this.dup.calls.males, collapse = '\t'), sep = '\t')
	# Add this to the vcf.txt object
	vcf.txt <- c(vcf.txt, vcf.line)
}

# Now write the vcf file to text
write(vcf.txt, file = 'CYP9K1_CNV_diploid_males.vcf', sep = '\n')


