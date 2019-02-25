# SEACR
## SEACR: Sparse Enrichment Analysis for CUT&RUN

SEACR is intended to call peaks and enriched regions from sparse CUT&RUN or chromatin profiling data in which background is dominated by "zeroes" (i.e. regions with no read coverage). It requires R (https://www.r-project.org) and Bedtools (https://bedtools.readthedocs.io/en/latest/) to be available in your path, and it requires bedgraphs as input, which can be generated from fragment BAM or BED files using bedtools genomecov with the "-bg" flag. 

## Usage: 

	bash SEACR_1.0.sh experimental bedgraph [control bedgraph | numeric threshold] ["norm" | "non"] ["union" | "AUC"]
	
## Description of input fields:

Field 1: Target data bedgraph file in UCSC bedgraph format (https://genome.ucsc.edu/goldenpath/help/bedgraph.html) that omits regions containing 0 signal.

Field 2: Control (IgG) data bedgraph file to generate an empirical threshold for peak calling. Alternatively, a numeric threshold n between 0 and 1 returns the top n fraction of peaks based on total signal within peaks. 

Field 3: “norm” denotes normalization of control to target data, “non” skips this behavior. "norm" is recommended unless experimental and control data are already rigorously normalized to each other (e.g. via spike-in).

Field 4: “union” forces implementation of a maximum signal threshold in addition to the total signal threshold, and corresponds to the “union” mode described in the text, whereas “AUC” avoids this behavior, and corresponds to “AUC only” mode.

## Output file:

	<experimental bedgraph>.auc.threshold.merge.bed (BED file of enriched regions)
## Output data structure: 
	
	<chr>	<start>	<end>	<total signal>	<max signal>	<max signal region>

## Description of output fields:

Field 1: Chromosome

Field 2: Start coordinate

Field 3: End coordinate

Field 4: Total signal contained within denoted coordinates

Field 5: Maximum begraph signal attained at any base pair within denoted coordinates

Field 6: Region representing the farthest upstream and farthest downstream bases within the denoted coordinates that are represented by the maximum bedgraph signal

## Examples:

	bash SEACR_1.0.sh target.bedgraph IgG.bedgraph norm AUC
Calls enriched regions in target data using normalized IgG control track with AUC threshold
	
	bash SEACR_1.0.sh target.bedgraph IgG.bedgraph non union
Calls enriched regions in target data using non-normalized IgG control track with AUC and max signal thresholds 

	bash SEACR_1.0.sh target.bedgraph 0.01 non AUC
Calls enriched regions in target data by selecting the top 1% of regions by AUC
