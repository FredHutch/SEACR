# SEACR
## SEACR: *S*parse *E*nrichment *A*nalysis for *C*UT&*R*UN

SEACR is intended to call peaks and enriched regions from sparse CUT&RUN or chromatin profiling data in which background is dominated by "zeroes" (i.e. regions with no read coverage). It requires R (https://www.r-project.org) and Bedtools (https://bedtools.readthedocs.io/en/latest/) to be available in your path, and it requires bedgraphs from paired-end sequencing as input, which can be generated from *read pair* BED files (i.e. BED coordinates reflecting the 5' and 3' termini of each read pair) using bedtools genomecov with the "-bg" flag, or alternatively from name-sorted paired-end BAM files as described in "Preparing input bedgraph files" below. 

A description of the method can be found in the following manuscript:

Meers MP, Bryson TD, Henikoff S. A streamlined protocol and analysis pipeline for CUT&RUN chromatin profiling. bioRxiv doi: https://doi.org/10.1101/569129

Direct link: https://www.biorxiv.org/content/10.1101/569129v2

## Recent changes

### v1.1
- Changed "union" and "AUC" modes to "relaxed" and "stringent" modes, respectively.
- Removed maximum signal threshold from "relaxed" mode and replaced it with an alternate total signal threshold that uses the point halfway between the knee and the peak of the total signal curve as described in the manuscript text. This change improves performance at high read depth.
- Implemented alternate threshold test that searches for any thresholds that come within 95% of the optimal threshold. This change avoids spurious thresholds that are overselective in some datasets.

## Usage: 

	bash SEACR_1.1.sh experimental bedgraph [control bedgraph | numeric threshold] ["norm" | "non"] ["relaxed" | "stringent"] output prefix

## Description of input fields:

Field 1: Target data bedgraph file in UCSC bedgraph format (https://genome.ucsc.edu/goldenpath/help/bedgraph.html) that omits regions containing 0 signal.

Field 2: Control (IgG) data bedgraph file to generate an empirical threshold for peak calling. Alternatively, a numeric threshold *n* between 0 and 1 returns the top *n* fraction of peaks based on total signal within peaks. 

Field 3: “norm” denotes normalization of control to target data, “non” skips this behavior. "norm" is recommended unless experimental and control data are already rigorously normalized to each other (e.g. via spike-in).

Field 4: “relaxed” uses a total signal threshold between the knee and peak of the total signal curve, and corresponds to the “relaxed” mode described in the text, whereas “stringent” uses the peak of the curve, and corresponds to “stringent” mode.

Field 5: Output prefix

## Preparing input bedgraph files

Bedgraph files should reflect density across *read pairs* rather than individual reads. If starting from BAM files, we recommend converting to paired end BED files using bedtools bamtobed with the -bedpe flag, then selecting the 5' and 3' coordinates of the read pair to generate a new BED3 file, and finally converting that file to a bedgraph using bedtools genomecov.

## Output file:

	<output prefix>.auc.threshold.merge.bed (BED file of enriched regions)
## Output data structure: 
	
	<chr>	<start>	<end>	<total signal>	<max signal>	<max signal region>

## Description of output fields:

Field 1: Chromosome

Field 2: Start coordinate

Field 3: End coordinate

Field 4: Total signal contained within denoted coordinates

Field 5: Maximum bedgraph signal attained at any base pair within denoted coordinates

Field 6: Region representing the farthest upstream and farthest downstream bases within the denoted coordinates that are represented by the maximum bedgraph signal

## Examples:

	bash SEACR_1.1.sh target.bedgraph IgG.bedgraph norm stringent output
Calls enriched regions in target data using normalized IgG control track with stringent threshold
	
	bash SEACR_1.1.sh target.bedgraph IgG.bedgraph non relaxed output
Calls enriched regions in target data using non-normalized IgG control track with relaxed threshold

	bash SEACR_1.1.sh target.bedgraph 0.01 non stringent output
Calls enriched regions in target data by selecting the top 1% of regions by AUC
