# SEACR
## SEACR: *S*parse *E*nrichment *A*nalysis for *C*UT&*R*UN

SEACR is intended to call peaks and enriched regions from sparse CUT&RUN or chromatin profiling data in which background is dominated by "zeroes" (i.e. regions with no read coverage). It requires R (https://www.r-project.org) and Bedtools (https://bedtools.readthedocs.io/en/latest/) to be available in your path, and it requires bedgraphs from paired-end sequencing as input, which can be generated from *read pair* BED files (i.e. BED coordinates reflecting the 5' and 3' termini of each read pair) using bedtools genomecov with the "-bg" flag, or alternatively from name-sorted paired-end BAM files as described in "Preparing input bedgraph files" below. 

A description of the method can be found in the following manuscript, which we respectfully request that you cite if you find SEACR useful in your research:

Meers MP, Tenenbaum D, Henikoff S. (2019). Peak calling by Sparse Enrichment Analysis for CUT&RUN chromatin profiling. *Epigenetics and Chromatin* 12(1):42. 

Direct link: https://doi.org/10.1186/s13072-019-0287-4

## SEACR web server

A web interface for SEACR analysis can be found at https://seacr.fredhutch.org

## Recent changes

### v1.4

- Added dashed command line parameter support: inputs 1, 2, 3, 4, and 5 now correspond to -b, -c, -n, -m, and -o, respectively
- Added a function to test a window of control scaling thresholds in order to render scaling more stable in "norm" mode
- Added an optional input for a numeric peak extension factor (-e) during peak merging, where all peaks will be extended by e*(mean(peak length)) (default = 0.1)
- Added an optional yes/no input to dictate whether peaks overlapping robust IgG peaks get filtered out or not (-r) (default = "yes")
- Fixed a bug in which empty peak files failed due to a division by zero error before completion
- Fixed a bug in which "NaN" is erroneously passed to quantile() in files with a low complexity of peaks

### v1.3

- Fixed a bug in which the bedgraph line thresholding added in v1.2 was failing for some datasets.
- Added a check to filter out any input bedgraph lines containing zero signal.

### v1.2

- Fixed a bug in lines 166 and 168 in which misplaced brackets caused the misreporting of the max signal region terminal coordinate for merged signal blocks.
- Added a counter to keep track of the number of component bedgraph lines that compose each signal block, and a function to calculate the minimum threshold of lines per signal block at which there is a smaller percentage of target signal blocks remaining than control. This is meant to be used as a filter for signal blocks that pass the total signal threshold despite being composed of very few bedgraph lines, which are unlikely to be true peaks.
- Changed how the dataframe for density plotting is truncated (previously a hard-coded 90% cutoff): a dataframe of list quantile (i.e. line #/max line#) vs. value quantile (i.e. value/max value) is derived, and the threshold is selected by finding the dataframe pair for which the orthogonal distance below the line defined by (0,0);(1,1) is maximized.

### v1.1
- Changed "union" and "AUC" modes to "relaxed" and "stringent" modes, respectively.
- Removed maximum signal threshold from "relaxed" mode and replaced it with an alternate total signal threshold that uses the point halfway between the knee and the peak of the total signal curve as described in the manuscript text. This change improves performance at high read depth.
- Implemented alternate threshold test that searches for any thresholds that come within 95% of the optimal threshold. This change avoids spurious thresholds that are overselective in some datasets.

## Usage: 

	bash SEACR_1.3.sh experimental bedgraph [control bedgraph | numeric threshold] ["norm" | "non"] ["relaxed" | "stringent"] output prefix

## Description of input fields:

Field 1: Target data bedgraph file in UCSC bedgraph format (https://genome.ucsc.edu/goldenpath/help/bedgraph.html) that omits regions containing 0 signal.

Field 2: Control (IgG) data bedgraph file to generate an empirical threshold for peak calling. Alternatively, a numeric threshold *n* between 0 and 1 returns the top *n* fraction of peaks based on total signal within peaks. 

Field 3: “norm” denotes normalization of control to target data, “non” skips this behavior. "norm" is recommended unless experimental and control data are already rigorously normalized to each other (e.g. via spike-in).

Field 4: “relaxed” uses a total signal threshold between the knee and peak of the total signal curve, and corresponds to the “relaxed” mode described in the text, whereas “stringent” uses the peak of the curve, and corresponds to “stringent” mode.

Field 5: Output prefix

## Preparing input bedgraph files

Bedgraph files should reflect density across *read pairs* rather than individual reads. If starting from BAM files, we recommend converting to paired end BED files using bedtools bamtobed with the -bedpe flag, then selecting the 5' and 3' coordinates of the read pair to generate a new BED3 file, and finally converting that file to a bedgraph using bedtools genomecov.

Here is some example code for converting from a paired-end BAM to a fragment bedgraph file as described above:

	bedtools bamtobed -bedpe -i $sample.bam > $sample.bed
	awk '$1==$4 && $6-$2 < 1000 {print $0}' $sample.bed > $sample.clean.bed
	cut -f 1,2,6 $sample.clean.bed | sort -k1,1 -k2,2n -k3,3n > $sample.fragments.bed
	bedtools genomecov -bg -i $sample.fragments.bed -g my.genome > $sample.fragments.bedgraph

## Output file:

	<output prefix>.stringent.bed OR <output prefix>.relaxed.bed (BED file of enriched regions)
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

	bash SEACR_1.4.sh -b target.bedgraph -c IgG.bedgraph -n norm -m stringent -o output
Calls enriched regions in target data using normalized IgG control track with stringent threshold and default peak extension factor of 0.1
	
	bash SEACR_1.4.sh -b target.bedgraph -c IgG.bedgraph -n non -m relaxed -o output -e 0.5
Calls enriched regions in target data using non-normalized IgG control track with relaxed threshold and peak extension factor of 0.5

	bash SEACR_1.4.sh -b target.bedgraph -c 0.01 -n non -m stringent -o output
Calls enriched regions in target data by selecting the top 1% of regions by total signal
