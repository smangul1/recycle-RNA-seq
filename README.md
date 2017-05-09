# recycle-RNA-seq-reads
Scripts and commands we used in our study : Comprehensive analysis of RNA-sequencing to find the source of 1 trillion reads across diverse adult human tissues


# VJ recombinations across GTEx samples

Inferred VJ recombinations are freely available at 
- [https://github.com/smangul1/recycle-RNA-seq/blob/master/IGK_VJ_recomb.csv](https://github.com/smangul1/recycle-RNA-seq/blob/master/IGK_VJ_recomb.csv)
- [https://github.com/smangul1/recycle-RNA-seq/blob/master/IGL_VJ_recomb.csv](https://github.com/smangul1/recycle-RNA-seq/blob/master/IGL_VJ_recomb.csv)

# Uncategorized reads  

Uncategorized reads (i.e. RNA-seq reads not categoried by ROP) from SRA samples are freely available. Curently we have deposited uncategorized reads for 71 SRA RNA-Seq samples. 

- [https://github.com/smangul1/recycle-RNA-seq/tree/master/unReadsSRA](https://github.com/smangul1/recycle-RNA-seq/tree/master/unReadsSRA)


# Workflow to categorize the mapped reads
Map reads onto human genome and transcriptome  
We mapped reads onto the human transcriptome (Ensembl GRCh37) and genome reference (Ensembl hg19) using tophat2 (v 2.0.13) with the default parameters. Tophat2 was supplied with a set of known transcripts (as a GTF formatted file, Ensembl GRCh37) using –G option.  The mapped reads of each sample are stored in a binary format (.bam).  

## Categorize mapped reads into genomic categories
ROP categorizes the reads into genomic categories based on the compatibility of each read from the pair with the features defined by Ensembl (GRCh37) gene annotations. First, we determined CDS, UTR3, UTR5 coordinates. We downloaded annotations for CDS, UTR3, UTR5 from UCSC Genome Browser (http://genome.ucsc.edu/cgi-bin/hgTables) in BED (browser extensible data) format. Next, we used gene annotations (a GTF formatted file, Ensembl GRCh37) to determine intron coordinates and inter-genic regions. We defined two types of inter-genic regions: ‘(proximate) inter-genic’ region (1Kb from the gene boundaries) and ‘deep inter-genic’ (beyond a proximity of 1Kb from the gene boundaries). 

Next, we checked the compatibility of the mapped reads with the defined genomic features, as follows:  

a.	Read mapped to multiple locations on the reference genome is categorized as a multi-mapped read.
b.	Read fully contained within the CDS, intron, UTR3, or UTR5 boundaries of a least one transcript is classified as a CDS, intronic, UTR3, or UTR5, respectively.
c.	Read simultaneously overlapping UTR3 and UTR5 regions is classified as a UTR read.
d.	Read spanning exon-exon boundary is defined as a junction read.
e.	Read mapped outside of gene boundaries and within a proximity of 1Kb is defined as a (proximal) inter-genic read.
f.	Read mapped outside of gene boundaries and beyond the proximity of 1Kb is defined as a deep inter-genic read.
g.	Read mapped to mitochondrial DNA (MT tag in hg19) is classified as a mitochondrial read.
h.	 Reads from a pair mapped to different chromosomes are classified as a fusion read.
