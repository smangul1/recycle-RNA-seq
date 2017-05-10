# recycle-RNA-seq-reads
Scripts and commands we used in our study : Comprehensive analysis of RNA-sequencing to find the source of 1 trillion reads across diverse adult human tissues


# VJ recombinations across GTEx samples

Inferred VJ recombinations are freely available at 
- [https://github.com/smangul1/recycle-RNA-seq/blob/master/IGK_VJ_recomb.csv](https://github.com/smangul1/recycle-RNA-seq/blob/master/IGK_VJ_recomb.csv)
- [https://github.com/smangul1/recycle-RNA-seq/blob/master/IGL_VJ_recomb.csv](https://github.com/smangul1/recycle-RNA-seq/blob/master/IGL_VJ_recomb.csv)

# Bacterial taxa assigned with Metaphlan2 

339 bacterial taxa were assigned with Metaphlan2  available at 
- https://github.com/smangul1/recycle-RNA-seq/blob/master/bacteriaTaxa.txt


# Uncategorized reads  

Uncategorized reads (i.e. RNA-seq reads not categoried by ROP) from SRA samples are freely available. Curently we have deposited uncategorized reads for 71 SRA RNA-Seq samples. 

- [https://github.com/smangul1/recycle-RNA-seq/tree/master/unReadsSRA](https://github.com/smangul1/recycle-RNA-seq/tree/master/unReadsSRA)


# Workflow to categorize the mapped reads

## Map reads onto human genome and transcriptome  
We mapped reads onto the human transcriptome (Ensembl GRCh37) and genome reference (Ensembl hg19) using tophat2 (v 2.0.13) with the default parameters. Tophat2 was supplied with a set of known transcripts (as a GTF formatted file, Ensembl GRCh37) using –G option.  The mapped reads of each sample are stored in a binary format (.bam).  

## Categorize mapped reads into genomic categories
ROP categorizes the reads into genomic categories based on the compatibility of each read from the pair with the features defined by Ensembl (GRCh37) gene annotations. First, we determined CDS, UTR3, UTR5 coordinates. We downloaded annotations for CDS, UTR3, UTR5 from UCSC Genome Browser (http://genome.ucsc.edu/cgi-bin/hgTables) in BED (browser extensible data) format. Next, we used gene annotations (a GTF formatted file, Ensembl GRCh37) to determine intron coordinates and inter-genic regions. We defined two types of inter-genic regions: ‘(proximate) inter-genic’ region (1Kb from the gene boundaries) and ‘deep inter-genic’ (beyond a proximity of 1Kb from the gene boundaries). 

Next, we checked the compatibility of the mapped reads with the defined genomic features, as follows:  

- Read mapped to multiple locations on the reference genome is categorized as a multi-mapped read.
- Read fully contained within the CDS, intron, UTR3, or UTR5 boundaries of a least one transcript is classified as a CDS, intronic, UTR3, or UTR5, respectively.
- c.	Read simultaneously overlapping UTR3 and UTR5 regions is classified as a UTR read.
- Read spanning exon-exon boundary is defined as a junction read.
- Read mapped outside of gene boundaries and within a proximity of 1Kb is defined as a (proximal) inter-genic read.
- Read mapped outside of gene boundaries and beyond the proximity of 1Kb is defined as a deep inter-genic read.
- Read mapped to mitochondrial DNA (MT tag in hg19) is classified as a mitochondrial read.
- Reads from a pair mapped to different chromosomes are classified as a fusion read.

More details are available at https://github.com/smangul1/rop/wiki/ROP-output-details


## Categorize mapped reads overlapping repeat instances 

Mapped reads were categorized based on the overlap with the repeat instances defined by RepeatMasker annotation (Repeatmasker v3.3, Repeat Library 20120124). RepeatMasker masks the repeats using the RepBase library: (http://www.girinst.org/repbase/update/index.html), which contains prototypic sequences representing repetitive DNA from different eukaryotic species. We use GTF files generated from the RepeatMasker annotations by Jin, Ying, et al. (Jin et al., 2015) and downloaded from: 
http://labshare.cshl.edu/shares/mhammelllab/www-data/TEToolkit/TE_GTF/hg19_rmsk_TE.gtf.gz 

Following  Melé, Marta, et al. (2015), repeat elements overlapping CDS regions are excluded from the analysis. We filtered out 6,873 repeat elements overlapping CDS regions. Prepared repeat annotations (bed formatted file) are available at https://drive.google.com/file/d/0Bx1fyWeQo3cORi1UNWhxOW9kYUk/view?pref=2&pli=1 


The prepared repeat annotations contain 8 Classes and 43 Families.  Number of elements per family and class are available  from the table below:

| classID |  N |
| --- | ---| 
| DNA | 458223 |
| LINE | 1478382 |
| LTR | 707384 | 
| RC | | 2226 | 
| SVA | 3582 |
| RNA | 717 | 
| Satellite | 8950 | 
| SINE | 1765403 | 


## Categorize mapped reads overlapping B cell receptor (BCR) and T cell receptor (TCR) loci

We used the gene annotations (Ensembl GRCh37) to extract BCR and TCR genes. We extracted gene annotations of the ‘constant’ (labeled as IG_C_gene, Ensembl GRCh37), ‘variable’ (labeled as IG_V_gene, Ensembl GRCh37), ‘diversity’ (labeled as IG_D_gene, Ensembl GRCh37), and ‘joining’ genes (labeled as IG_J_gene, Ensembl GRCh37) of BCR and TCR loci.  We excluded the BCR and TCR pseudogenes (labeled as IG_C_pseudogene, IG_V_pseudogene, IG_D_pseudogene, IG_J_pseudogene, TR_C_pseudogene, TR_V_pseudogene, TR_D_pseudogene, and TR_J_pseudogene). In addition, we excluded the patch contigs HG1592_PATCH and HG7_PATCH, as they are not part of the Ensembl hg19 reference, and reads are not mapped on the patch contigs by high throughput aligners.  After following the filtering steps described above, we extracted a total of 386 immune genes: 207 BCR genes and 179 TCR genes.  The gene annotations for antibody genes (GTF formatted file) are available at https://drive.google.com/file/d/0Bx1fyWeQo3cObFZNT3kyQlZUS1E/view?pref=2&pli=1

The number of VDJ genes per locus is reported in the Table  bellow: 


| | C domain|V domain|D domain| J  domain|
| ---| --- | --- | ---|  --- | 
| IGH locus| 8| 55 | 38| 6|
| IGK locus| 1| 46| -| 5|
|IGL locus| 4| 37| -| 7|
|TCRA locus| |1| 46|-| 57|
|TCRB locus| 1| 39| 0| 8|
|TRG locus| 2| 9| -| 5|
|TRD locus| 1| 3| 11| 4|

The number of reads mapping to each C-V-D-J genes was obtained by counting the number of sequencing reads that align, with high confidence, to each of the genes (HTSeq v0.6.1) (Anders et al., 2014). Script “htseq-count” is supplied with the gene annotations for BCR and TCR genes (genes_Ensembl_GRCh37_BCR_TCR.gtf) and a bam file. The bam file contains reads mapped to the human genome and transcriptome using TopHat2 (See Section “Map reads onto human genome and transcriptome” for details). The script generates individual gene counts by examining the read compatibility with BCR and TCR genes. We chose a conservative setting (--mode=intersection-strict) to handle reads overlapping more than one feature. Thus, a read overlapping several genes simultaneously is marked as a read with no feature and is excluded from the consideration. 

# Workflow for categorizing the unmapped reads

We first converted the unmapped reads saved by tophat2 from a BAM file into a FASTQ file (using bamtools). The FASTQ file of unmapped contain full read pairs (both ends of a read pair were unmapped) and discordant read pairs (one read end was mapped while the other end was unmapped). We disregarded the pairing information of the unmapped reads and categorize unmapped reads using the following steps:

## A. Quality Control
Low quality reads, defined as reads that have quality lower than 30 in at least 75% of their base pairs, were identified by [FASTX](http://hannonlab.cshl.edu/fastx_toolkit/) (v 0.0.13).  Low complexity reads, defined as reads with sequences of consecutive repetitive nucleotides, are identified by SEQCLEAN.  As a part of the quality control, we also excluded unmapped reads mapped onto the rRNA repeat sequence (HSU13369 Human ribosomal DNA complete repeating unit) (BLAST+ 2.2.30). We prepared the index from rRNA repeat sequence using makeblastdb and makembindex from BLAST+.  We used the following command for makeblastdb: 

```
makeblastdb -parse_seqids -dbtype nucl -in `<fasta file>`. 
```

We used the following command for makembindex: 

```
makembindex -input <fasta file> -output <index> -iformat blastdb
```

## B. Mapping unmapped reads onto the human references. 
We remapped the unmapped reads to the human reference sequences using Megablast (BLAST+ 2.2.30). We mapped reads onto the following references:

- Reference transcriptome (known transcripts), Ensembl GRCh37
- Reference genome, hg19 Ensembl

We prepared the index from each reference sequence using makeblastdb and makembindex. We mapped the reads separately onto each reference in the order listed above. Reads mapped to the reference genome and transcriptome were merged into a ‘lost human reads’ category. The following options were used to map the reads using Megablast: for each reference: 

```
task = megablast, use_index = true, perc_identity = 90, outfmt = 6, max_target_seqs =1, e-value = 1e-05.  
```



## C. Identification of hyper-edited reads 
We have used hyper-editing pipeline (HE-pipeline http://levanonlab.ls.biu.ac.il/resources/zip), which is capable of identifying hyper-edited reads.  When running HE-pipeline, additional changes can be made to parallelize the scripts for use with UCLA's Hoffman2 cluster. Before proceeding, follow the instructions in the README that is included with the scripts to prepare the reference and provide the necessary third-party tools. Ensure that the output directory is set correctly in config_file.sh (it is acceptable to use a single output directory), and check that the list of input files has been prepared correctly.

Details on how to run HE-pipeline are available here: 
https://github.com/smangul1/rop/wiki/How-to-run-hyper-editing-pipeline 



## D. Mapping unmapped reads onto the repeat sequences

We filtered out the reads that failed QC and lost human reads. The remaining reads were mapped to the reference repeat sequences.  The reference repeat sequences were downloaded from Repbase v20.07 (http://www.girinst.org/repbase/). Human repeat elements (humrep.ref and  humsub.ref) were merged into a single reference. We prepared the index from the merged repeat reference using makeblastdb and makembindex from BLAST+. In total, we obtained sequences for 1,117 repeat elements. The following options were used to map the reads using the Megablast: 

```
task = megablast, use_index = true,  perc_identity = 90, outfmt = 6, max_target_seqs = 1, e-value = 1e-05. 
```

Blast hits with alignment length shorter than 80% of the read length were discarded (corresponding to 80bp of the 100bp read). 

The repeat elements from humrep.ref and humsub.ref were classified into families and classes using RepeatMasker annotations (hg19_rmsk_TE_prepared_noCDS.bed). Repetitive reads identified from the unmapped reads were confirmed by directly applying Repeatmasker (Tarailo-Graovac & Chen, 2009).

## E. Workflow to detect ‘non-co-linear’ reads (trans-splicing, gene fusions, and circRNAs)

We divide non-co-linear reads into three categories: 

-	gene fusion characterized by reads that map on different chromosomes
-	trans-splicing events characterized by reads that map on the same chromosome, but are at least 1 Mb apart from each other 
-	circRNAs characterized by reads that map in a head-to-tail configuration on the same chromosome

To distinguish between these three categories, we make use of circExplorer2 (Zhang et al., 2016), which was recently identified as one of the best tools to detect circRNAs (Hansen et al., 2015). CircExplorer2 relies on Tophat-Fusion and thus allows also the monitoring NCL events in the same run. TopHat-Fusion (v2.0.13, bowtie1 v0.12.9) and circExplorer2 (v2.2.4) were invoked with the following commands:

```
tophat2 -o tophat-output-directory -p 4 --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search bowtie1-index fastq-file
````

```
python CIRCexplorer2 parse -t TopHat-Fusion -o circrna-output-folder  tophat-output-directory/accepted_hits.bam		
```
```
python CIRCexplorer2 annotate -r ensemble-reference.txt -g genome.fa circrna-output-folder  
```

To separate potential gene and trans-fusions from the TopHat-Fusion output, we ran a ruby custom script, which is part of the ROP pipeline.

## F. Mapping unmapped reads onto the V(D)J recombinations of B and T cell receptors

Gene segments of B cell receptors (BCR) and T cell receptors (TCR) were imported from IMGT (International ImMunoGeneTics information system): (http://www.imgt.org/vquest/refseqh.html#V-D-J-C-sets). 
IMGT database contains:

- Variable (V) gene segments
- Diversity (D) gene segments
- Joining (J) gene segments 

Unmapped reads categorized by step (A)-(D) were filtered out. We used IgBLAST (v. 1.4.0) with stringent e-value threshold (e-value < 10-20) to map the remaining high-quality unmapped reads onto the V(D)J regions of the of the BCR and TCR loci.  Reference files with BCR and TCR VDJ gene segments are distributed with ROP protocol and available at https://drive.google.com/folderview?id=0Bx1fyWeQo3cOTkhKdHFDb3c5MjA&usp=sharing

We prepared the index from each reference sequence using makeblastdb and makembindex from BLAST+. The following options were used to map the reads using IgBLAST: 

```
-germline_db_V; germline_db_D; -germline_db_J; -organism=human; -outfmt = 7; –evalue = 1e-20
```

## G. Identification of microbial reads
Unmapped reads mapping in step (A -(E) were filtered out. The remaining reads were high-quality non-human reads used to profile the taxonomic composition of the microbial communities. We used MetaPhlAn2 (Metagenomic Phylogenetic Analysis, v 2.0) to assign reads on microbial genes and to obtain a taxonomic profile. The database of the microbial marker genes is provided by MetaPhlAn. We run MetaPhlAn in two stages as follow: the first stage identifies the candidate microbial reads (i.e. reads hitting a marker), while the second stage profiles metagenomes in terms of relative abundances – the commands used are as follow:

```
metaphlan.py <fastq> <map> --input_type multifastq --bowtie2db bowtie2db/mpa -t reads_map --nproc 8 --bowtie2out 
metaphlan.py --input_type blastout <bowtie2out.txt> -t rel_ab <tsv>
```

The output of the first stage is a file containing a list of candidate microbial reads with the microbial taxa assigned (.map file). The second stage outputs the taxonomic profile (taxa detected and its relative abundance, in tab separated format (.tsv file). We used taxa detected from stage 2 to extract the reads associated with it in stage 1.  

In addition to MetaPhlAn2 we used to create the curated database of taxa-specific genes, we mapped the reads onto the entire reference genomes of microbial organisms. We used Megablast (BLAST+ 2.2.30) to align reads onto the collection of bacterial, viral, and eukaryotic pathogens reference genomes. Bacterial and viral genomes were downloaded from NCBI ftp://ftp.ncbi.nih.gov/ on February 1, 2015.  Genomes of eukaryotic pathogens were downloaded from EuPathDB database, which is available at: http://eupathdb.org/eupathdb/. 
The following parameters were used for the megablast alignment:  

```
e-value = 10-5, perc_identity = 90
```

The Megablast hits shorter than 80% of the input read sequence were removed (corresponding to 80bp of the 100bp read). 



