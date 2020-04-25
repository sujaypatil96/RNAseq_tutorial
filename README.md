# RNA-seq Data Analysis Pipeline
The following tutorial describes the various steps involved in the development/construction of a bioinformatics pipeline for the analysis of RNAseq data. The main topics/concepts covered as part of this tutorial are: Analysis of data using TopHat2, Isoform Analysis, Pathway Analysis.
> Note: Consider Module 0 as a pre-processing step, which includes collating all the metadata about the disease that has infected your samples under analysis, sequencing information and metrics and also some basic quality control steps.

<hr>

##### Module 0
### 1. **Biology/Introduction**
> In this step we will describe the biology/human disease under analysis.
- What is the disease / disorder?

The disease/disorder that we will be considering for the analysis pipeline is 'Alzheimer's Disease'. The RNA-seq dataset was obtained from GEO. You can find the [FASTQ files for analysis here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113524). There are 39 samples that the research team analyzed as part of an experimental study. You can [access the paper here](https://www.nature.com/articles/s41380-019-0563-5). For our analysis we will consider only 9 of these samples and we will split into 3 groups on the basis of their Braak stage. The 9 selected samples have been grouped as follows:

| Sample # | SRA record           | Sample| Age  | Sex  | Braak stage  | Group  |
| ------------- |:-------------:|:-----:|:-----:|:-----:|:-----:|-----:|
| 5 | SRR7056889 | Healthy5 RNA-Seq | 85 | M | 1 | Early (CTRL) |
| 6      | SRR7056890 |   Healthy6 RNA-Seq | 88 | M | 1 | Early (CTRL) |
| 38 | SRR8942876 |    Healthy38 RNA-Seq | 94 | M | 1 | Early (CTRL) |

| Sample # | SRA record           | Sample| Age  | Sex  | Braak stage  | Group  |
| ------------- |:-------------:|:-----:|:-----:|:-----:|:-----:|-----:|
| 20 | SRR7056904 | Patient10 RNA-Seq | 88 | M | 4 | Mid AD (MID_BRAAK) |
| 35      | SRR8942873 | Patient35 RNA-Seq | 87 | M | 4 | Mid AD (MID_BRAAK) |
| 37 | SRR8942875 | Patient37 RNA-Seq | 95 | M | 3 | Mid AD (MID_BRAAK) |

| Sample # | SRA record           | Sample| Age  | Sex  | Braak stage  | Group  |
| ------------- |:-------------:|:-----:|:-----:|:-----:|:-----:|-----:|
| 5 | SRR7056899 | Patient5 RNA-Seq | 87 | M | 6 | Late AD (LATE_BRAAK) |
| 6      | SRR8942871 | Patient33 RNA-Seq | 96 | M | 5 | Late AD (LATE_BRAAK) |
| 38 | SRR7056901 | Patient7 RNA-Seq | 94 | 82 | 6 | Late AD (LATE_BRAAK) |
- Describe an a priori hypothesis about a gene, pathway, etc. based on the disease/experiment? <span style="color:red">\[this needs to be tested]</span>

We know that the APOE gene mutation, speciifically, the E4 allele of the APOE gene is one of the most firmly established risk factors of late onset Alzheimer's Disease. So we will simply be testing the hypothesis that the APOE gene is a significantly differentially expressed gene in Alzheimer's related pathways.

![alt text](https://github.com/sujaypatil96/rnaseq-pipeline/blob/master/assets/images/TRGN515_Sujay_Patil-1.jpg)
![alt text](https://github.com/sujaypatil96/rnaseq-pipeline/blob/master/assets/images/TRGN515_Sujay_Patil-2.jpg)

### 2. **Sequencing type / source**
> In this step we will describe the sequencing data and source.

- Where did you get the data? Are the sequencing reads single end or paired end? Were the samples run in multiple lanes or in groups that required concatenation?

As mentioned in the above section, the RNA-seq dataset was obtained from GEO (find the link above). As part of the original experiment/study RNA-Seq of the olfactory bulb of 19 AD patients and 20 age-matched controls were generated by deep sequencing postmortem tissue.

For the purposes of our analysis, we considered 9 samples prior to alignment. As can been seen/summarized in the tables above, the 9 samples were split into 3 groups: CTRL, MID_BRAAK, LATE_BRAAK on the basis of Braak stage.

Sequencing platform: Illumina Hi-seq 2500.

    * Read type: single read
    * Length of oligomer (k-mer): 60 bp 
    * No. of lanes: 4 lanes of an Illumina Hi-Seq
    * Sequencing depth:~17 M reads/sample

### 3. **Pre-alignment sequencing metrics**
- Describe finalised number of samples before alignment.
- A graph (of your choice) and description of mean, minimum and maximum # of reads <span style="color:red">\[use terms ‘read’ or ‘read pairs’ correctly]
 
![alt text](https://github.com/sujaypatil96/rnaseq-pipeline/blob/master/assets/images/TRGN515_Sujay_Patil-4.jpg)
 
<hr>

##### Module 1A
### 4. **Tophat2 - alignment process description**
- Describe the alignment process/algorithm employed by Tophat2 for RNA-seq alignment.
> Note: Not necessary to understand ‘nitty gritty’ details about Tophat, but should demonstrate basic process: 
i. Starting with unaligned reads (FASTQ)
ii. Aligning to a reference genome (FASTA)
iii. Outputting a file of aligned reads (BAM) and associated statistics.

_Note: For our analysis we will use the human GRCh38 version of the genome from Ensembl._

From within (or outside, it is up to you) the directory that contains your fasta, and gtf files, use the command bowtie2-build to make index files.

Command usage is as follows:

    bowtie2-build FASTA-reference desired-prefix-of-index-files

First we will we be using bowtie2 to build index of the reference (in this human) genome.

In the end, your FASTA reference and bowtie index files must have the same prefix. So you can either name the indexes with the title already provided with the FASTA reference, or you can make a new name, and then rename the FASTA reference afterwards with mv to match the index files. Both example are below:

Keeping current FASTA naming convention:

    bowtie2-build Homo_sapiens.GRCh38.dna.primary_assembly.fa Homo_sapiens.GRCh38.dna.primary_assembly &

Changing naming convention (what I did):

    bowtie2-build Homo_sapiens.GRCh38.dna.primary_assembly.fa ensembl.GRCh38.99 &

    mv Homo_sapiens.GRCh38.dna.primary_assembly.fa ensembl.GRCh38.99.fa

You can start 6x TopHat2 alignments:

1) Use 2 single-end (SE) fastq files of your choosing

2) We will have you do 3x separate alignments with these two samples (for total of 6x alignments)- one where we align to the entire genome, one where we we only align to the whole transcriptom, and another where you strictly align to ribosomal RNA.. then you can compare alignment metrics in class. 

Command Line for SE alignemnt (whole genome):

    nohup tophat2 -p 4 -G /scratch/sujaysan/515/References/Homo_sapiens.GRCh38.99.chr.gtf -o name_of_desired_output_directory /scratch/sujaysan/515/References/ensembl.GRCh38.99 File1.fastq > NAME_genome.nohup.out &

Command Line for SE alignemnt (transcriptome only):

    nohup tophat2 -p 4 -G /scratch/sujaysan/515/References/Homo_sapiens.GRCh38.99.chr.gtf -T -o name_of_desired_output_directory /scratch/sujaysan/515/References/ensembl.GRCh38.99 File1.fastq > NAME_transcriptome.nohup.out &

Command Line for SE alignemnt (rRNA only):

    nohup tophat2 -p 4 -G /scratch/sujaysan/515/References/rRNA_Homo_sapiens.GRCh38.99.chr.gtf -T -o name_of_desired_output_directory /scratch/sujaysan/515/References/ensembl.GRCh38.99File1.fastq > NAME_rRNA.nohup.out &

List of command options and files (in order of appearance):

    nohup === don't hang up
    tophat2 === run tophat2 (v2.1.1)
    -p 4 === use 4 threads
    -G Homo_sapiens.GRCh38.99.chr.gtf === use this gtf file as a transcriptome guide and to make transcriptome index files (note I gave entire path to where this file is since you are not working in that directory)
    -T == align to transcriptome only (this option is only included in the last two commands as that will dictate that only things that align to your given GTF (transcriptome or rRNA) will be kept!!)
    -o OUTPUT_DIRECTORY === whatever you want the directory to be called where all of your alignment files to go (This must be unique for each run or else will overwrite the file!!!)
    ensembl.GRCh38.99 === prefix used for the FASTA reference and bowtie2 index files (Only give the prefix- do not give file extenstions like .fa or .bt2!!!) (note I gave entire path to where this file is since you are not working in that directory)
    File1.fastq === your chosen fastq (if SE)
    ** make sure to repeat for 2 samples (i.e., you should be doing 6x alignments)
    NAME.nohup.out & === Name of your nohup file (This must be unique for each run or else will overwrite the file!!! This will record all the steps and what happened during the alignment so is a very good things to have recorded)`

_Note: Building index files could take up to an hour- I highly suggest using nohup if you don't want to keep your terminal open or are worried about internet connectivity._

FASTQC is primarily for pre-alignment and it takes as input FASTQ or FASTA files. To make sure sequence content, sequence quality, sequence representation (no over-representation of adapters), and KMER representation are all adequate for alignment.
As far as pre-alignment quality control, we can manually load fastq files in FastQC to make sure the average sequence quality score is at least 25. I would like to trim any overrepresented sequences (usually from adapters) and long mononucleotide repeats (length threshold is not yet defined).

To trim sequence adapter from the read FASTQ files. The output of this step will be trimmed FASTQ files for each data set.
Once you have trimmed the reads, compare a pre- and post- trimming FastQ file using the FastQC tool.
My typical workflow is to run FastQC pre- and post-trimming with trimmomatic to compare, though it isn't always necessary to do any trimming depending on the usage/quality of the data.

Genes are the functional units of a reference genome and gene annotations describe the structure of transcripts expressed from those gene loci. GTF (Gene Transfer Format) is a common file format used to store gene and transcript annotation information. If annotation is available as a GTF file, TopHat will extract the transcript sequences and use Bowtie2 to align reads to this virtual transcriptome first. Only the reads that do not fully map to the transcriptome will then be mapped on the genome. The reads that remain unmapped are split into shorter segments, which are then aligned to the genome. Segment mappings are used to find potential splice sites.

![alt text](https://github.com/sujaypatil96/rnaseq-pipeline/blob/master/assets/images/TRGN515_Sujay_Patil-5.jpg)

This tool returns the following file: `tophat.bam`: BAM file containing the alignments

If your RNA-seq data was produced with a stranded/directional protocol, it is important that you select the correct strandedness option in the parameter "Library type": fr-unstranded (Standard Illumina).

### 5. **Tophat2 - alignment metrics and statistics**
- Information about the alignment process and output.
- Preapre a graph (of your choice) and description of mean, minimum and maximum for the following items:
i. % Alignment (# aligned reads / # initial reads)
ii. # of Aligned Reads
iii. % Discordance

_Note: $ discordance is a concept that applies only to paired-end (PE) data, we do not need to worry about it since our data is single-end (SE)._

![alt text](https://github.com/sujaypatil96/rnaseq-pipeline/blob/master/assets/images/TRGN515_Sujay_Patil-6.jpg)

##### Module 1B
### 6. **Cuffdiff results**
- Select a list of genes that were differentially expressed in your cuffdiff_exp file from the analysis you did with the stats test turned ON.
- Get the FPKM values for ALL of your samples in your genes_fpkm.tracking file from the analysis you did with the stats turned OFF.
- Limit FPKM list to genes selected from step 1.

- Make a Heatmap, Dendrogram, and PCA plot(s) of the filtered FPKM values – For dendrogram and PCA plot, plot SAMPLES not GENES. For Heatmap you can plot SAMPLES, GENES or BOTH.
\[Understand how the samples cluster using these plots. Do they group the way you would expect based on the conditions? What is your interpretation of this view of the data?]

<hr>

##### Module 2
> Note: The steps in this module are not necessary, but can be useful to help identify the best aligner for your pipeline.

### 7. **HiSat2 or STAR comparison**
- Choose 3-4 of your samples. Align with HiSat2 or STAR. Compare alignment metrics with TopHat2 output for these same files.
- Should include graphs (of your choice) comparing mean for both methods evaluating:
i. % Alignment (# aligned reads / # initial reads)
ii.	# of Aligned Reads
iii. % Discordance
- Based on any/all of the above results, which alignment tool performed better?
- Are there any command options you included in TopHat2 (or HiSat2/STAR) but not the other that make it difficult to compare the results directly or make a conclusion about performance?

<hr>

##### Module 3
### 8. **Pathway analysis**
- Information about upregulated/downregulated gene lists and source and affected pathways.
- Select a list of genes that were differentially expressed and were UPREGULATED in a treatment vs. control comparison. 
- Select a list of genes that were differentially expressed and were DOWNREGULATED in a treatment vs. control comparison.
- Submit each list to IPA for a simple pathway analysis (ORA - overrepresentation analysis). Describe affected pathways that interest you, a hypothesis about why this pathway might be UP or DOWN in the conditions you chose.
