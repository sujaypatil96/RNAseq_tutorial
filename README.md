# RNAseq Tutorial
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
| 5 | SRR7056889 | Healthy5 RNA-Seq | 85 | M | 1 | Early (Control) |
| 6      | SRR7056890 |   Healthy6 RNA-Seq | 88 | M | 1 | Early (Control) |
| 38 | SRR8942876 |    Healthy38 RNA-Seq | 94 | M | 1 | Early (Control) |

| Sample # | SRA record           | Sample| Age  | Sex  | Braak stage  | Group  |
| ------------- |:-------------:|:-----:|:-----:|:-----:|:-----:|-----:|
| 20 | SRR7056904 | Patient10 RNA-Seq | 88 | M | 4 | Mid AD |
| 35      | SRR8942873 | Patient35 RNA-Seq | 87 | M | 4 | Mid AD |
| 37 | SRR8942875 | Patient37 RNA-Seq | 95 | M | 3 | Mid AD |

| Sample # | SRA record           | Sample| Age  | Sex  | Braak stage  | Group  |
| ------------- |:-------------:|:-----:|:-----:|:-----:|:-----:|-----:|
| 5 | SRR7056899 | Patient5 RNA-Seq | 87 | M | 6 | Late AD |
| 6      | SRR8942871 | Patient33 RNA-Seq | 96 | M | 5 | Late AD |
| 38 | SRR7056901 | Patient7 RNA-Seq | 94 | 82 | 6 | Late AD |
- What type of cells / tissue was used?
- Describe an a priori hypothesis about a gene, pathway, etc. based on the disease/experiment? <span style="color:red">\[this needs to be tested]</span>

### 2. **Sequencing type / source**
> In this step we will describe the sequencing data and source.
- Where did you get the data?
- Are the sequencing reads single end or paired end?
- Were the samples run in multiple lanes or in groups that required concatenation?

### 3. **Pre-alignment sequencing metrics**
- Describe finalised number of samples before alignment.
- A graph (of your choice) and description of mean, minimum and maximum # of reads <span style="color:red">\[use terms ‘read’ or ‘read pairs’ correctly]
 
<hr>

##### Module 1A
### 4. **Tophat2 - alignment process description**
- Describe the alignment process/algorithm employed by Tophat2 for RNA-seq alignment.
> Note: Not necessary to understand ‘nitty gritty’ details about Tophat, but should demonstrate basic process: 
i. Starting with unaligned reads (FASTQ)
ii. Aligning to a reference genome (FASTA)
iii. Outputting a file of aligned reads (BAM) and associated statistics.

### 5. **Tophat2 - alignment metrics and statistics**
- Information about the alignment process and output.
- Preapre a graph (of your choice) and description of mean, minimum and maximum for the following items:
i. % Alignment (# aligned reads / # initial reads)
ii. # of Aligned Reads
iii. % Discordance

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
