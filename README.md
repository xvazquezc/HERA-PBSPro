# Introduction
This is a PBSPro (Protable Batch System) based pipeline to perform HERA (Highly Efficient Repeat Assembly) with the scripts "pipline.sh", "04-Qsub-Mapping2Ctg.pl", "08-qsub_job_index.pl", "09-Qsub-Pair_Alignment.pl", and "21-Daligner_New.pl" has been modified and re-saved as "PBS_pipline.sh", "04-PBS-Qsub-Mapping2Ctg.pl", "08-PBS-Qsub_job_index.pl", "09-PBS-Qsub-Pair_Alignment.pl", and "21-PBS-Daligner_New.pl"

# HERA (Highly Efficient Repeat Assembly)
HERA is a local assembly tool using assembled contigs and self-corrected long reads as input. HERA is highly efficient using SMS data to resolve repeats, which enables the assembly of highly contiguous genomes. With the help of BioNano genome maps and chromosomal anchoring information, HERA can generate ultra-long, even chromosome-scale, contigs.

It is important to note that even though HERA can be used to improve the sequence contiguity of highly heterozygous genomes, it require HiC data (and better also with BioNano data) to resolve the haplotype sequences.

# Installation

The running of HERA requires a few other software programs. I suggest to put them in the HERAFile.
1. Downloading and installing bwa-0.7.10
```   
   mkdir HERAFile && cd HERAFile
   git clone https://github.com/lh3/bwa.git  
   cd bwa; make.
```
2. Downloading and installing DALIGNER
```
   git clone https://github.com/thegenemyers/DALIGNER.git
   cd DALIGNER
   make
```
3. Downloading and installing DAZZ_DB
```
   git clone https://github.com/thegenemyers/DAZZ_DB.git
   cd DAZZ_DB
   make
```
4. Downloading and installing HERA
```
   git clone https://github.com/Github-Yilei/HERA.git
   cd HERA
   chmod 777 ./*
```
5. Notice

Because of the different between C99 and the older versions, there is a Error that "'for' loop initial declarations are only allowed in C99 mode".
We need to transform it from the C99 style to the C90 style and `make`.
```
    cd DALIGNER/DAZZ_DB
    grep -n  "for (int" ./*

    #C99
    for (int i = 0; i < n; ++i)
　　do();

    #C90
    int i ;
    for ( i = 0; i < n; ++i)
    do();
```

## With conda

It assumes you have the `bioconda` and `conda-forge` as part of your conda channels.

```
	conda create -n hera bwa daligner dazz_db
```

# Quick Start

### Step 0: Correct the noisy long reads by CANU and finish genome assembly by CANU or MECAT, or FALCON or other assemblers to generate contigs with high sequence accuracy.
The example data of running HERA program is included in HERA.

### Step 1: Create a config file

Before running HERA, you need to create a config file template. HERA provides two kinds of running patterns for connecting the whole-genome assembled contigs and filling the gaps between the paired contigs with or without the BioNano maps.  

The template looks like

```
############################### Variable information ##########################################
#Path of output folder (needs to be created before running)
outfolder=path-to/output

#Path to executables (e.g. conda bin/)
HERAFile=path-to/HERAFile

#the genome name(less 5 words)
genome_name=your_genome_name

#the whole genome assembled sequences with absolute path
genome_seq=path-to/work_shop/Test_Genome.fasta

#the corrected pacbio file with absolute path
Corrected_Pacbio=path-to/work_shop/Test_CorrectedPacbio.fasta
```
you need to fill and modify the relevant information, such as the whole genome assembled contigs or scaffold and the self-corrected long reads.

### Step 2: Running HERA

```Shell
sh pipeline.sh
```

# Results

After the successful submission of pipeline.sh, HERA will take a few steps to get the reassembled genome sequences with the name of "genome_name-Final_Genome_HERA.fasta". HERA mainly includes the following five parts:
1. Mapping the corrected pacbio long reads to the whole genome assembled contigs;
2. Filtering the corrected pacbio long reads which are used to assemble the contigs;
3. Constructing the Contig-Reads and Reads-Reads overlaping graph;
4. Traversing the overlapping graph taking the contig nodes as start and end to find the connecting paths;
5. Constructing and traversing the contig-to-contig path graph to define the order and orientation of the contigs;
6. Constructing the consensus sequence to fill the gap and produce the final genome.

Finally, the users can get the combined seq, the super-contig genome and the connection information by HERA in the work_shop/genome_name-Final_Genome_HERA.fasta, 06-Daligner/SuperContig.fasta and 06-Daligner/Ctg_Position.txt.

# Details of HERA pipeline

Please go to the mainpake of [liangclab](https://github.com/liangclab/HERA) for more information.

# Usage Limitations

HERA is highly efficient for generating highly contiguous and complete or nearly complete sequences for small genomes such as fungi as well as homozygous genomes. HERA may be applied to a genome for several rounds to get desired results. For highly heterozygous genomes, a lot of manual work may be required.

# Citing HERA

Du, H., Liang, C. (2018). Assembly of chromosome-scale contigs by efficiently resolving repetitive sequences with long reads. bioRxiv doi: https://doi.org/10.1101/345983

#  Seeking Help

The detailed usage is described in the man page available together with the source code. If you have questions about HERA, you may send the questions to cliang@genetics.ac.cn.
