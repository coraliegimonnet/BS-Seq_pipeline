<img src="logo/INRA_logo.jpg" align="right" width="350" height="70"/>

# BS-seq Pipeline

This pipeline was developed for fastq reads from WGBS (Whole Genome Bisulfite Sequencing) - not adapted for RRBS for the moment. 

This pipeline is based on different tools :

- [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for fastq analysis
- [Trim Galore](https://github.com/FelixKrueger/TrimGalore) for trimming (
Warning: Trim galore is used as recommended in manual for WGBS, so stringency to cut adapter was fixed at 1bp)
- [Bismark](https://github.com/FelixKrueger/Bismark) for alignment and bedGraph creation
- [MethylKit](http://bioconductor.org/packages/release/bioc/manuals/methylKit/man/methylKit.pdf) for methylation extraction

**Warnings** : This pipeline assumes all software are installed in PATH

It is writen in groovy language with [Nextflow framework](https://www.nextflow.io/).

It takes as input a fasta genome file and raw sequenced reads and make all bioinformatic step until methylation extraction.

Print help of BS-Seq_pipeline.nf :

	nextflow run BS-Seq_pipeline.nf --help


Parameters :

	--genome:	fasta file of genome 
	--reads:	paired end reads need to be write like this 'toto_R{1,2}.fq.gz' (surrounded by quotes)
	--outdir:	directory where all results are placed
	--cov:		minimum coverage for methylation extraction

Optional parameters :

	--index: path to bismark index files [default = FALSE]
	--se :		single-end data [default = paired-end]
	--fastqc:	make a fastQC analysis [default = FALSE]
	--notrim:	skip trimming step [default = FALSE]
	--bedGraph:	bedGraph creation [default = FALSE]
	--extract: methylation extraction with MethylKit [default = FALSE]


Example of usage :

    nextflow run BS-Seq_pipeline.nf --genome 'galGal6_all.fa' --reads '*_R{1,2}*' --outdir 'WGBS_results'


If an incident happens, you could rerun your command line with `-resume` option.
  
