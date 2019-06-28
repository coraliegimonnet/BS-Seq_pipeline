#!/usr/bin/env nextflow

/*
####################################################################
### BS-Seq Methylation Analysis : Trimming, Mapping and Cleaning ###
####################################################################
Script started in July 2018. 
Authors :
    Coralie Gimonnet <coralie.gimonnet@inra.fr>
This pipeline was adapted for WGBS analysis with paired-end data.
Version 1.6
*/


/* Default parameters */

params.genome = "genome.fasta"
params.reads = "*_R{1,2}.fq.gz" 
params.outdir ='results'
params.cov = 0

params.notrim = false
params.help = false
params.fastqc = false
params.index = false
params.se = false
params.extract = false


read_pairs = file(params.reads)
fasta = file(params.genome)
 
log.error """\
			BISMARK ALIGNMENT PIPELINE
		========================================
		default parameters :
			genome : ${params.genome}
			reads  : ${params.reads}
			outdir : ${params.outdir}
			cov : ${params.cov}
		optional parameters :
			index : path to bismark index files [default : FALSE]
			se : single-end data [default : paired-end] - se not optimized for the moment
			fastqc : make a FastQC analysis [default : FALSE]
			notrim : skip trimming step [default : FALSE]
			bedGraph : make a bedGraph of alignment - bedGraph creation not fonctionnal for the moment
			extract : extract methylation [default : FALSE]
		Warning : This pipeline assumes all software are on the PATH.
		"""
		.stripIndent()


if (params.help) exit 1

 
/* 
* Create the `read_pairs` channel that emits tuples containing three elements :
* the pair ID, the first read-pair file and the second read-pair file
*/

Channel
	.fromFilePairs( params.reads, size: params.se ? 1 : 2 )
	.ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}" }
	.into { read_files_fastqc ; read_files_trimming ; reads_mapping } 


/*
* Step 1. Builds the genome index required by the mapping process
*/
if (params.index)
{
	bismark_index = Channel
		.fromPath(params.index)
		.ifEmpty { exit 1, "Bismark index not found : ${params.index}" }
}
else
{
	process buildIndex {
		publishDir "${params.outdir}/reference_genome", mode : 'copy'
		input:
			file "genome.fasta" from fasta
		output:
			file "BismarkIndex" into bismark_index
		script:
		"""
			mkdir BismarkIndex
			cp ${fasta} BismarkIndex/
			bismark_genome_preparation --verbose BismarkIndex
		"""
	}
}

 
/*
* Step 2. FASTQC Analysis
*/
if (params.fastqc)
{
	process fastqc {
		tag "$pair_id"
		publishDir "${params.outdir}/fastqc", mode: 'copy',
			saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}
		input:
			set pair_id, file(reads) from read_files_fastqc
		output:
			file '*_fastqc.{zip,html}' into fastqc_results
		script:
		"""
			fastqc -q $reads
		"""
	}
}
else
{
	fastqc_results = Channel.from(false)
}


/*
* Step 3. Trim Galore!
*/
if (params.notrim){
	trimmed_reads = read_files_trimming
	trimgalore_results = Channel.from(false)
}
else
{
	if (params.fastqc) 
	{
		process trimming {
			tag "$pair_id"
			publishDir "${params.outdir}/trim_galore", mode : 'copy',
				saveAs: {filename -> 
					if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
					else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
				}
			input:
				set pair_id, file(reads) from read_files_trimming
			output:
				set pair_id, file('*.fq.gz') into trimmed_reads
				file "*trimming_report.txt" into trimgalore_results
				file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports
			script:
				if (params.se)
			    {
			        """
			            trim_galore --fastqc --gzip $reads
    			    """
	    		}
		    	else
			    {
					"""
						trim_galore --paired --fastqc --gzip $reads
					"""
				}
		}
	}
	else
	{
		process trimming {
			tag "$pair_id"
			publishDir "${params.outdir}/trim_galore",mode : 'copy',
				saveAs: {filename ->
					if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
				}
			input:
				set pair_id, file(reads) from read_files_trimming
			output:
				set pair_id, file('*.fq.gz') into trimmed_reads
				file "*trimming_report.txt" into trimgalore_results
			script:
				if (params.se)
			    {
			        """
			            trim_galore --gzip $reads
			        """
			    }
			    else
			    {
					"""
						trim_galore --paired --gzip $reads
					"""
				}
		}
	}
}	

/*
* Step 4. Mapping with Bismark
*/
process mapping {
	cpus 16
	tag "$pair_id"
	publishDir "${params.outdir}/mapping", mode : 'copy',
		saveAs: {filename ->
			if (filename.indexOf("report.txt") > 0) "reports/$filename"
		}
	input:
		file index from bismark_index.collect()
		set pair_id, file(reads) from trimmed_reads
	output:
		file "*.bam" into bam_aln
		file "*report.txt" into bismark_align_log
	script:
		if (params.se)
	    {
	        """
	            bismark --fastq --genome $index $reads -p ${task.cpus}
	        """
	    }
	    else
	    {
			"""
				bismark --fastq --genome $index -1 ${reads[0]} -2 ${reads[1]} -p ${task.cpus}
			"""
		}
}

/*
* Step 5. Sorting alignment and stats
*/
process sorting {
	cpus 8
	tag "${bam.simpleName}"
	publishDir "${params.outdir}/sort", mode : 'copy',
		saveAs: {filename -> 
			if (filename.indexOf(".txt") > 0) "logs/$filename" 
			}
	input:
		file bam from bam_aln
	output:
		file "${bam.simpleName}.sorted.bam" into bam_sort
		file "${bam.simpleName}_flagstat.txt" into flagstat_results
	script:
	"""
		samtools sort $bam -@ ${task.cpus} -o ${bam.simpleName}.sorted.bam
		samtools flagstat ${bam.simpleName}.sorted.bam > ${bam.simpleName}_flagstat.txt
	"""
}


/*
* Step 6. Deduplication
*/

process deduplicate {
	cpus 8
	tag "${bam.simpleName}"
	input:
		file bam from bam_sort
	output:
		file "${bam.simpleName}.dedup.bam" into bam_dedup, flag_dedup
		file "${bam.simpleName}.dedup.sortn.bam" into bam_dedup_sort
	script:
		if (params.se)
	    {
	        """
	            samtools rmdup -s $bam ${bam.simpleName}.dedup.bam
	            samtools sort -@ ${task.cpus} -n ${bam.simpleName}.dedup.bam -o ${bam.simpleName}.dedup.sortn.bam
	        """
	    }
	    else
	    {
			"""
				samtools rmdup $bam ${bam.simpleName}.dedup.bam
				samtools sort -@ ${task.cpus} -n ${bam.simpleName}.dedup.bam -o ${bam.simpleName}.dedup.sortn.bam
			"""
		}
}


/*
* Step 7. Cleaning 
*/

process cleaning {
	cpus 4
	tag "${bam.simpleName}"
	publishDir "${params.outdir}/cleaned", mode : 'copy',
		saveAs: { filename -> 
					if (filename.indexOf(".cleaned.bam") > 0) "final_bam/$filename"}
	input:
		file bam from bam_dedup_sort
	output:
		file "${bam.simpleName}.dedup.sortn.del.bam" into bam_dedup_sort_del
		file "${bam.simpleName}.dedup.cleaned.bam" into bam_dedup_cleaned, flag_dedup_cleaned, bam_extract, bam_bedgraph
	shell:
	'''
		samtools view -h !{bam} | awk '/^@/{print;next}$1==id{print l"\\n"$0;next}{id=$1;l=$0}' | samtools view -bS -o !{bam.simpleName}.dedup.sortn.del.bam
		samtools sort -@ !{task.cpus} !{bam.simpleName}.dedup.sortn.del.bam -o !{bam.simpleName}.dedup.cleaned.bam
	'''
}


/*
* Step 8. Flagstat of final bam
*/
process flagstat {
	tag "${bam.simpleName}"
	publishDir "${params.outdir}/cleaned", mode : 'copy'
	input :
		file bam from bam_dedup_cleaned
	output:
		file "${bam.simpleName}_flagstat.cleaned.txt" into flagstat_clean
	shell:
	"""
		samtools flagstat $bam > ${bam.simpleName}_flagstat.cleaned.txt
	"""
}
	
/*
* Step 9. Methylation extraction
*/
if (params.extract)
{
	process methylation_extraction {
    	tag "${bam.simpleName}"
	    publishDir "${params.outdir}/Methylation_extraction", mode : 'copy'
    	input:
        	file bam from bam_extract
	    output:
    	    file "${bam.simpleName}_MethylKit" into meth_extract
	    script:
    	"""
        	echo 'library (methylKit)' > script.R
	        echo 'my.methylRaw = processBismarkAln (location = \"$bam\" , sample.id = \"${bam.simpleName}\", assembly = \"Cotja\", read.context = \"CpG\", save.folder = \"${bam.simpleName}_MethylKit\", mincov =${params.cov})' >> script.R
	        Rscript script.R 
	    """
	}
}

/* 
* Step 10. BedGraphs creation
*/
/*if (params.bedGraph)
{
	process bedgraph {
	    cpus 8
	    tag "${bam.simpleName}"
	    publishDir "${params.outdir}/bismark_methylation_calls", mode : 'copy'
	    input:
	        file bam from bam_bedgraph
	    output:
	    	file "${bam.simpleName}.sortn.bam" into bam_sortn 
	        file "${bam.simpleName}_BedGraphs" into bedgraph_out
	    script:
	    	if (params.se)
	    	{
	        	"""
	        	    cat $bam | samtools sort -n -@ ${task.cpus} -o ${bam.simpleName}.sortn.bam
		            bismark_methylation_extractor -s --gzip --multicore ${task.cpus} --bedGraph --scaffolds ${bam.simpleName}.sortn.bam -o ${bam.simpleName}_BedGraphs
		        """
		    }
		    else
	    	{
		    	"""
	    	 	   cat $bam | samtools sort -n -@ ${task.cpus} -o ${bam.simpleName}.sortn.bam 
	    	 	   bismark_methylation_extractor -p --gzip --multicore ${task.cpus} --bedGraph --scaffolds ${bam.simpleName}.sortn.bam -o ${bam.simpleName}_BedGraphs
			    """
		    }
	}
}*/


/* 
* Step 11. Version of all tools used in this pipeline
*/

process software_version {
	publishDir "${params.outdir}/version", mode : 'copy'
	output:
		file 'software_version.txt' into software_version
	script:
	"""
		echo 'FastQC version:' > software_version.txt
		fastqc --version >> software_version.txt
		echo 'Trim Galore! version:' >> software_version.txt
		trim_galore --version >> software_version.txt
		echo 'Cutadapt version:' >> software_version.txt
		cutadapt --version >> software_version.txt
		echo 'Bismark version:' >> software_version.txt
		bismark --version >> software_version.txt
		echo 'Samtools version:' >> software_version.txt
		samtools --version >> software_version.txt
    """
}
