#!/usr/bin/env nextflow


/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

if(params.dir == '') {
	println "directory is null. --Please indicate a valid dir"
	exit 1
}

log.info """\
 =======================================
 R I B O M E T H S E Q   P I P E L I N E
 ======================================
 Usage:
    The typical command for running the pipeline is as follows:
	nextflow run main.nf --results_path 'path/to/results_directory' --dir '/path/to/fastq/files/*.gz' -profile {human or mouse}
 Mandatory arguments:
 	--results_path: 		path to where the results will be saved
	--dir				path of the sequenced fastq files
	-profile			human or mouse profile
 directory: ${params.dir}
 index        : ${params.bowtie_index}
 profile       : ${workflow.profile}
 result_path	: ${params.results_path}
 """
 
 

Channel
    .fromPath( params.dir )
    .map {file -> tuple(file.baseName.tokenize("_")[0], file)}	
    .into { read_pairs_ch; read_pairs2_ch; read_pairs3_ch }


results_path = params.results_path
if(results_path == '') {
	println "output directory is null. --Please indicate a valid results_path"
	exit 1
}

/*
 * process reformat {
 *	tag "$dataset_id"
 *       publishDir "$results_path/reformat"
 *       input:
 *	set datasetID, file(reads_3)     from read_pairs3_ch
 *
 *	output:
 *	set datasetID, file("${datasetID}_reformat.fastq") into reformated_files
 *	
 *	"""
 *	reformat.sh in1=$reads_3 out1=${datasetID}_reformat.fastq samplereadstarget=7000000
 *	"""
 * }
 */

process fastq_1 {
	conda params.conda + "riboFastQ/"
	tag "$dataset_id"
	publishDir "$results_path/fastqc_before"
	input:
	set datasetID, file(reads_2)	 from read_pairs2_ch
	
	output:
	file '*_fastqc.{zip,html}' into fastqc_1, fastqc_12
		
	"""
	fastqc -t 1 $reads_2
	"""
}



/*
 * trimmomatic: trim and crop Illumina (FASTQ) data and remove adapters
 * SE: signle end
 * phred33 : the base quality encoding
 * threads: indicates the number of threads to use
 * ILLUMINACLIP: <fastaWithAdaptersEtc: specifies the path to a fasta file containing all the Illumina adapters>:
 * < seedMismatches: specifies the maximum mismatch count which will still allow a full match to be performed>:
 * < palindromeClipThreshold: specifies how accurate the match between the two 'adapter ligated' reads must be for PE palindrome read alignment>:
 * < simpleClipThreshold: specifies how accurate the match between any adapter etc. sequence must be against a read.>
 * LEADING:<quality> Remove low quality bases from the beginning.
 * TRAILING:<quality> Remove low quality bases from the end
 * SLIDINGWINDOW:<windowSize>:<requiredQuality> Perform a sliding window trimming, cutting once the average quality within the window falls below a threshold.
 * AVGQUAL: Drop the read if the average quality is below the specified level
 * MINLEN:<length> This module removes reads that fall below the specified minimal length.
 */

process trimming {
	conda params.conda + "riboTrimm/"
	tag "$datasetID"
	publishDir "$results_path/$datasetID", pattern : "*.fastq"
	publishDir "$results_path/trim_logs", pattern : "*.trimmomatic.stats.log", mode: 'copy'
	input:
	set datasetID, file(datasetFile)	 from read_pairs_ch
	
	output:
	set datasetID, file("${datasetID}.fastq") into trimmed_files,fastq_files_2
	set datasetID, file("${datasetID}.trimmomatic.stats.log") into trimmomatic_logs, trimmomatic_logs2
	
    """
    trimmomatic SE -phred33 -threads ${params.trimmomatic_threads} ${datasetFile} ${datasetID}.fastq ILLUMINACLIP:${params.adapters}:2:30:10 LEADING:${params.LEADING} TRAILING:${params.TRAILING} SLIDINGWINDOW:${params.SLIDINGWINDOW} AVGQUAL:${params.AVGQUAL} MINLEN:${params.MINLEN} 2> ${datasetID}.trimmomatic.stats.log
    """
}

process fastq_2 {
	conda params.conda + "riboFastQ/"
	tag "$datasetID"
	publishDir "$results_path/fastqc_after"
	input:
	set datasetID, file(reads_qc)	 from fastq_files_2
	
	output:
	file '*_fastqc.{zip,html}' into fastqc_results
	
	"""
	fastqc -t 1 $reads_qc
	"""
}





process align_bowtie2 {
	conda params.conda + "riboBowtie2/"
	tag "$datasetID"
	publishDir "$results_path/$datasetID", pattern : "*.sam"
	publishDir "$results_path/bowtie2_logs", pattern : "*.bowtie2.stats.log", mode: 'copy'

	input:
	set datasetID, file(align_file)	 from trimmed_files
	output:
	set datasetID, file("${datasetID}.sam") into bowtie_files
	set datasetID, file("${datasetID}.bowtie2.stats.log") into bowtie_logs, bowtie_logs2
	
	"""
	PERL5LIB=$CONDA_PREFIX/lib/5.26.2
	bowtie2 -x ${params.bowtie_index} ${params.bowtie} $align_file -S ${datasetID}.sam 2>> ${datasetID}.bowtie2.stats.log 
	"""
}

process multiqc {
	conda params.conda + "riboMultiQC/"
	tag "$datasetID"
        publishDir "$results_path", mode: 'copy'
        input:
        file('*') from fastqc_1.collect()
        file('*') from trimmomatic_logs.collect()
		file('*') from bowtie_logs.collect()
        output:
        file('multiqc_report.html')

        """
        multiqc .
        """
}

process count5 {
	conda params.conda + "riboTools/"
	tag "$datasetID"
	publishDir "$results_path/$datasetID", mode:'copy'
	input:
	set datasetID, file(aligned_sam) from bowtie_files
	output:
	set datasetID, file("${datasetID}_3_counts.csv") into counts_3_end
	set datasetID, file("${datasetID}_5_counts.csv") into counts_5_end

	script:
	"""
	samtools view -bS ${datasetID}.sam | samtools sort -o tem_${datasetID}.bam
	samtools view -h -F 4 -b tem_${datasetID}.bam > mapped_RNA_${datasetID}.bam
	samtools view -h mapped_RNA_${datasetID}.bam > mapped_RNA_${datasetID}.sam
	rm ${datasetID}.sam
	rm mapped_RNA_${datasetID}.bam
	rm tem_${datasetID}.bam
	grep -E "@|NM:" mapped_RNA_${datasetID}.sam | grep -v "XS:" > unique_rRNA_${datasetID}.sam
	rm mapped_RNA_${datasetID}.sam
	samtools view -Sb  unique_rRNA_${datasetID}.sam > unique_rRNA_${datasetID}.bam
	bedtools genomecov -d -3 -ibam unique_rRNA_${datasetID}.bam > ${datasetID}_3_counts.csv
	bedtools genomecov -d -5 -ibam unique_rRNA_${datasetID}.bam > ${datasetID}_5_counts.csv
	rm unique_rRNA_${datasetID}.bam
	"""
}


process r_refine {
	tag "$datasetID"
	publishDir "$results_path/$datasetID", mode: 'copy'
	input:
	set datasetID, file(counts3) from counts_3_end
	set datasetID, file(counts5) from counts_5_end
	output:
	set datasetID, file("Treatment_5.8S_Sample_${datasetID}.csv") into final_58S
	set datasetID, file("Treatment_18S_Sample_${datasetID}.csv") into final_18S
	set datasetID, file("Treatment_28S_Sample_${datasetID}.csv") into final_28S
	set datasetID, file("Treatment_5S_Sample_${datasetID}.csv") into final_5S
		
	script:	
	"""
	Rscript ${params.r_script} $counts5 $counts3 ${params.fasta_58S} ${params.fasta_18S} ${params.fasta_28S} ${params.fasta_5S} ${datasetID}
	"""
}

process r_export_run {
	conda params.conda + "riboQC/"
	tag "$datasetID"
	publishDir "$results_path", mode: 'copy'
	input:
	file('*') from fastqc_12.collect()
	file('*') from trimmomatic_logs2.collect()
	file('*') from bowtie_logs2.collect()
	output:
	file("pipeline_report.html")
	file("QCtable.csv")
	file("overrep.fasta")
	
	script:
	"""
	Rscript ${params.r_export} ${params.install_path} .
	"""
	
}


