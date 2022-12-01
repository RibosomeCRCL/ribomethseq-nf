#!/usr/bin/env nextflow

// Copyright (C) 2022 CRCL

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

log.info ""
log.info "--------------------------------------------------------------------"
log.info "      $workflow.manifest.name: $workflow.manifest.description       "
log.info "--------------------------------------------------------------------"
log.info "Copyright (C) CRCL"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------------------"
log.info ""


/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

params.fqdir         = null
params.fastq_pattern = 'fastq.gz'

params.help          = false

// trimmomatic
params.adapters      = "$baseDir/data/adapters/TruSeq3-SE.fa"
params.leading       = 30
params.trailing      = 30
params.slidingwindow = '4:15'
params.avgqual       = 30
params.minlen        = 8
params.threads       = 3

// bowtie
params.bowtie_opts = "--sensitive -L 17"

// r utils
r_refine = "$baseDir/Rscripts/Refine/r_refine.R"
r_export = "$baseDir/Rscripts/QC/r_export.R"

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                             INPUT PARAMETERS CHECK
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

if (params.help) exit 0, helpMessage()

//
// Check and configure input/output parameters
//
if (!params.fqdir) {
	exitMessage '--fqdir option must be set'
}

def outdir = file(params.outdir, type: 'dir', checkIfExists: true).toAbsolutePath()
def logdir = file(params.logdir, type: 'dir', checkIfExists: true).toAbsolutePath()
def fqdir  = file(params.fqdir,  type: 'dir', checkIfExists: true).toAbsolutePath().toString() + '/'


/* Log execution parameters to stdout */
log.info "----------------------- Workflow properties ------------------------"
log.info "Workflow project directory       : $workflow.projectDir"
log.info "Workflow version                 : $workflow.manifest.version"
log.info "Workflow launched from           : $workflow.launchDir"
log.info "Workflow workDir                 : $workflow.workDir"
log.info "Workflow config file(s)          : $workflow.configFiles"
log.info "Workflow profile                 : $workflow.profile"
log.info "Executed command line            : $workflow.commandLine"
log.info ""
log.info "------------------------ Execution parameters ----------------------"
log.info "Fastq input directory            : ${fqdir}"
log.info "Fastq selection pattern          : ${params.fastq_pattern}"
log.info "Trimmomatic ILLUMINACLIP         : ${params.adapters}"
log.info "Trimmomatic LEADING              : ${params.leading}"
log.info "Trimmomatic TRAILING             : ${params.trailing}"
log.info "Trimmomatic SLIDINGWINDOW        : ${params.slidingwindow}"
log.info "Trimmomatic AVGQUAL              : ${params.avgqual}"
log.info "Trimmomatic MINLEN               : ${params.minlen}"
log.info "Bowtie options                   : ${params.bowtie_opts}"
log.info "Threads for fastqc               : ${params.fastqc_threads}"
log.info "Threads for bowtie               : ${params.bowtie_threads}"
log.info "Threads for trimmomatic          : ${params.trimmo_threads}"
log.info "Threads for samtools             : ${params.samtools_threads}"
log.info "Log directory                    : ${logdir}"
log.info "Output directory                 : ${outdir}"
log.info "Max parallel jobs                : ${params.qsize}"
log.info ""


/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                      MAIN
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

Channel
	.fromPath( fqdir + "/*${params.fastq_pattern}" )
	.ifEmpty { error "No fastq file matching \'${params.fastq_pattern}\' in: ${fqdir}" }
	.map { fq -> [
		file(fq).getName().replaceAll("${params.fastq_pattern}\$",''),
		file(fq)
	]}
	.view()
	.into { read_pairs_ch; read_pairs2_ch; read_pairs3_ch }

//exit 1

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

process fastqc_before {
	tag { datasetID }

	publishDir "$outdir/fastqc_before"

	input:
	set datasetID, file(reads_2) from read_pairs2_ch

	output:
	file '*_fastqc.{zip,html}' into fastqc_1, fastqc_12

	"""
	fastqc -t 1 $reads_2
	"""
}

process trim {
	tag { datasetID }

	publishDir "$outdir/$datasetID", pattern : "*.fastq"
	publishDir "$outdir/trim_logs", pattern : "*.trimmomatic.stats.log", mode: 'copy'

	input:
	set datasetID, file(datasetFile) from read_pairs_ch

	output:
	set datasetID, file("${datasetID}.fastq") into trimmed_files,fastq_files_2
	set datasetID, file("${datasetID}.trimmomatic.stats.log") into trimmomatic_logs, trimmomatic_logs2

	"""
	trimmomatic SE -phred33 -threads ${params.trimmo_threads} \\
	${datasetFile} ${datasetID}.fastq \\
	ILLUMINACLIP:${params.adapters}:2:30:10 \\
	LEADING:${params.leading} \\
	TRAILING:${params.trailing} \\
	SLIDINGWINDOW:${params.slidingwindow} \\
	AVGQUAL:${params.avgqual} \\
	MINLEN:${params.minlen} 2> ${datasetID}.trimmomatic.stats.log
	"""
}

process fastqc_after {
	tag { datasetID }

	publishDir "$outdir/fastqc_after"

	input:
	set datasetID, file(reads_qc) from fastq_files_2

	output:
	file '*_fastqc.{zip,html}' into fastqc_results

	"""
	fastqc -t 1 $reads_qc
	"""
}

process bowtie2 {
	tag { datasetID }

	publishDir "$outdir/$datasetID", pattern : "*.sam"
	publishDir "$outdir/bowtie2_logs", pattern : "*.bowtie2.stats.log", mode: 'copy'

	input:
	set datasetID, file(trimmed_file) from trimmed_files

	output:
	set datasetID, file("${datasetID}.sam") into bowtie_files
	set datasetID, file("${datasetID}.bowtie2.stats.log") into bowtie_logs, bowtie_logs2

	"""
	# PERL5LIB="/opt/conda/lib/5.34.0"
	bowtie2 -x ${params.bowtie_index} --threads ${params.bowtie_threads} \\
	${params.bowtie_opts} ${trimmed_file} -S ${datasetID}.sam 2>> ${datasetID}.bowtie2.stats.log
	"""
}

process multiqc {
	tag { datasetID }

	publishDir "$outdir", mode: 'copy'

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
	tag { datasetID }

	publishDir "$outdir/$datasetID", mode:'copy'

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
	tag { datasetID }

	publishDir "$outdir/$datasetID", mode: 'copy'

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
	Rscript ${params.r_refine} $counts5 $counts3 ${params.fasta_58S} ${params.fasta_18S} ${params.fasta_28S} ${params.fasta_5S} ${datasetID}
	"""
}

process r_export {
	tag { datasetID }

	publishDir "$outdir", mode: 'copy'

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
	Rscript ${params.r_export} $baseDir .
	"""
}

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                      UTILS
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

def exitMessage(msg) {
	log.error "${msg}"
	this.helpMessage()
	exit 1
}

def helpMessage() {
	log.info ""
	log.info "                       +++++++++++++++++++++++++"
	log.info "                       +  ribomethseq-nf help  +"
	log.info "                       +++++++++++++++++++++++++"
	log.info ""
	log.info "--fqdir              DIR    Fastq files location                       Required"
	log.info "--fastq_pattern      STR    Pattern for fastq file selection           Optional (.bam)"

	log.info "--qsize              INT    Max number of parallel jobs                Optional (20)"
	log.info "--outdir             DIR    Output directory                           Optional (.)"
	log.info "--logdir             DIR    Log directory                              Optional (\$outdir)"
	log.info "--help               FLAG   Displays this help"
	log.info ""
}

workflow.onComplete {
	if (workflow.success) {
		file(logdir + '/ribomethseq-nf.success').text = ''
	} else {
		def msg = "Launch directory : $workflow.launchDir\n"
		msg += "Work directory : $workflow.workDir\n"
		msg += "Exit status : $workflow.exitStatus\n"
		msg += "Error message : $workflow.errorMessage\n"
		msg += "Error report : $workflow.errorReport\n"
		file(logdir + '/ribomethseq-nf.failure').text = msg
	}
}
