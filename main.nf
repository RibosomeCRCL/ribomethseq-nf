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
params.fastq_pattern = '.fastq.gz'

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

// samtools
params.samtools_opts = "--no-PG -h -u -d 'NM' -e '![XS]'"

params.split         = false
params.help          = false

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
log.info "Bowtie index                     : ${params.bowtie_index}"
log.info "samtools options                 : ${params.samtools_opts}"
log.info "Threads for bowtie               : ${params.bowtie_threads}"
log.info "Threads for fastqc               : ${params.fastqc_threads}"
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
	.into { reads_ch; reads_ch2; reads_ch3 }


process fastqc {
	tag { sample_id }

	publishDir "${outdir}/fastqc", mode: 'copy', pattern: '*_fastqc.html'
	publishDir "${outdir}/fastqc/zip", mode: 'copy', pattern: '*_fastqc.zip'

	input:
	set sample_id, file(reads) from reads_ch

	output:
	file "*_fastqc.html"
	file "*_fastqc.zip" into fastqc_ch, fastqc_ch2

	"""
	fastqc --threads ${params.fastqc_threads} \\
	--outdir . ${reads}
	"""
}

process trim {
	tag { sample_id }

	publishDir "${outdir}/trimmomatic", mode: 'copy', pattern : "${sample_id}.trim.fastq.gz"
	publishDir "${outdir}/trimmomatic/logs", mode: 'copy', pattern : "${sample_id}.trimmomatic.stats.log"

	input:
	set sample_id, file(reads) from reads_ch2

	output:
	set sample_id, file("${sample_id}.trim.fastq.gz") into trimmed_files, trimmed_files2
	file("${sample_id}.trimmomatic.stats.log") into trimmomatic_logs, trimmomatic_logs2

	"""
	trimmomatic SE -phred33 -threads ${params.trimmo_threads} \\
	${reads} ${sample_id}.trim.fastq.gz \\
	ILLUMINACLIP:${params.adapters}:2:30:10 \\
	LEADING:${params.leading} \\
	TRAILING:${params.trailing} \\
	SLIDINGWINDOW:${params.slidingwindow} \\
	AVGQUAL:${params.avgqual} \\
	MINLEN:${params.minlen} 2> ${sample_id}.trimmomatic.stats.log
	"""
}

process bowtie2 {
	tag { sample_id }

	publishDir "${outdir}/bowtie2/logs", mode: 'copy', pattern : "${sample_id}.bowtie2.stats.log"

	input:
	set sample_id, file(reads) from trimmed_files

	output:
	set sample_id, file("${sample_id}.sam") into bowtie_files
	file("${sample_id}.bowtie2.stats.log") into bowtie_logs, bowtie_logs2

	"""
	bowtie2 -x ${params.bowtie_index} --threads ${params.bowtie_threads} \\
	${params.bowtie_opts} ${reads} -S ${sample_id}.sam 2>> ${sample_id}.bowtie2.stats.log
	"""
}

process filter {
	tag { sample_id }

	publishDir "${outdir}/bowtie2", mode: 'copy', pattern : "${sample_id}.uniq.{bam,bam.bai}"

	input:
	set sample_id, file(sam) from bowtie_files

	output:
	set sample_id, file("${sample_id}.uniq.bam") into filtered_bam_ch

	"""
	set -o pipefail
	samtools view ${params.samtools_opts} ${sam} | \\
	samtools sort -@${params.samtools_threads} -o ${sample_id}.uniq.bam
	samtools index ${sample_id}.uniq.bam -o ${sample_id}.uniq.bam.bai
	"""
}

process multiqc {

	publishDir "${outdir}", mode: 'copy', pattern: 'multiqc_report.html'

	input:
	file('*') from fastqc_ch.collect()
	file('*') from trimmomatic_logs.collect()
	file('*') from bowtie_logs.collect()

	output:
	file('multiqc_report.html')

	"""
	multiqc .
	"""
}

process counts {
	tag { sample_id }

	publishDir "${outdir}/counts", mode:'copy', pattern: "${sample_id}.*_counts.csv"

	input:
	set sample_id, file(filtered_bam) from filtered_bam_ch

	output:
	set sample_id, file("${sample_id}.5_counts.csv"), file("${sample_id}.3_counts.csv") into counts_ch
	set sample_id, file("${sample_id}.5_counts.csv") into counts_ch2
	script:
	"""
	bedtools genomecov -d -3 -ibam ${filtered_bam} > ${sample_id}.3_counts.csv
	bedtools genomecov -d -5 -ibam ${filtered_bam} > ${sample_id}.5_counts.csv
	"""
}

process split {
	tag { sample_id }

	publishDir "${outdir}/split", mode: 'move', pattern: "${sample_id}.*.csv"

	when: params.split

	input:
	set sample_id, file(counts5), file(counts3) from counts_ch

	output:
	file("${sample_id}.58S.csv")
	file("${sample_id}.18S.csv")
	file("${sample_id}.28S.csv")
	file("${sample_id}.5S.csv")

	script:
	"""
	Rscript ${baseDir}/Rscripts/Refine/r_refine.R ${counts5} ${counts3} \\
	${params.fasta_58S} ${params.fasta_18S} ${params.fasta_28S} ${params.fasta_5S} \\
	${sample_id}
	"""
}

process report {
	publishDir "${outdir}", mode: 'move'

	input:
	file('*') from counts_ch2.collect()

	output:
	file("rms_report.html")

	script:
	"""
	Rscript ${baseDir}/Rscripts/QC/generate_report.R .
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
	log.info "                           +++++++++++++++++++++++++"
	log.info "                           +  ribomethseq-nf help  +"
	log.info "                           +++++++++++++++++++++++++"
	log.info ""
	log.info "--fqdir              DIR    Fastq files location                       Required"
	log.info "--fastq_pattern      STR    Pattern for fastq file selection           Optional (.fastq.gz)"
	log.info ""
	log.info "--adapters           FILE   (Trimmomatic) Path to illumina adapters    Optional (\$baseDir/data/adapters/TruSeq3-SE.fa)"
	log.info "--leading            INT    (Trimmomatic) LEADING parameter            Optional (30)"
	log.info "--trailing           INT    (Trimmomatic) TRAILING parameter           Optional (30)"
	log.info "--slidingwindow      STR    (Trimmomatic) SLIDINGWINDOW parameter      Optional (4:15)"
	log.info "--avgqual            INT    (Trimmomatic) AVGQUAL parameter            Optional (30)"
	log.info "--minlen             INT    (Trimmomatic) MINLEN parameter             Optional (8)"
	log.info ""
	log.info "--bowtie_index       FILE   (Bowtie) Path to index                     Optional (\$baseDir/data/bowtie/human/human_index)"
	log.info "--bowtie_opts        STR    (Bowtie) additional options                Optional (--sensitive -L 17)"
	log.info ""
	log.info "--samtools_opts      STR    (samtools) options to view                 Optional (--no-PG -h -u -d 'NM' -e '![XS]')"
	log.info ""
	log.info "--bowtie_threads     INT    Threads for bowtie                         Optional (7)"
	log.info "--fastqc_threads     INT    Threads for fastqc                         Optional (2)"
	log.info "--trimmo_threads     INT    Threads for trimmomatic                    Optional (3)"
	log.info "--samtools_threads   INT    Threads for samtools                       Optional (4)"
	log.info ""
	log.info "--split              FLAG   Split count files by RNA                   Optional (false)"
	log.info "--scheduler          STR    Job scheduler                              Optional (slurm)"
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
