#!/usr/bin/env nextflow

/*
#####Usage#####
#
#

*/

params.output_dir = "/home/raghu/ExonCNV/viscap/doc"
//params.bam = "/home/raghu/CNVK/bams/Normals/all/*bam"
params.bam = "/home/raghu/CNVK/bams/PSG/KCS/*bam"
params.bed = "/home/raghu/ExonCNV/decon/ref_data/PST_decon.bed"

Channel
    .fromPath(params.bam)
    .ifEmpty{error "Cannot find any files matching ${params.bam}"}
//    .view()
    .set{bam_ch}


process callIntragenicEvents {
    input:
    path bam_file from bam_ch

    output:
    stdout into result

    """
    gatk3 -T DepthOfCoverage -R /home/raghu/PubData/hg19.fa -o '${params.output_dir}/${bam_file.baseName}' -I '${bam_file}' -L '${params.bed}'
    """
}

result.subscribe{println it}
