#!/usr/bin/env nextflow 

date = new Date().format( 'yyyyMMdd' )

params.help                   = null
params.config                 = null
params.cpu                    = "12"
params.mem                    = "100"
params.tmpdir                 = "tmp/"
params.bulk_fq_sheet          = null
params.reference              = "/media/theboocock/Data/Dropbox/Postdoc/projects/crispr_coupling/data/tn5_96/saccer3_plk88.fasta"
params.oligo_table            = "/media/theboocock/Data/Dropbox/Postdoc/projects/crispr_coupling/data/tn5_96/otableLong.RDS"

log.info ""
log.info "------------------------------------------------------------------------"
log.info "  Nextflow pipeline to call variants from bulk population FASTQ files "
log.info "------------------------------------------------------------------------"
log.info ""


if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "  USAGE                                                 "
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "------------------ BSA with strelka2 -------------------"
    log.info ""
    log.info "nextflow run bulk_variant_calling.nf --bulk_fq_sheet sample_sheet.csv"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--ref                  FILE                 Genome reference file" 
    log.info "--bulk_fq_sheet        FOLDER               Folder containing bulk FASTQ files"
    log.info ""
    log.info "--------------------------------------------------------"
    log.info "Optional arguments:"
    log.info "--cpu                  INTEGER              Number of cpu to use (default=2)"
    log.info "--output_folder        PATH                 Output directory for vcf files (default=strelka_ouptut)"
    log.info "--config               FILE                 Use custom configuration file"
    log.info "--email                STRING               email address for job notifications"
    log.info ""
    log.info "Flags:"
    log.info "--help                                      Display this message"
    log.info ""
    exit 1
} else {

params.output_folder = "BulkGWA-${date}"

/* Software information */
log.info ""
log.info "ref                     = ${params.reference}"
log.info "cpu                     = ${params.cpu}"
log.info "mem                     = ${params.mem}Gb"
log.info "output_folder           = ${params.output_folder}"
log.info "bulk_fq_sheet           = ${params.bulk_fq_sheet}"
log.info ""
}

// Variant Filtering
params.min_depth = 5
params.qual = 30.0
params.strand_odds_ratio = 5.0
params.quality_by_depth = 2.0
params.fisherstrand = 100.0
params.high_missing = 0.95



/*
    ===============
    Iniitalize input parameters for pipeline
    ===============

*/

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* ~~~~~~~~~~~~~~~~~~~ Read in bulk FASTQ information ~~~~~~~~~~~~~~~~~~ */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


workflow {

    File fq_file = new File(params.bulk_fq_sheet)

    in_fastas = Channel.from(fq_file.collect { it.tokenize( ',' ) })
                 .map { SM, ID, fq1, fq2 -> [SM, ID, file("${fq1}"), file("${fq2}")] }
    File reference = new File("${params.reference}")
    params.reference_handle = reference.getAbsolutePath()
 /*   plasmid_spades (in_fastas) */
    trim_reads(in_fastas)
    perform_bsa_alignment(trim_reads.out) | 
        tag_dups
    assemble_plasmid_and_filter(tag_dups.out.outs)
     
         
    scer_chrom=Channel
        .from("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")
    chrom_combos = scer_chrom.combine(assemble_plasmid_and_filter.out.outs)

    
    /*chrom_combos.view()*/
/*    call_variants_individual(chrom_combos) */

    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    /* ~~~~~~~~~~~~~~~~~~~ Initialize reference sequence ~~~~~~~~~~~~~~~~~~~ */
    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  /*  File reference = new File("${params.reference}")
    reference_handle = reference.getAbsolutePath()
    sample_map = fqs_to.map { "${it[0]}\t${it[0]}.g.vcf.gz" }.collectFile(name: "sample_map.tsv", newLine: true)
    sample_map.view()
    bulk_gvcf
        .groupTuple()
        .into{g_to_import; g_to_nothing}

    ce_chrom=Channel
         .from("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")
         .combine(strain_g)
         .groupTuple()
         .combine(sample_map)
         .into{strain_g_to_import; strain_g_to_nothing}
    Channel.from("BULK_POPS")
        .combine(cohort_chrom_vcf)
        .groupTuple()
        .into{cohort_vcf_to_concat; cohort_vcf_to_nothing}
*/

}
/* 
=====================================================================
=================================================================
                                                                   ===============
                                                                    ======================= BEGIN PIPELINE
                                                                   ===============
=================================================================
=====================================================================
*/ 

/*
   ce_chrom=Channel
     .from("I","II","III","IV","V","X","MtDNA")
     .combine(fqs_align)
     .println()
*/



/*
    ====================
    Bulk Fastq Alignment
    ====================
*/

process trim_reads {
    cpus 1
    memory '12 GB'
    tag { ID }

    input:
        tuple val(SM), val(ID), path(fq1), path(fq2)
        
    output: 
        tuple val(SM), val(ID), path("${SM}.r1.filt.fq.gz"), path("${SM}.r2.filt.fq.gz")

    """
        trimmomatic PE $fq1 $fq2 ${SM}.r1.filt.fq.gz ${SM}.r1.unpaired.fq.gz \\
        ${SM}.r2.filt.fq.gz ${SM}.r2.unpaired.fq.gz ILLUMINACLIP:NexteraPE-PE.fa:2:30:10:2:True \\
        LEADING:3 TRAILING:3 MINLEN:36
    """


}

process plasmid_spades {
    cpus 4
    memory '12 GB'
    tag { ID }

    input:
        tuple val(SM), val(ID), path(fq1), path(fq2)
        
    output: 
        tuple val(SM)
    """
        plasmidspades.py -1 ${fq1} -2 ${fq2} -o . 
        cp scaffold.fasta ${SM}.fasta 
    """

}

process perform_bsa_alignment {

    cpus 4
    memory '12 GB'
    tag { ID }

    input:
        tuple val(SM), val(ID), path(fq1), path(fq2) 

    output:
        tuple val(SM), path("${ID}.bam"), path("${ID}.bam.bai") 

    """
    bwa mem -t ${task.cpus} -R '@RG\\tID:${ID}\\tSM:${SM}\\tLB:${SM}_${ID}\\tPL:ILLUMINA' ${params.reference_handle} ${fq1} ${fq2} | \\
    sambamba view --nthreads=${task.cpus} --show-progress --sam-input --format=bam --with-header /dev/stdin | \\
    sambamba sort --nthreads=${task.cpus} --show-progress --tmpdir=${params.tmpdir} --out=${ID}.bam /dev/stdin
    sambamba index --nthreads=${task.cpus} ${ID}.bam

    if [[ ! \$(samtools view ${ID}.bam | head -n 10) ]]; then
        exit 1;
    fi
    """
}



process tag_dups {

    cpus 4
    memory '12 GB'

    publishDir "${params.output_folder}/alignments", mode: "copy"

    tag { SM }

    input:
        tuple val(SM), path(bam), path(bami) 

    output:
        tuple val(SM), path("${SM}.bam"), path("${SM}.bam.bai"), emit: outs
        path("${SM}.picard.sam.markduplicates") 

    """
    picard MarkDuplicates I=${bam} \\
        O=${SM}.bam \\
        M=${SM}.picard.sam.markduplicates \\
        VALIDATION_STRINGENCY=SILENT \\
        REMOVE_DUPLICATES=false \\
        TAGGING_POLICY=All \\
        REMOVE_SEQUENCING_DUPLICATES=TRUE \\
        SORTING_COLLECTION_SIZE_RATIO=0.1
        
    sambamba index --nthreads=${task.cpus} ${SM}.bam
    """
}

process assemble_plasmid_and_filter {
    
    cpus 4
    memory '12 GB'
    cache false
    publishDir "${params.output_folder}/filtered_bams_and_edits", mode: "copy"
    tag {SM}
    input:
        tuple val(SM), file(bam), file(bami) 
    output:
        path("${SM}_stats.txt"), emit: stats 
        tuple val(SM), path("${SM}_final.bam"), path("${SM}_final.bam.bai"), emit: outs, optional: true
        tuple path("${SM}_plasmid.txt"), path("${SM}_blast.txt"),path("${SM}_combo_match.txt"), optional: true
    """
        python ${baseDir}/scripts/spades_assembly_plk88.py \\
        -o ${SM} \\
        --oligo-table ${params.oligo_table} \\
        ${bam}
    """
}

/* 
   =================================
    Call Variants in Bulk Experiment
   =================================


/*=========================================
~ ~ ~ > * Generate Interval List  * < ~ ~ ~ 
=========================================*/


/*
===========================================
~ ~ ~ > * Run GATK4 HaplotypeCaller * < ~ ~ ~ 
===========================================
*/

process call_variants_individual {

    publishDir "${params.output_folder}/vcfs", mode: "copy"

    cpus params.cpu
    memory '40 GB'

    clusterOptions = '-V -l highp,h_data=10G'

    input:
        tuple val(CHROM), val(SM), path(bam), path(bambai) 

    output:
        tuple val(SM), path("${CHROM}_bulk.g.vcf.gz"), path("${CHROM}_bulk.g.vcf.gz.csi") 

    """
        gatk HaplotypeCaller --java-options "-Xmx${task.memory.toGiga()}g -Xms1g -XX:ConcGCThreads=${task.cpus}" \\
            --emit-ref-confidence GVCF \\
            --annotation DepthPerAlleleBySample \\
            --annotation Coverage \\
            --annotation GenotypeSummaries \\
            --annotation TandemRepeat \\
            --annotation StrandBiasBySample \\
            --annotation ChromosomeCounts \\
            --annotation ReadPosRankSumTest \\
            --annotation AS_ReadPosRankSumTest \\
            --annotation AS_QualByDepth \\
            --annotation AS_StrandOddsRatio \\
            --annotation AS_MappingQualityRankSumTest \\
            --annotation DepthPerSampleHC \\
            --annotation-group StandardAnnotation \\
            --annotation-group AS_StandardAnnotation \\
            --annotation-group StandardHCAnnotation \\
            --do-not-run-physical-phasing \\
            -R ${params.reference} \\
            -I ${bam} \\
            -L ${CHROM} \\
            -O ${SM}_${CHROM}_bulk.g.vcf   
        bcftools view -O z ${SM}_${CHROM}_bulk.g.vcf > ${SM}_${CHROM}_bulk.g.vcf.gz
        bcftools index ${SM}_${CHROM}_bulk.g.vcf.gz
        rm ${SM}_${CHROM}_bulk.g.vcf
    """
}


/*
=============================================
~ ~ ~ > *  Merge Sample gVCF files  * < ~ ~ ~ 
=============================================
*/

process concat_strain_gvcfs {

    tag { "${SM}" }

    input:
        set val(SM), file("*_bulk.g.vcf.gz"), file("*_bulk.g.vcf.gz.csi") from g_to_import

    output:
        set file("${SM}.g.vcf.gz"), file("${SM}.g.vcf.gz.tbi") into strain_g

    """
        ls | grep "bulk.g.vcf.gz\$" > contig_set.tsv
        bcftools concat  -O z --file-list contig_set.tsv |\\
        bcftools sort -Oz > ${SM}.g.vcf.gz
        bcftools index --tbi ${SM}.g.vcf.gz
    """
}


 process import_genomics_db {

    tag { "${chrom}" }
    
    cpus params.cpu
    memory '40 GB'

    clusterOptions = '-V -l highp,h_data=10G'

    input:
        set val(chrom), file(vcf), file(vcfi), path(sample_map) from strain_g_to_import

    output:
        set val(chrom), file("${chrom}.db") into strain_gdb

    """
        gatk  --java-options "-Xmx${task.memory.toGiga()-3}g -Xms${task.memory.toGiga()-4}g -XX:ConcGCThreads=${task.cpus}" \\
            GenomicsDBImport \\
            --genomicsdb-workspace-path ${chrom}.db \\
            --batch-size 16 \\
            -L ${chrom} \\
            --sample-name-map ${sample_map} \\
            --reader-threads ${task.cpus}
    """
}

/*==================================
~ ~ ~ > *  Genotype gVCFs  * < ~ ~ ~ 
==================================*/

process genotype_cohort_gvcf_db {
    // Heterozygous polarization is also performed here.

    publishDir "${params.output_folder}/vcfs", mode: "copy"

    tag { "${contig}" }

    cpus params.cpu
    memory '40 GB'

    clusterOptions = '-V -l highp,h_data=10G'

    input:
        set val(contig), file("${contig}.db") from strain_gdb

    output:
        set file("${contig}_cohort.vcf.gz"), file("${contig}_cohort.vcf.gz.csi") into cohort_chrom_vcf


    /*
        het_polarization polarizes het-variants to REF or ALT (except for mitochondria)
    */

    """
        gatk  --java-options "-Xmx${task.memory.toGiga()}g -Xms1g" \\
            GenotypeGVCFs \\
            -R ${reference_handle} \\
            -V gendb://${contig}.db \\
            -G StandardAnnotation \\
            -G AS_StandardAnnotation \\
            -G StandardHCAnnotation \\
            -L ${contig} \\
            -O ${contig}_cohort.vcf

        bcftools view -O z --threads=${task.cpus-1} ${contig}_cohort.vcf -o ${contig}_cohort.vcf.gz
        bcftools index ${contig}_cohort.vcf.gz
        
    """
}

/*===============================================
~ ~ ~ > *   Concatenate VCFs  * < ~ ~ ~
===============================================*/

process concatenate_vcf {

    publishDir "${params.output_folder}/vcfs", mode: "copy"

    cpus params.cpu
    memory '16 GB'

    clusterOptions = '-V -l highp,h_data=4G'

    input: 
        set val(pops), file("*_cohort.vcf.gz"), file("*_cohort.vcf.gz.csi") from cohort_vcf_to_concat

    output:
        set file("WI.raw.vcf.gz"), file("WI.raw.vcf.gz.tbi") into cohort_concat_vcf

    """
        ls *_cohort.vcf.gz > contig_set.tsv
        bcftools concat  -O z --file-list contig_set.tsv |\\
        bcftools sort -Oz > WI.raw.vcf.gz
        bcftools index --tbi WI.raw.vcf.gz
    """
}

/*===========================================
~ ~ ~ > *   Apply SNV Soft Filters  * < ~ ~ ~
===========================================*/

process soft_filter {

    publishDir "${params.output_folder}/variation", mode: 'copy'

    cpus params.cpu
    memory '16 GB'

    clusterOptions = '-V -l highp,h_data=4G'

    input:
        set file(vcf), file(index) from cohort_concat_vcf

    output:
        set file("WI.${date}.soft-filter.vcf.gz"), file("WI.${date}.soft-filter.vcf.gz.csi"), file("WI.${date}.soft-filter.tsv") into soft_vcf
    

    """
        function cleanup {
            rm out.vcf.gz
        }
        trap cleanup EXIT
        gatk --java-options "-Xmx${task.memory.toGiga()}g -Xms1g" \\
            VariantFiltration \\
            -R ${reference_handle} \\
            --variant ${vcf} \\
            --genotype-filter-expression "DP < ${params.min_depth}"    --genotype-filter-name "DP_min_depth" \\
            --filter-expression "QUAL < ${params.qual}"                --filter-name "QUAL_quality" \\
            --filter-expression "FS > ${params.fisherstrand}"          --filter-name "FS_fisherstrand" \\
            --filter-expression "QD < ${params.quality_by_depth}"      --filter-name "QD_quality_by_depth" \\
            --filter-expression "SOR > ${params.strand_odds_ratio}"    --filter-name "SOR_strand_odds_ratio" \\
            -O out.vcf
        
        bgzip out.vcf
        bcftools index --tbi out.vcf.gz
        
        # Get Stats
        bcftools view out.vcf.gz -Oz -o WI.${date}.soft-filter.vcf.gz
        bcftools index WI.${date}.soft-filter.vcf.gz
        bcftools index --tbi WI.${date}.soft-filter.vcf.gz

        bcftools query -f '[%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER\\t%SAMPLE\\t%GT\\t%AD\\t%DP\\t%FT\\t%GQ\\n]' WI.${date}.soft-filter.vcf.gz > WI.${date}.soft-filter.tsv
    """
}

/*
   ====================================
    Extract ALT counts for all variants
   ====================================
*/ 

/* 
   ================================================
    Run Josh's code to calculate strain frequencies
   ================================================
*/ 

/* 
   =====================================
    Summary plots of strain frequencies
   =====================================
*/ 
