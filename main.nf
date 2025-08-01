#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameters
params.genome = params.genome ?: 'hg38'  // Default to hg38 if not specified
params.platform = params.platform ?: 'illumina_pe'  // Default to Illumina paired-end
params.input_type = null
params.vcf_type = null
params.refhg38 = '/usr/src/app/ref38/hg38/hg381_22XYM/Homo_sapiens_assembly38cleaned.fasta' //chr1_22 X_Y_M only
params.refhg37 = '/usr/src/app/ref38/hg19/hg19122XYM/hg19122XYM.fa' //chr1_22 X_Y_M only
params.bedhg38 = params.bedhg38 ?:'/usr/src/app/ref38/BED/hg38_exome.bed'
params.bedhg37 = params.bedhg37 ?:'/usr/src/app/ref38/BED/hg37_exome.bed'
params.cnv_reference38 = '/usr/src/app/reference/reference/cnv_cnvkit/flat_referencehg38.cnn'
params.cnv_reference37 = '/usr/src/app/reference/reference/cnv_cnvkit/flat_referencehg19_cleaned.cnn'
params.mito_reference = '/usr/src/app/reference/reference/mitochondria/rCRS.fasta'
params.mito_annotation = '/usr/src/app/reference/reference/mitochondria/rCRS_annotation_2020-08-20.txt'
params.vusprize_script = '/usr/src/app/VusPrize/vusprize/VusPrize.py'
params.dependent_script = '/usr/src/app/VusPrize/vusprize/RF_niu.joblib'

//new
params.secondaryfindings_gene_file = '/usr/src/app/modules/ACMG_SF3_gene_list.xlsx'
params.clingen_file38 = '/usr/src/app/modules/ClinGen_gene_curation_list_GRCh38_cleaned.tsv'
params.clingen_file37 = '/usr/src/app/modules/ClinGen_gene_curation_list_GRCh37_cleaned.tsv'
params.gene_data_file = '/usr/src/app/modules/gene_data_with_lof_class.csv' //taken from clingen gene-disease association
//params.gender = params.gender ?: 'male'
params.moi_dk = '/usr/src/app/modules/gene.tsv'
params.omim = '/usr/src/app/modules/merge_9314_final_formatted.csv'

//// VEP input Databases Hg38
params.cache_dir = '/usr/src/app/vepC/cache'
params.dirPlugin = '/usr/src/app/ensembl-vep/Plugins'
params.dbnsfp38 = '/usr/src/app/vepDB/DBs/hg38/dbNSFP5.1a_grch38-004.gz'
params.dbnsfp4_38 = '/usr/src/app/vepDB/DBs/hg38/dbNSFP4.7a_grch38.gz'
params.LoFtool = '/usr/src/app/vepPlugins/LoFtool_scores.txt'
params.CADDsnv38 = '/usr/src/app/vepDB/DBs/hg38/whole_genome_SNVs.tsv.gz'
params.CADDindels38 = '/usr/src/app/vepDB/DBs/hg38/gnomad.genomes.r4.0.indel_inclAnno.tsv.gz'
params.dbscSNV38 = '/usr/src/app/vepDB/DBs/hg38/dbscSNV1.1_GRCh38.txt.gz'
params.dosage_collins = '/usr/src/app/vepPlugins/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz'
params.spliceaiindel38 = '/usr/src/app/vepDB/DBs/hg38/spliceai_scores.raw.indel.hg38.vcf.gz'
params.spliceaisnv38 = '/usr/src/app/vepDB/DBs/hg38/spliceai_scores.raw.snv.hg38.vcf.gz'
params.pli38 = '/usr/src/app/vepDB/DBs/hg38/plI_gene.txt'
params.phenotypesfile38 = '/usr/src/app/vepDB/DBs/hg38/Phenotypes.pm_homo_sapiens_114_GRCh38.gvf.gz'
params.loeuf38 = '/usr/src/app/vepDB/DBs/hg38/loeuf_dataset_grch38.tsv.gz'
params.clinpred38 = '/usr/src/app/vepDB/DBs/hg38/ClinPred_hg38_sorted_tabbed.tsv.gz'
params.alphamissense38 = '/usr/src/app/vepDB/DBs/hg38/AlphaMissense_hg38.tsv.gz'
params.primateai38 = '/usr/src/app/vepDB/DBs/hg38/PrimateAI_scores_v0.2_GRCh38_sorted.tsv.bgz'
params.revel38 = '/usr/src/app/vepDB/DBs/hg38/new_tabbed_revel_grch38.tsv.gz'
params.clinvar38 = '/usr/src/app/vepDB/DBs/hg38/clinvar.vcf.gz'

//// VEP input Databases Hg37
params.dbnsfp37 ='/usr/src/app/vepDB/DBs/hg19/dbNSFP5.1a_grch37-003.gz'
params.dbnsfp4_37 = '/usr/src/app/vepDB/DBs/hg19/dbNSFP4.7a_grch37.gz'
params.CADDsnv37 = '/usr/src/app/vepDB/DBs/hg19/whole_genome_SNVs_hg37.tsv.gz'
params.CADDindels37 = '/usr/src/app/vepDB/DBs/hg19/gnomad.genomes-exomes.r4.0.indel_inclAnno_hg37.tsv.gz'
params.dbscSNV37 = '/usr/src/app/vepDB/DBs/hg19/dbscSNV1.1_GRCh37.txt.gz'
params.LoFtool = '/usr/src/app/vepPlugins/LoFtool_scores.txt'
params.CADDsnv37 = '/usr/src/app/vepDB/DBs/hg19/whole_genome_SNVs_hg37.tsv.gz'
params.CADDindels37 = '/usr/src/app/vepDB/DBs/hg19/gnomad.genomes.r4.0.indel_inclAnno_hg37.tsv.gz'
params.dbscSNV37 = '/usr/src/app/vepDB/DBs/hg19/dbscSNV1.1_GRCh37.txt.gz'
params.dosage_collins = '/usr/src/app/vepPlugins/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz'
params.spliceaiindel37 = '/usr/src/app/vepDB/DBs/hg19/spliceai_scores.raw.indel.hg19.vcf.gz'
params.spliceaisnv37 = '/usr/src/app/vepDB/DBs/hg19/spliceai_scores.raw.snv.hg19.vcf.gz'
params.pli37 = '/usr/src/app/vepDB/DBs/hg19/plI_gene.txt'
params.phenotypesfile37 = '/usr/src/app/vepDB/DBs/hg19/Phenotypes.pm_homo_sapiens_114_GRCh37.gvf.gz'
params.loeuf37 = '/usr/src/app/vepDB/DBs/hg19/loeuf_dataset_hg37.tsv.gz'
params.clinpred37 = '/usr/src/app/vepDB/DBs/hg19/ClinPred_tabbed.tsv.gz'
params.alphamissense37 = '/usr/src/app/vepDB/DBs/hg19/AlphaMissense_hg19.tsv.gz'
params.primateai37 = '/usr/src/app/vepDB/DBs/hg19/PrimateAI_scores_v0.2_hg19.tsv.bgz'
params.revel37 = '/usr/src/app/vepDB/DBs/hg19/new_tabbed_revel_grch37.tsv.gz'
params.clinvar37 = '/usr/src/app/vepDB/DBs/hg19/clinvar.vcf.gz'

/// Exomiser parameters
params.template_yml = params.template_yml ?:'/usr/src/app/exomiser/db/exomiser-template.yml'  // Path inside the container
params.phenotypes = []                        // List of phenotypes
params.exomiser_jar = '/usr/src/app/exomiser/exomiser-cli-14.0.0/exomiser-cli-14.0.0.jar'   // Correct jar path in Docker
params.exomiser_hg38_properties = '/usr/src/app/vepDB/DBs/exomiser/hg38application.properties' // Properties path
params.exomiser_hg19_properties = '/usr/src/app/vepDB/DBs/exomiser/hg37application.properties'

/// CNV baseline sample dir
params.cnv_baseline = '/usr/src/app/ref38/CNV_baseline_samples/hg38/'


    // ACMG
params.cache_biashg38 = '/usr/src/app/ref38/ACMG/data/GRCh38/Cache/GRCh38/Both'
params.sd_biashg38 = '/usr/src/app/ref38/ACMG/data/GRCh38/SupplementaryAnnotation/GRCh38/'
params.ref_biashg38 = '/usr/src/app/ref38/ACMG/data/GRCh38/References/Homo_sapiens.GRCh38.Nirvana.dat'
params.parameter_biashg38 = '/usr/src/app/ref38/ACMG/hg38_required_paths.json'
    
params.cache_biashg19 = '/usr/src/app/ref38/ACMG/data/GRCh37/Cache/GRCh37/Both'
params.sd_biashg19 = '/usr/src/app/ref38/ACMG/data/GRCh37/SupplementaryAnnotation/GRCh37/'
params.ref_biashg19 = '/usr/src/app/ref38/ACMG/data/GRCh37/References/Homo_sapiens.GRCh37.Nirvana.dat'
params.parameter_biashg19 = '/usr/src/app/ref38/ACMG/hg19_required_paths.json'

// Genome reference selection
if (params.genome == 'hg38') {
    params.ref = params.refhg38
    params.ref_dir = params.ref_dirhg38
    params.reffai = params.refhg38fai
    params.bed = params.bedhg38
    params.gnomad = params.gnomad38
    params.db1000g = params.db1000g38
    params.msimodel = params.msimodelhg38
    params.assembly = 'GRCh38'
    params.cachee = params.cache_dir
    params.dirplugin = params.dirPlugin
    params.dbNSFP = params.dbnsfp38
    params.loftool = params.LoFtool
    params.CADDsnv = params.CADDsnv38
    params.CADDindel = params.CADDindels38
    params.dbscSNV = params.dbscSNV38
    params.cnv_reference = params.cnv_reference38
    params.mito_reference = params.mito_reference
    params.mito_annotation = params.mito_annotation
    params.exomiser_properties = params.exomiser_hg38_properties
    params.template_yml = params.template_yml
    params.vusprize_script = params.vusprize_script
    params.dependent_script = params.dependent_script
    params.spliceaiindel = params.spliceaiindel38
    params.spliceaisnv = params.spliceaisnv38
    params.phenotypesfile = params.phenotypesfile38
    params.loeuf = params.loeuf38
    params.clinpred = params.clinpred38
    params.alphamissense = params.alphamissense38
    params.primateai = params.primateai38
    params.revel = params.revel38
    params.clinvar = params.clinvar38
    params.clingen_file = params.clingen_file38
    params.pli = params.pli38
    params.dbnsfp4 = params.dbnsfp4_38
    // ACMG
    params.cache_bias = params.cache_biashg38
    params.sd_bias = params.sd_biashg38
    params.ref_bias = params.ref_biashg38
    params.parameter_bias = params.parameter_biashg38

} else if (params.genome == 'hg19') {
    params.ref = params.refhg37
    params.ref_dir = params.ref_dirhg37
    params.reffai = params.refhg37fai
    params.bed = params.bedhg37
    params.gnomad = params.gnomad37
    params.db1000g = params.db1000g37
    params.msimodel = params.msimodelhg19
    params.assembly = 'GRCh37'
    params.cachee = params.cache_dir
    params.dirplugin = params.dirPlugin
    params.dbNSFP = params.dbnsfp37
    params.loftool = params.LoFtool
    params.CADDsnv = params.CADDsnv37
    params.CADDindel = params.CADDindels37
    params.dbscSNV = params.dbscSNV37
    params.cnv_reference = params.cnv_reference37
    params.mito_reference = params.mito_reference
    params.mito_annotation = params.mito_annotation
    params.exomiser_properties = params.exomiser_hg19_properties
    params.template_yml = params.template_yml
    params.clingen_file = params.clingen_file37
    params.pli = params.pli37
    params.vusprize_script = params.vusprize_script
    params.dependent_script = params.dependent_script
    params.spliceaiindel = params.spliceaiindel37
    params.spliceaisnv = params.spliceaisnv37
    params.phenotypesfile = params.phenotypesfile37
    params.loeuf = params.loeuf37
    params.clinpred = params.clinpred37
    params.alphamissense = params.alphamissense37
    params.primateai = params.primateai37
    params.revel = params.revel37
    params.clinvar = params.clinvar37
    params.dbnsfp4 = params.dbnsfp4_37
    // ACMG
    params.cache_bias = params.cache_biashg19
    params.sd_bias = params.sd_biashg19
    params.ref_bias = params.ref_biashg19
    params.parameter_bias = params.parameter_biashg19

} else {
    error "Unsupported genome reference: ${params.genome}. Please use 'hg38' or 'hg19'."
}

params.maphg38 = '/usr/src/app/mapfileDelly/Hg38.map'
params.maphg37 = '/usr/src/app/mapfileDelly/'

// Module imports
include { fastqc_reads } from './modules/1.fastqc_module.nf'
include { trim_with_fastp_paired; trim_with_fastp_single } from './modules/2.fastp_module.nf'
include { align_bwa_paired; align_bwa_single } from './modules/3.bwa_module.nf'
include { markDup } from './modules/4.picard_module.nf'
include { samtool_index } from './modules/5.bamIndex_module.nf'
include { germVarCall } from './modules/6.Deepvariant_module.nf'
include { KeepPASS } from './modules/8.keeppass.nf'
include { annotVEP } from './modules/9.VEPannotation.nf'
include { vusVEP } from './modules/10.process_vus.nf'
include { processVEPOutput; runVusPrioritization } from './modules/11.process_vus_2.nf'
include { MODIFY_CONFIG; run_prioritizer } from './modules/16.exomiser_worked.nf'
include { SV_germline } from './modules/13.SVgermline.nf'
include { CNV_germline } from './modules/14.CNV.nf'
include { mito_variant } from './modules/15.mito_variant.nf'
include { acmgpre; acmgclass } from './modules/acmg.nf'
include { vcfvep_process } from './modules/vcfvep_process.nf'
include { mergeExomiserVEPResults } from './modules/merge_exomiser_vep_results.nf'
include { runSecondaryFindings } from './modules/secondary_findings.nf'
include { lofDetection } from './modules/lof_detection.nf'
include { runG2P } from './modules/G2P_filteration.nf'

// Helper function
def hasPhenotypes() {
    return !params.phenotypes.isEmpty()
}

/*
 * PREPROCESSING WORKFLOW
 * Handles input processing, QC, alignment, and variant calling
 */
workflow PREPROCESSING {
    take:
    input_ch

    main:
    // Initialize channels
    indexed_bams = Channel.empty()
    qc_reports = Channel.empty()
    snv_ch = Channel.empty()
    sv_ch = Channel.empty()
    cnv_ch = Channel.empty()
    mito_ch = Channel.empty()
    final_p2g_results = Channel.empty()    

    // Skip flags
    boolean SKIP_SNV = params.skip_snv ?: false
    boolean SKIP_SV = params.skip_sv ?: false
    boolean SKIP_CNV = params.skip_cnv ?: false
    boolean SKIP_MITO = params.skip_mito ?: false
    String VCF_TYPE = params.input_type == 'vcf' ? params.vcf_type?.toLowerCase() : null

    // Handle different input types
    if (params.input_type == 'vcf') {
        // Handle VCF input based on type
        if (VCF_TYPE == 'snv' && !SKIP_SNV) {
            snv_ch = input_ch
        } else if (VCF_TYPE == 'sv' && !SKIP_SV) {
            sv_ch = input_ch
        } else if (VCF_TYPE == 'cnv' && !SKIP_CNV) {
            cnv_ch = input_ch
        } else if (VCF_TYPE == 'mito' && !SKIP_MITO) {
            mito_ch = input_ch
        } else if (!VCF_TYPE) {
            // Assume mixed VCF, assign to SNV by default
            snv_ch = input_ch
        }
        indexed_bams = Channel.empty()
        
    } else if (params.input_type == 'bam') {
        // Index BAM files and call variants
        indexed_bams = samtool_index(input_ch)
        
        // Variant calling from BAM files
        if (!SKIP_SNV) {
            germVarCall(indexed_bams, params.ref, params.bed)
            raw_snv = germVarCall.out.raw_vcfs
            KeepPASS(raw_snv)
            snv_ch = KeepPASS.out.final_vcfs
        }
        
        if (!SKIP_SV) {
            SV_germline(indexed_bams, params.ref)
            sv_ch = SV_germline.out.sv_vcf
        }
        
        if (!SKIP_CNV) {
            CNV_germline(indexed_bams, params.cnv_baseline)
            cnv_ch = CNV_germline.out.cnv_vcf
        }
        
        if (!SKIP_MITO) {
            mito_variant(indexed_bams, params.mito_reference, params.mito_annotation)
            mito_ch = mito_variant.out.mito_vcf
        }
        
    } else { // fastq
        // Full preprocessing pipeline
        fastqc_reads(input_ch)
        qc_reports = fastqc_reads.out.fastqc_out
        
        // Handle paired-end vs single-end processing
        if (params.platform.endsWith('_pe')) {
            // Paired-end processing
            trim_with_fastp_paired(input_ch)
            trimmed_reads = trim_with_fastp_paired.out.fastp_out_pe
            
            align_bwa_paired(trimmed_reads, params.ref)
            aligned_bams = align_bwa_paired.out.sorted_bams_pe
        } else {
            // Single-end processing
            trim_with_fastp_single(input_ch)
            trimmed_reads = trim_with_fastp_single.out.fastp_out_se
            
            align_bwa_single(trimmed_reads, params.ref)
            aligned_bams = align_bwa_single.out.sorted_bams_se
        }
        
        markDup(aligned_bams)
        markdup_bams = markDup.out.markdup_bams
        
        samtool_index(markdup_bams)
        indexed_bams = samtool_index.out.markdup_bams
        
        // Variant calling from processed BAMs
        if (!SKIP_SNV) {
            germVarCall(indexed_bams, params.ref, params.bed)
            raw_snv = germVarCall.out.raw_vcfs
            KeepPASS(raw_snv)
            snv_ch = KeepPASS.out.final_vcfs
        }
        
        if (!SKIP_SV) {
            SV_germline(indexed_bams, params.ref)
            sv_ch = SV_germline.out.sv_germline
        }
        
        if (!SKIP_CNV) {
            CNV_germline(indexed_bams, params.cnv_baseline)
            cnv_ch = CNV_germline.out.cnv_Germline
        }
        
        if (!SKIP_MITO) {
            mito_variant(indexed_bams, params.mito_reference, params.mito_annotation)
            mito_ch = mito_variant.out.mito_vcf
        }
    }

    emit:
    snv_variants = snv_ch
    sv_variants = sv_ch
    cnv_variants = cnv_ch
    mito_variants = mito_ch
    qc_reports = qc_reports
    aligned_bams = indexed_bams
}

/*
 * P2G WORKFLOW
 * Phenotype-to-Genotype analysis workflow
 */
workflow P2G {
    take:
    snv_variants
    phenotypes

    main:
    // Initialize channels
    vep_ch = Channel.empty()
    acmg_class = Channel.empty() // added on 23082025
    //lof_ch = Channel.empty()
    pvs1_lof_ch = Channel.empty()
    vus_ch = Channel.empty()
    exomiser_ch = Channel.empty()
    sec_ch = Channel.empty()
    csv_output = Channel.empty()
    vcf_output = Channel.empty()
    final_p2g_results = Channel.empty()
    
    boolean SKIP_SNV = params.skip_snv ?: false
    boolean SEC_FIND = params.secondary_find ?: false
    String GENDER = params.gender?.toLowerCase() ?: 'male'
    String FREQUENCY_COLUMN = params.frequency_column ?: 'gnomAD_exomes_SAS_AF'

    log.info "Running P2G (Phenotype2Genotype) analysis..."

    // VEP annotation and LoF analysis (only on SNVs if available)
    if (!SKIP_SNV && snv_variants) {
        annotVEP(snv_variants, params.cachee, params.dirplugin, params.ref, 
                params.assembly, params.dbNSFP, params.loftool, 
                params.CADDsnv, params.CADDindel, params.dbscSNV,
                params.dosage_collins, params.spliceaiindel, params.spliceaisnv, params.pli, params.phenotypesfile, params.loeuf,
                params.clinpred, params.alphamissense, params.primateai, params.revel, params.clinvar)
        vep_ch = annotVEP.out.annotVEP_vcfs
        
        vcf_vep_input = vep_ch.join(snv_variants)
        // ACMG 
        acmgpre(snv_variants, params.cache_bias, params.sd_bias, params.ref_bias)
        acmg_pre = acmgpre.out.acmg_pre
        acmgclass(acmg_pre, params.parameter_bias) 
        acmg_class = acmgclass.out.acmg_class
        // ACMG
        
        vcfvep_process(vcf_vep_input, acmg_class, params.clingen_file,params.gender)
        csv_output = vcfvep_process.out.csv_output
        vcf_output = vcfvep_process.out.vcf_output
       
        //runLossOfFunctionAnalysis(vep_ch)
        //lof_ch = runLossOfFunctionAnalysis.out.lof_variants
        lofDetection(vep_ch, 
                    Channel.fromPath(params.clingen_file),
                    Channel.fromPath(params.gene_data_file))
        pvs1_lof_ch = lofDetection.out.lof_results
        
        // VusPrize for P2G
        vusVEP(snv_variants, params.cachee, params.dirplugin, params.ref, 
              params.assembly, params.dbnsfp4, params.loftool, 
              params.CADDsnv, params.CADDindel, params.dbscSNV)
        vus_vep = vusVEP.out.vus_vep
        
        processVEPOutput(vus_vep)
        processed_vus = processVEPOutput.out.processed_vcf
        
        vusprize_script_ch = Channel.value(file(params.vusprize_script))
        dependent_script_ch = Channel.value(file(params.dependent_script))
        runVusPrioritization(processed_vus, 
                           vusprize_script_ch,
                           dependent_script_ch)
        vus_ch = runVusPrioritization.out.vusprize_output
        
        // Exomiser for P2G
        hpoIds_ch = Channel.fromList(phenotypes)
        
        MODIFY_CONFIG(Channel.fromPath(params.template_yml),
                     hpoIds_ch.collect())
        modified_config = MODIFY_CONFIG.out.config
        
        vcf_files = snv_variants.map { sample_id, file -> tuple(sample_id, file) }
        
        run_prioritizer(vcf_files,
                      modified_config,
                      Channel.value(params.genome),
                      Channel.fromPath(params.exomiser_properties))
        //exomiser_ch = run_prioritizer.out.exomiser_results
        exomiser_tsv = run_prioritizer.out.variants_tsv
        exomiser_json = run_prioritizer.out.variants_json
        exomiser_ch = exomiser_tsv.join(exomiser_json)

        // Python script to process VEP + Exomiser results
        //runVEPExomiserScript(vep_ch, exomiser_ch)
        //vep_exomiser_results = runVEPExomiserScript.out.combined_results
        final_merge_input = vcf_output
            .join(run_prioritizer.out.variants_tsv)
            .join(run_prioritizer.out.variants_json)
            .join(snv_variants)
 
            //.map { sample_id, merged_vcf,  tsv_tuple, json_tuple, original_vcf -> 
                    // Assuming exomiser_results contains [tsv_file, json_file]
              //     tuple(sample_id, merged_vcf, tsv_tuple[1], json_tuple[1], original_vcf) }

        mergeExomiserVEPResults(final_merge_input)
        final_p2g_results = mergeExomiserVEPResults.out.final_p2g_results
        
        // Secondary findings for P2G
        if (SEC_FIND) {
            //acmg_sf = Channel.fromPath(params.secondaryfindings_gene_file)
            runSecondaryFindings(pvs1_lof_ch, exomiser_tsv, vcf_output)
            sec_ch = runSecondaryFindings.out.secondary_findings
        }
    }

    emit:
    annotated_vars = vep_ch
    //lof_variants = lof_ch
    csv_output = csv_output    
    vcf_output = vcf_output 
    pvs1_lof_results = pvs1_lof_ch
    vusprize_scores = vus_ch
    exomiser_results = exomiser_ch
    secondary_finds = sec_ch
    final_p2g_results = final_p2g_results 
}

/*
 * G2P WORKFLOW
 * Genotype-to-Phenotype analysis workflow
 */
workflow G2P {
    take:
    snv_variants
    //sv_variants
    //cnv_variants
    //mito_variants

    main:
    // Initialize channels
    vep_ch = Channel.empty()
    acmg_class = Channel.empty() // added on 23082025
    lof_ch = Channel.empty()
    //vus_ch = Channel.empty()
    g2p_ch = Channel.empty()
    pvs1_lof_ch = Channel.empty()
    g2p_primary = Channel.empty()
    csv_output = Channel.empty()
    vcf_output = Channel.empty()    


    boolean SKIP_SNV = params.skip_snv ?: false
    //String GENDER = params.gender?.toLowerCase() ?: 'unknown'
    //String FREQUENCY_COLUMN = params.frequency_column ?: 'gnomAD_AF'

    log.info "Running G2P (Genotype2Phenotype) analysis..."
    if (!SKIP_SNV && snv_variants) {
        annotVEP(snv_variants, params.cachee, params.dirplugin, params.ref,
                params.assembly, params.dbNSFP, params.loftool,
                params.CADDsnv, params.CADDindel, params.dbscSNV,
                params.dosage_collins, params.spliceaiindel, params.spliceaisnv, params.pli, params.phenotypesfile, params.loeuf,
                params.clinpred, params.alphamissense, params.primateai, params.revel, params.clinvar)
        vep_ch = annotVEP.out.annotVEP_vcfs

        vcf_vep_input = vep_ch.join(snv_variants)
        
        // ACMG 
        acmgpre(snv_variants, params.cache_bias, params.sd_bias, params.ref_bias)
        acmg_pre = acmgpre.out.acmg_pre
        acmgclass(acmg_pre, params.parameter_bias) 
        acmg_class = acmgclass.out.acmg_class
        // ACMG 

        vcfvep_process(vcf_vep_input, acmg_class, params.clingen_file, params.gender)
        csv_output = vcfvep_process.out.csv_output
        vcf_output = vcfvep_process.out.vcf_output
 
        vcf_vep_input = vep_ch.join(snv_variants)

        runG2P(vcf_output)
        g2p_ch = runG2P.out.final_g2p_score
        g2p_primary = runG2P.out.primary_g2p_score        

        //runLossOfFunctionAnalysis(vep_ch)
        //lof_ch = runLossOfFunctionAnalysis.out.lof_variants
        lofDetection(vep_ch,
                    Channel.fromPath(params.clingen_file),
                    Channel.fromPath(params.gene_data_file))
        pvs1_lof_ch = lofDetection.out.lof_results
        /*
        // VusPrize for G2P
        vusVEP(snv_variants, params.cachee, params.dirplugin, params.ref,
              params.assembly, params.dbnsfp4, params.loftool,
              params.CADDsnv, params.CADDindel, params.dbscSNV)
        vus_vep = vusVEP.out.vus_vep

        processVEPOutput(vus_vep)
        processed_vus = processVEPOutput.out.processed_vcf

        vusprize_script_ch = Channel.value(file(params.vusprize_script))
        dependent_script_ch = Channel.value(file(params.dependent_script))
        runVusPrioritization(processed_vus,
                           vusprize_script_ch,
                           dependent_script_ch)
        vus_ch = runVusPrioritization.out.vusprize_output
        */

    }

    emit:
    annotated_vars = vep_ch
    lof_variants = pvs1_lof_ch
    //vusprize_scores = vus_ch
    g2p_vals = g2p_ch
    g2p_primary = g2p_primary
    pvs1_lof_results = pvs1_lof_ch
    csv_output = csv_output
    vcf_output = vcf_output
    
}

/*
 * MAIN WORKFLOW
 * Entry point that coordinates preprocessing and mode-specific analysis
 */
workflow {
    // Parse phenotypes safely at the start
    def parsedPhenotypes = params.phenotypes instanceof String ?
        new groovy.json.JsonSlurper().parseText(params.phenotypes).flatten() :
        params.phenotypes.flatten()

    // Check for phenotypes and log warning at the start if needed
    if (parsedPhenotypes.isEmpty()) {
        log.warn "No phenotypes provided. Variant prioritization will not happen."
    }

    // Validation
    boolean SKIP_SNV = params.skip_snv ?: false
    boolean SKIP_SV = params.skip_sv ?: false
    boolean SKIP_CNV = params.skip_cnv ?: false
    boolean SKIP_MITO = params.skip_mito ?: false
    String MODE = params.mode?.toUpperCase() ?: 'P2G'
    boolean SEC_FIND = params.secondary_find ?: false
    String VCF_TYPE = params.input_type == 'vcf' ? params.vcf_type?.toLowerCase() : null
    String GENDER = params.gender?.toLowerCase() ?: 'unknown'
    String FREQUENCY_COLUMN = params.frequency_column ?: 'gnomAD_AF'

    // Validate mode
    if (!(MODE in ['P2G', 'G2P'])) {
        error "Invalid mode: ${MODE}. Must be 'P2G' or 'G2P'"
    }

    // Validate phenotypes for P2G mode
    if (MODE == 'P2G' && parsedPhenotypes.isEmpty()) {
        error "P2G mode requires phenotypes to be provided"
    }

    // Validate gender
    if (!(GENDER in ['male', 'female', 'unknown'])) {
        error "Invalid gender: ${GENDER}. Must be 'male', 'female', or 'unknown'"
    }

    // Validate platform
    def validPlatforms = ['illumina_pe', 'illumina_se', 'nanopore', 'bgi_pe', 'bgi_se', 'mgi_pe', 'mgi_se']
    if (params.input_type == 'fastq' && !(params.platform in validPlatforms)) {
        error "Unsupported platform: ${params.platform}. Please choose from ${validPlatforms.join(', ')}."
    }

    // Validate frequency column
    def validFreqColumns = [
        'gnomAD_exomes_AF', 'gnomAD_exomes_AFR_AF', 'gnomAD_exomes_AMR_AF', 
        'gnomAD_exomes_EAS_AF', 'gnomAD_exomes_FIN_AF', 'gnomAD_exomes_NFE_AF', 'gnomAD_exomes_SAS_AF',
    ]
    
    if (!(FREQUENCY_COLUMN in validFreqColumns)) {
        log.warn "Frequency column '${FREQUENCY_COLUMN}' is not in the standard list. Make sure it exists in your VEP annotations."
    }

    log.info "Analysis Mode: ${params.mode}"
    log.info "Platform: ${params.platform}"
    log.info "Gender: ${GENDER}"
    log.info "Frequency column: ${FREQUENCY_COLUMN}"
    log.info "Skip flags - SNV: ${SKIP_SNV}, SV: ${SKIP_SV}, CNV: ${SKIP_CNV}, MITO: ${SKIP_MITO}"

    // Input channel setup based on input type and platform
    def input_ch
    
    if (params.input_type == 'fastq') {
        if (params.platform == 'illumina_pe' || params.platform == 'bgi_pe' || params.platform == 'mgi_pe') {
            input_ch = Channel
                .fromFilePairs("${params.data_dir}/*_R{1,2}.fastq.gz", size: 2)
                .ifEmpty { error "No paired-end FASTQ files found in ${params.data_dir}" }
        } else if (params.platform == 'illumina_se' || params.platform == 'nanopore' || params.platform == 'bgi_se' || params.platform == 'mgi_se') {
            input_ch = Channel
                .fromFilePairs("${params.data_dir}/*.fastq.gz", size: 1)
                .ifEmpty { error "No single-end FASTQ files found in ${params.data_dir}" }
        } else {
            error "Unsupported platform: ${params.platform}. Please choose from 'illumina_pe', 'illumina_se', 'nanopore', 'bgi_pe', 'bgi_se', 'mgi_pe', or 'mgi_se'."
        }
    } else if (params.input_type == 'bam') {
        input_ch = Channel
            .fromPath("${params.data_dir}/*.bam")
            .map { file -> tuple(file.simpleName, file) }
            .ifEmpty { error "No BAM files found in ${params.data_dir}" }
    } else if (params.input_type == 'vcf') {
        input_ch = Channel
            .fromPath("${params.data_dir}/*.vcf.gz")
            .map { file -> tuple(file.simpleName, file) }
            .ifEmpty { error "No VCF files found in ${params.data_dir}" }
    } else {
        error "Invalid input_type: ${params.input_type}. Must be 'fastq', 'bam', or 'vcf'"
    }

    // Run preprocessing
    PREPROCESSING(input_ch)

    // Run mode-specific analysis
    if (MODE == 'P2G') {
        P2G(
            PREPROCESSING.out.snv_variants,
            //PREPROCESSING.out.sv_variants,
            //PREPROCESSING.out.cnv_variants,
            //PREPROCESSING.out.mito_variants,
            parsedPhenotypes
        )
        
        // Emit P2G results
        emit:
        snv_variants = PREPROCESSING.out.snv_variants
        annotated_vars = P2G.out.annotated_vars
        //lof_variants = P2G.out.lof_results
        pvs1_lof_results = P2G.out.pvs1_lof_results
        vusprize_scores = P2G.out.vusprize_scores
        exomiser_results = P2G.out.exomiser_results
        secondary_finds = P2G.out.secondary_finds
        qc_reports = PREPROCESSING.out.qc_reports
        aligned_bams = PREPROCESSING.out.aligned_bams
        final_p2g_results = P2G.out.final_p2g_results
        merged_vcf_vep = P2G.out.csv_output
        
    } else { // G2P mode
        G2P(
            PREPROCESSING.out.snv_variants,
            //PREPROCESSING.out.sv_variants,
            //PREPROCESSING.out.cnv_variants,
            //PREPROCESSING.out.mito_variants
        )
        
        // Emit G2P results
        emit:
        snv_variants = PREPROCESSING.out.snv_variants
        //sv_variants = PREPROCESSING.out.sv_variants
        //cnv_variants = PREPROCESSING.out.cnv_variants
        //mito_variants = PREPROCESSING.out.mito_variants
        annotated_vars = G2P.out.annotated_vars
        //lof_variants = G2P.out.lof_variants
        pvs1_lof_results = G2P.out.pvs1_lof_results
        ///vusprize_scores = G2P.out.vusprize_scores
        //exomiser_results = G2P.out.exomiser_results
        //secondary_finds = G2P.out.secondary_finds
        qc_reports = PREPROCESSING.out.qc_reports
        aligned_bams = PREPROCESSING.out.aligned_bams
        g2p_primary = G2P.out.g2p_primary
    }
}
