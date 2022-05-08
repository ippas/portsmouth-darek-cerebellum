import "https://gitlab.com/intelliseq/workflows/raw/fq-organize@2.1.1/src/main/wdl/tasks/fq-organize/fq-organize.wdl" as fq_organize_task
import "https://gitlab.com/intelliseq/workflows/-/raw/rna-seq-custom-hisat@1.0.1/src/main/wdl/tasks/rna-seq-custom-hisat/rna-seq-custom-hisat.wdl" as rna_seq_custom_hisat_task
import "https://gitlab.com/intelliseq/workflows/-/raw/rna-seq-ensembl-hisat@1.0.1/src/main/wdl/tasks/rna-seq-ensembl-hisat/rna-seq-ensembl-hisat.wdl" as rna_seq_ensembl_hisat_task
import "https://gitlab.com/intelliseq/workflows/raw/rna-seq-hisat@1.0.6/src/main/wdl/tasks/rna-seq-hisat/rna-seq-hisat.wdl" as rna_seq_hisat_task
import "https://gitlab.com/intelliseq/workflows/raw/rna-seq-cuffquant@1.0.5/src/main/wdl/tasks/rna-seq-cuffquant/rna-seq-cuffquant.wdl" as rna_seq_cuffquant_task
import "https://gitlab.com/intelliseq/workflows/raw/rna-seq-cuffnorm@1.0.6/src/main/wdl/tasks/rna-seq-cuffnorm/rna-seq-cuffnorm.wdl" as rna_seq_cuffnorm_task
import "https://gitlab.com/intelliseq/workflows/raw/rna-seq-concat-summary@1.0.3/src/main/wdl/tasks/rna-seq-concat-summary/rna-seq-concat-summary.wdl" as rna_seq_concat_summary_task
import "https://gitlab.com/intelliseq/workflows/raw/rna-seq-fastqc@1.0.9/src/main/wdl/tasks/rna-seq-fastqc/rna-seq-fastqc.wdl" as rna_seq_fastqc_task
import "https://gitlab.com/intelliseq/workflows/raw/rna-seq-qc-stats@1.0.3/src/main/wdl/tasks/rna-seq-qc-stats/rna-seq-qc-stats.wdl" as rna_seq_qc_stats_task
import "https://gitlab.com/intelliseq/workflows/raw/rna-seq-fastqc-overrep@1.0.3/src/main/wdl/tasks/rna-seq-fastqc-overrep/rna-seq-fastqc-overrep.wdl" as rna_seq_fastqc_overrep_task
import "https://gitlab.com/intelliseq/workflows/raw/bco@1.0.0/src/main/wdl/modules/bco/bco.wdl" as bco_module

workflow rna_seq_paired_end {

  Array[String] samples_names = select_first([samples_names])
  File? gtf
  Array[File] bam_file
  String organism_name = "Homo sapiens"
  String release_version = "100"
  String analysis_id = "no_id_provided"
  String genome_basename = sub(organism_name, " ", "_") + "_genome"
  String gtf_basename = sub(organism_name, " ", "_") + "_gtf"
  String chromosome_name = "primary_assembly"
  String summary_file_name = "summary"
  String pipeline_name = "rna_seq_paired_end"
  String pipeline_version = "1.12.9"

  scatter (index in range(length(bam_file))) {
    call rna_seq_cuffquant_task.rna_seq_cuffquant {
        input:
            gtf_file = gtf,
            bam_file = bam_file[index],
            sample_id = samples_names[index],
            index = index
    }
  }

  call rna_seq_cuffnorm_task.rna_seq_cuffnorm {
    input:
        gtf_file = gtf,
        abundances_files = rna_seq_cuffquant.abundances_file
  }

  output {

    Array[File] abundances_file = rna_seq_cuffquant.abundances_file
    Array[File] cuffnorm_output = rna_seq_cuffnorm.cuffnorm_output
  }
}
