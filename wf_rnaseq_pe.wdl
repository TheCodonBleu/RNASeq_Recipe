version 1.0

import "tasks/validation/validate_inputs.wdl" as Validation
import "tasks/qc/reads_qc.wdl" as QC
import "tasks/qc/trimming.wdl" as Trimming
import "tasks/feature_counts.wdl" as FC
import "tasks/utils/decompress.wdl" as Decompress
import "wf/wf_salmon.wdl" as WF_Salmon
import "wf/wf_star_align.wdl" as WF_Star_Align
import "wf/wf_kallisto.wdl" as WF_Kallisto

workflow RNASeq_pe {

  meta {
    description: "Perform RNASeq analysis."
    author: "David Maimoun (The Codon Bleu)"
  }

  input {
    String samplename
    File read1
    File read2
    String aligner = 'salmon'
    
    String? species
    File? transcripts
    File? index  
    File? ref_genome_fasta
    File? gtf_file                   
  }

  # call Validation.validate_aligner {
  #   input: aligner = aligner
  # }

  # call QC.fastqc as fastqc {
  #   input: 
  #     read1 = read1,
  #     read2 = read2
  # }

  call Trimming.trimmomatic_pe as trimming {
    input: 
      read1 = read1,
      read2 = read2,
      samplename = samplename
  }
  

  if (aligner == 'salmon') {
    call WF_Salmon.wf_salmon as wf_salmon {
      input:
        read1 = trimming.read1_trimmed,
        read2 = trimming.read2_trimmed, 
        index = index,
        transcripts = transcripts,
        samplename = samplename
    }
  }

  if (aligner == 'kallisto') {
    call WF_Kallisto.wf_kallisto as wf_kallisto {
      input:
        read1 = trimming.read1_trimmed,
        read2 = trimming.read2_trimmed, 
        index = index,
        transcripts = transcripts,
        samplename = samplename
    }
  }


  # if (aligner == "star") {
  #   call Decompress.decompress_file as decompress_file {
  #     input:
  #       file = gtf_file
  #   }
    
  #   call WF_Star_Align.star_align_wf as star_align {
  #     input:
  #       read1 = trimming.read1_trimmed,
  #       read2 = trimming.read2_trimmed,
  #       index_files = index_files,
  #       ref_genome_fasta = ref_genome_fasta,
  #       gtf_file = decompress_file.file_uncompressed
  #   }  

  #   call FC.feature_counts {
  #     input:
  #       bam = star_align.bam,
  #       gtf = gtf

  #   }
  # }

  # call QC.multiqc as multiqc {
  #   input:
  #     input_files = [fastqc.zip_report] 
  # }


  output {
    File? salmon_quant_sf = wf_salmon.salmon_quant_sf
    File? salmon_lib_format_counts = wf_salmon.salmon_lib_format_counts
    File? salmon_cmd_info = wf_salmon.salmon_cmd_info
    File? salmon_quant_log = wf_salmon.salmon_quant_log

    File? kallisto_quant = wf_kallisto.kallisto_quant_abundance
    File? kallisto_run_info = wf_kallisto.kallisto_run_info

    # File qc_html = fastqc.html_report
    # File multiqc = multiqc.report
    # File? quant = salmon_quant.quant_file
    # File? counts = feature_counts.counts
  }

}
