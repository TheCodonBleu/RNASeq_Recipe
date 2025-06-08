version 1.0

workflow RNASeq_Visualization {
    meta {
      description: "Create a heatmap from you RNASeq quantification file."
      author: "David Maimoun (The Codon Bleu)"
  }

  input {
    Array[File] quant_files
    Int top_n = 50
    String output_name = "counts.png"
  }

  call counts_visualization {
    input:
      input_files = quant_files,
      top_n = top_n,
      output_plot = output_name
  }

  output {
    File plot = counts_visualization.plot
  }
}

task counts_visualization {
    input {
        Array[File] input_files
        Int top_n
        String output_plot
    }

    command <<<
        mkdir -p input_dir

        for f in ~{sep=' ' input_files}; do
            filename=$(basename "$f")
            name="${filename%%.*}"  # Extract part before first '.'
            cp "$f" "input_dir/${name}.quant.sf"
        done

        rnaseq_visualizing_counts \
            --input_dir input_dir \
            --output_plot ~{output_plot} \
            --top_n ~{top_n}
    >>>

    output {
        File plot = "~{output_plot}"
    }

  runtime {
    docker: "thecodonbleu/rnaseq_counts_visualization:1"
    memory: "2G"
    cpu: 1
  }
}
