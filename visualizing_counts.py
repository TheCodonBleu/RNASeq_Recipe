import os
import pandas as pd
import glob
import seaborn as sns
import matplotlib.pyplot as plt

def load_quant_files(input_dir):
    quant_files = glob.glob(f"{input_dir}/*/*.sf")
    data = {}
    for f in quant_files:
        sample_name = os.path.basename(os.path.dirname(f))
        df = pd.read_csv(f, sep='\t')
        df = df.set_index('Name')
        data[sample_name] = df['TPM']
    return pd.DataFrame(data)

def main(input_dir, output_plot, top_n=50):
    df = load_quant_files(input_dir)

    # Filter for top N most expressed transcripts (average TPM)
    top_transcripts = df.mean(axis=1).sort_values(ascending=False).head(top_n).index
    df_top = df.loc[top_transcripts]

    # Log-transform for better visualization
    df_log = df_top.applymap(lambda x: np.log2(x + 1))

    # Plot heatmap
    sns.set(font_scale=0.8)
    plt.figure(figsize=(10, 12))
    sns.heatmap(df_log, cmap="viridis", linewidths=0.5)
    plt.title(f"Top {top_n} Expressed Transcripts (TPM)")
    plt.tight_layout()
    plt.savefig(output_plot)
    print(f"✅ Heatmap saved to {output_plot}")

if __name__ == "__main__":
    import argparse
    import numpy as np

    parser = argparse.ArgumentParser(description="Generate heatmap from Salmon quant.sf files")
    parser.add_argument("--input_dir", required=True, help="Directory containing per-sample subfolders with quant.sf files")
    parser.add_argument("--output_plot", default="heatmap.png", help="Output plot filename")
    parser.add_argument("--top_n", type=int, default=50, help="Top N transcripts by mean TPM to plot")

    args = parser.parse_args()
    main(args.input_dir, args.output_plot, args.top_n)
