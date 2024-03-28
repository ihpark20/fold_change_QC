"""
gene fc에 대한 summary statistics를 구하여 저장

- median and mad from log2 transformed fc values
"""

import pandas as pd
import numpy as np
from scipy import stats

if __name__ == "__main__":
    
    input = "data_raw/gene_fc_from_all_run_ids.tsv"
    df = pd.read_csv(input, sep="\t")
    target_genes = df["GeneName"].unique()

    all_stats = []
    for target_gene in target_genes:
        df_target = df[df["GeneName"]==target_gene]
        fc = list(df_target["FC"].values)
        fc = [x for x in fc if x>0]

        fc_log2 = np.log2(fc)
        fc_log2_median = np.median(fc_log2)
        fc_log2_mad = stats.median_abs_deviation(fc_log2)
        print(target_gene, np.round(fc_log2_median, 5), np.round(fc_log2_mad, 5))
        all_stats.append([target_gene, fc_log2_median, fc_log2_mad])
    
    all_stats_df = pd.DataFrame(all_stats)
    all_stats_df.columns = ["GeneName", "Log2FC_Median", "Log2FC_MAD"]
    all_stats_df.to_csv("example/target_genes_log2fc_median_mad.csv", sep="\t", index=False)