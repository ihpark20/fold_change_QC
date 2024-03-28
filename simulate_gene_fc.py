"""
Generate Gene FC data (simulate)

- Use median and mad from the log2 transformed

- CN propability (use config file) 

- We will have median and mad for each target genes.
"""

import configparser
import argparse
import numpy as np
import pandas as pd
import random
from scipy import stats
from anonymizedf.anonymizedf import anonymize



class SimulateCN:
    def __init__(self, p_hom_del, p_het_del, p_neutral, p_gain_3, p_gain_4, p_gain_5):
        prob_list = [p_hom_del, p_het_del, p_neutral, p_gain_3, p_gain_4, p_gain_5]
        self.cum_prob = np.cumsum(prob_list)

    # generate random copy number based on the probability
    def get_cn(self):
        rv = random.random()
        for i in range(len(self.cum_prob)):
            if rv < self.cum_prob[i]:
                return i            
        return 2


class SimulateFC:
    def __init__(self, sim_cn):
        self.sim_cn = sim_cn

    # output: use FC value
    def simulate_fc(self, log2_fc_median, log2_fc_mad, n_samples, tumor_percent=0.5):
        
        fc_neutral_values = self.simulate_neutral_fc(log2_fc_median, log2_fc_mad, n_samples)
        fc_all = []
        for fc_neutral in fc_neutral_values:
            target_cn = self.sim_cn.get_cn()
            sp = self.adjust_fc_cn(fc_neutral, target_cn, tumor_percent)
            fc_all.append([fc_neutral, target_cn, sp])
    
        fc_all_df = pd.DataFrame(fc_all)
        fc_all_df.columns = ["NeutralFC", "CN", "FC"]
        return fc_all_df

    def simulate_neutral_fc(self, log2_fc_median, log2_fc_mad, n_samples):        
        simulated_neutral_log2fc = stats.norm.rvs(loc=log2_fc_median, scale=log2_fc_mad, size=n_samples)
        simulated_neutral_fc = np.exp2(simulated_neutral_log2fc)
        return simulated_neutral_fc

    def adjust_fc_cn(self, fc_neutral, target_cn, tumor_percent):
        adjusted_fc = (target_cn * tumor_percent + 2 * (1 - tumor_percent)) * fc_neutral/2
        return adjusted_fc
    

def get_log2fc_median_mad():
    log2fc_median_mad = "example/target_genes_log2fc_median_mad.csv"
    log2fc_median_mad = pd.read_csv(log2fc_median_mad, sep="\t")
    #print(log2fc_median_mad)
    return log2fc_median_mad


if __name__=="__main__":

    log2fc_median_mad_df = get_log2fc_median_mad()
    template_gene = "CCND2"
    num_samples_in_a_batch = 8

    gene_df = log2fc_median_mad_df[log2fc_median_mad_df["GeneName"] == template_gene]
    log2_fc_median_original = gene_df["Log2FC_Median"].values[0]
    log2_fc_mad_original = gene_df["Log2FC_MAD"].values[0]
    
    config = configparser.ConfigParser()
    config.read("simulate_gene_fc_config.ini")    
    config = config["DEFAULT"]
    sim_cn = SimulateCN(float(config["p_hom_del"]), float(config["p_het_del"]), float(config["p_neutral"]), float(config["p_gain_3"]), float(config["p_gain_4"]), float(config["p_gain_5"]))
    
    sim_fc = SimulateFC(sim_cn)
        
    all_segments = [(480, 1.0), (80, 0.7)]
    all_dfs = []
    segment_no = 1


    for num_samples, fc_ratio in all_segments:

        log2_fc_median = log2_fc_median_original + np.log2(fc_ratio)
        log2_fc_mad = log2_fc_mad_original        

        print("Segment"+str(segment_no), log2_fc_median, log2_fc_mad)
        fc_df = sim_fc.simulate_fc(log2_fc_median, log2_fc_mad, num_samples, tumor_percent=0.5)
        fc_df["SegmentNo"] = segment_no        
        all_dfs.append(fc_df)
        segment_no += 1

    all_df = pd.concat(all_dfs)
    all_df.reset_index(drop=True, inplace=True)
    #print(all_df)

    #all_df["SampleId"] = 
    all_df["BatchId"] = list(all_df.index.values)
    all_df["BatchId"] = all_df["BatchId"].apply(lambda x:"T"+str(int(x/num_samples_in_a_batch)+1).zfill(5))

    all_df["SampleId"] = list(all_df.index.values)
    all_df["GeneName"] = "CCND2"
    an = anonymize(all_df)
    fake_df = an.fake_ids("SampleId", chaining=True).show_data_frame()
    fake_df = fake_df[["GeneName", "FC", "BatchId", "Fake_SampleId", "CN", "SegmentNo"]].copy(deep=True)
    fake_df.rename(columns={"Fake_SampleId": "SampleId"}, inplace=True)

    fake_df.to_csv("example/FC_" + template_gene + ".csv", sep="\t", index=False)
    print(fake_df)
  