import argparse
import pandas as pd
import numpy as np
import ruptures as rpt
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy import stats


# minimum 0.01
def get_data(min_fc=0.01):
    data = "/BiO2/codes/tso_abc/data/gene_fc_from_all_run_ids.tsv"
    df = pd.read_csv(data, sep="\t")
    
    df.loc[df["FC"] < 0.01, "FC"] = min_fc
    df["Log2FC"] = np.log2(df["FC"])
    return df



def get_penalties_for_PELT(data_length):
    p = np.linspace(0.01, 1, 25)
    return p * np.log(data_length)    



# find segments using PELT with a penalty
def find_segments_using_pelt(fold_change_data: np.array, penalty: float):

    algo = rpt.Pelt(model="l2").fit(fold_change_data)
    change_points = algo.predict(pen=penalty)        
    if len(change_points) < 2:
        return None

    # the change_points are the index of the first element of changed segments
    print("data legnth:", len(fold_change_data))
    print("penalty:", penalty)
    print(change_points)
    
    all_segments = []    
    start = 0

    # [start, end)
    # start - included, end - not included
    for ch in change_points:

        print(start, ch)
        num_points = ch - start

        median_log2fc = np.median(fold_change_data[start:ch])       
        median_fc = np.power(2, median_log2fc)
       
        data_items = [start, ch, num_points, median_log2fc, median_fc]
        all_segments.append(data_items)
        start = ch

    all_segments_df = pd.DataFrame(all_segments)
    all_segments_df.columns = ["Start", "End", "NumPoints", "Log2FC_Median", "FC_Median"]
    all_segments_df["Penalty"] = penalty
    return all_segments_df


def find_segments_using_pelt_across_all_penalties(fold_change_data, penalties, max_segments):

    df_all = []
    for pen in penalties:
        segments = find_segments_using_pelt(fold_change_data, pen)

        if segments is None:
            break

        if segments.shape[0] <= max_segments:
            df_all.append(segments)
        
        # if we find 2 segments (the last possible segments), then stop
        if segments.shape[0]==2:
            break

    if len(df_all)==0: return None

    df_all = pd.concat(df_all)
    df_all.reset_index(drop=True, inplace=True)
    return df_all


def reduce_segments(target_data_values, segments_df: pd.DataFrame):

    segments = segments_df[["Start", "End"]].drop_duplicates().copy(deep=True)
    segments.reset_index(drop=True, inplace=True)
    segments.loc[len(segments.index)] = [0, len(target_data_values)]
    segments["Start"] = segments["Start"].astype(int)
    segments["End"] = segments["End"].astype(int)

    segments["NumPoints"] = segments["End"] - segments["Start"]
    segments["Log2FC_Median"] = segments.apply(lambda x: np.median(target_data_values[x["Start"]:x["End"]]), axis=1)
    segments["FC_Median"] = segments["Log2FC_Median"].apply(lambda x: np.power(2, x)) 
    segments.sort_values(["Start", "End"], inplace=True)
    segments.reset_index(drop=True, inplace=True)
    return segments
    

# outliers detection
# the samples and the batches are ordered
# num_batches: the number of batches for the outlier detection
# w: weight 
# lo limit: median - w * mad
# hi limit: median + w * mad
# output column
# NumOutliersDetected: 몇 번의 test에서 outlier로 detection 되었는지
# NumTests: 모두 몇 번의 test를 했는지,
# Outlier: Outlier 여부 (True or False), 한번이라도 outlier로 detection되었으면 outlier로 표시
def mark_outliers_from_each_batch(df: pd.DataFrame, data_column: str, batch_column: str, num_batches: int = 3, w: float = 3, out_column: str = "Outlier"):

    # unqiue batches
    batches = df[batch_column].drop_duplicates().values
    
    #sliding batches
    is_outlier_all = []
    for i in range(len(batches) - num_batches + 1):            
        selected_batches = batches[i:i+num_batches]        
        batch_df = df[df[batch_column].isin(selected_batches)]        
        data_values = np.array(batch_df[data_column])
        data_median = np.median(data_values)
        
        q25 = np.quantile(data_values, 0.25)
        q75 = np.quantile(data_values, 0.75)

        iqr = q75 - q25
        lo = data_median - w * iqr
        hi = data_median + w * iqr

        outside_lo = data_values < lo 
        outside_hi = data_values > hi
        outside = outside_lo | outside_hi
        is_outlier = pd.DataFrame({"ID": batch_df.index, "IsOutlier": outside})
        is_outlier_all.append(is_outlier)

    is_outlier_all = pd.concat(is_outlier_all)
    num_outlier_detected = is_outlier_all.groupby("ID").sum()
    num_outlier_detected.columns = ["NumOutliersDetected"]

    test_count = is_outlier_all["ID"].value_counts()
    test_count.name = "NumTests"
    
    outlier_df = pd.merge(num_outlier_detected, test_count, left_index=True, right_index=True, how="inner")
    df_with_outliers_info = pd.merge(df, outlier_df, left_index=True, right_index=True, how="inner")    
    df_with_outliers_info[out_column] = df_with_outliers_info["NumOutliersDetected"] > 0

    return df_with_outliers_info
   

# the ends of each segment are both-closed (left-included, right-included)
def detect_change_points(gene_fc_df: pd.DataFrame):
    
    gene_fc_outlier_marked = mark_outliers_from_each_batch(gene_fc_df, "Log2FC", "BatchId", num_batches=3, w=3, out_column="Outlier")
    gene_fc_normal = gene_fc_outlier_marked[gene_fc_outlier_marked["Outlier"] == False].copy(deep=True)    

    target_values = np.array(gene_fc_normal["Log2FC"].values)    

    # PELT 알고리즘에서 사용할 penalties 값을 구한다.
    penalties = get_penalties_for_PELT(len(target_values))    
    max_segments = 5

    all_segments = find_segments_using_pelt_across_all_penalties(target_values, penalties, max_segments)    
    original_index = list(gene_fc_normal.index.values)
    
    #print(gene_fc_outlier_marked)
    batchId_from_original = gene_fc_outlier_marked["BatchId"]

    segments_batch_corrected = []
    if all_segments is not None and not all_segments.empty:   

        all_segments = reduce_segments(target_values, all_segments)      
        for _, segment_info in all_segments.iterrows():
            start = int(segment_info["Start"])
            end = int(segment_info["End"]-1)
            start_original = original_index[start]
            end_original = original_index[end]
            
            selected_batches = batchId_from_original[start_original:end_original+1]
            selected_batches = selected_batches.value_counts()

            final_batches = []
            for batch_id, n_samples in selected_batches.items():
                all_samples_in_batch = (gene_fc_outlier_marked["BatchId"]==batch_id).sum()
                
                if n_samples / all_samples_in_batch > 0.5:
                    final_batches.append(batch_id)
            
            v_df = gene_fc_outlier_marked[gene_fc_outlier_marked["BatchId"].isin(final_batches)]
            findex = v_df.index[0]
            lindex = v_df.index[-1]

            start_batch = v_df["BatchId"].values[0]
            last_batch = v_df["BatchId"].values[-1]

            segments_batch_corrected.append([findex, lindex, len(final_batches), start_batch, last_batch, segment_info["NumPoints"], segment_info["Log2FC_Median"], segment_info["FC_Median"]])

        segments_batch_corrected_df = pd.DataFrame(segments_batch_corrected)
        segments_batch_corrected_df.columns = ["Start", "End", "NumBatches", "StartBatchId", "LastBatchId", "Num_Non_Outliers", "Log2FC_Median", "FC_Median"]
        return segments_batch_corrected_df


    return None
    


def create_argparser():
    parser = argparse.ArgumentParser(
                    prog='Detector',
                    description='Detect change point from the FC values')
    
    parser.add_argument("-i", help="input fc values (with GeneName, FC, BatchId, SampleId)")
    parser.add_argument("-o", help="out file")
    return parser

# how to adjust sex of samples.
# hmm... ahha, fold change values are already sex adjusted.
if __name__ == "__main__":

    parser = create_argparser()
    args = parser.parse_args()

    input_fc = args.i
    out = args.o
 
    gene_fc_df = pd.read_csv(input_fc, sep="\t")
    gene_name = gene_fc_df["GeneName"].unique()[0]
    
    min_segment_size = 16
    lo_fc_ratio = 0.85
    hi_fc_ratio = 1.15

    # we will use log2fc for the change detection algorithm
    gene_fc_passed = gene_fc_df[gene_fc_df["FC"]>0.01].copy(deep=True)
    gene_fc_passed["Log2FC"] = gene_fc_passed["FC"].apply(np.log2)

    fc_median = np.median(gene_fc_passed["FC"])
    print("Global Median:", str(fc_median))


    segments_df = detect_change_points(gene_fc_passed)
    
   
    if segments_df is not None:
        segments_df["GeneName"] = gene_name
    
        segments_df["SegmentSize"] = segments_df.apply(lambda x: x["End"] - x["Start"] + 1, axis=1)
        segment_size_condition = segments_df["SegmentSize"] >= min_segment_size
        segments_df["FC_Median_Ratio"] = segments_df["FC_Median"] / fc_median        

        fc_condition = (segments_df["FC_Median_Ratio"] < lo_fc_ratio) | (segments_df["FC_Median"] > hi_fc_ratio)
        segments_df["Call"] = fc_condition

        print(segments_df)

        segments_df.to_csv(out, sep="\t", index=False)


    bks = list(segments_df["Start"].values)
    bks.append(gene_fc_df.shape[0])
    bks = [x for x in bks if x!=0]

    rpt.display(gene_fc_df["FC"], bks, )
    plt.suptitle(gene_name, x=0.1, fontsize=9)
    plt.savefig(fname=out+".png", dpi=150)