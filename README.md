##### 1. 유전자별 log2FC의 median과 MAD값 저장 (example/target_genes_log2fc_median_mad.csv)

##### 2. 유전자에 대한 fold change 값을 simulation 함
###### - 위에 저장되어 있는 특정 유전자의 log2fc median과 MAD 이용
     
```{sh}
python simulate_gene_fc.py
```

##### Deletion및 Gain에 대한 확률값 설정: simulate_gene_fc_config.ini 
```{}
p_hom_del = 0.08
p_het_del = 0.02
p_neutral = 0.8
p_gain_3 = 0.03
p_gain_4 = 0.03
p_gain_5 = 0.02
p_gain_6 = 0.02
```
##### simulation 결과: example/FC_CCND2.csv

##### change point detection
```{bash}
python change_detector -i example/FC_CCND2.csv -o out/CCND2_segments.txt
```

##### 결과예시
```
Start   End     NumBatches      StartBatchId    LastBatchId     Num_Non_Outliers        Log2FC_Median   FC_Median       GeneName        SegmentSize     FC_Median_Ratio Call
0       479     60      T00001  T00060  405.0   -0.10121952822259782    0.932244620297326       CCND2   480     1.0165706469548015      False
0       559     70      T00001  T00070  473.0   -0.11618258555651556    0.9226257189864954      CCND2   560     1.006081669581645       False
480     559     10      T00061  T00070  68.0    -0.6248066851663631     0.6485066685646945      CCND2   80      0.7071672276393027      True
```

##### 결과그림
