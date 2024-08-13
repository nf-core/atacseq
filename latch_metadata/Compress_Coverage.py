from os import listdir

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pyBigWig as bw


def Return_Chromosome_Level_Coverages(bigwig, chromosome, chr_len, batch=10000000):
    cov_vec = np.array([])
    for i in range(1, chr_len, batch):
        v = bigwig.values(chromosome, i, min(i + batch, chr_len))
        v = np.array(v)
        v = np.nan_to_num(v)
        cov_vec = np.append(cov_vec, v)
    return pa.array(cov_vec)


def Compute_Coverage_Across_Samples(data_path, outPath):
    samples = listdir(data_path)
    chroms_list = []
    for s in samples:
        if s.endswith(".bigWig"):
            fp = bw.open(data_path + s)
            chrom = fp.chroms()
            break
    chroms_list = list(chrom.keys())

    for c in chroms_list[0:]:
        d = {}
        for s in samples:
            if s.endswith(".bigWig"):
                print(s)
                fp = bw.open(data_path + s)
                cov = Return_Chromosome_Level_Coverages(fp, c, chrom[c])
                d[s.replace(".mLb.clN.bigWig", "")] = cov
        pq.write_table(pa.table(d), outPath + "/" + c + ".parquet")
    return outPath
