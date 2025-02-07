#!/usr/bin/env python3
import re, pickle, sys

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import scipy.stats as stats
import statistics

species = sys.argv[1]
term = sys.argv[2]
intervals = [0, 10, 20, 30, 40]
pattern = r"CG"

def readFasta(path):
    count = 0
    seq_dict = {}
    with open(path) as FILE:
        seq = ""
        for line in FILE:
            line = line.strip()
            if line.find(">") != -1:
                if seq != "":
                    seq_dict.setdefault(seqname[1:].split(" ")[0], seq)
                    seq = ""
                seqname = line
                count += 1
            else:
                seq += line

    seq_dict.setdefault(seqname[1:].split(" ")[0], seq)

    return seq_dict

sp_dict = {}
for interval in intervals:
    with open(f"{species}.count.{term}.{interval}.bed") as FILE:
        sp_dict.setdefault(interval, FILE.readlines())

outputs_dict = {}
for interval in intervals:
    outputs_dict.setdefault(interval, [])
outputs_dict.setdefault("no_insert", [])

mix_count = 0
print(f"Total: {len(sp_dict[0])}")
for i in range( len(sp_dict[0]) ):
    counts = [ int(sp_dict[interval][i].strip().split("\t")[-1]) for interval in intervals]
    gene_region = "\t".join(sp_dict[0][i].strip().split("\t")[:-1])

    if sum(counts) == 0: # 加起來是0，那就是都沒有insert
        outputs_dict["no_insert"].append(gene_region)
    else:
        PURE_TAG = False
        for interval in intervals:
            index_interval = intervals.index(interval)

            main_count = counts[index_interval] #例如: 0-15區間DNA count
            sum_count = sum(counts) # 全部加起來的count

            if main_count == sum_count: #如果「區間DNA count」 == 全部count, 那表示這個基因裡面只有「該區間」
                outputs_dict[interval].append(gene_region)
                PURE_TAG = True
                break
        if not PURE_TAG:
            mix_count += 1

### Print Statistics and Beds
with open(f"{species}.only.{term}.statistics.txt", 'w') as S_OUT:
    for catelog in intervals + ["no_insert"]:
        print(catelog, len(outputs_dict[catelog]))
        S_OUT.write(f"{catelog}:\t{len(outputs_dict[catelog])}\n")
        with open(f"{species}.only.{term}.{catelog}.bed", 'w') as OUT:
            for line in outputs_dict[catelog]:
                OUT.write(line + "\n")
    S_OUT.write(f"mix:\t{mix_count}\n")
    print("mix", mix_count)




### Calculate GC percentage
insert_dict = {}
seq_dict = readFasta(f"/home/JuiHung/RA/genome/{species}.genome.fa")
for catelog in intervals + ["no_insert"]:
    for line in outputs_dict[catelog]:
        line = line.split("\t")
        scaffold, start, end = line[0], int(line[1]), int(line[2])

        subseq = seq_dict[scaffold][start-1:end].upper()
        matches = re.findall(pattern, subseq)
        GC_content = (subseq.count("C") + subseq.count("G"))/len(subseq)

        observed_CpG = len(matches)/len(subseq)
        expected_CpG = (GC_content/2)**2

        if len(subseq) > 0 and expected_CpG > 0:
            insert_dict.setdefault(catelog, [])
            insert_dict[catelog].append((observed_CpG, expected_CpG, observed_CpG/expected_CpG))

### Store dict
with open(f"{species}.insertDict.{term}.pkl", 'wb') as file:
    pickle.dump(insert_dict, file)


# values = []
observed_CpGs = []
expected_CpGs = []
Normalized_CpGs = []
categories = []
for category, value_data in insert_dict.items():
    # values.extend(value_data)
    observed_CpGs.extend([CpG[0] for CpG in value_data])
    expected_CpGs.extend([CpG[1] for CpG in value_data])
    Normalized_CpGs.extend([CpG[2] for CpG in value_data])
    categories.extend([category] * len(value_data))

data = pd.DataFrame({
    #'Value': values,
    'observed_CpG': observed_CpGs,
    'expected_CpG': expected_CpGs,
    'Normalized_CpG': Normalized_CpGs,
    'Category': categories
})

CpG_type = ['observed_CpG', 'expected_CpG', 'Normalized_CpG']

with open(f"{species}.only.{term}.mean.txt", 'w') as OUT:
    for i in range(len(CpG_type)):
        print(f"=== {CpG_type[i]} ===")
        OUT.write(f"=== {CpG_type[i]} ===\n")

        no_insert_CpGs = [CpG[i] for CpG in insert_dict['no_insert']]
        no_insert_mean = round(sum(no_insert_CpGs)/len(no_insert_CpGs), 4)
        no_insert_median = round(statistics.median(no_insert_CpGs),4)
        for catelog in intervals:
            if catelog in insert_dict:
                ct_CpGs = [CpG[i] for CpG in insert_dict[catelog]]
                ct_mean = round(sum(ct_CpGs)/len(ct_CpGs), 4)
                ct_median = round(statistics.median(ct_CpGs), 4)
            else:
                ct_mean = "NA"
                ct_median = "NA"
            
            if ct_mean != "NA":
                print(f"{catelog:>2}-{catelog+10} : mean={ct_mean}, median={ct_median} ({ct_mean - no_insert_mean:.4f}, {ct_median - no_insert_median:.4f})")
                OUT.write(f"{catelog:>2}-{catelog+10} : mean={ct_mean}, median={ct_median} ({ct_mean - no_insert_mean:.4f}, {ct_median - no_insert_median:.4f})\n")
            else:
                print(f"{catelog:>2}-{catelog+10} : {ct_mean}")
                OUT.write(f"{catelog:>2}-{catelog+10} : {ct_mean}\n")
        print(f"No Insert : mean={no_insert_mean}, median={no_insert_median}")        
        OUT.write(f"No Insert : mean={no_insert_mean}, median={no_insert_median}\n")
