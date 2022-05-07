
"""
count.txt and interesting_hets.txt
chrm	snppos 	ref	mat_gtyp	pat_gtyp	c_gtyp	phase	mat_all	pat_all	cA	cC	cG	cT	winning	SymCls	SymPval	BindingSite	cnv
"""

import sys
import os
import pandas
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from collections import OrderedDict

program = sys.argv[0]
results = sys.argv[1]
bed = sys.argv[2]

path_to_results = os.getcwd() + "/" + results
path_to_bed = os.getcwd() + "/" + bed


def ImportData(filename):
    df = []
    df = pandas.read_csv(filename, sep = "\t", lineterminator = "\n")
    print("ImportData()")
    print(df.head())
    return df

def BreakDataFrame(df):
    matrix = df.to_numpy()
    num_rows = len(df)
    datasets = np.array_split(matrix, num_rows)
    return datasets

def get_snp_pos(dataset):
    snp_pos = dataset[1]
    print(snp_pos)
    return snp_pos

def get_maternal_allele(dataset):
    mat_allele = dataset[7]
    print("Maternal Allele: %s " % mat_allele)
    return mat_allele

def get_paternal_allele(dataset):
    pat_allele = dataset[8]
    print("Paternal Allele: %s" % pat_allele )
    return pat_allele

def get_allele_index(allele):
    allele_dict = {"A": 9, "C": 10, "G": 11, "T": 12}
    return allele_dict[allele]

def get_maternal_count(dataset):
    mat_allele = get_maternal_allele(dataset)
    allele_index = get_allele_index(mat_allele)
    count = dataset[allele_index]
    print("Maternal count: " + str(count))
    return count

def get_paternal_count(dataset):
    pat_allele = get_paternal_allele(dataset)
    allele_index = get_allele_index(pat_allele)
    count = dataset[allele_index]
    print("Paternal count: " + str(count))
    return count

def extract_gene_ranges(bed_row):
    start_pos = bed_row[1]
    end_pos = bed_row[2]
    return (start_pos, end_pos)

def extract_gene_name(bed_row):
    gene_name = bed_row[3]
    return gene_name

def create_gene_dictionary(datasets):
    gene_dict = OrderedDict()
    for row in datasets:
        gene_name = extract_gene_name(row[0])
        range = extract_gene_ranges(row[0])
        gene_dict[range] = gene_name
    return gene_dict

def return_genes_covering_snp(snppos, gene_dict):
    genes_covering_snp = []
    for key in gene_dict.keys():
        if snppos in range(key[0], key[1]):
            print("gene snp and match found!")
            gene_range_and_name = (key, gene_dict[key])
            genes_covering_snp.append(gene_range_and_name)
    print(genes_covering_snp)
    return genes_covering_snp

def get_gene_matches(datasets, gene_dict):
    for row in datasets:
        snppos = get_snp_pos(row[0])
        mat_count = get_maternal_count(row[0])
        pat_count = get_paternal_count(row[0])
        genes_covering_snp = return_genes_covering_snp(snppos, gene_dict)
        print("*****************")
        graph_single_snp(mat_count, pat_count, genes_covering_snp, False)

def snp_counts_to_percent(mat_count, pat_count):
    sum = mat_count + pat_count
    mat_percent = mat_count / float(sum) * 100
    pat_percent = pat_count / float(sum) * 100
    return [mat_percent, pat_percent]

def return_first_gene_covering_snp(genes_covering_snp):
    if len(genes_covering_snp) != 0:
        return "No Gene Classified"
    else:
        print(genes_covering_snp)
        return genes_covering_snp[1]

def graph_single_snp(mat_count, pat_count, genes_covering_snp, do_percentage):
    mat_color = "#70114f"
    pat_color = "#248af0"
    if do_percentage:
        counts = snp_counts_to_percent(mat_count, pat_count)
    else:
        counts = [mat_count, pat_count]
    print("****")
    print(genes_covering_snp)
    if len(genes_covering_snp)!= 0:
        gene_name = genes_covering_snp[1]
    else:
        gene_name = "Gene not classified"
    title = "Number of SNP Alignments\n Per Maternal and Paternal Allele\n of the X Chromosome.\n " + gene_name
    width = 0.25
    x = np.arange(2)
    fig, ax = plt.subplots()
    bars = ax.bar(x, counts, color = [mat_color, pat_color])
    ax.set_title(title)
    ax.set_xticks(x)
    ax.set_xticklabels(["maternal", "paternal"])
    ax.set_ylabel('SNP Read Counts')
    ax.set_xlabel('Paternal and Maternal Chromosome')
    ax.bar_label(bars, padding = 1, fontsize = 8)
    ax.set_xlim(-0.5, 2)
    plt.show()




df = ImportData(results)
datasets = BreakDataFrame(df)
first_dataset = datasets[0][0]
print(first_dataset)
get_paternal_count(first_dataset)
get_maternal_count(first_dataset)


bed_df = ImportData(bed)
bed_datasets = BreakDataFrame(bed_df)
gene_dict = create_gene_dictionary(bed_datasets)



return_genes_covering_snp(228498, gene_dict)
get_gene_matches(datasets, gene_dict)
