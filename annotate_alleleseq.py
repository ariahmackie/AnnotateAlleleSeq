"""
This module takes a count.txt or interesting_hets.txt file and graphs the count of maternal and paternal allele alignments for all
single nucleotide polymorphisms found. It also reads a corresponding BED file and stores the genome positions for all genes in the file.
If a SNP lies within a genome region, the gene is reported in the graph.

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

# Get Command Line Arguments
program = sys.argv[0]  # annotate_alleleseq.py
results = sys.argv[1]  # counts.txt
bed = sys.argv[2]      # file.bed

path_to_results = os.getcwd() + "/" + results
path_to_bed = os.getcwd() + "/" + bed
os.mkdir("saved")
def main():
    datasets = process_results_file()
    gene_dictionary = process_bed_file()
    get_gene_matches(datasets, gene_dictionary)

def process_results_file():
    df = ImportData(results)
    datasets = BreakDataFrame(df)
    return datasets

def process_bed_file():
    bed_dataframe = ImportData(bed)
    bed_rows = BreakDataFrame(bed_dataframe)
    gene_dictionary = create_gene_dictionary(bed_rows)
    return gene_dictionary

def ImportData(filename):
    """Import tab seperated files and store as dataframe."""
    df = []
    df = pandas.read_csv(filename, sep = "\t", lineterminator = "\n")
    PrintDataFrameHead(df)
    return df

def PrintDataFrameHead(df):
    """Print beginning of Dataframe"""
    print("Dataframe Head:")
    print(df.head())

def BreakDataFrame(df):
    """Split Dataframe into seperate rows."""
    matrix = df.to_numpy()
    num_rows = len(df)
    datasets = np.array_split(matrix, num_rows)
    datasets_single_bracket = strip_brackets(datasets)
    return datasets_single_bracket

def strip_brackets(datasets):
    dataset2 = []
    for i in datasets:
        dataset2.append(i[0])
    return dataset2

def get_snp_pos(dataset):
    """Return the SNP position (column 2) in the count.txt/interesting_het.txt file"""
    snp_pos = dataset[1]
    print("SNP Position " + str(snp_pos))
    return snp_pos

def get_maternal_allele(dataset):
    """Return the maternal allele (letter) for a SNP"""
    mat_allele = dataset[7]
    print("Maternal Allele: %s " % mat_allele)
    return mat_allele

def get_paternal_allele(dataset):
    """Return the paternal allele (letter) for a SNP"""
    pat_allele = dataset[8]
    print("Paternal Allele: %s" % pat_allele )
    return pat_allele

def get_allele_index(allele):
    """return the index in count.txt file for a given nucleotide"""
    allele_dict = {"A": 9, "C": 10, "G": 11, "T": 12}
    return allele_dict[allele]

def get_maternal_count(dataset):
    """return the number of alignments that match with the maternal allele"""
    mat_allele = get_maternal_allele(dataset)
    allele_index = get_allele_index(mat_allele)
    count = dataset[allele_index]
    print("Maternal count: " + str(count))
    return count

def get_paternal_count(dataset):
    """return the number of alignments that match the paternal allele"""
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
    """return the name of a gene in a bedfile row"""
    gene_name = bed_row[3]
    return gene_name

def create_gene_dictionary(datasets):
    gene_dict = OrderedDict()
    for row in datasets:
        gene_name = extract_gene_name(row)
        range = extract_gene_ranges(row)
        gene_dict[range] = gene_name
    return gene_dict

def return_genes_covering_snp(snppos, gene_dict):
    genes_covering_snp = []
    for key in gene_dict.keys():
        if snppos in range(key[0], key[1]):
            print("gene snp and match found!")
            gene_range_and_name = (key, gene_dict[key])
            genes_covering_snp.append(gene_range_and_name)
    return genes_covering_snp

def get_gene_matches(datasets, gene_dict):
    for row in datasets:
        snppos = get_snp_pos(row)
        mat_count = get_maternal_count(row)
        pat_count = get_paternal_count(row)
        genes_covering_snp = return_genes_covering_snp(snppos, gene_dict)
        gene_name = get_gene_names(genes_covering_snp)
        filename = create_file_name(snppos, gene_name)
        graph_single_snp(mat_count, pat_count, gene_name, filename, True)

def snp_counts_to_percent(mat_count, pat_count):
    sum = mat_count + pat_count
    mat_percent = mat_count / float(sum) * 100
    pat_percent = pat_count / float(sum) * 100
    return [mat_percent, pat_percent]

def get_gene_names(genes_covering_snp):
    if len(genes_covering_snp) == 0:
        return "No Gene Match Found."
    gene_names = "Matches: "
    for i in range(len(genes_covering_snp)):
            if i == len(genes_covering_snp) - 1:
                gene_names = gene_names + genes_covering_snp[i][1]
            else:
                gene_names = gene_names +  genes_covering_snp[i][1] + ", "
    return gene_names

def create_file_name(snp_pos, gene_name):
    if gene_name == "No Gene Classified":
        gene_name = "unclassified"
    print(snp_pos)
    filename = gene_name + "_" + str(snp_pos) + ".png"
    return filename

def graph_single_snp(mat_count, pat_count, gene_name, filename, do_percentage):
    mat_color = "#70114f"
    pat_color = "#248af0"
    if do_percentage:
        counts = snp_counts_to_percent(mat_count, pat_count)
    else:
        counts = [mat_count, pat_count]
    title = "Number of SNP Alignments Per Maternal and Paternal Allele \n of the X Chromosome. " + gene_name
    width = 0.25
    x = np.arange(2)
    fig, ax = plt.subplots()
    bars = ax.bar(x, counts, color = [mat_color, pat_color])
    ax.set_title(title, fontsize = 12)
    ax.set_xticks(x)
    ax.set_xticklabels(["maternal", "paternal"])
    ax.set_ylabel('SNP Read Counts')
    ax.set_xlabel('Paternal and Maternal Chromosome')
    ax.bar_label(bars, padding = 1, fontsize = 8)
    ax.set_xlim(-0.5, 2)
#    plt.show()

    plt.savefig("saved/"+ filename)

#
main()
