import pandas as pd
import matplotlib.pyplot as plt
import os

# Set the working directory to the script's location
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

# Before y2t Pipeline
df = pd.read_csv("../data/gene_drug_interaction_inchikey.csv")

gen_to_drugs = df.groupby('orf')['inchiKey'].nunique()
drug_to_genes = df.groupby('inchiKey')['orf'].nunique()

def categorize(values):
    categories = values.apply(lambda x: (
        "1" if x == 1 else
        "2-3" if 2 <= x <= 3 else
        "4-5" if 4 <= x <= 5 else
        "6-10" if 6 <= x <= 10 else
        ">10"
    ))
    return categories

cat_gen = categorize(gen_to_drugs)
cat_drug = categorize(drug_to_genes)

count_gen = cat_gen.value_counts().sort_index()
count_drug = cat_drug.value_counts().sort_index()

fig, axs = plt.subplots(1, 2, figsize=(14, 5))

axs[0].bar(count_gen.index, count_gen.values, color='green')
axs[0].set_title('Number of Drugs per Gene - Before y2t - sce genes')
axs[0].set_xlabel('Number of Associated Drugs')
axs[0].set_ylabel('Number of Genes')

axs[1].bar(count_drug.index, count_drug.values, color='orange')
axs[1].set_title('Number of Genes per Drug - Before y2t - sce genes')
axs[1].set_xlabel('Number of Associated Genes')
axs[1].set_ylabel('Number of Drugs')

plt.tight_layout()

# Save plots
fig.savefig("cardinality_analysis_Before_y2t_Pipeline_sce.png", dpi=300)
fig.savefig("cardinality_analysis_Before_y2t_Pipeline_sce.svg")

plt.show()

## After y2t Pipeline
df = pd.read_csv("../results/gdi_results.csv")

gen_to_drugs = df.groupby('gene_name_sce')['inchiKey'].nunique()
drug_to_genes = df.groupby('inchiKey')['gene_name_sce'].nunique()

def categorize(values):
    categories = values.apply(lambda x: (
        "1" if x == 1 else
        "2-3" if 2 <= x <= 3 else
        "4-5" if 4 <= x <= 5 else
        "6-10" if 6 <= x <= 10 else
        ">10"
    ))
    return categories

cat_gen = categorize(gen_to_drugs)
cat_drug = categorize(drug_to_genes)

count_gen = cat_gen.value_counts().sort_index()
count_drug = cat_drug.value_counts().sort_index()

fig, axs = plt.subplots(1, 2, figsize=(14, 5))

axs[0].bar(count_gen.index, count_gen.values, color='orange')
axs[0].set_title('Number of Drugs per Gene - After y2t - sce genes')
axs[0].set_xlabel('Number of Associated Drugs')
axs[0].set_ylabel('Number of Genes')

axs[1].bar(count_drug.index, count_drug.values, color='green')
axs[1].set_title('Number of Genes per Drug- After y2t - sce genes')
axs[1].set_xlabel('Number of Associated Genes')
axs[1].set_ylabel('Number of Drugs')

plt.tight_layout()

# Save plots
fig.savefig("cardinality_analysis_After_y2t_Pipeline_sce.png", dpi=300)
fig.savefig("cardinality_analysis_After_y2t_Pipeline_sce.svg")

plt.show()

## After y2t Pipeline - tcr genes
df = pd.read_csv("../results/gdi_results.csv")

gen_to_drugs = df.groupby('gene_name_tcr')['inchiKey'].nunique()
drug_to_genes = df.groupby('inchiKey')['gene_name_tcr'].nunique()

def categorize(values):
    categories = values.apply(lambda x: (
        "1" if x == 1 else
        "2-3" if 2 <= x <= 3 else
        "4-5" if 4 <= x <= 5 else
        "6-10" if 6 <= x <= 10 else
        ">10"
    ))
    return categories

cat_gen = categorize(gen_to_drugs)
cat_drug = categorize(drug_to_genes)

count_gen = cat_gen.value_counts().sort_index()
count_drug = cat_drug.value_counts().sort_index()

fig, axs = plt.subplots(1, 2, figsize=(14, 5))

axs[0].bar(count_gen.index, count_gen.values, color='orange')
axs[0].set_title('Number of Drugs per Gene - After y2t - tcr genes')
axs[0].set_xlabel('Number of Associated Drugs')
axs[0].set_ylabel('Number of Genes')

axs[1].bar(count_drug.index, count_drug.values, color='green')
axs[1].set_title('Number of Genes per Drug - After y2t - tcr genes')
axs[1].set_xlabel('Number of Associated Genes')
axs[1].set_ylabel('Number of Drugs')

plt.tight_layout()

# Save plots
fig.savefig("cardinality_analysis_After_y2t_Pipeline_tcr.png", dpi=300)
fig.savefig("cardinality_analysis_After_y2t_Pipeline_tcr.svg")

plt.show()
