import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind

# Load matrices
counts = pd.read_csv("counts_matrix.csv")
meta = pd.read_csv("sample_metadata.csv")

# Set gene names as index
counts = counts.set_index("gene")

# Split samples by condition
control_samples = meta[meta["condition"] == "Control"]["sample"]
lesion_samples = meta[meta["condition"] == "Lesion"]["sample"]

# Extract expression matrices
control = counts[control_samples]
lesion = counts[lesion_samples]

# Compute mean per condition
mean_control = control.mean(axis=1)
mean_lesion = lesion.mean(axis=1)

# Compute log2 fold change
log2fc = np.log2((mean_lesion + 1) / (mean_control + 1))

# Compute p-values using simple t-test
p_values = []
for gene in counts.index:
    t, p = ttest_ind(control.loc[gene], lesion.loc[gene])
    p_values.append(p)

p_values = np.array(p_values)

# Volcano plot
plt.figure(figsize=(8,6))
plt.scatter(log2fc, -np.log10(p_values), color="gray")

# highlight key inflammatory genes
highlight_genes = ["IL1A","IL6","CXCL8","TNF","TSLP","CCL5","MMP9","MPO","ELANE","LTF","LCN2","CAMP"]
for gene in highlight_genes:
    if gene in log2fc.index:
        plt.scatter(log2fc.loc[gene], -np.log10(p_values[counts.index == gene]), color="red")
        plt.text(log2fc.loc[gene], -np.log10(p_values[counts.index == gene]) + 0.2, gene)

plt.xlabel("log2 Fold Change (Lesion vs Control)")
plt.ylabel("-log10(p-value)")
plt.title("Volcano Plot of Differential Expression")
plt.tight_layout()
plt.savefig("volcano.png")
plt.close()

# Heatmap of selected inflammatory genes
selected_genes = ["IL1A","IL6","CXCL8","TNF","TSLP","CCL5","MMP9","MPO","ELANE","LTF","LCN2","CAMP"]
heatmap_data = counts.loc[selected_genes]

plt.figure(figsize=(8,6))
sns.heatmap(np.log2(heatmap_data+1), cmap="viridis")
plt.title("Expression Heatmap (log2 CPM)")
plt.tight_layout()
plt.savefig("heatmap.png")
plt.close()

print("Finished! Generated volcano.png and heatmap.png")
