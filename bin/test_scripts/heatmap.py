import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy

# Read data from scores.txt and counts.txt
scores_data = pd.read_csv("Microviridae/results/scores.txt", sep='\t', index_col=0)
counts_data = pd.read_csv("Microviridae/results/counts.txt", sep='\t', index_col=0)

# Perform hierarchical clustering on rows with Bray-Curtis distance and average method
row_linkage = hierarchy.linkage(scores_data, method='average', metric='braycurtis')

# Perform hierarchical clustering on columns with Euclidean distance and average method
col_linkage = hierarchy.linkage(scores_data.T, method='average', metric='euclidean')

# Truncate row labels by extracting the first field before "|"
counts_data.index = counts_data.index.str.split('|').str[0]

########################################################
# Calculate the font size based on the number of rows and columns
nrows, ncols = counts_data.shape
fontsize = max(8, min(12, int(100 / max(nrows, ncols))))

# Calculate label size based on font size
label_size = max(fontsize - 7, 1)
#print(f"label size = {label_size}")
########################################################

# Create a custom colormap (yellow to orange to red)
custom_cmap = sns.color_palette("YlOrRd", as_cmap=True)

# Adjust the figure size, cell size, and font size for better annotation display
plt.figure(figsize=(12, 10))  # Increase the figure size
#sns.set(font_scale=0.9)  # Adjust font size
sns.set(font_scale=0.5)  # Adjust font size

# Plot heatmap with row and column clustering, dendrograms, and color gradient
heatmap = sns.clustermap(counts_data, annot=True , row_linkage=row_linkage, col_linkage=col_linkage,
                         cmap=custom_cmap, dendrogram_ratio=(0.2, 0.2), linewidths=0.5, linecolor='black',
                         cbar_pos=None, annot_kws={"size": fontsize},
                         xticklabels=label_size, yticklabels=label_size)  # Add a color bar label

# Rotate row labels for better visibility
plt.setp(heatmap.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)

# Save the heatmap as a PNG file at 300dpi resolution
plt.savefig("heatmap.png", dpi=300, bbox_inches='tight')

## Display the plot (optional)
##plt.show()