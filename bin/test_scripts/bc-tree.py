import argparse
import numpy as np
import pandas as pd
from scipy.spatial import distance
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import to_tree
import matplotlib.pyplot as plt
import seaborn as sns

def parse_args():
    parser = argparse.ArgumentParser(description="Compute Bray-Curtis distances and generate a distance tree.")
    parser.add_argument("-i", "--input", required=True, help="Input tab-delimited file with scores")
    parser.add_argument("-o", "--output", required=True, help="Output tree file in Newick format")
    return parser.parse_args()

def read_data(input_file):
    data = {}
    ids = []
    with open(input_file, 'r') as file:
        lines = file.readlines()
        header = lines[0].strip().split('\t')[1:]  # Extract column headers
        for line in lines[1:]:
            fields = line.strip().split('\t')
            id = fields[0]
            scores = list(map(float, fields[1:]))
            data[id] = scores
            ids.append(id)
    return ids, header, data

def bray_curtis_distance(v1, v2):
    numerator = sum(abs(x - y) for x, y in zip(v1, v2))
    denominator = sum(v1) + sum(v2)
    return numerator / denominator

def calculate_bray_curtis_matrix(ids, data):
    num_samples = len(ids)
    matrix = np.zeros((num_samples, num_samples))

    for i in range(num_samples):
        for j in range(i, num_samples):
            distance = bray_curtis_distance(data[ids[i]], data[ids[j]])
            matrix[i][j] = distance
            matrix[j][i] = distance  # Matrix is symmetric, so set the mirrored value

    return matrix

def compute_distance_tree(distance_matrix):
    condensed_distance = distance.squareform(distance_matrix.values)
    linkage = hierarchy.linkage(condensed_distance, method='average')
    dendrogram = hierarchy.dendrogram(linkage, labels=distance_matrix.index, orientation='right', leaf_font_size=8)
    return linkage, dendrogram

def generate_newick_from_tree(tree, leaf_names):
    if tree.is_leaf():
        return leaf_names[tree.id]
    else:
        left_node = generate_newick_from_tree(tree.left, leaf_names)
        right_node = generate_newick_from_tree(tree.right, leaf_names)
        branch_length = tree.dist
        return f"({left_node}:{branch_length},{right_node}:{branch_length})"

def main():
    args = parse_args()
    input_file = args.input
    output_file = args.output

    # Step 1: Read data and calculate Bray-Curtis matrix
    ids, header, data = read_data(input_file)
    matrix = calculate_bray_curtis_matrix(ids, data)
    distance_matrix = pd.DataFrame(matrix, index=ids, columns=ids)
    leaf_names = list(distance_matrix.index)

    # Step 2: Compute the distance tree
    linkage, dendrogram = compute_distance_tree(distance_matrix)

    # Step 3: Convert linkage matrix to a tree structure
    tree = to_tree(linkage)

    # Step 4: Generate Newick format from the tree structure
    newick_tree = generate_newick_from_tree(tree, leaf_names)

    # Step 5: Save the Newick tree to the output file
    with open(output_file, "w") as f:
        f.write(newick_tree + ";")

if __name__ == "__main__":
    main()
