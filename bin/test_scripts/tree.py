import argparse
import pandas as pd
from scipy.spatial import distance
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import to_tree

def parse_args():
    parser = argparse.ArgumentParser(description="Compute a distance tree from a symmetric distance matrix.")
    parser.add_argument("-i", "--input", required=True, help="Input distance matrix file")
    parser.add_argument("-o", "--output", required=True, help="Output tree file in Newick format")
    return parser.parse_args()

def read_distance_matrix(input_file):
    df = pd.read_csv(input_file, sep='\t', index_col=0)
    return df

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

    # Step 1: Read the symmetric distance matrix
    distance_matrix = read_distance_matrix(input_file)
    leaf_names = list(distance_matrix.index) ######test

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
