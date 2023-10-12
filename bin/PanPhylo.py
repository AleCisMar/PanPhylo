#!/usr/bin/env python
import os
import argparse
import shutil
import subprocess
from Bio import SeqIO, SearchIO
import numpy as np
import pandas as pd
from scipy.spatial import distance
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import to_tree
import seaborn as sns
import matplotlib.pyplot as plt

def create_directory(directory_name):
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)

def search_and_extract_records(database_file, search_string, output_directory):
    #create_directory(output_directory)
    #print(f"{output_directory} directory created")
    with open(database_file, "r") as db_file:
        header = ""
        record = ""
        found_match = False  # Flag to track if a match was found
        for line in db_file:
            if line.startswith(">"):
                if found_match and header:
                    save_record(header, record, output_directory)
                header = line.strip()
                record = ""
                # Check if the search_string is in the header
                found_match = search_string in header
            else:
                if found_match:
                    record += line
        # Check the last record in the file
        if found_match and header:
            save_record(header, record, output_directory)

def save_record(header, record, output_directory):
    # Extract the first field from the header
    first_field = header.split("|")[0].lstrip(">")
    output_filename = os.path.join(output_directory, f"{first_field}.fasta")
    with open(output_filename, "w") as output_file:
        # Write the header and the sequence
        output_file.write(header + "\n" + record)
    print(f"{first_field}.fasta saved at {output_directory}")

def move_files_to_directory(file_list_file, output_directory):
    with open(file_list_file, 'r') as list_file:
        for line in list_file:
            file_path = line.strip()
            if file_path.endswith('.fasta'):
                # Move .fasta files from the list to the created directory
                shutil.move(file_path, os.path.join(output_directory, os.path.basename(file_path)))
                print(f"{file_path} moved to {output_directory}")

def run_prodigal(output_directory):
    prodigal_out_path = os.path.join(output_directory, "prodigal.out")
    with open(prodigal_out_path, 'w') as prodigal_out_file:
        for fasta_file in os.listdir(output_directory):
            if fasta_file.endswith('.fasta'):
                fasta_path = os.path.join(output_directory, fasta_file)
                cds_file = os.path.splitext(fasta_path)[0] + '.cds'
                print(f"Running prodigal for {fasta_path}")
                subprocess.run(['prodigal', '-i', fasta_path, '-a', cds_file, '-p', 'meta', '-q'], stdout=prodigal_out_file)
                print(f"{cds_file} created")
                modify_cds_header(cds_file)

def modify_cds_header(cds_file):
    # Read the content of the .cds file
    with open(cds_file, 'r') as file:
        lines = file.readlines()
    header_count = 0  # Initialize a counter for headers
    modified_lines = []  # Create a list to store modified lines
    for line in lines:
        if line.startswith(">"):
            header_count += 1
            # Split the header by "#" and keep the first field
            header_fields = line.strip().split("#")
            header = header_fields[0].lstrip(">")

            # Check if there is a "|" in the header
            if "|" in header:
                # Extract the accession from the header
                accession = header.split("|")[0]
                # Append the sequence number to the accession
                modified_accession = f"{accession}_{header_count}"
                modified_lines.append(f">{modified_accession}\n")
                print(f"{line} modified to >{modified_accession}")
            else:
                # If no "|" is present, keep the header as is
                modified_lines.append(f">{header}\n")
                print(f"{line} modified to >{header}")
        else:
            modified_lines.append(line)
    # Write the modified content back to the .cds file
    with open(cds_file, 'w') as file:
        file.writelines(modified_lines)

def split_cds_into_faa(output_directory):
    #create_directory(os.path.join(cds_directory, "CDS"))
    cds_directory = os.path.join(output_directory, "PRODIGAL")
    for cds_file in os.listdir(cds_directory):
        if cds_file.endswith('.cds'):
            cds_path = os.path.join(cds_directory, cds_file)
            diamond_directory = os.path.join(output_directory, "DIAMOND")
            print(f"Splitting {cds_path} into individual .faa files")
            split_fasta(cds_path, diamond_directory)

def split_fasta(input_file, output_directory):
    # Parse the input file and split records into separate files
    for record in SeqIO.parse(input_file, "fasta"):
        # Get the first field of the identifier (before the first space)
        identifier = record.id.split()[0]
        # Create a new file with the identifier as the filename
        output_file_name = os.path.join(output_directory, f"{identifier}.faa")        
        # Write the record to the output file
        with open(output_file_name, "w") as output_file:
            SeqIO.write(record, output_file, "fasta")
        print(f"{output_file_name} created")

def create_diamond_databases(faa_directory):
    output = os.path.join(faa_directory, "diamond_makedb.out")
    with open(output, 'w') as output_file:
        for faa_file in os.listdir(faa_directory):
            if faa_file.endswith('.faa'):
                faa_path = os.path.join(faa_directory, faa_file)
                db_name = os.path.splitext(faa_path)[0]
                subprocess.run(['diamond', 'makedb', '--in', faa_path, '--db', db_name], stdout=output_file)
                print(f"{db_name}.dmnd created")

def run_diamond_blastp(faa_directory):
    faa_files = sorted([f for f in os.listdir(faa_directory) if f.endswith('.faa')])
    dmnd_files = sorted([d for d in os.listdir(faa_directory) if d.endswith('.dmnd')])
    
    result_file = os.path.join(faa_directory, "blastp.abc")
    #diamond_err = os.path.join(faa_directory, "diamond.err")

    with open(result_file, 'w') as result_out_file:
        for i in range(len(faa_files)):
            for j in range(i + 1, len(dmnd_files)):
                query = os.path.join(faa_directory, faa_files[i])
                database = os.path.join(faa_directory, dmnd_files[j])

                print(f"Comparing: {query} - {database}")

                subprocess.run(['diamond', 'blastp', '--query', query, '--db', database, '--outfmt', '6', 'qseqid', 'sseqid', 'bitscore'], stdout=result_out_file)

def run_mcl(output_directory):
    blastp_abc = os.path.join(output_directory, "blastp.abc")
    mcl_out = os.path.join(output_directory, "mcl.out")
    
    subprocess.run(['mcl', blastp_abc, '--abc', '-o', mcl_out])

def process_mcl_clusters(output_directory):
    mcl_out = os.path.join(output_directory, "mcl.out")
    nclusters = sum(1 for line in open(mcl_out))

    for i in range(1, nclusters + 1):
        cluster_members = []
        with open(mcl_out, 'r') as mcl_file:
            for line_num, line in enumerate(mcl_file, 1):
                if line_num == i:
                    cluster_members = line.strip().split('\t')

        # Add .faa extension to cluster members
        cluster_members = [f"{member}.faa" for member in cluster_members]

        # Concatenate cluster members into a single file
        cluster_filename = os.path.join(output_directory, f"{i}.cl")
        with open(cluster_filename, 'w') as cluster_file:
            for member in cluster_members:
                with open(os.path.join(output_directory, member), 'r') as member_file:
                    shutil.copyfileobj(member_file, cluster_file)
        print(f"{cluster_members} form a cluster. Concatenating files into {cluster_filename}")

def find_and_copy_orphan_faa_files(output_directory):
    mcl_out = os.path.join(output_directory, "DIAMOND", "mcl.out")
    faa_directory = os.path.join(output_directory, "DIAMOND")
    hmmer_directory = os.path.join(output_directory, "HMMER")

    # Create PPHMMS directory if it doesn't exist
    #create_directory(ppmmhs_directory)

    # Get the list of cluster members from mcl.out
    cluster_members = set()
    with open(mcl_out, 'r') as mcl_file:
        for line in mcl_file:
            members = line.strip().split('\t')
            cluster_members.update(members)

    # Find orphan .faa files and copy them to HMMER directory with .aln extension
    for faa_file in os.listdir(faa_directory):
        if faa_file.endswith('.faa') and faa_file.split('.faa')[0] not in cluster_members:
            source_path = os.path.join(faa_directory, faa_file)
            target_path = os.path.join(hmmer_directory, f"{faa_file.split('.faa')[0]}.aln")
            shutil.copyfile(source_path, target_path)
            print(f"{source_path} does not belong to any cluster. {target_path} created")

def run_mafft_for_clusters(faa_directory, hmmer_directory):
    cl_files = [f for f in os.listdir(faa_directory) if f.endswith('.cl')]

    for cl_file in cl_files:
        input_cl = os.path.join(faa_directory, cl_file)
        output_aln = os.path.join(hmmer_directory, f"{os.path.splitext(cl_file)[0]}.aln")
        mafft = os.path.join(hmmer_directory, "mafft.err")
        with open(mafft, 'w') as mafft_err:
            # Run MAFFT with the input .cl file and redirect the output to the .aln file
            subprocess.run(['mafft', input_cl], stdout=open(output_aln, 'w'), stderr=mafft_err)
        print(f"{output_aln} created")

def create_hmms(hmmer_directory):
    aln_files = [f for f in os.listdir(hmmer_directory) if f.endswith('.aln')]

    for aln_file in aln_files:
        input_aln = os.path.join(hmmer_directory, aln_file)
        hmm_file = os.path.join(hmmer_directory, f"{os.path.splitext(aln_file)[0]}.hmm")
        hmmbuild_out = os.path.join(hmmer_directory, "hmmbuild.out")
        with open(hmmbuild_out, 'w') as hmmbuild_out_file:
            # Create an HMM from the input .aln file
            subprocess.run(['hmmbuild', '--amino', hmm_file, input_aln], stdout=hmmbuild_out_file)
        print(f"{hmm_file} created")

def index_hmms(hmmer_directory):
    hmm_files = [f for f in os.listdir(hmmer_directory) if f.endswith('.hmm')]

    for hmm_file in hmm_files:
        input_hmm = os.path.join(hmmer_directory, hmm_file)

        # Index the HMM file using hmmpress
        subprocess.run(['hmmpress', input_hmm])

def annotate_hmms_using_hmmsearch(hmmer_directory, reference_db):
    hmm_files = [f for f in os.listdir(hmmer_directory) if f.endswith('.hmm')]

    for hmm_file in hmm_files:
        input_hmm = os.path.join(hmmer_directory, hmm_file)
        tblout_file = os.path.join(hmmer_directory, f"{os.path.splitext(hmm_file)[0]}.tblout")
        hmmsearch_out = os.path.join(hmmer_directory, "hmmsearch.out")
        with open(hmmsearch_out, 'w') as hmmsearch_out_file:
            # Run hmmsearch to annotate the HMMs and save the output in the .tblout file
            print(f"Searching {input_hmm} against {reference_db}")
            subprocess.run(['hmmsearch', '--noali', '--tblout', tblout_file, input_hmm, reference_db], stdout=hmmsearch_out_file)
            print(f"{tblout_file} created")

        # Process the .tblout file to extract the annotation
        print(f"Processing {tblout_file} to extract the annotation")
        extract_and_rename_annotation(hmmer_directory, hmm_file, tblout_file)


def extract_and_rename_annotation(hmmer_directory, hmm_file, tblout_file):
    # Load the HMMER3 tabular output file using SearchIO
    search_results = list(SearchIO.parse(tblout_file, 'hmmer3-tab'))

    annotation_count = {}  # Dictionary to keep track of annotation counts

    for result in search_results:
        if len(result) > 0:
            hit = result[0]  # Take the first hit (best match)

            # Check if the Hit object evalue attribute is greater than 1e-05
            if hit.evalue > 1e-05:
                annotation = "Uncharacterized_protein"  # Annotate as "Uncharacterized protein"
            else:
                # Read the QueryResult object accession attribute
                #accession = result.accession

                # Extract annotation before "OS=" from the best hit's description
                description = hit.description
                annotation = description.split("OS=")[0].strip()

                # Check if the annotation is "Uncharacterized protein" and look for a better annotation
                if annotation == "Uncharacterized protein":
                    # Iterate through remaining hits to find a better annotation
                    for next_hit in result[1:]:
                        next_description = next_hit.description
                        next_annotation = next_description.split("OS=")[0].strip()
                        if next_annotation != "Uncharacterized protein":
                            annotation = next_annotation
                            break  # Use the first non-"Uncharacterized protein" annotation

            # Replace spaces with underscores and remove other special characters
            annotation = "_".join(annotation.split())
            annotation = ''.join(e for e in annotation if e.isalnum() or e == '_')

            # Count the number of times this annotation has appeared
            if annotation not in annotation_count:
                annotation_count[annotation] = 1
            else:
                annotation_count[annotation] += 1

            # Append a number if the annotation has appeared more than once
            if annotation_count[annotation] > 1:
                annotation = f"{annotation}_{annotation_count[annotation]}"

            # Rename the HMM and associated files using the annotation
            annotation_base = os.path.splitext(annotation)[0]
            new_hmm_file = os.path.join(hmmer_directory, f"{annotation_base}.hmm")
            new_aln_file = os.path.join(hmmer_directory, f"{annotation_base}.aln")
            new_tblout_file = os.path.join(hmmer_directory, f"{annotation_base}.tblout")
            shutil.move(os.path.join(hmmer_directory, hmm_file), new_hmm_file)
            print(f"{os.path.join(hmmer_directory, hmm_file)} renamed as {new_hmm_file}")
            shutil.move(os.path.join(hmmer_directory, f"{os.path.splitext(hmm_file)[0]}.aln"), new_aln_file)
            old_aln_file = os.path.join(hmmer_directory, f"{os.path.splitext(hmm_file)[0]}.aln")
            print(f"{old_aln_file} renamed as {new_aln_file}")
            shutil.move(os.path.join(hmmer_directory, f"{os.path.splitext(hmm_file)[0]}.tblout"), new_tblout_file)
            old_tblout_file = os.path.join(hmmer_directory, f"{os.path.splitext(hmm_file)[0]}.tblout")
            print(f"{old_tblout_file} renamed as {new_tblout_file}")
            for ext in ['.hmm.h3f', '.hmm.h3i', '.hmm.h3m', '.hmm.h3p']:
                old_ext_file = os.path.join(hmmer_directory, f"{os.path.splitext(hmm_file)[0]}{ext}")
                new_ext_file = os.path.join(hmmer_directory, f"{annotation_base}{ext}")
                shutil.move(old_ext_file, new_ext_file)
                print(f"{old_ext_file} renamed as {new_ext_file}")

def compare_hmms_to_cds(output_directory):
    hmm_directory = os.path.join(output_directory, "HMMER")
    cds_directory = os.path.join(output_directory, "PRODIGAL")
    hmm_files = sorted([f for f in os.listdir(hmm_directory) if f.endswith('.hmm')])
    cds_files = sorted([f for f in os.listdir(cds_directory) if f.endswith('.cds')])

    # Initialize dictionaries to store hit counts and score sums for each .hmm-.cds comparison
    hit_counts = {}
    score_sums = {}

    # Numbering system for .hmm and .cds files
    hmm_numbering = {hmm_file: f"hmm{index + 1}" for index, hmm_file in enumerate(hmm_files)}
    cds_numbering = {cds_file: f"cds{index + 1}" for index, cds_file in enumerate(cds_files)}

    # Iterate through .hmm files
    for hmm_file in hmm_files:
        hmm_path = os.path.join(hmm_directory, hmm_file)
        #hmm_name = hmm_file  # Keep the file extension

        # Initialize dictionaries for this .hmm file
        hit_counts[hmm_file] = {}
        score_sums[hmm_file] = {}

        # Iterate through .cds files
        for cds_file in cds_files:
            cds_path = os.path.join(cds_directory, cds_file)
            #cds_name = cds_file  # Keep the file extension

            # Define a unique temporary file name
            tmp_output_file = os.path.join(output_directory, "HMMER", f"{hmm_numbering[hmm_file]}VS{cds_numbering[cds_file]}.tmp")
            hmmscan_out = os.path.join(output_directory, "HMMER", "hmmscan.out")
            with open(hmmscan_out, 'w') as hmmscan_out_file:
                # Run hmmscan and save the output to the temporary file
                print(f"Comparing {hmm_path} with {cds_path}")
                subprocess.run(['hmmscan', '--noali', '--tblout', tmp_output_file, hmm_path, cds_path], stdout=hmmscan_out_file)

            # Load the HMMER3 tabular output file using SearchIO
            search_results = list(SearchIO.parse(tmp_output_file, 'hmmer3-tab'))

            # Initialize variables to store hit counts and score sums
            hit_count = 0
            score_sum = 0.0

            # Check if there are matches
            if search_results:
                for result in search_results:
                    for hit in result.hits:
                        hit_count += 1
                        score_sum += hit.bitscore

            # Store the hit count and score sum in the dictionaries
            hit_counts[hmm_file][cds_file] = hit_count
            score_sums[hmm_file][cds_file] = score_sum

            # Remove the temporary file
            os.remove(tmp_output_file)

    # Initialize output files
    scores_output_file = os.path.join(output_directory, "results", "scores.txt")
    counts_output_file = os.path.join(output_directory, "results", "counts.txt")

    # Write header row to scores.txt and counts.txt
    with open(scores_output_file, 'a') as scores_file:
        hmm_names_no_ext = [os.path.splitext(hmm_file)[0] for hmm_file in hmm_files]
        scores_file.write("\t" + "\t".join(hmm_names_no_ext) + "\n")

    with open(counts_output_file, 'a') as counts_file:
        hmm_names_no_ext = [os.path.splitext(hmm_file)[0] for hmm_file in hmm_files]
        counts_file.write("\t" + "\t".join(hmm_names_no_ext) + "\n")

    # Open output files for writing data rows
    with open(scores_output_file, 'a') as scores_file, open(counts_output_file, 'a') as counts_file:
        for cds_file in cds_files:
            #cds_name = os.path.splitext(cds_file)[0]  # Remove file extension
            score_values = [str(score_sums[hmm_file][cds_file]) for hmm_file in hmm_files]
            count_values = [str(hit_counts[hmm_file][cds_file]) for hmm_file in hmm_files]
    
            scores_file.write(f"{cds_file}\t" + "\t".join(score_values) + "\n")
            counts_file.write(f"{cds_file}\t" + "\t".join(count_values) + "\n")

    print(f"PCPs written to {scores_output_file} and {counts_output_file}")

def annotate_profiles(output_directory, results_directory):
    for fasta_file in os.listdir(output_directory):
        if fasta_file.endswith('.fasta'):
            fasta_name = os.path.splitext(fasta_file)[0]
            cds_name = f"{fasta_name}.cds"
            fasta_path = os.path.join(output_directory, fasta_file)
            scores_path = os.path.join(results_directory, "scores.txt")
            counts_path = os.path.join(results_directory, "counts.txt")

            # Read the header from the corresponding .fasta file
            with open(fasta_path, 'r') as fasta_file:
                header = fasta_file.readline().strip()[1:]  # Skip the ">" character

            # Process scores.txt
            with open(scores_path, 'r') as scores_file:
                scores_content = scores_file.read()
                scores_content = scores_content.replace(cds_name, header)

            # Process counts.txt
            with open(counts_path, 'r') as counts_file:
                counts_content = counts_file.read()
                counts_content = counts_content.replace(cds_name, header)

            # Write back the updated contents to scores.txt and counts.txt
            with open(scores_path, 'w') as scores_file:
                scores_file.write(scores_content)

            with open(counts_path, 'w') as counts_file:
                counts_file.write(counts_content)

def read_data(input_file):
    # Modify this function to read data from the specified path (scores.txt)
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

def create_heatmap(counts_file_path, scores_file_path, output_heatmap_path):
    # Read data from scores.txt and counts.txt
    scores_data = pd.read_csv(scores_file_path, sep='\t', index_col=0)
    counts_data = pd.read_csv(counts_file_path, sep='\t', index_col=0)

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
    ########################################################

    # Create a custom colormap (yellow to orange to red)
    custom_cmap = sns.color_palette("YlOrRd", as_cmap=True)
    
    # Adjust the figure size, cell size, and font size for better annotation display
    plt.figure(figsize=(12, 10))  # Increase the figure size
    sns.set(font_scale=0.5)  # Adjust font size
    
    # Plot heatmap with row and column clustering, dendrograms, and color gradient
    heatmap = sns.clustermap(counts_data, annot=True , row_linkage=row_linkage, col_linkage=col_linkage,
                             cmap=custom_cmap, dendrogram_ratio=(0.2, 0.2), linewidths=0.5, linecolor='black',
                             cbar_pos=None, annot_kws={"size": fontsize},
                             xticklabels=label_size, yticklabels=label_size)  # Adjust label size
    
    # Rotate row labels for better visibility
    plt.setp(heatmap.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    
    # Save the heatmap as a PNG file at 300dpi resolution    
    plt.savefig(output_heatmap_path, dpi=300, bbox_inches='tight')

def main():
    parser = argparse.ArgumentParser(description="Compare .fasta files specified in an input list with FASTA records extracted from the VMR complete genomes database through string search to calculate Protein Clusters Profiles (PCPs), and compute a dendrogram, and a heatmap in order to know how are they related and which Protein Clusters support such relationships, respectively")
    parser.add_argument('-l', '--list', required=True, help='Input list of query .fasta files')
    parser.add_argument('-d', '--database', required=True, help='Database to search the string and extract fasta records (VMR complete genomes)')
    parser.add_argument('-s', '--string', required=True, help='String to search in the VMR complete genomes database. For example: Coronaviridae')
    parser.add_argument('-r', '--reference_db', required=True, help='Reference database to annotate PPHMMS (UniProtKB reference viral proteomes)')
    args = parser.parse_args()

    output_directory = args.string
    create_directory(output_directory)
    print(f"{output_directory} directory created")
    print(f"All data and results will be stored in {output_directory}/")

    prodigal_directory = os.path.join(output_directory, "PRODIGAL")
    create_directory(prodigal_directory)
    print(f"{prodigal_directory} directory created")

    print(f"Searching and extracting fasta records from {args.database}. Records whose header match the word '{args.string}'")
    search_and_extract_records(args.database, args.string, prodigal_directory)

    # Move .fasta files from the list to the created directory
    print(f"Moving .fasta files listed in {args.list} to {prodigal_directory}")
    move_files_to_directory(args.list, prodigal_directory)

    # Run Prodigal on the extracted .fasta files
    run_prodigal(prodigal_directory)

    diamond_directory = os.path.join(output_directory, "DIAMOND")
    create_directory(diamond_directory)
    print(f"{diamond_directory} directory created")

    # Split .cds files into individual .faa files
    print(f"Spliting .cds files into individual .faa files")
    split_cds_into_faa(output_directory)

    # Create diamond databases from .faa files
    print(f"Creating diamond databases from .faa files")
    create_diamond_databases(diamond_directory)

    # Run Diamond blastp
    print(f"Running diamond blastp")
    run_diamond_blastp(diamond_directory)

    # Run MCL
    print(f"Running MCL clustering")
    run_mcl(diamond_directory)

    # Process MCL clusters
    print(f"Processing MCL clusters")
    process_mcl_clusters(diamond_directory)

    hmmer_directory = os.path.join(output_directory, "HMMER")
    create_directory(hmmer_directory)
    print(f"{hmmer_directory} directory created")

    # Find and copy orphan .faa files to HMMER directory
    print(f"Looking for orphan .faa files")
    find_and_copy_orphan_faa_files(output_directory)

    # Run MAFFT for each cluster
    print(f"Performing MSA with MAFFT for .cl files")
    run_mafft_for_clusters(diamond_directory, hmmer_directory)

    # Create HMMs from .aln files
    print(f"Creating HMMs from .aln files with hmmbuild")
    create_hmms(hmmer_directory)

    # Index HMMs using hmmpress
    print(f"Creating index files with hmmpress")
    index_hmms(hmmer_directory)

    # Annotate HMMs using hmmsearch and process .tblout files
    print(f"Annotating HMMs using hmmsearch")
    annotate_hmms_using_hmmsearch(hmmer_directory, args.reference_db)

    results_directory = os.path.join(output_directory, "results")
    create_directory(results_directory)
    print(f"{results_directory} directory created")

    # Compare .hmm files against .cds files using hmmscan and collect hit counts and score sums
    print(f"Comparing HMMs agains .cds files with hmmscan to build Protein Clusters Profiles (PCPs)")
    compare_hmms_to_cds(output_directory)

    # Annotate PCPs scores.txt and counts.txt
    print(f"Annotating PCPs scores.txt and counts.txt")
    annotate_profiles(prodigal_directory, results_directory)

    # Step 1: Read data and calculate Bray-Curtis matrix
    scores_file_path = os.path.join(output_directory, "results", "scores.txt")
    print(f"Computing Bray-Curtis matrix from {scores_file_path}")
    ids, header, data = read_data(scores_file_path)
    matrix = calculate_bray_curtis_matrix(ids, data)
    distance_matrix = pd.DataFrame(matrix, index=ids, columns=ids)
    leaf_names = list(distance_matrix.index)

    # Step 2: Compute the distance tree
    print(f"Computing tree from distance matrix")
    linkage, dendrogram = compute_distance_tree(distance_matrix)

    # Step 3: Convert linkage matrix to a tree structure
    tree = to_tree(linkage)

    # Step 4: Generate Newick format from the tree structure
    newick_tree = generate_newick_from_tree(tree, leaf_names)

    # Step 5: Save the Newick tree to the output file
    tree_output_path = os.path.join(output_directory, "results", f"{args.string}.tree")
    with open(tree_output_path, "w") as f:
        f.write(newick_tree + ";")
    print(f"Tree written to {tree_output_path}")

    counts_file_path = os.path.join(output_directory, "results", "counts.txt")

    output_heatmap_path = os.path.join(output_directory, "results", f"{args.string}.png")
    print(f"Creating heatmap from {counts_file_path}")
    create_heatmap(counts_file_path, scores_file_path, output_heatmap_path)
    print(f"Heatmap saved in {output_heatmap_path}")

if __name__ == "__main__":
    main()