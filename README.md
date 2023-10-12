# PanPhylo (Pangenomic and Phylogenomic analysis)

PanPhylo is intended to reconstruct the phylogenomic relationships between a group of assembled viruses (metagenomic) and a group of known viruses selected from previous knowledge of taxonomic membership. It predicts coding sequences (CDS) for every genome, performs an all vs all comparison to obtain protein clusters (PCs) which are annotated and used to build Protein Clusters Profiles (PCPs). A dendrogram is built using a Bray-Curtis dissimilarity matrix calculated from bitscore sums of per-genome CDS matching each PC (protein content dissimlarity). A heat map is built using counts of per-genome CDS matching a PC. The heat map allows to visually evaluate which proteins are shared (or not) and explain the dendrogram topology. With PanPhylo we want to answer, for metagenomic assembled viruses: how are they related to other viruses?, are they novel endemic viruses?, and how many and which genes they share with their closest relatives?

## Dependencies

PanPhylo was developed and tested with Python 3.9.18 using the following libraries:
* biopython 1.78
* diamond 0.9.14
* hmmer 3.3.2
* mafft 7.520
* matplotlib 3.7.2
* mcl 22.282
* numpy 1.26.0
* pandas 2.0.3
* prodigal 2.6.3
* scipy 1.11.3
* seaborn 0.13.0

## Installation

You can download the repository as a ZIP archive or clone the repository with:

```{bash, eval=FALSE, echo=TRUE}
git clone https://github.com/AleCisMar/PanPhylo.git
``` 
Once unpacked, within the VirTak directory:
### Option 1 (recommended): create environment with all dependencies from environment.yml file:
```{bash, eval=FALSE, echo=TRUE}
conda env create -f environment.yml
```

### Option 2: create environment with all dependencies from a single command:

```{bash, eval=FALSE, echo=TRUE}
conda create -n PanPhylo python=3.9 biopython diamond=0.9.14 hmmer mafft matplotlib mcl numpy pandas prodigal scipy seaborn=0.13.0
```

### Option 3: create a conda environment and install further dependencies one at a time:

```{bash, eval=FALSE, echo=TRUE}
conda create -n PanPhylo python=3.9
```

```{bash, eval=FALSE, echo=TRUE}
conda activate PanPhylo
```
Example:
```{bash, eval=FALSE, echo=TRUE}
conda install biopython=1.78

### Download the databases
The VMR complete genomes database is available at https://github.com/AleCisMar/VirTaK/blob/master/db/VMR_MSL38_v1_complete_genomes.fasta. For updating the database see build_kmer_database.py in https://github.com/AleCisMar/VirTaK

The reference database for annotating PCs can be downloaded from https://www.uniprot.org/uniprotkb/?query=taxonomy%5Fname:%22Viruses%22+AND+keyword:%22Reference%20proteome%22

## Execution

Make sure to activate PanPhylo environment before executing the code:

```{bash, eval=FALSE, echo=TRUE}
conda activate PanPhylo
```

```{bash, eval=FALSE, echo=TRUE}
usage: PanPhylo.py [-h] -l LIST -d DATABASE -s STRING -r REFERENCE_DB

Compare .fasta files specified in an input list with FASTA records extracted from the VMR complete genomes database through string search to calculate Protein Clusters Profiles (PCPs), and compute a
dendrogram, and a heatmap in order to know how are they related and which Protein Clusters support such relationships, respectively

optional arguments:
  -h, --help            show this help message and exit
  -l LIST, --list LIST  Input list of query .fasta files
  -d DATABASE, --database DATABASE
                        Database to search the string and extract fasta records (VMR complete genomes)
  -s STRING, --string STRING
                        String to search in the VMR complete genomes database. For example: Coronaviridae
  -r REFERENCE_DB, --reference_db REFERENCE_DB
                        Reference database to annotate PPHMMS (UniProtKB reference viral proteomes)
```