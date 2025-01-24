# Some required python packages
import pandas as pd
import sys
import os
import anndata
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

pd.set_option('display.max_columns', 100)

# Define our base directory for the analysis
# os.chdir('/home/jovyan/cpdb_tutorial')

# Step01 : set path
cpdb_file_path = 'cellphonedb-data-5.0.0/cellphonedb.zip'
meta_file_path = 'BrM.meta.tsv'
counts_file_path = 'BrM.norm.log.h5ad'

# Step02 : check files
metadata = pd.read_csv(meta_file_path, sep = '\t')
metadata.head(3)

adata = anndata.read_h5ad(counts_file_path)
adata.shape
list(adata.obs.index).sort() == list(metadata['Cell']).sort() # should be True

# Step03 : run cpdb - subsample=False
cpdb_results = cpdb_statistical_analysis_method.call(
    cpdb_file_path = cpdb_file_path,                 # mandatory: CellphoneDB database zip file.
    meta_file_path = meta_file_path,                 # mandatory: tsv file defining barcodes to cell label.
    counts_file_path = counts_file_path,             # mandatory: normalized count matrix - a path to the counts file, or an in-memory AnnData object
    counts_data = 'hgnc_symbol',                     # defines the gene annotation in counts matrix.
#    active_tfs_file_path = active_tf_path,           # optional: defines cell types and their active TFs.
#    microenvs_file_path = microenvs_file_path,       # optional (default: None): defines cells per microenvironment.
    score_interactions = True,                       # optional: whether to score interactions or not. 
    iterations = 1000,                               # denotes the number of shufflings performed in the analysis.
    threshold = 0.1,                                 # defines the min % of cells expressing a gene for this to be employed in the analysis.
    threads = 10,                                     # number of threads to use in the analysis.
    debug_seed = 42,                                 # debug randome seed. To disable >=0.
    result_precision = 3,                            # Sets the rounding for the mean values in significan_means.
    pvalue = 0.05,                                   # P-value threshold to employ for significance.
    subsampling = False,                             # To enable subsampling the data (geometri sketching).
    subsampling_log = False,                         # (mandatory) enable subsampling log1p for non log-transformed data inputs.
    subsampling_num_pc = 100,                        # Number of componets to subsample via geometric skectching (dafault: 100).
#    subsampling_num_cells = 1000,                    # Number of cells to subsample (integer) (default: 1/3 of the dataset).
    separator = '|',                                 # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
    debug = False,                                   # Saves all intermediate tables employed during the analysis in pkl format.
    output_path = './results/BrM_normlog_stat/',                          # Path to save results.
    output_suffix = None                             # Replaces the timestamp in the output files by a user defined string in the  (default: None).
    )
