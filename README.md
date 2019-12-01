# himc_helper_functions
This repo will contain a script with functions that are used for 10x single cell data. These functions include:

* Input/Output of single cell data (e.g. CytoBank)
* Filtering sparse matrices (e.g. GEX UMI)
* De-hashing samples
* Read/write from/to AWS S3
* Visualize data for QC purposes
* Generate cell and feature meta-data
* Harmonize TCR/BCR Clonotypes across lanes (treat as cell meta-data)

# Data Flow Diagram

```
Disk             In-Memory                Disk
----     ------------------------       ---------
MTX  ->  feature_data  ->   df      ->   parquets    (-> Database)
           (sparse)       (dense)   ->   MTX
                                    ->   CytoBank
                                    ->   Zarr
```

# himc_helper_functions Roadmap
Warning -- the `himc_helper_functions.py` script is currently a grab bag full of stuff, but will become more organized over time. This repo may eventually become a pip installable library or may (at least partially) become absorbed by [Clustergrammer2](https://github.com/ismms-himc/clustergrammer2).

# Single-Cell Data Read/Write
This set of functions will be used to read and write single cell data between several commonly used and custom data formats. This section of this README will describe these files and some schemas for what we want to save and in what data formats we want to use.

## Metadata Overview
This section outlines the metadata we would like to store for different entities (e.g. cells). For the time being we are thinking of storing this metadata as separate parquet files (with metadata type as column and entities as rows). The advantage to this workflow is that we can add cell-level metadata as it becomes available (e.g. cell type prediction) without having to re-write the underlying data (e.g. GEX data). We should also look into [Anndata](https://scanpy.readthedocs.io/en/stable/basic_usage.html#anndata) used by Scanpy.

Embedding cell metadata into the MTX file format could be done (by treating numeric metadata as a 'Custom' feature, but the same cannot be done for feature-based metadata (e.g. we would not want to have a metadata barcode).

* cell_meta_data

  * Gene Expression Level Meta-data
     * gex_umi_sum
     * num_genes_meas
     * fraction_mito_umi
     * sum_mito_umi
     * mean_mito_umi (can be used in place of multiple genes)
     * mean_ribo_umi (can be used in place of multiple genes)

  * Feature Level Meta-Data
     * adt_umi_sum (if applicable)
     * hto_umi_sum (if applicable)
     * hto_first_vs_secont_highest_ratio (if applicable)

  * Cell Type Level Meta-Data
     * cell_type_broad (if applicable)
     * cell_type_narrow (if applicable)
     * cell_state (if applicable)
     * t_cell_clnonotye (if applicable)
     * b_cell_clnonotye (if applicable)
     * Subject (if applicable)
     * Sample-Metadata (e.g. timepoint, if applicable)

* gex_meta_data
   * ensemble ID (for the purpose of exporting to Cell Ranger format)
   * freaction_of_cells_expressing (divide by total number of cells)
   * mean_umi (across all cells)
   * var_umi (across all cells)

* adt_meta_data
   * freaction_of_cells_expressing (divide by total number of cells)
   * mean_umi_level (across all cells)

* hto_meta_data
   * freaction_of_cells_expressing (divide by total number of cells)
   * mean_umi_level (across all cells)

## Data Formats
* Cell Ranger Version 2 Sparse Matrix MTX Format (uncompressed, read-only format)
    * barcodes.tsv
    * genes.tsv
    * matrix.mtx
* Cell Ranger Version 3 Sparse Matrix MTX Format (compressed, read and write format)
    * barcodes.tsv.gz
    * features.tsv.gz
    * matrix.mtx.gz
* Parquet Format (read and write format)
    * gex.parquet
    * adt.parquet
    * hto.parquet
    * meta_cell.parquet
    * meta_gex.parquet
    * meta_adt.parquet
    * meta_hto.parquet

* CytoBank Format
    * merge of top var gex (w/o ribo/mito)

### Cytobank Upload Format

1. General matrix format is 1 row per cell and 1 column per feature
2. First column should be a cell index, with the column header "cell_index", incremented by 0.01
3. Features should be indicated by the prefix `ADT_`, `HTO_` or `GEX_`
4. Any derived/calculated data will include `der` after after the feature type. By default, always include `GEX_der_umi_sum`, `GEX_der_unique_gene_count`, `GEX_der_mito_proportion` (proprtion of all gex umi that derive from mitochondrial genes)  `GEX_der_mito_avg` (avg umi expression of all mitochondrial genes), `GEX_der_ribo_avg` (avg umi expression of all ribosomal genes), `HTO_der_umi_sum`, `ADT_der_umi_sum`
5. Add random noise to ADT and HTO data to 2 decimal places for visualization purposes
6. Include top 500 differentially expressed genes, excluding mitochondiral or ribosomal genes
7. By default, merge data from replicate lanes but add `lane` as a feature to indicate which cells came from which lane
