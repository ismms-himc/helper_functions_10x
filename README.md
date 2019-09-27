# helper_functions_10x
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

# helper_functions_10x Roadmap
Warning -- the `helper_functions_10x.py` script is currently a grab bag full of stuff, but will become more organized over time. This repo may eventually become a pip installable library or may (at least partially) become absorbed by [Clustergrammer2](https://github.com/ismms-himc/clustergrammer2). 

## Private Test Data
This repo contains pivate HIMC data that is used for demonstrating and testing functionality - so this specific repo will probably always remain private. 

# Single-Cell Data Read/Write
This set of functions will be used to read and write single cell data between several commonly used and custom data formats. This section of this README will describe these files and some schemas for what we want to save and in what data formats we want to use.

## MetaData
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
     * t_cell_clnonotye (if applicable) 
     * b_cell_clnonotye (if applicable) 
     * cell_type_broad (if applicable)
     * cell_type_narrow (if applicable)
     * cell_state (if applicable)
   
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
