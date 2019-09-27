# helper_functions_10x
This repo will contain a script with functions that are used for 10x single cell data. These functions include: 

* Input/Output of single cell data (e.g. CytoBank)
* Filtering sparse matrices (e.g. GEX UMI)
* De-hashing samples
* Read/write from/to AWS S3
* Visualize data for QC purposes
* Generate cell and feature meta-data 

# helper_functions_10x Roadmap
Warning -- the `helper_functions_10x.py` script is currently a grab bag full of stuff, but will become more organized over time. This repo may eventually become a pip installable library or may (at least partially) become absorbed by [Clustergrammer2](https://github.com/ismms-himc/clustergrammer2). 

## Private Test Data
This repo contains pivate HIMC data that is used for demonstrating and testing functionality - so this specific repo will probably always remain private. 

# Single-Cell Data Read/Write
This set of functions will be used to read and write single cell data between several commonly used and custom data formats. This section of this README will describe these files and some schemas for what we want to save and in what data formats we want to use.

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
 
