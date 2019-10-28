
# Processing Modules

## Load Sparse MTX FBM
This gives something like

```gex, custom```

Here we should divide up custom into `HTO` and `ADT` (delete `custom`). 

### Sparse
```
feat_data
  gex
    mat (scipy sparse)
    barcodes
    features
  adt
    mat
    barcodes
    features    
  hto
    mat
    barcodes
    features        
```

## Initialize Cell Metadata
Calc UMI-sum and number of unique measurements for GEX, ADT, HTO

```
df_meta = calc_feat_sum_and_meas_across_cells(feat_data, inst_feat)
```

## Filter Debris
If GEX only, define debris based GEX UMI Sum (manually set based on plot). 

If HTO and GEX, defin debris based on not meeting GEX or HTO threshold. 

```
keep_gex ...
keep_hto ...

keep_cells = keep_gex + keep_hto
```

Histogram red/blue thresholding (later make scatterplot). 

The debris cells are permanently filtered. Reduces from 700K/6M -> 10-20K cells.

Clean up meta_cell by dropping metadata for debris.

## Dehash (if necessary)
Will make function for this 
```
# dense copy of HTO data
df_hto = deepcopy(df['hto'])
```

### HTO Spillover Correction (if necessary)


### Dehash Threshold
Set threshold for positive hto on a per-hto basis. Define debris, singlets, multiplets based on threshoolding. 


### 








### Dense (later)
```
df
  gex
  hto
  adt
 ```
