
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

### HTO Compensation Correction (if necessary)
HTO26 and HTO27. Subtract some fraction of one HTO from the other. This is only an issue for two HTOs from 3'.

```
# corr is the correction factor
HTO1 = HTO1 - HTO2 * corr
```

This could be an issue. This is basically a contamination issue. May save the HTO data.

### Dehash Threshold
Set threshold for positive HTO on a per-hto basis. Define debris, singlets, multiplets based on threshoolding. 

Signal to noise (sn, log2 signal of highest threshold vs second-highest threshold). This signal-to-noise data can be used to adjust our previous threshold de-hashing results. We allow for three sn-thresholds that can rescue debris and multiplets or throw out singlets.  

```
sn_thresh
  debris: rescue to singlet if sn above threshold
  singlet: convert to multiplet if sn below threshold
  multiplet: rescue to singlet if sn above thresh (~1)
```

We give two de-hashing opinions for each cell (thresh, thresh-sn). Finally, the sample is assigned based on the dehashing results.  

### 








### Dense (later)
```
df
  gex
  hto
  adt
 ```
