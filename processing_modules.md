
# Processing Modules

## Load Sparse MTX FBM
This gives something like

```gex, custom```
Here we should divide up custom into `HTO` and `ADT`. 

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

### Dehash (if necessary)

```
# dense copy of HTO data
df_hto = deepcopy(df['hto'])
```
