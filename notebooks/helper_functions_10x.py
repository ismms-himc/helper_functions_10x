# Version: 0.3.0
# This is a set of scripts that are used in processing 10x single cell data

import gzip
from scipy import io
from scipy.sparse import csc_matrix
from ast import literal_eval as make_tuple
import pandas as pd
import numpy as np
from copy import deepcopy
import os
import matplotlib.pyplot as plt

def get_version():
    print('0.3.0', 'De-hash viz and meta-data')

def make_dir(directory):
    if not os.path.exists(directory):
        os.mkdir(directory)

def load_crv3_feature_matrix(inst_path):
    # Read Barcodes
    ###########################
    # need to check whether we have tuples
    barcodes_cats = False

    # barcodes
    filename = inst_path + 'barcodes.tsv.gz'
    f = gzip.open(filename, 'rt')
    lines = f.readlines()
    f.close()

    barcodes = []
    for inst_bc in lines:
        inst_bc = inst_bc.strip().split('\t')

        if barcodes_cats == False:
            # remove dash from barcodes if necessary
            if '-' in inst_bc[0]:
                inst_bc[0] = inst_bc[0].split('-')[0]

        barcodes.append(inst_bc[0])

    # parse tuples if necessary
    if barcodes_cats:
        try:
            barcodes = [make_tuple(x) for x in barcodes]
        except:
            pass

    # Load Matrix
    #################
    mat = io.mmread(inst_path + 'matrix.mtx.gz')
    mat_csr = mat.tocsr()

    # Get Indexes of Feature Types
    ##################################
    filename = inst_path + 'features.tsv.gz'
    f = gzip.open(filename, 'rt')
    lines = f.readlines()
    f.close()

    feature_indexes = {}
    feature_lines = {}
    for index in range(len(lines)):

        inst_line = lines[index].strip().split('\t')
        inst_feat = inst_line[2].replace('Gene Expression', 'gex').replace('Antibody Capture', 'adt').replace('Custom', 'custom')


        if inst_feat not in feature_indexes:
            feature_indexes[inst_feat] = []

        feature_indexes[inst_feat].append(index)

    feature_data = {}

    for inst_feat in feature_indexes:
        feature_data[inst_feat] = {}

        feature_data[inst_feat]['barcodes'] = barcodes

        inst_indexes = feature_indexes[inst_feat]

        # Separate feature lists
        ser_lines = pd.Series(lines)
        ser_lines_found = ser_lines[inst_indexes]
        lines_found = ser_lines_found.get_values().tolist()

        # save feature lines
        feature_lines[inst_feat] = lines_found

        # save as compressed sparse column matrix (for barcode filtering)
        mat_filt = mat_csr[inst_indexes, :].tocsc()

        feature_data[inst_feat]['mat'] = mat_filt

    # Make unique feature names
    for inst_feat in feature_lines:
        feat_lines = feature_lines[inst_feat]
        feat_lines = [x.strip().split('\t') for x in feat_lines]

        # find non-unique initial feature names (add id later if necessary)
        ini_names = [x[1] for x in feat_lines]

        ini_name_count = pd.Series(ini_names).value_counts()
        duplicate_names = ini_name_count[ini_name_count > 1].index.tolist()

        new_names = [x[1] if x[1] not in duplicate_names else x[1] + '_' + x[0] for x in feat_lines]

        # quick hack to clean up names
        new_names = [x.replace('_TotalSeqB', '') for x in new_names]

        feature_data[inst_feat]['features'] = new_names

    return feature_data

def load_crv2_gene_matrix(inst_path):
    '''
    Loads gene expression data from 10x in sparse matrix format and returns a
    Pandas dataframe
    '''

    import pandas as pd
    from scipy import io
    from scipy import sparse
    from ast import literal_eval as make_tuple

    # matrix
    mat = io.mmread( inst_path + 'matrix.mtx').tocsc()


    # genes
    filename = inst_path + 'genes.tsv'
    f = open(filename, 'r')
    lines_genes = f.readlines()
    f.close()

    # make unique gene names
    #############################
    gene_list = [x.strip().split('\t') for x in lines_genes]

    # find non-unique initial gene names
    ini_names = [x[1] for x in gene_list]

    ini_name_count = pd.Series(ini_names).value_counts()
    duplicate_names = ini_name_count[ini_name_count > 1].index.tolist()
    genes = [x[1] if x[1] not in duplicate_names else x[1] + '_' + x[0] for x in gene_list]


    # barcodes
    filename = inst_path + 'barcodes.tsv'
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    cell_barcodes = []
    for inst_bc in lines:
        inst_bc = inst_bc.strip().split('\t')

        # remove dash from barcodes if necessary
        if '-' in inst_bc[0]:
          inst_bc[0] = inst_bc[0].split('-')[0]

        cell_barcodes.append(inst_bc[0])

    # parse tuples if necessary
    try:
        cell_barcodes = [make_tuple(x) for x in cell_barcodes]
    except:
        pass

    try:
        genes = [make_tuple(x) for x in genes]
    except:
        pass

    # generate feature_data
    feature_data = {}
    feature_data['gex'] = {}
    feature_data['gex']['features'] = genes
    feature_data['gex']['barcodes'] = cell_barcodes
    feature_data['gex']['mat'] = mat

    return feature_data

def plot_umi_levels(feature_data, feature_type='gex', logy=True, logx=False,
                    figsize=(10,5), min_umi=0, max_umi=1e8, zscore_features=False):
    '''
    This function takes a feature data format or dictionary of DataFrames and plots
    UMI levels
    '''

    if 'mat' in feature_data[feature_type]:
        mat_csc = feature_data[feature_type]['mat']

        if zscore_features:
            print('z-scoring feature_data')
            inst_df = pd.DataFrame(data=mat_csc.todense(), columns=feature_data[feature_type]['barcodes'])

            net.load_df(inst_df)
            net.normalize(axis='row', norm_type='zscore')
            inst_df = net.export_df()

            # sort
            ser_sum = inst_df.sum(axis=0).sort_values(ascending=False)

        else:
            # drop cells with fewer than threshold events
            ser_sum = mat_csc.sum(axis=0)
            arr_sum = np.asarray(ser_sum[0,:])

            # sort
            ser_sum = pd.Series(arr_sum[0], index=feature_data[feature_type]['barcodes']).sort_values(ascending=False)

        # filter
        ser_sum = ser_sum[ser_sum >= min_umi]
        ser_sum = ser_sum[ser_sum <= max_umi]

    else:
        inst_df = feature_data[feature_type]

        if zscore_features:
            print('zscore features')
            net.load_df(inst_df)
            net.normalize(axis='row', norm_type='zscore')
            inst_df = net.export_df()

        # sort
        ser_sum = inst_df.sum(axis=0).sort_values(ascending=False)

        # filter
        ser_sum = ser_sum[ser_sum >= min_umi]
        ser_sum = ser_sum[ser_sum <= max_umi]

    ser_sum.plot(logy=logy, logx=logx, figsize=figsize)
    return ser_sum

def filter_barcodes_by_umi(feature_data, feature_type, min_umi=0, max_umi=1e8,
                                      make_sparse=True, zscore_features=False):

    # feature data format
    ########################
    if 'mat' in feature_data[feature_type]:
        mat_csc = feature_data[feature_type]['mat']

        if zscore_features:
            print('*** warning, z-scoring not supported in feature_data format')

        # drop barcodes with fewer than threshold UMI
        ser_sum = mat_csc.sum(axis=0)
        arr_sum = np.asarray(ser_sum[0,:])
        ser_sum = pd.Series(arr_sum[0])
        ser_keep = ser_sum[ser_sum >= min_umi]
        ser_keep = ser_keep[ser_keep <= max_umi]

        # these indexes will be used to filter all features
        keep_indexes = ser_keep.index.tolist()

        # filter barcodes
        barcodes = feature_data[feature_type]['barcodes']
        ser_barcodes = pd.Series(barcodes)
        barcodes_filt = ser_barcodes[keep_indexes].get_values()

        # return Dictionary of DataFrames
        filtered_data = {}
        for inst_feat in feature_data:

            inst_mat = feature_data[inst_feat]['mat']
            mat_filt = inst_mat[:, keep_indexes]
            feature_names = feature_data[inst_feat]['features']

            inst_data = {}
            inst_data['mat'] = mat_filt
            inst_data['barcodes'] = barcodes_filt
            inst_data['features'] = feature_names

            filtered_data[inst_feat] = inst_data

    else:
        # drop barcodes with fewer than threshold UMI
        inst_df = feature_data[feature_type]

        if zscore_features:
            print('z-scoring features')
            net.load_df(inst_df)
            net.normalize(axis='row', norm_type='zscore')
            inst_df = net.export_df()

        ser_sum = inst_df.sum(axis=0)
        ser_keep = ser_sum[ser_sum >= min_umi]
        ser_keep = ser_keep[ser_keep <= max_umi]
        keep_cols = ser_keep.index.tolist()

        # filter data
        filtered_data = {}
        for inst_feat in feature_data:

            filtered_data[inst_feat] = feature_data[inst_feat][keep_cols]

    return filtered_data

def convert_feature_data_to_df_dict(feature_data, make_sparse=True):

    # return Dictionary of DataFrames
    df = {}
    for inst_feat in feature_data:

        inst_mat = feature_data[inst_feat]['mat']
        feature_names = feature_data[inst_feat]['features']
        barcodes = feature_data[inst_feat]['barcodes']

        if make_sparse:
            inst_data = pd.SparseDataFrame(data=inst_mat, index=feature_names, columns=barcodes, default_fill_value=0)
        else:
            inst_data = pd.DataFrame(data=inst_mat.todense(), index=feature_names, columns=barcodes)

        df[inst_feat] = inst_data

    return df

def check_feature_data_size(feature_data):
    for inst_feat in feature_data:
        print(inst_feat)
        print(len(feature_data[inst_feat]['features']), len(feature_data[inst_feat]['barcodes']))
        print(feature_data[inst_feat]['mat'].shape, '\n')

def calc_mito_gene_umi_fraction(df_gex, meta_cell, plot_mito=False, mito_thresh=0.9):

    # Removing Mitochondrial Genes
    list_mito_genes = ['MTRNR2L11', 'MTRF1', 'MTRNR2L12', 'MTRNR2L13', 'MTRF1L', 'MTRNR2L6', 'MTRNR2L7',
                    'MTRNR2L10', 'MTRNR2L8', 'MTRNR2L5', 'MTRNR2L1', 'MTRNR2L3', 'MTRNR2L4']


    all_genes = df_gex.index.tolist()
    mito_genes = [x for x in all_genes if 'MT-' == x[:3] or
                 x.split('_')[0] in list_mito_genes]


    mito_sum = df_gex.loc[mito_genes].sum(axis=0)
    gex_sum = df_gex.sum(axis=0)

    mito_fraction = mito_sum/gex_sum

    if plot_mito:
        mito_fraction.sort_values(ascending=False).plot()

    list_mito_dead = []
    cells = mito_fraction.index.tolist()
    for inst_cell in cells:
        inst_mito = mito_fraction[inst_cell]
        if inst_mito >= mito_thresh:
            inst_state = 'dead-cell'
        else:
            inst_state = 'live-cell'
        list_mito_dead.append(inst_state)

    ser_dead = pd.Series(list_mito_dead, index=cells)

    meta_cell['mito-fraction-umi'] = mito_fraction
    meta_cell['dead-cell-mito'] = ser_dead

    return meta_cell

def set_hto_thresh(df_hto, meta_hto, hto_name, max_plot_hto=7, thresh=1, ylim =100):

    if 'hto-threshold' not in meta_hto.columns.tolist():
        ser_thresh = pd.Series(np.nan, index=meta_hto.index)
        meta_hto['hto-threshold'] = ser_thresh


    ser_hto = df_hto.loc[hto_name]

    n, bins, patches = plt.hist(ser_hto, bins=100, range=(0, max_plot_hto))

    colors = []
    for inst_bin in bins:
        if inst_bin <= thresh:
            colors.append('red')
        else:
            colors.append('blue')

    # apply the same color for each class to match the map
    for patch, color in zip(patches, colors):
        patch.set_facecolor(color)

    meta_hto.loc[hto_name, 'hto-threshold'] = thresh
    plt.ylim((0,100))

def ini_meta_cell(df):
    cells = df['gex'].columns.tolist()
    list_ser = []

    # look for available data types
    found_types = list(set(['gex', 'adt', 'hto']).intersection(df.keys()))
    for inst_type in found_types:

        # calc umi sum
        inst_ser = df[inst_type].sum(axis=0)
        inst_ser.name = inst_type + '-umi-sum'

        list_ser.append(inst_ser)

    df['meta_cell'] = pd.DataFrame(data=list_ser).transpose()

    df_gex = deepcopy(df['gex'])
    df_gex[df_gex >= 1] = 1
    ser_gene_num = df_gex.sum(axis=0)
    df['meta_cell']['num_expressed_genes'] = ser_gene_num

    return df

def meta_cell_gex_wo_mito_ribo(df_gex_ini, meta_cell):

    df_gex = deepcopy(df_gex_ini)

    # calc umi sum
    ser_umi_sum = df_gex.sum(axis=0)

    meta_cell['gex-umi-sum-no-ribo-mito'] = ser_umi_sum

    # count number of measured genes

    df_gex[df_gex >= 1] = 1
    ser_gene_num = df_gex.sum(axis=0)

    meta_cell['num_expressed_genes_no-ribo-mito'] = ser_gene_num

    return meta_cell

def ini_meta_gene(df_gex_ini):

    df_gex = deepcopy(df_gex_ini)

    # Mean UMI
    ser_gene_mean = df_gex.mean(axis=1)
    ser_gene_mean.name = 'mean'

    # Variance UMI
    ser_gene_var = df_gex.mean(axis=1)
    ser_gene_var.name = 'variance'

    # fraction of cells measured
    df_gex[df_gex >= 1] = 1
    ser_gene_meas = df_gex.sum(axis=1)/df_gex.shape[1]
    ser_gene_meas.name = 'fraction of cells measured'

    meta_gene = pd.concat([ser_gene_mean, ser_gene_var, ser_gene_meas], axis=1)

    return meta_gene

def assign_htos(df_hto_ini, meta_hto, meta_cell, sn_thresh):

    ser_list = []
    df_hto = deepcopy(df_hto_ini)

    for inst_row in df_hto.index.tolist():

        # get data for a HTO
        inst_ser = deepcopy(df_hto.loc[inst_row])

        # load threshold level for this HTO
        inst_thresh = meta_hto.loc[inst_row, 'hto-threshold']

        # binarize HTO values about threshold
        inst_ser[inst_ser < inst_thresh] = 0
        inst_ser[inst_ser >= inst_thresh] = 1

        # assemble list of series to make dataframe later
        ser_list.append(inst_ser)

        # find cells that are positive for this HTO
        pos_hto = inst_ser[inst_ser==1].index.tolist()


    # generate binarized dataframe
    df_binary = pd.concat(ser_list, axis=1).transpose()

    # find singlets
    ser_sum = df_binary.sum(axis=0)
    ct_list = {}
    ct_list['debris']    = ser_sum[ser_sum == 0].index.tolist()
    ct_list['singlet']   = ser_sum[ser_sum == 1].index.tolist()
    ct_list['multiplet'] = ser_sum[ser_sum > 1].index.tolist()

    # initialize dehash-thresh
    if 'dehash-thresh' not in meta_cell.columns.tolist():
        ser_type = pd.Series(np.nan, index=meta_cell.index)
        meta_cell['dehash-thresh'] = ser_type

    # save dehash-thresh
    for inst_type in ct_list:
        meta_cell.loc[ct_list[inst_type], 'dehash-thresh'] = inst_type

    # find the highest hto
    ser_max_hto = df_hto.idxmax(axis=0)
    meta_cell['hto-max'] = ser_max_hto

    # calc signal vs noise
    list_first = []
    list_second = []
    list_cells = []

    for inst_cell in df_hto.columns.tolist():
        inst_ser = df_hto[inst_cell].sort_values(ascending=False)
        inst_first  = inst_ser.get_values()[0]
        inst_second = inst_ser.get_values()[1]

        list_first.append(inst_first)
        list_second.append(inst_second)
        list_cells.append(inst_cell)

    ser_first  = pd.Series(data=list_first,  index=list_cells, name='first highest HTO')
    ser_second = pd.Series(data=list_second, index=list_cells, name='second highest HTO')

    df_comp = pd.concat([ser_first, ser_second], axis=1).transpose()

    sn_ratio = np.log2(df_comp.loc['first highest HTO']/df_comp.loc['second highest HTO'])

    # save to metadata
    meta_cell['hto-log2-sn'] = sn_ratio

    # assign tentative sample
    list_samples = []
    for inst_cell in meta_cell.index.tolist():

        inst_type = meta_cell.loc[inst_cell, 'dehash-thresh']

        if inst_type == 'singlet':
            inst_hto = meta_cell.loc[inst_cell, 'hto-max']
            inst_sample = meta_hto.loc[inst_hto]['Sample']
        else:
            inst_sample = 'N.A.'

        list_samples.append(inst_sample)

    ser_sample = pd.Series(list_samples, index=meta_cell.index.tolist())
    meta_cell['Sample-thresh'] = ser_sample

    # signal-to-noise adjustment
    cells = meta_cell.index.tolist()
    for inst_cell in cells:
        inst_type = meta_cell.loc[inst_cell, 'dehash-thresh']
        inst_sn = meta_cell.loc[inst_cell, 'hto-log2-sn']
        inst_max_hto = meta_cell.loc[inst_cell, 'hto-max']

        # change singlet to multiplet if low sn
        if inst_type == 'singlet':
            if inst_sn < sn_thresh['singlets']:
                # convert to multiplet
                meta_cell.loc[inst_cell, 'dehash-thresh-sn'] = 'multiplet'
                meta_cell.loc[inst_cell, 'Sample-thresh-sn'] = 'N.A.'
            else:
                meta_cell.loc[inst_cell, 'dehash-thresh-sn'] = 'singlet'
                meta_cell.loc[inst_cell, 'Sample-thresh-sn'] = meta_hto.loc[inst_max_hto]['Sample']
        elif inst_type == 'debris':
            if inst_sn >= sn_thresh['debris']:
                meta_cell.loc[inst_cell, 'dehash-thresh-sn'] = 'singlet'
                meta_cell.loc[inst_cell, 'Sample-thresh-sn'] = meta_hto.loc[inst_max_hto]['Sample']
            else:
                meta_cell.loc[inst_cell, 'dehash-thresh-sn'] = 'debris'
                meta_cell.loc[inst_cell, 'Sample-thresh-sn'] = 'N.A.'
        elif inst_type == 'multiplet':
            if inst_sn >= sn_thresh['multiplets']:
                meta_cell.loc[inst_cell, 'dehash-thresh-sn'] = 'singlet'
                meta_cell.loc[inst_cell, 'Sample-thresh-sn'] = meta_hto.loc[inst_max_hto]['Sample']
            else:
                meta_cell.loc[inst_cell, 'dehash-thresh-sn'] = 'multiplet'
                meta_cell.loc[inst_cell, 'Sample-thresh-sn'] = 'N.A.'

    return meta_cell

def plot_signal_vs_noise(df, alpha=0.25, s=10, hto_range=7, inf_replace=5):

    fig, axes = plt.subplots(nrows=1, ncols=2)

    list_first = []
    list_second = []
    list_cells = []
    for inst_cell in df.columns.tolist():
        inst_ser = df[inst_cell].sort_values(ascending=False)
        inst_first  = inst_ser.get_values()[0]
        inst_second = inst_ser.get_values()[1]

        list_first.append(inst_first)
        list_second.append(inst_second)
        list_cells.append(inst_cell)

    ser_first  = pd.Series(data=list_first,  index=list_cells, name='first highest HTO')
    ser_second = pd.Series(data=list_second, index=list_cells, name='second highest HTO')

    df_comp = pd.concat([ser_first, ser_second], axis=1).transpose()

    df_comp.transpose().plot(kind='scatter', figsize=(5,5),
                             x='first highest HTO', y='second highest HTO',
                             ylim=(0,hto_range), xlim=(0,hto_range), alpha=alpha, s=s, ax=axes[0])

    sn_ratio = np.log2(df_comp.loc['first highest HTO']/df_comp.loc['second highest HTO'])


    # replace positive infinities with set value
    sn_ratio = sn_ratio.replace(np.Inf, inf_replace)
    sn_ratio.hist(bins=100, ax=axes[1], figsize=(15,7))

    return df_comp, sn_ratio

def filter_ribo_mito_from_gex(df_ini, meta_cell):

    df = deepcopy(df_ini)
    all_genes = df.index.tolist()

    ribo_rpl = [x for x in all_genes if 'RPL' in x]
    ribo_rps = [x for x in all_genes if 'RPS' in x]
    ribo_genes = ribo_rpl + ribo_rps

    # calculate average ribo gene expression
    ser_ribo = df_ini.loc[ribo_genes].mean(axis=0)
    ser_ribo.name = 'Ribosomal-Avg'

    keep_genes = [x for x in all_genes if x not in ribo_genes]

    df = df.loc[keep_genes]

    # Removing Mitochondrial Genes
    list_mito_genes = ['MTRNR2L11', 'MTRF1', 'MTRNR2L12', 'MTRNR2L13', 'MTRF1L', 'MTRNR2L6', 'MTRNR2L7',
                    'MTRNR2L10', 'MTRNR2L8', 'MTRNR2L5', 'MTRNR2L1', 'MTRNR2L3', 'MTRNR2L4']


    all_genes = df.index.tolist()

    mito_genes = [x for x in all_genes if 'MT-' == x[:3] or
                 x.split('_')[0] in list_mito_genes]


    # calculate average ribo gene expression
    ser_mito = df_ini.loc[mito_genes].mean(axis=0)
    ser_mito.name = 'Mitochondrial-Avg'

    keep_genes = [x for x in all_genes if x not in mito_genes]

    df = df.loc[keep_genes]

    # calculate average ribo gene expression
    df_meta = pd.concat([ser_ribo, ser_mito], axis=1)

    meta_cell = pd.concat([meta_cell, df_meta], axis=1)

    return df, meta_cell

def add_cats_from_meta(barcodes, df_meta, add_cat_list):
    '''
    Add categories from df_meta.
    '''

    # get metadata of interest (add_cat_list) from barcodes of interest
    df_cats = df_meta.loc[barcodes][add_cat_list]

    # get list of cats
    list_cat_ini = [list(x) for x in df_cats.values]

    # add titles to cats
    list_cat_titles = [ list([str(x) + ': ' + str(y) for x,y in zip(add_cat_list, a)]) for a in list_cat_ini]

    # add barcodes to new columns
    new_cols = [tuple([x] + y) for x,y in zip(barcodes, list_cat_titles)]

    return new_cols


def make_dehash_meta_cell(df, ct_list, hto_names, ct_max_hto, sn_ratio_all, sn_singlet=1, sn_multiplet=1, sn_debris=1):

    inst_list = {}
    for inst_type in ['sample', 'hto-assigned', 'hto-max', 'sn-ratio']:
        inst_list[inst_type] = []

    rows = df['gex'].columns.tolist()
    for inst_col in df['gex'].columns.tolist():

        hto_cat = ''
        sample_cat = ''
        max_hto = ct_max_hto[inst_col]

        # signal-to-noise threshold
        inst_sn = sn_ratio_all[inst_col]

        # singlet
        if inst_col in ct_list['singlets']:

            if inst_sn >= sn_singlet:
                hto_cat = max_hto
                sample_cat = hto_names[hto_cat]

            # deemed a doublet based on signal-to-noise threshold
            else:
                hto_cat = 'multiplet-sn-ratio'
                sample_cat = 'N.A.'

        elif inst_col in ct_list['debris']:

            if inst_sn >= sn_debris:
                hto_cat = max_hto
                sample_cat = hto_names[hto_cat]
            else:
                hto_cat = 'debris'
                sample_cat = 'N.A.'

        elif inst_col in ct_list['multiplets']:

            # rescue multiplet to singlet if above signal-to-noise threshold
            if inst_sn >= sn_multiplet:
                hto_cat = max_hto
                sample_cat = hto_names[hto_cat]
            else:
                hto_cat = 'multiplet'
                sample_cat = 'N.A.'


        # save to lists
        inst_list['hto-assigned'].append(hto_cat)
        inst_list['hto-max'].append(max_hto)
        inst_list['sample'].append(sample_cat)
        inst_list['sn-ratio'].append(sn_ratio_all[inst_col])

    ser_list = []
    for inst_type in ['sample', 'hto-assigned', 'hto-max', 'sn-ratio']:
        inst_ser = pd.Series(inst_list[inst_type], name=inst_type, index=rows)
        ser_list.append(inst_ser)

    df_meta = pd.concat(ser_list, axis=1)


    return df_meta

def make_cyto_export(df, num_var_genes=500):

    keep_meta = ['hto-umi-sum',
                 'gex-umi-sum',
                 'adt-umi-sum',
                 'num_expressed_genes',
                 'hto-log2-sn',
                 'mito-fraction-umi',
                 'Ribosomal-Avg',
                 'Mitochondrial-Avg',
                 'gex-umi-sum-no-ribo-mito',
                 'num_expressed_genes_no-ribo-mito']

    df_cyto = None

    for inst_type in ['gex', 'adt', 'hto', 'meta_cell']:
        inst_df = deepcopy(df[inst_type])

        # filter for top var genes
        if inst_type == 'gex':
            keep_var_genes = inst_df.var(axis=1).sort_values(ascending=False).index.tolist()[:num_var_genes]
            inst_df = inst_df.loc[keep_var_genes]

        if 'meta' not in inst_type:
            inst_df.index = [inst_type.upper() + '_' + x for x in inst_df.index.tolist()]

        else:
            inst_df = inst_df[keep_meta].transpose()
            inst_df.index = ['DER_' + x for x in inst_df.index.tolist()]

        print(inst_type, inst_df.shape)

        if df_cyto is None:
            df_cyto = inst_df
        else:
            df_cyto = df_cyto.append(inst_df)

    df_export = df_cyto.transpose()

    cells = df_export.index.tolist()
    index_cells = [str(x/100) for x in range(len(cells))]
    df_export.index = index_cells

    ser_index = pd.Series(data=index_cells, index=cells)
    df['meta_cell']['Cytobank-Index'] = ser_index

    df['cyto-export'] = df_export

    return df

# # alternate lambda function
# def sum_field(dataframe, field):
#     return dataframe[field].sum(axis=0)

# list_ser_functions = {**{inst_type+'-umi-sum':(lambda y,inst_type=inst_type: sum_field(y,inst_type))\
#                          for inst_type in ['gex', 'adt', 'hto']},
#                      }

# for key,value in list_ser_functions.items():
#     list_ser.append(value(df))
# df['meta_cell'] = pd.DataFrame(data=list_ser).transpose()