# Version: 0.1.0
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
    print('0.1.0')

def make_dir(directory):
    if not os.path.exists(directory):
        os.mkdir(directory)

def load_v3_comp_sparse_feat_matrix(inst_path):
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

def calc_mito_gene_umi_fraction(df, plot_mito=False):

    # Removing Mitochondrial Genes
    list_mito_genes = ['MTRNR2L11', 'MTRF1', 'MTRNR2L12', 'MTRNR2L13', 'MTRF1L', 'MTRNR2L6', 'MTRNR2L7',
                    'MTRNR2L10', 'MTRNR2L8', 'MTRNR2L5', 'MTRNR2L1', 'MTRNR2L3', 'MTRNR2L4']


    all_genes = df.index.tolist()
    mito_genes = [x for x in all_genes if 'MT-' == x[:3] or
                 x.split('_')[0] in list_mito_genes]


    mito_sum = df.loc[mito_genes].sum(axis=0)
    gex_sum = df.sum(axis=0)

    mito_fraction = mito_sum/gex_sum

    if plot_mito:
        mito_fraction.sort_values(ascending=False).plot()

    return mito_fraction

def plot_hto(df_hto, thresh_levels, hto_name, max_plot_hto=7, thresh=1, ylim =100):
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

    thresh_levels[hto_name] = thresh
    plt.ylim((0,100))

def assign_htos(df_hto_ini, thresh_levels):
    hto_key = {}
    ser_list = []
    df_hto = deepcopy(df_hto_ini)

    for inst_row in df_hto.index.tolist():

        # get data for a HTO
        inst_ser = df_hto.loc[inst_row]

        # load threshold level for this HTO
        inst_thresh = thresh_levels[inst_row]
        print(inst_row, 'thresh level: ', inst_thresh)

        # binarize HTO values about threshold
        inst_ser[inst_ser < inst_thresh] = 0
        inst_ser[inst_ser >= inst_thresh] = 1

        # assemble list of series to make dataframe later
        ser_list.append(inst_ser)

        # find cells that are positive for this HTO
        pos_hto = inst_ser[inst_ser==1].index.tolist()
        hto_key[inst_row] = pos_hto
        print('num hto pos: ', len(hto_key[inst_row]))


    # generate binarized dataframe
    df_binary = pd.concat(ser_list, axis=1).transpose()

    # find singlets
    ser_sum = df_binary.sum(axis=0)
    ct_list = {}
    ct_list['debris']     = ser_sum[ser_sum == 0].index.tolist()
    ct_list['singlets']   = ser_sum[ser_sum == 1].index.tolist()
    ct_list['multiplets'] = ser_sum[ser_sum > 1].index.tolist()


    # generate cell type dict
    ct_dict = {}
    print('\noverview')
    for inst_type in ct_list:
        print(inst_type, ': ', len(ct_list[inst_type]))
        for inst_cell in ct_list[inst_type]:
            ct_dict[inst_cell] = inst_type

    return hto_key, ct_list, ct_dict

def calc_first_vs_second(df, alpha=0.25, s=10, hto_range=7):
    list_first = []
    list_second = []
    list_cells = []
    for inst_cell in df.columns.tolist():
        inst_ser = df[inst_cell].sort_values(ascending=False)
        # print(inst_ser)
        inst_first  = inst_ser.get_values()[0]
        inst_second = inst_ser.get_values()[1]
        # print(inst_first, inst_second)

        list_first.append(inst_first)
        list_second.append(inst_second)
        list_cells.append(inst_cell)

    ser_first  = pd.Series(data=list_first,  index=list_cells, name='first highest HTO')
    ser_second = pd.Series(data=list_second, index=list_cells, name='second highest HTO')

    df_comp = pd.concat([ser_first, ser_second], axis=1).transpose()

    df_comp.transpose().plot(kind='scatter', figsize=(5,5),
                             x='first highest HTO', y='second highest HTO',
                             ylim=(0,hto_range), xlim=(0,hto_range), alpha=alpha, s=s)

    return df_comp


def filter_ribo_mito_from_gex(df_ini):

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
    df_meta = pd.concat([ser_ribo, ser_mito], axis=1).transpose()


    return df, df_meta

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
