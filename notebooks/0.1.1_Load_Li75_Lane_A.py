#!/usr/bin/env python
# coding: utf-8

# # 0.1.0 Load Lw5 Lane 1

# In[1]:


import numpy as np
import pandas as pd
from clustergrammer2 import net
from copy import deepcopy
df = {}


# ### Lane

# In[258]:





# In[259]:


import sys
inst_lane = (sys.argv[1] if len(sys.argv)>1 else '1')
base_dir='..'
base_dir='/home/manu/cellranger_base_reps/liver_lw5'
sys.path.insert(0, base_dir+'/himc_helper_functions/notebooks/')
import himc_helper_functions as hf


# In[260]:


feature_data = hf.load_crv3_feature_matrix( base_dir+'/data/primary_data/LW5-2-'+ inst_lane +'_raw_feature_bc_matrix/')


# In[261]:


hf.check_feature_data_size(feature_data)


# ### GEX

# In[262]:


ser_sum = hf.plot_umi_levels(feature_data, feature_type='adt', min_umi=1)


# In[263]:


feature_filtered = hf.filter_barcodes_by_umi(feature_data, 'adt', min_umi=1)


# ### Filter barcodes based on min GEX UMI level

# ### Convert to Dense

# In[264]:


df = hf.convert_feature_data_to_df_dict(feature_filtered, make_sparse=False)
print('adt', df['adt'].shape)

# add lane to barcode
df['adt'].columns = [x + '-' + inst_lane for x in df['adt'].columns.tolist()]


# In[265]:


df['adt']


# ### Separate ADT and HTO

# In[266]:


rows = df['adt'].index.tolist()
adt_rows = [x for x in rows if 'adt_' in x]
hto_rows = [x for x in rows if 'HTO' in x]
initial_adt=deepcopy(df['adt'])
df['adt'] = initial_adt.loc[adt_rows]
df['hto'] = initial_adt.loc[hto_rows]

print(df['adt'].shape, df['hto'].shape)


# ### Make ADT Metadata

# In[267]:


rows = df['adt'].index.tolist()
rows = [x.replace('adt_5', '') for x in rows]
df['adt'].index = rows


# In[268]:


ser_chem = pd.Series(['5-prime']*len(rows), name='Chemistry', index=rows)
df['meta_adt'] = pd.concat([ser_chem], axis=1)
df['meta_adt'].shape


# ### Make HTO Metadata

# In[269]:


rows = df['hto'].index.tolist()
rows = [x.split('_')[-1].replace('HTO', 'hto-') for x in rows]
df['hto'].index = rows


# In[270]:


rows = df['hto'].index.tolist()
ser_chem = pd.Series(['5-prime']*len(rows), name='Chemistry', index=rows)
df['meta_hto'] = pd.concat([ser_chem], axis=1)
df['meta_hto'].shape


# In[271]:


df['meta_hto']


# ### Metadata Cells

# In[272]:


import importlib
importlib.reload(hf)
df = hf.ini_meta_cell(df)


# In[273]:


df['meta_cell'].head()


# ### Load Hashtag Sample Data

# In[274]:


df['hto-key'] = pd.read_csv( base_dir+'/data/primary_data/dehash_key.txt', sep='\t', index_col=0)


# In[275]:


df['hto-key'].index


# In[276]:


df['hto-key'].index = ['hto-' + x.split('HTO')[-1] for x in df['hto-key'].index]


# In[277]:


df['hto-key']


# In[278]:


df['meta_hto']


# In[279]:


df['meta_hto'] = pd.concat([df['meta_hto'], df['hto-key']], axis=1)


# In[280]:


df['meta_hto']


# In[281]:


df_hto = deepcopy(df['hto'])
df_hto = np.arcsinh(df_hto/5)
df_hto.index.tolist()


# In[282]:


# net.load_df(df_hto)
# net.widget()


# In[283]:


df['meta_hto']


# ### Threshold HTOs

# In[284]:


hf.set_hto_thresh(df_hto, df['meta_hto'], 'hto-21', thresh=2.2)


# In[285]:


hf.set_hto_thresh(df_hto, df['meta_hto'], 'hto-22', thresh=1.6)


# In[286]:


hf.set_hto_thresh(df_hto, df['meta_hto'], 'hto-23', thresh=1.6)


# In[287]:


hf.set_hto_thresh(df_hto, df['meta_hto'], 'hto-25', thresh=2)


# In[288]:


hf.set_hto_thresh(df_hto, df['meta_hto'], 'hto-26', thresh=2)


# In[289]:


hf.set_hto_thresh(df_hto, df['meta_hto'], 'hto-27', thresh=2.3)


# In[290]:


hf.set_hto_thresh(df_hto, df['meta_hto'], 'hto-28', thresh=2)


# In[291]:


hf.set_hto_thresh(df_hto, df['meta_hto'], 'hto-29', thresh=2.2)


# In[292]:


hf.set_hto_thresh(df_hto, df['meta_hto'], 'hto-30', thresh=2.2)


# In[293]:


hf.set_hto_thresh(df_hto, df['meta_hto'], 'hto-44', thresh=2.4)


# In[294]:


df['meta_hto']


# In[295]:


sn_thresh = {}
sn_thresh['singlets'] = 1
sn_thresh['debris'] = 1
sn_thresh['multiplets'] = 1


# In[296]:


# df['meta_cell'] = hf.assign_htos(df_hto, df['meta_hto'], df['meta_cell'], sn_thresh)


# In[297]:


df['meta_cell'].head()


# ### Signal-to-Noise Visuals
# Here we're defining signal as the highest hashtag and the 'noise' as the second highest hashtag. I prefer this name to something like "first highest hashtag vs second highest hashtag".

# ### All Cells

# In[298]:


# df_comp, sn_ratio_all = hf.plot_signal_vs_noise(df_hto)


# ### Debris

# In[299]:


# df_comp, sn_ratio = plot_signal_vs_noise(df_hto[ct_list['debris']])


# ### Multiplets

# In[300]:


# df_comp, sn_ratio = plot_signal_vs_noise(df_hto[ct_list['multiplet']])


# ### Singlets

# In[301]:


# df_comp, sn_ratio = plot_signal_vs_noise(df_hto[ct_list['singlet']])


# In[302]:


df.keys()


# ### Overview of Dehashing

# In[303]:


# df['meta_cell']['dehash-thresh'].value_counts()


# In[304]:


# df['meta_cell']['dehash-thresh-sn'].value_counts()


# In[305]:


# df['meta_cell']['Sample-thresh-sn'].value_counts()


# In[306]:


df['meta_cell'].head()


# ### Filter out Dead (high mito fraction) cells

# In[307]:


df['meta_cell'].head()


# ### Drop Ribo and Mito Genes and Calc Avg Exp of Ribo and Mito
# Permanently filter

# In[308]:


df['meta_cell'].head()


# ### More GEX Metadata
# After we have dropped mito and ribo genes

# ### Save Parquets

# In[312]:


df.keys()


# In[310]:


hf.make_dir( base_dir+'/data/processed_data/')
hf.make_dir( base_dir+'/data/processed_data/individual_lanes/')
hf.make_dir( base_dir+'/data/processed_data/individual_lanes/lane-' + inst_lane)


# In[311]:

df = hf.make_cyto_export(df)
hf.make_dir( base_dir+'/data/processed_data/cytobank_export/lane-'+inst_lane+'/')
df['cyto-export'].to_csv( base_dir+'/data/processed_data/cytobank_export/lane-'+inst_lane+'/LW5-2-'+inst_lane+'_cytobank_export.csv')


# In[ ]:




