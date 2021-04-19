#!/usr/bin/python
# Hernando M Vergara Mar 2021
# batch_process_CP_data.py

# processes several animals.
# See Inmuno_4channels_analysis notebook in CellProfiler_AnalysisPipelines repo

import os
import sys
import numpy as np
sys.path.append('../')
import CellProfiler_AnalysisPipelines.utils.data_reader as data_reader
import CellProfiler_AnalysisPipelines.utils.generic_functions as gf


# Indicate the animals ID you want to load
animal_list_to_load = ['PH301', 'PH302', 'PH303', 'PH304', 'PH305', 'PH306']
# Indicate the main path where all processed data is
processed_data_path = '/mnt/c/Users/herny/Desktop/SWC/Data/Microscopy_Data/Plasticity/PH3_inmuno/Processed_data/'
# Remove some manual ROIs that are not good
g_names_to_remove = ['PH301-1-5-R-Tail']
# define parameters for selection of ph3 cells
c4_thr = 0.75
c4_dif = 0.8
# select a criteria for what is an SPN
c2_thr = .7
# define parameters for selection of d2 cells
c3_thr_for_d1 = 0.43
c3_thr_for_d2 = 0.65
c3_dif_for_d2 = 0.95

# loop through
for animal_to_load in animal_list_to_load:
    print('---> Analysing animal {}'.format(animal_to_load))
    # get paths
    data_path = os.path.join(processed_data_path, animal_to_load)
    # path where the dataframe is (the output of CellProfiler)
    CPoutput_path = os.path.join(data_path, 'Cell_profiler_output/')
    # Import the data
    df = data_reader.PH3_data_reader(CPoutput_path)
    # Create a unique identifier for every instance of measure (manual ROI)
    df['group_name'] = df.apply(gf.group_name, axis=1)
    df['manual_roi_name'] = df.apply(gf.manual_roi_name, axis=1)
    # Add the data path as an attribute of the dataframe
    df.attrs['datapath'] = data_path
    # remove ROIs
    df = df[~df.manual_roi_name.isin(g_names_to_remove)]
    # find indexes of PH3 cells
    PH3_indexes = df[np.logical_and(df.I_cell_C4 >= c4_thr, (df.I_cell_C4 * c4_dif) >= df.I_surround_C4)].index.values
    # find indexes of spn cells
    SPN_indexes = df[df.I_cell_C2 >= c2_thr].index.values
    # find d1 and d2 indexes
    d1_indexes = df[df.I_cell_C3 <= c3_thr_for_d1].index.values
    d2_indexes = df[np.logical_and(df.I_cell_C3 >= c3_thr_for_d2,
                 (df.I_cell_C3 * c3_dif_for_d2) >= df.I_surround_C3)].index.values
    # combine them with the SPNs
    d1_SPN_indexes = np.intersect1d(SPN_indexes, d1_indexes)
    d2_SPN_indexes = np.intersect1d(SPN_indexes, d2_indexes)
    # percentage of d2 cells
    d2_num = len(d2_SPN_indexes)
    d1_num = len(d1_SPN_indexes)
    d2_perc = 100 * d2_num / (d1_num + d2_num)
    print('Percentage of d2 cells {}'.format(d2_perc))
    # total number of PH3 cells
    print("There are {} PH3+ cells".format(len(PH3_indexes)))
    # total number of d1+ PH3 cells
    d1_SPN_PH3_indexes = np.intersect1d(PH3_indexes, d1_SPN_indexes)
    print("There are {} d1 + cells".format(len(d1_SPN_PH3_indexes)))
    # total number of d2+ PH3 cells
    d2_SPN_PH3_indexes = np.intersect1d(PH3_indexes, d2_SPN_indexes)
    print("There are {} d2 + cells".format(len(d2_SPN_PH3_indexes)))
    # total number of spn PH3 cells
    SPN_PH3_indexes = np.intersect1d(PH3_indexes, SPN_indexes)
    print("There are {} spn + cells".format(len(SPN_PH3_indexes)))

    # save a dataframe with only the ph3 for the gradient
    df_out = df.loc[PH3_indexes]
    # create new column
    df_out['cell_label'] = 'undetermined'
    df_out.at[SPN_PH3_indexes, 'cell_label'] = 'spn'
    df_out.at[d1_SPN_PH3_indexes, 'cell_label'] = 'd1'
    df_out.at[d2_SPN_PH3_indexes, 'cell_label'] = 'd2'
    df_out.to_pickle(os.path.join(data_path, animal_to_load + '_ph3_cells.pkl'))

    # save the dataframe with the SPN information to calculate, grouping mice together,
    # the ph3 positive cells that are d1 or d2 (normalization would have to be done per mouse)
    # create columns
    df['SPN'] = 0
    df.at[SPN_indexes, 'SPN'] = 1
    df['d1'] = 0
    df.at[d1_SPN_indexes, 'd1'] = 1
    df['d2'] = 0
    df.at[d2_SPN_indexes, 'd2'] = 1
    df['PH3'] = 0
    df.at[PH3_indexes, 'PH3'] = 1
    # save indexes for going back if necessary
    df['df_indexes'] = df.index.values
    df.to_pickle(os.path.join(data_path, animal_to_load + '_all_cells.pkl'))

print('#### DONE')
