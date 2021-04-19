#!/usr/bin/python
# Hernando M Vergara Mar 2021
# merge_csvs.py

import sys
import glob
import os
import pandas as pd


def merge_csv_tables(main_path):
    # get the list of the processed mice in the folder
    csv_files = glob.glob(os.path.join(main_path, '**', '*_ARA_coordinates.csv'))
    # create a list to concatenate dataframes
    dfs_list = []

    for csv_file in csv_files:
        # read them as dataframe
        df_sub = pd.read_csv(csv_file)
        # add the name of the mouse
        df_sub['animal_id'] = csv_file.split(os.sep)[-2]
        # append to list
        dfs_list.append(df_sub)

    # concatenate to a dataframe
    df = pd.concat(dfs_list, ignore_index=True)
    # write main df as csv
    file_out_path = os.path.join(main_path, 'ph3_cells_ARA_coordinates.csv')
    df.to_csv(file_out_path, index=False)


if __name__ == '__main__':
    # check input
    if len(sys.argv) not in [2]:
        sys.exit('Incorrect number of arguments, please run like this:\
            python merge_csvs.py path_to_processed_data_folder')
    # catch input
    inpath = sys.argv[1]
    # run function
    merge_csv_tables(main_path=inpath)
