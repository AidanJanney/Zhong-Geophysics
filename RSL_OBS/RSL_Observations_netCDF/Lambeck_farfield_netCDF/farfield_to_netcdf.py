import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd

## Open File as a dataframe
def file_to_dataframe(file_name, sheet_num=0):
    """Open file as dataframe and check execution

    Args:
        file_name (string): file name
        sheet_num (integer): index of individual sheet to retrieve
    """
    try:
        df = pd.read_excel(file_name, sheet_name=sheet_num)
    except:
        raise Exception(f"EXCEPTION: File Not Found -> Dataframe not created for {file_name}")
    else:
        return df
        
#Establish all dataframes
RSL_Data_df = file_to_dataframe('Lambeck_Farfield_Adjusted.xlsx', 0)
References_df = file_to_dataframe('Lambeck_Farfield_Adjusted.xlsx', 1)
Additional_Info_df = file_to_dataframe('Lambeck_Farfield_Adjusted.xlsx', 2)
print("Dataframes Opened")

## Cleaning up RSL_Data_df, had some duplicate rows, don't need observation numbers
## See notebook for removed/merged items
RSL_Data_df.drop(['Observation Number'], axis = 1, inplace=True)
RSL_Data_df.drop_duplicates(inplace=True)
print("Initial Processing: Duplicates Dropped")

## Group different sets togehter so we can remove duplicate indices -> this lets us provide duplicate coordinate
## Entries as a list of results to access later
RSL_Adj = RSL_Data_df.groupby(['Site Location', 'Latitude (degree)', 'Longitude (degree)', 'Age (ka)', 'Sample Type'])['RSL (m)'].apply(list).reset_index(name = 'RSL (m)')
Sigma_RSL_Adj = RSL_Data_df.groupby(['Site Location', 'Latitude (degree)', 'Longitude (degree)', 'Age (ka)', 'Sample Type'])['Sigma RSL (m)'].apply(list).reset_index(name = 'Sigma RSL (m)')
References_Adj = RSL_Data_df.groupby(['Site Location', 'Latitude (degree)', 'Longitude (degree)', 'Age (ka)', 'Sample Type'])['Reference Number'].apply(list).reset_index(name = 'Reference Number')
# Combine them all together now that each column RSL, Sigma RSL, and Reference Number has been reconciled
RSL_Data_Adj = (RSL_Adj.merge(Sigma_RSL_Adj).merge(References_Adj))[RSL_Data_df.columns.to_list()]
print("Grouping of RSL and Uncertanties Complete")


## Cleanup Reference Number
def adj_ref_list(ref_nums):
    ref_nums_adj = ref_nums.to_list().copy()
    for i, elem in enumerate(ref_nums):
        elem_temp = elem.copy()
        for j, entry in enumerate(elem):
            if isinstance(entry, str):
                elem_temp[j] = [int(num) for num in entry.split('/')]
        ref_nums_adj[i] = elem_temp
    return ref_nums_adj
# Apply the adj_ref_list function, does not edit in place (though that was a pain to figure out haha!)
RSL_Data_Adj['Reference Number'] = adj_ref_list(RSL_Data_Adj['Reference Number'])
print("Reference Numbers Cleaned Up")

## Set the index to the elements that we would like and save as netcdf file
print("Converting to Xarray")
RSL_Data = RSL_Data_Adj.set_index(['Site Location', 'Latitude (degree)', 'Longitude (degree)', 'Age (ka)']).to_xarray()
print("Converted to Xarray")
RSL_Data.to_netcdf('testing_this.nc')
print("Saved to netCDF")