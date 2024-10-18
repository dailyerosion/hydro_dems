# hydro_dems
a variety of arcpy scripts that create digital elevation models designed for enhancing surface water flow

The .pyt scripts are designed to be run with ArcGIS Pro as Python toolboxes. They can also be run via the command line by removing the comments from the 'if __name__ == __main__' code block (from there to the end). This code breaks a Python toolbox.

# First
Start with the script 'get_merge_lidar_datasets.pyt'. This will download a couple datasets and intersect them for future work. The EPT/WESM output feature class from this is used in the next step.

# Second
Next, use the cmd_build_DEM.pyt with the EPT/WESM output from above and a polygon of the area you want to download. You also need the pdal.exe as in put (below). It will build a 'pit-filled' DEM and other optional output if you desire. 
## Required
You MUST also install pdal to make this script work. I suggest installing Anaconda and then adding PDAL to it as described here.  
https://pdal.io/en/2.7-maintenance/quickstart.html#install-conda
Then you need to find the path to the pdal.exe file to include in your arguments. That should look somethig like:
C:\Users\bkgelder\Anaconda3\envs\pdal\Library\bin\pdal.exe
but that could vary a lot depending on how you installed Anaconda (all users or just you).

# Optional next steps
A couple optional DEM re-interpolation processes follow that. First is cmd_flattener_DEM.pyt and then cmd_cleaner_DEM.pyt

# Third (or Fourth or Fifth)
Next is 'cmd_puncher_DEM.pyt'. This will 'punch' holes in a DEM based on the depression depth and/or area, leaving only depressions that are deeper/larger than the criteria.
