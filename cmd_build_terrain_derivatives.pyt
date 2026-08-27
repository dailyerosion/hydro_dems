# -*- coding: utf-8 -*-
"""Creates DEMs and derivatives (return counts, intensity) from lidar datasets
available via the Entwine Point cloud format for a given polygon (usually 
buffered HUC12 watershed boundary). Uses a previously created merge of
the USGS WESM data and the EPT GeoJSON dataset to locate the AWS bucket for each
project, then joins multiple projects and creates a lidar DEM and derivatives
at one or multiple resolutions. Also creates a feature class of lidar datasets
used to generate the DEM and derivatives."""

import arcpy
from arcpy.sa import *
import arcpy.metadata as md
import sys
import os
import time
import subprocess
import pdal
import platform
import getpass
import glob
import traceback
import re
import shutil
import difflib
from os.path import join as opj
from pathlib import Path
from math import ceil

sys.path.append("C:\\DEP\\Scripts\\basics")

import dem_functions as df

#https://stackoverflow.com/questions/7006238/how-do-i-hide-the-console-when-i-use-os-system-or-subprocess-call
CREATE_NO_WINDOW = 0x08000000

class msgStub:
    def addMessage(self,text):
        arcpy.AddMessage(text)
    def addErrorMessage(self,text):
        arcpy.AddErrorMessage(text)
    def addWarningMessage(self,text):
        arcpy.AddWarningMessage(text)

class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Toolbox"
        self.alias = "toolbox"

        # List of tool classes associated with this toolbox
        self.tools = [Tool]


class Tool(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Build DEMs from LiDAR by EPT download"
        self.description = "Must install PDAL and run get_merge_lidar_datasets first to merge EPT and WESM. Build a DEM for a polygon by downloading the data through EPT format via PDAL."
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""

        param0 = arcpy.Parameter(
            name="monthly_ept_wesm_mashup",
            displayName="EPT WESM Merged Features From Previous",
            datatype="DEFeatureClass",
            parameterType='Required',
            direction="Input")
        
        param1 = arcpy.Parameter(
            name="dem_polygon",
            displayName="Polygon Feature (single) of DEM Area",
            datatype="DEFeatureClass",
            parameterType='Required',
            direction="Input")
        
        param2 = arcpy.Parameter(
            name = "pdal_exe",
            displayName="PDAL.exe Location",
            datatype="DEFile",
            parameterType='Required',
            direction="Input")
        
        param3 = arcpy.Parameter(
            name = "gsds",
            displayName="Integer Resolution/Ground Sample Distance (in meters) of rasters, multiples joined by comma",
            datatype="GPString",
            parameterType='Required',
            direction="Input")
        param3.values = "3,2,1"#default gsds value to create 3, 2, and 1 meter rasters
        
        param4 = arcpy.Parameter(
            name = "procDir",
            displayName="Local Processing Directory",
            datatype="DEFolder",
            parameterType='Optional',
            direction="Input")
        
        param5 = arcpy.Parameter(
            name="snap",
            displayName="Snap Raster",
            datatype="DERasterDataset",
            parameterType='Optional',
            direction="Input")
        
        param6 = arcpy.Parameter(
            name = "breakpolys",
            displayName="Input HUC12 Merged Breakline Polygon Features",
            datatype="DEFeatureClass",
            parameterType='Optional',
            direction="Input")
        
        param7 = arcpy.Parameter(
            name = "breaklines",
            displayName="Input HUC12 Merged Breakline Polyline Features",
            datatype="DEFeatureClass",
            parameterType='Optional',
            direction="Input")
        
        param8 = arcpy.Parameter(
            name = "tElevFile",
            displayName="Output Pit-Filled Elevation Model",
            datatype="DERasterDataset",
            parameterType='Required',
            direction="Output")
        
        param9 = arcpy.Parameter(
            name = "bareEarthReturnMinFile",
            displayName="Output Bare Earth Minimum Elevation Model",
            datatype="DERasterDataset",
            parameterType='Optional',
            direction="Output")
        
        param10 = arcpy.Parameter(
            name = "firstReturnMaxFile",
            displayName="Output First Return Maximum Elevation/Surface Model",
            datatype="DERasterDataset",
            parameterType='Optional',
            direction="Output")
        
        param11 = arcpy.Parameter(
            name = "cntBeFile",
            displayName="Output Bare Earth Return Count Raster",
            datatype="DERasterDataset",
            parameterType='Optional',
            direction="Output")
        
        param12 = arcpy.Parameter(
            name = "cnt1rFile",
            displayName="Output First Return Count Raster",
            datatype="DERasterDataset",
            parameterType='Optional',
            direction="Output")
        
        param13 = arcpy.Parameter(
            name = "cntPlsFile",
            displayName="Output Pulse Count Raster",
            datatype="DERasterDataset",
            parameterType='Optional',
            direction="Output")
        
        param14 = arcpy.Parameter(
            name = "int1rMinFile",
            displayName="Output Intensity First Return Minimum Raster",
            datatype="DERasterDataset",
            parameterType='Optional',
            direction="Output")
        
        param15 = arcpy.Parameter(
            name = "int1rMaxFile",
            displayName="Output Intensity First Return Maximum Raster",
            datatype="DERasterDataset",
            parameterType='Optional',
            direction="Output")
        
        param16 = arcpy.Parameter(
            name = "intBeMaxFile",
            displayName="Output Intensity Bare Earth Maximum Raster",
            datatype="DERasterDataset",
            parameterType='Optional',
            direction="Output")
        
        param17 = arcpy.Parameter(
            name = "ept_wesm_project_file",
            displayName="EPT WESM Feature for AOI",
            datatype="DEFeatureClass",
            parameterType='Optional',
            direction="Output")
        
        param18 = arcpy.Parameter(
            name = "lidar_download_directory",
            displayName="Lidar Data Download Directory",
            datatype="DEFolder",
            parameterType='Optional',
            direction="Output")
        
        param19 = arcpy.Parameter(
            name = "cleanup",
            displayName="end of run data deletion",
            datatype="GPBoolean",
            parameterType='Optional',
            direction="Output")

        # param18.filter.type = "ValueList"
        # param18.filter.list = [True, False]
        # param19.value = True
                        
        params = [param0, param1, param2, param3,
                  param4, param5, param6, param7,
                  param8, param9, param10, param11,
                  param12, param13, param14, param15,
                  param16, param17, param18, param19]
        return params


    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        cleanup = False
        # doLidarDEMs(parameters[0].valueAsText, parameters[1].valueAsText, parameters[2].valueAsText, parameters[3].valueAsText, parameters[4].valueAsText, 
        #             parameters[5].valueAsText, parameters[6].valueAsText, parameters[7].valueAsText, parameters[8].valueAsText, parameters[9].valueAsText, 
        #             parameters[10].valueAsText, parameters[11].valueAsText, parameters[12].valueAsText, parameters[13].valueAsText, parameters[14].valueAsText, 
        #             parameters[15].valueAsText, parameters[16].valueAsText, parameters[17].valueAsText, parameters[18].valueAsText, parameters[19].valueAsText, messages)
        return

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return

##----------------------------------------------------------------------
## Set environments and begin

# def fillOCSinks(inDEM, log):
#     # Return a raster will all one-cell-sinks filled
# ##    arcpy.AddMessage("-----Find Pits...")
#     log.info('finding flow direction')
#     sinkFDir = FlowDirection(inDEM)
#     log.info('finding sinks')
#     allSinks = Sink(sinkFDir)
#     # arcpy.BuildRasterAttributeTable_management(allSinks)
#     # log.info('sinks for ' + str(inDEM) + ' is ' + str(int(arcpy.GetCount_management(allSinks).getOutput(0))))
# ##    arcpy.AddMessage("-----Fill everything else...")

#     ## Make a No-one-cell-sink DEM
#     log.info('finding all but sinks')
#     AllButSinks_DEM = Con(IsNull(allSinks), inDEM)

#     log.info('filling pits')
#     ## Fill the No-one-cell-sink DEM
#     absDEM_fill = Fill(AllButSinks_DEM)
#     log.info('filled DEM')

#     ## Add the Original 'real' sinks back into the filled DEM
#     fill_DEM = Con(IsNull(absDEM_fill), inDEM, absDEM_fill)
#     log.info('fixed pits')

#     return(fill_DEM, allSinks)


def getfields(infc, fieldString = '', fieldType = ''):
    flds = arcpy.ListFields(infc, fieldString, fieldType)
    if len(flds) > 0:
        fieldNames = [fld.name for fld in flds]
    else:
        fieldNames = []
    return fieldNames


def processEptLas(sgdb, sfldr, srOutCode, fixedFolder, geom, ept_las, srOut, inm, FDSet, procDir, allTilesList, log, time, work_id):
    '''process a cursor row of data by creating a suitable las file from the input las/laz/zlas dataset
    This inlcudes project and clipping las data files into output dataset and also creating a multipoint
    file from the las data if there is any within the extent'''
    try:
        ept_las_base = os.path.splitext(os.path.basename(ept_las))[0]
        sfx = arcpy.ValidateTableName('_' + ept_las_base, sgdb)
        log.debug('lidar file suffix is: ' + sfx)

        allLasd = arcpy.CreateLasDataset_management(ept_las, opj(sfldr, 'all' + sfx))

        # extract to tile geometry and project if necessary
        nameSfx = '_' + str(srOutCode)
        fixedLasBasename = os.path.basename(ept_las)[:-4] + nameSfx + '.las'
    ##            log.debug('fixedLasBasename: ' + fixedLasBasename)
        # some old 3DEP projects don't alway have boundaries and data lining up...
        log.debug(f'ExtractLas arguments are {allLasd}, {fixedFolder}, name_suffix = {nameSfx}, rearrange_points = {"MAINTAIN_POINTS"}, out_las_dataset = {opj(fixedFolder, "fixed" + sfx + ".lasd")}')
        if work_id < 0:
            tileGeomBuffer5 = geom.buffer(5) #tile Geometry column
            log.debug('--- Use ExtractLas, boundary option,')
            fixedLasd = arcpy.ExtractLas_3d(allLasd, fixedFolder, name_suffix = nameSfx, rearrange_points = "MAINTAIN_POINTS", out_las_dataset = opj(fixedFolder, 'fixed' + sfx + '.lasd'), boundary = tileGeomBuffer5)
        else:
            log.debug('--- Use ExtractLas, no boundary option,')
            fixedLasd = arcpy.ExtractLas_3d(allLasd, fixedFolder, name_suffix = nameSfx, rearrange_points = "MAINTAIN_POINTS", out_las_dataset = opj(fixedFolder, 'fixed' + sfx + '.lasd'))#, boundary = tileGeomBuffer5)
        log.debug(fixedLasd.getMessages())
        fixedLasdDescDa = arcpy.da.Describe(fixedLasd)
        fixedLasPath = opj(fixedLasdDescDa['path'], fixedLasBasename)

        log.debug('--- Done creating LAS dataset and extracting LAS at ')

        if fixedLasdDescDa['pointCount'] > 0:
            allTilesList.append(fixedLasPath)

        if fixedLasPath in allTilesList:#non 0 amount of lidar points in las
            allLAZ = fixedLasPath
            if allLAZ.endswith('.laz') or allLAZ.endswith('.las'):
                """Filters LAS points to class 2 and creates multipoints in FDSet"""
                lasBase = os.path.splitext(os.path.basename(allLAZ))[0]
                if allLAZ.endswith('.laz'):
                    log.debug('--- Using ConvertLas to decompress LAZ')
                    las_from_laz = arcpy.ConvertLas_conversion(allLAZ, target_folder=procDir, compression=None, las_options=None)
                else:
                    las_from_laz = allLAZ

                log.debug('--- Create las non-Minnesota Multipoint')
                lasMP = arcpy.LASToMultipoint_3d(las_from_laz, inm + "\\pts" + sfx, "1", class_code = [2,8], input_coordinate_system = srOut)

            elif allLAZ.endswith('.zlas'):
                log.debug('--- Create zlas non-Minnesota Multipoint')
                lasMP = arcpy.LASToMultipoint_3d(allLAZ, inm + "\\pts" + sfx, "1", class_code = [2,8], input_coordinate_system = srOut)
                las_from_laz = allLAZ

            if lasMP:
                ptsName = arcpy.ValidateTableName('pts_' + lasBase, os.path.join(str(FDSet)))
                ptOut = projIfNeeded(lasMP, os.path.join(str(FDSet), ptsName), srOut)

            else:
                log.warning('no ptOut created, setting to None')
                ptOut = None
                
            cl2Las = las_from_laz

        # ready so there is something to return
        if 'cl2Las' not in locals():
            cl2Las = None
        if 'ptOut' not in locals():
            ptOut = None
        log.info('ptOut: ' + str(ptOut) + ' and cl2Las ' + str(cl2Las))

        return cl2Las

    except Exception as e:
        print('handling as exception')
##        log.debug(e.message)
        if sys.version_info.major == 2:
            arcpy.AddError(e.message)
            print(e.message)
        elif sys.version_info.major == 3:
            arcpy.AddError(e)
            print(e)

        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]

        # Concatenate information together concerning the error into a message string
        pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
        # Return python error messages for use in script tool or Python Window
        arcpy.AddError(pymsg)
        # Print Python error messages for use in Python / Python Window
        print(pymsg + "\n")

        if arcpy.GetMessages(2) not in pymsg:
            msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
            arcpy.AddError(msgs)
            print(msgs)

    except:
        print('handling as except')
        # Get the traceback object
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]

        # Concatenate information together concerning the error into a message string
        pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
        # Return python error messages for use in script tool or Python Window
        arcpy.AddError(pymsg)
        # Print Python error messages for use in Python / Python Window
        print(pymsg + "\n")

        if arcpy.GetMessages(2) not in pymsg:
            msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
            arcpy.AddError(msgs)
            print(msgs)

def prepPolygonBoundary(dem_polygon, log, sgdb, srOut, srSfx, maskRastBase, demLists):

    try:
        assert int(arcpy.GetCount_management(dem_polygon).getOutput(0)) < 2, 'multiple features in polygon feature class'
        assert int(arcpy.GetCount_management(dem_polygon).getOutput(0)) > 0, 'no features in polygon feature class'
        maskFc = arcpy.CopyFeatures_management(dem_polygon, opj(sgdb, 'maskFc'))
        maskFc_area = [s[0] for s in arcpy.da.SearchCursor(maskFc, ['SHAPE@AREA'])][0]

        # # geom_copy = arcpy.management.CopyFeatures(huc12fc, opj(sgdb, 'huc' + huc12))
        # geom_copy = arcpy.Buffer_analysis(maskFc, buffer_distance_or_field = '-1000 METERS')
        # if df.testForZero(geom_copy):
        #     geom_copy = arcpy.Buffer_analysis(maskFc, buffer_distance_or_field = '-500 METERS')
        #     if df.testForZero(geom_copy):
        #         geom_copy = arcpy.CopyFeatures_management(dem_polygon)

        if 'id' not in df.getfields(maskFc):
            arcpy.AddField_management(maskFc, 'id', 'LONG')
            arcpy.CalculateField_management(maskFc, 'id', 1, 'PYTHON')

        log.info("maskFcPrelim complete at: ")

        ## Set up geodatabase to store the multipoint files and terrains (necessary all inputs be in feature dataset
        # Vertical units are in meters (float) so use a meter-based reference
        FDSet = arcpy.CreateFeatureDataset_management(sgdb, "Lidar_pts", srOut)
        maskFcOut = projIfNeeded(maskFc, os.path.join(str(FDSet), 'buf_huc' + srSfx), srOut)
        log.info("maskFcOut complete at: ")

        for demList in demLists:
            maskRastOut = arcpy.PolygonToRaster_conversion(maskFcOut, 'id', opj(sgdb, maskRastBase + str(demList[0])), cellsize = demList[0])
            # huc_rast_out = arcpy.conversion.PolygonToRaster(geom_copy, 'OBJECTID', opj(sgdb, 'huc_rast' + str(demList[0])), cellsize = demList[0])

        return maskFc, maskFc_area, maskFcOut, maskRastOut, None, FDSet
        # return maskFc, maskFc_area, maskFcOut, maskRastOut, huc_rast_out, FDSet

    except Exception as e:
        print('handling as exception')
##        log.debug(e.message)
        if sys.version_info.major == 2:
            print('handling as 2 exception')
            arcpy.AddError(e.message)
            print(e.message)
        elif sys.version_info.major == 3:
            print('handling as 3 exception')
            arcpy.AddError(e)
            print(e)

        # Get the traceback object
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        # Concatenate information together concerning the error into a message string
        pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
        # Return python error messages for use in script tool or Python Window
        arcpy.AddError(pymsg)
        # Print Python error messages for use in Python / Python Window
        print(pymsg + "\n")

        if arcpy.GetMessages(2) not in pymsg:
            msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
            arcpy.AddError(msgs)
            print(msgs)

    except:
        print('handling as except')
        # Get the traceback object
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]

        # Concatenate information together concerning the error into a message string
        pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
        # Return python error messages for use in script tool or Python Window
        arcpy.AddError(pymsg)
        # Print Python error messages for use in Python / Python Window
        print(pymsg + "\n")

        if arcpy.GetMessages(2) not in pymsg:
            msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
            arcpy.AddError(msgs)
            print(msgs)


def getLidarTimeframes(merged):#, tilesClip_local):

    try:
    # if merged is not None:#in locals():
        collect_starts = [s[0] for s in arcpy.da.SearchCursor(merged, ['collect_start'])]
        collect_starts_min = min(collect_starts).strftime('%Y %b %d')
        collect_ends = [s[0] for s in arcpy.da.SearchCursor(merged, ['collect_end'])]
        collect_ends_max = max(collect_ends).strftime('%Y %b %d')
        first_by_area = [s[0] for s in arcpy.da.SearchCursor(merged, ['collect_start', 'area_field'], sql_clause = (None, 'ORDER BY area_field DESC'))][0]
        collect_majority = first_by_area.strftime('%Y %b %d')
    except:
    # else:
        collect_ends_max = 'Unknown'
        collect_starts_min = 'Unknown'
        collect_majority = 'Unknown'

####    if df.testForZero(tilesClip_local):
####        tiles_t_or_f = 'True'
####    else:
####        tiles_t_or_f = 'False'

    return collect_ends_max, collect_starts_min, collect_majority#, tiles_t_or_f


def getLasRasterArguments(lasToRaster):
    raster_arguments_list = [lasToRaster.getInput(t) for t in range(0, 7)]
    raster_arguments = 'LasToRaster arguments: ' + ', '.join(raster_arguments_list)

    return raster_arguments

def createCountsFromMultipoints(sgdb, maskRastOut, demListVal, demPtString, huc12, finalMPinm, finalMP, log, cntBeFile, named_cell_size, pattern):#locDict):
    try:

        cntBeFile_sized = updateResolution(cntBeFile, named_cell_size, demListVal, pattern, log)
        log.debug('---Counting Bare Earth Returns for ' + demListVal)
        terrCountName = arcpy.ValidateTableName("cnt" + demPtString + "m_fl_" + huc12, sgdb)
        terrCount = arcpy.PointToRaster_conversion(finalMPinm, arcpy.Describe(finalMP).OIDFieldName, os.path.join(sgdb, terrCountName), "COUNT", "NONE", demListVal)
        cntBeFileRasterObj = clipCountRaster(terrCount, maskRastOut, cntBeFile_sized)
# def createCountsFromMultipoints2(sgdb, maskRastBase, demList, huc12, finalMPinm, finalMP, log, cntBeFile, named_cell_size, pattern):#locDict):
#     try:
#         maskRastOut = opj(sgdb, maskRastBase + demList)

#         cntBeFile_sized = updateResolution(cntBeFile, named_cell_size, demList, pattern, log)
#         log.debug('---Counting Bare Earth Returns for ' + demList)
#         terrCountName = arcpy.ValidateTableName("cnt" + demList + "m_fl_" + huc12, sgdb)
#         terrCount = arcpy.PointToRaster_conversion(finalMPinm, arcpy.Describe(finalMP).OIDFieldName, os.path.join(sgdb, terrCountName), "COUNT", "NONE", demList)
#         cntBeFileRasterObj = clipCountRaster(terrCount, maskRastOut, cntBeFile_sized)

    except Exception as e:
        print('handling as exception')
##        log.debug(e.message)
        if sys.version_info.major == 2:
            print('handling as 2 exception')
            arcpy.AddError(e.message)
            print(e.message)
        elif sys.version_info.major == 3:
            print('handling as 3 exception')
            arcpy.AddError(e)
            print(e)

        # Get the traceback object
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        # Concatenate information together concerning the error into a message string
        pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
        # Return python error messages for use in script tool or Python Window
        arcpy.AddError(pymsg)
        # Print Python error messages for use in Python / Python Window
        print(pymsg + "\n")

        if arcpy.GetMessages(2) not in pymsg:
            msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
            arcpy.AddError(msgs)
            print(msgs)

    except:
        print('handling as except')
        # Get the traceback object
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]

        # Concatenate information together concerning the error into a message string
        pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
        # Return python error messages for use in script tool or Python Window
        arcpy.AddError(pymsg)
        # Print Python error messages for use in Python / Python Window
        print(pymsg + "\n")

        if arcpy.GetMessages(2) not in pymsg:
            msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
            arcpy.AddError(msgs)
            print(msgs)

    return cntBeFileRasterObj

def buildTerrainsUSGS(finalMP, FDSet, tcdFdSet, BREAKLINES, log, windows, ql):
    terrains = []
    # create one terrain with ZMinMax option
    if ql == 'QL 1' or ql == 'QL 0':
        spacing = "0.35"
    elif ql == 'QL 2':
        spacing = "0.7"
    else:
        spacing = "1.4"

    tp = [spacing, "", "", "WINDOWSIZE", "", "MILD", 0.18]
    pyrmd_str = "2 1000;4 2500;8 5000;16 10000;32 20000;64 40000"

    for window in windows:
        log.info(f'creating and setting up terrain: {window}')
        tp[4] = window

        LTrrn = arcpy.CreateTerrain_3d(FDSet, "Lidar_Trn" + "_" + tp[4], tp[0], tp[1], tp[2], tp[3], tp[4], tp[5], tp[6])

        pyrmd_str = "2 1000;4 2500;8 5000;16 10000;32 20000;64 40000"
        pyramids = arcpy.AddTerrainPyramidLevel_3d(LTrrn, "WINDOWSIZE", pyrmd_str)

        pyramids_arguments = 'WINDOWSIZE, ' + pyrmd_str

        # Value table rows: [feature_class, height_field, SF_type, group,
        #                     min_resolution, max_resolution, embed, embed_name]

        data_sources = [[finalMP, "SHAPE", "MASSPOINTS", 1, "", "", "", ""]]
    
        for key, info in BREAKLINES.items():
            data_sources.append([
                breakline_paths[key],
                info["height_field"],
                info["sf_type"],
                info["group"],
                "", "", "", "",
            ])
    
        arcpy.ddd.AddFeatureClassToTerrain(terrain_path, data_sources)
        print("Added mass points and breaklines to terrain")
    
# #        tf = setupTerrain(LTrrn, tcdFdSet, finalHb, finalHl, finalMP, finalNoZHb, poorZHb, log)#, badHb)
# # def setupTerrain(terrain, mask, breakpolys, breaklines, points, noZbreakpolys, poorZ, log):#, badHb):
#         group = 1#0
#         # list of terrain features
#         tf = []
#         ret = arcpy.AddFeatureClassToTerrain_3d(LTrrn, str(points) + " SHAPE masspoints " + str(group) + " 0 64 true true LASmerge_emb <None>")
#         tf.append(ret)
#         group += 1
#         if breakpolys:
#             ret = arcpy.AddFeatureClassToTerrain_3d(LTrrn, str(breakpolys) + " SHAPE hardreplace " + str(group) + " 0 32 true false <None> <None>")
#             tf.append(ret)
#             group += 1
#         if breaklines:
#             ret = arcpy.AddFeatureClassToTerrain_3d(LTrrn, str(breaklines) + " SHAPE hardline " + str(group) + " 0 32 true false <None> <None>")
#             tf.append(ret)
#             group += 1
#         if poorZ:
#             ret = arcpy.AddFeatureClassToTerrain_3d(LTrrn, str(poorZ) + " Z_Max hardline " + str(group) + " 0 32 true false <None> <None>")
#             tf.append(ret)
#             group += 1
#         if noZbreakpolys:
#             ret = arcpy.AddFeatureClassToTerrain_3d(LTrrn, str(noZbreakpolys) + " <None> hardline " + str(group) + " 0 32 true false <None> <None>")
#             tf.append(ret)
#             group += 1
#         ret = arcpy.AddFeatureClassToTerrain_3d(LTrrn, str(mask) + " <None> hardclip " + str(group) + " 0 32 true false <None> <None>")
#         tf.append(ret)
#         # return tf

        arcpy.BuildTerrain_3d(LTrrn)
        terrains.append(LTrrn)

        terrain_arguments_list = [LTrrn.getInput(t) for t in range(0,10)]
        terrain_arguments = ', '.join(terrain_arguments_list)

    return terrains, tf, terrain_arguments, pyramids_arguments


def buildTerrains(finalMP, FDSet, tcdFdSet, finalHb, finalHl, finalNoZHb, poorZHb, log, windows, time):
    terrains = []
    # create one terrain with ZMinMax option
    tp = ["1", "", "", "WINDOWSIZE", "ZMEAN", "MILD", 0.18]
    pyrmd_str = "2 1000;4 2500;8 5000;16 10000;32 20000;64 40000"

####    windowsizeMethods = windows#['ZMEAN', 'ZMINMAX']
    for window in windows:
        log.info(f'creating and setting up terrain: {window}')
        tp[4] = window

        LTrrn = arcpy.CreateTerrain_3d(FDSet, "Lidar_Trn" + "_" + tp[4], tp[0], tp[1], tp[2], tp[3], tp[4], tp[5], tp[6])

        pyrmd_str = "2 1000;4 2500;8 5000;16 10000;32 20000;64 40000"
        pyramids = arcpy.AddTerrainPyramidLevel_3d(LTrrn, "WINDOWSIZE", pyrmd_str)

        pyramids_arguments = 'WINDOWSIZE, ' + pyrmd_str

        tf = setupTerrain(LTrrn, tcdFdSet, finalHb, finalHl, finalMP, finalNoZHb, poorZHb, log)#, badHb)

        arcpy.BuildTerrain_3d(LTrrn)
        terrains.append(LTrrn)

        terrain_arguments_list = [LTrrn.getInput(t) for t in range(0,10)]
        terrain_arguments = ', '.join(terrain_arguments_list)

    return terrains, tf, terrain_arguments, pyramids_arguments

def terrain_args_from_inputs(terrain):
    terrain_arguments_list = [terrain.getInput(t) for t in range(0,10)]
    terrain_arguments = ', '.join(terrain_arguments_list)
    return terrain_arguments

def createCmDemRastersFromTerrains(log, demListVal, demPtString, maskRastOut, procDir, terrains, huc12, lidar_metadata_info, pyramid_args, dem_metadata_template, tElevFile_initial, named_cell_size, pattern22, interpDict, srOutNoVCS):
    try:
        # log.debug('snapRaster for Terrain to raster: ' + arcpy.env.snapRaster)
        interpTechnique = 'NATURAL_NEIGHBORS'
        pyramidLevel = '4'
        # dem_cellSize = demList[0]
        arcpy.env.cellSize = demListVal#_cellSize
        # windows are types of terrain (ZMEAN, ZMINMAX, etc.)
        for terrain in terrains:#window in windows:
            window = terrain.getInput(6)
            log.debug(f"processing window: {window}")

            tempTerrName = generateTempTerrName(procDir, window, demListVal, huc12)
            demOut = arcpy.TerrainToRaster_3d(terrain, tempTerrName, "FLOAT", interpTechnique, "CELLSIZE " + str(demListVal), pyramidLevel)

            ttr_arguments_list = [demOut.getInput(t) for t in range(0,5)]
            ttr_args = ', '.join(ttr_arguments_list)
            log.debug(f'TerrainToRaster args: {ttr_args}')

            nowYmd, collect_starts_min, collect_ends_max, collect_majority = [i for i in lidar_metadata_info]

            terrain_args = terrain_args_from_inputs(terrain)

            log.debug(f"terrain building args: {terrain_args}")

            paraDict = {
                '\n\nACPF: DEM Generation and Pit Fill Tool     ' : '\nRun Date: %s' % nowYmd,
                # '\nUnknown Vintage Lidar Data: ' : False,#tiles_t_or_f,
                '\nEarliest 3DEP Lidar Data: ' : collect_starts_min,
                '\nLatest 3DEP Lidar Data: ' : collect_ends_max,
                '\nLatest 3DEP Lidar Data: ' : collect_majority,
                '\nOutput conditioned DEM raster: ' : tempTerrName,#pfFileTemp,#fElevFile_interp,
                '\nTerrain Interpolation Arguments: ' : terrain_args,
                '\nTerrain Pyramid Arguments: ' : pyramid_args,
                '\nTerrain To Raster Arguments: ' : ttr_args
                }

            # ## update metadata
            # log.debug('---Adding metadata')
            # addMetadata(tempTerrName, paraDict, dem_metadata_template, log)

            tElevFile_internal = updateResolution(tElevFile_initial, named_cell_size, demListVal, pattern22, log)

            interpType = interpDict[window]
            # default interpolation type is mean18
            if interpType != 'mean18':
                tElevFile_interp = tElevFile_internal.replace('mean18', interpType)
            else:
                tElevFile_interp = tElevFile_internal

            log.debug('---Creating Masked Integer Centimeter DEM for ' + str(demListVal))
            tempDEM = Raster(tempTerrName)#rastr)
            intCmDEM = Int(tempDEM*100)
            maskedDEMint = ExtractByMask(intCmDEM, maskRastOut)

            # clear output VCS from here on out, creating centimeter rasters
            arcpy.env.outputCoordinateSystem = None
            arcpy.env.outputCoordinateSystem = srOutNoVCS

            log.debug('tElevFile name will be ' + tElevFile_interp)
            arcpy.BuildPyramids_management(maskedDEMint)
            # arcpy.CopyRaster_management(maskedDEMint, tElevFile_interp, pixel_type = '16_BIT_SIGNED', format = 'TIFF')
            maskedDEMint.save(tElevFile_interp)
            log.debug(f'Saved DEM for: {demPtString}')

            # f_dict = copy_md_summary_args(rastr)

            # ## update the metadata to reflect the output DEM just created
            # sep = ': '
            # key = "Output DEM raster from terrain"
            # # '\nOutput conditioned DEM raster: ' : pfFileTemp,tElevFile_interp,
            # f_dict.update({'\n' + key + sep: tElevFile_interp})

            ## update metadata
            # log.debug(f'---Skipping metadata for {tElevFile_interp}')
            log.debug(f'---Adding metadata to {tElevFile_interp}')
            addMetadata(tElevFile_interp, paraDict, dem_metadata_template, log)
            log.debug('---Added metadata')


    except Exception as e:
        print('handling as exception')
##        log.debug(e.message)
        if sys.version_info.major == 2:
            print('handling as 2 exception')
            arcpy.AddError(e.message)
            print(e.message)
        elif sys.version_info.major == 3:
            print('handling as 3 exception')
            arcpy.AddError(e)
            print(e)

        # Get the traceback object
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        # Concatenate information together concerning the error into a message string
        pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
        # Return python error messages for use in script tool or Python Window
        arcpy.AddError(pymsg)
        # Print Python error messages for use in Python / Python Window
        print(pymsg + "\n")

        if arcpy.GetMessages(2) not in pymsg:
            msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
            arcpy.AddError(msgs)
            print(msgs)

    except:
        print('handling as except')
        # Get the traceback object
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]

        # Concatenate information together concerning the error into a message string
        pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
        # Return python error messages for use in script tool or Python Window
        arcpy.AddError(pymsg)
        # Print Python error messages for use in Python / Python Window
        print(pymsg + "\n")

        if arcpy.GetMessages(2) not in pymsg:
            msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
            arcpy.AddError(msgs)
            print(msgs)



def setupLasDataset(lasIn, mask, procDir, breakpolys, breaklines, counterId, badHb, log, time, proj = ''):
    try:
        surfcons = ''
        if breakpolys:
            surfcons += str(breakpolys) + " Shape.Z Hard_Replace;"# " + str(mask) + " <None> Hard_Clip"
        if breaklines:
            surfcons += str(breaklines) + " Shape.Z Hard_Line;"# " + str(mask) + " <None> Hard_Clip"
        if badHb:
            surfcons += str(badHb) + " Z_Max Hard_Line;"# " + str(mask) + " <None> Hard_Clip"
        else:
            surfcons += str(mask) + " <None> Hard_Clip"
        log.info(f'surfcons are: {surfcons}')
##        log.debug('surfcons are ' + surfcons)
        lasd = arcpy.CreateLasDataset_management(lasIn, os.path.join(procDir, 'huc_cl2' + counterId + '.lasd'), in_surface_constraints = surfcons, spatial_reference = proj)

        return lasd
    except:
        # Get the traceback object
        #
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]

        # Concatenate information together concerning the error into a message string
        #
        pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
        msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"

        log.warning(pymsg)
        log.warning(msgs)


def setupTerrain(terrain, mask, breakpolys, breaklines, points, noZbreakpolys, poorZ, log):#, badHb):
    group = 1#0
    # list of terrain features
    tf = []
    ret = arcpy.AddFeatureClassToTerrain_3d(terrain, str(points) + " SHAPE masspoints " + str(group) + " 0 64 true true LASmerge_emb <None>")
    tf.append(ret)
    group += 1
    if breakpolys:
        ret = arcpy.AddFeatureClassToTerrain_3d(terrain, str(breakpolys) + " SHAPE hardreplace " + str(group) + " 0 32 true false <None> <None>")
        tf.append(ret)
        group += 1
    if breaklines:
        ret = arcpy.AddFeatureClassToTerrain_3d(terrain, str(breaklines) + " SHAPE hardline " + str(group) + " 0 32 true false <None> <None>")
        tf.append(ret)
        group += 1
    if poorZ:
        ret = arcpy.AddFeatureClassToTerrain_3d(terrain, str(poorZ) + " Z_Max hardline " + str(group) + " 0 32 true false <None> <None>")
        tf.append(ret)
        group += 1
    if noZbreakpolys:
        ret = arcpy.AddFeatureClassToTerrain_3d(terrain, str(noZbreakpolys) + " <None> hardline " + str(group) + " 0 32 true false <None> <None>")
        tf.append(ret)
        group += 1
    ret = arcpy.AddFeatureClassToTerrain_3d(terrain, str(mask) + " <None> hardclip " + str(group) + " 0 32 true false <None> <None>")
    tf.append(ret)
    return tf


def clipCountRaster(cnt3mFull, maskRastUTM, cntBeFileName):

    cntUTMpre = Con(IsNull(cnt3mFull), 0, cnt3mFull)
    cntUTMter = Con(maskRastUTM, cntUTMpre)
    cntUTMter.save(cntBeFileName)
    arcpy.BuildPyramids_management(cntUTMter)

    return cntUTMter



def testForZero(dataset):
    if type(dataset) == 'Raster':
        if dataset.hasRAT != True:
            arcpy.BuildRasterAttributeTable_management(dataset)

    try:
        fcount = int(arcpy.GetCount_management(dataset).getOutput(0))
    except:
        fcount = 0
    return fcount


def projIfNeeded(input2, output, srOut):
    srInput = arcpy.Describe(input2).spatialReference#maybe use projectionCode?
    try:
        if srInput.PCSCode != srOut.PCSCode:
    ##        log.debug('projecting ' + str(input2))
            projInput = arcpy.Project_management(input2, output, srOut)
        else:
            projInput = arcpy.CopyFeatures_management(input2, output)
    except arcpy.ExecuteError:
        ## assume there is a terrain created for the input, so must copy point feature class before projecting
        ptCopy = arcpy.CopyFeatures_management(input2, os.path.join(os.path.dirname(os.path.dirname(output)), os.path.basename(output) + '_copy'))
        projInput = arcpy.CopyFeatures_management(ptCopy, output)
    return projInput


def fill_donut_slow(fc):
    '''Edits a layer in-place and fills all donut holes or gaps in the selected
    features. Will operate on entire layer if there are no features selected.
    Requires layer to honor selected features.
    '''
    desc = arcpy.Describe(fc)
    shapefield = desc.ShapeFieldName
    rows = arcpy.UpdateCursor(fc)
    n = 0
    polyGeo = arcpy.Array() # to hold edited shape
    polyOuter = arcpy.Array() # to hold outer ring
    for row in rows:
        feat = row.getValue(shapefield)
        qInterior = False
        for partNum in range(feat.partCount) :
            part = feat.getPart(partNum)
            for pt in iter(lambda:part.next(),None) : # iter stops on null pt
                polyOuter.append(pt)
            if part.count > len(polyOuter) :
                qInterior = True
            polyGeo.append(polyOuter) # reassemble each part
            polyOuter.removeAll() # ready for next part
        if qInterior : # in any part of this feature, otherwise skip update
            n+=1
            row.setValue(shapefield,polyGeo)
            rows.updateRow(row)
        polyGeo.removeAll()
    del rows,row

def genClass2AndMultiPoints(allLAZ, sfx, srOut, inm, FDSet, procDir, log):
    try:
        if allLAZ.endswith('.laz') or allLAZ.endswith('.las'):
            """Filters LAS points to class 2 and creates multipoints in FDSet"""
            lasBase = os.path.splitext(os.path.basename(allLAZ))[0]
            cl2LAS = os.path.join(procDir, arcpy.ValidateTableName(lasBase + '_cl2', procDir) + '.las')
        ## Filter LAS to class 2 and 8 (key points), LAS 1.4 now has class 8 as reserved
            # LASTools requires " around file names with spaces, ' not allowed, (Windows Command Line too?)
            # Use pdal for this? https://pdal.io/en/stable/apps/translate.html#example-1
            # log.debug('--- Use las2las')
            if allLAZ.endswith('.laz'):
                log.debug('--- Using ConvertLas to decompress LAZ')
                las_from_laz = arcpy.ConvertLas_conversion(allLAZ, target_folder=procDir, compression=None, las_options=None)
            else:
                las_from_laz = allLAZ

            # rc = subprocess.call(os.path.join(softwareDir, 'LASTools', 'bin', 'las2las') + ' -i "' + allLAZ + '" -keep_class 2 8 -o ' + cl2LAS, shell=True)
            # if rc == 0 and not os.path.isfile(cl2LAS):
            #     log.warning('las2las did not create class 2 points file: ' + cl2LAS)
            #     lasMP = None
            # elif rc == 0:
            log.debug('--- Create las non-Minnesota Multipoint')
            lasMP = arcpy.LASToMultipoint_3d(cl2LAS, inm + "\\pts" + sfx, "1", class_code = [2,8], input_coordinate_system = srOut)
            # else:
                # log.warning('las2las did not execute successfully')

        ##        os.remove(cl2LAS)

        elif allLAZ.endswith('.zlas'):
            log.debug('--- Create zlas non-Minnesota Multipoint')
            lasMP = arcpy.LASToMultipoint_3d(allLAZ, inm + "\\pts" + sfx, "1", class_code = [2,8], input_coordinate_system = srOut)

        if lasMP:#allLAZ.endswith('.laz') and os.path.isfile(cl2LAS) or allLAZ.endswith('.zlas'):
        ## Clip multipoints
    ####        if untiledByLas:
    ####            lasMP = arcpy.Clip_analysis(lasMP, untiledByLas, inm + "\\pts_clp_" + str(rowCounter))
            ptsName = arcpy.ValidateTableName('pts_' + lasBase, os.path.join(str(FDSet)))
            ptOut = projIfNeeded(lasMP, os.path.join(str(FDSet), ptsName), srOut)

        else:
            log.warning('no ptOut created, setting to None')
            ptOut = None

    ##    return cl2LAS, ptOut
        return ptOut, cl2LAS

    except Exception as e:
        print('handling as exception')
##        log.debug(e.message)
        if sys.version_info.major == 2:
            arcpy.AddError(e.message)
            print(e.message)
        elif sys.version_info.major == 3:
            arcpy.AddError(e)
            print(e)

        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]

        # Concatenate information together concerning the error into a message string
        pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
        # Return python error messages for use in script tool or Python Window
        arcpy.AddError(pymsg)
        # Print Python error messages for use in Python / Python Window
        print(pymsg + "\n")

        if arcpy.GetMessages(2) not in pymsg:
            msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
            arcpy.AddError(msgs)
            print(msgs)

    except:
        print('handling as except')
        # Get the traceback object
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]

        # Concatenate information together concerning the error into a message string
        pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
        # Return python error messages for use in script tool or Python Window
        arcpy.AddError(pymsg)
        # Print Python error messages for use in Python / Python Window
        print(pymsg + "\n")

        if arcpy.GetMessages(2) not in pymsg:
            msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
            arcpy.AddError(msgs)
            print(msgs)




def buildSelection(inList, field):
    sel = ''
    for index, item in enumerate(inList):
        if index == 0:
            sel = field + ' = ' + str(item)
        else:
            sel += ' OR ' + field + ' = ' + str(item)
    return sel


def setupPointsAndBreaklines(finalMP, inm, FDSet, breakpolys, breaklines, log):
    try:
        fd_string = FDSet.getOutput(0)
##        with arcpy.EnvManager(workspace = str(FDSet)):
        log.info('setting up points and breaklines')
        # an in-memory version of the feature class for faster generation of count statistics
        finalMPinm = arcpy.CopyFeatures_management(finalMP, os.path.join(inm, 'mp_merge'))

        # a list of breakpoint feature classes to merge together at the end, saved for later
        breakpolyList = []

        # list of polygon breakline feature classes, some are 'better' than others, Minnesota you kill me!
        hbList = arcpy.ListFeatureClasses('hb_sel_*', feature_type = 'POLYGON', feature_dataset = os.path.basename(fd_string))
        fd_hbList = [opj(os.path.basename(fd_string), c) for c in hbList]
        hbListM = arcpy.ListFeatureClasses('hbm_sel_*', feature_type = 'POLYGON', feature_dataset = os.path.basename(fd_string))
        fd_hbListM = [opj(os.path.basename(fd_string), c) for c in hbListM]
        # put those with M values first (largest number of characters in text field)
####            hbListM += hbList
        fd_hbListM += fd_hbList
        if len(fd_hbListM) > 0:
            finalHb = arcpy.Merge_management(fd_hbListM, os.path.join(str(FDSet), 'hb_merge'))
            breakpolyList.append(finalHb)
        else:
            finalHb = None

        poorZHbList = arcpy.ListFeatureClasses('hb_poor_*', feature_type = 'POLYGON', feature_dataset = os.path.basename(fd_string))
        if len(poorZHbList) > 0:
            poorZHb = arcpy.Merge_management(poorZHbList, os.path.join(str(FDSet), 'poor_z_hb_merge'))
            breakpolyList.append(poorZHb)
        else:
            poorZHb = None

        noZHbList = arcpy.ListFeatureClasses('bad_hb_*', feature_type = 'POLYGON', feature_dataset = os.path.basename(fd_string))
        if len(noZHbList) > 0:
            finalNoZHb = arcpy.Merge_management(noZHbList, os.path.join(str(FDSet), 'no_z_hb_merge'))
            breakpolyList.append(finalNoZHb)
        else:
            finalNoZHb = None

        # merge the breakline feature classes together
        if breaklines is not None:
            df.create_needed_dirs_and_gdbs(breaklines, log)
        # breakGdb = os.path.dirname(breaklines)
        # if not arcpy.Exists(breakGdb):
        #     if not os.path.isdir(os.path.basename(breakGdb)):
        #         os.makedirs(os.path.dirname(breakGdb))
        #     breakGdbResult = arcpy.CreateFileGDB_management(os.path.dirname(breakGdb), os.path.basename(breakGdb))

        log.debug(f'breakpolyList: {breakpolyList}')
        if len(breakpolyList) > 0:
            # mergedBreakpolys = arcpy.Merge_management(breakpolyList, os.path.join(breakGdb, 'break_polys_' + huc12))
            mergedBreakpolys = arcpy.Merge_management(breakpolyList, breakpolys)

        hlList = arcpy.ListFeatureClasses('hl_*', feature_type = 'POLYLINE', feature_dataset = os.path.basename(fd_string))
        log.debug(f'hlList: {hlList}')
        if len(hlList) > 0:
            finalHl = arcpy.Merge_management(hlList, os.path.join(str(FDSet), 'hl_merge'))

            copiedBreaklines = arcpy.CopyFeatures_management(finalHl, breaklines)
            # copiedBreaklines = arcpy.CopyFeatures_management(finalHl, os.path.join(breakGdb, 'break_lines_' + huc12))
            log.debug(f'copiedBreaklines: {copiedBreaklines}')
        else:
            finalHl = None

        return finalMPinm, finalHb, poorZHb, finalNoZHb, finalHl

    except Exception as e:
        print('handling as exception')
##        log.debug(e.message)
        if sys.version_info.major == 2:
            arcpy.AddError(e.message)
            print(e.message)
        elif sys.version_info.major == 3:
            arcpy.AddError(e)
            print(e)

        errorhandle(sys.exc_info(), arcpy, traceback)#[2])

    except:
        errorhandle(sys.exc_info(), arcpy, traceback)#[2])


def setupPointsAndUSGSBreaklines(FDSet, breaks_ponds_lakes, breaks_streams_rivers, breaks_islands, breaks_bridges, BREAKLINES, log):
    try:
        fd_string = str(FDSet)#.getOutput(0)
##        with arcpy.EnvManager(workspace = str(FDSet)):
        log.info('setting up breaklines')

        # a list of breakpoint feature classes to merge together at the end, saved for later
        breakpolyList = []

        if len(breaks_ponds_lakes) > 0:
            final_ponds_lakes = arcpy.Merge_management(breaks_ponds_lakes, os.path.join(str(FDSet), 'breaks_ponds_lakes_merge'))
        else:
            final_ponds_lakes = None
        BREAKLINES['InlandPondsLakes']['path'] = final_ponds_lakes

        if len(breaks_streams_rivers) > 0:
            final_streams_rivers = arcpy.Merge_management(breaks_streams_rivers, os.path.join(str(FDSet), 'breaks_streams_rivers_merge'))
        else:
            final_streams_rivers = None
        BREAKLINES['InlandStreamsRivers']['path'] = final_streams_rivers

        if len(breaks_islands) > 0:
            final_islands = arcpy.Merge_management(breaks_islands, os.path.join(str(FDSet), 'breaks_islands_merge'))
        else:
            final_islands = None
        BREAKLINES['Islands']['path'] = final_islands

        if len(breaks_bridges) > 0:
            final_bridges = arcpy.Merge_management(breaks_bridges, os.path.join(str(FDSet), 'breaks_bridges_merge'))
        else:
            final_bridges = None
        BREAKLINES['Bridges']['path'] = final_bridges

        return final_ponds_lakes, final_streams_rivers, final_islands, final_bridges

    except Exception as e:
        print('handling as exception')
##        log.debug(e.message)
        if sys.version_info.major == 2:
            arcpy.AddError(e.message)
            print(e.message)
        elif sys.version_info.major == 3:
            arcpy.AddError(e)
            print(e)

        errorhandle(sys.exc_info(), arcpy, traceback)#[2])

    except:
        errorhandle(sys.exc_info(), arcpy, traceback)#[2])


def errorhandle(sei, arcpy, traceback):
    '''Try to handle the errors and output information about them'''
    tb = sys.exc_info()[2]
    tbinfo = traceback.format_tb(tb)[0]
    # Concatenate information together concerning the error into a message string
    pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
    # Return python error messages for use in script tool or Python Window
    arcpy.AddError(pymsg)
    # Print Python error messages for use in Python / Python Window
    print(pymsg + "\n")

    if arcpy.GetMessages(2) not in pymsg:
        msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
        arcpy.AddError(msgs)
        print(msgs)


# def generateLasArea(tilesClip, FDSet):
# ##    try:
#     # buffer/debuffer this by 2 meters to get rid of some gaps in lidar file characterization
#     tilesClipBuffer = arcpy.Buffer_analysis(tilesClip, buffer_distance_or_field = '2 METERS', dissolve_option = 'ALL')
#     tilesClipDeBuffer = arcpy.Buffer_analysis(tilesClipBuffer, buffer_distance_or_field = '-2 METERS')
#     tilesClipDslv = arcpy.Dissolve_management(tilesClipDeBuffer)
#     tilesClipDslvElim = arcpy.EliminatePolygonPart_management(tilesClipDslv, condition = 'PERCENT', part_area_percent = 50)
#     tcdFdSet = arcpy.CopyFeatures_management(tilesClipDslvElim, os.path.join(str(FDSet), 'local_las_area'))

#     return tcdFdSet

def updateResolution(filepath, init_res, new_res, pattern, log):
    """Take a filename with a specified resolution and alter it to the current processing resolution.
    This is done to reduce the number of arguments that are passed to the program."""
    # try:
    if init_res != new_res:
        if '.' in str(new_res):
            new_res = str(new_res).replace('.', 'pt')
        filename_path = Path(filepath)
        # check to see if it follows HUC DEM naming procedure
        st = filename_path.stem
        if re.search(pattern, st):#match(pattern, st):
            # replace within the name
            updated_filename = str(filename_path.name).replace(str(init_res) + 'm', str(new_res) + 'm')
            updated_filepath = str(filename_path.parent.joinpath(updated_filename))
        else:
            # append to the name
            updated_filename = str(filename_path.stem + "_" + str(new_res) + 'm' + filename_path.suffix)
            updated_filepath = str(filename_path.parent.joinpath(updated_filename))

        # updated_filename = filename.replace(str(init_res) + 'm' + huc12, str(new_res) + 'm' + huc12)
        log.debug(f"filename was: {filepath}; updated filename: {updated_filepath}")
    else:
        updated_filepath = filepath
    return updated_filepath


def buildLASRasters(lasdAll, beLayer, log, demListVal, demPtString, huc12, srSfx, maskRastOut, sgdb, procDir, int1rMaxFile, int1rMinFile, surfaceElevFile, intBeMaxFile, bareEarthReturnMinFile, cnt1rFile, cntPlsFile, cntBeFile, named_cell_size, internal_regions, lidar_metadata_info, derivative_metadata, pattern22):
##def buildLASRasters(lasdAll, lasdGround, log, demList, huc12, srSfx, maskRastBase, sgdb, procDir, int1rMaxFile, int1rMinFile, surfaceElevFile, frMinFile, intBeMaxFile, intBeMinFile, lastReturnMinFile, bareEarthReturnMinFile, cnt1rFile, named_cell_size, int_regions, ptr):
    '''creates multiple rasters from a las dataset, including min/max intensity of
    first return and bare earth surfaces, first return max and min surface, and z_range'''
    try:
        nowYmd, collect_starts_min, collect_ends_max, collect_majority = [i for i in lidar_metadata_info]

        paraDict = {
                '\n\nACPF: DEM Generation and Pit Fill Tool     ' : '\nRun Date: %s' % nowYmd,
                # '\nUnknown Vintage Lidar Data: ' : False,#tiles_t_or_f,
                '\nEarliest 3DEP Lidar Data: ' : collect_starts_min,
                '\nLatest 3DEP Lidar Data: ' : collect_ends_max,
                '\nLatest 3DEP Lidar Data: ' : collect_majority
                }

        if bareEarthReturnMinFile is not None or intBeMaxFile is not None:
            log.debug('---Creating LR Min layer')

            beReturnsMinTempFile = os.path.join(procDir, '_'.join(['tmp_bemin', demPtString + 'm', huc12, 'out.tif']))
            log.debug('---Creating LR Min raster')
            beReturnsMin = arcpy.LasDatasetToRaster_conversion(beLayer, beReturnsMinTempFile, interpolation_type = 'BINNING MINIMUM NONE', sampling_type = 'CELLSIZE', sampling_value = float(demListVal), data_type = 'FLOAT')
            log.debug('---Creating LR Min cm raster')
            beReturnsMinCm = Int(Times(beReturnsMin, 100))
            if bareEarthReturnMinFile is not None:
                bareEarthReturnMinFile_sized = updateResolution(bareEarthReturnMinFile, named_cell_size, demListVal, pattern22, log)
                log.debug('---Saving LR Min cm raster')
                beReturnsMinCm.save(bareEarthReturnMinFile_sized)#locDict['bareEarthReturnMinFile'])#.replace('fr', 'be'))
                log.debug('---Adding metadata to LR Min cm raster')
                addMetadata(bareEarthReturnMinFile_sized, paraDict, derivative_metadata, log)

        if int1rMaxFile is not None or int1rMinFile is not None or intBeMaxFile is not None:
            if demListVal == '2' or demListVal == '1':# should only run for 2m, otherwise too slow, but 1m needed for IA
                log.debug('---Creating intensity rasters')
                if int1rMaxFile is None and int1rMinFile is not None: 
                    log.warning('Faking int1rMaxFile value due to requested int1rMinFile')
                    int1rMaxFile_faked = int1rMinFile.replace('fr_int_min', 'fr_int_max')
                    int1rMaxFile = int1rMaxFile_faked
                elif int1rMaxFile is None and intBeMaxFile is not None:
                    log.warning('Faking int1rMaxFile value due to requested intBeMaxFile')
                    int1rMaxFile_faked = intBeMaxFile.replace('be_int_max', 'fr_int_max')
                    int1rMaxFile = int1rMaxFile_faked
                log.debug('---Creating FR Max Intensity')
                recode_tf = False
                log.debug(f'ir.max: {internal_regions.maximum},ir.min: {internal_regions.minimum}')
                int1rMaxFile_sized = updateResolution(int1rMaxFile, named_cell_size, demListVal, pattern22, log)
                if internal_regions.maximum - internal_regions.minimum != 0:
                    log.info('multiple regions')
                    # int1rMaxFile_sized_temp = opj(os.path.dirname(int1rMaxFile_sized), 'temp_' + os.path.basename(int1rMaxFile_sized))
                    int1rMaxFile_sized_temp = os.path.join(procDir, '_'.join(['tmp_frmax', demPtString + 'm', huc12, 'out.tif']))
                    lasd1rMaxIntensity = arcpy.LasDatasetToRaster_conversion(lasdAll, int1rMaxFile_sized_temp, 'INTENSITY', 'BINNING MAXIMUM NONE', sampling_type = 'CELLSIZE', sampling_value = float(demListVal), data_type = 'INT')
                    int_zs_max = ZonalStatistics(internal_regions, 'VALUE', int1rMaxFile_sized_temp, 'MAXIMUM')
                    if int_zs_max.minimum < 256 and int_zs_max.maximum > 256:
                        int_lt_256 = LessThan(int_zs_max, 256)
                        recode_areas = ZonalStatistics(internal_regions, 'VALUE', int_lt_256, 'MAXIMUM')
                        multiplied_intensities = Raster(int1rMaxFile_sized_temp) * 256
                        recoded_intensities = Con(recode_areas, multiplied_intensities, int1rMaxFile_sized_temp)
                        if int1rMaxFile is not None:
                            recoded_intensities.save(int1rMaxFile_sized)
                            arcpy.Delete_management(int1rMaxFile_sized_temp)
                        recode_tf = True
                    else:
                        if int1rMaxFile is not None:
                            log.info('all regions equal max intensity')
                            arcpy.CopyRaster_management(int1rMaxFile_sized_temp, int1rMaxFile_sized)
                            arcpy.Delete_management(int1rMaxFile_sized_temp)
                else:
                    log.info('one region')
                    if int1rMaxFile is not None:
                        lasd1rMaxIntensity = arcpy.LasDatasetToRaster_conversion(lasdAll, int1rMaxFile_sized, 'INTENSITY', 'BINNING MAXIMUM NONE', sampling_type = 'CELLSIZE', sampling_value = float(demListVal), data_type = 'INT')

                if int1rMaxFile is not None:
                    addMetadata(int1rMaxFile_sized, paraDict, derivative_metadata, log)

                if int1rMinFile is not None:
                    log.debug('---Creating FR Min Intensity')
                    int1rMinFile_sized = updateResolution(int1rMinFile, named_cell_size, demListVal, pattern22, log)
                    if recode_tf:
                        # int1rMinFile_sized_temp = opj(os.path.dirname(int1rMinFile_sized), 'temp_' + os.path.basename(int1rMaxFile_sized))
                        int1rMinFile_sized_temp = os.path.join(procDir, '_'.join(['tmp_frmin', demPtString + 'm', huc12, 'out.tif']))
                        lasd1rMinIntensity = arcpy.LasDatasetToRaster_conversion(lasdAll, int1rMinFile_sized_temp, 'INTENSITY', 'BINNING MINIMUM NONE', sampling_type = 'CELLSIZE', sampling_value = float(demListVal), data_type = 'INT')
                        multiplied_intensities = Raster(int1rMinFile_sized_temp) * 256
                        recoded_intensities = Con(recode_areas, multiplied_intensities, int1rMinFile_sized_temp)
                        recoded_intensities.save(int1rMinFile_sized)
                        arcpy.Delete_management(int1rMinFile_sized_temp)
                    else:
                        lasd1rMinIntensity = arcpy.LasDatasetToRaster_conversion(lasdAll, int1rMinFile_sized, 'INTENSITY', 'BINNING MINIMUM NONE', sampling_type = 'CELLSIZE', sampling_value = float(demListVal), data_type = 'INT')
                    addMetadata(int1rMinFile_sized, paraDict, derivative_metadata, log)

                if intBeMaxFile is not None:
                    log.debug('---Creating BE Max Intensity')
                    intBeMaxFile_sized = updateResolution(intBeMaxFile, named_cell_size, demListVal, pattern22, log)
                    if recode_tf:
                        # intBeMaxFile_sized_temp = opj(os.path.dirname(intBeMaxFile_sized), 'temp_' + os.path.basename(intBeMaxFile_sized))
                        intBeMaxFile_sized_temp = os.path.join(procDir, '_'.join(['tmp_bemax', demPtString + 'm', huc12, 'out.tif']))
                        lasdBeMaxIntensity = arcpy.LasDatasetToRaster_conversion(beLayer, intBeMaxFile_sized_temp, 'INTENSITY', 'BINNING MAXIMUM NONE', sampling_type = 'CELLSIZE', sampling_value = float(demListVal), data_type = 'INT')
                        multiplied_intensities = Raster(intBeMaxFile_sized_temp) * 256
                        recoded_intensities = Con(recode_areas, multiplied_intensities, intBeMaxFile_sized_temp)
                        recoded_intensities.save(intBeMaxFile_sized)
                    else:
                        lasdBeMaxIntensity = arcpy.LasDatasetToRaster_conversion(beLayer, intBeMaxFile_sized, 'INTENSITY', 'BINNING MAXIMUM NONE', sampling_type = 'CELLSIZE', sampling_value = float(demListVal), data_type = 'INT')
                    addMetadata(intBeMaxFile_sized, paraDict, derivative_metadata, log)

        if surfaceElevFile is not None:
            log.debug('---Creating FR Max surface')
            frMaxFile_sized = updateResolution(surfaceElevFile, named_cell_size, demListVal, pattern22, log)
            allReturnsMaxTempFile = os.path.join(procDir, '_'.join(['tmp_frmax', demPtString + 'm', huc12, 'out.tif']))
            # allReturnsMax = arcpy.LasDatasetToRaster_conversion(lasdAll, allReturnsMaxTempFile, interpolation_type = 'BINNING MAXIMUM SIMPLE', sampling_type = 'CELLSIZE', sampling_value = float(demListVal), data_type = 'FLOAT')
            allReturnsMax = arcpy.LasDatasetToRaster_conversion(lasdAll, allReturnsMaxTempFile, interpolation_type = 'BINNING MAXIMUM NONE', sampling_type = 'CELLSIZE', sampling_value = float(demListVal), data_type = 'FLOAT')
            allReturnsMaxCm = Int(Times(allReturnsMax, 100))
            allReturnsMaxCm.save(frMaxFile_sized)#locDict['surfaceElevFile'])#allReturnsMaxFile)
            addMetadata(frMaxFile_sized, paraDict, derivative_metadata, log)

        # log.debug('---Creating FR Min surface')
        # allReturnsMinTempFile = os.path.join(procDir, '_'.join(['tmp_frmin', demListVal + 'm', huc12, 'out.tif']))
        # allReturnsMin = arcpy.LasDatasetToRaster_conversion(lasdAll, allReturnsMinTempFile, interpolation_type = 'BINNING MINIMUM SIMPLE', sampling_type = 'CELLSIZE', sampling_value = float(demListVal), data_type = 'FLOAT')
        # allReturnsMinCm = Int(Times(allReturnsMin, 100))
        # allReturnsMinCm.save(frMinFile_sized)#locDict['firstReturnMinFile'])#allReturnsMinFile)

        if demListVal == '2' or demListVal == '1':# should only run for 2m, otherwise too slow, but 1m needed for IA
            if cnt1rFile is not None:
                log.debug('---Counting First Returns')
                cnt1rFile_sized = updateResolution(cnt1rFile, named_cell_size, demListVal, pattern22, log)
                cfrFileTemp = 'cnt_fr_' + demPtString + "m_" + huc12 + srSfx + '.tif'
                lasdCount = arcpy.LasPointStatsAsRaster_management(lasdAll, os.path.join(procDir, cfrFileTemp), 'POINT_COUNT', 'CELLSIZE', sampling_value = demListVal)
                cfrFileRasterObj = clipCountRaster(lasdCount, maskRastOut, cnt1rFile_sized)
                addMetadata(cnt1rFile_sized, paraDict, derivative_metadata, log)

            if cntPlsFile is not None:
                log.debug('---Counting All Returns')
                cntPlsFile_sized = updateResolution(cntPlsFile, named_cell_size, demListVal, pattern22, log)
                cntPlsFileTemp = 'cnt_pls_' + demPtString + "m_" + huc12 + srSfx + '.tif'
                lasdCount = arcpy.LasPointStatsAsRaster_management(lasdAll, os.path.join(procDir, cntPlsFileTemp), 'PULSE_COUNT', 'CELLSIZE', sampling_value = demListVal)
                cntPlsFileRasterObj = clipCountRaster(lasdCount, maskRastOut, cntPlsFile_sized)
                addMetadata(cntPlsFile_sized, paraDict, derivative_metadata, log)

            if cntBeFile is not None:
                log.debug('---Counting BE Returns')
                cntBeFile_sized = updateResolution(cntBeFile, named_cell_size, demListVal, pattern22, log)#.replace('_be_', '_belas_')
                cntBeFileTempSize = 'cnt_be_laspsr_' + demPtString + "m_" + huc12 + srSfx + '.tif'
                be_lasdCount = arcpy.LasPointStatsAsRaster_management(beLayer, os.path.join(procDir, cntBeFileTempSize), 'PULSE_COUNT', 'CELLSIZE', sampling_value = demListVal)
                # save the count raster with nulls converted to zeros, so that metadata can be added
                be_nulls = IsNull(be_lasdCount)#cntPlsFileRasterObj)
                be_nulls_as_zero = Con(be_nulls, 0, be_lasdCount)#cntPlsFileRasterObj)
                be_cntPlsFileRasterObj = clipCountRaster(be_lasdCount, maskRastOut, cntBeFile_sized)
##                        be_nulls_as_zero.save(cntBeFile_sized) 
                addMetadata(cntBeFile_sized, paraDict, derivative_metadata, log)

        # log.debug('---Counting Z Range')
        # zrangeFileTemp = 'zrng_all_' + demPtString + "m_" + huc12 + srSfx + '.tif'
        # lasdCount = arcpy.LasPointStatsAsRaster_management(lasdAll, os.path.join(procDir, zrangeFileTemp), 'Z_RANGE', 'CELLSIZE', sampling_value = demListVal)

##        lastLayer = arcpy.MakeLasDatasetLayer_management(lasdOut, 'ground_layer', [2,8], 'Last Return')
        # if sys.version_info.minor < 9:
        #     lastLayer = arcpy.MakeLasDatasetLayer_management(lasdGround, 'last_layer', return_values = 'Last Return')
        # else:
        #     lastLayer = arcpy.MakeLasDatasetLayer_management(lasdGround, 'last_layer', return_values = 'LAST')
        # log.debug('---Creating BE Max Intensity')
        # # minimum be seems to be a little less 'noisy' than maximum
        # lasdBeMaxIntensity = arcpy.LasDatasetToRaster_conversion(lastLayer, intBeMaxFile_sized, 'INTENSITY', 'BINNING MAXIMUM NONE', sampling_type = 'CELLSIZE', sampling_value = float(demListVal), data_type = 'INT')

        # log.debug('---Creating BE Min Intensity')
        # # minimum be seems to be a little less 'noisy' than maximum
        # lasdBeMinIntensity = arcpy.LasDatasetToRaster_conversion(lastLayer, intBeMinFile_sized, 'INTENSITY', 'BINNING MINIMUM NONE', sampling_type = 'CELLSIZE', sampling_value = float(demListVal), data_type = 'INT')

        # log.debug('---Creating LR Min surface')
        # lastReturnsMinTempFile = os.path.join(procDir, '_'.join(['tmp_lrmin', demPtString + 'm', huc12, 'out.tif']))
        # lastReturnsMin = arcpy.LasDatasetToRaster_conversion(lastLayer, lastReturnsMinTempFile, interpolation_type = 'BINNING MINIMUM SIMPLE', sampling_type = 'CELLSIZE', sampling_value = float(demListVal), data_type = 'FLOAT')
        # lastReturnsMinCm = Int(Times(lastReturnsMin, 100))
        # lastReturnsMinCm.save(lastReturnMinFile_sized)#locDict['lastReturnMinFile'])#.replace('fr', 'lr'))

##        # my current favorites  - minmax mild 18 terrain, binning min simple (almost too detailed),
##        interps = ['BINNING MINIMUM SIMPLE', 'BINNING AVERAGE SIMPLE', 'BINNING IDW SIMPLE', 'TRIANGULATION NATURAL_NEIGHBOR WINDOW_SIZE MINIMUM ' + demListVal, 'TRIANGULATION LINEAR WINDOW_SIZE MINIMUM ' + demListVal, 'TRIANGULATION NATURAL_NEIGHBOR WINDOW_SIZE CLOSEST_TO_MEAN ' + demListVal, 'TRIANGULATION LINEAR WINDOW_SIZE CLOSEST_TO_MEAN ' + demListVal]
##
##        log.debug('---Creating extra LR Min surfaces')
##        for interp in interps:
##            log.debug('interp is ' + str(interp))
######        interp = 'BINNING MINIMUM SIMPLE'
##            if interp[:3] == 'BIN':
##                selection = interp.split()[1][:3]
##                void = interp.split()[2][:3]
##                interpString = '_'.join([interp[:3], selection, void])
##            if interp[:3] == 'TRI':
##                selection = interp.split()[3][:3]
##                void = interp.split()[1][:3]
##                interpString = '_'.join([interp[:3], selection, void])
##            beReturnsMinTempFile = os.path.join(procDir, '_'.join(['tmp_be', interpString, demPtString + 'm', huc12, 'out.tif']))
##            beReturnsMin = arcpy.LasDatasetToRaster_conversion(beLayer, beReturnsMinTempFile, interpolation_type = interp, sampling_type = 'CELLSIZE', sampling_value = float(demListVal), data_type = 'FLOAT')
##
        arcpy.env.cellSize = None

    except Exception as e:
        print('handling as exception')
##        log.debug(e.message)
        if sys.version_info.major == 2:
            arcpy.AddError(e.message)
            print(e.message)
        elif sys.version_info.major == 3:
            arcpy.AddError(e)
            print(e)

        # Get the traceback object
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]

        # Concatenate information together concerning the error into a message string
        pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
        # Return python error messages for use in script tool or Python Window
        arcpy.AddError(pymsg)
        # Print Python error messages for use in Python / Python Window
        print(pymsg + "\n")

        if arcpy.GetMessages(2) not in pymsg:
            msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
            arcpy.AddError(msgs)
            print(msgs)

    except:
        print('handling as except')
        # Get the traceback object
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]

        # Concatenate information together concerning the error into a message string
        pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
        # Return python error messages for use in script tool or Python Window
        arcpy.AddError(pymsg)
        # Print Python error messages for use in Python / Python Window
        print(pymsg + "\n")

        if arcpy.GetMessages(2) not in pymsg:
            msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
            arcpy.AddError(msgs)
            print(msgs)


def generateTempTerrName(procDir, window, cellsize, huc12):
    if '.' not in str(cellsize):#not a decimal resolution raster
        tempTerrName = os.path.join(procDir, '_'.join(['tmp_ter', window, str(cellsize) + 'm', huc12, 'out.tif']))
    else:
        dec_res = str(cellsize)
        pt_res_name = dec_res.replace('.', 'pt')
        tempTerrName = os.path.join(procDir, '_'.join(['tmp_ter', window, pt_res_name + 'm', huc12, 'out.tif']))

    return tempTerrName


def dismantleTerrains(terrainList, finalHb, finalNoZHb, badHb, finalHl, tcdFdSet, log):
    '''Cleans up after building terrains and DEMs by deleting all input data'''

    for terrain in terrainList:
        if arcpy.Exists(terrain):
##        if 'LTrrn' in locals():
            if 'finalHb' in locals():
                if finalHb is not None:#terrain == True and finalHb is not None:#'finalHb' in locals():
                    arcpy.RemoveFeatureClassFromTerrain_3d(terrain, os.path.basename(str(finalHb)))
            if 'finalNoZHb' in locals():
                if finalNoZHb is not None:#terrain == True and finalHb is not None:#'finalHb' in locals():
                    arcpy.RemoveFeatureClassFromTerrain_3d(terrain, os.path.basename(str(finalNoZHb)))
            if 'badHb' in locals():
                if badHb is not None:#terrain == True and finalHb is not None:#'finalHb' in locals():
                    arcpy.RemoveFeatureClassFromTerrain_3d(terrain, os.path.basename(str(badHb)))
            if 'finalHl' in locals():
                if finalHl is not None:#terrain == True and finalHb is not None:#'finalHb' in locals():
                    arcpy.RemoveFeatureClassFromTerrain_3d(terrain, os.path.basename(str(finalHl)))
            if 'tcdFdSet' in locals():
                    arcpy.RemoveFeatureClassFromTerrain_3d(terrain, os.path.basename(str(tcdFdSet)))
##            arcpy.RemoveFeatureClassFromTerrain_3d(terrain, 'LASmerge_emb')
##            arcpy.Delete_management(terrain)

def copy_md_summary_args(md_item):
    """
    Copy an ACPF formatted metadata dictionary from one raster to another.
    Designed for use copying inputs from a pre-pitfilled DEM to a pit-filled one.
    md_item - dataset with metadata summary
    f_dict - dictionary with line breaks (as ACPF style) and arguments
    """    
    f_dict = {}
    sep = ': '
    f_metadata = md.Metadata(md_item)
    if f_metadata.summary is not None and f_metadata.summary != '':
        f_met_split = f_metadata.summary.split('\n')
        for i in f_met_split[2:]:
            key, value = i.split(sep, 1)
            # print(f'key {key} and value {value}')
            if i == f_met_split[2]: #prev row (f_met_split[1]) was \n
                f_dict.update({'\n\n' + key + sep: value})
            else:
                f_dict.update({'\n' + key + sep: value})
    return f_dict

def addMetadata(outDEM, paraDict, template_file_path, log=None):
    """ Claude AI contributed revision of addMetadata function to fix some bugs and improve functionality.
    Set the standard-format metadata XML file's path
    need to load metadata editor via 'import arcpy.metadata as md'
    outDEM = raster to receive updated metadata
    paraDict = dictionary of key/value pairs to be stored in metadata
      values stored include things like analyst, lidar acquisition date, etc.
    template_file_path = a template to load a basic summary from
    log = otional logging of error messages to a log file
    scriptPath = sys.path[0]
    """
    def _log(msg):
        """Safe logger that works whether log is provided or not."""
        if log:
            log.info(msg)

    try:
        # Get the target item's Metadata object
        tgt_item_md = md.Metadata(outDEM)

        if tgt_item_md is None:
            _log(f'Could not retrieve metadata object for {outDEM}')
            return

        if not tgt_item_md.isReadOnly:
            _log(f'Adding metadata to {outDEM}')

            # Import template metadata and save before further edits
            tgt_item_md.importMetadata(template_file_path)
            tgt_item_md.save()  # FIX: Save after import before modifying

            # Set title
            tgt_item_md.title = os.path.split(outDEM)[1]

            # FIX: Use getpass.getuser() instead of os.getlogin()
            tgt_item_md.credits = 'Analyst: %s' % getpass.getuser()

            # Build and set summary
            src_desc = tgt_item_md.summary or ''  # FIX: cleaner None handling
            for key, value in paraDict.items():
                src_desc += f'{key} {value}\n'   # Added newline for readability
            tgt_item_md.summary = src_desc

            # FIX: Removed unused tgt_item_md.copy() call before save
            tgt_item_md.save()
            _log(f'Metadata saved successfully for {outDEM}')

        else:
            _log(f'Metadata is read-only for {outDEM}, skipping.')

    except Exception as e:
        _log(f'Error adding metadata to {outDEM}: {e}')
        raise  # Re-raise so the caller knows something went wrong


def try_to_delete(rasRes, log):
    if arcpy.Exists(rasRes):
        try:
            arcpy.Delete_management(rasRes)
        except arcpy.ExecuteError:
            log.warning('could not remove using arcpy.Delete, trying os.remove')
            os.remove(rasRes)


BREAKLINES = {
    "Islands": {
        "reference_name": "Islands",
        "path": r"C:\gis\project\source.gdb\Island",
        "sf_type": "hardclip",
        "height_field": "SHAPE",
        "group": 1,
    },
    "Bridges": {
        "reference_name": "Bridges",
        "path": r"C:\gis\project\source.gdb\Bridge",
        "sf_type": "hardline",
        "height_field": "SHAPE",
        "group": 1,
    },
    "InlandPondsLakes": {
        "reference_name": "Inland_Ponds_Lakes",
        "path": r"C:\gis\project\source.gdb\InlandPondLake",
        "sf_type": "hardreplace",
        "height_field": "SHAPE",
        "group": 1,
    },
    "InlandStreamsRivers": {
        "reference_name": "Inland_Streams_Rivers",
        "path": r"C:\gis\project\source.gdb\InlandStreamRiver",
        "sf_type": "hardline",
        "height_field": "SHAPE",
        "group": 1,
    },
}

def doLidarDEMs(dem_boundary, wesm_huc12_tiles, laz_download_dir, 
        pdal_exe, gsds, procDir, snap, breakpolys, breaklines, 
        tElevFile, bareEarthReturnMinFile, firstReturnMaxFile, cntBeFile, cnt1rFile, cntPlsFile,
        int1rMinFile, int1rMaxFile, intBeMaxFile, cleanup, messages):
    arguments = [dem_boundary, wesm_huc12_tiles, laz_download_dir, 
        pdal_exe, gsds, procDir, snap, breakpolys, breaklines, 
        tElevFile, bareEarthReturnMinFile, firstReturnMaxFile, cntBeFile, cnt1rFile, cntPlsFile,
        int1rMinFile, int1rMaxFile, intBeMaxFile, cleanup]
                                                       
    for a in arguments:
        if a == arguments[0]:
            arg_str = str(a) + '\n'
        else:
            arg_str += str(a) + '\n'

    messages.addMessage("Tool: Executing with parameters:\n" + arg_str)

    arcpy.env.overwriteOutput = True

    arcpy.CheckOutExtension("Spatial")
    arcpy.CheckOutExtension("3D")

    arcpy.env.ZResolution = "0.01"

    try:
        huc12, huc8 = df.figureItOut(tElevFile)
        # the DEP huc DEM naming convention
        pattern27 = 'e[cfpxv][0-9]m\\d{10,16}'
        pattern26 = '[0-9]m\\d{10,16}'
        pattern25 = '[0-9]m_d{10,16}'
        pattern22 = '[+_][0-9]m[+_]'
        # huc12, huc8, named_cell_size = df.figureItOut(tElevFile)

        if procDir is not None:
            # if os.path.isdir(procDir):
            #     # log.info('nuking: ' + procDir)
            #     df.nukedir(procDir)

            if not os.path.isdir(procDir):
                os.makedirs(procDir)

            arcpy.env.scratchWorkspace = procDir
            sfldr = arcpy.env.scratchFolder
        else:
            sfldr = arcpy.env.scratchFolder
            procDir = sfldr

        # if lidar_download_directory is None:
        #     lidar_download_directory = procDir

        if not os.path.isdir(laz_download_dir):
            os.makedirs(laz_download_dir)

        sgdb = arcpy.env.scratchGDB
        arcpy.env.scratchWorkspace = sgdb
        arcpy.env.workspace = sgdb

        if procDir is None:
            procDir = sfldr

        # if lidar_download_directory is None:
        #     lidar_download_directory = sfldr

        #figure out where to create log files
        node = platform.node()
        logProc = df.defineLocalProc(node)
        if not os.path.isdir(logProc):
            logProc = sfldr

        if cleanup:
            log, nowYmd, logName, startTime = df.setupLoggingNoCh(logProc, sys.argv[0], huc12)
            arcpy.SetLogHistory = False
        else:
            # log to file and console
            log, nowYmd, logName, startTime = df.setupLoggingNew(logProc, sys.argv[0], huc12)
            arcpy.SetLogHistory = True

        # if not os.path.isfile(flib_metadata_template):
        #     log.warning('flib_metadata does not exist')
        # if not os.path.isfile(derivative_metadata):
        #     log.warning('derivative_metadata does not exist')
        log.info("Beginning execution: ")
        log.debug('sys.argv is: ' + str(sys.argv) + '\n')
        log.info("Processing HUC: " + huc12)
        log.info(f"procDir: {procDir}")
        # log.info(f"lidar_download_directory: {lidar_download_directory}")

        fElevDesc = arcpy.da.Describe(wesm_huc12_tiles)#dem_polygon)
        srOut = fElevDesc['spatialReference']
        srOutCode = srOut.PCSCode

        assert srOutCode < 32768, "EPSG spatial reference code too large, PDAL will not recognize"

        log.info("Output will be in EPSG Code (spatial reference): " + str(srOutCode))#sys.argv[9])
        log.info("Log file at " + logName)
        messages.addMessage("Log file at " + logName)

        flib_metadata_template, derivative_metadata = df.getMetadata(['flib', 'deriv'], procDir, log)

        ## store a list of all DEMs (lidar based, others) that must be joined to create HUC12
        ## Now a list of lists to facilitate creating two DEM resolutions easily (2 and 3 meter)
        rezes = gsds.split(",")
        log.info(f'Resolutions: {rezes}')
        ordered_rezes = []
        # for r in rezes:
        #     filename_path = Path(fElevFile)
        #     # check to see if it follows HUC DEM naming procedure
        #     stem = filename_path.stem
        #     pattern28 = '[0-9]m'
            # if re.search(pattern28, stem):
            #     log.debug(f"found match in {stem} of resolution {r}m")
            #     ordered_rezes.append(r)
            #     rezes.pop(r)

        for r in rezes:
            ordered_rezes.append(r)
        log.debug(f"ordered_rezes is {ordered_rezes}")

        # do lower to higher resolution
        # rezes.sort(reverse = True)
        demLists = [r for r in ordered_rezes]
        log.debug(f"demLists is {demLists}")
        # named_cell_size = demLists[0]
        init_res = demLists[0][0]
        log.debug(f"init_res is {init_res}")

        if init_res + 'm' not in tElevFile:
            tElevFile = os.path.splitext(tElevFile)[0] + '_' + str(init_res) + 'm' + os.path.splitext(tElevFile)[1]

        ## windowsizeMethods are the criterion used to select which point(s) in the window define the terrain
        interpDict = df.loadInterpDict()
        keys = interpDict.keys()
        windowsizeMethods = list(keys)[:2]# just get ['ZMEAN', 'ZMINMAX'] - extend to get 'ZMIN' and 'ZMAX'
        vals = list(interpDict.values())
        vals_lower = [v.lower() for v in vals]
        for v in vals_lower:
            # print(v)
            if v in tElevFile:
                interpFound = True
                interpType = v
                break
            else:
                interpFound = False
        else:
            if not interpFound:
                interpType = interpDict[windowsizeMethods[0]]
            else:
                log.warning('interpolation type not set')
        if not interpFound:# not in tElevFile:
            tElevFile = os.path.splitext(tElevFile)[0] + '_' + interpType + os.path.splitext(tElevFile)[1]
        log.debug(f'Revised tElevFile: {tElevFile}')

        # create output directories
        for filename in [tElevFile, cntBeFile, int1rMaxFile, firstReturnMaxFile]:
            if filename is not None:
                df.create_needed_dirs_and_gdbs(filename, log)

            # if not os.path.isdir(os.path.dirname(filename)):
            #     os.makedirs(os.path.dirname(filename))
            # create directories for alternate interpolation types/windowsize methods
                if interpType in filename:
                    for window in windowsizeMethods:
                        abbrev = interpDict[window]
                        if abbrev != interpType:
                            altfilename = filename.replace(interpType, abbrev)
                            df.create_needed_dirs_and_gdbs(altfilename, log)
                            # if not os.path.isdir(os.path.dirname(altfilename)):
                            #     os.makedirs(os.path.dirname(altfilename))

    ## If you set a scratch workspace first you can control where the scratchGDB or scratchFolder are created
    ## otherwise it defaults to a user's temp folder
    ## if you don't set anything it will go to 'in_memory'
        inm = 'in_memory'
        if snap is not None:#!= "":
            arcpy.env.snapRaster = snap

        # also set output to VCS 5703, NAVD88 Meters
        srOut = arcpy.SpatialReference(int(srOutCode), 5703)
        # since final output is in centimeters, create one without VCS
        srOutNoVCS = arcpy.SpatialReference(int(srOutCode))
        arcpy.env.outputCoordinateSystem = srOut
        srSfx = '_'+str(srOutCode)
        work_id_name = 'workunit_id'

        # maskFc, maskFc_area, maskFcOut, maskRastOut, hucRastOut, FDSet = prepPolygonBoundary(dem_polygon, log, sgdb, srOut, srSfx, maskRastBase, demLists)

        # assert int(arcpy.GetCount_management(dem_polygon).getOutput(0)) < 2, 'multiple features in polygon feature class'
        # assert int(arcpy.GetCount_management(dem_polygon).getOutput(0)) > 0, 'no features in polygon feature class'
        # tiles_dissolve = arcpy.Dissolve_management(wesm_huc12_tiles, 'in_memory/dem_poly_dissolve')
        maskFc = arcpy.CopyFeatures_management(dem_boundary)
        merged_area = [s[0] for s in arcpy.da.SearchCursor(maskFc, ['SHAPE@AREA'])][0]

        if 'id' not in df.getfields(maskFc):
            arcpy.AddField_management(maskFc, 'id', 'LONG')
            arcpy.CalculateField_management(maskFc, 'id', 1, 'PYTHON')
        log.info("maskFcPrelim complete")

        ## Set up geodatabase to store the multipoint files and terrains (necessary all inputs be in feature dataset
        # Vertical units are in meters (float) so use a meter-based reference
        FDSet = arcpy.CreateFeatureDataset_management(sgdb, "Lidar_pts", srOut)
        maskFcOut = projIfNeeded(maskFc, os.path.join(str(FDSet), 'buf_huc' + srSfx), srOut)
        log.info("maskFcOut complete at: ")

        maskRastBase = 'mask_rast_'
        for demListVal in demLists:
            if '.' in demListVal:
                demPtString = demListVal.replace('.', 'pt')
            else:
                demPtString = demListVal
            maskRastOut = arcpy.PolygonToRaster_conversion(maskFcOut, 'id', opj(sgdb, maskRastBase + demPtString), cellsize = float(demListVal))
            # huc_rast_out = arcpy.conversion.PolygonToRaster(geom_copy, 'OBJECTID', opj(sgdb, 'huc_rast' + str(demList[0])), cellsize = demList[0])

##----------------------------------------------------------------------

    # # check for collection change (different priorities) to restrict further data
        fixedFolder = str(arcpy.CreateFolder_management(sfldr, 'fixed_crs'))


##----------------------------------------------------------------------


##----------------------------------------------------------------------
        build_threshold = 0.99#99 - reduced to 0.99 for Lake Erie at Toledo (04100009)
        if df.testForZero(wesm_huc12_tiles):

            arcpy.env.outputCoordinateSystem = srOut

            allTilesList = []

            # with arcpy.da.SearchCursor(wesm_huc12_tiles, ['SHAPE@', 'laz_file', work_id_name]) as scur:#, sql_clause = [None, 'ORDER BY order DESC']) as scur:
            with arcpy.da.SearchCursor(wesm_huc12_tiles, ['SHAPE@', 'laz_file', work_id_name, 'ql', 'dem_gsd_meters', 'horiz_crs', 'vert_crs']) as scur:
                for row_counter, srow in enumerate(scur):
                    work_id = srow[2]
                    storage_laz = srow[1]
                    geom = srow[0]
                    ql = srow[3]
                    vert_crs = srow[-1]
                    hor_crs = srow[-2]
                    # if 'DEP\\laz' in storage_laz:
                    #     log.info('found DEP\\laz in storage_laz path, replacing with DEP\\USGS_LPC')
                    #     if 'DEP\\laz' in storage_laz:
                    #         log.info('found DEP\\laz in storage_laz path, replacing with DEP\\USGS_LPC')
                    #         ept_las = storage_laz.replace('DEP\\laz', 'DEP\\USGS_LPC')
                    # elif 'E:\\DEP_Checkout' in storage_laz:
                    #     log.info('found E:\\DEP_Checkout in storage_laz path, replacing with M:\\DEP')
                    #     ept_las = storage_laz.replace('E:\\DEP_Checkout', 'M:\\DEP')
                    # else:
                    #     ept_las = storage_laz
                    if os.path.exists(storage_laz):#ept_las):

                        storage_laz_altsep = storage_laz.replace(os.path.sep, os.path.altsep)
                        if arcpy.Exists(storage_laz_altsep.replace('M:/DEP/USGS_LPC/IL_10CountyNRCS_D23/IL_10CoNRCS_1_D23/LAZ', 'E:/DEP_Proc/DEMProc/LAS_dem2013_1m_071300060105/laz_for_huc12')):
                            storage_laz_altsep = storage_laz_altsep.replace('M:/DEP/USGS_LPC/IL_10CountyNRCS_D23/IL_10CoNRCS_1_D23/LAZ', 'E:/DEP_Proc/DEMProc/LAS_dem2013_1m_071300060105/laz_for_huc12')
                        fixedLasPath = opj(fixedFolder, os.path.basename(storage_laz).replace('.laz','.las'))
                        fixedLasPath_altsep = fixedLasPath.replace(os.path.sep, os.path.altsep)##
                        if ql == 'QL 0':
                            pipeline_laz_las = '''{
"pipeline": [
    {
    "type": "readers.las",
    "filename": "''' + storage_laz_altsep + '''",
    "tag": "source1"
    },
    {
      "type": "filters.decimation",
      "step": 4,
      "offset": 0,
      "tag": "decimated"
    },
    {
    "type": "filters.reprojection",
    "in_srs": "''' + 'EPSG:' + hor_crs + '''+''' + vert_crs + '''",
    "out_srs": "''' + 'EPSG:' + str(srOutCode) + '''+5703",
    "tag": "reprojected_local_NAVD88m"
    },
    {
    "type": "writers.las",
    "filename": "''' + fixedLasPath_altsep + '''",
    "tag": "writerslas"
    }
]}'''
                        elif ql == 'QL 1':
                            pipeline_laz_las = '''{
"pipeline": [
    {
    "type": "readers.las",
    "filename": "''' + storage_laz_altsep + '''",
    "tag": "source1"
    },
    {
      "type": "filters.decimation",
      "step": 4,
      "offset": 0,
      "tag": "decimated"
    },
    {
    "type": "filters.reprojection",
    "in_srs": "''' + 'EPSG:' + hor_crs + '''+''' + vert_crs + '''",
    "out_srs": "''' + 'EPSG:' + str(srOutCode) + '''+5703",
    "tag": "reprojected_local_NAVD88m"
    },
    {
    "type": "writers.las",
    "filename": "''' + fixedLasPath_altsep + '''",
    "tag": "writerslas"
    }
]}'''
                        else:
                            pipeline_laz_las = '''{
"pipeline": [
    {
    "type": "readers.las",
    "filename": "''' + storage_laz_altsep + '''",
    "tag": "source1"
    },
    {
    "type": "filters.reprojection",
    "in_srs": "''' + 'EPSG:' + hor_crs + '''+''' + vert_crs + '''",
    "out_srs": "''' + 'EPSG:' + str(srOutCode) + '''+5703",
    "tag": "reprojected_local_NAVD88m"
    },
    {
    "type": "writers.las",
    "filename": "''' + fixedLasPath_altsep + '''",
    "tag": "writerslas"
    }
]}'''
                        if ql == 'QL 0' or ql == 'QL 1':
                            log.debug(f'using decimating pipeline_laz_las')
                        pl_laz_las = pdal.Pipeline(pipeline_laz_las)
                        ex = pl_laz_las.execute()

                        ept_las_base = os.path.splitext(os.path.basename(storage_laz))[0]


            # # with arcpy.da.SearchCursor(wesm_huc12, ['SHAPE@', work_id_name] + url_list, sql_clause = [None, 'ORDER BY ' + addOrderField.getInput(1) + ' DESC']) as scur:#work_id_name, 'SHAPE@AREA', 'lpc_link']) as scur:
            #     for srow in scur:
            #         work_id = srow[2]
            #         storage_laz = srow[1]
            #         geom = srow[0]

            #         # print(srow)
            #         # ept_las_full_filename = laz
            #         # handle some inconsistent paths
            #         if 'DEP\\laz' in storage_laz:
            #             log.info('found DEP\\laz in storage_laz path, replacing with DEP\\USGS_LPC')
            #             if 'DEP\\laz' in storage_laz:
            #                 log.info('found DEP\\laz in storage_laz path, replacing with DEP\\USGS_LPC')
            #                 ept_las = storage_laz.replace('DEP\\laz', 'DEP\\USGS_LPC')
            #         elif 'E:\\DEP_Checkout' in storage_laz:
            #             log.info('found E:\\DEP_Checkout in storage_laz path, replacing with M:\\DEP')
            #             ept_las = storage_laz.replace('E:\\DEP_Checkout', 'M:\\DEP')
            #         else:
            #             ept_las = storage_laz
            #         if os.path.exists(ept_las):#_full_filename):# and stats.st_size > las_size_threshold:
            #         #     cl2Las = processEptLas(sgdb, sfldr, srOutCode, fixedFolder, geom_srOut, ept_las_full_filename, srOut, inm, FDSet, procDir, allTilesList, log, time, work_id)

            #         # '''process a cursor row of data by creating a suitable las file from the input las/laz/zlas dataset
            #         # This inlcudes project and clipping las data files into output dataset and also creating a multipoint
            #         # file from the las data if there is any within the extent'''
            #         # try:
            #             ept_las_base = os.path.splitext(os.path.basename(ept_las))[0]
                        sfx = arcpy.ValidateTableName('_' + ept_las_base, sgdb)
                        # log.debug('lidar file suffix is: ' + sfx)

                        # allLasd = arcpy.CreateLasDataset_management(ept_las, opj(sfldr, 'all' + sfx))
                        log.debug(f'lidar file suffix is: {sfx}')
                        lasdAll = arcpy.CreateLasDataset_management(fixedFolder, os.path.join(procDir, 'huc_all.lasd'), spatial_reference=arcpy.SpatialReference(int(srOutCode)))

                    #     # extract to tile geometry and project if necessary
                    #     nameSfx = '_' + str(srOutCode)
                    #     fixedLasBasename = os.path.basename(ept_las)[:-4] + nameSfx + '.las'
                    # ##            log.debug('fixedLasBasename: ' + fixedLasBasename)
                    #     # some old 3DEP projects don't alway have boundaries and data lining up...
                    #     log.debug(f'ExtractLas arguments are {allLasd}, {fixedFolder}, name_suffix = {nameSfx}, rearrange_points = {"MAINTAIN_POINTS"}, out_las_dataset = {opj(fixedFolder, "fixed" + sfx + ".lasd")}')
                    #     # if work_id < 0:
                    #     #     tileGeomBuffer5 = geom.buffer(5) #tile Geometry column
                    #     #     log.debug('--- Use ExtractLas, boundary option,')
                    #     #     fixedLasd = arcpy.ExtractLas_3d(allLasd, fixedFolder, name_suffix = nameSfx, rearrange_points = "MAINTAIN_POINTS", out_las_dataset = opj(fixedFolder, 'fixed' + sfx + '.lasd'), compression="NO_COMPRESSION", boundary = tileGeomBuffer5)
                    #     # else:
                    #     #     log.debug('--- Use ExtractLas, no boundary option,')
                    #     #     fixedLasd = arcpy.ExtractLas_3d(allLasd, fixedFolder, name_suffix = nameSfx, rearrange_points = "MAINTAIN_POINTS", out_las_dataset = opj(fixedFolder, 'fixed' + sfx + '.lasd'), compression="NO_COMPRESSION")#, boundary = tileGeomBuffer5)
                    #     # log.debug(fixedLasd.getMessages())
                    #     # fixedLasdDescDa = arcpy.da.Describe(fixedLasd)
                    #     # fixedLasPath = opj(fixedLasdDescDa['path'], fixedLasBasename)

                    #     # log.debug('--- Done creating LAS dataset and extracting LAS at ')

                    #     # if fixedLasdDescDa['pointCount'] > 0:
                    #     #     allTilesList.append(fixedLasPath)

                        if os.path.exists(fixedLasPath_altsep):
                            # 'Filters LAS points to class 2 and creates multipoints in FDSet'
                            # lasBase = os.path.splitext(os.path.basename(allLAZ))[0]
                            log.debug('--- Create las non-Minnesota Multipoint')
                            if ql == 'QL 0' or ql == 'QL 1':
                                spacing = '0.125' # meters (UTM)
                            else:
                                spacing = '1'
                            lasMP = arcpy.LASToMultipoint_3d(fixedLasPath_altsep, inm + '\\pts' + sfx, spacing, class_code=[2, 8], input_coordinate_system=srOut)

                        # if fixedLasPath in allTilesList:#non 0 amount of lidar points in las
                        #     allLAZ = fixedLasPath
                        #     if allLAZ.endswith('.laz') or allLAZ.endswith('.las'):
                        #         """Filters LAS points to class 2 and creates multipoints in FDSet"""
                        #         lasBase = os.path.splitext(os.path.basename(allLAZ))[0]
                        #         if allLAZ.endswith('.laz'):
                        #             log.debug('--- Using ConvertLas to decompress LAZ')
                        #             las_from_laz = arcpy.ConvertLas_conversion(allLAZ, target_folder=procDir, compression=None, las_options=None)
                        #         else:
                        #             las_from_laz = allLAZ

                        #         log.debug('--- Create las non-Minnesota Multipoint')
                        #         lasMP = arcpy.LASToMultipoint_3d(las_from_laz, inm + "\\pts" + sfx, "1", class_code = [2,8], input_coordinate_system = srOut)

                        #     elif allLAZ.endswith('.zlas'):
                        #         log.debug('--- Create zlas non-Minnesota Multipoint')
                        #         lasMP = arcpy.LASToMultipoint_3d(allLAZ, inm + "\\pts" + sfx, "1", class_code = [2,8], input_coordinate_system = srOut)
                        #         las_from_laz = allLAZ

                            if lasMP:
                                ptsName = arcpy.ValidateTableName('pts' + sfx + '_' + str(srOutCode), os.path.join(str(FDSet)))
                                ptOut = projIfNeeded(lasMP, os.path.join(str(FDSet), ptsName), srOut)
                                arcpy.Delete_management(lasMP)
                            else:
                                log.warning('no ptOut created, setting to None')
                                ptOut = None
                        #     cl2Las = las_from_laz
                        # if 'cl2Las' not in locals():
                        #     cl2Las = None
                        if 'ptOut' not in locals():
                            ptOut = None
                        log.info('ptOut: ' + str(ptOut))# + ' and cl2Las ' + str(cl2Las))

                        # now find the USGS standard breaklines and breakpolygons for this work_id
                        if row_counter == 0:# set things up so the below code can be run for the first row and then only when work_id changes
                            prev_work_id = work_id

                        if row_counter == 0 or work_id != prev_work_id:
                            breaks_md_dir = opj(os.path.dirname(os.path.dirname(storage_laz)), 'breaks_md')
                            arcpy.env.workspace = breaks_md_dir
                            breaks_gdb = arcpy.ListWorkspaces()[0]

                            # populate breaklines and breakpolys with the first row's project
                            for k in BREAKLINES.keys():
                                ref = BREAKLINES[k]['reference_name']
                                candidate = arcpy.ListFeatureClasses(ref + '*')[0]
                                if not candidate:
                                    possibilities = arcpy.ListFeatureClasses()
                                    match = difflib.get_close_matches(ref, possibilities, n=1, cutoff=0.6)
                                    if match:
                                        log.warning(f"No exact match for breaks found. Closest match found: '{match}'")
                                        candidate = match
                                    else:
                                        log.info("No close match found above the threshold.")
                                        candidate = None
                                print(candidate)
                                if candidate and row_counter == 0:
                                    breakline_fc = opj(breaks_gdb, candidate)
                                    if ref == 'Islands':
                                        breaks_islands = [breakline_fc]
                                    elif ref == 'Bridges':
                                        breaks_bridges = [breakline_fc]
                                    elif ref == 'Inland_Ponds_Lakes':
                                        breaks_ponds_lakes = [breakline_fc]
                                    elif ref == 'Inland_Streams_Rivers':
                                        breaks_streams_rivers = [breakline_fc]

                                elif candidate:
                                    breakline_fc = opj(breaks_gdb, candidate)
                                    if ref == 'Islands':
                                        breaks_islands = [breakline_fc]
                                    elif ref == 'Bridges':
                                        breaks_bridges = [breakline_fc]
                                    elif ref == 'Inland_Ponds_Lakes':
                                        breaks_ponds_lakes = [breakline_fc]
                                    elif ref == 'Inland_Streams_Rivers':
                                        breaks_streams_rivers = [breakline_fc]

                                else:
                                    log.warning(f"No matches for breaklines found for {breaks_gdb}")

                            prev_work_id = work_id

            wesm_huc12_tiles_buffer = arcpy.analysis.Buffer(wesm_huc12_tiles, opj(sgdb, 'wesm_huc12_tiles_buffer'), '0.1 Meters')#, dissolve_option = 'ALL')    
            wesm_huc12_tiles_buffer_dissolve = arcpy.management.Dissolve(wesm_huc12_tiles_buffer, opj(sgdb, 'wesm_huc12_tiles_dissolve'), work_id_name)

            ptr_poly = arcpy.management.CopyFeatures(wesm_huc12_tiles_buffer_dissolve, opj(sgdb, 'ptr_poly'))
            arcpy.AddField_management(ptr_poly, 'area_field', 'DOUBLE')
            arcpy.CalculateField_management(ptr_poly, 'area_field', '!shape.area!', 'Python 3')

            # find internal regions where work_id_name differs for checking for different intensity maximums
            ptr = arcpy.conversion.PolygonToRaster(ptr_poly, work_id_name, opj(sgdb, 'ptr'))
            int_ptr = Int(ptr)
            fs1 = FocalStatistics(int_ptr, NbrRectangle(5,5), 'RANGE')
            internal_regions = Con(fs1 == 0, int_ptr)
            internal_regions.save(opj(sgdb, 'intrnl_rgns'))
            
            ##----------------------------------------------------------------------
            # now build datasets
            arcpy.env.workspace = sgdb
            mpList = arcpy.ListFeatureClasses('pts_*', feature_type = 'POINT', feature_dataset = os.path.basename(FDSet.getOutput(0)))
            if len(mpList) > 0:
                finalMP = arcpy.Merge_management(mpList, os.path.join(str(FDSet), 'mp_merge'))
                if df.testForZero(finalMP):
                    final_breaks_ponds_lakes, final_breaks_streams_rivers, final_breaks_islands, final_breaks_bridges = setupPointsAndUSGSBreaklines(FDSet, breaks_ponds_lakes, breaks_streams_rivers, breaks_islands, breaks_bridges, BREAKLINES, log)

                    tcdFdSet = arcpy.management.Dissolve(wesm_huc12_tiles_buffer_dissolve, os.path.join(str(FDSet), 'ept_and_local_las'))
                    fill_donut_slow(tcdFdSet)

                    terrains, tf, terrain_args, pyramid_args = buildTerrainsUSGS(finalMP, FDSet, tcdFdSet, BREAKLINES, log, windowsizeMethods, ql)

                    if sys.version_info.minor < 9:
                        beLayer = arcpy.MakeLasDatasetLayer_management(lasdAll, 'ground_layer', [2,8], 'Last Return')
                    else:
                        beLayer = arcpy.MakeLasDatasetLayer_management(lasdAll, 'ground_layer', [2,8], 'LAST')

                    collect_ends_max, collect_starts_min, collect_majority = getLidarTimeframes(ptr_poly)#prev_merged)#merged_copy)#, tilesClip_local)

                    lidar_metadata_info = [nowYmd, collect_starts_min, collect_ends_max, collect_majority]

                    # lasdAll = arcpy.CreateLasDataset_management(allTilesList, os.path.join(procDir, 'huc_all.lasd'), spatial_reference = arcpy.SpatialReference(int(srOutCode)))
                    # ## Following code runs slowly at times and is not being used further 2023.12.21
                    # # classify overlap in lasdAll

        arcpy.env.cellSize = None

##----------------------------------------------------------------------
        dem_boundary_area = [s[0] for s in arcpy.da.SearchCursor(dem_boundary, ['SHAPE@AREA'])][0]

        if merged_area / dem_boundary_area >= build_threshold:
            for demListVal in demLists:
                log.debug(f'---Processing resolution {demListVal}')
                if '.' in demListVal:
                    demPtString = demListVal.replace('.', 'pt')
                else:
                    demPtString = demListVal
##                    maskRastOut = arcpy.PolygonToRaster_conversion(maskFcOut, 'id', opj(sgdb, maskRastBase + demPtString), cellsize = float(demListVal))
                maskRastOutName = opj(sgdb, maskRastBase + demPtString)#demListVal)
                # if cntBeFile is not None:
                #     if int(demListVal) >= 1: #creating point counts for high resolution is very slow (slower than from LAS Datasets)
                #         maskRastOutName = opj(sgdb, maskRastBase + demPtString)#demListVal)
                #         cntBeFileRasterObj = createCountsFromMultipoints(sgdb, maskRastOutName, demListVal, demPtString, huc12, finalMPinm, finalMP, log, cntBeFile, init_res, pattern22)

                # if not arcpy.Exists(tElevFile):
                terrainList = createCmDemRastersFromTerrains(log, demListVal, demPtString, maskRastOutName, procDir, terrains, huc12, lidar_metadata_info, pyramid_args, flib_metadata_template, tElevFile, init_res, pattern22, interpDict, srOutNoVCS)

                buildLASRasters(lasdAll, beLayer, log, demListVal, demPtString, huc12, srSfx, maskRastOutName, sgdb, procDir, int1rMaxFile, int1rMinFile, firstReturnMaxFile, intBeMaxFile, bareEarthReturnMinFile, cnt1rFile, cntPlsFile, cntBeFile, init_res, internal_regions, lidar_metadata_info, derivative_metadata, pattern22)
#                 lasdAll, beLayer, log, demListVal, demPtString, huc12, srSfx, maskRastOut, sgdb, procDir, int1rMaxFile, int1rMinFile, surfaceElevFile, intBeMaxFile, bareEarthReturnMinFile, cnt1rFile, cntPlsFile, cntBeFile, named_cell_size, internal_regions, lidar_metadata_info, derivative_metadata, pattern22 = lasdAll, beLayer, log, demListVal, demPtString, huc12, srSfx, maskRastOutName, sgdb, procDir, int1rMaxFile, int1rMinFile, firstReturnMaxFile, intBeMaxFile, bareEarthReturnMinFile, cnt1rFile, cntPlsFile, cntBeFile, init_res, internal_regions, lidar_metadata_info, derivative_metadata, pattern22
#                 nowYmd, collect_starts_min, collect_ends_max, collect_majority = [i for i in lidar_metadata_info]

#                 paraDict = {
#                         '\n\nACPF: DEM Generation and Pit Fill Tool     ' : '\nRun Date: %s' % nowYmd,
#                         # '\nUnknown Vintage Lidar Data: ' : False,#tiles_t_or_f,
#                         '\nEarliest 3DEP Lidar Data: ' : collect_starts_min,
#                         '\nLatest 3DEP Lidar Data: ' : collect_ends_max,
#                         '\nLatest 3DEP Lidar Data: ' : collect_majority
#                         }

#                 if bareEarthReturnMinFile is not None or intBeMaxFile is not None:
#                     log.debug('---Creating LR Min layer')

#                     beReturnsMinTempFile = os.path.join(procDir, '_'.join(['tmp_bemin', demPtString + 'm', huc12, 'out.tif']))
#                     log.debug('---Creating LR Min raster')
#                     beReturnsMin = arcpy.LasDatasetToRaster_conversion(beLayer, beReturnsMinTempFile, interpolation_type = 'BINNING MINIMUM NONE', sampling_type = 'CELLSIZE', sampling_value = float(demListVal), data_type = 'FLOAT')
#                     log.debug('---Creating LR Min cm raster')
#                     beReturnsMinCm = Int(Times(beReturnsMin, 100))
#                     if bareEarthReturnMinFile is not None:
#                         bareEarthReturnMinFile_sized = updateResolution(bareEarthReturnMinFile, named_cell_size, demListVal, pattern22, log)
#                         log.debug('---Saving LR Min cm raster')
#                         beReturnsMinCm.save(bareEarthReturnMinFile_sized)#locDict['bareEarthReturnMinFile'])#.replace('fr', 'be'))
#                         log.debug('---Adding metadata to LR Min cm raster')
#                         addMetadata(bareEarthReturnMinFile_sized, paraDict, derivative_metadata, log)

#                 if int1rMaxFile is not None or int1rMinFile is not None or intBeMaxFile is not None:
#                     if demListVal == '2':# only run for 2m DEMs, otherwise too slow
#                         log.debug('---Creating intensity rasters')
#                         if int1rMaxFile is None and int1rMinFile is not None: 
#                             log.warning('Faking int1rMaxFile value due to requested int1rMinFile')
#                             int1rMaxFile_faked = int1rMinFile.replace('fr_int_min', 'fr_int_max')
#                             int1rMaxFile = int1rMaxFile_faked
#                         elif int1rMaxFile is None and intBeMaxFile is not None:
#                             log.warning('Faking int1rMaxFile value due to requested intBeMaxFile')
#                             int1rMaxFile_faked = intBeMaxFile.replace('be_int_max', 'fr_int_max')
#                             int1rMaxFile = int1rMaxFile_faked
#                         log.debug('---Creating FR Max Intensity')
#                         recode_tf = False
#                         log.debug(f'ir.max: {internal_regions.maximum},ir.min: {internal_regions.minimum}')
#                         int1rMaxFile_sized = updateResolution(int1rMaxFile, named_cell_size, demListVal, pattern22, log)
#                         if internal_regions.maximum - internal_regions.minimum != 0:
#                             log.info('multiple regions')
#                             # int1rMaxFile_sized_temp = opj(os.path.dirname(int1rMaxFile_sized), 'temp_' + os.path.basename(int1rMaxFile_sized))
#                             int1rMaxFile_sized_temp = os.path.join(procDir, '_'.join(['tmp_frmax', demPtString + 'm', huc12, 'out.tif']))
#                             lasd1rMaxIntensity = arcpy.LasDatasetToRaster_conversion(lasdAll, int1rMaxFile_sized_temp, 'INTENSITY', 'BINNING MAXIMUM NONE', sampling_type = 'CELLSIZE', sampling_value = float(demListVal), data_type = 'INT')
#                             int_zs_max = ZonalStatistics(internal_regions, 'VALUE', int1rMaxFile_sized_temp, 'MAXIMUM')
#                             if int_zs_max.minimum < 256 and int_zs_max.maximum > 256:
#                                 int_lt_256 = LessThan(int_zs_max, 256)
#                                 recode_areas = ZonalStatistics(internal_regions, 'VALUE', int_lt_256, 'MAXIMUM')
#                                 multiplied_intensities = Raster(int1rMaxFile_sized_temp) * 256
#                                 recoded_intensities = Con(recode_areas, multiplied_intensities, int1rMaxFile_sized_temp)
#                                 if int1rMaxFile is not None:
#                                     recoded_intensities.save(int1rMaxFile_sized)
#                                     arcpy.Delete_management(int1rMaxFile_sized_temp)
#                                 recode_tf = True
#                             else:
#                                 if int1rMaxFile is not None:
#                                     log.info('all regions equal max intensity')
#                                     arcpy.CopyRaster_management(int1rMaxFile_sized_temp, int1rMaxFile_sized)
#                                     arcpy.Delete_management(int1rMaxFile_sized_temp)
#                         else:
#                             log.info('one region')
#                             if int1rMaxFile is not None:
#                                 lasd1rMaxIntensity = arcpy.LasDatasetToRaster_conversion(lasdAll, int1rMaxFile_sized, 'INTENSITY', 'BINNING MAXIMUM NONE', sampling_type = 'CELLSIZE', sampling_value = float(demListVal), data_type = 'INT')

#                         if int1rMaxFile is not None:
#                             addMetadata(int1rMaxFile_sized, paraDict, derivative_metadata, log)

#                         if int1rMinFile is not None:
#                             log.debug('---Creating FR Min Intensity')
#                             int1rMinFile_sized = updateResolution(int1rMinFile, named_cell_size, demListVal, pattern22, log)
#                             if recode_tf:
#                                 # int1rMinFile_sized_temp = opj(os.path.dirname(int1rMinFile_sized), 'temp_' + os.path.basename(int1rMaxFile_sized))
#                                 int1rMinFile_sized_temp = os.path.join(procDir, '_'.join(['tmp_frmin', demPtString + 'm', huc12, 'out.tif']))
#                                 lasd1rMinIntensity = arcpy.LasDatasetToRaster_conversion(lasdAll, int1rMinFile_sized_temp, 'INTENSITY', 'BINNING MINIMUM NONE', sampling_type = 'CELLSIZE', sampling_value = float(demListVal), data_type = 'INT')
#                                 multiplied_intensities = Raster(int1rMinFile_sized_temp) * 256
#                                 recoded_intensities = Con(recode_areas, multiplied_intensities, int1rMinFile_sized_temp)
#                                 recoded_intensities.save(int1rMinFile_sized)
#                                 arcpy.Delete_management(int1rMinFile_sized_temp)
#                             else:
#                                 lasd1rMinIntensity = arcpy.LasDatasetToRaster_conversion(lasdAll, int1rMinFile_sized, 'INTENSITY', 'BINNING MINIMUM NONE', sampling_type = 'CELLSIZE', sampling_value = float(demListVal), data_type = 'INT')
#                             addMetadata(int1rMinFile_sized, paraDict, derivative_metadata, log)

#                         if intBeMaxFile is not None:
#                             log.debug('---Creating BE Max Intensity')
#                             intBeMaxFile_sized = updateResolution(intBeMaxFile, named_cell_size, demListVal, pattern22, log)
#                             if recode_tf:
#                                 # intBeMaxFile_sized_temp = opj(os.path.dirname(intBeMaxFile_sized), 'temp_' + os.path.basename(intBeMaxFile_sized))
#                                 intBeMaxFile_sized_temp = os.path.join(procDir, '_'.join(['tmp_bemax', demPtString + 'm', huc12, 'out.tif']))
#                                 lasdBeMaxIntensity = arcpy.LasDatasetToRaster_conversion(beLayer, intBeMaxFile_sized_temp, 'INTENSITY', 'BINNING MAXIMUM NONE', sampling_type = 'CELLSIZE', sampling_value = float(demListVal), data_type = 'INT')
#                                 multiplied_intensities = Raster(intBeMaxFile_sized_temp) * 256
#                                 recoded_intensities = Con(recode_areas, multiplied_intensities, intBeMaxFile_sized_temp)
#                                 recoded_intensities.save(intBeMaxFile_sized)
#                             else:
#                                 lasdBeMaxIntensity = arcpy.LasDatasetToRaster_conversion(beLayer, intBeMaxFile_sized, 'INTENSITY', 'BINNING MAXIMUM NONE', sampling_type = 'CELLSIZE', sampling_value = float(demListVal), data_type = 'INT')
#                             addMetadata(intBeMaxFile_sized, paraDict, derivative_metadata, log)

#                 if surfaceElevFile is not None:
#                     log.debug('---Creating FR Max surface')
#                     frMaxFile_sized = updateResolution(surfaceElevFile, named_cell_size, demListVal, pattern22, log)
#                     allReturnsMaxTempFile = os.path.join(procDir, '_'.join(['tmp_frmax', demPtString + 'm', huc12, 'out.tif']))
#                     # allReturnsMax = arcpy.LasDatasetToRaster_conversion(lasdAll, allReturnsMaxTempFile, interpolation_type = 'BINNING MAXIMUM SIMPLE', sampling_type = 'CELLSIZE', sampling_value = float(demListVal), data_type = 'FLOAT')
#                     allReturnsMax = arcpy.LasDatasetToRaster_conversion(lasdAll, allReturnsMaxTempFile, interpolation_type = 'BINNING MAXIMUM NONE', sampling_type = 'CELLSIZE', sampling_value = float(demListVal), data_type = 'FLOAT')
#                     allReturnsMaxCm = Int(Times(allReturnsMax, 100))
#                     allReturnsMaxCm.save(frMaxFile_sized)#locDict['surfaceElevFile'])#allReturnsMaxFile)
#                     addMetadata(frMaxFile_sized, paraDict, derivative_metadata, log)

#                 # log.debug('---Creating FR Min surface')
#                 # allReturnsMinTempFile = os.path.join(procDir, '_'.join(['tmp_frmin', demListVal + 'm', huc12, 'out.tif']))
#                 # allReturnsMin = arcpy.LasDatasetToRaster_conversion(lasdAll, allReturnsMinTempFile, interpolation_type = 'BINNING MINIMUM SIMPLE', sampling_type = 'CELLSIZE', sampling_value = float(demListVal), data_type = 'FLOAT')
#                 # allReturnsMinCm = Int(Times(allReturnsMin, 100))
#                 # allReturnsMinCm.save(frMinFile_sized)#locDict['firstReturnMinFile'])#allReturnsMinFile)

#                 if demListVal == '2':# only run for 2m DEMs, otherwise too slow
#                     if cnt1rFile is not None:
#                         log.debug('---Counting First Returns')
#                         cnt1rFile_sized = updateResolution(cnt1rFile, named_cell_size, demListVal, pattern22, log)
#                         cfrFileTemp = 'cnt_fr_' + demPtString + "m_" + huc12 + srSfx + '.tif'
#                         lasdCount = arcpy.LasPointStatsAsRaster_management(lasdAll, os.path.join(procDir, cfrFileTemp), 'POINT_COUNT', 'CELLSIZE', sampling_value = demListVal)
#                         cfrFileRasterObj = clipCountRaster(lasdCount, maskRastOut, cnt1rFile_sized)
#                         addMetadata(cnt1rFile_sized, paraDict, derivative_metadata, log)

#                     if cntPlsFile is not None:
#                         log.debug('---Counting All Returns')
#                         cntPlsFile_sized = updateResolution(cntPlsFile, named_cell_size, demListVal, pattern22, log)
#                         cntPlsFileTemp = 'cnt_pls_' + demPtString + "m_" + huc12 + srSfx + '.tif'
#                         lasdCount = arcpy.LasPointStatsAsRaster_management(lasdAll, os.path.join(procDir, cntPlsFileTemp), 'PULSE_COUNT', 'CELLSIZE', sampling_value = demListVal)
#                         cntPlsFileRasterObj = clipCountRaster(lasdCount, maskRastOut, cntPlsFile_sized)
#                         addMetadata(cntPlsFile_sized, paraDict, derivative_metadata, log)

#                     if cntBeFile is not None:
#                         log.debug('---Counting BE Returns')
#                         cntBeFile_sized = updateResolution(cntBeFile, named_cell_size, demListVal, pattern22, log)#.replace('_be_', '_belas_')
#                         cntBeFileTempSize = 'cnt_be_laspsr_' + demPtString + "m_" + huc12 + srSfx + '.tif'
#                         be_lasdCount = arcpy.LasPointStatsAsRaster_management(beLayer, os.path.join(procDir, cntBeFileTempSize), 'PULSE_COUNT', 'CELLSIZE', sampling_value = demListVal)
#                         # save the count raster with nulls converted to zeros, so that metadata can be added
#                         be_nulls = IsNull(be_lasdCount)#cntPlsFileRasterObj)
#                         be_nulls_as_zero = Con(be_nulls, 0, be_lasdCount)#cntPlsFileRasterObj)
#                         be_cntPlsFileRasterObj = clipCountRaster(be_lasdCount, maskRastOut, cntBeFile_sized)
# ##                        be_nulls_as_zero.save(cntBeFile_sized) 
#                         addMetadata(cntBeFile_sized, paraDict, derivative_metadata, log)

        else:
            log.warning('lidar data area does not exist or does not exceed build threshold; DEM was not built')

##----------------------------------------------------------------------

        # cleanup
        if cleanup:
##            if 'terrains' in locals():
##                dismantleTerrains(terrains, finalHb, finalNoZHb, poorZHb, finalHl, tcdFdSet, log)
            df.cleanupOther(procDir, log, sgdb, inm)

        return tElevFile

    except AssertionError:
        log.warning('assertion failure on: ' + huc12)
        sys.exit(1)

    except:
        # Get the traceback object
        #
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]

        # Concatenate information together concerning the error into a message string
        #
        pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
        msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"

        # Print Python error messages for use in Python / Python Window
        #
        log.warning(pymsg)
        log.warning(msgs)

        log.warning('failure on: ' + huc12)
        sys.exit(1)

    finally:
        log.info("Finished at ")
        handlers = log.handlers
        for h in handlers:
            log.info('shutting it down!')
            log.removeHandler(h)
            h.close()


##----------------------------------------------------------------------
## below should be commented out when using as a Python Toolbox (.pyt) - in 2025, .pyt cannot handle running code in the main block
## remove the comments below for use from the windows command line

# if __name__ == "__main__":
#     if len(sys.argv) == 1:
#         #Paste arguments into here for use within Python Window
#         arcpy.AddMessage("Whoo, hoo! Running from Python Window!")
#         # cleanup = False
#         parameters = 	["C:/Program Files/ArcGIS/Pro/bin/Python/envs/arcgispro-py3/pythonw.exe",
#     "C:/Temp/cmd_build_terrain_derivatives_main.pyt",
#     "E:/DEP_Checkout/Man_Data_ACPF/dep_ACPF2023/07080103/idepACPF070801030408.gdb/buf070801030408",
#     "E:/DEP_Checkout/Man_Data_ACPF/dep_ACPF2023/07080103/idepACPF070801030408.gdb/wesm_tiles_2026_02_10_070801030408",
#     "E:/DEP_Proc/DEMProc/LAS_dem2013_1m_070801030408/laz_all_points",
#     "C:/Users/bkgelder/.conda/envs/pdal_python/Library/bin/pdal.exe",
#     "1,2,3",
#     "E:/DEP_Proc/DEMProc/LAS_dem2013_1m_070801030408",
#     "E:/DEP_Checkout/Basedata_Summaries/Basedata_26915.gdb/Snap1m",
#     "",
#     "",
#     "E:/DEP_Checkout/LiDAR_Current/elev_TLib_mean18/07080103/et_1m_070801030408.tif",
#     "E:/DEP_Checkout/LiDAR_Current/surf_el_Lib/07080103/be_min_1m_070801030408.tif",
#     "E:/DEP_Checkout/LiDAR_Current/surf_el_Lib/07080103/fr_max_1m_070801030408.tif",
#     "E:/DEP_Checkout/LiDAR_Current/count_Lib/07080103/cnt_be_1m_070801030408.tif",
#     "E:/DEP_Checkout/LiDAR_Current/count_Lib/07080103/cnt_fr_1m_070801030408.tif",
#     "E:/DEP_Checkout/LiDAR_Current/count_Lib/07080103/cnt_pls_1m_070801030408.tif",
#     "E:/DEP_Checkout/LiDAR_Current/int_Lib/07080103/fr_int_min_1m_070801030408.tif",
#     "E:/DEP_Checkout/LiDAR_Current/int_Lib/07080103/fr_int_max_1m_070801030408.tif",
#     "E:/DEP_Checkout/LiDAR_Current/int_Lib/07080103/be_int_max_1m_070801030408.tif",
#     "False"]
#         for i in parameters[2:]:
#             sys.argv.append(i)

#     else:
#         #For use via Windows Command Line
#         #above 'parameters' come in via command line arguments, nothing else needed
#         arcpy.AddMessage("Whoo, hoo! Command-line enabled!")
#         # clean up the folder after done processing
#         # cleanup = True

#     # inputs then outputs, change "" to Python None
#     (dem_boundary, wesm_huc12_tiles, laz_download_dir,
#          pdal_exe, gsds, procDir, snap, breakpolys, breaklines, 
#          tElevFile, bareEarthReturnMinFile, firstReturnMaxFile, cntBeFile, cnt1rFile, cntPlsFile,
#          int1rMinFile, int1rMaxFile, intBeMaxFile, cleanup
#         ) = [i if i != "" else None for i in sys.argv[1:]]

#     # switch a text 'True' into a real Python True
#     cleanup = True if cleanup == "True" else False

#     messages = msgStub()

#     doLidarDEMs(dem_boundary, wesm_huc12_tiles, laz_download_dir,
#          pdal_exe, gsds, procDir, snap, breakpolys, breaklines, 
#          tElevFile, bareEarthReturnMinFile, firstReturnMaxFile, cntBeFile, cnt1rFile, cntPlsFile,
#          int1rMinFile, int1rMaxFile, intBeMaxFile, cleanup, messages)

#     arcpy.AddMessage("Back from doEPT!")

