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
import platform
import glob
import traceback
import re
from os.path import join as opj
from pathlib import Path

sys.path.append("C:\\DEP\\Scripts\\basics")

import dem_functions as df



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
            displayName="EPT WESM Merged Features",
            datatype="DEFeatureClass",
            parameterType='Required',
            direction="Input")
        
        param1 = arcpy.Parameter(
            name="dem_polygon",
            displayName="Buffered HUC12 Feature",
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
            displayName="Integer Resolution/Ground Sample Distance (meters) of output rasters, multiples joined by comma",
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
            name = "fElevFile",
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
        param19.value = True
                        
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
        doLidarDEMs(parameters[0].valueAsText, parameters[1].valueAsText, parameters[2].valueAsText, parameters[3].valueAsText, parameters[4].valueAsText, 
                    parameters[5].valueAsText, parameters[6].valueAsText, parameters[7].valueAsText, parameters[8].valueAsText, parameters[9].valueAsText, 
                    parameters[10].valueAsText, parameters[11].valueAsText, parameters[12].valueAsText, parameters[13].valueAsText, parameters[14].valueAsText, 
                    parameters[15].valueAsText, parameters[16].valueAsText, parameters[17].valueAsText, parameters[18].valueAsText, parameters[19].valueAsText, messages)
        return

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return

##----------------------------------------------------------------------
## Set environments and begin

def fillOCSinks(inDEM):
    # Return a raster will all one-cell-sinks filled
##    arcpy.AddMessage("-----Find Pits...")
    sinkFDir = FlowDirection(inDEM)
    allSinks = Sink(sinkFDir)
    arcpy.BuildRasterAttributeTable_management(allSinks)
##    log.warning('sinks for ' + str(inDEM) + ' is ' + str(int(arcpy.GetCount_management(allSinks).getOutput(0))))
##    arcpy.AddMessage("-----Fill everything else...")

    ## Make a No-one-cell-sink DEM
    AllButSinks_DEM = Con(IsNull(allSinks), inDEM)

    ## Fill the No-one-cell-sink DEM
    absDEM_fill = Fill(AllButSinks_DEM)

    ## Add the Original 'real' sinks back into the filled DEM
    fill_DEM = Con(IsNull(absDEM_fill), inDEM, absDEM_fill)

    return(fill_DEM, allSinks)


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

        log.debug('--- Done creating LAS dataset and extracting LAS at ' + time.asctime())

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
        assert int(arcpy.GetCount_management(dem_polygon).getOutput(0)) < 2, 'multiple features in feature class'
        maskFc = arcpy.CopyFeatures_management(dem_polygon)
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

        log.info("maskFcPrelim complete at: " + time.asctime())

        ## Set up geodatabase to store the multipoint files and terrains (necessary all inputs be in feature dataset
        # Vertical units are in meters (float) so use a meter-based reference
        FDSet = arcpy.CreateFeatureDataset_management(sgdb, "Lidar_pts", srOut)
        maskFcOut = projIfNeeded(maskFc, os.path.join(str(FDSet), 'buf_huc' + srSfx), srOut)
        log.info("maskFcOut complete at: " + time.asctime())

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

def createCountsFromMultipoints(sgdb, maskRastBase, demList, huc12, finalMPinm, finalMP, log, cntBeFile, named_cell_size, pattern):#locDict):
    try:
        maskRastOut = opj(sgdb, maskRastBase + str(demList[0]))

        cntBeFile_sized = updateResolution(cntBeFile, named_cell_size, demList[0], pattern, log)
        log.debug('---Counting Bare Earth Returns for ' + str(demList[0]))
        terrCountName = arcpy.ValidateTableName("cnt" + str(demList[0]) + "m_fl_" + huc12, sgdb)
        terrCount = arcpy.PointToRaster_conversion(finalMPinm, arcpy.Describe(finalMP).OIDFieldName, os.path.join(sgdb, terrCountName), "COUNT", "NONE", str(demList[0]))
        cntBeFileRasterObj = clipCountRaster(terrCount, maskRastOut, cntBeFile_sized)

        return cntBeFileRasterObj

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


def buildTerrains(finalMP, FDSet, tcdFdSet, finalHb, finalHl, finalNoZHb, poorZHb, log, windows, time):
    terrains = []
    # create one terrain with ZMinMax option
    tp = ["1", "", "", "WINDOWSIZE", "ZMEAN", "MILD", 0.18]
    pyrmd_str = "2 1000;4 2500;8 5000;16 10000;32 20000;64 40000"

####    windowsizeMethods = windows#['ZMEAN', 'ZMINMAX']
    for window in windows:
        log.info(f'creating and setting up terrain: {window} at {time.asctime()}')
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

def createRastersFromTerrains(log, demList, procDir, terrains, huc12, lidar_metadata_info, pyramid_args, dem_metadata_template):
    try:
        # log.debug('snapRaster for Terrain to raster: ' + arcpy.env.snapRaster)
        interpTechnique = 'NATURAL_NEIGHBORS'
        pyramidLevel = '4'
        dem_cellSize = demList[0]
        arcpy.env.cellSize = dem_cellSize
        for terrain in terrains:
            pfFileTemp = os.path.join(procDir, '_'.join(['tmp_ter', terrain.getInput(6), str(dem_cellSize) + 'm', huc12, 'out.tif']))
            log.debug('---Creating Raster from Terrain for ' + str(dem_cellSize) + ' using ' + terrain.getInput(6))
            log.debug(f'Creating Raster from Terrain at {pfFileTemp}')
            # tempTerrName = generateTempTerrName(procDir, terrain.getInput(6), dem_cellSize, huc12)
            demOut = arcpy.TerrainToRaster_3d(terrain, pfFileTemp, "FLOAT", interpTechnique, "CELLSIZE " + str(dem_cellSize), pyramidLevel)

            ttr_arguments_list = [demOut.getInput(t) for t in range(0,5)]
            ttr_args = ', '.join(ttr_arguments_list)
            log.debug(f'TerrainToRaster args: {ttr_args}')
####            demList.append(demOut.getOutput(0))

            nowYmd, collect_starts_min, collect_ends_max, collect_majority = [i for i in lidar_metadata_info]

            terrain_args = terrain_args_from_inputs(terrain)

            log.debug(f"terrain building args: {terrain_args}")

            paraDict = {
                '\n\nACPF: DEM Generation and Pit Fill Tool     ' : '\nRun Date: %s' % nowYmd,
                # '\nUnknown Vintage Lidar Data: ' : False,#tiles_t_or_f,
                '\nEarliest 3DEP Lidar Data: ' : collect_starts_min,
                '\nLatest 3DEP Lidar Data: ' : collect_ends_max,
                '\nLatest 3DEP Lidar Data: ' : collect_majority,
                '\nOutput conditioned DEM raster: ' : pfFileTemp,#fElevFile_interp,
                '\nTerrain Interpolation Arguments: ' : terrain_args,
                '\nTerrain Pyramid Arguments: ' : pyramid_args,
                '\nTerrain To Raster Arguments: ' : ttr_args
                }

            ## update metadata
            log.debug('---Adding metadata')
            addMetadata(pfFileTemp, paraDict, dem_metadata_template, log)



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
        log.info(f'surfcons are: {surfcons} at {time.asctime()}')
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


def buildLASRasters(lasdAll, lasdGround, log, demList, huc12, srSfx, maskRastBase, sgdb, procDir, int1rMaxFile, int1rMinFile, surfaceElevFile, intBeMaxFile, bareEarthReturnMinFile, cnt1rFile, cntPlsFile, named_cell_size, internal_regions, lidar_metadata_info, derivative_metadata, pattern22):
##def buildLASRasters(lasdAll, lasdGround, log, demList, huc12, srSfx, maskRastBase, sgdb, procDir, int1rMaxFile, int1rMinFile, surfaceElevFile, frMinFile, intBeMaxFile, intBeMinFile, lastReturnMinFile, bareEarthReturnMinFile, cnt1rFile, named_cell_size, int_regions, ptr):
    '''creates multiple rasters from a las dataset, including min/max intensity of
    first return and bare earth surfaces, first return max and min surface, and z_range'''
    try:
        maskRastOut = opj(sgdb, maskRastBase + str(demList[0]))

        # log.debug('snapRaster for LAS Dataset to raster: ' + arcpy.env.snapRaster)

        # procDir = locDict['fProcDir']
##        frMinFile_sized = updateResolution(frMinFile, named_cell_size, demList[0], huc12, log)
##        intBeMinFile_sized = updateResolution(intBeMinFile, named_cell_size, demList[0], huc12, log)
##        lastReturnMinFile_sized = updateResolution(lastReturnMinFile, named_cell_size, demList[0], huc12, log)


        nowYmd, collect_starts_min, collect_ends_max, collect_majority = [i for i in lidar_metadata_info]

        paraDict = {
                '\n\nACPF: DEM Generation and Pit Fill Tool     ' : '\nRun Date: %s' % nowYmd,
                # '\nUnknown Vintage Lidar Data: ' : False,#tiles_t_or_f,
                '\nEarliest 3DEP Lidar Data: ' : collect_starts_min,
                '\nLatest 3DEP Lidar Data: ' : collect_ends_max,
                '\nLatest 3DEP Lidar Data: ' : collect_majority
                }

        if bareEarthReturnMinFile is not None or intBeMaxFile is not None:
            log.debug('---Creating LR Min surface')
            if sys.version_info.minor < 9:
                beLayer = arcpy.MakeLasDatasetLayer_management(lasdGround, 'ground_layer', [2,8], 'Last Return')
            else:
                beLayer = arcpy.MakeLasDatasetLayer_management(lasdGround, 'ground_layer', [2,8], 'LAST')

            beReturnsMinTempFile = os.path.join(procDir, '_'.join(['tmp_bemin', str(demList[0]) + 'm', huc12, 'out.tif']))
            beReturnsMin = arcpy.LasDatasetToRaster_conversion(beLayer, beReturnsMinTempFile, interpolation_type = 'BINNING MINIMUM NONE', sampling_type = 'CELLSIZE', sampling_value = demList[0], data_type = 'FLOAT')
            beReturnsMinCm = Int(Times(beReturnsMin, 100))
            if bareEarthReturnMinFile is not None:
                bareEarthReturnMinFile_sized = updateResolution(bareEarthReturnMinFile, named_cell_size, demList[0], pattern22, log)
                beReturnsMinCm.save(bareEarthReturnMinFile_sized)#locDict['bareEarthReturnMinFile'])#.replace('fr', 'be'))
                addMetadata(bareEarthReturnMinFile_sized, paraDict, derivative_metadata, log)

        if int1rMaxFile is not None or int1rMinFile is not None or intBeMaxFile is not None:
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
            int1rMaxFile_sized = updateResolution(int1rMaxFile, named_cell_size, demList[0], pattern22, log)
            if internal_regions.maximum - internal_regions.minimum != 0:
                log.info('multiple regions')
                # int1rMaxFile_sized_temp = opj(os.path.dirname(int1rMaxFile_sized), 'temp_' + os.path.basename(int1rMaxFile_sized))
                int1rMaxFile_sized_temp = os.path.join(procDir, '_'.join(['tmp_frmax', str(demList[0]) + 'm', huc12, 'out.tif']))
                lasd1rMaxIntensity = arcpy.LasDatasetToRaster_conversion(lasdAll, int1rMaxFile_sized_temp, 'INTENSITY', 'BINNING MAXIMUM NONE', sampling_type = 'CELLSIZE', sampling_value = demList[0], data_type = 'INT')
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
                    lasd1rMaxIntensity = arcpy.LasDatasetToRaster_conversion(lasdAll, int1rMaxFile_sized, 'INTENSITY', 'BINNING MAXIMUM NONE', sampling_type = 'CELLSIZE', sampling_value = demList[0], data_type = 'INT')

            if int1rMaxFile is not None:
                addMetadata(int1rMaxFile_sized, paraDict, derivative_metadata, log)

            if int1rMinFile is not None:
                log.debug('---Creating FR Min Intensity')
                int1rMinFile_sized = updateResolution(int1rMinFile, named_cell_size, demList[0], pattern22, log)
                if recode_tf:
                    # int1rMinFile_sized_temp = opj(os.path.dirname(int1rMinFile_sized), 'temp_' + os.path.basename(int1rMaxFile_sized))
                    int1rMinFile_sized_temp = os.path.join(procDir, '_'.join(['tmp_frmin', str(demList[0]) + 'm', huc12, 'out.tif']))
                    lasd1rMinIntensity = arcpy.LasDatasetToRaster_conversion(lasdAll, int1rMinFile_sized_temp, 'INTENSITY', 'BINNING MINIMUM NONE', sampling_type = 'CELLSIZE', sampling_value = demList[0], data_type = 'INT')
                    multiplied_intensities = Raster(int1rMinFile_sized_temp) * 256
                    recoded_intensities = Con(recode_areas, multiplied_intensities, int1rMinFile_sized_temp)
                    recoded_intensities.save(int1rMinFile_sized)
                    arcpy.Delete_management(int1rMinFile_sized_temp)
                else:
                    lasd1rMinIntensity = arcpy.LasDatasetToRaster_conversion(lasdAll, int1rMinFile_sized, 'INTENSITY', 'BINNING MINIMUM NONE', sampling_type = 'CELLSIZE', sampling_value = demList[0], data_type = 'INT')
                addMetadata(int1rMinFile_sized, paraDict, derivative_metadata, log)

            if intBeMaxFile is not None:
                log.debug('---Creating BE Max Intensity')
                intBeMaxFile_sized = updateResolution(intBeMaxFile, named_cell_size, demList[0], pattern22, log)
                if recode_tf:
                    # intBeMaxFile_sized_temp = opj(os.path.dirname(intBeMaxFile_sized), 'temp_' + os.path.basename(intBeMaxFile_sized))
                    intBeMaxFile_sized_temp = os.path.join(procDir, '_'.join(['tmp_bemax', str(demList[0]) + 'm', huc12, 'out.tif']))
                    lasdBeMaxIntensity = arcpy.LasDatasetToRaster_conversion(beLayer, intBeMaxFile_sized_temp, 'INTENSITY', 'BINNING MAXIMUM NONE', sampling_type = 'CELLSIZE', sampling_value = demList[0], data_type = 'INT')
                    multiplied_intensities = Raster(intBeMaxFile_sized_temp) * 256
                    recoded_intensities = Con(recode_areas, multiplied_intensities, intBeMaxFile_sized_temp)
                    recoded_intensities.save(intBeMaxFile_sized)
                else:
                    lasdBeMaxIntensity = arcpy.LasDatasetToRaster_conversion(beLayer, intBeMaxFile_sized, 'INTENSITY', 'BINNING MAXIMUM NONE', sampling_type = 'CELLSIZE', sampling_value = demList[0], data_type = 'INT')
                addMetadata(intBeMaxFile_sized, paraDict, derivative_metadata, log)

        if surfaceElevFile is not None:
            log.debug('---Creating FR Max surface')
            frMaxFile_sized = updateResolution(surfaceElevFile, named_cell_size, demList[0], pattern22, log)
            allReturnsMaxTempFile = os.path.join(procDir, '_'.join(['tmp_frmax', str(demList[0]) + 'm', huc12, 'out.tif']))
            allReturnsMax = arcpy.LasDatasetToRaster_conversion(lasdAll, allReturnsMaxTempFile, interpolation_type = 'BINNING MAXIMUM SIMPLE', sampling_type = 'CELLSIZE', sampling_value = demList[0], data_type = 'FLOAT')
            allReturnsMaxCm = Int(Times(allReturnsMax, 100))
            allReturnsMaxCm.save(frMaxFile_sized)#locDict['surfaceElevFile'])#allReturnsMaxFile)
            addMetadata(frMaxFile_sized, paraDict, derivative_metadata, log)

        # log.debug('---Creating FR Min surface')
        # allReturnsMinTempFile = os.path.join(procDir, '_'.join(['tmp_frmin', str(demList[0]) + 'm', huc12, 'out.tif']))
        # allReturnsMin = arcpy.LasDatasetToRaster_conversion(lasdAll, allReturnsMinTempFile, interpolation_type = 'BINNING MINIMUM SIMPLE', sampling_type = 'CELLSIZE', sampling_value = demList[0], data_type = 'FLOAT')
        # allReturnsMinCm = Int(Times(allReturnsMin, 100))
        # allReturnsMinCm.save(frMinFile_sized)#locDict['firstReturnMinFile'])#allReturnsMinFile)

        if cnt1rFile is not None:
            log.debug('---Counting First Returns')
            cnt1rFile_sized = updateResolution(cnt1rFile, named_cell_size, demList[0], pattern22, log)
            cfrFileTemp = 'cnt_fr_' + str(demList[0]) + "m_" + huc12 + srSfx + '.tif'
            lasdCount = arcpy.LasPointStatsAsRaster_management(lasdAll, os.path.join(procDir, cfrFileTemp), 'POINT_COUNT', 'CELLSIZE', demList[0])
            cfrFileRasterObj = clipCountRaster(lasdCount, maskRastOut, cnt1rFile_sized)
            addMetadata(cnt1rFile_sized, paraDict, derivative_metadata, log)

        if cntPlsFile is not None:
            log.debug('---Counting All Returns')
            cntPlsFile_sized = updateResolution(cntPlsFile, named_cell_size, demList[0], pattern22, log)
            cntPlsFileTemp = 'cnt_pls_' + str(demList[0]) + "m_" + huc12 + srSfx + '.tif'
            lasdCount = arcpy.LasPointStatsAsRaster_management(lasdAll, os.path.join(procDir, cntPlsFileTemp), 'PULSE_COUNT', 'CELLSIZE', demList[0])
            cntPlsFileRasterObj = clipCountRaster(lasdCount, maskRastOut, cntPlsFile_sized)
            addMetadata(cntPlsFile_sized, paraDict, derivative_metadata, log)

        # log.debug('---Counting Z Range')
        # zrangeFileTemp = 'zrng_all_' + str(demList[0]) + "m_" + huc12 + srSfx + '.tif'
        # lasdCount = arcpy.LasPointStatsAsRaster_management(lasdAll, os.path.join(procDir, zrangeFileTemp), 'Z_RANGE', 'CELLSIZE', demList[0])

##        lastLayer = arcpy.MakeLasDatasetLayer_management(lasdOut, 'ground_layer', [2,8], 'Last Return')
        # if sys.version_info.minor < 9:
        #     lastLayer = arcpy.MakeLasDatasetLayer_management(lasdGround, 'last_layer', return_values = 'Last Return')
        # else:
        #     lastLayer = arcpy.MakeLasDatasetLayer_management(lasdGround, 'last_layer', return_values = 'LAST')
        # log.debug('---Creating BE Max Intensity')
        # # minimum be seems to be a little less 'noisy' than maximum
        # lasdBeMaxIntensity = arcpy.LasDatasetToRaster_conversion(lastLayer, intBeMaxFile_sized, 'INTENSITY', 'BINNING MAXIMUM NONE', sampling_type = 'CELLSIZE', sampling_value = demList[0], data_type = 'INT')

        # log.debug('---Creating BE Min Intensity')
        # # minimum be seems to be a little less 'noisy' than maximum
        # lasdBeMinIntensity = arcpy.LasDatasetToRaster_conversion(lastLayer, intBeMinFile_sized, 'INTENSITY', 'BINNING MINIMUM NONE', sampling_type = 'CELLSIZE', sampling_value = demList[0], data_type = 'INT')

        # log.debug('---Creating LR Min surface')
        # lastReturnsMinTempFile = os.path.join(procDir, '_'.join(['tmp_lrmin', str(demList[0]) + 'm', huc12, 'out.tif']))
        # lastReturnsMin = arcpy.LasDatasetToRaster_conversion(lastLayer, lastReturnsMinTempFile, interpolation_type = 'BINNING MINIMUM SIMPLE', sampling_type = 'CELLSIZE', sampling_value = demList[0], data_type = 'FLOAT')
        # lastReturnsMinCm = Int(Times(lastReturnsMin, 100))
        # lastReturnsMinCm.save(lastReturnMinFile_sized)#locDict['lastReturnMinFile'])#.replace('fr', 'lr'))

##        # my current favorites  - minmax mild 18 terrain, binning min simple (almost too detailed),
##        interps = ['BINNING MINIMUM SIMPLE', 'BINNING AVERAGE SIMPLE', 'BINNING IDW SIMPLE', 'TRIANGULATION NATURAL_NEIGHBOR WINDOW_SIZE MINIMUM ' + str(demList[0]), 'TRIANGULATION LINEAR WINDOW_SIZE MINIMUM ' + str(demList[0]), 'TRIANGULATION NATURAL_NEIGHBOR WINDOW_SIZE CLOSEST_TO_MEAN ' + str(demList[0]), 'TRIANGULATION LINEAR WINDOW_SIZE CLOSEST_TO_MEAN ' + str(demList[0])]
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
##            beReturnsMinTempFile = os.path.join(procDir, '_'.join(['tmp_be', interpString, str(demList[0]) + 'm', huc12, 'out.tif']))
##            beReturnsMin = arcpy.LasDatasetToRaster_conversion(beLayer, beReturnsMinTempFile, interpolation_type = interp, sampling_type = 'CELLSIZE', sampling_value = demList[0], data_type = 'FLOAT')
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
    tempTerrName = os.path.join(procDir, '_'.join(['tmp_ter', window, str(cellsize) + 'm', huc12, 'out.tif']))

    return tempTerrName


def mosaicDEMsAndPitfill(demList, maskRastBase, huc12, log, sgdb, terrains, procDir, fElevFile, interpDict, named_cell_size,
                          srOutNoVCS, dem_metadata_template, lidar_metadata_info, pattern22, pyramid_args):
    '''Takes whole or partial DEMs from the demList and mosaics them together
    if there are multiple DEMs. Then pit-fills the result (fills all one cell
    sinks). Also processes 'ZMEAN' and 'ZMINMAX' (or other terrain->raster
    options if they're added to the code) to separate directories.'''

    try:
        arcpy.env.cellSize = demList[0]
        maskRastOut = opj(sgdb, maskRastBase + str(demList[0]))

        noLASdem = demList[1:]#[]
        log.debug('noLASdem (if present) is: ' + str(noLASdem))

        # windows are types of terrain (ZMEAN, ZMINMAX, etc.)
        # if len(windows):
        for terrain in terrains:#window in windows:
            window = terrain.getInput(6)
            log.debug(f"processing window: {window}")

            fElevFile = updateResolution(fElevFile, named_cell_size, demList[0], pattern22, log)
            log.debug(f"pit filling for {fElevFile}")

            interpType = interpDict[window]
            # default interpolation type is mean18
            if interpType != 'mean18':
                fElevFile_interp = fElevFile.replace('mean18', interpType)
            else:
                fElevFile_interp = fElevFile

            tempTerrName = generateTempTerrName(procDir, window, demList[0], huc12)
            # get rasters created from Terrain, should only be one
            arcpy.env.workspace = os.path.dirname(tempTerrName)
            rastrList = arcpy.ListRasters(os.path.basename(tempTerrName))#[0]

            # if len(noLASdem) > 0:
            #     mosaicList = rastrList + noLASdem
            #     if True:#len(mosaicList) > 1:
            #         log.debug('resampling method for mosaic to new raster: ' + arcpy.env.resamplingMethod)
            #         log.debug('cellsize for mosaic to new raster: ' + arcpy.env.cellSize)
            #         log.debug('---Mosaicing DEMs for ' + str(demList[0]) + ')
            #         mosDEM = arcpy.MosaicToNewRaster_management(mosaicList, procDir, 'mos_' + str(demList[0]) + 'm_' + huc12 + '.tif', number_of_bands = '1', pixel_type = '32_BIT_FLOAT', mosaic_method = 'MINIMUM')
            #         maskedDEM = Con(Plus(maskRastOut, 2), mosDEM)
            #         maskedDEMintPre = Int(maskedDEM*100)

            #         # filter any null values in mask (buffered HUC12) using Nibble
            #         isnCmDemInt = IsNull(maskedDEM)
            #         ndToNibble = SetNull(isnCmDemInt == 1, 1)
            #         cmDEMintNoNulls = Con(isnCmDemInt == 0, maskedDEMintPre, 1)
            #         nibbleInt = Nibble(cmDEMintNoNulls, ndToNibble)
            #         maskedDEMint = Con(Plus(maskRastOut, 2), nibbleInt)

            #     else:
            #         log.debug('resampling method for Resample: ' + arcpy.env.resamplingMethod)
            #         log.debug('cellsize for resample: ' + arcpy.env.cellSize)
            #         log.debug('---Resampling DEMs for ' + str(demList[0]) + ')
            #         mosDEM = arcpy.Resample_management(mosaicList[0], os.path.join(procDir, 'mos_' + str(demList[0]) + 'm_' + huc12 + '.tif'), cell_size = arcpy.env.cellSize, resampling_type = arcpy.env.resamplingMethod)
            #         maskedDEM = Con(Plus(maskRastOut, 2), mosDEM)
            #         maskedDEMintPre = Int(maskedDEM*100)

            #         # filter any null values in mask (buffered HUC12) using Nibble
            #         isnCmDemInt = IsNull(maskedDEM)
            #         ndToNibble = SetNull(isnCmDemInt == 1, 1)
            #         cmDEMintNoNulls = Con(isnCmDemInt == 0, maskedDEMintPre, 1)
            #         nibbleInt = Nibble(cmDEMintNoNulls, ndToNibble)
            #         maskedDEMint = Con(Plus(maskRastOut, 2), nibbleInt)

            # else:
            rastr = rastrList[0]
            mosDEM = Raster(rastr)
            intCmDEM = Int(mosDEM*100)
            maskedDEMint = ExtractByMask(intCmDEM, maskRastOut)

            # clear output VCS from here on out, creating centimeter rasters
            arcpy.env.outputCoordinateSystem = None
            arcpy.env.outputCoordinateSystem = srOutNoVCS

            log.debug('msk_dem to be saved now')
            maskedDEMint.save('msk_dem_' + str(demList[0]) + 'm_' + huc12 + '.tif')


            cmDEMnocs, cmDEMsinks = fillOCSinks(maskedDEMint)
            log.debug('pfFile name will be ' + fElevFile_interp)#paths['fElevFile'])
            cmDEMnocs.save(fElevFile_interp)#paths['fElevFile'])
            log.debug('Saved DEM for ' + str(demList[0]))# + ' at ' + str(time.clock()))
            arcpy.BuildPyramids_management(cmDEMnocs)

            f_dict = {}
            sep = ': '
            f_metadata = md.Metadata(rastr)
            if f_metadata.summary is not None and f_metadata.summary != '':
                f_met_split = f_metadata.summary.split('\n')
                for i in f_met_split[2:]:
                    key, value = i.split(sep, 1)
                    print(f'key {key} and value {value}')
                    f_dict.update({key: value})

            log.debug(f"terrain building args: {f_dict}")

            ## update metadata
            log.debug('---Adding metadata')
            addMetadata(fElevFile_interp, f_dict, dem_metadata_template, log)

    except Exception as e:
        print('handling as exception')
##        log.debug(e.message)
        if sys.version_info.major == 2:
            arcpy.AddError(e.message)
            print(e.message)
            log.warning(e.message)
        elif sys.version_info.major == 3:
            arcpy.AddError(e)
            print(e)
            log.warning(e)

        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]

        # Concatenate information together concerning the error into a message string
        pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
        # Return python error messages for use in script tool or Python Window
        arcpy.AddError(pymsg)
        # Print Python error messages for use in Python / Python Window
        print(pymsg + "\n")
        log.warning(pymsg)

        if arcpy.GetMessages(2) not in pymsg:
            msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
            arcpy.AddError(msgs)
            print(msgs)
            log.warning(msgs)

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
        log.warning(pymsg)

        if arcpy.GetMessages(2) not in pymsg:
            msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
            arcpy.AddError(msgs)
            print(msgs)
            log.warning(msgs)


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


def addMetadata(outDEM, paraDict, template_file_path, log = None):
    # Set the standard-format metadata XML file's path
    # need to load metadata editor via 'import arcpy.metadata as md'
    # outDEM = raster to receive updated metadata
    # paraDict = dictionary of key/value pairs to be stored in metadata
    #   values stored include things like analyst, lidar acquisition date, etc.
    # template_file_path = a template to load a basic summary from
    # log = otional logging of error messages to a log file
    # scriptPath = sys.path[0]
    try:
        src_file_path = template_file_path

        # Get the target item's Metadata object
        tgt_item_md = md.Metadata(outDEM)    

        # Import the ACPF metadata content to the target item
        if not tgt_item_md.isReadOnly:
            tgt_item_md.importMetadata(src_file_path)
            tgt_item_md.title = os.path.split(outDEM)[1]
            tgt_item_md.credits = 'Analyst: %s' % os.getlogin()#getpass.getuser()

            src_desc = tgt_item_md.summary
            if src_desc == None:
                src_desc = ''
            for key, value in paraDict.items():  
                src_desc = src_desc + ('%s %s' % (key, value))
            tgt_item_md.summary = src_desc
            
            tgt_item_md.save()

    except TypeError as e:
        print('handling as exception')
##        log.debug(e.message)
        if sys.version_info.major == 2:
            arcpy.AddError(e.message)
            print(e.message)
            log.warning(e.message)
        elif sys.version_info.major == 3:
            arcpy.AddError(e)
            print(e)
            if log is not None:
                log.warning(e)

        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]

        # Concatenate information together concerning the error into a message string
        pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
        # Return python error messages for use in script tool or Python Window
        arcpy.AddError(pymsg)
        # Print Python error messages for use in Python / Python Window
        print(pymsg + "\n")
        if log is not None:
            log.warning(pymsg)

        if arcpy.GetMessages(2) not in pymsg:
            msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
            arcpy.AddError(msgs)
            print(msgs)
            if log is not None:
                log.warning(msgs)

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
        if log is not None:
            log.warning(pymsg)

        if arcpy.GetMessages(2) not in pymsg:
            msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
            arcpy.AddError(msgs)
            print(msgs)
            if log is not None:
                log.warning(msgs)

def create_cl2_json_pipeline(cl2_json_filename, eleDir, all_las_file, cl2_las_full_filename):
    '''Writes a json pipeline for use by pdal (point data abstraction library)'''

    cl2_json_full_filename = os.altsep.join([eleDir, cl2_json_filename])

    json_str = '''{
"pipeline": [
{
    "filename": "''' + all_las_file + '''",
    "type": "readers.las",
    "tag": "readdata"
},
{
    "tag": "getclass2",
    "type": "filters.range",
    "limits": "Classification[2:2]"
},
{
    "filename": "''' + cl2_las_full_filename + '''",
    "tag": "writerslas",
    "type": "writers.las"
}]}'''

    json_file_obj = open(cl2_json_full_filename, 'w')
    json_file_obj.write(json_str)
    json_file_obj.close()

    return cl2_json_full_filename


def create_laz_json_pipeline(laz_json_filename, saveDir, all_las_file, laz_full_filename):
    '''Writes a json pipeline for use by pdal (point data abstraction library)'''

    laz_json_full_filename = os.altsep.join([saveDir, laz_json_filename])

    json_str = '''{
"pipeline": [
{
    "filename": "''' + all_las_file + '''",
    "type": "readers.las",
    "tag": "readdata"
},
{
    "filename": "''' + laz_full_filename + '''",
    "tag": "writerslas",
    "type": "writers.las"
}]}'''

    json_file_obj = open(laz_json_full_filename, 'w')
    json_file_obj.write(json_str)
    json_file_obj.close()

    return laz_json_full_filename


def create_ept_json_pipeline(ept_json_filename, eleDir, ept_las_full_filename, extent_request, ept_address, srOutCode):
    '''Writes a json pipeline for use by pdal (point data abstraction library)'''

    ept_json_full_filename = os.altsep.join([eleDir, ept_json_filename])

    json_str = '''{
"pipeline": [
{
    "bounds": "([''' + extent_request + '''])",
    "filename": "''' + ept_address + '''",
    "type": "readers.ept",
    "tag": "readdata"
},
{
    "out_srs": "EPSG:''' + str(srOutCode) + '''",
    "tag": "reprojectUTM",
    "type": "filters.reprojection"
},
{
    "filename": "''' + ept_las_full_filename + '''",
    "tag": "writerslas",
    "type": "writers.las"
}]}'''

    json_file_obj = open(ept_json_full_filename, 'w')
    json_file_obj.write(json_str)
    json_file_obj.close()

    return ept_json_full_filename

def organizeProjectsByDate(wesm_huc12, work_id_name, maskFc_area, build_threshold, log):
    """Organize the 3DEP projects by acquisition date so we use the
    most recent data first, then fill in with older data"""
    valid_order = arcpy.ValidateFieldName('order', os.path.dirname(wesm_huc12.getOutput(0)))
    addOrderField = arcpy.AddField_management(wesm_huc12, valid_order, 'SHORT')

    # wesm_fields =
    ordered_work_ids = [s[0] for s in arcpy.da.SearchCursor(wesm_huc12, [work_id_name], sql_clause = [None, 'ORDER BY collect_start DESC'])]
    merged_area = 0
    for cnt, o in enumerate(ordered_work_ids):
        if merged_area <= maskFc_area * build_threshold:
            log.debug(o)
            if o < 0:
                selected = arcpy.Select_analysis(wesm_huc12, 'select_wkid_neg' + str(abs(o)), work_id_name + ' = ' + str(o))
            else:
                selected = arcpy.Select_analysis(wesm_huc12, 'select_wkid_' + str(o), work_id_name + ' = ' + str(o))
            arcpy.CalculateField_management(selected, valid_order, cnt+1, 'PYTHON3')
            if cnt > 0:
                erased = arcpy.Erase_analysis(selected, prev_merged)
                if df.testForZero(erased):
                    merged = arcpy.Merge_management([erased, prev_merged])
                    merged_area += [s[0] for s in arcpy.da.SearchCursor(erased, ['SHAPE@AREA'])][0]
            else:
                merged = selected
                merged_area += [s[0] for s in arcpy.da.SearchCursor(selected, ['SHAPE@AREA'])][0]
            prev_merged = merged
            log.debug(f"merged_area: {merged_area}")
    log.info(f'merged_area was: {merged_area} and maskFc_area was: {maskFc_area}')

    # if merged_area >= maskFc_area * 0.9999:
    #     if 'tcdFdSet' not in locals():
    #         tcdFdSet = arcpy.CopyFeatures_management(maskFc, os.path.join(str(FDSet), 'ept_las_area'))

#             bigVoids = arcpy.Erase_analysis(maskFcOut, tcdFdSet_local)
#             wesm_clipped = arcpy.analysis.Clip(wesm_huc12, bigVoids)

    return prev_merged, merged_area, addOrderField

def queryParts(geom, geom_extent, maskFcOut, srOut, sgdb, log):#maskFc_3857, maskFc_3857_desc)
    """Subdivide the EPT request into one or four parts based on size"""
    parts = []
        # if more than 500 sq km, split into 4
    x_range = geom_extent.XMax - geom_extent.XMin
    y_range = geom_extent.YMax - geom_extent.YMin
    square_area = (x_range)*(y_range)
    log.debug(f'Square area in km^2 is: {round(square_area/pow(1000,2), 1)}')
    log.debug(f'Geometry area in km^2 is: {round(geom.area/pow(1000,2), 1)}')

    if square_area / (1000**2) > 500.0:
        splits = 4
        # switch to do intersect/clip in 3857
        arcpy.env.outputCoordinateSystem = 3857
        for p in range(0,splits):
            p_div = p // pow(splits, 0.5)
            p_mod = p % pow(splits, 0.5)
            x_start = geom_extent.XMin + x_range * (p_div + 0) / pow(splits, 0.5)
            x_end = geom_extent.XMin + x_range * (p_div + 1) /pow(splits, 0.5)
            y_start = geom_extent.YMin + y_range * (p_mod + 0) /pow(splits, 0.5)
            y_end = geom_extent.YMin + y_range * (p_mod + 1) /pow(splits, 0.5)
            log.debug(f'p_div: {p_div}, p_mod: {p_mod}')
            log.debug(f'x_start: {x_start} and x_end: {x_end}')
            log.debug(f'y_start: {y_start} and y_end: {y_end}')
            # refine request by intersecting preliminary extent polygon with HUC12 buffered boundary - helps in NE-SW and SE-NW watersheds
            ext1 = arcpy.Extent(x_start, y_start, x_end, y_end)
            poly1 = ext1.polygon

            maskFc_3857 = arcpy.management.Project(maskFcOut, opj(sgdb, 'maskfc_3857'), 3857)
            maskFc_3857_desc = arcpy.da.Describe(maskFc_3857)
            maskFc_3857_extent = maskFc_3857_desc['extent']
            # Should always be False
            #poly1.disjoint(maskFc_3857_extent)
            poly1_3857_int = poly1.intersect(maskFc_3857_extent, 4)

            clip3 = arcpy.analysis.Clip(maskFc_3857, poly1_3857_int)
            clip3_extent = arcpy.da.Describe(clip3)['extent']

            x_start = clip3_extent.XMin
            x_end = clip3_extent.XMax
            y_start = clip3_extent.YMin
            y_end = clip3_extent.YMax

            log.debug(f'final_x_start: {x_start} and x_end: {x_end}')
            log.debug(f'final_y_start: {y_start} and y_end: {y_end}')
            ept_extent = str(x_start) + ', ' + str(x_end) + '], [' + str(y_start) + ', ' + str(y_end)

            parts.append(['_pt' + str(p), ept_extent])

        arcpy.env.outputCoordinateSystem = srOut

    else:
        splits = 1
        ept_extent = str(geom.extent.XMin) + ', ' + str(geom.extent.XMax) + '], [' + str(geom.extent.YMin) + ', ' + str(geom.extent.YMax)
        parts.append(['', ept_extent])

    return parts, square_area


def getLidarFiles(wesm_huc12, work_id_name, pdal_exe, prev_merged, addOrderField, log, sgdb, sfldr, srOut, srOutCode, huc12, eleDir, maskFcOut, fixedFolder, inm, FDSet, allTilesList, procDir):
    try:
        if df.testForZero(prev_merged):
            # requests to EPT must be in 3857
            prev_merged_projected_3857 = arcpy.management.Project(prev_merged, 'proj_trial', 3857)
            # max_area = 0
            url_list = df.getfields(wesm_huc12, 'url*')
            with arcpy.da.SearchCursor(prev_merged_projected_3857, ['SHAPE@', work_id_name] + url_list, sql_clause = [None, 'ORDER BY ' + addOrderField.getInput(1) + ' DESC']) as scur:#work_id_name, 'SHAPE@AREA', 'lpc_link']) as scur:
            # with arcpy.da.SearchCursor(prev_merged_projected_3857, ['SHAPE@', work_id_name] + url_list, work_id_name + ' = -1117', sql_clause = [None, 'ORDER BY ' + addOrderField.getInput(1) + ' DESC']) as scur:#work_id_name, 'SHAPE@AREA', 'lpc_link']) as scur:
            # with arcpy.da.SearchCursor(prev_merged_projected_3857, ['SHAPE@', work_id_name, 'SHAPE@AREA', 'lpc_link']) as scur:
                for srow in scur:
                    log.debug(f'{work_id_name} is: {srow[1]}')
                    geom = srow[0]
                    geom_extent = geom.extent
                    las_size_threshold = 750 #bytes, then assume .las file has points
                    parts, square_area = queryParts(geom, geom_extent, maskFcOut, srOut, sgdb, log)#maskFc_3857, maskFc_3857_desc)

                    arcpy.env.outputCoordinateSystem = srOut
                    for part in parts:
                        work_id = srow[1]
                        work_id_part = str(srow[1]) + part[0]
                        extent_request = part[1]

                        geom_srOut = geom.projectAs(arcpy.SpatialReference(srOutCode))
                        valid_geom_name = arcpy.ValidateTableName('geom_proj_' + work_id_part, sgdb)
                        geom_srOut_copy = arcpy.CopyFeatures_management(geom_srOut, valid_geom_name)

                        # get the first non-None url, that will tell us address of EPT.JSON
                        urls = srow[1:]
                        for u in urls:
                            if u is not None:
                                url = u
                        ept_address = url

                        ept_zlas_filename = "_".join(["ept", huc12, str(work_id_part) + ".zlas"])
                        ept_zlas_full_filename = os.altsep.join([eleDir.replace(os.path.sep, os.path.altsep), ept_zlas_filename])

                        ept_las_filename = ept_zlas_filename.replace(".zlas", ".las")

                        ept_laz_filename = ept_zlas_filename.replace(".zlas", ".laz")
                        ept_laz_full_filename = os.altsep.join([eleDir.replace(os.path.sep, os.path.altsep), ept_laz_filename])
                        local_ept_laz_full_filename = os.altsep.join([procDir.replace(os.path.sep, os.path.altsep), ept_laz_filename])

                        # pipeline json requires / not \ for path separator
                        ept_las_full_filename = os.altsep.join([procDir.replace(os.path.sep, os.path.altsep), ept_las_filename])
                        # NEEDS UPDATE - check for files in lidar_download_directory as well - getting ready below
                        # dl_ept_laz_full_filename = ept_laz_full_filename.replace(eptDir, eleDir)
                        # dl_ept_zlas_full_filename = ept_zlas_full_filename.replace(eptDir, eleDir)
                        if os.path.isfile(ept_laz_full_filename) and not os.path.isfile(ept_las_full_filename):
                            log.info('converting laz to las')
                            log.info(f"arguments: {ept_laz_full_filename}, {procDir}")
                            las_result = arcpy.conversion.ConvertLas(ept_laz_full_filename, procDir)#, compression = 'ZLAS')
                            log.info(las_result)
                        elif os.path.isfile(ept_zlas_full_filename) and not os.path.isfile(ept_las_full_filename):
                            log.info('converting zlas to las')
                            log.info(f"arguments: {ept_zlas_full_filename}, {procDir}")
                            las_result = arcpy.conversion.ConvertLas(ept_zlas_full_filename, procDir)#, compression = 'ZLAS')
                            log.info(las_result)
                        # if zlas does not exist, get las then convert to zlas
                        if not os.path.isfile(ept_zlas_full_filename):
                            log.info('Getting LAS from EPT')
                            log.info(ept_zlas_filename)

                            ept_json_filename = "_".join(["get", "ept", huc12, str(work_id_part) + ".json"])

                            df.create_needed_dirs_and_gdbs(ept_las_full_filename, log)
                            ept_json_full_filename = create_ept_json_pipeline(ept_json_filename, eleDir, ept_las_full_filename, extent_request, ept_address, srOutCode)
                            df.create_needed_dirs_and_gdbs(ept_json_full_filename, log)

                            if not os.path.exists(ept_las_full_filename):
                                run_string = " ".join([pdal_exe, "pipeline", ept_json_full_filename])
                                # estimate download time based on 102500040309 (area 1175 km2) in 4 parts
                                m2_per_sec = 1175.2*1000**2/len(parts)/2200
                                log.debug(f'pdal run_string: {run_string}')
                                log.info(f'Estimated pdal download time (for QL2 lidar): {round(square_area/(m2_per_sec * len(parts) * 60), 2)} minutes for {ept_json_filename}')
                                co = subprocess.run(run_string)
                                log.debug(f'completed pdal run_string')

                            # archive as zlas for use later in this script and re-use
                            stats = os.stat(ept_las_full_filename)
                            if stats.st_size > las_size_threshold:
                                # log.info('converting las to zlas for archive')
                                # zlas_result = arcpy.conversion.ConvertLas(ept_las_full_filename, ele, compression = 'ZLAS', las_options = None)
                                # log.info(zlas_result)

                                # log.debug('converting las to laz for archive')
                                # laz_json_full_filename = create_laz_json_pipeline(ept_laz_filename, procDir, ept_las_full_filename, ept_laz_full_filename)
                                # laz_run_string = " ".join([pdal_exe, "pipeline", laz_json_full_filename])
                                # log.debug(f'pdal run_string: {laz_run_string}')
                                # co = subprocess.run(laz_run_string)

                                laz_result = arcpy.conversion.ConvertLas(ept_las_full_filename, eleDir, compression = 'LAZ', las_options = None)
                                log.debug(laz_result)

                                # arcpy.Delete_management(ept_las_full_filename)
                            else:
                                log.warning(f"{ept_las_full_filename} has very small file size; plotting extent as poly")
                                poly35 = geom_extent.polygon
                                p35 = arcpy.management.CopyFeatures(poly35, opj(sgdb, 'failed_' + valid_geom_name))
                                ## Use to get requested bounds
                                # geom_extent.JSON
                                # '{"xmin":-11386571.4549,"ymin":5310909.4501999989,"xmax":-11370238.307800001,"ymax":5314697.4098999985,"spatialReference":{"wkid":102100,"latestWkid":3857}}'

                                ## Use to plot query bounds from PDAL Debug or from EPT.JSON file header
                                # arcpy.env.outputCoordinateSystem = 3857
                                # ept_json_extent = arcpy.Extent(-11583422,5262814,-11396830,5449406)
                                # ept_json_extent_polygon = ept_json_extent.polygon
                                # ept_json_polygon = arcpy.management.CopyFeatures(ept_json_extent_polygon, opj(sgdb, 'poly37'))
                        else:
                            if os.path.isfile(ept_las_full_filename):
                                stats = os.stat(ept_las_full_filename)
                            else:
                                log.warning("can't get good LAS data, skipping to next project")
                                continue

                        # if co.returncode == 0:
                        if os.path.exists(ept_las_full_filename) and stats.st_size > las_size_threshold:
                            cl2Las = processEptLas(sgdb, sfldr, srOutCode, fixedFolder, geom_srOut, ept_las_full_filename, srOut, inm, FDSet, procDir, allTilesList, log, time, work_id)
                            #remove invalid geometry
                            if cl2Las is None:
                                log.warning('deleting ' + str(geom_srOut_copy))
                                delete = arcpy.Delete_management(geom_srOut_copy)
                                # remove from dates
                                prev_merged = arcpy.Select_analysis(prev_merged, where_clause = work_id_name + ' <> ' + str(work_id))
                            else:
                                cl2LasNotNone = cl2Las
                        elif os.path.exists(ept_las_full_filename):
                            log.warning('EPT LAS file too small, no points')
                            run_string += ' --debug'
                            log.warning(f'Invalid call was {run_string}')

                        else:
                            log.warning('No valid output from ept.json request')
                            log.warning(f'Invalid call was {run_string}')
                            sys.exit(1)

        #revert to previous good cl2Las if last is None
        if cl2Las is None:
            cl2Las = cl2LasNotNone

        return cl2Las, geom_srOut_copy

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


def try_to_delete(rasRes, log):
    if arcpy.Exists(rasRes):
        try:
            arcpy.Delete_management(rasRes)
        except arcpy.ExecuteError:
            log.warning('could not remove using arcpy.Delete, trying os.remove')
            os.remove(rasRes)


def doLidarDEMs(monthly_wesm_ept_mashup, dem_polygon, 
         pdal_exe, gsds, procDir, snap, breakpolys, breaklines, 
         fElevFile, bareEarthReturnMinFile, firstReturnMaxFile, cntBeFile, cnt1rFile, cntPlsFile,
         int1rMinFile, int1rMaxFile, intBeMaxFile, ept_wesm_project_file, lidar_download_directory, cleanup, messages):
    
    arguments = [monthly_wesm_ept_mashup, dem_polygon, 
        pdal_exe, gsds, procDir, snap, breakpolys, breaklines, 
        fElevFile, bareEarthReturnMinFile, firstReturnMaxFile, cntBeFile, cnt1rFile, cntPlsFile,
        int1rMinFile, int1rMaxFile, intBeMaxFile, ept_wesm_project_file, lidar_download_directory, cleanup]
                                                       
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
        huc12, huc8 = df.figureItOut(fElevFile)
        # the DEP huc DEM naming convention
        pattern27 = 'e[cfpxv][0-9]m\\d{10,16}'
        pattern26 = '[0-9]m\\d{10,16}'
        pattern25 = '[0-9]m_d{10,16}'
        pattern22 = '[+_][0-9]m[+_]'
        # huc12, huc8, named_cell_size = df.figureItOut(fElevFile)

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

        if lidar_download_directory is None:
            lidar_download_directory = procDir

        sgdb = arcpy.env.scratchGDB
        arcpy.env.scratchWorkspace = sgdb
        arcpy.env.workspace = sgdb

        if procDir is None:
            procDir = sfldr

        if lidar_download_directory is None:
            lidar_download_directory = sfldr

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
        log.info("Beginning execution: " + time.asctime())
        log.debug('sys.argv is: ' + str(sys.argv) + '\n')
        log.info("Processing HUC: " + huc12)
        log.info(f"procDir: {procDir}")
        log.info(f"lidar_download_directory: {lidar_download_directory}")

        log.info(f"procDir: {procDir}")
        log.info(f"lidar_download_directory: {lidar_download_directory}")

        fElevDesc = arcpy.da.Describe(dem_polygon)
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

        if init_res + 'm' not in fElevFile:
            fElevFile = os.path.splitext(fElevFile)[0] + '_' + str(init_res) + 'm' + os.path.splitext(fElevFile)[1]

        ## windowsizeMethods are the criterion used to select which point(s) in the window define the terrain
        interpDict = df.loadInterpDict()
        windowsizeMethods = ['ZMEAN', 'ZMINMAX']#, 'ZMIN']
        interpType = interpDict[windowsizeMethods[0]]

        if interpType not in fElevFile:
            fElevFile = os.path.splitext(fElevFile)[0] + '_' + interpType + os.path.splitext(fElevFile)[1]
        log.debug(f'Revised fElevFile: {fElevFile}')

        # delete any pre-existing inputs
        for ras in [fElevFile, cntBeFile, int1rMaxFile, int1rMinFile, cnt1rFile, firstReturnMaxFile]:
            if ras is not None:
                filename_path = Path(ras)
                # check to see if it follows HUC DEM naming procedure
                stem = filename_path.stem
                if re.search(pattern22, stem):
                    if len(demLists) > 0:
                        for demList in demLists:
                            rasRes = updateResolution(ras, init_res, demList[0], pattern22, log)
                            try_to_delete(rasRes, log)
                    # if str(named_cell_size) + 'm' in ras:
                    #         # updateResolution(filepath, init_res, new_res, huc12, log)
                    #         rasRes = ras.replace(str(named_cell_size) + 'm', str(demList[0]) + 'm')
                    #         try_to_delete(rasRes, log)
                    # else:
                    #     try_to_delete(rasRes, log)

        # create output directories
        for filename in [fElevFile, cntBeFile, int1rMaxFile, firstReturnMaxFile]:
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

        maskRastBase = 'mask_rast_'
        maskFc, maskFc_area, maskFcOut, maskRastOut, hucRastOut, FDSet = prepPolygonBoundary(dem_polygon, log, sgdb, srOut, srSfx, maskRastBase, demLists)
        
# Build the DEM using one of two ways
# First see if there is any lidar LAS data (preferred)
# Second, use lidar derived DEMs by county/tile to build out the rest

##----------------------------------------------------------------------

    # # check for collection change (different priorities) to restrict further data
        fixedFolder = str(arcpy.CreateFolder_management(sfldr, 'fixed'))
        localLidarFolder = str(arcpy.CreateFolder_management(sfldr, 'localLidar'))
        clippedFolder = str(arcpy.CreateFolder_management(sfldr, 'clipped'))
        projectedFolder = str(arcpy.CreateFolder_management(sfldr, 'projected'))

        allTilesList = []

        #eptDir = os.path.dirname(os.path.dirname(monthly_wesm_ept_mashup))
        wesm_huc12 = arcpy.analysis.Clip(monthly_wesm_ept_mashup, maskFcOut, 'wesm_' + huc12)

##----------------------------------------------------------------------


##----------------------------------------------------------------------
        work_id_name = 'workunit_id'
        build_threshold = 0.9999
        if df.testForZero(wesm_huc12):
            prev_merged, merged_area, addOrderField = organizeProjectsByDate(wesm_huc12, work_id_name, maskFc_area, build_threshold, log)

            cl2Las, geom_srOut_copy = getLidarFiles(wesm_huc12, work_id_name, pdal_exe, prev_merged, addOrderField, log, sgdb, sfldr, srOut, srOutCode, huc12, lidar_download_directory, maskFcOut, fixedFolder, inm, FDSet, allTilesList, procDir)

            arcpy.env.outputCoordinateSystem = srOut

            log.debug('tileList is ' + str(allTilesList))
            # save the merged wesm_las_dates

            arcpy.AddField_management(prev_merged, 'area_field', 'DOUBLE')
            arcpy.CalculateField_management(prev_merged, 'area_field', '!shape.area!', 'Python 3')

            ptr = arcpy.conversion.PolygonToRaster(prev_merged, work_id_name, opj(sgdb, 'ptr'))
            int_ptr = Int(ptr)
            fs1 = FocalStatistics(int_ptr, NbrRectangle(5,5), 'RANGE')
            internal_regions = Con(fs1 == 0, int_ptr)
            internal_regions.save(opj(sgdb, 'intrnl_rgns'))
            
            if ept_wesm_project_file is not None:
                #assume not in a feature dataset...
                df.create_needed_dirs_and_gdbs(ept_wesm_project_file, log)
                # ept_wesm_project_gdb = os.path.dirname(ept_wesm_project_file)
                # if not arcpy.Exists(ept_wesm_project_gdb):
                #     if not os.path.isdir(os.path.dirname(ept_wesm_project_gdb)):
                #         os.makedirs(os.path.dirname(ept_wesm_project_gdb))
                #     log.debug(f"making gdb: {ept_wesm_project_gdb}")
                #     ept_gdb = arcpy.CreateFileGDB_management(os.path.dirname(ept_wesm_project_gdb), os.path.basename(ept_wesm_project_gdb))
                merged_copy = arcpy.CopyFeatures_management(prev_merged, ept_wesm_project_file)

            ept_lidar_fcs = arcpy.ListFeatureClasses(os.path.basename(geom_srOut_copy.getOutput(0))[:10] + '*')
            tcdFdSet_ept = arcpy.Union_analysis(ept_lidar_fcs)

            ##----------------------------------------------------------------------
            # now build datasets
            arcpy.env.workspace = sgdb
            mpList = arcpy.ListFeatureClasses('pts_*', feature_type = 'POINT', feature_dataset = os.path.basename(FDSet.getOutput(0)))
            if len(mpList) > 0:
                finalMP = arcpy.Merge_management(mpList, os.path.join(str(FDSet), 'mp_merge'))
                if df.testForZero(finalMP):
                    finalMPinm, finalHb, poorZHb, finalNoZHb, finalHl = setupPointsAndBreaklines(finalMP, inm, FDSet, breakpolys, breaklines, log)

                    # if 'tcdFdSet_local' in locals() and 'tcdFdSet_ept' in locals():
                    #     tcdFdSet_union = arcpy.Union_analysis([tcdFdSet_local, tcdFdSet_ept], os.path.join(str(FDSet), 'ept_and_local_las_union'))
                    # elif 'tcdFdSet_local' in locals():
                    #     tcdFdSet_union = tcdFdSet_local
                    # else:
                    tcdFdSet_union = tcdFdSet_ept

                    tcdFdSet = arcpy.management.Dissolve(tcdFdSet_union, os.path.join(str(FDSet), 'ept_and_local_las'))
                    fill_donut_slow(tcdFdSet)

                    terrains, tf, terrain_args, pyramid_args = buildTerrains(finalMP, FDSet, tcdFdSet, finalHb, finalHl, finalNoZHb, poorZHb, log, windowsizeMethods, time)

                    # cl2_tiles_list = filter_to_cl2()
                    check_4_all = os.altsep.join([os.path.dirname(cl2Las), '*' + cl2Las[-7:]])
                    all_tiles_list = glob.glob(check_4_all)
                    cl2_tiles_list = []
                    for all_tile in all_tiles_list:
                        # JSON likes / not \
                        all_tile_alt = all_tile.replace(os.sep, os.altsep)
                        work_id = os.path.basename(all_tile_alt).split('_')[2]

                        cl2_json_filename = "_".join(["filter", "fixed_cl2", huc12, work_id + ".json"])

                        cl2_las_full_filename = all_tile_alt.replace('.las', '_cl2.las')
                        cl2_tiles_list.append(cl2_las_full_filename)
                        
                        cl2_las_json = create_cl2_json_pipeline(cl2_json_filename, lidar_download_directory, all_tile_alt, cl2_las_full_filename)

                        cl2_run_string = " ".join([pdal_exe, "pipeline", cl2_las_json])
                        log.info(f"running: {cl2_run_string}")
                        co = subprocess.run(cl2_run_string)
                        if co.returncode != 0:
                            log.warning(f"Failed on {cl2_run_string}")

                    log.debug(f"cl2_tiles_list: {cl2_tiles_list}")
                    # assume lidar data in same spatial reference as output, ExtractLAS should handle that
                    lasdGround = setupLasDataset(cl2_tiles_list, tcdFdSet, procDir, None, None, srSfx, None, log, time, arcpy.SpatialReference(int(srOutCode)))
                    log.info('finished setting up las dataset')

                    collect_ends_max, collect_starts_min, collect_majority = getLidarTimeframes(prev_merged)#merged_copy)#, tilesClip_local)

                    lidar_metadata_info = [nowYmd, collect_starts_min, collect_ends_max, collect_majority]

                    lasdAll = arcpy.CreateLasDataset_management(allTilesList, os.path.join(procDir, 'huc_all.lasd'), spatial_reference = arcpy.SpatialReference(int(srOutCode)))
                    # ## Following code runs slowly at times and is not being used further 2023.12.21
                    # # classify overlap in lasdAll
                    # arcpy.ddd.ClassifyLasOverlap(lasdAll)
                    # overlapLayer = arcpy.MakeLasDatasetLayer_management(lasdAll, 'overlap_layer', [12], )
                    # overlapMaxIntensity = arcpy.LasDatasetToRaster_conversion(overlapLayer, opj(procDir, 'overlap' + str(demList[0])), 'INTENSITY', 'BINNING MAXIMUM NONE', sampling_type = 'CELLSIZE', sampling_value = demList[0], data_type = 'INT')
                    # # overlapCD = arcpy.

                    # overlap_args = getLasRasterArguments(overlapMaxIntensity)
                    # if 'merged_copy' not in locals():
                    #     merged_copy = None

                    # paraDict = {
                    #     '\n\nDEP: DEM Overlap Intensity Output     ' : '\nRun Date: %s' % nowYmd,
                    #     '\nUnknown Vintage Lidar Data' : False,
                    #     '\nEarliest 3DEP Lidar Data: ' : collect_starts_min,
                    #     '\nLatest 3DEP Lidar Data: ' : collect_ends_max,
                    #     '\nLas Dataset To Raster Arguments: ' : overlap_args
                    #     }

                    # template_interp = derivative_metadata

                    # ## update metadata
                    # log.debug('---Adding metadata')
                    # addMetadata(overlapMaxIntensity.getOutput(0), paraDict, template_interp, log)

                    # try:
                    #     overlap_cost_surface = Con(IsNull(overlapMaxIntensity), 1, 0.001)
                    #     arcpy.env.mask = Raster(maskRastOut)
                    #     overlap_cost_dist2 = DistanceAccumulation(overlapMaxIntensity, '', '', overlap_cost_surface)
                    #     gaps = Con(overlap_cost_dist2 > 2*demList[0], 1)
                    #     gaps_regions = RegionGroup(gaps)
                    #     arcpy.env.mask = Raster(hucRastOut)
                    #     gaps_regions_max2 = ZonalStatistics(gaps_regions, 'Value', overlap_cost_dist2, 'MAXIMUM')
                    #     below_flight_lines = Con(overlap_cost_dist2 > 0.75*gaps_regions_max2.maximum, 1)
                    # except:
                    #     log.warning('Failure on calculating flight lines from returns')
                    # arcpy.env.mask = None

        # indicate there are no datasets to process further
        else:
            prev_merged = None
            merged_area = 0.00

        arcpy.env.cellSize = None

##----------------------------------------------------------------------
        if merged_area / maskFc_area >= build_threshold:
            for demList in demLists:
                log.debug(f'---Processing resolution {demList}')
                # arcpy.env.snapRaster
                if cntBeFile is not None:
                    cntBeFileRasterObj = createCountsFromMultipoints(sgdb, maskRastBase, demList, huc12, finalMPinm, finalMP, log, cntBeFile, init_res, pattern22)

                terrainList = createRastersFromTerrains(log, demList, procDir, terrains, huc12, lidar_metadata_info, pyramid_args, flib_metadata_template)

                buildLASRasters(lasdAll, lasdGround, log, demList, huc12, srSfx, maskRastBase, sgdb, procDir, int1rMaxFile, int1rMinFile, firstReturnMaxFile, intBeMaxFile, bareEarthReturnMinFile, cnt1rFile, cntPlsFile, init_res, internal_regions, lidar_metadata_info, derivative_metadata, pattern22)

                mosaicDEMsAndPitfill(demList, maskRastBase, huc12, log, sgdb, terrains, procDir, fElevFile, interpDict, init_res, srOutNoVCS, flib_metadata_template, lidar_metadata_info, pattern22, pyramid_args)
        else:
            log.warning('lidar data area does not exist or does not exceed build threshold; DEM was not built')

##----------------------------------------------------------------------

        # cleanup
        if cleanup:
##            if 'terrains' in locals():
##                dismantleTerrains(terrains, finalHb, finalNoZHb, poorZHb, finalHl, tcdFdSet, log)
            df.cleanupOther(procDir, log, sgdb, inm)

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
        log.warning("Finished at " + time.asctime())
        handlers = log.handlers
        for h in handlers:
            log.info('shutting it down!')
            log.removeHandler(h)
            h.close()

    return fElevFile


##----------------------------------------------------------------------



# if __name__ == "__main__":
#     if len(sys.argv) == 1:
#         #Paste arguments into here for use within Python Window
#         arcpy.AddMessage("Whoo, hoo! Running from Python Window!")
#         # cleanup = False

#         parameters = 	["C:/Program Files/ArcGIS/Pro/bin/Python/envs/arcgispro-py3/pythonw.exe",
#     "C:/DEP/Scripts/basics/cmd_build_DEM_main.pyt",
#     "//dep2.ae.iastate.edu/D$/DEP/Elevation_databases/ept.gdb/ept_resources_2024_09_01",
#     "//dep2.ae.iastate.edu/D$/DEP/Man_Data_ACPF/dep_ACPF2022/07080105/idepACPF070801050307.gdb/buf_070801050307",
#     "C:/Users/bkgelder/.conda/envs/pdal/Library/bin/pdal.exe",
#     "5,3,2,1",
#     "D:/DEP_Proc/DEMProc/LAS_dem2013_2m_070801050307",
#     "//dep2.ae.iastate.edu/D$/DEP/Basedata_Summaries/Basedata_26915.gdb/Snap1m",
#     "",
#     "",
#     "//dep2.ae.iastate.edu/M$/DEP/LiDAR_Current/elev_FLib_mean18/07080105/ef_2m_070801050307.tif",
#     "//dep2.ae.iastate.edu/M$/DEP/LiDAR_Current/surf_el_Lib/07080105/be_min_2m_070801050307.tif",
#     "//dep2.ae.iastate.edu/M$/DEP/LiDAR_Current/surf_el_Lib/07080105/fr_max_2m_070801050307.tif",
#     "//dep2.ae.iastate.edu/M$/DEP/LiDAR_Current/count_Lib/07080105/cnt_be_2m_070801050307.tif",
#     "//dep2.ae.iastate.edu/M$/DEP/LiDAR_Current/count_Lib/07080105/cnt_fr_2m_070801050307.tif",
#     "//dep2.ae.iastate.edu/M$/DEP/LiDAR_Current/count_Lib/07080105/cnt_pls_2m_070801050307.tif",
#     "//dep2.ae.iastate.edu/M$/DEP/LiDAR_Current/int_Lib/07080105/fr_int_min_2m_070801050307.tif",
#     "//dep2.ae.iastate.edu/M$/DEP/LiDAR_Current/int_Lib/07080105/fr_int_max_2m_070801050307.tif",
#     "//dep2.ae.iastate.edu/M$/DEP/LiDAR_Current/int_Lib/07080105/be_int_max_2m_070801050307.tif",
#     "//dep2.ae.iastate.edu/D$/DEP/Man_Data_ACPF/dep_ACPF2022/07080105/idepACPF070801050307.gdb/wesm_ept_resources_2024_09_01_070801050307",
#     "//dep2.ae.iastate.edu/M$/DEP/Elev_Base_Data",
#     "True"]
#         for i in parameters[2:]:
#             sys.argv.append(i)

#     else:
#         #For use via Windows Command Line
#         #above 'parameters' come in via command line arguments, nothing else needed
#         arcpy.AddMessage("Whoo, hoo! Command-line enabled!")
#         # clean up the folder after done processing
#         # cleanup = True

#     # inputs then outputs, change "" to Python None
#     (monthly_wesm_ept_mashup, dem_polygon, 
#          pdal_exe, gsds, procDir, snap, breakpolys, breaklines, 
#          fElevFile, bareEarthReturnMinFile, firstReturnMaxFile, cntBeFile, cnt1rFile, cntPlsFile,
#          int1rMinFile, int1rMaxFile, intBeMaxFile, ept_wesm_project_file, lidar_download_directory, cleanup
#         ) = [i if i != "" else None for i in sys.argv[1:]]

#     # switch a text 'True' into a real Python True
#     cleanup = True if cleanup == "True" else False

#     messages = msgStub()

#     doLidarDEMs(monthly_wesm_ept_mashup, dem_polygon, 
#          pdal_exe, gsds, procDir, snap, breakpolys, breaklines, 
#          fElevFile, bareEarthReturnMinFile, firstReturnMaxFile, cntBeFile, cnt1rFile, cntPlsFile,
#          int1rMinFile, int1rMaxFile, intBeMaxFile, ept_wesm_project_file, lidar_download_directory, cleanup, messages)#msgStub())

#     arcpy.AddMessage("Back from doEPT!")