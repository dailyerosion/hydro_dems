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
import shutil
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
            displayName="Output Terrain-based Elevation Model",
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
        assert int(arcpy.GetCount_management(dem_polygon).getOutput(0)) < 2, 'multiple features in polygon feature class'
        assert int(arcpy.GetCount_management(dem_polygon).getOutput(0)) > 0, 'no features in polygon feature class'
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

        for demListVal in demLists:
            maskRastOut = arcpy.PolygonToRaster_conversion(maskFcOut, 'id', opj(sgdb, maskRastBase + demListVal), cellsize = int(demListVal))
            # huc_rast_out = arcpy.conversion.PolygonToRaster(geom_copy, 'OBJECTID', opj(sgdb, 'huc_rast' + demListVal), cellsize = demList[0])

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

def queryParts(geom, geom_extent, maskFcOut, srOut, sgdb, log, ql1):#maskFc_3857, maskFc_3857_desc)
    """Subdivide the EPT request into one or four parts based on size"""
    parts = []
        # if more than 500 sq km, split into 4
    x_range = geom_extent.XMax - geom_extent.XMin
    y_range = geom_extent.YMax - geom_extent.YMin
    square_area = (x_range)*(y_range)
    log.debug(f'Square area in km^2 is: {round(square_area/pow(1000,2), 1)}')
    log.debug(f'Geometry area in km^2 is: {round(geom.area/pow(1000,2), 1)}')

    if square_area / (1000**2) > 150.0 or ql1:
        if ql1:
            x_net_size = 2000
        else:
            x_net_size = 5000
        y_net_size = x_net_size
        log.debug(f'x/y_net_size for fishnet: {x_net_size}')

        maskfc_5070 = arcpy.Project_management(maskFcOut, opj(sgdb, 'maskfc_5070'), 5070)
        mask_5070_extent = arcpy.da.Describe(maskfc_5070)['extent']

        origin_coords = f"{mask_5070_extent.XMin} {mask_5070_extent.YMin}"
        y_axis_coords = f"{mask_5070_extent.XMin} {mask_5070_extent.YMin + y_net_size}"
        rows = ceil((mask_5070_extent.YMax - mask_5070_extent.YMin) / y_net_size)
        cols = ceil((mask_5070_extent.XMax - mask_5070_extent.XMin) / x_net_size)

        # the geometry extents are coming in 5070
        sr_5070 = arcpy.SpatialReference(5070)
        arcpy.env.outputCoordinateSystem = sr_5070

        fishnet_name = opj(sgdb, 'fishnet')
        log.debug(f"fishnet args are: {fishnet_name, [geom_extent.XMin, geom_extent.YMin], None, x_net_size, y_net_size}")
        # fishnet = arcpy.CreateFishnet_management(opj(sgdb, 'fishnet'), origin_coord =  + str(geom_extent.XMin, geom_extent.YMin], None, x_net_size, y_net_size)
        fishnet2 = arcpy.CreateFishnet_management(fishnet_name, origin_coords, y_axis_coords, x_net_size, y_net_size, rows, cols, geometry_type = "POLYGON")

        fishnet_int_mask = arcpy.Intersect_analysis([fishnet2, maskFcOut])

        pt_field_name = 'PT_Identifier'
        arcpy.AddField_management(fishnet_int_mask, pt_field_name, 'TEXT', field_length='10')
        arcpy.CalculateField_management(fishnet_int_mask, pt_field_name, '!FID_' + os.path.basename(fishnet_name) + '!', 'PYTHON3')

        # # PDAL requests must be in 3857
        # arcpy.env.outputCoordinateSystem = 3857
        fishnet_int_mask_3857 = arcpy.management.Project(fishnet_int_mask, opj(sgdb, 'fishnet_3857'), 3857)

        with arcpy.da.SearchCursor(fishnet_int_mask_3857, ['OID@', 'SHAPE@', pt_field_name]) as scur:
            for p, srow in enumerate(scur):
                oid = srow[0]
                geom_fish = srow[1]
                pt_code = srow[2]
                # splits.append(parts, geom_fish)

                clip3_extent = geom_fish.extent
                x_start = clip3_extent.XMin
                x_end = clip3_extent.XMax
                y_start = clip3_extent.YMin
                y_end = clip3_extent.YMax

                log.debug(f'fishnet_int_mask oid: {oid}')
                log.debug(f'final_x_start: {x_start} and x_end: {x_end}')
                log.debug(f'final_y_start: {y_start} and y_end: {y_end}')
                ept_extent = str(x_start) + ', ' + str(x_end) + '], [' + str(y_start) + ', ' + str(y_end)


                parts.append(['_pt' + str(pt_code), ept_extent])

        # # splits = 4 # for splitting manually, not fishnet
        # # switch to do intersect/clip in 3857
        # arcpy.env.outputCoordinateSystem = 3857
        # for p in range(0,splits):
        #     p_div = p // pow(splits, 0.5)
        #     p_mod = p % pow(splits, 0.5)
        #     x_start = geom_extent.XMin + x_range * (p_div + 0) / pow(splits, 0.5)
        #     x_end = geom_extent.XMin + x_range * (p_div + 1) /pow(splits, 0.5)
        #     y_start = geom_extent.YMin + y_range * (p_mod + 0) /pow(splits, 0.5)
        #     y_end = geom_extent.YMin + y_range * (p_mod + 1) /pow(splits, 0.5)
        #     log.debug(f'p_div: {p_div}, p_mod: {p_mod}')
        #     log.debug(f'x_start: {x_start} and x_end: {x_end}')
        #     log.debug(f'y_start: {y_start} and y_end: {y_end}')
        #     # refine request by intersecting preliminary extent polygon with HUC12 buffered boundary - helps in NE-SW and SE-NW watersheds
        #     ext1 = arcpy.Extent(x_start, y_start, x_end, y_end)
        #     poly1 = ext1.polygon

        #     maskFc_3857 = arcpy.management.Project(maskFcOut, opj(sgdb, 'maskfc_3857'), 3857)
        #     maskFc_3857_desc = arcpy.da.Describe(maskFc_3857)
        #     maskFc_3857_extent = maskFc_3857_desc['extent']
        #     # Should always be False
        #     #poly1.disjoint(maskFc_3857_extent)
        #     poly1_3857_int = poly1.intersect(maskFc_3857_extent, 4)

        #     clip3 = arcpy.analysis.Clip(maskFc_3857, poly1_3857_int)
        #     clip3_extent = arcpy.da.Describe(clip3)['extent']

        #     x_start = clip3_extent.XMin
        #     x_end = clip3_extent.XMax
        #     y_start = clip3_extent.YMin
        #     y_end = clip3_extent.YMax

        #     log.debug(f'final_x_start: {x_start} and x_end: {x_end}')
        #     log.debug(f'final_y_start: {y_start} and y_end: {y_end}')
        #     ept_extent = str(x_start) + ', ' + str(x_end) + '], [' + str(y_start) + ', ' + str(y_end)

        #     parts.append(['_pt' + str(p), ept_extent])

        # arcpy.env.outputCoordinateSystem =\ srOut

    else:
        splits = 1
        ept_extent = str(geom.extent.XMin) + ', ' + str(geom.extent.XMax) + '], [' + str(geom.extent.YMin) + ', ' + str(geom.extent.YMax)
        parts.append(['_pt1', ept_extent])

    return parts, square_area



def try_to_delete(rasRes, log):
    if arcpy.Exists(rasRes):
        try:
            arcpy.Delete_management(rasRes)
        except arcpy.ExecuteError:
            log.warning('could not remove using arcpy.Delete, trying os.remove')
            os.remove(rasRes)






import os
import re
import requests
from bs4 import BeautifulSoup
from urllib.parse import urljoin, urlparse
from time import sleep
import logging

def download_usgs_laz(
    page_url: str,
    output_dir: str,
    alt_output_dir: str,
    log: logging.Logger,
    timeout: int = 60,
    max_retries: int = 5,
    retry_delay: int = 5,
    user_agent: str = "Mozilla/5.0"
):
    """
    Download LAZ files from a USGS Apache directory listing.
    Skips existing files with matching size.
    Retries failed downloads up to max_retries times.

    Parameters
    ----------
    page_url : str
        Apache directory URL containing LAZ files
    output_dir : str
        Local directory for downloads
    alt_output_dir : str
        Local directory for downloads
    timeout : int
        HTTP timeout (seconds)
    max_retries : int
        Maximum download retry attempts per file
    retry_delay : int
        Seconds to wait between retries
    user_agent : str
        HTTP user-agent
    """
    try:
        os.makedirs(output_dir, exist_ok=True)

        session = requests.Session()
        session.headers.update({"User-Agent": user_agent})

        log.info(f"Scanning: {page_url}")
        resp = session.get(page_url, timeout=timeout)
        resp.raise_for_status()

        soup = BeautifulSoup(resp.text, "html.parser")

        laz_files = []

        for link in soup.find_all("a", href=True):
            href = link["href"]
            if href.lower().endswith(".laz"):
                filename = os.path.basename(urlparse(href).path)

                # Extract file size from Apache listing line
                line_text = link.parent.get_text(" ", strip=True)
                size_match = re.search(
                    rf"{re.escape(filename)}\s+\d+\-\w+\-\d+\s+\d+:\d+\s+(\d+)",
                    line_text
                )

                expected_size = int(size_match.group(1)) if size_match else None
                laz_files.append((filename, urljoin(page_url, href), expected_size))

        log.info(f"Found {len(laz_files)} LAZ files")

        for filename, url, expected_size in laz_files:
            out_path = os.path.join(output_dir, filename)
            alt_out_path = os.path.join(alt_output_dir, filename.lower())

            # Skip valid existing file
            if (os.path.exists(out_path) and expected_size) or (os.path.exists(alt_out_path) and expected_size):
                if os.path.exists(out_path):
                    if os.path.getsize(out_path) == expected_size:
                        log.info(f"Exists (size OK): {filename}")
                        return_path = out_path
                        continue
                elif os.path.exists(alt_out_path):
                    if os.path.getsize(alt_out_path) == expected_size:
                        log.info(f"Exists (size OK): {filename}")
                        return_path = alt_out_path
                        continue
                else:
                    log.info(f"Re-downloading (size mismatch): {filename}") 

            success = False

            for attempt in range(1, max_retries + 1):
                try:
                    log.info(f"Downloading ({attempt}/{max_retries}): {filename}")  

                    with session.get(url, stream=True, timeout=timeout) as r:
                        r.raise_for_status()
                        with open(out_path, "wb") as f:
                            for chunk in r.iter_content(chunk_size=1024 * 1024):
                                if chunk:
                                    f.write(chunk)

                    downloaded_size = os.path.getsize(out_path)

                    if expected_size and downloaded_size != expected_size:
                        raise ValueError(
                            f"Size mismatch (expected {expected_size}, got {downloaded_size})"
                        )

                    log.info(f"Download complete: {filename}")
                    return_path = out_path
                    success = True
                    break

                except Exception as e:
                    log.warning(f"Attempt {attempt} failed: {e}")

                    if os.path.exists(out_path):
                        os.remove(out_path)

                    if attempt < max_retries:
                        sleep(retry_delay)

            if not success:
                log.warning(f"FAILED after {max_retries} attempts: {filename}")

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

    finally:
        return return_path


import json
import os
import pickle
import pdal


def get_laz_bounds_and_crs(laz_file, pkl_path, write_pickle=True):
    """
    Extract X/Y/Z bounds and CRS from a LAZ file.
    Optionally writes results to a pickle file with the same basename.

    Parameters
    ----------
    laz_file : str
        Path to LAS/LAZ file
    pkl_path : str
        Path to pickle file
    write_pickle : bool, optional
        If True, write results to pkl_path

    Returns
    -------
    dict
        Dictionary containing bounds and CRS information
    """

    pipeline_def = {
        "pipeline": [
            {
                "type": "readers.las",
                "filename": laz_file
            }
        ]
    }

    pipeline = pdal.Pipeline(json.dumps(pipeline_def))
    pipeline.execute()

    metadata = pipeline.metadata
    readers_meta = metadata["metadata"]["readers.las"]

    result = {
        "source_file": os.path.abspath(laz_file),
        "bounds": {
            "minx": readers_meta.get("minx"),
            "maxx": readers_meta.get("maxx"),
            "miny": readers_meta.get("miny"),
            "maxy": readers_meta.get("maxy"),
            "minz": readers_meta.get("minz"),
            "maxz": readers_meta.get("maxz")
        },
        "crs": {
            "authority": readers_meta.get("srs", {}).get("authority"),
            "horizontal": readers_meta.get("srs", {}).get("horizontal"),
            "vertical": readers_meta.get("srs", {}).get("vertical"),
            "wkt": readers_meta.get("srs", {}).get("wkt")
        }
    }

    if write_pickle:

        with open(pkl_path, "wb") as f:
            pickle.dump(result, f, protocol=pickle.HIGHEST_PROTOCOL)

        result["pickle_file"] = os.path.abspath(pkl_path)

    return result


##if __name__ == "__main__":
##    laz_path = "input.laz"
##    info = get_laz_bounds_and_crs(laz_path)
##    print(info)



def doLidarDEMs(monthly_wesm_ept_mashup, dem_polygon, 
         pdal_exe, gsds, procDir, snap, breakpolys, breaklines, 
         tElevFile, bareEarthReturnMinFile, firstReturnMaxFile, cntBeFile, cnt1rFile, cntPlsFile,
         int1rMinFile, int1rMaxFile, intBeMaxFile, ept_wesm_project_file, lidar_download_directory, cleanup, messages):
    
    arguments = [monthly_wesm_ept_mashup, dem_polygon, 
        pdal_exe, gsds, procDir, snap, breakpolys, breaklines, 
        tElevFile, bareEarthReturnMinFile, firstReturnMaxFile, cntBeFile, cnt1rFile, cntPlsFile,
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
        huc12, huc8 = df.figureItOut(tElevFile)
        # the DEP huc DEM naming convention
        pattern27 = 'e[cfpxv][0-9]m\\d{10,16}'
        pattern26 = '[0-9]m\\d{10,16}'
        pattern25 = '[0-9]m_d{10,16}'
        pattern22 = '[+_][0-9]m[+_]'

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

        if not os.path.isdir(lidar_download_directory):
            os.makedirs(lidar_download_directory)

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
            print(v)
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
        if not interpFound:# not in fElevFile:
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

        maskRastBase = 'mask_rast_'
        maskFc, maskFc_area, maskFcOut, maskRastOut, hucRastOut, FDSet = prepPolygonBoundary(dem_polygon, log, sgdb, srOut, srSfx, maskRastBase, demLists)
        
##----------------------------------------------------------------------

    # # check for collection change (different priorities) to restrict further data
        fixedFolder = str(arcpy.CreateFolder_management(sfldr, 'fixed'))

        allTilesList = []

        #eptDir = os.path.dirname(os.path.dirname(monthly_wesm_ept_mashup))
        wesm_huc12_all = arcpy.analysis.Clip(monthly_wesm_ept_mashup, maskFcOut, opj('in_memory', 'meets_3dep'))
        wesm_huc12 = arcpy.Select_analysis(wesm_huc12_all, 'wesm_' + huc12, where_clause= "lpc_category = 'Meets'")
        select_ql0_ql1 = arcpy.Select_analysis(wesm_huc12, where_clause = "ql = 'QL 1' OR ql = 'QL 0'")
        count_ql0_ql1 = int(str(arcpy.GetCount_management(select_ql0_ql1)))
        if count_ql0_ql1 >= 1:
            ql1 = True
        else:
            ql1 = False
        # code fails on QL1 data for 071200030402, downloaded LAS for 79951 work id was 143 GB and caused ExtractLas to fail
        # assert count_ql0_ql1 < 1, 'DEM builder not yet configured for QL1 or QL0 density data, email bkgelder@iastate.edu to request upgrade'
        if ql1:
            log.warning('DEM builder not yet configured for QL1 or QL0 density data, email bkgelder@iastate.edu to request upgrade')

##----------------------------------------------------------------------


##----------------------------------------------------------------------
        work_id_name = 'workunit_id'
        build_threshold = 0.99#99 - reduced to 0.99 for Lake Erie at Toledo (04100009)
        if df.testForZero(wesm_huc12):
            prev_merged, merged_area, addOrderField = organizeProjectsByDate(wesm_huc12, work_id_name, maskFc_area, build_threshold, log)

        if df.testForZero(prev_merged):
            get_USGS_LAZ = True # otherwise use PDAL to get data via EPT 
            if get_USGS_LAZ:
                with arcpy.da.SearchCursor(prev_merged, ['SHAPE@', work_id_name, 'lpc_link'] + url_list, sql_clause = [None, 'ORDER BY ' + addOrderField.getInput(1) + ' DESC']) as scur:
                    for srow in scur:
                        log.debug(f'{work_id_name} is: {srow[1]} and lpc_link is: {srow[2]} ')
                
##            geom_srOut_copy = getLidarFiles(wesm_huc12, work_id_name, pdal_exe, prev_merged, addOrderField, log, sgdb, sfldr, srOut, srOutCode, huc12, lidar_download_directory, maskFcOut, fixedFolder, inm, FDSet, allTilesList, procDir, ql1, getEPT)
##
##def getLidarFiles(wesm_huc12, work_id_name, pdal_exe, prev_merged, addOrderField, log, sgdb, sfldr, srOut, srOutCode, huc12, eleDir, maskFcOut, fixedFolder, inm, FDSet, allTilesList, procDir, ql1, getEPT):
##    try:
##            if not getEPT:
##                geom_srOut_copy = None
            
            else:
                eleDir = lidar_download_directory
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
                        parts, square_area = queryParts(geom, geom_extent, maskFcOut, srOut, sgdb, log, ql1)#maskFc_3857, maskFc_3857_desc)

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
                            ept_laz_path = Path(ept_laz_full_filename)
                            # check for a local backup copy, if exists, update path to avoid re-downloading
                            new_path = Path('M:/', *ept_laz_path.parts[1:])
                            alt_ept_laz_full_filename = new_path#ept_laz_path.replace(anchor = 'M:')
                            if os.path.isfile(alt_ept_laz_full_filename):
                                log.debug(f'using alt laz file {alt_ept_laz_full_filename}')
                                ept_laz_full_filename = str(alt_ept_laz_full_filename)
                            # local_ept_laz_full_filename = os.altsep.join([procDir.replace(os.path.sep, os.path.altsep), ept_laz_filename])

                            # pipeline json requires / not \ for path separator
                            ept_las_full_filename = os.altsep.join([procDir.replace(os.path.sep, os.path.altsep), ept_las_filename])
                            # NEEDS UPDATE - check for files in lidar_download_directory as well - getting ready below
                            # dl_ept_laz_full_filename = ept_laz_full_filename.replace(eptDir, eleDir)
                            # dl_ept_zlas_full_filename = ept_zlas_full_filename.replace(eptDir, eleDir)
                            if os.path.isfile(ept_laz_full_filename) and not os.path.isfile(ept_las_full_filename):
                                log.info('converting laz to las from archive')
                                log.info(f"arguments: {ept_laz_full_filename}, {procDir}")
                                las_result = arcpy.conversion.ConvertLas(ept_laz_full_filename, procDir)#, compression = 'ZLAS')
                                log.info(las_result)
                            elif os.path.isfile(ept_zlas_full_filename) and not os.path.isfile(ept_las_full_filename):
                                log.info('converting zlas to las')
                                log.info(f"arguments: {ept_zlas_full_filename}, {procDir}")
                                las_result = arcpy.conversion.ConvertLas(ept_zlas_full_filename, procDir)#, compression = 'ZLAS')
                                log.info(las_result)
                            # if zlas does not exist, get las then convert to zlas
                            if not os.path.isfile(ept_zlas_full_filename) and not os.path.isfile(ept_laz_full_filename):
                                log.info('Getting LAS from EPT')
                                log.info(ept_zlas_filename)

                                ept_json_filename = "_".join(["get", "ept", huc12, str(work_id_part) + ".json"])

                                df.create_needed_dirs_and_gdbs(ept_las_full_filename, log)
                                df.create_needed_dirs_and_gdbs(eleDir, log)
                                ept_json_full_filename = create_ept_json_pipeline(ept_json_filename, eleDir, ept_las_full_filename, extent_request, ept_address, srOutCode)
                                df.create_needed_dirs_and_gdbs(ept_json_full_filename, log)

                                if not os.path.exists(ept_las_full_filename):
                                    run_string = " ".join([pdal_exe, "pipeline", ept_json_full_filename])
                                    # estimate download time based on 102500040309 (area 1175 km2) in 4 parts
                                    m2_per_sec = 1175.2*1000**2/len(parts)/2200
                                    log.debug(f'pdal run_string: {run_string}')
                                    log.info(f'Estimated pdal download time (for QL2 lidar): {round(square_area/(m2_per_sec * len(parts) * 60), 2)} minutes for {ept_json_filename}')
                                    co = subprocess.call(run_string, creationflags=CREATE_NO_WINDOW)
                                    # co = subprocess.run(run_string)
                                    log.debug(f'completed pdal run_string')

                                # archive as zlas for use later in this script and re-use
                                stats = os.stat(ept_las_full_filename)
                                if stats.st_size > las_size_threshold:
                                    if not os.path.isfile(ept_laz_full_filename) and not os.path.isfile(ept_zlas_full_filename):
                                        # takes too long, laz is much faster to create
                                        # log.info('converting las to zlas for archive')
                                        # zlas_result = arcpy.conversion.ConvertLas(ept_las_full_filename, ele, compression = 'ZLAS', las_options = None)
                                        # log.info(zlas_result)

                                        log.debug('converting las to laz for archive')
                                        laz_json_filename = "_".join(["laz", "ept", huc12, str(work_id_part) + ".json"])
                                        laz_json_full_filename = create_laz_json_pipeline(laz_json_filename, eleDir, ept_las_full_filename, ept_laz_full_filename)
                                        laz_run_string = " ".join([pdal_exe, "pipeline", laz_json_full_filename])
                                        log.debug(f'pdal run_string: {laz_run_string}')
                                        co = subprocess.call(laz_run_string, creationflags=CREATE_NO_WINDOW)
                                        # co = subprocess.run(laz_run_string)
                                        log.debug(f'ran pdal run_string')

                                        # ADD CODE to do local then copy, SLOW on network drives 
                                        # laz_result = arcpy.conversion.ConvertLas(ept_las_full_filename, eleDir, compression = 'LAZ', las_options = None)
                                        # log.debug(laz_result)

                                    # arcpy.Delete_management(ept_las_full_filename)
                                else:
                                    poly35 = geom_extent.polygon
                                    p35 = arcpy.management.CopyFeatures(poly35, opj(sgdb, 'failed_' + valid_geom_name))
                                    log.warning(f"{ept_las_full_filename} has very small file size; plotting extent as poly: {p35}")
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

                        # # if co.returncode == 0:
                        # if os.path.exists(ept_las_full_filename) and stats.st_size > las_size_threshold:
                        #     cl2Las = processEptLas(sgdb, sfldr, srOutCode, fixedFolder, geom_srOut, ept_las_full_filename, srOut, inm, FDSet, procDir, allTilesList, log, time, work_id)
                        #     #remove invalid geometry
                        #     if cl2Las is None:
                        #         log.warning('deleting ' + str(geom_srOut_copy))
                        #         delete = arcpy.Delete_management(geom_srOut_copy)
                        #         # remove from dates
                        #         prev_merged = arcpy.Select_analysis(prev_merged, where_clause = work_id_name + ' <> ' + str(work_id))
                        #     else:
                        #         cl2LasNotNone = cl2Las
                        # elif os.path.exists(ept_las_full_filename):
                        #     log.warning('EPT LAS file too small, no points')
                        #     run_string += ' --debug'
                        #     log.warning(f'Invalid call was {run_string}')

                        # else:
                        #     log.warning('No valid output from ept.json request')
                        #     log.warning(f'Invalid call was {run_string}')
                        #     sys.exit(1)

        # #revert to previous good cl2Las if last is None
        # if cl2Las is None:
        #     cl2Las = cl2LasNotNone

        # return cl2Las, geom_srOut_copy
##        return geom_srOut_copy
##
##    except Exception as e:
##        print('handling as exception')
####        log.debug(e.message)
##        if sys.version_info.major == 2:
##            arcpy.AddError(e.message)
##            print(e.message)
##        elif sys.version_info.major == 3:
##            arcpy.AddError(e)
##            print(e)
##
##        tb = sys.exc_info()[2]
##        tbinfo = traceback.format_tb(tb)[0]
##
##        # Concatenate information together concerning the error into a message string
##        pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
##        # Return python error messages for use in script tool or Python Window
##        arcpy.AddError(pymsg)
##        # Print Python error messages for use in Python / Python Window
##        print(pymsg + "\n")
##
##        if arcpy.GetMessages(2) not in pymsg:
##            msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
##            arcpy.AddError(msgs)
##            print(msgs)
##
##    except:
##        print('handling as except')
##        # Get the traceback object
##        tb = sys.exc_info()[2]
##        tbinfo = traceback.format_tb(tb)[0]
##
##        # Concatenate information together concerning the error into a message string
##        pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
##        # Return python error messages for use in script tool or Python Window
##        arcpy.AddError(pymsg)
##        # Print Python error messages for use in Python / Python Window
##        print(pymsg + "\n")
##
##        if arcpy.GetMessages(2) not in pymsg:
##            msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
##            arcpy.AddError(msgs)
##            print(msgs)


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

        # indicate there are no datasets to process further
        else:
            prev_merged = None
            merged_area = 0.00


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

    return tElevFile



import arcpy
import pickle
import os
import glob

def build_bounds_of_pkls(pkl_dir, out_gdb, out_fc_name, out_fc, work_id_name, work_id):
    '''
    # --------------------------------
    # Inputs
    # --------------------------------
    pkl_dir = r"E:\\DEP\\USGS_LPC\\IA_NorthCentral_2020_D20\\IA_NorthCentral_2_2020\\pkl"
    out_gdb = r"E:\\DEP\\USGS_LPC\\IA_NorthCentral_2020_D20\\IA_NorthCentral_2_2020\\bounds\\laz_bounds.gdb"
    out_fc_name = "laz_bounds"
    out_fc = os.path.join(out_gdb, out_fc_name)
    '''
    mem_fc = "in_memory\\laz_bounds_tmp"

    # --------------------------------
    # Create file geodatabase if needed
    # --------------------------------
    arcpy.env.overwriteOutput = True
    gdb_dir, gdb_name = os.path.split(out_gdb)

    if not arcpy.Exists(out_gdb):
        arcpy.management.CreateFileGDB(gdb_dir, gdb_name)

    # --------------------------------
    # Create in-memory feature class
    # --------------------------------
    if arcpy.Exists(mem_fc):
        arcpy.management.Delete(mem_fc)

    spatial_ref = None  # will be set from first valid PKL

    arcpy.management.CreateFeatureclass(
        out_path="in_memory",
        out_name="laz_bounds_tmp",
        geometry_type="POLYGON"
    )

    arcpy.management.AddField(
        mem_fc,
        "laz_file",
        "TEXT",
        field_length=200
    )

    arcpy.management.AddField(
        mem_fc,
        work_id_name,
        "LONG"
    )

    # --------------------------------
    # Loop over PKL files
    # --------------------------------
    pkl_files = glob.glob(os.path.join(pkl_dir, "*.pkl"))

    inserted = 0
    skipped = 0

    with arcpy.da.InsertCursor(mem_fc, ["SHAPE@", "laz_file", work_id_name]) as cursor:
        for pkl_path in pkl_files:
            try:
                with open(pkl_path, "rb") as f:
                    data = pickle.load(f)
                print(pkl_path)

                bounds = data["bounds"]
                laz_path = data["source_file"]
                wkt = data["crs"]["wkt"]

                minx = bounds["minx"]
                miny = bounds["miny"]
                maxx = bounds["maxx"]
                maxy = bounds["maxy"]

                # Initialize spatial reference from first PKL
                if spatial_ref is None:
                    spatial_ref = arcpy.SpatialReference()
                    spatial_ref.loadFromString(wkt)
                    arcpy.management.DefineProjection(mem_fc, spatial_ref)

                # Build polygon
                array = arcpy.Array([
                    arcpy.Point(minx, miny),
                    arcpy.Point(maxx, miny),
                    arcpy.Point(maxx, maxy),
                    arcpy.Point(minx, maxy),
                    arcpy.Point(minx, miny)
                ])

                polygon = arcpy.Polygon(array, spatial_ref)

                cursor.insertRow([polygon, laz_path, work_id])
                inserted += 1

            except Exception as e:
                skipped += 1
                arcpy.AddWarning(f"Skipped {os.path.basename(pkl_path)}: {e}")

    # --------------------------------
    # Copy to file GDB
    # --------------------------------
    if arcpy.Exists(out_fc):
        arcpy.management.Delete(out_fc)

    arcpy.management.CopyFeatures(mem_fc, out_fc)

    # --------------------------------
    # Cleanup
    # --------------------------------
    arcpy.management.Delete(mem_fc)

    print(f"Done. Inserted {inserted} features. Skipped {skipped} files.")

##----------------------------------------------------------------------



if __name__ == "__main__":
    if len(sys.argv) == 1:
        #Paste arguments into here for use within Python Window
        arcpy.AddMessage("Whoo, hoo! Running from Python Window!")
        # cleanup = False
        parameters = 	["C:/Program Files/ArcGIS/Pro/bin/Python/envs/arcgispro-py3/pythonw.exe",
    "C:/Temp/cmd_build_terrain_derivatives_main.pyt",
    "E:/DEP/Elevation_databases/ept.gdb/ept_resources_2026_01_01",
    "M:/DEP/Man_Data_ACPF/dep_ACPF2023/07080105/idepACPF0708010509010101.gdb/buf0708010509010101",
    "C:/Users/bkgelder/.conda/envs/pdal/Library/bin/pdal.exe",
    "E:/DEP_Proc/DEMProc/LAS_dem2013_1m_0708010509010101",
    "get_USGS_LAZ"#alternative - "get_PDAL_LAZ",
    "E:/DEP_Checkout/Man_Data_ACPF/dep_ACPF2023/07080105/idepACPF0708010509010101.gdb/wesm_ept_resources_2026_01_01_0708010509010101",
    "E:/DEP_Checkout/Man_Data_ACPF/dep_ACPF2023/07080105/idepACPF0708010509010101.gdb/wesm_tiles_2026_01_23_0708010509010101",
    "E:/DEP_Checkout/USGS_LPC",
    "False"]
        for i in parameters[2:]:
            sys.argv.append(i)

    else:
        #For use via Windows Command Line
        #above 'parameters' come in via command line arguments, nothing else needed
        arcpy.AddMessage("Whoo, hoo! Command-line enabled!")
        # clean up the folder after done processing
        # cleanup = True

    # inputs then outputs, change "" to Python None
    (monthly_wesm_ept_mashup, dem_polygon, 
         pdal_exe, procDir, get_lidar_method,
         ept_wesm_project_file, wesm_huc12_tiles, lidar_download_directory, cleanup
        ) = [i if i != "" else None for i in sys.argv[1:]]

    # switch a text 'True' into a real Python True
    cleanup = True if cleanup == "True" else False

    messages = msgStub()

    arguments = [monthly_wesm_ept_mashup, dem_polygon, pdal_exe, procDir, get_lidar_method, ept_wesm_project_file, wesm_huc12_tiles, lidar_download_directory, cleanup]
    for a in arguments:
        if a == arguments[0]:
            arg_str = str(a) + '\n'
        else:
            arg_str += str(a) + '\n'
    messages.addMessage('Tool: Executing with parameters:\n' + arg_str)
    arcpy.env.overwriteOutput = True
    arcpy.CheckOutExtension('Spatial')
    arcpy.CheckOutExtension('3D')
    arcpy.env.ZResolution = '0.01'
    try:
        huc12, huc8 = df.figureItOut(dem_polygon)#tElevFile)
        pattern27 = 'e[cfpxv][0-9]m\\d{10,16}'
        pattern26 = '[0-9]m\\d{10,16}'
        pattern25 = '[0-9]m_d{10,16}'
        pattern22 = '[+_][0-9]m[+_]'
        if procDir is not None:
            if not os.path.isdir(procDir):
                os.makedirs(procDir)
            arcpy.env.scratchWorkspace = procDir
            sfldr = arcpy.env.scratchFolder
        else:
            sfldr = arcpy.env.scratchFolder
            procDir = sfldr
        if lidar_download_directory is None:
            lidar_download_directory = procDir
        if not os.path.isdir(lidar_download_directory):
            os.makedirs(lidar_download_directory)
        sgdb = arcpy.env.scratchGDB
        arcpy.env.scratchWorkspace = sgdb
        arcpy.env.workspace = sgdb
        if procDir is None:
            procDir = sfldr
        if lidar_download_directory is None:
            lidar_download_directory = sfldr
        node = platform.node()
        logProc = df.defineLocalProc(node)
        if not os.path.isdir(logProc):
            logProc = sfldr
        if cleanup:
            log, nowYmd, logName, startTime = df.setupLoggingNoCh(logProc, sys.argv[0], huc12)
            arcpy.SetLogHistory = False
        else:
            log, nowYmd, logName, startTime = df.setupLoggingNew(logProc, sys.argv[0], huc12)
            arcpy.SetLogHistory = True
        log.info('Beginning execution: ' + time.asctime())
        log.debug('sys.argv is: ' + str(sys.argv) + '\n')
        log.info('Processing HUC: ' + huc12)
        log.info(f'procDir: {procDir}')
        log.info(f'lidar_download_directory: {lidar_download_directory}')
        fElevDesc = arcpy.da.Describe(dem_polygon)
        srOut = fElevDesc['spatialReference']
        srOutCode = srOut.PCSCode
        assert srOutCode < 32768, 'EPSG spatial reference code too large, PDAL will not recognize'
        log.info('Output will be in EPSG Code (spatial reference): ' + str(srOutCode))
        log.info('Log file at ' + logName)
        messages.addMessage('Log file at ' + logName)
        flib_metadata_template, derivative_metadata = df.getMetadata(['flib', 'deriv'], procDir, log)
        inm = 'in_memory'
        srOut = arcpy.SpatialReference(int(srOutCode), 5703)
        srOutNoVCS = arcpy.SpatialReference(int(srOutCode))
        arcpy.env.outputCoordinateSystem = srOut
        srSfx = '_' + str(srOutCode)
        # maskRastBase = 'mask_rast_'
        # maskFc, maskFc_area, maskFcOut, maskRastOut, hucRastOut, FDSet = prepPolygonBoundary(dem_polygon, log, sgdb, srOut, srSfx, maskRastBase, demLists)

        assert int(arcpy.GetCount_management(dem_polygon).getOutput(0)) < 2, 'multiple features in polygon feature class'
        assert int(arcpy.GetCount_management(dem_polygon).getOutput(0)) > 0, 'no features in polygon feature class'
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

        # ## Set up geodatabase to store the multipoint files and terrains (necessary all inputs be in feature dataset
        # # Vertical units are in meters (float) so use a meter-based reference
        # FDSet = arcpy.CreateFeatureDataset_management(sgdb, "Lidar_pts", srOut)
        maskFcOut = projIfNeeded(maskFc, os.path.join(str(sgdb), 'get_laz_bnds' + srSfx), srOut)
        log.info("maskFcOut complete at: " + time.asctime())

        fixedFolder = str(arcpy.CreateFolder_management(sfldr, 'fixed'))
        allTilesList = []
        wesm_huc12_all = arcpy.analysis.Clip(monthly_wesm_ept_mashup, maskFcOut, opj('in_memory', 'meets_3dep'))
        wesm_huc12 = arcpy.Select_analysis(wesm_huc12_all, 'wesm_' + huc12, where_clause="lpc_category = 'Meets'")
        select_ql0_ql1 = arcpy.Select_analysis(wesm_huc12, where_clause="ql = 'QL 1' OR ql = 'QL 0'")
        count_ql0_ql1 = int(str(arcpy.GetCount_management(select_ql0_ql1)))
        if count_ql0_ql1 >= 1:
            ql1 = True
        else:
            ql1 = False
        if ql1:
            log.warning('DEM builder not yet configured for QL1 or QL0 density data, email bkgelder@iastate.edu to request upgrade')
        work_id_name = 'workunit_id'
        build_threshold = 0.99
        if df.testForZero(wesm_huc12):
            prev_merged, merged_area, addOrderField = organizeProjectsByDate(wesm_huc12, work_id_name, maskFc_area, build_threshold, log)
        if df.testForZero(prev_merged):
            if get_lidar_method == "get_USGS_LAZ":
            # if get_USGS_LAZ:
                bounds_list = []
                with arcpy.da.SearchCursor(prev_merged, ['SHAPE@', work_id_name, 'lpc_link'], sql_clause=[None, 'ORDER BY ' + addOrderField.getInput(1) + ' DESC']) as scur:
                    for srow in scur:
                        log.debug(f'{work_id_name} is: {srow[1]} and lpc_link is: {srow[2]} ')
                        # dl_dir = os.path.join('E:\\DEP\\USGS_LPC', os.path.basename(os.path.dirname(srow[2])), os.path.basename(srow[2]), 'LAZ')
                        dl_dir = os.path.join(lidar_download_directory, os.path.basename(os.path.dirname(srow[2])), os.path.basename(srow[2]), 'LAZ')
                        alt_dl_dir = opj('M:/DEP/USGS_LPC', os.path.basename(os.path.dirname(srow[2])), os.path.basename(srow[2]), 'LAZ')
                        page_url = srow[2] + '/LAZ/'
                        try:
                            return_path = download_usgs_laz(page_url = page_url, output_dir = dl_dir, alt_output_dir = alt_dl_dir, log = log)
                        except:
                            log.warning(f'failed download {page_url}')
                            return_path = None

                        if return_path is not None:
                            # create pkl dir and get bounds from each laz file
                            pkl_dir = os.path.dirname(return_path).replace('LAZ', 'pkl')
                            log.debug(f'pkl_dir is: {pkl_dir}')
                            df.create_needed_dirs_and_gdbs(pkl_dir, log)
                            # pkl_dir = dl_dir.replace('LAZ', 'pkl')

                            return_path_po = Path(return_path)
                            laz_dir = return_path_po.parent
                            cwd = os.getcwd()
                            os.chdir(laz_dir)
                            lazs = glob.glob('*.laz')
                            for laz_file in lazs:
                                write_pickle = True
                                basename = os.path.basename(laz_file)
                                base, _ = os.path.splitext(basename)
                                pkl_file = base + ".pkl"
                                pkl_path = os.path.join(pkl_dir, pkl_file)
                                if not os.path.isfile(pkl_path):
                                    log.info(f'writing pickle file: {pkl_path}')
                                    get_laz_bounds_and_crs(laz_file, pkl_path, write_pickle=write_pickle)
                            
                            bounds_dir = pkl_dir.replace('pkl', 'bounds')
                            df.create_needed_dirs_and_gdbs(bounds_dir, log)
                            out_gdb = opj(bounds_dir, 'laz_bounds_' + str(srow[1]) + '.gdb')
                            out_fc_name = 'laz_bounds_' + str(srow[1])
                            out_fc = os.path.join(out_gdb, out_fc_name)
                            df.create_needed_dirs_and_gdbs(out_fc, log)
                            if not arcpy.Exists(out_fc):
                                build_bounds_of_pkls(pkl_dir, out_gdb, out_fc_name, out_fc, work_id_name, srow[1])
                            bounds_list.append(out_fc)

                out_fc_merge = arcpy.Merge_management(bounds_list, os.path.join('in_memory', 'out_fc_merge'))

                out_fc_clip = arcpy.Clip_analysis(out_fc_merge, maskFc)

                po_procDir = Path(procDir)
                po_all_laz_dir = po_procDir.joinpath('laz_all_points')
                if not os.path.isdir(po_all_laz_dir):
                    os.makedirs(po_all_laz_dir)

                with arcpy.da.SearchCursor(out_fc_clip, ['OBJECTID', 'laz_file']) as scur:
                    for srow in scur:
                        laz_file = srow[1].replace('\\laz\\', '\\USGS_LPC\\')
                        if not os.path.isfile(laz_file):
                            laz_file = laz_file.replace('E:\\DEP\\USGS_LPC', 'M:\\DEP\\USGS_LPC')
                        po_laz_file = Path(laz_file)
                        dest = po_all_laz_dir.joinpath(po_laz_file.name)
                        if not os.path.isfile(dest):
                            shutil.copy(po_laz_file, dest)

                df.create_needed_dirs_and_gdbs(wesm_huc12_tiles, log)
                wesm_huc12_tiles = arcpy.CopyFeatures_management(out_fc_clip, wesm_huc12_tiles)
                                

##                        if not os.path.isfile(ept_laz_full_filename) and (not os.path.isfile(ept_zlas_full_filename)):
##                            log.debug('converting las to laz for archive')
##                            laz_json_filename = '_'.join(['laz', 'ept', huc12, str(work_id_part) + '.json'])
##                            laz_json_full_filename = create_laz_json_pipeline(laz_json_filename, eleDir, ept_las_full_filename, ept_laz_full_filename)
##                            laz_run_string = ' '.join([pdal_exe, 'pipeline', laz_json_full_filename])
##                            log.debug(f'pdal run_string: {laz_run_string}')
##                            co = subprocess.call(laz_run_string, creationflags=CREATE_NO_WINDOW)
##                            log.debug(f'ran pdal run_string')



                        
            elif get_lidar_method == "get_PDAL_LAZ":
                eleDir = lidar_download_directory
                prev_merged_projected_3857 = arcpy.management.Project(prev_merged, 'proj_trial', 3857)
                url_list = df.getfields(wesm_huc12, 'url*')
                with arcpy.da.SearchCursor(prev_merged_projected_3857, ['SHAPE@', work_id_name] + url_list, sql_clause=[None, 'ORDER BY ' + addOrderField.getInput(1) + ' DESC']) as scur:
                    for srow in scur:
                        log.debug(f'{work_id_name} is: {srow[1]}')
                        geom = srow[0]
                        geom_extent = geom.extent
                        las_size_threshold = 750
                        parts, square_area = queryParts(geom, geom_extent, maskFcOut, srOut, sgdb, log, ql1)
                        arcpy.env.outputCoordinateSystem = srOut
                        for part in parts:
                            work_id = srow[1]
                            work_id_part = str(srow[1]) + part[0]
                            extent_request = part[1]
                            geom_srOut = geom.projectAs(arcpy.SpatialReference(srOutCode))
                            valid_geom_name = arcpy.ValidateTableName('geom_proj_' + work_id_part, sgdb)
                            geom_srOut_copy = arcpy.CopyFeatures_management(geom_srOut, valid_geom_name)
                            arcpy.management.AddField(geom_srOut_copy, work_id_name, "LONG")
                            arcpy.management.CalculateField(geom_srOut_copy, work_id_name, work_id, "PYTHON")
                            urls = srow[1:]
                            for u in urls:
                                if u is not None:
                                    url = u
                            ept_address = url
                            ept_zlas_filename = '_'.join(['ept', huc12, str(work_id_part) + '.zlas'])
                            ept_zlas_full_filename = os.altsep.join([eleDir.replace(os.path.sep, os.path.altsep), ept_zlas_filename])
                            ept_las_filename = ept_zlas_filename.replace('.zlas', '.las')
                            ept_laz_filename = ept_zlas_filename.replace('.zlas', '.laz')
                            ept_laz_full_filename = os.altsep.join([eleDir.replace(os.path.sep, os.path.altsep), ept_laz_filename])
                            ept_laz_path = Path(ept_laz_full_filename)
                            new_path = Path('M:/', *ept_laz_path.parts[1:])
                            alt_ept_laz_full_filename = new_path
                            if os.path.isfile(alt_ept_laz_full_filename):
                                log.debug(f'using alt laz file {alt_ept_laz_full_filename}')
                                ept_laz_full_filename = str(alt_ept_laz_full_filename)
                            ept_las_full_filename = os.altsep.join([procDir.replace(os.path.sep, os.path.altsep), ept_las_filename])
                            if os.path.isfile(ept_laz_full_filename) and (not os.path.isfile(ept_las_full_filename)):
                                log.info('converting laz to las from archive')
                                log.info(f'arguments: {ept_laz_full_filename}, {procDir}')
                                las_result = arcpy.conversion.ConvertLas(ept_laz_full_filename, procDir)
                                log.info(las_result)
                            elif os.path.isfile(ept_zlas_full_filename) and (not os.path.isfile(ept_las_full_filename)):
                                log.info('converting zlas to las')
                                log.info(f'arguments: {ept_zlas_full_filename}, {procDir}')
                                las_result = arcpy.conversion.ConvertLas(ept_zlas_full_filename, procDir)
                                log.info(las_result)
                            if not os.path.isfile(ept_zlas_full_filename) and (not os.path.isfile(ept_laz_full_filename)):
                                log.info('Getting LAS from EPT')
                                log.info(ept_zlas_filename)
                                ept_json_filename = '_'.join(['get', 'ept', huc12, str(work_id_part) + '.json'])
                                df.create_needed_dirs_and_gdbs(ept_las_full_filename, log)
                                df.create_needed_dirs_and_gdbs(eleDir, log)
                                ept_json_full_filename = create_ept_json_pipeline(ept_json_filename, eleDir, ept_las_full_filename, extent_request, ept_address, srOutCode)
                                df.create_needed_dirs_and_gdbs(ept_json_full_filename, log)
                                if not os.path.exists(ept_las_full_filename):
                                    run_string = ' '.join([pdal_exe, 'pipeline', ept_json_full_filename])
                                    m2_per_sec = 1175.2 * 1000 ** 2 / len(parts) / 2200
                                    log.debug(f'pdal run_string: {run_string}')
                                    log.info(f'Estimated pdal download time (for QL2 lidar): {round(square_area / (m2_per_sec * len(parts) * 60), 2)} minutes for {ept_json_filename}')
                                    co = subprocess.call(run_string, creationflags=CREATE_NO_WINDOW)
                                    log.debug(f'completed pdal run_string')
                                stats = os.stat(ept_las_full_filename)
                                if stats.st_size > las_size_threshold:
                                    if not os.path.isfile(ept_laz_full_filename) and (not os.path.isfile(ept_zlas_full_filename)):
                                        log.debug('converting las to laz for archive')
                                        laz_json_filename = '_'.join(['laz', 'ept', huc12, str(work_id_part) + '.json'])
                                        laz_json_full_filename = create_laz_json_pipeline(laz_json_filename, eleDir, ept_las_full_filename, ept_laz_full_filename)
                                        laz_run_string = ' '.join([pdal_exe, 'pipeline', laz_json_full_filename])
                                        log.debug(f'pdal run_string: {laz_run_string}')
                                        co = subprocess.call(laz_run_string, creationflags=CREATE_NO_WINDOW)
                                        log.debug(f'ran pdal run_string')
                                else:
                                    poly35 = geom_extent.polygon
                                    p35 = arcpy.management.CopyFeatures(poly35, opj(sgdb, 'failed_' + valid_geom_name))
                                    log.warning(f'{ept_las_full_filename} has very small file size; plotting extent as poly: {p35}')
                            elif os.path.isfile(ept_las_full_filename):
                                stats = os.stat(ept_las_full_filename)
                            else:
                                log.warning("can't get good LAS data, skipping to next project")
                                continue

                log.info('processing EPT lidar feature classes into unioned polygon feature class')
                ept_lidar_fcs = arcpy.ListFeatureClasses(os.path.basename(geom_srOut_copy.getOutput(0))[:10] + '*')
                wesm_huc12_tiles = arcpy.Union_analysis(ept_lidar_fcs)

            df.joinDict(wesm_huc12_tiles, work_id_name, wesm_huc12, work_id_name, ['collect_start', 'collect_end'])

            
    except AssertionError:
        log.warning('assertion failure on: ' + huc12)
        sys.exit(1)
    except:
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = 'PYTHON ERRORS:\nTraceback info:\n' + tbinfo + '\nError Info:\n' + str(sys.exc_info()[1])
        msgs = 'ArcPy ERRORS:\n' + arcpy.GetMessages(2) + '\n'
        log.warning(pymsg)
        log.warning(msgs)
        log.warning('failure on: ' + huc12)
        sys.exit(1)
    finally:
        log.warning('Finished at ' + time.asctime())
        handlers = log.handlers
        for h in handlers:
            log.info('shutting it down!')
            log.removeHandler(h)
            h.close()

#     doLidarDEMs(monthly_wesm_ept_mashup, dem_polygon, 
#          pdal_exe, gsds, procDir, snap, breakpolys, breaklines, 
#          fElevFile, bareEarthReturnMinFile, firstReturnMaxFile, cntBeFile, cnt1rFile, cntPlsFile,
#          int1rMinFile, int1rMaxFile, intBeMaxFile, ept_wesm_project_file, lidar_download_directory, cleanup, messages)#msgStub())

#     arcpy.AddMessage("Back from doEPT!")
