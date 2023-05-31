# -*- coding: utf-8 -*-

import arcpy
from arcpy import env
from arcpy.sa import *
import arcpy.metadata as md
import os
from os.path import join as opj
import traceback
import time
from subprocess import call
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
        self.label = "Build DEM from Lidar via EPT"
        self.description = ""
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
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
        return

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return


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


def processEptLas(sgdb, softwareDir, sfldr, srOutCode, fixedFolder, geom, ept_las, srOut, inm, FDSet, procDir, allTilesList, log, time, work_id):
    '''process a cursor row of data by creating a suitable las file from the input las/laz/zlas dataset
    This inlcudes project and clipping las data files into output dataset and also creating a multipoint
    file from the las data if there is any within the extent'''
    try:
        ept_las_base = os.path.splitext(os.path.basename(ept_las))[0]
        sfx = arcpy.ValidateTableName('_' + ept_las_base, sgdb)
        log.debug('lidar file suffix is: ' + sfx + ' at ' + time.asctime())

        # put in las or zlas format


        # fixedLasPath = ept_las#lasCopy
    ##                                        allTilesList.append(fixedLasPath)
        allLasd = arcpy.CreateLasDataset_management(ept_las, opj(sfldr, 'all' + sfx))
        # allLasdDescDa = arcpy.da.Describe(allLasd)

        # extract to tile geometry and project if necessary
        nameSfx = '_' + str(srOutCode)
        fixedLasBasename = os.path.basename(ept_las)[:-4] + nameSfx + '.las'
    ##            log.debug('fixedLasBasename: ' + fixedLasBasename)
        # some old 3DEP projects don't alway have boundaries and data lining up...
        if work_id < 0:
            tileGeomBuffer5 = geom.buffer(5) #tile Geometry column
            log.debug('--- Use ExtractLas, boundary option, at ' + time.asctime())
            fixedLasd = arcpy.ExtractLas_3d(allLasd, fixedFolder, name_suffix = nameSfx, out_las_dataset = opj(fixedFolder, 'fixed' + sfx + '.lasd'), boundary = tileGeomBuffer5)
        else:
            log.debug('--- Use ExtractLas, no boundary option, at ' + time.asctime())
            fixedLasd = arcpy.ExtractLas_3d(allLasd, fixedFolder, name_suffix = nameSfx, out_las_dataset = opj(fixedFolder, 'fixed' + sfx + '.lasd'))#, boundary = tileGeomBuffer5)
        log.debug(fixedLasd.getMessages())
        fixedLasdDescDa = arcpy.da.Describe(fixedLasd)
        fixedLasPath = opj(fixedLasdDescDa['path'], fixedLasBasename)

        log.debug('--- Done creating LAS dataset and extracting LAS at ' + time.asctime())

        if fixedLasdDescDa['pointCount'] > 0:
            allTilesList.append(fixedLasPath)

        if fixedLasPath in allTilesList:#non 0 amount of lidar points in las
            ## generate multipoint feature class from lidar class 2 (BE) points
            ptOut, cl2Las = genClass2AndMultiPoints(fixedLasPath, sfx, srOut, inm, FDSet, procDir, softwareDir, log)
        # ready so there is something to return
        if 'cl2Las' not in locals():
            cl2Las = None
        if 'ptOut' not in locals():
            ptOut = None
        log.warning('ptOut: ' + str(ptOut) + ' and cl2Las ' + str(cl2Las))

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

def prepPolygonBoundary(huc12fc, huc12, log, sgdb, srOut, srSfx, maskRastBase, demLists):

    try:
        if int(arcpy.GetCount_management(huc12fc).getOutput(0)) > 1:
                # Buffer the HUC12 watershed to the appropriate amount
            with arcpy.da.SearchCursor(huc12fc, ['SHAPE@', 'HUC12'], '"HUC12" = \'' + huc12 + "'") as scur:
                for srow in scur:
                    geom = srow[0]
                    log.info("got geometry at: " + time.asctime())
                    buf = srow[0].buffer(1000)
                    log.info("did buffer at: " + time.asctime())
                    maskFc = arcpy.CopyFeatures_management(buf)

                    geom_copy = arcpy.management.CopyFeatures(geom, opj(sgdb, 'huc' + huc12))

        else:
            maskFc = arcpy.Buffer_analysis(huc12fc, buffer_distance_or_field = '1000 METERS')
            log.info("did buffer at: " + time.asctime())
            # maskFc = arcpy.CopyFeatures_management(buf)
            # log.warning("copied buffer at: " + time.asctime())

            geom_copy = arcpy.management.CopyFeatures(huc12fc, opj(sgdb, 'huc' + huc12))

        arcpy.AddField_management(maskFc, 'id', 'LONG')
        arcpy.CalculateField_management(maskFc, 'id', 1, 'PYTHON')

        log.warning("maskFcPrelim complete at: " + time.asctime())

        ## Set up geodatabase to store the multipoint files and terrains (necessary all inputs be in feature dataset
        # Vertical units are in meters (float) so use a meter-based reference
        FDSet = arcpy.CreateFeatureDataset_management(sgdb, "Lidar_pts", srOut)
        maskFcOut = projIfNeeded(maskFc, os.path.join(str(FDSet), 'buf_huc' + srSfx), srOut)
        log.warning("maskFcOut complete at: " + time.asctime())

        for demList in demLists:
            maskRastOut = arcpy.PolygonToRaster_conversion(maskFcOut, 'id', opj(sgdb, maskRastBase + str(demList[0])), cellsize = demList[0])
            huc_rast_out = arcpy.conversion.PolygonToRaster(geom_copy, 'OBJECTID', opj(sgdb, 'huc_rast' + str(demList[0])), cellsize = demList[0])

        return maskFc, FDSet, maskFcOut, maskRastOut, huc_rast_out

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

def createCountsFromMultipoints(sgdb, maskRastBase, demList, huc12, finalMPinm, finalMP, log, cntFile):#locDict):
    try:
        maskRastOut = opj(sgdb, maskRastBase + str(demList[0]))

        log.warning('---Counting Bare Earth Returns for ' + str(demList[0]) + ' at ' + time.asctime())
        terrCountName = arcpy.ValidateTableName("cnt" + str(demList[0]) + "m_fl_" + huc12, sgdb)
        terrCount = arcpy.PointToRaster_conversion(finalMPinm, arcpy.Describe(finalMP).OIDFieldName, os.path.join(sgdb, terrCountName), "COUNT", "NONE", str(demList[0]))
        cntFileRasterObj = clipCountRaster(terrCount, maskRastOut, cntFile)

        return cntFileRasterObj

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



def createRastersFromTerrains(log, demList, procDir, terrains, huc12):
    try:
        log.debug('snapRaster for Terrain to raster: ' + arcpy.env.snapRaster)
        interpTechnique = 'NATURAL_NEIGHBORS'
        pyramidLevel = '4'
        dem_cellSize = demList[0]
        arcpy.env.cellSize = dem_cellSize
        for terrain in terrains:
            log.warning('---Creating Raster from Terrain for ' + str(demList[0]) + ' using ' + terrain.getInput(6) + ' at ' + time.asctime())
            tempTerrName = generateTempTerrName(procDir, terrain.getInput(6), dem_cellSize, huc12)
            pfFileTemp = os.path.join(procDir, '_'.join(['tmp_ter', terrain.getInput(6), str(demList[0]) + 'm', huc12, 'out.tif']))
            demOut = arcpy.TerrainToRaster_3d(terrain, pfFileTemp, "FLOAT", interpTechnique, "CELLSIZE " + str(demList[0]), pyramidLevel)
####            demList.append(demOut.getOutput(0))

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
##        log.warning('surfcons are ' + surfcons)
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


def clipCountRaster(cnt3mFull, maskRastUTM, cntFileName):

####    maskRastUTM = arcpy.PolygonToRaster_conversion(mask, 'id', cellsize = cellSize)
    cntUTMpre = Con(IsNull(cnt3mFull), 0, cnt3mFull)
    cntUTMter = Con(maskRastUTM, cntUTMpre)
    cntUTMter.save(cntFileName)
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
    ##        log.warning('projecting ' + str(input2))
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

def genClass2AndMultiPoints(allLAZ, sfx, srOut, inm, FDSet, procDir, softwareDir, log):
    try:
        if allLAZ.endswith('.laz') or allLAZ.endswith('.las'):
            """Filters LAS points to class 2 and creates multipoints in FDSet"""
            lasBase = os.path.splitext(os.path.basename(allLAZ))[0]
            cl2LAS = os.path.join(procDir, arcpy.ValidateTableName(lasBase + '_cl2', procDir) + '.las')
        ## Filter LAS to class 2 and 8 (key points), LAS 1.4 now has class 8 as reserved
            # LASTools requires " around file names with spaces, ' not allowed, (Windows Command Line too?)
            # Use pdal for this? https://pdal.io/en/stable/apps/translate.html#example-1
            log.debug('--- Use las2las at ' + time.asctime())
            rc = call(os.path.join(softwareDir, 'LASTools', 'bin', 'las2las') + ' -i "' + allLAZ + '" -keep_class 2 8 -o ' + cl2LAS, shell=True)
            if rc == 0 and not os.path.isfile(cl2LAS):
                log.warning('las2las did not create class 2 points file: ' + cl2LAS)
                lasMP = None
            elif rc == 0:
                log.debug('--- Create las non-Minnesota Multipoint at ' + time.asctime())
                lasMP = arcpy.LASToMultipoint_3d(cl2LAS, inm + "\\pts" + sfx, "1", input_coordinate_system = srOut)
            else:
                log.warning('las2las did not execute successfully')

        ##        os.remove(cl2LAS)

        elif allLAZ.endswith('.zlas'):
            log.debug('--- Create zlas non-Minnesota Multipoint at ' + time.asctime())
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

##def convertDEMunitsIfNeeded(srVCS, coDEM):
def convertDEMunitsIfNeeded(srVCS, coDEM, state, log):
    '''convert units, but not to centimeters'''
    log.warning('coDEM is: ' + coDEM)
    if srVCS is not None:
        log.warning('srVCS: ' + str(srVCS.factoryCode) + ', ' + str(srVCS.name) + ', ' + srVCS.linearUnitName)
        unitsLower = srVCS.linearUnitName.lower()
        if 'foot' in unitsLower or 'feet' in unitsLower:
            log.warning('Converting DEM')# at ' + str(time.clock()))
            if 'survey' in unitsLower and 'foot' in unitsLower:
                regionDEM = Raster(coDEM)* 1200.0/3937.0
            elif 'foot' in unitsLower or 'feet' in unitsLower:
                regionDEM = Raster(coDEM)*0.3048
        elif 'meter' in unitsLower:
            regionDEM = Raster(coDEM)

    elif state == 'Illinois' or state == 'Wisconsin':
        log.warning('srVCS is None for ' + coDEM)
        # assume data in feet
        regionDEM = Raster(coDEM)*0.3048
    elif state == 'Missouri':
        log.warning('srVCS is None for ' + coDEM)
        # assume data in meters
        regionDEM = Raster(coDEM)

    return regionDEM

def setupPointsAndBreaklines(finalMP, inm, FDSet, breakpolys, breaklines, huc8, huc12, log):
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
        # breakGdb = os.path.join(os.path.dirname(blFile), "breaks_" + huc8 + ".gdb")
        breakGdb = os.path.dirname(breaklines)
        if not arcpy.Exists(breakGdb):
            if not os.path.isdir(os.path.basename(breakGdb)):
                os.makedirs(os.path.dirname(breakGdb))
            breakGdbResult = arcpy.CreateFileGDB_management(os.path.dirname(breakGdb), os.path.basename(breakGdb))

        log.warning(f'breakpolyList: {breakpolyList}')
        if len(breakpolyList) > 0:
            # mergedBreakpolys = arcpy.Merge_management(breakpolyList, os.path.join(breakGdb, 'break_polys_' + huc12))
            mergedBreakpolys = arcpy.Merge_management(breakpolyList, breakpolys)

        hlList = arcpy.ListFeatureClasses('hl_*', feature_type = 'POLYLINE', feature_dataset = os.path.basename(fd_string))
        log.warning(f'hlList: {hlList}')
        if len(hlList) > 0:
            finalHl = arcpy.Merge_management(hlList, os.path.join(str(FDSet), 'hl_merge'))

            copiedBreaklines = arcpy.CopyFeatures_management(finalHl, breaklines)
            # copiedBreaklines = arcpy.CopyFeatures_management(finalHl, os.path.join(breakGdb, 'break_lines_' + huc12))
            log.warning(f'copiedBreaklines: {copiedBreaklines}')
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


def generateLasArea(tilesClip, FDSet):
##    try:
    # buffer/debuffer this by 2 meters to get rid of some gaps in lidar file characterization
    tilesClipBuffer = arcpy.Buffer_analysis(tilesClip, buffer_distance_or_field = '2 METERS', dissolve_option = 'ALL')
    tilesClipDeBuffer = arcpy.Buffer_analysis(tilesClipBuffer, buffer_distance_or_field = '-2 METERS')
    tilesClipDslv = arcpy.Dissolve_management(tilesClipDeBuffer)
    tilesClipDslvElim = arcpy.EliminatePolygonPart_management(tilesClipDslv, condition = 'PERCENT', part_area_percent = 50)
    tcdFdSet = arcpy.CopyFeatures_management(tilesClipDslvElim, os.path.join(str(FDSet), 'local_las_area'))

    return tcdFdSet

def updateResolution(filename, init_res, new_res, huc12, log):
    """Take a filename with a specified resolution and alter it to the current processing resolution.
    This is done to reduce the number of arguments that are passed to the program."""
    # try:
    if init_res != new_res:
        updated_filename = filename.replace(str(init_res) + 'm' + huc12, str(new_res) + 'm' + huc12)
        log.debug(f"filename was: {filename}; updated filename: {updated_filename}")
    else:
        updated_filename = filename
    return updated_filename
    # except:



def buildLASRasters(lasdAll, lasdGround, log, demList, huc12, srSfx, maskRastBase, sgdb, procDir, int1rMaxFile, int1rMinFile, surfaceElevFile, intBeMaxFile, bareEarthReturnMinFile, cnt1rFile, init_cellSize, internal_regions, ptr):
##def buildLASRasters(lasdAll, lasdGround, log, demList, huc12, srSfx, maskRastBase, sgdb, procDir, int1rMaxFile, int1rMinFile, surfaceElevFile, frMinFile, intBeMaxFile, intBeMinFile, lastReturnMinFile, bareEarthReturnMinFile, cnt1rFile, init_cellSize, int_regions, ptr):
    '''creates multiple rasters from a las dataset, including min/max intensity of
    first return and bare earth surfaces, first return max and min surface, and z_range'''
    try:
        maskRastOut = opj(sgdb, maskRastBase + str(demList[0]))

        log.debug('snapRaster for LAS Dataset to raster: ' + arcpy.env.snapRaster)

        # procDir = locDict['fProcDir']
        int1rMaxFile_sized = updateResolution(int1rMaxFile, init_cellSize, demList[0], huc12, log)
        int1rMinFile_sized = updateResolution(int1rMinFile, init_cellSize, demList[0], huc12, log)
        frMaxFile_sized = updateResolution(surfaceElevFile, init_cellSize, demList[0], huc12, log)
##        frMinFile_sized = updateResolution(frMinFile, init_cellSize, demList[0], huc12, log)
        intBeMaxFile_sized = updateResolution(intBeMaxFile, init_cellSize, demList[0], huc12, log)
##        intBeMinFile_sized = updateResolution(intBeMinFile, init_cellSize, demList[0], huc12, log)
##        lastReturnMinFile_sized = updateResolution(lastReturnMinFile, init_cellSize, demList[0], huc12, log)
        bareEarthReturnMinFile_sized = updateResolution(bareEarthReturnMinFile, init_cellSize, demList[0], huc12, log)
        cnt1rFile_sized = updateResolution(cnt1rFile, init_cellSize, demList[0], huc12, log)


        log.warning('---Creating FR Max Intensity at ' + time.asctime())
        lasd1rMaxIntensity = arcpy.LasDatasetToRaster_conversion(lasdAll, int1rMaxFile_sized, 'INTENSITY', 'BINNING MAXIMUM NONE', sampling_type = 'CELLSIZE', sampling_value = demList[0], data_type = 'INT')
        recode_tf = False
        log.debug(f'ir.max: {internal_regions.maximum},ir.min: {internal_regions.minimum}')
        if internal_regions.maximum - internal_regions.minimum != 0:
            log.info('multiple regions')
            int_zs_max = ZonalStatistics(internal_regions, 'VALUE', int1rMaxFile_sized, 'MAXIMUM')
            if int_zs_max.minimum < 256 and int_zs_max.maximum > 256:
                int_lt_256 = LessThan(int_zs_max, 256)
                recode_areas = ZonalStatistics(ptr, 'VALUE', int_lt_256, 'MAXIMUM')
                multiplied_intensities = Raster(int1rMaxFile_sized) * 256
                recoded_intensities = Con(recode_areas, multiplied_intensities, int1rMaxFile_sized)
                recoded_intensities.save(int1rMaxFile_sized)
                recode_tf = True
            else:
                log.info('all regions equal max intensity')
        else:
            log.info('one region')

        log.warning('---Creating FR Min Intensity at ' + time.asctime())
        lasd1rMinIntensity = arcpy.LasDatasetToRaster_conversion(lasdAll, int1rMinFile_sized, 'INTENSITY', 'BINNING MINIMUM NONE', sampling_type = 'CELLSIZE', sampling_value = demList[0], data_type = 'INT')
        if recode_tf:
            multiplied_intensities = Raster(int1rMinFile_sized) * 256
            recoded_intensities = Con(recode_areas, multiplied_intensities, int1rMinFile_sized)
            recoded_intensities.save(int1rMinFile_sized)



        log.warning('---Creating FR Max surface at ' + time.asctime())
        allReturnsMaxTempFile = os.path.join(procDir, '_'.join(['tmp_frmax', str(demList[0]) + 'm', huc12, 'out.tif']))
        allReturnsMax = arcpy.LasDatasetToRaster_conversion(lasdAll, allReturnsMaxTempFile, interpolation_type = 'BINNING MAXIMUM SIMPLE', sampling_type = 'CELLSIZE', sampling_value = demList[0], data_type = 'FLOAT')
        allReturnsMaxCm = Int(Times(allReturnsMax, 100))
        allReturnsMaxCm.save(frMaxFile_sized)#locDict['surfaceElevFile'])#allReturnsMaxFile)

        # log.warning('---Creating FR Min surface at ' + time.asctime())
        # allReturnsMinTempFile = os.path.join(procDir, '_'.join(['tmp_frmin', str(demList[0]) + 'm', huc12, 'out.tif']))
        # allReturnsMin = arcpy.LasDatasetToRaster_conversion(lasdAll, allReturnsMinTempFile, interpolation_type = 'BINNING MINIMUM SIMPLE', sampling_type = 'CELLSIZE', sampling_value = demList[0], data_type = 'FLOAT')
        # allReturnsMinCm = Int(Times(allReturnsMin, 100))
        # allReturnsMinCm.save(frMinFile_sized)#locDict['firstReturnMinFile'])#allReturnsMinFile)

        log.warning('---Counting All Returns at ' + time.asctime())
        cfrFileTemp = 'cnt_fr_' + str(demList[0]) + "m_" + huc12 + srSfx + '.tif'
        lasdCount = arcpy.LasPointStatsAsRaster_management(lasdAll, os.path.join(procDir, cfrFileTemp), 'POINT_COUNT', 'CELLSIZE', demList[0])
        cfrFileRasterObj = clipCountRaster(lasdCount, maskRastOut, cnt1rFile_sized)

        # log.warning('---Counting Z Range at ' + time.asctime())
        # zrangeFileTemp = 'zrng_all_' + str(demList[0]) + "m_" + huc12 + srSfx + '.tif'
        # lasdCount = arcpy.LasPointStatsAsRaster_management(lasdAll, os.path.join(procDir, zrangeFileTemp), 'Z_RANGE', 'CELLSIZE', demList[0])

##        lastLayer = arcpy.MakeLasDatasetLayer_management(lasdOut, 'ground_layer', [2,8], 'Last Return')
        # if sys.version_info.minor < 9:
        #     lastLayer = arcpy.MakeLasDatasetLayer_management(lasdGround, 'last_layer', return_values = 'Last Return')
        # else:
        #     lastLayer = arcpy.MakeLasDatasetLayer_management(lasdGround, 'last_layer', return_values = 'LAST')
        # log.warning('---Creating BE Max Intensity at ' + time.asctime())
        # # minimum be seems to be a little less 'noisy' than maximum
        # lasdBeMaxIntensity = arcpy.LasDatasetToRaster_conversion(lastLayer, intBeMaxFile_sized, 'INTENSITY', 'BINNING MAXIMUM NONE', sampling_type = 'CELLSIZE', sampling_value = demList[0], data_type = 'INT')

        # log.warning('---Creating BE Min Intensity at ' + time.asctime())
        # # minimum be seems to be a little less 'noisy' than maximum
        # lasdBeMinIntensity = arcpy.LasDatasetToRaster_conversion(lastLayer, intBeMinFile_sized, 'INTENSITY', 'BINNING MINIMUM NONE', sampling_type = 'CELLSIZE', sampling_value = demList[0], data_type = 'INT')

        # log.warning('---Creating LR Min surface at ' + time.asctime())
        # lastReturnsMinTempFile = os.path.join(procDir, '_'.join(['tmp_lrmin', str(demList[0]) + 'm', huc12, 'out.tif']))
        # lastReturnsMin = arcpy.LasDatasetToRaster_conversion(lastLayer, lastReturnsMinTempFile, interpolation_type = 'BINNING MINIMUM SIMPLE', sampling_type = 'CELLSIZE', sampling_value = demList[0], data_type = 'FLOAT')
        # lastReturnsMinCm = Int(Times(lastReturnsMin, 100))
        # lastReturnsMinCm.save(lastReturnMinFile_sized)#locDict['lastReturnMinFile'])#.replace('fr', 'lr'))

        log.warning('---Creating LR Min surface at ' + time.asctime())
        if sys.version_info.minor < 9:
            beLayer = arcpy.MakeLasDatasetLayer_management(lasdGround, 'ground_layer', [2,8], 'Last Return')
        else:
            beLayer = arcpy.MakeLasDatasetLayer_management(lasdGround, 'ground_layer', [2,8], 'LAST')

        beReturnsMinTempFile = os.path.join(procDir, '_'.join(['tmp_bemin', str(demList[0]) + 'm', huc12, 'out.tif']))
        beReturnsMin = arcpy.LasDatasetToRaster_conversion(beLayer, beReturnsMinTempFile, interpolation_type = 'BINNING MINIMUM NONE', sampling_type = 'CELLSIZE', sampling_value = demList[0], data_type = 'FLOAT')
        beReturnsMinCm = Int(Times(beReturnsMin, 100))
        beReturnsMinCm.save(bareEarthReturnMinFile_sized)#locDict['bareEarthReturnMinFile'])#.replace('fr', 'be'))

        lasdBeMaxIntensity = arcpy.LasDatasetToRaster_conversion(beLayer, intBeMaxFile_sized, 'INTENSITY', 'BINNING MAXIMUM NONE', sampling_type = 'CELLSIZE', sampling_value = demList[0], data_type = 'INT')
        if recode_tf:
            multiplied_intensities = Raster(intBeMaxFile_sized) * 256
            recoded_intensities = Con(recode_areas, multiplied_intensities, intBeMaxFile_sized)
            recoded_intensities.save(intBeMaxFile_sized)

##        # my current favorites  - minmax mild 18 terrain, binning min simple (almost too detailed),
##        interps = ['BINNING MINIMUM SIMPLE', 'BINNING AVERAGE SIMPLE', 'BINNING IDW SIMPLE', 'TRIANGULATION NATURAL_NEIGHBOR WINDOW_SIZE MINIMUM ' + str(demList[0]), 'TRIANGULATION LINEAR WINDOW_SIZE MINIMUM ' + str(demList[0]), 'TRIANGULATION NATURAL_NEIGHBOR WINDOW_SIZE CLOSEST_TO_MEAN ' + str(demList[0]), 'TRIANGULATION LINEAR WINDOW_SIZE CLOSEST_TO_MEAN ' + str(demList[0])]
##
##        log.warning('---Creating extra LR Min surfaces at ' + time.asctime())
##        for interp in interps:
##            log.warning('interp is ' + str(interp))
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


def loadStateAbbrevs(state):

    us_state_abbrev = {
        'Alabama': 'AL',
        'Alaska': 'AK',
        'Arizona': 'AZ',
        'Arkansas': 'AR',
        'California': 'CA',
        'Colorado': 'CO',
        'Connecticut': 'CT',
        'Delaware': 'DE',
        'Florida': 'FL',
        'Georgia': 'GA',
        'Hawaii': 'HI',
        'Idaho': 'ID',
        'Illinois': 'IL',
        'Indiana': 'IN',
        'Iowa': 'IA',
        'Kansas': 'KS',
        'Kentucky': 'KY',
        'Louisiana': 'LA',
        'Maine': 'ME',
        'Maryland': 'MD',
        'Massachusetts': 'MA',
        'Michigan': 'MI',
        'Minnesota': 'MN',
        'Mississippi': 'MS',
        'Missouri': 'MO',
        'Montana': 'MT',
        'Nebraska': 'NE',
        'Nevada': 'NV',
        'New Hampshire': 'NH',
        'New Jersey': 'NJ',
        'New Mexico': 'NM',
        'New York': 'NY',
        'North Carolina': 'NC',
        'North Dakota': 'ND',
        'Ohio': 'OH',
        'Oklahoma': 'OK',
        'Oregon': 'OR',
        'Pennsylvania': 'PA',
        'Rhode Island': 'RI',
        'South Carolina': 'SC',
        'South Dakota': 'SD',
        'Tennessee': 'TN',
        'Texas': 'TX',
        'Utah': 'UT',
        'Vermont': 'VT',
        'Virginia': 'VA',
        'Washington': 'WA',
        'West Virginia': 'WV',
        'Wisconsin': 'WI',
        'Wyoming': 'WY',
        'District of Columbia': 'DC',
        'Northern Mariana Islands':'MP',
        'Palau': 'PW',
        'Puerto Rico': 'PR',
        'Virgin Islands': 'VI'
    }
    abbrev = us_state_abbrev.get(state)
    return abbrev

def generateTempTerrName(procDir, window, cellsize, huc12):
    tempTerrName = os.path.join(procDir, '_'.join(['tmp_ter', window, str(cellsize) + 'm', huc12, 'out.tif']))

    return tempTerrName


def mosaicDEMsAndPitfill(demList, maskRastBase, huc12, log, sgdb, windows, procDir, fElevFile, interpDict, init_cellSize, srOutNoVCS, LTrrn, pyramids, dem_metadata_template, lidar_metadata_info):
    '''Takes whole or partial DEMs from the demList and mosaics them together
    if there are multiple DEMs. Then pit-fills the result (fills all one cell
    sinks). Also processes 'ZMEAN' and 'ZMINMAX' (or other terrain->raster
    options if they're added to the code) to separate directories.'''

    try:
        arcpy.env.cellSize = demList[0]
        maskRastOut = opj(sgdb, maskRastBase + str(demList[0]))

        noLASdem = demList[1:]#[]
        log.warning('noLASdem (if present) is: ' + str(noLASdem))

        # windows are types of terrain (ZMEAN, ZMINMAX, etc.)
        if len(windows):
            for window in windows:

                fElevFile = updateResolution(fElevFile, init_cellSize, demList[0], huc12, log)

                interpType = interpDict[window]
                # default interpolation type is mean18
                if interpType != 'mean18':
                    fElevFile_interp = fElevFile.replace('mean18', interpType)
                else:
                    fElevFile_interp = fElevFile


                # paths = df.loadVariablesDict(node, ACPFyear, huc12, srOutCode, interpType, demList[0], nowYmd, version)

                # procDir = paths['fProcDir']

                tempTerrName = generateTempTerrName(procDir, window, demList[0], huc12)
                # get rasters created from Terrain, should only be one
                arcpy.env.workspace = os.path.dirname(tempTerrName)
                rastrList = arcpy.ListRasters(os.path.basename(tempTerrName))#[0]

                if len(noLASdem) > 0:
                    mosaicList = rastrList + noLASdem
                    if True:#len(mosaicList) > 1:
                        log.debug('resampling method for mosaic to new raster: ' + arcpy.env.resamplingMethod)
                        log.debug('cellsize for mosaic to new raster: ' + arcpy.env.cellSize)
                        log.warning('---Mosaicing DEMs for ' + str(demList[0]) + ' at ' + time.asctime())
                        mosDEM = arcpy.MosaicToNewRaster_management(mosaicList, procDir, 'mos_' + str(demList[0]) + 'm_' + huc12 + '.tif', number_of_bands = '1', pixel_type = '32_BIT_FLOAT', mosaic_method = 'MINIMUM')
                        maskedDEM = Con(Plus(maskRastOut, 2), mosDEM)
                        maskedDEMintPre = Int(maskedDEM*100)

                        # filter any null values in mask (buffered HUC12) using Nibble
                        isnCmDemInt = IsNull(maskedDEM)
                        ndToNibble = SetNull(isnCmDemInt == 1, 1)
                        cmDEMintNoNulls = Con(isnCmDemInt == 0, maskedDEMintPre, 1)
                        nibbleInt = Nibble(cmDEMintNoNulls, ndToNibble)
                        maskedDEMint = Con(Plus(maskRastOut, 2), nibbleInt)

                    else:
                        log.debug('resampling method for Resample: ' + arcpy.env.resamplingMethod)
                        log.debug('cellsize for resample: ' + arcpy.env.cellSize)
                        log.warning('---Resampling DEMs for ' + str(demList[0]) + ' at ' + time.asctime())
                        mosDEM = arcpy.Resample_management(mosaicList[0], os.path.join(procDir, 'mos_' + str(demList[0]) + 'm_' + huc12 + '.tif'), cell_size = arcpy.env.cellSize, resampling_type = arcpy.env.resamplingMethod)
                        maskedDEM = Con(Plus(maskRastOut, 2), mosDEM)
                        maskedDEMintPre = Int(maskedDEM*100)

                        # filter any null values in mask (buffered HUC12) using Nibble
                        isnCmDemInt = IsNull(maskedDEM)
                        ndToNibble = SetNull(isnCmDemInt == 1, 1)
                        cmDEMintNoNulls = Con(isnCmDemInt == 0, maskedDEMintPre, 1)
                        nibbleInt = Nibble(cmDEMintNoNulls, ndToNibble)
                        maskedDEMint = Con(Plus(maskRastOut, 2), nibbleInt)

                else:
                    rastr = rastrList[0]
                    mosDEM = Raster(rastr)#tempTerrName)
                    intCmDEM = Int(mosDEM*100)
                    maskedDEMint = ExtractByMask(intCmDEM, maskRastOut)

                # clear output VCS from here on out, creating centimeter rasters
                arcpy.env.outputCoordinateSystem = None
                arcpy.env.outputCoordinateSystem = srOutNoVCS

                log.warning('msk_dem to be saved now')
                maskedDEMint.save('msk_dem_' + str(demList[0]) + 'm_' + huc12 + '.tif')


                cmDEMnocs, cmDEMsinks = fillOCSinks(maskedDEMint)
                log.warning('pfFile name will be ' + fElevFile_interp)#paths['fElevFile'])
                cmDEMnocs.save(fElevFile_interp)#paths['fElevFile'])
                log.warning('Saved DEM for ' + str(demList[0]))# + ' at ' + str(time.clock()))
                arcpy.BuildPyramids_management(cmDEMnocs)

                terrain_args, nowYmd, collect_starts_min, collect_ends_max, collect_majority, pyramid_args = [
                    lidar_metadata_info[i] for i in lidar_metadata_info]

                if 'ZMINMAX'.lower() in fElevFile_interp:
                    terrain_args_updated = terrain_args.replace('MEAN', 'MINMAX')
                else:
                    terrain_args_updated = terrain_args

                paraDict = {
                    '\n\nACPF: DEM Generation and Pit Fill Tool     ' : '\nRun Date: %s' % nowYmd,
                    '\nUnknown Vintage Lidar Data' : False,#tiles_t_or_f,
                    '\nEarliest 3DEP Lidar Data: ' : collect_starts_min,
                    '\nLatest 3DEP Lidar Data: ' : collect_ends_max,
                    '\nLatest 3DEP Lidar Data: ' : collect_majority,
                    '\nOutput conditioned DEM raster: ' : fElevFile_interp,
                    '\nTerrain Interpolation Arguments: ' : terrain_args_updated,
                    '\nTerrain Pyramid Arguments: ' : pyramid_args
                    }

##                template_interp = r'\\EL3354-02\O$\DEP\toolMetadata\F-Lib1mDEM2020_metadata.xml'

                ## update metadata
                log.warning('---Adding metadata at ' + time.asctime())
                addMetadata(fElevFile_interp, paraDict, dem_metadata_template, log)

                # ## update metadata
                # log.warning('---Updating metadata at ' + time.asctime())

                # if 'basics' in os.getcwd():
                #     updateScript = os.path.join(os.getcwd(), 'update_metadata.py')
                # else:
                #     updateScript = os.path.join(os.getcwd(), 'basics', 'update_metadata.py')

                # mdversion = '10.3'
                # callList = ['C:\\Python27\\ArcGIS' + mdversion + '\\pythonw.exe', updateScript, fElevFile_interp, os.path.basename(sys.argv[0])]
                # callString = ' '.join(callList)

                # if sys.winver.startswith('2'):
                #     try:
                #         stdout=subprocess.check_output(callString, stderr=subprocess.STDOUT, universal_newlines=True)#, env=env)
                #     except subprocess.CalledProcessError as e:
                #         log.warning('stdout output:\n' + e.output)
                # elif sys.winver.startswith('3'):

                #     subprocess.run(callString)

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


def addMetadata(outDEM, paraDict, template_file_path, log):
    # Set the standard-format metadata XML file's path
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




def lidarDEM(ept_wesm_file, cleanup, messages):
    import sys, os, time
    import subprocess
    from subprocess import call
    import platform, glob, traceback
    from os.path import join as opj
    if sys.version_info.major == 2:
        import getpass
        login = getpass.getuser()
    else:
        login = os.getlogin()
        
    if login == 'bkgelder':
        boxes = ['C:\\Users\\bkgelder\\Box\\Data_Sharing\\Scripts\\basics', 'O:\\DEP\\Scripts\\basics']
    else:
        boxes = ['C:\\Users\\idep2\\Box\\Scripts\\basics', 'O:\\DEP\\Scripts\\basics']

    for box in boxes:
        if os.path.isdir(box):
            sys.path.append(box)

    import dem_functions2 as df




    return

if __name__ == "__main__":
    import sys


    if len(sys.argv) == 1:
        arcpy.AddMessage("Whoo, hoo! Running from Python Window!")
        cleanup = False

        parameters = ["C:/Program Files/ArcGIS/Pro/bin/Python/envs/arcgispro-py3/pythonw.exe",
'C:\\DEP\\Scripts\\basics\\cmd_BuildDEM_newArgs.py', '\\\\EL3354-02\\O$\\DEP\\Basedata_Summaries\\Basedata_26916.gdb\\MW_HUC12_v2019', '\\\\EL3354-02\\O$\\DEP\\Basedata_Summaries\\Basedata_26916.gdb/Snap1m', '\\\\EL3354-02\\O$\\DEP\\Elev_Base_Data\\ept\\ept.gdb\\ept_resources_2023_04_01', '\\\\EL3354-02\\O$\\DEP\\toolMetadata\\FLib_DEMs2022_mTemplate.xml', '\\\\EL3354-02\\O$\\DEP\\toolMetadata\\FLib_Derivatives2022_mTemplate.xml', 'D:\\DEP_Proc\\DEMProc\\LAS_dem2013_3m_040302020105', '\\\\EL3354-02\\O$\\DEP\\Elev_Base_Data', '\\\\EL3354-02\\O$\\DEP\\LiDAR_Current\\elev_FLib_mean18\\04030202\\ef3m040302020105.tif', '\\\\EL3354-02\\O$\\DEP\\LiDAR_Current\\surf_el_Lib\\04030202\\bemin3m040302020105.tif', '\\\\EL3354-02\\O$\\DEP\\LiDAR_Current\\surf_el_Lib\\04030202\\frmax3m040302020105.tif', '\\\\EL3354-02\\O$\\DEP\\LiDAR_Current\\count_Lib\\04030202\\cbe3m040302020105.tif', '\\\\EL3354-02\\O$\\DEP\\LiDAR_Current\\count_Lib\\04030202\\cfr3m040302020105.tif', '\\\\EL3354-02\\O$\\DEP\\LiDAR_Current\\int_Lib\\04030202\\fr_int_min3m040302020105.tif', '\\\\EL3354-02\\O$\\DEP\\LiDAR_Current\\int_Lib\\04030202\\fr_int_max3m040302020105.tif', '\\\\EL3354-02\\O$\\DEP\\LiDAR_Current\\int_Lib\\04030202\\be_int_max3m040302020105.tif', '\\\\EL3354-02\\O$\\DEP\\LiDAR_Current\\bl_Lib\\04030202\\breaks_04030202.gdb\\break_polys_040302020105', '\\\\EL3354-02\\O$\\DEP\\LiDAR_Current\\bl_Lib\\04030202\\breaks_04030202.gdb\\break_lines_040302020105', '\\\\EL3354-02\\D$\\DEP\\Man_Data_ACPF\\dep_ACPF2021\\04030202\\idepACPF040302020105.gdb\\wesm_ept_resources_2023_04_01_040302020105']

        for i in parameters[2:]:
            sys.argv.append(i)

    else:
        arcpy.AddMessage("Whoo, hoo! Command-line enabled!")
        # clean up the folder after done processing
        cleanup = True

    # ept_fc = "C:/DEP/Elev_Base_Data/ept/ept.gdb/ept_resources_2023_05_20"
    ept_fc = sys.argv[1]

    # inputs then outputs
    (dem_fc, snap, ept_wesm_file, flib_metadata_template, derivative_metadata_template,
        procDir, eleBaseDir,
        fElevFile, bareEarthReturnMinFile, firstReturnMaxFile, cntFile,cnt1rFile,
        int1rMinFile, int1rMaxFile, intBeMaxFile, breakpolys, breaklines, wesm_project_file
        ) = (sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5],
            sys.argv[6], sys.argv[7], 
            sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11], sys.argv[12],
            sys.argv[13], sys.argv[14], sys.argv[15], sys.argv[16], sys.argv[17], sys.argv[18])


    # doLidarDEM(ept_fc, cleanup, msgStub())

    # arcpy.AddMessage("Back from doEPT!")