##This program was written to help remove issues in lidar DEMs caused
##by poor performance of interpolation algorithms in specific lcocations
##It is part of a series of DEM processing algorithms designed to improve
##flow of water across an elevaiton surface.
##
##Primarily written by Brian Gelder, bkgelder@iastate.edu
##Substantial help/advice from David James
##
##Separate code section created on 2019 February 26.
##for Python 2.7 and ArcGIS 10.3.1
##
# removal watershed calculations for distingushing between two fill regions with
#   same max depth. now uses a second region group to code them uniquely. 2019/04/26 bkgelder
##
## MAJOR bug - punched DEMs still had holes after IsNull code
## fixed to test null values appropriately on 2019/08/02 bkgelder
## added code to save 'holes' out as separate dataset 2019/08/13 bkgelder
## 2019.09.12 - added log.warning('Invalid Topology [Incomplete void poly] for rtp, attempting repair'), bkgelder
## 2020-03-03 cleanup of unused arguments, added comments bkgelder
## 2021-01-20 paths now loaded internall via arguments passed via command line - bkgelder
## 2021.10.22 - now can use void fixed DEMs as potential inputs
## 2022.06.03 - added additional region group to make separate sinks at same elevation in one fill region
##              this was causing errors in cutting as they were both at same cut level
## 2024.02.11 - rewritten for ArcGIS Pro and Python 3

# Import system modules
import arcpy
import sys
import os 
import traceback
import time
import platform

# from shapely import area

sys.path.append("C:\\DEP\\Scripts\\basics")
import dem_functions as df

from arcpy.sa import *

from os.path import join as opj

class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "CMD_Puncher"
        self.alias = "CMD_Puncher"
        # List of tool classes associated with this toolbox
        self.tools = [Tool]


class Tool(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Hole_Puncher"
        self.description = "Punches holes in a DEM and fills remaining to remove depressions shallower than a criteria"
        self.canRunInBackground = False
        self.category = "DEM Puncher"

    def getParameterInfo(self):
        """Define parameter definitions"""

        param0 = arcpy.Parameter(
            name="input_dem",
            displayName="Input Elevation Model",
            datatype="DERasterDataset",
            parameterType='Required',
            direction="Input")
        
        param1 = arcpy.Parameter(
            name="output_dem",
            displayName="Output Punched Elevation Model",
            datatype="DERasterDataset",
            parameterType='Required',
            direction="Output")
        
        param2 = arcpy.Parameter(
            name="plib_metadata",
            displayName="Punched DEM metadata template",
            datatype="DEFile",
            parameterType='Required',
            direction="Input")
        
        param3 = arcpy.Parameter(
            name="depressions_fc",
            displayName="Punched depressions feature class",
            datatype="DEFeatureClass",
            parameterType='Required',
            direction="Output")
        
        param4 = arcpy.Parameter(
            name="depth_threshold",
            displayName="Depression Punch Depth Threshold",
            datatype="GPString",
            parameterType='Required',
            direction="Input")
        
        param5 = arcpy.Parameter(
            name="area_threshold",
            displayName="Depression Punch Area Threshold",
            datatype="GPString",
            parameterType='Optional',
            direction="Input")
        
        param6 = arcpy.Parameter(
            name = "procDir",
            displayName="Local Processing Directory",
            datatype="DEFolder",
            parameterType='Optional',
            direction="Input")
        
        parameters = [param0, param1, param2, param3, param4, param5, param6]
        return parameters
        

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
        params = parameters
        doPuncher(params[0].valueAsText, params[1].valueAsText, params[2].valueAsText, params[3].valueAsText, params[4].valueAsText, params[5].valueAsText, params[6].valueAsText, cleanup, messages)
        return

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return



def doPuncher(input_dem, output_dem, plib_metadata, depressions_fc, depth_threshold, area_threshold, procDir, cleanup, messages):

    try:
        arguments = [input_dem, output_dem, plib_metadata, depressions_fc, depth_threshold, area_threshold, procDir, cleanup]

        for a in arguments:
            if a == arguments[0]:
                arg_str = str(a) + '\n'
            else:
                arg_str += str(a) + '\n'

        messages.addMessage("Tool: Executing with parameters:\n" + arg_str)
        print(arg_str)

        huc12, huc8 = df.figureItOut(input_dem)
        print(huc12)

        # control where scratchFolder and GDB are created
        if procDir != "":
            if os.path.isdir(procDir):
                # log.info('nuking: ' + procDir)
                df.nukedir(procDir)

            if not os.path.isdir(procDir):
                os.makedirs(procDir)

            arcpy.env.scratchWorkspace = procDir
        print(procDir)

        # if no procDir, next step will create in users directory
        # otherwise, previously set scratchworkspace will control location
        sfldr = arcpy.env.scratchFolder
        sgdb = arcpy.env.scratchGDB
        arcpy.env.scratchWorkspace = sgdb
        arcpy.env.workspace = sgdb
        print(sfldr)

        #figure out where to create log files
        node = platform.node()
        logProc = df.defineLocalProc(node)
        if not os.path.isdir(logProc):
            logProc = sfldr
        print(logProc)

        if cleanup:
            # log to file only
            log, nowYmd, logName, startTime = df.setupLoggingNoCh(logProc, sys.argv[0], huc12)
            verbose = False
            arcpy.SetLogHistory = False
        else:
            # log to file and console
            log, nowYmd, logName, startTime = df.setupLoggingNew(logProc, sys.argv[0], huc12)
            verbose = True
            arcpy.SetLogHistory = True

        startTime = time.time()
        log.info("Beginning execution:")
        log.info("Tool: Executing with parameters:\n" + arg_str)
        messages.addMessage("Log file at " + logName)

        ## Make sure depression output locations exist
        output_gdb = os.path.dirname(depressions_fc)
        folder_dir = os.path.dirname(output_gdb)
        if not os.path.isdir(folder_dir):
            os.makedirs(folder_dir)
        if not arcpy.Exists(output_gdb):
            arcpy.CreateFileGDB_management(folder_dir, os.path.basename(output_gdb))

        df.create_needed_dirs_and_gdbs(output_dem, log)

        # if os.path.isdir(procDir):
        #     log.warning('nuking: ' + procDir)
        #     df.nukedir(procDir)

        # if not os.path.isdir(procDir):
        #     os.makedirs(procDir)

        # ## Set the environments
        # arcpy.env.scratchWorkspace = procDir

        # sfldr = arcpy.env.scratchFolder
        # sgdb = arcpy.env.scratchGDB
        # arcpy.env.scratchWorkspace = sfldr#
        # arcpy.env.workspace = sgdb

        inm = 'in_memory'

        arcpy.env.snapRaster = input_dem

        arcpy.env.cellSize = input_dem#

        arcpy.CheckOutExtension("Spatial")
        arcpy.env.overwriteOutput = True

        frFld = 'FILL_RGN'
        fillLvlFld = 'FILL_LVL'

        minElFld = 'FR_MIN_EL'

        frDepthFld = 'FR_DEPTH'
        frVolFld = 'FR_VOLUME'
        frMaxSlopeFld = 'FR_MAX_SLP'
        frMeanSlopeFld = 'FR_MEAN_SLP'

        wsLvlFld = 'ws_lvl'

        gridfield = 'gridcode'

        dfsList = []
        wsList = []
        fr0List = []

        ######------------------------------------------------------------------------------

        ## Begin running the fast hole punching algorithm for first time
        index = 0
        sfx = '_' + str(index)

        meterNdDEM = 0.01 * Raster(input_dem)
        ProcSize = meterNdDEM.meanCellHeight
        slopePct = Slope(meterNdDEM, 'PERCENT_RISE')
        slopePct.save(opj(procDir, 'slope_pct'))

        ## Ideas for filtering noise from the DEM
        # also try Perona Malik filter???
        # fst5x5Std = FocalStatistics(bestestDEM, NbrRectangle(5, 5, 'CELL'), 'STD')

        ## calculate relative depth raster, maximum depth value, and filled DEM for the input DEM
        rDepth, maxDepth, fillLvl = df.findMaxDepth(input_dem)
        rDepth.save(opj(procDir, 'rdpth' + sfx))
        log.info('Initial maximum depth: ' + str(maxDepth))
        fillLvl.save(opj(procDir, 'fill_lvl' + sfx))

        # refers to current holed out DEM
        inDEM = input_dem

        prevMaxValue = 0
        indexMax = 9 # stop punching after this many times
        while maxDepth > float(depth_threshold) and index < indexMax:

            log.info("Punch holes for " + str(sfx) + ' at ' + str(time.asctime()))
        ## Define fill regions (and cost surface)
            fillGT0 = Con(rDepth > 0, 1)

            rgMinFilledEl = RegionGroup(fillGT0, 'EIGHT')

        ## Calculate depth statistics for fill regions
            # could just use ZonalStatistics here if only depth punching desired
            zstMaxDepth = ZonalStatisticsAsTable(rgMinFilledEl, "VALUE", rDepth, 'in_memory\\zst_max_dpth')
            joinMax = arcpy.JoinField_management(rgMinFilledEl, "VALUE", zstMaxDepth, "VALUE", "MAX")

        ## Select the appropriate fill regions for enforcement (deep enough or large enough)
        ##        # fldName is need to fix boggling of max field name for some versions (10.5, 10.6, others?)
        ##        if version.find('10.5') > -1 or version.find('10.6') > -1:
            fldName = arcpy.ListFields(rgMinFilledEl, 'MAX*')[0].name
            rgToPunch = Con(rgMinFilledEl, rgMinFilledEl, "", '"' + fldName + '" > ' + str(depth_threshold) + ' OR "COUNT" * ' + str(float(ProcSize)**2) + ' > ' + str(area_threshold))
        ##        else:
        ##            rgToPunch = Con(rgMinFilledEl, rgMinFilledEl, "", '"MAX" > ' + str(depth_threshold) + ' OR "COUNT" * ' + str(ProcSize**2) + ' > ' + str(area_threshold))
            arcpy.BuildRasterAttributeTable_management(rgToPunch)
            rgToPunch.save(opj(procDir, 'rg2punch' + sfx))

        ## Define sink (and thus watershed) number by minimum elevation sink
            #(minimum value if tie for minimum to make one sink/hole the 'dominant' fill region
            rgMinSinkEl = ZonalStatistics(rgToPunch, 'VALUE', inDEM, 'MINIMUM')
            sinksAtMin = Con(rgMinSinkEl == inDEM, rgToPunch)
            #Unique numbering for each iteration by adding prevMaxValue, helps cutting process
            if index > 0:
                holes2PunchMulti = RegionGroup(sinksAtMin, "EIGHT") + prevMaxValue
            else:
                holes2PunchMulti = RegionGroup(sinksAtMin, "EIGHT")

            minHoles2Punch = ZonalStatistics(rgToPunch, 'VALUE', holes2PunchMulti, 'MINIMUM')
            holes2Punch = Con(minHoles2Punch == holes2PunchMulti, holes2PunchMulti)
            holes2Punch.save(opj(procDir, 'snkunq' + sfx))
            prevMaxValue = int(holes2Punch.maximum)

        ## Fill in the little fill regions (force flow to deeper/larger depressions), flow is not forced elsewhere
            rgToFillFilled = Con(IsNull(rgToPunch), fillLvl, inDEM)
            log.debug("After rgToFillFilled for " + str(sfx) + ' at ' + str(time.asctime()))

        ## Make holes in the DEM where we want them (from above)
            demWithHoles = Con(IsNull(holes2Punch), rgToFillFilled, '')
##            demWithHoles.save(opj(procDir, "nwdm" + sfx))



        ## Now do the next depth check because we're going to use the DEM for water flow
            ## The water in the next depth check flows to the holes just created 
            ## We'll use that to define the watersheds for that fill region
            log.debug("Before df.findMaxDepth( for " + str(sfx) + ' at ' + str(time.asctime()))
            ## re-do the calculations for the new DEM
            rDepthNext, maxDepthNext, fillLvlNext = df.findMaxDepth(demWithHoles)
            log.info('Maximum depth for ' + str(sfx) + ': ' + str(maxDepth) + ' at ' + str(time.asctime()))
            # log.debug("After df.findMaxDepth( for " + str(sfx) + ' at ' + str(time.asctime()))
            ## save fillLvl, rDepth with next iteration number since that's when they would be used
            fillLvl.save(opj(procDir, 'fill_lvl_' + str(index+1)))
            # rDepth.save(opj(procDir, 'rdpth_' + str(index+1)))


            
            fdTemp = FlowDirection(fillLvlNext)
            log.debug("After flowDirection for " + str(sfx) + ' at ' + str(time.asctime()))
        ##        fdTemp.save(opj(cp, 'fd_lvl0_' + str(index + 1))#' + sfx)#

            wsLvl = Watershed(fdTemp, holes2Punch)
            ## wsLvl is watersheds from fdTemp and cumSinks
            wsLvl.save(opj(procDir, wsLvlFld + sfx))
            log.debug("After wsLvl for " + str(sfx) + ' at ' + str(time.asctime()))

            ## Needed for matcher, not for puncher
            wsPolysMultiFr = arcpy.RasterToPolygon_conversion(wsLvl, opj(inm, 'ws_polys_multi_features' + sfx), 'SIMPLIFY')
            # 'dissolve' the features with the same gridcode (watershed number) into one single multipart feature in attribute table
            wsPolys = arcpy.Dissolve_management(wsPolysMultiFr, opj(inm, 'ws_polys' + sfx), gridfield)
            df.tryAddField(wsPolys, frFld, 'LONG')
            arcpy.CalculateField_management(wsPolys, frFld, '!' + gridfield + '!', 'python')
            wsList.append(wsPolys)
            df.copyfc(verbose, wsPolys, sgdb)

            wsMinFilledEl = ZonalStatistics(wsLvl, 'VALUE', fillLvl, 'MINIMUM')
            # wsMinFilledEl = ZonalStatistics(wsLvl, 'VALUE', fillLvlPrev, 'MINIMUM')
            log.debug("After wsMinFilledEl for " + str(sfx) + ' at ' + str(time.asctime()))
            # where fill greater than or equal zero, code to wsLvl/fill Region
            fr0 = Con(fillLvl == wsMinFilledEl, wsLvl)
            # fr0 = Con(fillLvlPrev == wsMinFilledEl, wsLvl)
            fr0.save(opj(procDir, 'fr0' + sfx))
            fr0List.append(fr0)

        ## Convert all fill regions to polygons to store search data
            try:
                rtp = arcpy.RasterToPolygon_conversion(fr0, opj(inm, 'rtp' + sfx), 'SIMPLIFY')
                dfsFC = arcpy.Dissolve_management(rtp, opj(inm, "dfs_frToPoly" + sfx), gridfield)
            except:
                log.warning('Invalid Topology [Incomplete void poly] for rtp, attempting repair')
                arcpy.RepairGeometry_management(rtp)
                dfsFC = arcpy.Dissolve_management(rtp, opj(inm, "dfs_frToPoly" + sfx), gridfield)
        ##        rtp = arcpy.RasterToPolygon_conversion(fr0, opj(inm, 'rtp' + sfx, 'SIMPLIFY')
        ##        dfsFC = arcpy.Dissolve_management(rtp, opj(inm, "dfs_frToPoly" + sfx, gridfield)
            log.debug("After dfsFC for " + str(sfx) + ' at ' + str(time.asctime()))
            arcpy.Delete_management(rtp)

        ## Calculate Fill Region statistics
            zstFrMinEl = ZonalStatisticsAsTable(fr0, 'value', input_dem, opj(inm, 'zst_' + minElFld + sfx), '', 'MINIMUM')
            log.debug('start loading of dfs at ' + time.asctime())
            valueDict1 = {r[0]:(r[1:]) for r in arcpy.da.SearchCursor(zstFrMinEl, ['VALUE', 'MIN'])}
            df.tryAddField(dfsFC, frFld, 'LONG')
            df.tryAddField(dfsFC, fillLvlFld, 'LONG')
            df.tryAddField(dfsFC, minElFld, 'LONG')
            df.tryAddField(dfsFC, 'dfs_area', 'DOUBLE')
            df.tryAddField(dfsFC, 'compactness', 'DOUBLE')
            
            with arcpy.da.UpdateCursor(dfsFC, [frFld, fillLvlFld, minElFld, gridfield, 'dfs_area', "SHAPE@AREA", 'compactness', "SHAPE@LENGTH"]) as ucur:
                for urow in ucur:
                    urow[0] = urow[3]
                    urow[1] = index
                    urow[2] = valueDict1[urow[0]][0]
                    urow[4] = urow[5]
                    urow[6] = 4 * 3.14159 * urow[5] / (urow[7]**2)
                    ucur.updateRow(urow)

            df.joinDict(wsPolys, frFld, dfsFC, frFld, [minElFld])
            log.debug('done with loading dfs at ' + time.asctime())

            ## calc max FR depth and volume
            zstFrDepth = ZonalStatisticsAsTable(fr0, 'value', rDepth, opj(inm, 'zst_' + frDepthFld + sfx), '', 'ALL')
            df.joinDict(dfsFC, frFld, zstFrDepth, 'value', ['MAX', 'SUM', 'MEAN', 'MEDIAN', 'MAJORITY'], [frDepthFld, frVolFld, 'FR0_MEAN_DPTH', 'FR0_MEDIAN', 'FR0_MAJORITY'])
            df.joinDict(wsPolys, frFld, zstFrDepth, 'value', ['MAX', 'SUM', 'MEAN', 'MEDIAN', 'MAJORITY'], [frDepthFld, frVolFld, 'FR0_MEAN_DPTH', 'FR0_MEDIAN', 'FR0_MAJORITY'])
            df.joinDict(wsPolys, frFld, dfsFC, frFld, ['dfs_area'])
            df.copytbl(verbose, zstFrDepth, sgdb)

            wsDpthFld = 'WS_MEAN_DPTH'
            dfsWsFractionFld = 'DFS_WS_FRACT'
            df.tryAddField(wsPolys, wsDpthFld, 'DOUBLE')
            df.tryAddField(wsPolys, dfsWsFractionFld, 'DOUBLE')
            with arcpy.da.UpdateCursor(wsPolys, [wsDpthFld, frVolFld, 'SHAPE@AREA', dfsWsFractionFld, 'dfs_area']) as ucur:
                for urow in ucur:
                    volume = urow[1] * float(ProcSize)**2
                    ws_area = urow[2]
                    urow[0] = volume/ws_area
                    urow[3] = urow[4]/ws_area
                    ucur.updateRow(urow)

            # nonZeroDepth = Con(rDepth > 0, rDepth)
            # zstFrNzDepth = ZonalStatisticsAsTable(fr0, 'value', nonZeroDepth, opj(inm, 'zst_nz_' + frDepthFld + sfx), '', 'ALL')
            # df.joinDict(dfsFC, frFld, zstFrNzDepth, 'value', ['MEAN', 'MEDIAN', 'MAJORITY'], ['FR_NZ_MEAN', 'FR_NZ_MEDIAN', 'FR_NZ_MAJORITY'])
            # df.copytbl(verbose, zstFrNzDepth, sgdb)

            zstFrSlope = ZonalStatisticsAsTable(fr0, 'value', slopePct, opj(inm, 'zst_fr_slope' + sfx), '', 'ALL')
            df.joinDict(dfsFC, frFld, zstFrSlope, 'value', ['MAX', 'MEAN'], [frMaxSlopeFld, frMeanSlopeFld])
            log.debug('done with adding stats to dfs at ' + time.asctime())
            dfsList.append(dfsFC)

            df.copyfc(verbose, dfsFC, sgdb)
            log.debug('done with dfs copy if needed')



            ## Get things ready for next iteration
            inDEM = demWithHoles

            rDepth = rDepthNext
            fillLvl = fillLvlNext
            maxDepth = maxDepthNext
        

            index += 1
            sfx = '_' + str(index)



        ## End additional stuff
        ## copy the fill region database
        arcpy.env.workspace = sgdb#inm
        log.debug(f"dfsList: {dfsList}")
        merged = arcpy.Merge_management(dfsList, opj(inm, 'merged'))

        depressionsCopy = arcpy.CopyFeatures_management(merged, depressions_fc)

        merged_ws = arcpy.Merge_management(wsList, opj(inm, 'merged_ws'))
        ws_output = opj(os.path.dirname(str(depressionsCopy)), 'ws_fr0_' + huc12)
        wsCopy = arcpy.CopyFeatures_management(merged_ws, ws_output)

        # fr0_merged = arcpy.Merge_management(fr0List, opj(inm, 'fr0_merged'))
        fr0_output = opj("M:\\DEP\\LiDAR_Current\\fr0", huc8, 'fr0_' + huc12 + '.tif')
        fr0_stack_string = ";".join([str(f) for f in fr0List])
        log.info(f"Stacking {fr0_stack_string} to {fr0_output}")
        df.create_needed_dirs_and_gdbs(fr0_output, log)
        fr0_stacked = arcpy.CompositeBands_management(fr0_stack_string, fr0_output)

        dontDrainWhere = "WS_MEAN_DPTH >= 30 AND DFS_WS_FRACT >= 0.7 AND dfs_area >= 100"
        dontDrainWhere2 = "Shape_Area < 2000 AND Shape_Area > 200 AND FR0_MEAN_DPTH >= 50 AND FR0_MEDIAN > FR0_MEAN_DPTH"
        dontDrain = arcpy.Select_analysis(merged_ws, opj(inm, 'dontDrain'), dontDrainWhere)
        df.copyfc(verbose, dontDrain, sgdb)
        dontDrainList = [row[0] for row in arcpy.da.SearchCursor(dontDrain, [frFld])]
        # ptr = arcpy.PolygonToRaster_management(dontDrain, frFld, opj(inm, 'dontDrainRaster'), 'MAXIMUM', '', ProcSize)

        ##--------------------------------------------------------
        ## Create the punchedDEM (drains to all sinks greater than min depth or min area)
        punchedDEM = Con(IsNull(demWithHoles) == 1, input_dem, demWithHoles)

        allow_internal_drains = True
        if allow_internal_drains:
            # ## generate a raster of all hole locations for easy recreation of punched DEM
            arcpy.env.workspace = os.path.dirname(str(holes2Punch))
            holeRasters = arcpy.ListRasters('snkunq_*')
            punchSinks = CellStatistics(holeRasters, 'MINIMUM')
            arcpy.env.workspace = sgdb#inm
            # punchSinks.save(holesTif)
            places2Drain = df.conByList(dontDrainList, "VALUE", punchSinks, procDir)
            isn = IsNull(places2Drain)
            drainedPunchedDEM = Con(isn, punchedDEM, '', 'VALUE = 1')
            punchedDEM = drainedPunchedDEM

        punchedDEM.save(output_dem)

    except:
        # Get the traceback object
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]

        # Concatenate information together concerning the error into a message string
        pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
        # Return python error messages for use in script tool or Python Window
        arcpy.AddError(pymsg)
        # Print Python error messages for use in Python / Python Window
        log.warning(pymsg + "\n")

        if arcpy.GetMessages(2) not in pymsg:
            msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
            arcpy.AddError(msgs)
            log.warning(msgs)

        sys.exit(1)

    finally:
        
        if 'logName' in locals():
            log.info("Ending script execution at " + time.asctime())
            log.info("Script execution lasted " + str(time.time()-startTime) + " seconds or " + str((time.time()-startTime)/60) + " minutes\n")




class msgStub:
    def addMessage(self,text):
        arcpy.AddMessage(text)
    def addErrorMessage(self,text):
        arcpy.AddErrorMessage(text)
    def addWarningMessage(self,text):
        arcpy.AddWarningMessage(text)

##----------------------------------------------------------------------
## below should be commented out when using as a Python Toolbox (.pyt) - in 2025, .pyt cannot handle running code in the main block
## remove the comments below for use from the windows command line

# if __name__ == "__main__":

#     if len(sys.argv) == 1:
#         arcpy.AddMessage("Whoo, hoo! Running from Python Window!")
#         cleanup = False

#         parameters = ["C:/Program Files/ArcGIS/Pro/bin/Python/envs/arcgispro-py3/pythonw.exe",
# 	"C:/DEP/Scripts/basics/cmd_puncher_DEM.pyt",
# 	"//10.27.15.155/M$/DEP/LiDAR_Current/elev_FLib_mean18/11030011/ef_3m_110300110104.tif",
# 	"//10.27.15.155/M$/DEP/LiDAR_Current/elev_PLib_mean18/11030011/ep_3m_110300110104.tif",
# 	"//10.27.15.155/M$/DEP/toolMetadata/PLib_DEMs2022_mTemplate.xml",
# 	"//10.27.15.155/D$/DEP/Man_Data_ACPF/dep_ACPF2022/11030011/idepACPF110300110104.gdb/dprsns_mean18_dem2013_3m_110300110104",
# 	"9.0",
# 	"500",
# 	"D:/DEP_Proc/DEMProc/Cut_dem2013_3m_110300110104"]
#   #  "False"]
#         for i in parameters[2:]:
#             sys.argv.append(i)
#     else:
#         #For use via Windows Command Line
#         #above 'parameters' come in via command line arguments, nothing else needed
#         arcpy.AddMessage("Whoo, hoo! Command-line enabled!")
#         # DO NOT clean up the folder after done processing - matcher needs this data
#         cleanup = False

#     # switch a text 'True' into a real Python True
#     # make sure puncher cleanup is False until we log the outputs elsewhere so they aren't deleted
#     cleanup = False#cleanup = True if cleanup == "True" else False

#     messages = msgStub()

#     # inputs then outputs, change "" to Python None
#     (input_dem, output_dem, plib_metadata, depressions_fc, depth_threshold, 
#      area_threshold, procDir
#      ) = [i if i != "" else None for i in sys.argv[1:]]#[i for i in sys.argv[1:]]

#     doPuncher(input_dem, output_dem, plib_metadata, depressions_fc, depth_threshold, area_threshold, procDir, cleanup, messages)
#     arcpy.AddMessage("Back from doing!")
