# -*- coding: utf-8 -*-
'''A Python program to get match upstream and downstream features from high resolution
DEMs facilitating flow improvements to the DEM. This table of match points is then used 
in later programs to cut and 'improve' the DEM.'''
# 2024.04.03 - moved to AG Pro 3.2, Python 3. Seems to still be some bugs with in_memory
#               visualization in Pro 3.2 interface. Pro saves in_memory to the project 
#               database and Select Layer sometimes bombs during this. As does JoinField

import arcpy
import sys
import os
import platform
import time
from os.path import join as opj
import traceback
from arcpy.sa import *
import pickle

sys.path.append("C:\\DEP\\Scripts\\basics")
import dem_functions as df


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
        self.label = "EPT_WESM_download"
        self.description = "Creates a feature class to enable EPT downloads"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""

        param0 = arcpy.Parameter(
            name="buffered_field_boundaries",
            displayName="Buffered ACPF Field Boundaries",
            datatype="DEFeatureClass",
            parameterType='Required',
            direction="Input")
        
        param1 = arcpy.Parameter(
            name="fill_tif",
            displayName="Pit Filled DEM",
            datatype="DERasterDataset",
            parameterType='Required',
            direction="Input")
        
        param2 = arcpy.Parameter(
            name="punch_tif",
            displayName="Punched DEM",
            datatype="DERasterDataset",
            parameterType='Required',
            direction="Input")
               
        param3 = arcpy.Parameter(
            name="void_tif",
            displayName="Void Filled DEM",
            datatype="DERasterDataset",
            parameterType='Optional',
            direction="Input")
        
        param4 = arcpy.Parameter(
            name="merged_medians",
            displayName="Merged Medians Feature Class",
            datatype="DEFeatureClass",
            parameterType='Required',
            direction="Input")
        
        param5 = arcpy.Parameter(
            displayName="Local Processing Directory",
            datatype="DEFolder",
            parameterType='Required',
            direction="Input")
                
        param6 = arcpy.Parameter(
            name="pickle_distance_file",
            displayName="Pickle Distance Output File",
            datatype="DEFile",
            parameterType='Required',
            direction="Output")
        
        param7 = arcpy.Parameter(
            displayName="Depression Punch Area Threshold",
            datatype="GPString",
            parameterType='Optional',
            direction="Input")
               
        parameters = [param0, param1, param2, param3,
                  param4, param5, param6, param7]

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
        doMatcher(parameters[0].valueAsText, parameters[1].valueAsText, parameters[2].valueAsText, parameters[3].valueAsText, parameters[4].valueAsText, parameters[5].valueAsText, parameters[6].valueAsText, parameters[7].valueAsText, cleanup, messages)
        return

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return


def setupNoMedians(upPtsCmb, medianFrFld, dfs2cut4step, minElFld, frFld, maxWsMdnFld, mdnFracFld, sgdb):
    arcpy.AddField_management(upPtsCmb, 'bfr_dist', 'DOUBLE')
    arcpy.AddField_management(upPtsCmb, medianFrFld, 'LONG')#'bfr_dist', 'DOUBLE')
    arcpy.CalculateField_management(upPtsCmb, medianFrFld, '0', 'PYTHON')
    try:
        maxWsMedian = arcpy.Statistics_analysis(dfs2cut4step, 'ws_max_median', [[minElFld, 'COUNT']], frFld)
    except:
        maxWsMedian = arcpy.Statistics_analysis(dfs2cut4step, opj(sgdb, 'ws_max_median'), [[minElFld, 'COUNT']], frFld)
    arcpy.AddField_management(maxWsMedian, maxWsMdnFld, 'DOUBLE')
    arcpy.AddField_management(maxWsMedian, mdnFracFld, 'DOUBLE')
    return maxWsMedian


def selectFrAndWsDpPts(gdbTable, dnPts, frFld, sfx, log, inm = 'in_memory'):
    dnFrWsList = []
    with arcpy.da.SearchCursor(gdbTable, ['FILL_RGN', 'ws_lvl'], sql_clause = (None, 'GROUP BY ' + frFld + ', ws_lvl')) as scur:
        for srow in scur:
            dnFrWsList.append(srow)
    log.debug("list Lenth is " + str(len(dnFrWsList)))

    stringLength = len(dnFrWsList) * (len(frFld) + len(' = AND = ') + len('ws_lvl ') + len(str(dnFrWsList[-1])) + len(' OR '))
    log.debug("stringLenth is " + str(stringLength))
    
##    str_threshold = 150000 # Arc Map 10.3
    str_threshold = 65000 # AG Pro
    if stringLength > str_threshold:
        iterNeeded = stringLength//str_threshold + 1
        appendList = []
        sliceDenom = int(len(dnFrWsList)/iterNeeded)
        for i in range(iterNeeded):
            ofSfx = '_' + str(i)
            sublist = dnFrWsList[i*sliceDenom:(i+1)*sliceDenom]
            if i == iterNeeded - 1:
                sublist = dnFrWsList[i*sliceDenom:]
            
##            gcSel4step = buildSelection(sublist, in_field)
            sel = ''
            for index, item in enumerate(sublist):
                if index == 0:
                    sel = frFld + ' = ' + str(int(item[0])) + ' AND ws_lvl = ' + str(int(item[1]))
                else:
                    sel += ' OR ' + frFld + ' = ' + str(int(item[0])) + ' AND ws_lvl = ' + str(int(item[1]))
                    
            if i == 1:
                selInit = arcpy.Select_analysis(dnPts, opj(inm, 'dn_pts_bst2' + sfx), sel)
            else:
                selLater = arcpy.Select_analysis(dnPts, opj(inm, 'dn_pts_bst2' + sfx + ofSfx), sel)
                appendList.append(selLater)
        selOut = arcpy.Append_management(appendList, selInit)
    else:
        sel = ''
        for index, item in enumerate(dnFrWsList):
            if index == 0:
                sel = frFld + ' = ' + str(int(item[0])) + ' AND ws_lvl = ' + str(int(item[1]))
            else:
                sel += ' OR ' + frFld + ' = ' + str(int(item[0])) + ' AND ws_lvl = ' + str(int(item[1]))
##        gcSel4step = buildSelection(in_list, in_field)
        selOut = arcpy.Select_analysis(dnPts, opj(inm, 'dn_pts_bst2' + sfx), sel)
    return selOut



def setAllSearchDistances(inFC, wsSearchDistFld, minFrDistFld, frFld, minElFld, ofElFld, srchMult, crestMin, maxWsMedian, maxWsMdnFld, wsMdnSearchDistFld, pntMdnSearchDistFld, maxExtraSearchDist, mdnFracFld, frDepthFld):
    try:
        ## bring in the max median buffer for each fill region
        arcpy.JoinField_management(inFC, frFld, maxWsMedian, frFld, [maxWsMdnFld, mdnFracFld])
        df.tryAddField(inFC, wsSearchDistFld, "DOUBLE")
        df.tryAddField(inFC, wsMdnSearchDistFld, "DOUBLE")
        df.tryAddField(inFC, pntMdnSearchDistFld, "DOUBLE")

        maxSearchDist = 0
        pntMdnMaxSearchDist = 0
        baseErrorList = []
        maxErrorList = []
        ## Calculate base search distance for all fill regions
        ## Calculate new buffer distance for channelzed areas and areas where fill region is steep
        with arcpy.da.UpdateCursor(inFC, [wsSearchDistFld, minFrDistFld, frFld, minElFld, ofElFld, maxWsMdnFld, wsMdnSearchDistFld, pntMdnSearchDistFld, mdnFracFld, frDepthFld]) as ucur:#'MIN_NEAR_DIST']) as ucur:
            for row in ucur:
                crestWidth = 1.41*crestMin/2.0
                frDepth = row[-1]#row[4] - row[3]
                frDepth = min(frDepth, 2000)
            # set minimum and maximum search distances based on backslope (%)
                maxBackslope = 0.50#fraction
                if frDepth < 150:
                    minBackslope = 0.05#fraction
                else:
                    minBackslope = 0.10#fraction
                minFrSearchDist = crestWidth + (frDepth)*0.01*(1/maxBackslope)*srchMult
                maxFrSearchDist = crestWidth + (frDepth)*0.01*(1/minBackslope)*srchMult
            # set base distance is srchMult * minDist3 * 1.25 (minDist3 could be ~80% depth, so extrapolate to 100% depth)
                # fix all None by setting equal to else
                if row[ucur.fields.index(maxWsMdnFld)] is not None:
                    if row[1] == None:
                        baseDist = minFrSearchDist
                        baseErrorList.append(str(row[2]))
                        ## errors seem to be due to no minFrDist data
        ##row is [None, None, 4572, 25795, 25800, 0.0, None, None, 0.0]
                    ## if upstream points are farther to median than fill region distance, use greater distance #why?
                    elif row[ucur.fields.index(maxWsMdnFld)] > row[1]:
                        actFrSearchDist = crestWidth + 1.25*row[ucur.fields.index(maxWsMdnFld)]*srchMult
                        baseDistPre = min(maxFrSearchDist, actFrSearchDist)
                        baseDist = max(minFrSearchDist, baseDistPre)
                    else:
                        actFrSearchDist = crestWidth + 1.25*row[1]*srchMult
                        baseDistPre = min(maxFrSearchDist, actFrSearchDist)
                        baseDist = max(minFrSearchDist, baseDistPre)
                else:
                    actFrSearchDist = crestWidth + 1.25*row[1]*srchMult
                    baseDistPre = min(maxFrSearchDist, actFrSearchDist)
                    baseDist = max(minFrSearchDist, baseDistPre)
                    
                row[ucur.fields.index(wsSearchDistFld)] = baseDist

            # calculate ws median search distance, if crossing median add in highway width`
                ## fr with all points in median only get half the median width
                if row[ucur.fields.index(maxWsMdnFld)] is not None:
                    if row[ucur.fields.index(maxWsMdnFld)] > 0 and row[ucur.fields.index(mdnFracFld)] == 1.0:
                        wsTotalDist = baseDist + 0.5*row[ucur.fields.index(maxWsMdnFld)]
                    elif row[ucur.fields.index(maxWsMdnFld)] > 0:
                        wsTotalDist = baseDist + row[ucur.fields.index(maxWsMdnFld)]
                        # AG Pro 3.2 addition - set mdnFracFld to zero
                        row[ucur.fields.index(mdnFracFld)] = 0.0
                    else:
                        wsTotalDist = baseDist
                        ## set nulls to zero
                        row[ucur.fields.index(maxWsMdnFld)] = 0.0
                        # AG Pro 3.2 addition - set mdnFracFld to zero
                        row[ucur.fields.index(mdnFracFld)] = 0.0
                else:
                    wsTotalDist = baseDist
                    ## set nulls to zero
                    row[ucur.fields.index(maxWsMdnFld)] = 0.0
                    # AG Pro 3.2 addition - set mdnFracFld to zero
                    row[ucur.fields.index(mdnFracFld)] = 0.0

                row[ucur.fields.index(wsMdnSearchDistFld)] = wsTotalDist

            # calculate point median search distance, if crossing median add in highway width
                ## fr with all points in median only get half the median width
                if row[ucur.fields.index(maxWsMdnFld)] is not None:
                    if row[ucur.fields.index(maxWsMdnFld)] > 0 and row[ucur.fields.index(mdnFracFld)] == 1.0:
                        pntTotalDist = baseDist*2 + 0.5*row[ucur.fields.index(maxWsMdnFld)]
                    elif row[ucur.fields.index(maxWsMdnFld)] > 0:
                        pntTotalDist = baseDist*2 + row[ucur.fields.index(maxWsMdnFld)]
                    else:
                        pntTotalDist = baseDist*2
                else:
                    pntTotalDist = baseDist*2

                if pntTotalDist >= maxExtraSearchDist:
                    pntTotalDist= maxExtraSearchDist
    ##                log.warn('\nMaximum Extra Search Distance Exceeded!\nFill region was ' + str(row[ucur.fields.index(frFld)]) + '\n')
                    maxErrorList.append(row[ucur.fields.index(frFld)])

                row[ucur.fields.index(pntMdnSearchDistFld)] = pntTotalDist

                ucur.updateRow(row)

                pntMdnMaxSearchDist = max(pntMdnMaxSearchDist, pntTotalDist)

    ##    log.debug('baseDist ERROR for ' + str(baseErrorList))
    ##    log.debug('maxDist ERROR for ' + str(maxErrorList))
    ##
    ##    log.debug("Done creating buffer distance at " + time.asctime())
        return pntMdnMaxSearchDist

    except Exception as e:
##        log.debug(e.message)
##        arcpy.AddError(e.message)
        print('row is ' + str(row))

        # Get the traceback object
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]

        # Concatenate information together concerning the error into a message string
        pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
        print(pymsg)
##        log.warn(pymsg)




def idMedianPtsUp(upPtsAll, medianFrFld, mdnsOnly):
####------------------------------------------------------------------------------
## Define transportation areas for additional search - roads, railroads, and runways
    arcpy.AddField_management(upPtsAll, medianFrFld, 'LONG')
    arcpy.CalculateField_management(upPtsAll, medianFrFld, '0', 'PYTHON')
## Indicate FRs that are in the median to reduce search distance
    nrBstUp3Lyr = arcpy.MakeFeatureLayer_management(upPtsAll, 'nrBstUpPts3Lyr')
    if int(arcpy.GetCount_management(mdnsOnly).getOutput(0)) > 0:
        arcpy.SelectLayerByLocation_management(nrBstUp3Lyr, "WITHIN", mdnsOnly)

        arcpy.CalculateField_management(nrBstUp3Lyr, medianFrFld, '1', 'PYTHON')


def idMedianPtsDn(dnPtsAll, medianFrFld, mdnsOnly):
####------------------------------------------------------------------------------
## Define transportation areas for additional search - roads, railroads, and runways
    arcpy.AddField_management(dnPtsAll, medianFrFld, 'LONG')
    arcpy.CalculateField_management(dnPtsAll, medianFrFld, '0', 'PYTHON')
## Indicate FRs that are in the median to reduce search distance
    nrBstDn3Lyr = arcpy.MakeFeatureLayer_management(dnPtsAll, 'nrBstDnPts3Lyr')
    if int(arcpy.GetCount_management(mdnsOnly).getOutput(0)) > 0:
        arcpy.SelectLayerByLocation_management(nrBstDn3Lyr, "WITHIN", mdnsOnly)

        arcpy.CalculateField_management(nrBstDn3Lyr, medianFrFld, '2', 'PYTHON')


			
def condenseStats(inZone, inValue, inm, outTable, statType, sfx):
    if arcpy.Exists(opj(inm, outTable)) == False:
        zst = ZonalStatisticsAsTable(inZone, "VALUE", inValue, opj(inm, outTable), "", statType)
    else:
        copy_rows= arcpy.CopyRows_management(opj(inm, outTable), opj(inm, 'temp_fr0_stats'))
        zst1 = ZonalStatisticsAsTable(inZone, "VALUE", inValue, opj(inm, outTable + sfx), "", statType)
        zst = arcpy.Append_management([copy_rows, zst1], opj(inm, outTable), "NO_TEST")
        arcpy.Delete_management(zst1)

    return zst


def conByList(in_list, in_field, in_rast):#, cp):
    initWs = arcpy.env.workspace
    arcpy.env.workspace = os.path.dirname(in_rast)

    list_length = len(in_list)
        
    if list_length == 1:
        gcSel4step = df.buildSelection(in_list, in_field)
        sel_out = Con(in_rast, in_rast, '', gcSel4step)
    elif list_length > 14999:
        iter_needed = list_length//9999 + 1
        append_list = []
        slice_denom = int(len(in_list)/iter_needed)
        for i in range(iter_needed):
            ofSfx = '_' + str(i)
            sub_list = in_list[i*slice_denom:(i+1)*slice_denom]
            if i == iter_needed - 1:
                sub_list = in_list[i*slice_denom:]

            gc_sel4step = in_field + " IN " + str(tuple(sub_list))
            con4later = Con(in_rast, in_rast, '', gc_sel4step)
            con4later.save('c4l_' + str(i))
            append_list.append(con4later.name)

        sel_out_float = CellStatistics(append_list)
        sel_out = Int(sel_out_float)
    else:
        gc_sel4step = in_field + " IN " + str(tuple(in_list))
        sel_out = Con(in_rast, in_rast, '', gc_sel4step)

    arcpy.env.workspace = initWs
    return sel_out



def doMatcher(fill_or_void_tif, punch_tif, buffered_fc, merged_medians, fr0_rasters, ws_polys, dfs_polys, proc_dir, search_distance_file, match_depth, cleanup, messages):

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
    ##updated to Python 3.9 for AG Pro 3.2 - 2024.01.05
    # added dfs2CutList to store dfs2Cut as well as dfs feature classes 2019/04/26
    # 2019.09.10 - added try/except to wsLineBfrUnion to handle Invalide Topology error
    # 2019.09.13 - added arcpy.Exists(mergedMdnsHuc8) to remove problems when no merged medians in the HUC8
    # 2020.03.10 - fixed invalid topology error in area2search HUC12 090201060106 - bkgelder
    # 2020.03.11 - set extent to DEM to avoid block statistics out of memory errors HUC12 070102050704 - bkgelder
    # 2021.02.02 - updated mergedMdnsHuc8 to mergedMdnsHuc8FC, now load paths via dictionary - bkgelder
    # 2022.03.22 - many updates to run on Python 3.7 with minimal functionality change
    # 2022.03.29 - added SplitRaster to split raster in areas where there were too many combination for Combine
    # 2022.05.25 - fixed issue from above testing where frRaster did not advance past the first one
    # 2022.06.06 - fixed of above issue revealed that split raster was not correctly set up, failing when trying to merge down cells from all iterations
    # 2022.06.08 - above solution required adding a copy of raster data to build statistics for Combining
    # 2024.01.05 - added multiple try/except statements to handle upgrade to Python 3.9/AG Pro 3.2, mostly using join by dictionary instead of JoinField
    # 2024.01.10 - moved selectFrAndWsDpPts into script and lowered threshold for splitting into subselections due to change in AG Pro

    try:
        arguments = [fill_or_void_tif, punch_tif, buffered_fc, merged_medians, fr0_rasters, ws_polys, dfs_polys, proc_dir, search_distance_file, match_depth, cleanup]

        for a in arguments:
            if a == arguments[0]:
                arg_str = str(a) + '\n'
            else:
                arg_str += str(a) + '\n'

        messages.addMessage("Tool: Executing with parameters:\n" + arg_str)

        huc12, huc8 = df.figureItOut(fill_or_void_tif)

        if cleanup:
            # log to file only
            log, nowYmd, logName, startTime = df.setupLoggingNoCh(platform.node(), sys.argv[0], huc12)
            verbose = False
            arcpy.SetLogHistory = False
        else:
            # log to file and console
            log, nowYmd, logName, startTime = df.setupLoggingNew(platform.node(), sys.argv[0], huc12)
            verbose = True
            arcpy.SetLogHistory = True

        startTime = time.time()
        log.info("Beginning execution: " + time.asctime())
        messages.addMessage("Log file at " + logName)

        log.warning('log file is ' + logName)

        log.info("Tool: Executing with parameters:\n" + arg_str)

    # Create output directories
    ## Set the environments
        arcpy.env.scratchWorkspace = proc_dir

        sgdb = arcpy.env.scratchGDB
        sfldr = arcpy.env.scratchFolder
        arcpy.env.scratchWorkspace = sgdb
        arcpy.env.workspace = sgdb

        inm = 'in_memory'

        arcpy.env.snapRaster = fill_or_void_tif

        arcpy.env.extent = fill_or_void_tif

        arcpy.env.cellSize = int(arcpy.Raster(fill_or_void_tif).meanCellHeight)
        ProcSize = arcpy.env.cellSize

        arcpy.CheckOutExtension("Spatial")
        arcpy.env.overwriteOutput = True

    # set up output raster names

        ## factor to increase search distance by
        srchMult = 3

        ofElList = ["FR_OF_EL", "LONG"]
        ofElFld = ofElList[0]
        frFld = 'FILL_RGN'
        fillLvlFld = 'FILL_LVL'

        minElFld = 'FR_MIN_EL'
        cutElFld = 'FR_CUT_EL'
        maxFillFld = 'MAX_FILL'

        wsSearchDistFld = 'WS_SRCH_DST'

        minFrDistList3 = ['MIN_WS_DST3', 'DOUBLE']
        minFrDistFld = minFrDistList3[0]

        medianFrFld = "Median_FR"        

        frAllCrvFrac = ['ALL_FR_CRV_FRAC', 'DOUBLE']
        frThknsFld = 'FR_THKNS'

        frDepthFld = 'FR_DEPTH'

        minElDifFld = 'MIN_EL_DIF'
        minElDif = [minElDifFld, 'long']

        maxElDifFld = 'MAX_EL_DIF'
        maxElDif = [maxElDifFld, 'long']

        bsNbr = NbrRectangle(5, 5, 'CELL')
        
        minSummaryFld = 'MIN_MIN'
        revCutElFld = 'REV_CUT_EL'

        bsMinDnElFld = 'bs_min_el_dn'

        allDnCellsFld = 'all_dn'
        wsLvlFld = 'ws_lvl'

        crestMin = 9.0 # Meters

        gridfield = 'gridcode'
        gridfield2 = 'grid_code'

        frMaxSlopeFld = 'FR_MAX_SLP'

        ## wsSearchDistFld is an adjusted search distance used distance from deep point to ws boundary and a minimum width
        ## it is not adjusted for crossing a median, just used for figuring out if upstream points are close
        ## wsMdnSearchDistField is ws search dist + max fr median search dist
        ## Used for buffered watershed searches
        ## pntSearchDistField is inner and outer ws search dist (i.e. x2)
        ## pntMdnSearchDistField is inner and outer ws search dist (i.e. x2) + median search
        ## Used for buffered point searches
        maxWsMdnFld = 'MAX_bfr_dist'
        pntMdnSearchDistFld = 'PNT_MDN_DIST'
        wsMdnSearchDistFld = 'WS_MDN_DIST'
        mdnFracFld = 'up_mdn_frac'

        ## cutOneMinCrit is ratio of FR depth + min FR elevaiton that external elevation must be less than to be considered for cutting
        cutOneMinCrit = 0.5#0.4

        maxSearchDistList = []

        dfsList = []
        dfs2CutList = []
        srchList = []
        upPtsList = []
        dnPtsList = []

        # punchGdb = os.path.join(sgdb, os.path.splitext(os.path.basename(punch_tif))[0])
        # if not arcpy.Exists(punchGdb):
        #     punchGdbResult = arcpy.CopyRaster_management(punch_tif, punchGdb)
        # punchedDEMNoHoles = Con(IsNull(punchGdb) == 1, fill_or_void_tif, punch_tif)
        # bestestDEM = punchedDEMNoHoles#Raster(punchGdb)

        slopePct = Raster(opj(proc_dir, 'slope_pct'))
        flats2 = Con(slopePct == 0.0, 1, 0)
        noInteriorFlatsDEM = Con(flats2 == 0, fill_or_void_tif, '')

        fsPfdMdn = FocalStatistics(fill_or_void_tif, "RECTANGLE 3 3 CELL", "MEDIAN")#pitFilledDEM

        arcpy.env.workspace = inm

        match_stats = arcpy.CreateTable_management(inm, 'match_stats')
        df.tryAddField(match_stats, 'HUC12', 'TEXT')
        df.tryAddField(match_stats, 'LEVEL', 'SHORT')
        df.tryAddField(match_stats, 'UP_PTS_COUNT', 'LONG')
        df.tryAddField(match_stats, 'DN_PTS_COUNT', 'LONG')
        df.tryAddField(match_stats, 'UP_PTS_RATIO', 'DOUBLE')
        df.tryAddField(match_stats, 'DN_PTS_RATIO', 'DOUBLE')

        fr0_dir = os.path.dirname(fr0_rasters)
        arcpy.env.workspace = fr0_dir
        fr0_name = os.path.basename(fr0_rasters)
        frRasters = arcpy.ListRasters(fr0_name[:4] + '*')#'fr0_*')
        frRasters.sort()
        for i, item in enumerate(frRasters[:]):
            log.info('working on frRaster: ' + str(item))
            arcpy.env.workspace = proc_dir
            sfx = '_' + item.split('_')[-1]
            dfsFC = dfs_polys[:-2] + sfx
            # dfsFC = os.path.join(sgdb, 'dfs_frToPoly' + sfx)
            dfsList.append(dfsFC)
            # add identifier for medians
            df.tryAddField(dfsFC, medianFrFld, "SHORT")
            # wsPolys = os.path.join(sgdb, 'ws_polys' + sfx)
            wsPolys = ws_polys[:-2] + sfx

    ##                selFull = '(' + ofElFld + ' - ' + minElFld + ') > 67 OR ' + maxFrOfDistFld + ' >= ' + str(3 * ProcSize)
    ##        selFull = '(' + ofElFld + ' - ' + minElFld + ') > ' + str(3.0 * RMSE) + ' AND (' + ofElFld + ' - ' + minElFld + ')*0.01 / ' + minFrDistFld + ' >= 0.05'
            selFull = frDepthFld + ' > ' + str(2.0 * float(match_depth)) + ' AND ' + frMaxSlopeFld + ' > 5.0'
            ## instead of frPctDrop could analyze by water 'piling' up near edge

    ####        df.condDelete(verbose, wsPolys)
            dfs2cut4step = arcpy.Select_analysis(dfsFC, opj(inm, 'dfs2cut' + sfx), selFull)
            frs2cut = df.getFrsAsList(dfs2cut4step, frFld, '')
            if len(frs2cut) > 0:

                ws2cut4step = df.selectByList(frs2cut, gridfield, wsPolys, 'ws_2_cut', sfx, inm)

                if int(arcpy.GetCount_management(ws2cut4step).getOutput(0)) > 0:
                    df.copyfc(verbose, ws2cut4step, sgdb)
                    df.copyfc(verbose, dfs2cut4step, sgdb)
                ## Create list of fill regions to select upstream points for appropriate fill regions

                fr0 = conByList(frs2cut, 'VALUE', opj(proc_dir,'fr0' + sfx))
                wsLvl = Raster(wsLvlFld + sfx)
    ####            fdPrev = Raster('fd_lvl0' + sfx)
                # inDEM = Raster('nwdm' + sfx)
                # snkUnique = Raster('snkunq' + sfx)
                fillPrevLvl = Raster(fillLvlFld.lower() + sfx)

            ## Calculate maximum thickness in each fill region (i.e. how far off can the flowpath be?) for later processing
            ## Create cost surface - fill regions are 1, others 0
                if sys.version_info.major == 2:
                    fillRgnThickCost = Con(IsNull(fr0), 0, 1)
                ## Create cost distance raster of thickness
                    thknsCostDistance = CostDistance(ExtractByAttributes(fillRgnThickCost, 'VALUE = 0'), fillRgnThickCost)
                else:
                    # fillRgnThick = Con(IsNull(fr0), 0, 1)
                    fillRgnThickCost = Con(IsNull(fr0), 0.001, 1)
                    thknsCostDistance = CostDistance(ExtractByAttributes(fillRgnThickCost, 'VALUE = ' + str(fillRgnThickCost.minimum)), fillRgnThickCost)
                ##            thknsCostDistance.save(gdb + "thkns_cd" + sfx)
                
            ## Calculate maximum thickness in each fill region for later processing
                flrgThkns = condenseStats(fr0, thknsCostDistance, inm, "flrg_thkns", "ALL", sfx)
                df.addCalcJoin(dfs2cut4step, frFld, flrgThkns, 'value', [frThknsFld, 'DOUBLE'], '!MAX!')

    ## Fix ND road intersections by cutting up and down, then fix ND FR exits that are + curvature by replacing + curvature areas

            ## Calculate the overflow flowpath and calculate distance from edge of fill region0 to path
                ## We can use this distance to quantitatively evaluate if overflow compomises accuracy

                zstFrOfEl = ZonalStatisticsAsTable(fr0, 'value', fillPrevLvl, opj(inm, 'zst_' + ofElFld + sfx), '', 'MINIMUM')
                df.addCalcJoin(dfs2cut4step, frFld, zstFrOfEl, 'value', [ofElFld, 'LONG'], '!MIN!')
                zstMaxFilEl = ZonalStatisticsAsTable(fr0, 'value', opj(proc_dir, 'fill_lvl_0'), opj(inm, 'zst_' + maxFillFld + sfx), '', 'MINIMUM')
                df.addCalcJoin(dfs2cut4step, frFld, zstMaxFilEl, 'value', [maxFillFld, 'LONG'], '!MIN!')

                wsOfEl = ZonalStatistics(wsLvl, 'value', fillPrevLvl, 'MINIMUM')
                wsOfEl.save(opj(proc_dir, ofElFld))
                log.debug('done with ofElFld at ' + time.asctime())

            ## Define the quartermost deepest parts of the watershed (deeper or equal to 1/4 fr depth)
                wsMinEl = ZonalStatistics(wsLvl, 'value', fill_or_void_tif, 'MINIMUM')#pitFilledDEM, 'MINIMUM')
                qrtrDeepEl = wsMinEl + 0.25*(wsOfEl-wsMinEl)
                quarterDeepFr0 = Con(fill_or_void_tif <= qrtrDeepEl, fr0)
                quarterDeepFr1 = Con(quarterDeepFr0, 1)

            ## Define the deepest parts of the watershed (deeper or equal to fr median minimum elevation)
                minMdnEl = ZonalStatistics(fr0, 'value', fsPfdMdn, 'MINIMUM')
                cutEl = Con(minMdnEl < qrtrDeepEl, minMdnEl, qrtrDeepEl)
                cutElFr0 = Con(fill_or_void_tif <= cutEl, fr0)
                cutElFr1 = Con(cutElFr0, 1)

                zstCutElAlt = ZonalStatisticsAsTable(fr0, 'value', cutEl, opj(inm, 'zst_alt_' + minFrDistFld + sfx), '', 'MINIMUM')
        ####                df.copyTbl(verbose, zstCutElAlt, sgdb)
                df.addCalcJoin(dfs2cut4step, frFld, zstCutElAlt, 'value', [cutElFld, 'LONG'], '!MIN!')# + cutElFld + '!')'ALT_CUT_EL'

            ## Calculate maximum thickness in each watershed (i.e. distance to ws boundary) for later processing
                if sys.version_info.major == 2:
                ## Create cost surface - fill regions are 1, others 0
                    wsThickCost = Con(wsLvl, 1)
                else:
                    wsThickCost = Con(wsLvl, 1, 0.001)
        ##                        wsThickCost.save("ws_thk_cst" + sfx)

            ## Create cost distance raster of thickness from ws boundaries (range > 0)
                wsBnd = Con(FocalStatistics(wsLvl, "RECTANGLE 3 3 CELL", "RANGE") > 0, wsLvl, '')
                maxWsBndThick = 120
                wsCostDistance100 = CostDistance(wsBnd, wsThickCost, maxWsBndThick)

            ## Set null values to maximum
            ## Add 1.5 to all distances to remove 0 values
                wsCostDistance = Con(IsNull(wsCostDistance100), maxWsBndThick, wsCostDistance100) + 1.5
                wsCostDistance.save(opj(proc_dir, "ws_thkns_cd" + sfx))

            ## For cut el
            ## Preliminary calculations to determiine how many low points near each other
                fstCutFr1Sum = FocalStatistics(cutElFr1, "RECTANGLE 3 3 CELL", "SUM")

                zsCutFr1SumMax = ZonalStatistics(cutElFr0, 'value', fstCutFr1Sum, 'MAXIMUM')

                cdCutEl = Con(cutElFr0, wsCostDistance)
                fstSum3CutEl = Con(fstCutFr1Sum >=3, cdCutEl)

                cutElWsCD = Con(zsCutFr1SumMax >= 3, fstSum3CutEl, cdCutEl + ProcSize)
                wsCDinCutEl = ZonalStatistics(fr0, 'value', cutElWsCD, 'MINIMUM')

            ## For qrtr el
            ## Preliminary calculations to determiine how many low points near each other
                fstQrtrFr1Sum = FocalStatistics(quarterDeepFr1, "RECTANGLE 3 3 CELL", "SUM")

                zsQrtrFr1SumMax = ZonalStatistics(quarterDeepFr0, 'value', fstQrtrFr1Sum, 'MAXIMUM')

                cdQrtrEl = Con(quarterDeepFr0, wsCostDistance)
                fstSum3QrtrEl = Con(fstQrtrFr1Sum >=3, cdQrtrEl)

                qrtrElWsCD = Con(zsQrtrFr1SumMax >= 3, fstSum3QrtrEl, cdQrtrEl + ProcSize)
                wsCDinQrtrEl = ZonalStatistics(fr0, 'value', qrtrElWsCD, 'MINIMUM')

            ## Compare ratios of minimum distances and select, assuming cut el is minimum, qrtr should be about 25% less
                cutQrtrDistRatio = wsCDinQrtrEl/wsCDinCutEl

                deepEnoughFr0 = Con(cutQrtrDistRatio < 0.75, quarterDeepFr0, cutElFr0)
                deepEnoughFr0.save(opj(proc_dir, 'dp_fr' + sfx))
                zstFrMinDist = ZonalStatisticsAsTable(fr0, 'value', Con(cutQrtrDistRatio < 0.75, wsCDinQrtrEl, wsCDinCutEl), opj(inm, 'zst_' + minFrDistFld + sfx), '', 'MINIMUM')
                df.addCalcJoin(dfs2cut4step, frFld, zstFrMinDist, 'value', [minFrDistFld, 'LONG'], '!MIN!')

                zsFrMinDist = ZonalStatistics(fr0, 'value', Con(cutQrtrDistRatio < 0.75, wsCDinQrtrEl, wsCDinCutEl), 'MINIMUM')
                nearDeepEnoughFr0 = Con(wsCostDistance < srchMult*zsFrMinDist + crestMin/2.0, deepEnoughFr0)

                deepEnoughFr1 = Con(IsNull(deepEnoughFr0), 0, 1)
                if sys.version_info.major == 2:
                    deepFrCostDistance = CostDistance(ExtractByAttributes(deepEnoughFr1, 'VALUE = 0'), deepEnoughFr1)
                else:        
                    deepEnoughFr1Cost = Con(IsNull(deepEnoughFr0), 0.001, 1)
                ## Calculate maximum thickness for deep part of each fill region (i.e. how far off can the flowpath be?)
                ## Use to jest get outer cells of deep fill region in areas of flattened water
                ## Create cost distance raster of thickness
                    deepFrCostDistance = CostDistance(ExtractByAttributes(deepEnoughFr1, 'VALUE = 0'), deepEnoughFr1Cost)
                
            ## Calculate maximum thickness in each fill region for later processing
                deepFrCDMax = ZonalStatistics(deepEnoughFr0, 'value', deepFrCostDistance, 'MAXIMUM')
                if sys.version_info.major == 2:
                    outerDeepEnoughFr0 = Con(deepFrCDMax > ProcSize * 3.0, Con(deepFrCostDistance < ProcSize *2.0, nearDeepEnoughFr0), nearDeepEnoughFr0)
                
                else:
                    ProcSize2 = ProcSize *2.0
                    test1 = LessThan(deepFrCostDistance, ProcSize2)
                    outerDeepEnoughFr0Pre = Con(test1, nearDeepEnoughFr0)
        ##            outerDeepEnoughFr0.save(opj(proc_dir, 'dp_fr' + sfx))
                    ProcSize3 = ProcSize *3.0
                    test2 = GreaterThan(deepFrCDMax, ProcSize3)
                    outerDeepEnoughFr0 = Con(test2, outerDeepEnoughFr0Pre, nearDeepEnoughFr0)
        ##            outerDeepEnoughFr0.save(opj(proc_dir, 'dp_fr' + sfx))

        ######------------------------------------------------------------------------------


        ## Create upstream points by selective filtering
            ## Just get those FRs to be cut
                deepEnoughFr0ToCut = df.conByList(frs2cut, 'VALUE', outerDeepEnoughFr0, proc_dir)
    ##            deepEnoughFr0ToCut.save('dp_fr2cut' + sfx)

    ####            in_list, in_field, in_rast = frs2cut, 'VALUE', outerDeepEnoughFr0
    ####            gc_sel4step = in_field + " IN " + str(tuple(in_list))
    ####            sel0_out = Con(fr0, 1, '', gc_sel4step)
    ####            deepEnoughFr0ToCut = Con(sel0_out, outerDeepEnoughFr0)

            ## Just get those at the minimum block elevation
                upPtsDEM = Raster(fill_or_void_tif)#ndAllFixedDEM)#bigNdFixed
                bsMinUpEl = BlockStatistics(Con(deepEnoughFr0ToCut, upPtsDEM), bsNbr, 'MINIMUM')

                # rgBsMinUpEl = RegionGroup(bsMinUpEl)

                deepEnoughFr0Bs = Con(upPtsDEM == bsMinUpEl, deepEnoughFr0ToCut)
                deepEnoughFr0Bs.save(opj(proc_dir, 'up_bs' + sfx))
                deepBsMinUpElName = os.path.basename(str(deepEnoughFr0Bs))

                upstreamCombine = Combine([deepEnoughFr0Bs, upPtsDEM, wsCostDistance])#bestestDEM, wsCostDistance])
                upstreamCombine.save(opj(proc_dir, 'up_cmb' + sfx))
                log.debug(f"upPtsCombine count is: " + str(arcpy.GetCount_management(upstreamCombine)))
                deepEnoughFr0Name = df.getfields(upstreamCombine)[3]
                upPtsDemName = df.getfields(upstreamCombine)[4]
                wsCostDistName = df.getfields(upstreamCombine)[5]

                upPtsCmb = arcpy.RasterToPoint_conversion(upstreamCombine, opj(inm, 'cb_up' + sfx))
                log.debug(f"upPtsCmb count is: " + str(arcpy.GetCount_management(upPtsCmb).getOutput(0)))

                try:
                    arcpy.JoinField_management(upPtsCmb, gridfield2, upstreamCombine, 'value', [deepEnoughFr0Name, upPtsDemName, wsCostDistName])
                except:
                    copiedRows = arcpy.CopyRows_management(upstreamCombine, opj(inm, 'upstream_table'))
                    df.joinDict(upPtsCmb, gridfield2, copiedRows, 'value', [deepEnoughFr0Name, upPtsDemName, wsCostDistName])
                ## create consistent field names for all levels by altering name
                arcpy.AlterField_management(upPtsCmb, deepEnoughFr0Name, frFld)
                # if version.find('10.5') > -1 or version.find('10.6') > -1:
                #     arcpy.AlterField_management(upPtsCmb, df.getfields(upPtsCmb)[-1], 'WS_THKNS_CD')#frFld)
                # else:
                arcpy.AlterField_management(upPtsCmb, wsCostDistance.name, 'WS_THKNS_CD')#frFld)

                if arcpy.Exists(merged_medians):
                    mergedMdnsFc = arcpy.Clip_analysis(merged_medians, buffered_fc)

                    if df.testForZero(mergedMdnsFc):
        ####                if mergedMdns is not None:
                            # indicate fill regions that are in the median of a multilane road

                    ## Identify upstream points in median to enable differential processing (medians get half the extra distance)
                        idMedianPtsUp(upPtsCmb, medianFrFld, mergedMdnsFc)

                        ## Define if a transportation median should be used to increase search distance, and if so, which one
                        ## Calculate distance from upstream points that are near medians
                        ## use the minimum of this value to possibly increase search distance in these areas, if greater than minFrDistFld
        ####                df.conByList(
                        mergedMdnsCostMax = max([r[0] for r in arcpy.da.SearchCursor(mergedMdnsFc, ['bfr_dist'])])/2.0
                        upPtsMergedMdnsGnt = arcpy.GenerateNearTable_analysis(upPtsCmb, mergedMdnsFc, opj(sgdb, 'up_pts_mrg_mdn' + sfx), search_radius = mergedMdnsCostMax, closest = 'ALL')
                        arcpy.JoinField_management(upPtsMergedMdnsGnt, 'NEAR_FID', mergedMdnsFc, arcpy.Describe(mergedMdnsFc).OIDFieldName, 'bfr_dist')
                        arcpy.JoinField_management(upPtsMergedMdnsGnt, 'IN_FID', upPtsCmb, arcpy.Describe(upPtsCmb).OIDFieldName, frFld)
                        maxWsMedian = arcpy.Statistics_analysis(upPtsMergedMdnsGnt, opj(sgdb, 'ws_max_median_gnt' + sfx), [["bfr_dist", "MAX"], [frFld, "COUNT"]], frFld)

                    ## Calculate what fraction of the upstream fr points are in a median, if all, then fr is a 'median fr' and only gets half median distance                                        
                        upPtsCmbStats = arcpy.Statistics_analysis(upPtsCmb, opj(sgdb, 'up_pts_cmb_stats' + sfx), [[frFld, "COUNT"]], frFld)
                        df.tryAddField(maxWsMedian, mdnFracFld, 'DOUBLE')
                        arcpy.JoinField_management(maxWsMedian, frFld, upPtsCmbStats, frFld, 'COUNT_' + frFld)
                        arcpy.CalculateField_management(maxWsMedian, mdnFracFld, '!COUNT_' + frFld + '! / !COUNT_' + frFld + '_1!', 'PYTHON')

                        maxWsMdnFld1 = maxWsMedian.getInput(1).split(';')[0].split()[1] + '_' + maxWsMedian.getInput(1).split()[0]
                        log.debug('maxWsMdnFld from input query is ' + maxWsMdnFld1)
                    else:
                        ## set up dummy values
                        maxWsMedian = setupNoMedians(upPtsCmb, medianFrFld, dfs2cut4step, minElFld, frFld, maxWsMdnFld, mdnFracFld, sgdb)

                else:
                    ## set up dummy values
                    log.warning('No Medians to increase search distance for this watershed')
                    maxWsMedian = setupNoMedians(upPtsCmb, medianFrFld, dfs2cut4step, minElFld, frFld, maxWsMdnFld, mdnFracFld, sgdb)

                log.debug('done with maxWsMedian at ' + time.asctime())
                ## Set search distance for fill region,

                pntMdnMaxSearchDist = setAllSearchDistances(dfs2cut4step, wsSearchDistFld, minFrDistFld, frFld, minElFld, ofElFld, srchMult, crestMin, maxWsMedian, maxWsMdnFld, wsMdnSearchDistFld, pntMdnSearchDistFld, 500, mdnFracFld, frDepthFld)
                arcpy.JoinField_management(upPtsCmb, frFld, dfs2cut4step, frFld, pntMdnSearchDistFld)
                maxSearchDistList.append(pntMdnMaxSearchDist)

                df.copyfc(verbose, upPtsCmb, sgdb)

                arcpy.env.workspace = sgdb
                if df.testForZero(upPtsCmb):
                    log.debug('matching for ' + sfx)

                ## Create upstream point search radius
                    upPtsMbg = arcpy.MinimumBoundingGeometry_management(upPtsCmb, opj(inm, 'up_pts_mbg' + sfx), 'CONVEX_HULL', 'LIST', frFld)
                    arcpy.JoinField_management(upPtsMbg, frFld, dfs2cut4step, frFld, [wsSearchDistFld, pntMdnSearchDistFld, ofElFld, minElFld])#searchX2Fld)
                    df.copyfc(verbose, upPtsMbg, sgdb)

                    log.debug('matching Buffer for ' + sfx)
                ## Buffer by wsSearchDistFld (1/2 of point search distance) to find areas where points reach ws boundary
                    upMbgBfr = arcpy.Buffer_analysis(upPtsMbg, opj(inm, 'up_mbg_bfr' + sfx), wsSearchDistFld)#searchX2Fld)
                    df.copyfc(verbose, upMbgBfr, sgdb)

                    wsfcLine = arcpy.FeatureToLine_management(ws2cut4step, 'ws_2_line' + sfx)

                    intWsLines = arcpy.Intersect_analysis([wsfcLine, upMbgBfr], opj(inm, 'up_int' + sfx))#, output_type = 'LINE')
                    df.copyfc(verbose, intWsLines, sgdb)

                    selWsLines = arcpy.Select_analysis(intWsLines, opj(inm, 'sel_up' + sfx), frFld + ' = ' + frFld + '_1')
                    df.copyfc(verbose, selWsLines, sgdb)

                    log.debug('matching Dissolve for ' + sfx)
                    dslvWsLines = arcpy.Dissolve_management(selWsLines, opj(inm, 'dslv_ws_lines' + sfx), frFld)
                    arcpy.JoinField_management(dslvWsLines, frFld, dfs2cut4step, frFld, [wsSearchDistFld, wsMdnSearchDistFld])
                    df.copyfc(verbose, dslvWsLines, sgdb)

                    wsLineBfr = arcpy.Buffer_analysis(dslvWsLines, opj(inm, 'ws_line_bfr' + sfx), wsMdnSearchDistFld)#wsSearchDistFld)
                    df.copyfc(verbose, wsLineBfr, sgdb)

                    try:
                        wsLineBfrUnion = arcpy.Union_analysis([wsLineBfr, ws2cut4step], opj(inm, 'ws_line_bfr_unn' + sfx))
                        df.copyfc(verbose, wsLineBfrUnion, sgdb)
                    except:
                        log.warning("Invalid Topology, can't create wsLineBfrUnion, attempting repair")
                        arcpy.RepairGeometry_management(wsLineBfr)
                        arcpy.RepairGeometry_management(ws2cut4step)
                        wsLineBfrUnion = arcpy.Union_analysis([wsLineBfr, ws2cut4step], opj(inm, 'ws_line_bfr_unn' + sfx))
                        df.copyfc(verbose, wsLineBfrUnion, sgdb)
                        

                    if arcpy.Exists(merged_medians):
                        if df.testForZero(mergedMdnsFc):# is not None:#len(mergedMdnList) > 0:
                            wsLineBfrUnionErase = arcpy.Erase_analysis(wsLineBfrUnion, mergedMdnsFc, opj(inm, 'ws_line_bfr_union_no_mdns' + sfx))
                            if arcpy.Exists(wsLineBfrUnionErase):
                                df.copyfc(verbose, wsLineBfrUnionErase, sgdb)
                            else:
                                arcpy.RepairGeometry_management(wsLineBfrUnion)
                                wsLineBfrUnionErase = arcpy.Erase_analysis(wsLineBfrUnion, mergedMdnsFc, opj(inm, 'ws_line_bfr_union_no_mdns' + sfx))
                                df.copyfc(verbose, wsLineBfrUnionErase, sgdb)

                            area2search = arcpy.Select_analysis(wsLineBfrUnionErase, opj(inm, 'good_ds5' + sfx), frFld + ' <> ' + frFld + '_1 and FID_' + wsLineBfr[0].split('\\')[-1] + ' <> -1')
                            df.copyfc(verbose, area2search, sgdb)

                            area2searchWithMedians = arcpy.Select_analysis(wsLineBfrUnion, opj(inm, 'good_ds6' + sfx), frFld + ' <> ' + frFld + '_1 and FID_' + wsLineBfr[0].split('\\')[-1] + ' <> -1')
                            df.copyfc(verbose, area2searchWithMedians, sgdb)
                        else:
                            area2search = arcpy.Select_analysis(wsLineBfrUnion, opj(inm, 'good_ds5' + sfx), frFld + ' <> ' + frFld + '_1 and FID_' + wsLineBfr[0].split('\\')[-1] + ' <> -1')
                            df.copyfc(verbose, area2search, sgdb)

                            area2searchWithMedians = area2search
                    else:
                        area2search = arcpy.Select_analysis(wsLineBfrUnion, opj(inm, 'good_ds5' + sfx), frFld + ' <> ' + frFld + '_1 and FID_' + wsLineBfr[0].split('\\')[-1] + ' <> -1')
                        df.copyfc(verbose, area2search, sgdb)

                        area2searchWithMedians = area2search

            ## Convert the area2search to non-overlapping polygons that are 'rasterized'
                    area2searchCopy = arcpy.CopyFeatures_management(area2search, opj(inm, 'good_ds5_copy' + sfx))
                    df.copyfc(verbose, area2searchCopy, sgdb)
                    try:
                        area2searchSingles = arcpy.DeleteIdentical_management(area2searchCopy, 'SHAPE')
                        df.copyfc(verbose, area2searchSingles, sgdb)
                    except:
                        log.warning('Likely error in DeleteIdentical for sfx ' + sfx)
                        chckGeom = arcpy.CheckGeometry_management(area2searchCopy, opj(inm, 'chck_geom' + sfx))
    ####                                df.copyTbl(verbose, chckGeom, sgdb)
                        arcpy.JoinField_management(area2searchCopy, 'OBJECTID', chckGeom, 'FEATURE_ID', 'PROBLEM')

                        noProb = arcpy.Select_analysis(area2searchCopy, opj(inm, 'no_prob' + sfx), 'PROBLEM IS NULL')
                        df.copyfc(verbose, noProb, sgdb)
                        toFix = arcpy.Select_analysis(area2searchCopy, opj(inm, 'geom_prob' + sfx), 'PROBLEM IS NOT NULL')
                        df.copyfc(verbose, toFix, sgdb)
                        fixed = arcpy.RepairGeometry_management(toFix)
                        chckGeom2 = arcpy.CheckGeometry_management(fixed, opj(inm, 'chck_geom2' + sfx))
    ####                                df.copyTbl(verbose, chckGeom2, sgdb)
                        area2searchSingles = arcpy.Merge_management([noProb, fixed], opj(inm, 'good_ds5_fixed' + sfx))
                        df.copyfc(verbose, area2searchSingles, sgdb)


                    if df.testForZero(area2searchSingles):#int(arcpy.GetCount_management(area2searchSingles).getOutput(0)) > 0:
                ## Define linkages between overlapping and non-overlapping polygons (FID_area2searchSingles -> FID_area2search (this has multiple, overlapping polygons))
                        area2searchSinglesIntBack = arcpy.Intersect_analysis([area2search, area2searchSingles], opj(inm, 'a2s_singles_int_back' + sfx))
                        # fix 'Invalid Topology [Incomplete void poly.]' error
                        if not arcpy.Exists(area2searchSinglesIntBack):
                            log.warning('area2searchSinglesIntBack not completed successfully; attempting repair')
                            arcpy.RepairGeometry_management(area2search)
                            arcpy.RepairGeometry_management(area2searchSingles)
                            area2searchSinglesIntBack = arcpy.Intersect_analysis([area2search, area2searchSingles], opj(inm, 'a2s_singles_int_back' + sfx))
                        df.copyfc(verbose, area2searchSinglesIntBack, sgdb)

                ## Convert non-overlapping polygon to raster (FID_area2searchSingles -> VALUE (unique))
                        area2searchRasterResult = arcpy.PolygonToRaster_conversion(area2searchSingles, arcpy.Describe(area2searchSingles).OIDFieldName, opj(proc_dir, 'gd_ds5_rst'+ sfx), cellsize = ProcSize)
                ## Convert back to polygons (VALUE -> gridfield)
    ##                    zsBackToArea2Search = ZonalStatistics(wsLvl, 'value', fsPfdMdn, 'MINIMUM')
                        zstBackToArea2Search = ZonalStatisticsAsTable(area2searchRasterResult, 'value', fsPfdMdn, opj(inm, 'zst_el_srch_area_rast' + sfx))
                        arcpy.JoinField_management(area2searchSinglesIntBack, 'FID_' + area2searchSingles[0].split('\\')[-1], zstBackToArea2Search, 'value', 'MIN')
                        backToAreaSummaryStats = arcpy.Statistics_analysis(area2searchSinglesIntBack, opj(inm, 'bk2area_smry' + sfx), [['MIN', 'MIN']], frFld)
    ####                                df.copyTbl(verbose, backToAreaSummaryStats, sgdb)

                ## Define the minimum elevation in the search area
                        arcpy.JoinField_management(dfs2cut4step, frFld, backToAreaSummaryStats, frFld, minSummaryFld)
                    ## Only revise ones to cut
                        df.tryAddField(dfs2cut4step, revCutElFld, 'LONG')
                        with arcpy.da.UpdateCursor(dfs2cut4step, [cutElFld, revCutElFld, minSummaryFld]) as ucur:
                            for urow in ucur:
                                if urow[2] is not None:
                                    if urow[2] > urow[0]:
                                        urow[1] = urow[2]
                                    else:
                                        urow[1] = urow[0]
                                else:
                                    urow[1] = urow[0]
                                ucur.updateRow(urow)
                        df.copyfc(verbose, dfs2cut4step, sgdb)

                        area2searchDescribe = arcpy.da.Describe(area2searchWithMedians)
                        startDelIndex = df.getfields(area2searchWithMedians).index(frFld + '_1')
                        for fld in df.getfields(area2searchWithMedians)[startDelIndex:]:
                            if fld != area2searchDescribe['areaFieldName'] and fld != area2searchDescribe['lengthFieldName']:
                                arcpy.DeleteField_management(area2searchWithMedians, fld)

                        arcpy.JoinField_management(area2searchWithMedians, frFld, dfs2cut4step, frFld, [ofElFld, minElFld, revCutElFld, minSummaryFld])

                ## Select only parts of those FRs that are deeper than 50 cm or  revCutElFld < cutOneMinCrit * FR depth
                        area2search2 = arcpy.Select_analysis(area2searchWithMedians, opj(inm, 'good_ds8' + sfx), minSummaryFld + ' < (' + ofElFld + ' - 50) OR ' + minSummaryFld + ' < (' + minElFld + ' + ' + str(cutOneMinCrit) + '*(' + ofElFld + ' - ' + minElFld + '))')
                        if df.testForZero(area2search2):
                            df.copyfc(verbose, area2search2, sgdb)
                            maxEl4Dn = Raster(arcpy.PolygonToRaster_conversion(area2search2, revCutElFld, opj(proc_dir, "cut_el_dn" + sfx), "",revCutElFld, ProcSize))

                            areas2SearchEl = Raster(arcpy.PolygonToRaster_conversion(area2search2, minSummaryFld, opj(proc_dir, "srch_el_dn" + sfx), "",minSummaryFld, ProcSize))
                            log.debug("areas2SearchEl at " + time.asctime())


                        ## Thin potential downstream points by restricting to edges of large flats
                            allDnCells = Con(noInteriorFlatsDEM <= (maxEl4Dn), noInteriorFlatsDEM)
        ##                                allDnCells.save(opj(proc_dir, "all_dn"))# + sfx)
                            ## handle errors when the extent of allDnCells is less than blockstatistics height or width, results in 'Out of Memory' error?
                            if allDnCells.height < bsNbr.height or allDnCells.width < bsNbr.width:
                                arcpy.env.extent = noInteriorFlatsDEM
                                allDnCells = Con(noInteriorFlatsDEM <= (maxEl4Dn), noInteriorFlatsDEM)
                                arcpy.ClearEnvironment('extent')

                            ## Getting ERROR 000049: Failed to build attribute table on some areas
                            allDnCellsCopy = arcpy.CopyRaster_management(allDnCells, opj(sfldr, 'all_dn_cells' + sfx + '.tif'))
                            arcpy.BuildRasterAttributeTable_management(allDnCellsCopy)

                        ## Thin potential downstream points by calculating block statistics
                            bsMinDnEl = BlockStatistics(allDnCellsCopy, bsNbr, 'MINIMUM')
                            bsMinDnElName = str(bsMinDnEl).split('\\')[-1]
                            allDnCellsName =  str(allDnCells).split('\\')[-1]

                    ## Either need to rasterize a polygon dataset (to avoid losing small sliver polygons at intersections) accounting for overlap
                        ## or work with points and dissolve resulting buffers

                            if allDnCells.maximum != None:
                                arcpy.env.workspace = proc_dir
    # ------------------------------------------------------------------------------
    ##                            try:
    ##
    ##                                dnstreamCombine = Combine([allDnCells, bsMinDnEl, wsLvl])
    ##                                dnstreamCombine.save(opj(proc_dir, 'dn_cmb' + sfx))
    ##                                if df.testForZero(dnstreamCombine):
    ##                                    dnPtsCmb = arcpy.RasterToPoint_conversion(dnstreamCombine, opj(inm, "cb_dn" + sfx))
    ##
    ##                                if arcpy.Exists(merged_median):
    ##                                    if df.testForZero(mergedMdnsFc):# is not None:
    ##                                        idMedianPtsDn(dnPtsCmb, medianFrFld, mergedMdnsFc)
    ##                                    else:
    ##                                        arcpy.AddField_management(dnPtsCmb, medianFrFld, 'LONG')
    ##                                        arcpy.CalculateField_management(dnPtsCmb, medianFrFld, '0', 'PYTHON')
    ##                                else:
    ##                                    arcpy.AddField_management(dnPtsCmb, medianFrFld, 'LONG')
    ##                                    arcpy.CalculateField_management(dnPtsCmb, medianFrFld, '0', 'PYTHON')
    ##
    ##                                df.copyfc(verbose, dnPtsCmb, sgdb)
    ##                            ## Get the upstream region it might be a match for
    ##                                dnPts4step = arcpy.SpatialJoin_analysis(dnPtsCmb, area2search2, opj(inm, 'dn_pts_jn_bfr' + sfx), 'JOIN_ONE_TO_MANY', 'KEEP_ALL', gridfield2 + ' "' + gridfield2 + '" true true false 4 Long 0 0 , First, #, ' + dnPtsCmb[0] + ',' + gridfield2 + ',-1,-1; ' + frFld + ' "' + frFld + '" true true false 4 Long 0 0 , First, #, ' + area2search[0] + ',' + frFld + ',-1,-1')
    ##                                df.condDelete(verbose, dnPtsCmb)
    ##                                arcpy.JoinField_management(dnPts4step, gridfield2, dnstreamCombine, 'value', [allDnCellsName, bsMinDnElName, wsLvl.name])
    ##
    ##                            except arcpy.ExecuteError:
    ##                                # Since there are too many unique combinations we need to split the raster up and do more work
    ##                                # to handle all the different segments one at a time
    ##                                if u'Too many unique values' in arcpy.GetMessages():
    ##                                    log.warning('handling too many unique values error via SplitRaster')
    ##
    ##                                    splitname = 'ws_splt'
    ##                                    splitout = arcpy.SplitRaster_management(wsLvl, sfldr, splitname)
    ##                                    arcpy.env.workspace = sfldr
    ##                                    splits = arcpy.ListRasters(splitname + '*')
    ##                                    merge_tbls = []
    ##                                    log.info('done splitting raster into ' + str(splits))
    ##                                    for split in splits:
    ##                                        split_ext = os.path.splitext(split)[0]
    ##                                        add_on = split_ext.strip(splitname)
    ##                                        log.info('combining split ' + split + ' at ' + time.asctime())
    ####                                        splitDnRslt = arcpy.Clip_management(allDnCells, split.extent)
    ####                                        splitDn = ExtractByMask(allDnCells, split)
    ####                                        splitBs = ExtractByMask(bsMinDnEl, split)
    ####                                        dnstreamCombineSplit = Combine([splitDn, splitBs, split])
    ##                                        dnstreamCombineSplit = Combine([allDnCells, bsMinDnEl, split])
    ##                                        log.info('done combining split ' + split + ' at ' + time.asctime())
    ##                                        if df.testForZero(dnstreamCombineSplit):
    ##                                            dnPtsCmbSplit = arcpy.RasterToPoint_conversion(dnstreamCombineSplit, opj(inm, "cb_dn_split" + add_on + sfx))
    ##                                            log.info('done raster to point split ' + split + ' at ' + time.asctime())
    ####                                        split_tbl = arcpy.CopyRows_management(dnstreamCombineSplit, opj(sgdb, splitname + 'tbl' + add_on))
    ##                                            dnPts4stepSplit = arcpy.SpatialJoin_analysis(dnPtsCmbSplit, area2search2, opj(inm, 'dn_pts_jn_bfr' + add_on + sfx), 'JOIN_ONE_TO_MANY', 'KEEP_ALL')#,
    ##                                            log.info('done spatial join split ' + split + ' at ' + time.asctime())
    ####                                                gridfield2 + ' "' + gridfield2 + '" true true false 4 Long 0 0 , First, #, ' + dnPtsCmb[0] + ',' + gridfield2 + ',-1,-1; ' + frFld + ' "' + frFld + '" true true false 4 Long 0 0 , First, #, ' + area2search[0] + ',' + frFld + ',-1,-1')
    ##                                            arcpy.JoinField_management(dnPts4stepSplit, gridfield2, dnstreamCombineSplit, 'Value', [allDnCellsName, bsMinDnElName, split_ext])
    ##                                            arcpy.AlterField_management(dnPts4stepSplit, split_ext, wsLvl.name)
    ##                                            merge_tbls.append(dnPts4stepSplit)
    ##                                            df.condDelete(verbose, dnPtsCmbSplit)
    ##                                    dnPts4step = arcpy.Merge_management(merge_tbls, opj(inm, 'dn_pts_jn_bfr' + sfx))
    ##
    ######                                    dnPts4step = arcpy.SpatialJoin_analysis(dnPtsCmb, area2search2, opj(inm, 'dn_pts_jn_bfr' + sfx), 'JOIN_ONE_TO_MANY', 'KEEP_ALL',
    ######                                        gridfield2 + ' "' + gridfield2 + '" true true false 4 Long 0 0 , First, #, ' + dnPtsCmb[0] + ',' + gridfield2 + ',-1,-1; ' + frFld + ' "' + frFld + '" true true false 4 Long 0 0 , First, #, ' + area2search[0] + ',' + frFld + ',-1,-1')
    ##                                else:
    ##                                    log.warning('unhandled error')
    ##
    ##                            if df.testForZero(dnPts4step):#dnstreamCombine):
    ##                                arcpy.env.workspace = inm
    ##
    ######                                dnPtsCmb = arcpy.RasterToPoint_conversion(dnstreamCombine, "in_memory/cb_dn" + sfx)
    ##        ##                                    if len(mergedMdnList) > 0:
    ######                                if arcpy.Exists(merged_median):
    ######                                    if df.testForZero(mergedMdnsFc):# is not None:
    ######                                        idMedianPtsDn(dnPtsCmb, medianFrFld, mergedMdnsFc)
    ######                                    else:
    ######                                        arcpy.AddField_management(dnPtsCmb, medianFrFld, 'LONG')
    ######                                        arcpy.CalculateField_management(dnPtsCmb, medianFrFld, '0', 'PYTHON')
    ######                                else:
    ######                                    arcpy.AddField_management(dnPtsCmb, medianFrFld, 'LONG')
    ######                                    arcpy.CalculateField_management(dnPtsCmb, medianFrFld, '0', 'PYTHON')
    ######
    ######                                df.copyfc(verbose, dnPtsCmb, sgdb)
    ######                            ## Get the upstream region it might be a match for
    ######                                dnPts4step = arcpy.SpatialJoin_analysis(dnPtsCmb, area2search2, inm + 'dn_pts_jn_bfr' + sfx, 'JOIN_ONE_TO_MANY', 'KEEP_ALL', gridfield2 + ' "' + gridfield2 + '" true true false 4 Long 0 0 , First, #, ' + dnPtsCmb[0] + ',' + gridfield2 + ',-1,-1; ' + frFld + ' "' + frFld + '" true true false 4 Long 0 0 , First, #, ' + area2search[0] + ',' + frFld + ',-1,-1')
    ######                                df.condDelete(verbose, dnPtsCmb)
    ######                                arcpy.JoinField_management(dnPts4step, gridfield2, dnstreamCombine, 'value', [allDnCellsName, bsMinDnElName, wsLvl.name])
    ##
    ##                                ## generalize field names to they play nicely
    ##                                arcpy.AlterField_management(dnPts4step, allDnCellsName, allDnCellsFld)
    ##                                arcpy.AlterField_management(dnPts4step, wsLvl.name, wsLvlFld)
    ##                                arcpy.AlterField_management(dnPts4step, bsMinDnElName, bsMinDnElFld)
    ##                            ## Get the ds points that are deep enough for match
    ##                                arcpy.JoinField_management(dnPts4step, frFld, dfs2cut4step, frFld, [minElFld, ofElFld, cutElFld, revCutElFld, wsSearchDistFld, pntMdnSearchDistFld, minSummaryFld])
    ##                                dnPts4step2 = arcpy.Select_analysis(dnPts4step, inm + 'dn_pts_blk' + sfx, allDnCellsFld + ' = ' + bsMinDnElFld + ' OR ' + allDnCellsFld + ' = ' + minSummaryFld)
    ##                                df.condDelete(verbose, dnPts4step)
    ##                                dnPtsDeepEnough = arcpy.Select_analysis(dnPts4step2, inm + 'dn_dp_enuf' + sfx, allDnCellsFld + ' <= ' + revCutElFld)# ' + minElFld)
    ##                                df.copyfc(verbose, dnPts4step, sgdb)
    ##                                dnPtsDeepEnoughGdb = df.copyfc(verbose, dnPtsDeepEnough, sgdb)
    ##                                df.condDelete(verbose, dnPts4step2)
    ##
    ##                                dnPtsBst = dnPtsDeepEnough#Clip
    ##                                ptsMbg = arcpy.MinimumBoundingGeometry_management(dnPtsBst, inm + 'dp_pts_mbg' + sfx, 'CONVEX_HULL', 'LIST', frFld)
    ##                                arcpy.JoinField_management(ptsMbg, frFld, dfs2cut4step, frFld, pntMdnSearchDistFld)
    ##                                df.copyfc(verbose, ptsMbg, sgdb)
    ##
    ##                                dnMbgBfrOld = arcpy.Buffer_analysis(ptsMbg, inm + 'dn_mbg_bfr_old' + sfx, pntMdnSearchDistFld)
    ##                                df.copyfc(verbose, dnMbgBfrOld, sgdb)
    ##                                log.debug("dnMbgBfrOld at " + time.asctime())
    ##

    #-----------------------------------------------------------------------------------------------------
                                if True:
                                    # replace Combine with two Equal To and Boolean Or to make things work faster/more area in raster land
                                        # restrict downstream cells to those in 'watershed', not needed anymore?
                                        # this currently keeps it from finding flowing outwards downstream matches...
                                    allDnCellsWs = Con(wsLvl, allDnCells)
                                    bsMinDnElTest = allDnCellsWs == bsMinDnEl
                                    areas2SearchElTest = allDnCells == areas2SearchEl
                                    passedFirstElevationTest4Minimums = bsMinDnElTest | areas2SearchElTest

                                    passedFirstTestEl = Con(passedFirstElevationTest4Minimums, allDnCells)
                                    passedSecondElevationTest4Max = passedFirstTestEl <= maxEl4Dn
                                    passedSecondTestEl = Con(passedSecondElevationTest4Max, allDnCellsWs)
    ##                                passedSecondTestElName =  str(passedSecondTestEl).split('\\')[-1]
                                    # copy to a tif file to avoid the following error
    ##'''WARNING - ArcPy ERRORS:
    ##ERROR 999999: Error executing function.
    ##Class not registered

    ##ERROR 010423: ifthe_ras24.RASTER.1(Band_1) does not have valid statistics as required by the operation.
    ##ERROR 010069: Unable to open input raster(s).
    ##ERROR 010067: Error in executing grid expression.
    ##Failed to execute (Combine).'''
                                    secondDnCellsCopy = arcpy.CopyRaster_management(passedSecondTestEl, opj(sfldr, 'dn_cells2' + sfx + '.tif'))
    ##                                secondDnCellsCopy = arcpy.CopyRaster_management(allDnCells, opj(sfldr, 'all_dn_cells' + sfx + '.tif'))
                                    secondDnCellsCopyName =  os.path.splitext(os.path.basename(str(secondDnCellsCopy)))[0]
                                    arcpy.BuildRasterAttributeTable_management(secondDnCellsCopy)
                                    if df.testForZero(secondDnCellsCopy):

        ##                                dnPtsDeepEnough = arcpy.RasterToPoint_conversion(passedSecondTestEl, opj('in_memory', 'dn_pts_cnvrt' + sfx))

                                        dnstreamCombine = Combine([secondDnCellsCopy, bsMinDnEl, wsLvl])
                                        dnstreamCombine.save(opj(proc_dir, 'dn_cmb' + sfx))
                                        if df.testForZero(dnstreamCombine):
                                            dnPtsConversion = arcpy.RasterToPoint_conversion(dnstreamCombine, opj(inm, "cb_dn" + sfx))
                                            log.debug(f"dnPtsConversion count is: " + str(arcpy.GetCount_management(dnPtsConversion).getOutput(0)))

                                        if arcpy.Exists(merged_medians):
                                            if df.testForZero(mergedMdnsFc):# is not None:
                                                idMedianPtsDn(dnPtsConversion, medianFrFld, mergedMdnsFc)
                                            else:
                                                arcpy.AddField_management(dnPtsConversion, medianFrFld, 'LONG')
                                                arcpy.CalculateField_management(dnPtsConversion, medianFrFld, '0', 'PYTHON')
                                        else:
                                            arcpy.AddField_management(dnPtsConversion, medianFrFld, 'LONG')
                                            arcpy.CalculateField_management(dnPtsConversion, medianFrFld, '0', 'PYTHON')

                                    ## Get the upstream region it might be a match for
        ##arcpy.analysis.SpatialJoin("dn_pts_cnvrt_0", "good_ds8_0", r"O:\DEP\temp\Matcher_070801050901\Matcher_070801050901.gdb\dn_pts_cnvrt_0_SpatialJoin", "JOIN_ONE_TO_MANY", "KEEP_ALL", 'pointid "pointid" true true false 4 Long 0 0,First,#,dn_pts_cnvrt_0,pointid,-1,-1;grid_code "grid_code" true true false 4 Long 0 0,First,#,dn_pts_cnvrt_0,grid_code,-1,-1;FID_ws_line_bfr_0 "FID_ws_line_bfr_0" true true false 4 Long 0 0,First,#,good_ds8_0,FID_ws_line_bfr_0,-1,-1;FILL_RGN "FILL_RGN" true true false 4 Long 0 0,First,#,good_ds8_0,FILL_RGN,-1,-1;WS_SRCH_DST "WS_SRCH_DST" true true false 8 Double 0 0,First,#,good_ds8_0,WS_SRCH_DST,-1,-1;WS_MDN_DIST "WS_MDN_DIST" true true false 8 Double 0 0,First,#,good_ds8_0,WS_MDN_DIST,-1,-1;BUFF_DIST "BUFF_DIST" true true false 8 Double 0 0,First,#,good_ds8_0,BUFF_DIST,-1,-1;ORIG_FID "ORIG_FID" true true false 4 Long 0 0,First,#,good_ds8_0,ORIG_FID,-1,-1;FID_ws_2_cut_0 "FID_ws_2_cut_0" true true false 4 Long 0 0,First,#,good_ds8_0,FID_ws_2_cut_0,-1,-1;Id "Id" true true false 4 Long 0 0,First,#,good_ds8_0,Id,-1,-1;gridcode "gridcode" true true false 4 Long 0 0,First,#,good_ds8_0,gridcode,-1,-1;FR_MIN_EL "FR_MIN_EL" true true false 4 Long 0 0,First,#,good_ds8_0,FR_MIN_EL,-1,-1;FR_OF_EL "FR_OF_EL" true true false 4 Long 0 0,First,#,good_ds8_0,FR_OF_EL,-1,-1;MIN_MIN "MIN_MIN" true true false 8 Double 0 0,First,#,good_ds8_0,MIN_MIN,-1,-1;REV_CUT_EL "REV_CUT_EL" true true false 4 Long 0 0,First,#,good_ds8_0,REV_CUT_EL,-1,-1;Shape_Length "Shape_Length" false true true 8 Double 0 0,First,#,good_ds8_0,Shape_Length,-1,-1;Shape_Area "Shape_Area" false true true 8 Double 0 0,First,#,good_ds8_0,Shape_Area,-1,-1', "INTERSECT", None, '')
        ##arcpy.SpatialJoin_analysis(target_features="D:/DEP_Proc/DEMProc/Cut_dem2013_070801050901/scratch.gdb/dn_pts_cnvrt_0", join_features="D:/DEP_Proc/DEMProc/Cut_dem2013_070801050901/scratch.gdb/good_ds8_0", out_feature_class="D:/DEP_Proc/DEMProc/Cut_dem2013_070801050901/scratch.gdb/dn_pts_cnvrt_0_SpatialJoin2", join_operation="JOIN_ONE_TO_MANY", join_type="KEEP_ALL", field_mapping="""pointid "pointid" true true false 4 Long 0 0 ,First,#,D:\DEP_Proc\DEMProc\Cut_dem2013_070801050901\scratch.gdb\dn_pts_cnvrt_0,pointid,-1,-1;grid_code "grid_code" true true false 4 Long 0 0 ,First,#,D:\DEP_Proc\DEMProc\Cut_dem2013_070801050901\scratch.gdb\dn_pts_cnvrt_0,grid_code,-1,-1;FID_ws_line_bfr_0 "FID_ws_line_bfr_0" true true false 4 Long 0 0 ,First,#,D:\DEP_Proc\DEMProc\Cut_dem2013_070801050901\scratch.gdb\good_ds8_0,FID_ws_line_bfr_0,-1,-1;FILL_RGN "FILL_RGN" true true false 4 Long 0 0 ,First,#,D:\DEP_Proc\DEMProc\Cut_dem2013_070801050901\scratch.gdb\good_ds8_0,FILL_RGN,-1,-1;WS_SRCH_DST "WS_SRCH_DST" true true false 8 Double 0 0 ,First,#,D:\DEP_Proc\DEMProc\Cut_dem2013_070801050901\scratch.gdb\good_ds8_0,WS_SRCH_DST,-1,-1;WS_MDN_DIST "WS_MDN_DIST" true true false 8 Double 0 0 ,First,#,D:\DEP_Proc\DEMProc\Cut_dem2013_070801050901\scratch.gdb\good_ds8_0,WS_MDN_DIST,-1,-1;BUFF_DIST "BUFF_DIST" true true false 8 Double 0 0 ,First,#,D:\DEP_Proc\DEMProc\Cut_dem2013_070801050901\scratch.gdb\good_ds8_0,BUFF_DIST,-1,-1;ORIG_FID "ORIG_FID" true true false 4 Long 0 0 ,First,#,D:\DEP_Proc\DEMProc\Cut_dem2013_070801050901\scratch.gdb\good_ds8_0,ORIG_FID,-1,-1;FID_ws_2_cut_0 "FID_ws_2_cut_0" true true false 4 Long 0 0 ,First,#,D:\DEP_Proc\DEMProc\Cut_dem2013_070801050901\scratch.gdb\good_ds8_0,FID_ws_2_cut_0,-1,-1;Id "Id" true true false 4 Long 0 0 ,First,#,D:\DEP_Proc\DEMProc\Cut_dem2013_070801050901\scratch.gdb\good_ds8_0,Id,-1,-1;gridcode "gridcode" true true false 4 Long 0 0 ,First,#,D:\DEP_Proc\DEMProc\Cut_dem2013_070801050901\scratch.gdb\good_ds8_0,gridcode,-1,-1;FR_MIN_EL "FR_MIN_EL" true true false 4 Long 0 0 ,First,#,D:\DEP_Proc\DEMProc\Cut_dem2013_070801050901\scratch.gdb\good_ds8_0,FR_MIN_EL,-1,-1;FR_OF_EL "FR_OF_EL" true true false 4 Long 0 0 ,First,#,D:\DEP_Proc\DEMProc\Cut_dem2013_070801050901\scratch.gdb\good_ds8_0,FR_OF_EL,-1,-1;MIN_MIN "MIN_MIN" true true false 8 Double 0 0 ,First,#,D:\DEP_Proc\DEMProc\Cut_dem2013_070801050901\scratch.gdb\good_ds8_0,MIN_MIN,-1,-1;REV_CUT_EL "REV_CUT_EL" true true false 4 Long 0 0 ,First,#,D:\DEP_Proc\DEMProc\Cut_dem2013_070801050901\scratch.gdb\good_ds8_0,REV_CUT_EL,-1,-1;Shape_Length "Shape_Length" false true true 8 Double 0 0 ,First,#,D:\DEP_Proc\DEMProc\Cut_dem2013_070801050901\scratch.gdb\good_ds8_0,Shape_Length,-1,-1;Shape_Area "Shape_Area" false true true 8 Double 0 0 ,First,#,D:\DEP_Proc\DEMProc\Cut_dem2013_070801050901\scratch.gdb\good_ds8_0,Shape_Area,-1,-1""", match_option="INTERSECT", search_radius="", distance_field_name="")
                                        dnPtsDeepEnough = arcpy.SpatialJoin_analysis(dnPtsConversion, area2search2, opj(inm, 'dn_pts_sj_rstr' + sfx), 'JOIN_ONE_TO_MANY', 'KEEP_ALL', gridfield2 + ' "' + gridfield2 + '" true true false 4 Long 0 0 , First, #, ' + dnPtsConversion[0] + ',' + gridfield2 + ',-1,-1; ' + frFld + ' "' + frFld + '" true true false 4 Long 0 0 , First, #, ' + area2search[0] + ',' + frFld + ',-1,-1')
                                        df.copyfc(verbose, dnPtsDeepEnough, sgdb)

                                        try:
                                            arcpy.JoinField_management(dnPtsDeepEnough, gridfield2, dnstreamCombine, 'value', [secondDnCellsCopyName, bsMinDnElName, wsLvl.name])
                                        except:
                                            dnStreamRows = arcpy.CopyRows_management(dnstreamCombine, opj(inm, 'dnstreamCombine'))
                                            df.joinDict(dnPtsDeepEnough, gridfield2, dnStreamRows, 'value', [secondDnCellsCopyName, bsMinDnElName, wsLvl.name])


                                        ## generalize field names to they play nicely
                                        arcpy.AlterField_management(dnPtsDeepEnough, secondDnCellsCopyName, allDnCellsFld)
                                        arcpy.AlterField_management(dnPtsDeepEnough, wsLvl.name, wsLvlFld)
                                        arcpy.AlterField_management(dnPtsDeepEnough, bsMinDnElName, bsMinDnElFld)

                                        ptsMbgRstr = arcpy.MinimumBoundingGeometry_management(dnPtsDeepEnough, opj(inm, 'dp_pts_mbg_rstr' + sfx), 'CONVEX_HULL', 'LIST', frFld)

                                        arcpy.JoinField_management(ptsMbgRstr, frFld, dfs2cut4step, frFld, pntMdnSearchDistFld)
                                        df.copyfc(verbose, ptsMbgRstr, sgdb)

                                        dnMbgBfr = arcpy.Buffer_analysis(ptsMbgRstr, opj(inm, 'dn_mbg_bfr' + sfx), pntMdnSearchDistFld)
                                        df.copyfc(verbose, dnMbgBfr, sgdb)
                                        log.debug("dnMbgBfr at " + time.asctime())

        #-----------------------------------------------------------------------------------------------------


                                    ## Reduce the upstream points to just those in search area
                                        upPtsClip2 = arcpy.Clip_analysis(upPtsCmb.getOutput(0)[:-2] + sfx, dnMbgBfr, opj(inm, 'up_pts_clip2' + sfx))
                                        df.copyfc(verbose, upPtsClip2, sgdb)
                                        df.condDelete(verbose, upPtsCmb)

                                    ## Must use the DEM we built the up points with
                                        deepFrName = frFld
                                        distField = pntMdnSearchDistFld
                                        dn_demName = allDnCellsFld#str(allDnCells).split('\\')[-1]

                                    ## to track uniquely the 
                ##                                    gnt = arcpy.GenerateNearTable_analysis(upPts, dnPts, 'in_memory\\' + gnt1Name, maxSearchDist, "LOCATION", "NO_ANGLE", "ALL")
                ##                                    log.debug("started gnt at " + time.asctime())
                                        gnt = "gnt_1" + sfx
                                        inmgnt = arcpy.GenerateNearTable_analysis(upPtsClip2, dnPtsDeepEnough, opj(inm, gnt), maxSearchDistList[i], "LOCATION", "NO_ANGLE", "ALL")#, method = 'GEODESIC')
            ##                            SSgnt = arcpy.GenerateNearTable_analysis(upPtsClip2, dnPtsDeepEnough, SS + gnt, maxSearchDistList[i], "LOCATION", "NO_ANGLE", "ALL")#, method = 'GEODESIC')
                ##                                    log.debug("finished gnt at " + time.asctime())
                ####                                    df.copyTbl(verbose, SSgnt, sgdb)
                                        log.debug(f"inmgnt count is: " + str(arcpy.GetCount_management(inmgnt).getOutput(0)))

            ##                            upPts = copyfc2(upPtsClip2, SS)#Bst, SS)
            ##                            upPtsCount = float(arcpy.GetCount_management(upPts).getOutput(0))
                                        upPtsCount = float(arcpy.GetCount_management(upPtsClip2).getOutput(0))
                                        upPtsRatio = upPtsCount/float(arcpy.GetCount_management(dfs2cut4step).getOutput(0))
            ##                            dnPts = copyfc2(dnPtsDeepEnough, SS)#Bst, SS)
            ##                            dnPtsCount = float(arcpy.GetCount_management(dnPts).getOutput(0))
                                        dnPtsCount = float(arcpy.GetCount_management(dnPtsDeepEnough).getOutput(0))
                                        dnPtsRatio = dnPtsCount/float(arcpy.GetCount_management(dfs2cut4step).getOutput(0))
                ####                                    upDnProduct = upPtsCount * dnPtsCount
                                        log.debug('upPtsRatio is ' + str(upPtsRatio) + ' and dnPtsRatio is ' + str(dnPtsRatio))

                                        iCurPts = arcpy.da.InsertCursor(match_stats, ['HUC12', 'LEVEL', 'UP_PTS_COUNT', 'DN_PTS_COUNT', 'UP_PTS_RATIO', 'DN_PTS_RATIO'])#, 'UP_DN_PRODUCT'])
                                        iRowPts = iCurPts.insertRow([huc12, i, upPtsCount, dnPtsCount, upPtsRatio, dnPtsRatio])#, upDnProduct])
                                        del iRowPts, iCurPts
                ##                                    assert 
            ##                            SSdfs2cut4step = copyfc2(dfs2cut4step, SS)
                                        dfs2CutList.append(dfs2cut4step)
                                        
            ##                            upName = tno(upPts)
            ##                            df.condDelete(verbose, upPts)
            ##                            dnName = tno(dnPts)
            ##                            dfs = tno(SSdfs2cut4step)

            ##                            ## Select initial matches and put associated data into new table
            ####                            DB + upName + "." + str(upPtsDEM.name) + ", " + DB + upName + "." + frFld + ", " + DB + upName + "." + medianFrFld
            ####                            DB + dfs + "." + minElFld + ", " + DB + dfs + "." + ofElFld + ", " + DB + dfs + "." + revCutElFld + ", " + DB + dfs + "." + pntMdnSearchDistFld + ", " + DB + dfs + "." + wsSearchDistFld + ", " + DB + dfs + "." + mdnFracFld + ", " + DB + dfs + "." + maxWsMdnFld
            ####                            DB + dnName + "." + wsLvlFld + ", " + DB + dnName + "." + str(dn_demName)
            ##                                                                                          
            ##                            sql2pt1 = "SELECT " + DB + gnt + ".OBJECTID, " + DB + gnt + ".IN_FID, " + DB + gnt + ".NEAR_FID, " + DB + gnt + ".NEAR_DIST, NEAR_X, NEAR_Y, FROM_X, FROM_Y, " + DB + upName + "." + str(upPtsDEM.name) + ", " + DB + upName + "." + frFld + ", " + DB + upName + "." + medianFrFld + ", " + DB + dfs + "." + minElFld + ", " + DB + dfs + "." + ofElFld + ", " + DB + dfs + "." + revCutElFld + ", " + DB + dfs + "." + pntMdnSearchDistFld + ", " + DB + dfs + "." + wsSearchDistFld + ", " + DB + dfs + "." + mdnFracFld + ", " + DB + dfs + "." + maxWsMdnFld + ", " + DB + dnName + "." + wsLvlFld + ", " + DB + dnName + "." + str(dn_demName) + " INTO " + DB + gnt2 + " FROM " + DB + gnt
            ##    ####                                    sql2pt1 = "SELECT " + DB + gnt + ".OBJECTID, " + DB + gnt + ".IN_FID, " + DB + gnt + ".NEAR_FID, " + DB + gnt + ".NEAR_DIST, NEAR_X, NEAR_Y, FROM_X, FROM_Y, " + DB + upName + "." + str(upPtsDEM.name) + ", " + DB + upName + "." + frFld + ", " + DB + upName + "." + medianFrFld + ", " + DB + dfs + "." + minElFld + ", " + DB + dfs + "." + ofElFld + ", " + DB + dfs + "." + revCutElFld + ", " + DB + dfs + "." + pntMdnSearchDistFld + ", " + DB + dfs + "." + wsSearchDistFld + ", " + DB + dfs + "." + mdnFracFld + ", " + DB + dfs + "." + maxWsMdnFld + ", " + DB + dfs + "." + 'DS_BEARING' + ", " + DB + dnName + "." + wsLvlFld + ", " + DB + dnName + "." + str(dn_demName) + " INTO " + DB + gnt2 + " FROM " + DB + gnt
            ##                            sql2Join = " INNER JOIN " + DB + upName + " ON IN_FID = " + upName + ".OBJECTID INNER JOIN " + DB + dfs + " ON " + DB + upName + "." + frFld + " = " + DB + dfs + "." + frFld + " INNER JOIN " + DB + dnName + " ON " + DB + gnt + ".NEAR_FID = " + DB + dnName + ".OBJECTID"
            ##                            ## Build statements so matches must be nearer that total search distance, appropriate elevation, and not in same watershed as fill region
            ##                            where3 = " AND " + DB + gnt + ".NEAR_DIST <= " + DB + dfs + "." + pntMdnSearchDistFld + " AND " + DB + gnt + ".NEAR_DIST > " + str(ProcSize + ProcSize/2.0) + " AND " + DB + dnName + "." + str(dn_demName) + " <= " + DB + dfs + "." + revCutElFld + " AND " + DB + dfs + "." + frFld + " <> " + DB + dnName + "." + wsLvlFld
            ##
            ##                            sql2 = sql2pt1 + sql2Join + where3
            ##                            print(time.clock());sdeConn.execute(sql2);print(time.clock())
            ##                            sdeConn.commitTransaction()

                                        # inmgnt = arcpy.GenerateNearTable_analysis(upPtsClip2, dnPtsDeepEnough, inm + gnt, maxSearchDistList[i], "LOCATION", "NO_ANGLE", "ALL")

    ##                                    print(time.clock())
                                        log.debug('starting upPtsFlds at ' + time.asctime())
                                        # upPtsFlds = [deepBsMinUpElName, frFld, medianFrFld]
                                        # upPtsFlds = [df.getfields(upPtsCmb)[5], frFld, medianFrFld]#upPtsDEM.name
                                        upPtsFlds = [upPtsDemName, frFld, medianFrFld]
                                        dfsFlds = [minElFld, ofElFld, revCutElFld, pntMdnSearchDistFld, wsSearchDistFld, mdnFracFld, maxWsMdnFld]
                                        dnPtsFlds = [wsLvlFld, dn_demName]
                                        allFlds = upPtsFlds + dnPtsFlds + dfsFlds
                                        for field in allFlds:
                                            if field in allFlds[:-4]:
                                                arcpy.AddField_management(inmgnt, field, 'DOUBLE')
                                            else:
                                                arcpy.AddField_management(inmgnt, field, 'LONG')
    ##                                    print(time.clock())
                                        log.debug('done with upPtsFlds at ' + time.asctime())
                                        dictValue1 = {r[0]:(r[1:]) for r in arcpy.da.SearchCursor(upPtsClip2, ['OBJECTID'] + upPtsFlds)}
                                        #AG Pro 3.2
                                        with arcpy.da.UpdateCursor(dfs2cut4step, [mdnFracFld]) as ucur_dfs:
                                            for urow_dfs in ucur_dfs:
                                                # print(urow_dfs)
                                                if urow_dfs[0] is None:
                                                    urow_dfs [0] = 0
                                                    ucur_dfs.updateRow(urow_dfs)

                                        dictValue2 = {r[0]:(r[1:]) for r in arcpy.da.SearchCursor(dfs2cut4step, [frFld] + dfsFlds)}
        ##                                dnPtsDesc = arcpy.Describe(dnPts4step2)
                                        dnPtsDesc = arcpy.Describe(dnPtsDeepEnough)
                                        dictValue3 = {r[0]:(r[1:]) for r in arcpy.da.SearchCursor(dnPtsDeepEnough, [dnPtsDesc.OIDFieldName] + dnPtsFlds)}
                                        with arcpy.da.UpdateCursor(inmgnt, ['IN_FID'] + ['NEAR_FID'] + upPtsFlds + dfsFlds + dnPtsFlds) as ucur:
                                            for urow in ucur:
                                                ## Time comparison on 070801050901 - all using .index() 2.8609s, removing frFld .index() 2.2345s, no .index() 0.6103s
            ##                                    for i, field1 in enumerate(upPtsFlds):
                                                    urow[2] = dictValue1[urow[0]][0]
                                                    urow[3] = dictValue1[urow[0]][1]
                                                    urow[4] = dictValue1[urow[0]][2]
            ##                                        urow[ucur.fields.index(field1)] = dictValue1[urow[0]][i]
            ##                                    for i, field2 in enumerate(dfsFlds):
                                                    urow[5] = dictValue2[urow[3]][0]
                                                    urow[6] = dictValue2[urow[3]][1]
                                                    urow[7] = dictValue2[urow[3]][2]
                                                    urow[8] = dictValue2[urow[3]][3]
                                                    urow[9] = dictValue2[urow[3]][4]
                                                    urow[10] = dictValue2[urow[3]][5]
                                                    urow[11] = dictValue2[urow[3]][6]
            ##                                        urow[ucur.fields.index(field2)] = dictValue2[urow[ucur.fields.index(frFld)]][i]
            ##                                    for i, field3 in enumerate(dnPtsFlds):
                                                    urow[12] = dictValue3[urow[1]][0]
                                                    urow[13] = dictValue3[urow[1]][1]
            ##                                        urow[ucur.fields.index(field3)] = dictValue3[urow[1]][i]
                                                    ucur.updateRow(urow)
                                        log.debug('done with updating inmgnt at ' + time.asctime())

                                        where4 = "NEAR_DIST <= " + pntMdnSearchDistFld + " AND NEAR_DIST > " + str(ProcSize + ProcSize/2.0) + " AND " + str(dn_demName) + " <= " + revCutElFld + " AND " + frFld + " <> " + wsLvlFld

                                        gnt2 = "gnt_2" + sfx
                                        inmgnt2 = arcpy.TableSelect_analysis(inmgnt, opj(inm, gnt2), where4)
                                                    
                                        if arcpy.Exists(merged_medians):
                                            if df.testForZero(mergedMdnsFc):# is not None:#len(mergedMdnList) > 0:

                ##                                gntFields = df.getfields(SS + gnt2)[1:]
                                                gntFields = df.getfields(inmgnt2)[1:]
                                                log.debug(f"inmgnt2[1:] fields are: {gntFields}")
                                                # inFidIndx = gntFields.index('IN_FID')
                                                log.debug('processing in merged medians at ' + time.asctime())

                                            ## Figure out which point matches need to cross an extra buffer feature and select those that do and beyond base search distance
                                            ## Then create a polyline fc of unique crossings
                ##                                gntFrs2Test4X = arcpy.TableSelect_analysis(SS + gnt2, inm + 'gnt_3' + sfx, 'NEAR_DIST > (' + pntMdnSearchDistFld + ' - ' + maxWsMdnFld + ')')
                                                gntFrs2Test4X = arcpy.TableSelect_analysis(inmgnt2, opj(inm, 'gnt_3' + sfx), 'NEAR_DIST > (' + pntMdnSearchDistFld + ' - ' + maxWsMdnFld + ')')
                                                log.debug(f"gntFrs2Test4X fields are: " + str(df.getfields(gntFrs2Test4X)))

                                                # need data in a File GDB to use sql_clause...
                                                gdbGntFrs2Test4X = arcpy.CopyRows_management(gntFrs2Test4X, opj(sgdb, gntFrs2Test4X[0].split('\\')[1]))
                                                frWsList = []
                                                with arcpy.da.SearchCursor(gdbGntFrs2Test4X, ['FILL_RGN', 'ws_lvl'], sql_clause = (None, 'GROUP BY ' + frFld + ', ' + wsLvlFld + ' ORDER BY ' + frFld + ', ' + wsLvlFld)) as scur:
                                                    for srow in scur:
                                                        frWsList.append(srow)

                                            ## Just get the first, nearest, match for each upstream point to check for crossing extra buffer feature
                                                gntXLines = arcpy.CreateTable_management(inm, 'gnt_tr_lines_tbl' + sfx, inmgnt2)
                ##                                gntXLines = arcpy.CreateTable_management(inm, 'gnt_tr_lines_tbl' + sfx, SS + gnt2)
                                                frWsId = 'FR_WS_ID'
                                                df.tryAddField(gntXLines, frWsId, 'LONG')
                                                icur = arcpy.da.InsertCursor(gntXLines, gntFields + [frWsId])

                                                for frWsIndex, item in enumerate(frWsList):
                                                    scur = arcpy.da.SearchCursor(gntFrs2Test4X, gntFields, frFld + ' = ' + str(item[0]) + ' AND ws_lvl = ' + str(item[1]))
                                                    srow = scur.next()
                                                    irow = list(srow)
                                                    irow.append(frWsIndex)
                                                    icur.insertRow(irow)
                                                    del srow, scur
                                                del icur

                                                sr = arcpy.Describe(fill_or_void_tif).spatialReference
                                                trXLines = arcpy.XYToLine_management(gntXLines, opj(inm, 'tr_x_lines' + sfx), 'FROM_X', 'FROM_Y', 'NEAR_X', "NEAR_Y", id_field = frWsId, spatial_reference = sr)
                                                df.copyfc(verbose, trXLines, sgdb)

                                                trXInt = arcpy.Intersect_analysis([trXLines, mergedMdnsFc], opj(inm, 'tr_x_lines_int' + sfx), output_type = 'POINT')
                                                df.copyfc(verbose, trXInt, sgdb)

                                                trXIntSummary = arcpy.Statistics_analysis(trXInt, opj(inm, 'tr_x_int_smry' + sfx), [[frWsId, 'COUNT']], frWsId)
                                                df.copytbl(verbose, trXIntSummary, sgdb)

                                                arcpy.JoinField_management(trXIntSummary, frWsId, gntXLines, frWsId, [frFld, 'ws_lvl'])

                                                selStr = ''
                                                with arcpy.da.SearchCursor(trXIntSummary, [frFld, 'ws_lvl']) as scur:
                                                    for srow in scur:
                                                        selStr += frFld + ' = ' +str(srow[0]) + ' AND ws_lvl = ' + str(srow[1]) + ' OR '

                                                gntFrsX = arcpy.TableSelect_analysis(gntFrs2Test4X, opj(inm, 'gnt_4a_fr_x' + sfx), selStr[:-4])#'COUNT_' + frWsId + ' >= 1')
                                                df.copytbl(verbose, gntFrsX, sgdb)

                ##                                gntFrsNoX = arcpy.TableSelect_analysis(SS + gnt2, inm + 'gnt_4b_fr_no_x' + sfx, 'NEAR_DIST <= (' + pntMdnSearchDistFld + ' - ' + maxWsMdnFld + ') OR ' + maxWsMdnFld + ' IS NULL')
                                                gntFrsNoX = arcpy.TableSelect_analysis(inmgnt2, opj(inm, 'gnt_4b_fr_no_x' + sfx), 'NEAR_DIST <= (' + pntMdnSearchDistFld + ' - ' + maxWsMdnFld + ') OR ' + maxWsMdnFld + ' IS NULL')
                                                df.copytbl(verbose, gntFrsNoX, sgdb)

                                                gntXNoX = arcpy.Merge_management([gntFrsNoX, gntFrsX], opj(inm, 'gnt_5_x_no_x' + sfx))
                    ####                                        gdbGntXNoX = df.copyTbl(verbose, gntXNoX, sgdb)

                                            ## prefer those points not in median for downstream matches, if there are any
                                                frXSummary = arcpy.Statistics_analysis(gntXNoX, opj(inm, 'mdn_x_smry' + sfx), [[medianFrFld, 'MEAN']], frFld)
                                                df.copytbl(verbose, frXSummary, sgdb)
                    ##                                                    arcpy.JoinField_management(gntXNoX, frFld, frXSummary, frFld, 'MEAN_' + medianFrFld)
                    ##                                                    gntRemovedMedianMatches = arcpy.TableSelect_analysis(gntXNoX, inm + 'gnt_no_mdn_match' + sfx, 'MEAN_' + medianFrFld + ' = 1.0 AND ' + medianFrFld + ' = 1 OR MEAN_' + medianFrFld

                                            else:
                                            ## copy stuff over from pre-median testing
                ##                                gntXNoX = arcpy.CopyRows_management(SS + gnt2, inm + 'gnt_5_x_no_x' + sfx)
                                                gntXNoX = arcpy.CopyRows_management(inmgnt2, opj(inm, 'gnt_5_x_no_x' + sfx))
                    ####                                        gdbGntXNoX = df.copyTbl(verbose, gntXNoX, sgdb)
                                        else:
                                        ## copy stuff over from pre-median testing
            ##                                gntXNoX = arcpy.CopyRows_management(SS + gnt2, inm + 'gnt_5_x_no_x' + sfx)
                                            gntXNoX = arcpy.CopyRows_management(inmgnt2, opj(inm, 'gnt_5_x_no_x' + sfx))
                ####                                        gdbGntXNoX = df.copyTbl(verbose, gntXNoX, sgdb)

                                        log.debug('done with medians at ' + time.asctime())

                                    ## further winnow upstream points for lowest elevation
                                        statsUpCells = arcpy.Statistics_analysis(gntXNoX, opj(inm, 'stat_up_pts_dem' + sfx), [[upPtsDemName, 'MIN'], [upPtsDemName, 'MEAN']], frFld)#, [bestestDEM.name, 'STD']#str(upPtsDEM.name)
                                        minShort = str('MIN_' + upPtsDemName)
                                        meanShort = str('MEAN_' + upPtsDemName)
                                        df.joinDict(gntXNoX, frFld, statsUpCells, frFld, [meanShort, minShort])
                                        # arcpy.JoinField_management(gntXNoX, frFld, statsUpCells, frFld, ['MEAN_' + upPtsDemName, 'MIN_' + upPtsDemName])#str(upPtsDEM.name)
                                        gntXNoXUpBstDeep = arcpy.TableSelect_analysis(gntXNoX, opj(inm, 'gnt_7_up_bst_dp' + sfx), upPtsDemName + ' <= ' + revCutElFld + ' OR (' + upPtsDemName + ' <= ' + meanShort + ' AND ' + minShort + ' > ' + revCutElFld + ')')

                                        log.debug('done with gntXNoXUpBstDeep at ' + time.asctime())

                                        ## removed gntAfterDsBearing due to slow execution of calcDsFd function in production environment
                ####                                ## winnow downstream points for best downstream bearing
                ####                                    df.tryAddField(gntXNoXUpBstDeep, 'NA_BEARING', 'DOUBLE')
                ####                                    df.tryAddField(gntXNoXUpBstDeep, 'DS_BEAR_DIF', 'DOUBLE')
                ####                                    with arcpy.da.UpdateCursor(gntXNoXUpBstDeep, ['DS_BEAR_DIF', 'NA_BEARING', 'DS_BEARING', 'FROM_X', 'FROM_Y', 'NEAR_X', "NEAR_Y"]) as ucur:
                ####                                        for urow in ucur:
                ####                                            if urow[2] != None:
                ####
                ####                                                pnt_dx = urow[ucur.fields.index('NEAR_X')] - urow[ucur.fields.index('FROM_X')]
                ####                                                pnt_dy = urow[ucur.fields.index('NEAR_Y')] - urow[ucur.fields.index('FROM_Y')]
                ####                                                bearingAngle = math.atan2(pnt_dx, pnt_dy)*(360.0/(2*math.pi))
                ####                                                urow[1] = bearingAngle
                ####
                ####                                                urow[0] = angleDif(urow[1], urow[2])
                ####                                                ucur.updateRow(urow)
                ####                                    statsDsBearing = arcpy.Statistics_analysis(gntXNoXUpBstDeep, inm + 'stat_ds_bear_dif' + sfx, [['DS_BEAR_DIF', 'MIN'], [wsLvlFld, 'RANGE']], frFld)
                ####                                    arcpy.JoinField_management(gntXNoXUpBstDeep, frFld, statsDsBearing, frFld, ['MIN_DS_BEAR_DIF', 'RANGE_' + wsLvlFld])
                ####                                    gntAfterDsBearing = arcpy.TableSelect_analysis(gntXNoXUpBstDeep, inm + 'gnt_7b_ds_bearing' + sfx, 'DS_BEARING IS NULL OR RANGE_' + wsLvlFld + ' = 0 OR ((MIN_DS_BEAR_DIF < 180.0 AND DS_BEAR_DIF < 180.0) OR MIN_DS_BEAR_DIF >= 180)')


                                    ## further winnow downstream points by equal or less than mean slope
                                        udSlpFld = 'up_dn_slp_pct'
                                        df.tryAddField(gntXNoXUpBstDeep, udSlpFld, 'double')

                                    ## Calculate slope between points
                                        with arcpy.da.UpdateCursor(gntXNoXUpBstDeep, ['up_dn_slp_pct', 'NEAR_DIST', revCutElFld, dn_demName]) as ucur:
                                            for urow in ucur:
                                                if urow[1] < ProcSize/2:
                                                    urow[0] = -999
                                                else:
                                                    urow[0] = (urow[2] - urow[3]) / urow[1]
                                                ucur.updateRow(urow)

                                        ptsSlpSummary = arcpy.Statistics_analysis(gntXNoXUpBstDeep, opj(inm, 'pts_slp_sum_' + sfx), [['up_dn_slp_pct', 'MAX']], deepFrName)#'NEAR_FID')

                                        arcpy.JoinField_management(gntXNoXUpBstDeep, deepFrName, ptsSlpSummary, deepFrName, ['MAX_' + udSlpFld])
                                        df.copytbl(verbose, gntXNoXUpBstDeep, sgdb)

                                        dpSelection = '(' + udSlpFld + ' >= MAX_' + udSlpFld + '/2.0 AND MAX_' + udSlpFld + ' > 2.0) OR ((' + udSlpFld + ' + 1.0) >= MAX_' + udSlpFld + ' AND MAX_' + udSlpFld + ' <= 2.0)'
                                        gntAfterDsBearing2 = arcpy.TableSelect_analysis(gntXNoXUpBstDeep, opj(inm, 'gnt_8_up_bst_dp2' + sfx), dpSelection)
                                        df.copytbl(verbose, gntAfterDsBearing2, sgdb)

                                        gdbgntAfterDsBearing2 = arcpy.CopyRows_management(gntAfterDsBearing2, opj(sgdb, gntAfterDsBearing2[0].split('\\')[1]))
                                        log.debug('done with gdbgntAfterDsBearing2 at ' + time.asctime())

                                        if df.testForZero(gdbgntAfterDsBearing2):#int(arcpy.GetCount_management(gdbgntAfterDsBearing2).getOutput(0)) > 0:

                                        ## Get just those downstream point matches that remain (still not suitable for matching with all fill regions as a ds match may work for some FRs and not others)
                                            dnPtsBstDeepSummary = arcpy.Statistics_analysis(gntAfterDsBearing2, opj(inm, 'dn_pts_cnt_smry' + sfx), [['NEAR_FID', 'COUNT']], 'NEAR_FID')
                                            df.copytbl(verbose, dnPtsBstDeepSummary, sgdb)

            ##                                df.addCalcJoin(dnPts, 'OBJECTID', dnPtsBstDeepSummary, 'NEAR_FID', ['CNT_NR_FID_1', 'LONG'], '!COUNT_NEAR_FID!')
            ##                                dnPtsBst1pt5 = arcpy.Select_analysis(dnPts, inm + 'dn_pts_bst1pt5' + sfx, 'CNT_NR_FID_1 >= 1')
                                            dnPtsDeepEnoughDesc = arcpy.Describe(dnPtsDeepEnough)
                                            df.addCalcJoin(dnPtsDeepEnough, dnPtsDeepEnoughDesc.OIDFieldName, dnPtsBstDeepSummary, 'NEAR_FID', ['CNT_NR_FID_1', 'LONG'], '!COUNT_NEAR_FID!')
                                            dnPtsBst1pt5 = arcpy.Select_analysis(dnPtsDeepEnough, opj(inm, 'dn_pts_bst1pt5' + sfx), 'CNT_NR_FID_1 >= 1')
                                            df.copyfc(verbose, dnPtsBst1pt5, sgdb)

                                        ## Select unique FR and WS pairs for downstream points 
                                            gdbgntAfterDsBearing2 = arcpy.CopyRows_management(gntAfterDsBearing2, opj(sgdb, gntAfterDsBearing2[0].split('\\')[1]))
                                            dnPtsBst2 = selectFrAndWsDpPts(gdbgntAfterDsBearing2, dnPtsBst1pt5, frFld, sfx, log)
                                            dnPtsBst2Gdb = df.copyfc(verbose, dnPtsBst2, sgdb)

                                        ## Get just those upstream point matches that remain
                                            upPtsBstDeepSummary = arcpy.Statistics_analysis(gntAfterDsBearing2, opj(inm, 'up_pts_cnt_smry' + sfx), [['IN_FID', 'COUNT']], 'IN_FID')
                                            df.copytbl(verbose, upPtsBstDeepSummary, sgdb)
            ##                                arcpy.JoinField_management(upPts, 'OBJECTID', upPtsBstDeepSummary, 'IN_FID', 'COUNT_IN_FID')
            ##                                upPtsBst = arcpy.Select_analysis(upPts, inm + 'up_pts_bst' + sfx, 'COUNT_IN_FID >= 1')
                                            arcpy.JoinField_management(upPtsClip2, 'OBJECTID', upPtsBstDeepSummary, 'IN_FID', 'COUNT_IN_FID')
                                            upPtsBst = arcpy.Select_analysis(upPtsClip2, opj(inm, 'up_pts_bst' + sfx), 'COUNT_IN_FID >= 1')
                                            upPtsBstGdb = df.copyfc(verbose, upPtsBst, sgdb)

                                            # upPtsBstDistSummary = arcpy.Statistics_analysis(gntAfterDsBearing2, inm + 'pts_dist_smry' + sfx, [['NEAR_DIST', 'MIN'],['NEAR_DIST', 'MEAN']], frFld)
                                            # log.debug('done with upPtsBstDistSummary at ' + time.asctime())

                                            upPtsMbg2 = arcpy.MinimumBoundingGeometry_management(upPtsBst, opj(inm, 'up_pts_mbg2' + sfx), 'CONVEX_HULL', 'LIST', frFld)
                                            arcpy.JoinField_management(upPtsMbg2, frFld, dfs2cut4step, frFld, pntMdnSearchDistFld)
                                            df.copyfc(verbose, upPtsMbg2, sgdb)                                            

                                            upMbgBfr2 = arcpy.Buffer_analysis(upPtsMbg2, opj(inm, 'up_mbg_bfr2' + sfx), pntMdnSearchDistFld)
                                            df.copyfc(verbose, upMbgBfr2, sgdb)

                                        ## Create a connection search area by drawing buffers around upstream and downstrea points
                                            ptsMbg2 = arcpy.MinimumBoundingGeometry_management(dnPtsBst2, opj(inm, 'dp_pts_mbg2' + sfx), 'CONVEX_HULL', 'LIST', frFld)
                                            arcpy.JoinField_management(ptsMbg2, frFld, dfs2cut4step, frFld, pntMdnSearchDistFld)
                                            df.copyfc(verbose, ptsMbg2, sgdb)

                                            dnMbgBfr2 = arcpy.Buffer_analysis(ptsMbg2, opj(inm, 'dn_mbg_bfr2' + sfx), pntMdnSearchDistFld)
                                            df.copyfc(verbose, dnMbgBfr2, sgdb)
                                            log.debug('done with dnMbgBfr2 at ' + time.asctime())

                                        ## Distill the connection search area...
                                            try:
                                                allSearchIntPrelim = arcpy.Intersect_analysis([dnMbgBfr2, upMbgBfr2], opj(inm, 'all_srch_int_prlm' + sfx))
                                            except:
                                                arcpy.RepairGeometry_management(upMbgBfr2)
                                                sgdbDnMbgBfr2 = df.copyfc(True, dnMbgBfr2, sgdb)
                                                arcpy.RepairGeometry_management(sgdbDnMbgBfr2)
                                                allSearchIntPrelim = arcpy.Intersect_analysis([sgdbDnMbgBfr2, upMbgBfr2], opj(inm, 'all_srch_int_prlm' + sfx))
                                            allSearchInt = arcpy.Select_analysis(allSearchIntPrelim, opj(inm, 'all_srch_int' + sfx), frFld + ' = ' + frFld + '_1')
                                            allSearchIntDslv = arcpy.Dissolve_management(allSearchInt, opj(inm, 'all_srch_int_dslv' + sfx), frFld, multi_part = 'SINGLE_PART')
                                            df.copyfc(verbose, allSearchIntPrelim, sgdb)
                                            df.copyfc(verbose, allSearchInt, sgdb)
                                            df.copyfc(verbose, allSearchIntDslv, sgdb)

                                            searchLayerName = "srch_Lyr" + sfx
                                            goodSrchLayer = arcpy.MakeFeatureLayer_management(allSearchIntDslv, searchLayerName)
                                            goodDslvSrchLayerInDn = arcpy.SelectLayerByLocation_management(searchLayerName, 'INTERSECT', dnPtsBst2)
                                            goodDslvSrchInDn = arcpy.CopyFeatures_management(goodDslvSrchLayerInDn, opj(inm, 'gd_dslv_in_dn' + sfx))
                                            df.copyfc(verbose, goodDslvSrchInDn, sgdb)

                                            upPtsBst2Prelim = arcpy.Intersect_analysis([upPtsBst, goodDslvSrchInDn], opj(inm, 'up_pts_bst_int_srch' + sfx), output_type = "POINT")
                                            upPtsBst2PrelimFields = df.getfields(upPtsBst2Prelim)
                                            for fld in upPtsBst2PrelimFields:
                                                if fld.startswith('FID_'):
                                                    arcpy.DeleteField_management(upPtsBst2Prelim, fld)
                                            df.copyfc(verbose, upPtsBst2Prelim, sgdb)
                                            upPtsBst2 = arcpy.Select_analysis(upPtsBst2Prelim, opj(inm, 'up_pts_bst2' + sfx), frFld + ' = ' + frFld + '_1')# ' + minElFld)
                                            df.copyfc(verbose, upPtsBst2, sgdb)
                                            upPtsList.append(upPtsBst2)
                                            dnPtsList.append(dnPtsBst2)
                                            srchList.append(goodDslvSrchInDn)

                    df.tryAddField(dfs2cut4step, minSummaryFld, 'DOUBLE')
                    df.tryAddField(dfs2cut4step, revCutElFld, 'LONG')
                    df.copyfc(verbose, dfs2cut4step, sgdb)

    
        inmDfs = df.condenseDataLvls(dfsList, opj(inm, 'dfs_' + huc12))
        arcpy.DeleteField_management(inmDfs, gridfield)
        gdbDfs = df.copyfc(True, inmDfs, sgdb)

        inmDfs2Cut = df.condenseDataLvls(dfs2CutList, opj(inm, 'dfs2cut_' + huc12))
        arcpy.DeleteField_management(inmDfs2Cut, gridfield)
        gdbDfs2Cut = df.copyfc(True, inmDfs2Cut, sgdb)
    ####        SSdfs = copyfc2(inmDfs, SSF)
    ####        dfs = tno(SSdfs)

        log.debug("time check SSDfs completed at " + time.asctime())
        goodDslvAll = df.condenseDataLvls(srchList, opj(inm, 'good_int_dslv_all'))
        df.copyfc(True, goodDslvAll, sgdb)
        goodDnDslvAll = df.condenseDataLvls(dnPtsList, opj(inm, 'good_dn_pts_all'))
        df.copyfc(True, goodDnDslvAll, sgdb)
        goodUpDslvAll = df.condenseDataLvls(upPtsList, opj(inm, 'good_up_pts_all'))
        df.copyfc(True, goodUpDslvAll, sgdb)

        arcpy.env.workspace = proc_dir

        searchPickleFileObject = open(search_distance_file, 'w+b')
        pickle.dump(maxSearchDistList, searchPickleFileObject, protocol=2)
        searchPickleFileObject.close()


    except:
        # Get the traceback object
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]

        # Concatenate information together concerning the error into a message string
        pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
        msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"

        # Return python error messages for use in script tool or Python Window
        arcpy.AddError(pymsg)
        arcpy.AddError(msgs)

        # Print Python error messages for use in Python / Python Window
        log.warning(pymsg + "\n")
        log.warning(msgs)

        sys.exit(1)

    finally:
        
        if 'logName' in locals():
            log.warning("Ending script execution at " + time.asctime())
            log.warning("Script execution lasted " + str(time.time()-startTime) + " seconds or " + str((time.time()-startTime)/60) + " minutes\n")


class msgStub:
    def addMessage(self,text):
        arcpy.AddMessage(text)
    def addErrorMessage(self,text):
        arcpy.AddErrorMessage(text)
    def addWarningMessage(self,text):
        arcpy.AddWarningMessage(text)

if __name__ == "__main__":
    import sys

    if len(sys.argv) == 1:
        arcpy.AddMessage("Whoo, hoo! Running from Python Window!")
        cleanup = False

        parameters = ["C:/Program Files/ArcGIS/Pro/bin/Python/envs/arcgispro-py3/pythonw.exe",
	"C:/DEP/Scripts/basics/cmd_matcher_DEM.pyt",
	"C:/DEP/LiDAR_Current/elev_FLib_mean18/07080105/ef3m070801050901.tif",
	"C:/DEP/LiDAR_Current/elev_PLib_mean18/07080105/ep3m070801050901.tif",
	"C:/DEP/Man_Data_ACPF/dep_ACPF2022/07080105/idepACPF070801050901.gdb/buf_070801050901",
	"C:/DEP/LiDAR_Current/huc8_26915/huc_07080105.gdb/rd_rr_rd_rw_mrg_07080105",
	"C:/DEP_Proc/DEMProc/Cut_dem2013_3m_070801050901/fr0_0",
	"C:/DEP_Proc/DEMProc/Cut_dem2013_3m_070801050901/scratch.gdb/ws_polys_0",
	"C:/DEP_Proc/DEMProc/Cut_dem2013_3m_070801050901/scratch.gdb/dfs_frToPoly_0",
	"C:/DEP_Proc/DEMProc/Cut_dem2013_3m_070801050901",
	"C:/DEP_Proc/DEMProc/Cut_dem2013_3m_070801050901/search_070801050901_mean18.pkl",
	"9.0"]

        for i in parameters[2:]:
            sys.argv.append(i)
    else:
        arcpy.AddMessage("Whoo, hoo! Command-line enabled!")
        # clean up the folder after done processing
        cleanup = True

    fill_or_void_tif, punch_tif, buffered_fc, merged_medians, fr0_rasters, ws_polys, dfs_polys, proc_dir, search_distance_file, match_depth = [i for i in sys.argv[1:]]
    messages = msgStub()

    doMatcher(fill_or_void_tif, punch_tif, buffered_fc, merged_medians, fr0_rasters, ws_polys, dfs_polys, proc_dir, search_distance_file, match_depth, cleanup, messages)
    arcpy.AddMessage("Back from doMatcher!")