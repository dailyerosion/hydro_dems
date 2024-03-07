##This program was written to help remove issues in lidar DEMs caused
##by poor performance of interpolation algorithms in specific lcocations
##It is part of a series of DEM processing algorithms designed to improve
##flow of water across an elevaiton surface.
##
##Primarily written by Brian Gelder, bkgelder@iastate.edu
##Substantial help/advice from David James
##
##Separate code section created on 2019 February 26.
##for Python 3.7 and ArcGIS 10.3.1
##
## v2 2020.05.20 - added analysis of more count and intensity files to better discern water voids
## 2021.06.28 - shrinking and expanding of nodata is now controlled by shrinkAndExpand, set to 2 from 3
## 2021.07.09 - DEM processing adapted so lakes better flow to center after flattening
## 2022.06.02 - Extended bridge fix to the fill region minimum elevation and made the fix elevation fr minimum +1
## 2022.06.06 - added more error handling where mspName is created but 
# coding: utf-8
import arcpy
from arcpy.sa import *
import os
import time
import traceback
import sys
import platform
import winsound
import math
sys.path.append(os.getcwd())
# sys.path.append('O:\\DEP\\Scripts\\basics')
# import cmd_v9g_shallow_crv_plus_compact as vf

# fix for embedded code (using 'load code' in ArcGIS Pro) fails when testing len(sys.argv)
if not hasattr(sys, 'argv'):
    sys.argv = ['']
import dem_functions as df
##import dem_functionsNew as dfNew
from os.path import join as opj



class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "CMD_Flattener"
        self.alias = "CMD_Flattener"
        # List of tool classes associated with this toolbox
        self.tools = [Tool]


class Tool(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Flattener"
        self.description = "Flattens the void areas in a DEM and tries to stair step the rivers so they flow downstream"
        self.canRunInBackground = False
        self.category = "DEM Flattener"

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
            displayName="Output Flattened Elevation Model",
            datatype="DERasterDataset",
            parameterType='Required',
            direction="Output")
        
        param2 = arcpy.Parameter(
            name="vlib_metadata",
            displayName="Flattened DEM metadata template",
            datatype="DEFile",
            parameterType='Optional',
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
        doFlattener(params[0].valueAsText, params[1].valueAsText, params[2].valueAsText, params[3].valueAsText, params[4].valueAsText, params[5].valueAsText, params[6].valueAsText, cleanup, messages)
        return

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return





def smoothByShrink(regions, shrinkAmount, zoneValue):
    """Shrinks then expands regions by a given amount"""
    shrunk = Shrink(regions, shrinkAmount, zoneValue)
    expanded = Expand(shrunk, shrinkAmount, zoneValue)

    return expanded

def fence_and_flow_region(region, basins, full_dem):
    """Utilizes basins to select an area that intersects a flowing feature and
    makes every basin that touches the region part of an area that will flow"""
    partial_basins = Con(region, basins)
    full_basins = df.fullZoneByZs(partial_basins, basins)
    region_dem = Con(full_basins, full_dem)
    region_dem_fs = FocalStatistics(region_dem, NbrRectangle(5,5), 'MAXIMUM')
    region_dem_fs_plus1k = region_dem_fs + 1000
    fenced_dem = Con(IsNull(region_dem), region_dem_fs_plus1k, region_dem)
##    fenced_dem = Con(IsNull(region_dem), region_dem_fs, region_dem)
    fenced_dem_with_hole = Con(fenced_dem != fenced_dem.minimum, fenced_dem)
    filled_dem_with_hole = Fill(fenced_dem_with_hole)
    filled_no_hole = Con(IsNull(filled_dem_with_hole), fenced_dem, filled_dem_with_hole)
    filled_no_fence = Con(region_dem, filled_no_hole)

    return filled_no_fence

def startingBasinsByStreamLink(start, dem_for_cp, fd_for_cp, basin, initial_dem):
    '''Determines the starting stream link for a cost path.
    This is used for ensuring the inverted DEM correction starts with the correct
    elevation instead of skipping the first local minima.
    REVISED - 2021.07.27 - bkgelder
    Filter the basin minimum to the starting point because some basins
    along the main stream channel would send minimum values up into the
    higher basins'''

    cpForDoubleFilteredStarts = CostPath(start, dem_for_cp, fd_for_cp)
    slForDoubleFilteredStarts = StreamLink(cpForDoubleFilteredStarts, fd_for_cp)

    rgSlForDoubleFilteredStarts = RegionGroup(slForDoubleFilteredStarts, 'EIGHT')
    rgSegment = Con(basin, rgSlForDoubleFilteredStarts)
    rgForDoubleFilteredStartsWhole = df.fullZoneByZs(rgSegment, rgSlForDoubleFilteredStarts)

    fsMin = FocalStatistics(initial_dem, statistics_type = 'MINIMUM')
    doubleFilteredStartsFsMin = Con(start, fsMin)
    ##startingBasinsMinSl = ZonalStatistics(slForDoubleFilteredStartsWhole, 'Value', doubleFilteredStartsFsMin, 'MINIMUM')
    startingBasinsMinRg = ZonalStatistics(rgForDoubleFilteredStartsWhole, 'Value', doubleFilteredStartsFsMin, 'MINIMUM')

    startingBasinsMinSl = ZonalStatistics(slForDoubleFilteredStarts, 'Value', doubleFilteredStartsFsMin, 'MINIMUM')

    # check to see if the StreamLink/Region has a doubledFileteredStartsFsMin, if not, create one so it gets an elevation
    
##    log.warning('now back to your regularly scheduled programming')
    # now back to your regularly scheduled programming
    refiltered_starting_basins_min = Con(startingBasinsMinSl > startingBasinsMinRg, startingBasinsMinSl, startingBasinsMinRg)
    refiltered_starting_basins_min_replace_nulls = Con(IsNull(refiltered_starting_basins_min), startingBasinsMinRg, refiltered_starting_basins_min)

    ##    startingBasinsPartial = Con(start, basin)
    ##    startingBasins = df.fullZoneByZs(startingBasinsPartial, basin)
    ##    startingBasinsMin = ZonalStatistics(startingBasins, 'VALUE', initial_dem, 'MINIMUM')
    ##    doubleFilteredStartsMin = Con(start, startingBasinsMin)
    ##    startingBasinsMinSl = ZonalStatistics(slForDoubleFilteredStarts, 'Value', startingBasinsMin, 'MINIMUM')
    ##

####    extract1582 = ExtractByAttributes(slForDoubleFilteredStarts, 'Value = 1582')
####    rtp1582 = arcpy.conversion.RasterToPolygon(extract1582, opj(sgdb, 'ext_1582'))

    ### Doesn't work because some stream links are in two disconnected areas, that shouldn't happen
    ##slSegment = Con(basin, slForDoubleFilteredStarts)
    ##slForDoubleFilteredStartsWhole = df.fullZoneByZs(slSegment, slForDoubleFilteredStarts)

    ### Doesn't work because you can't do a count selection on a Con
    ##rgSlVar = ZonalStatistics(slForDoubleFilteredStarts, 'Value', rgSlForDoubleFilteredStarts, 'VARIETY')
    ##goofedUpStreamLinks = Con(rgSlVar > 1, rgSlForDoubleFilteredStarts)
    ##count1Goofups = Con(goofedUpStreamLinks, goofedUpStreamLinks, '', 'COUNT = 1')
    ##isnl_count1_goofups = IsNull(count1Goofups)
    ##slMaj = FocalStatistics(slForDoubleFilteredStarts, statistics_type = 'MAJORITY')
    ##fixedSl = Con(isnl_count1_goofups == 1, slForDoubleFilteredStarts, slMaj)
    ##
    ### This would work because using ExtractByMask allows the count field to move over - Argh ESRI and your .afr rasters!
    ##rgSlVarGT1 = Con(rgSlVar > 1, rgSlVar)
    ##extByVarGT1 = ExtractByMask(slForDoubleFilteredStarts, rgSlVarGT1)
    ##extRgByVarGT1 = ExtractByMask(rgSlForDoubleFilteredStarts, rgSlVarGT1)
    ##count1Goofups = Con(extRgByVarGT1, extRgByVarGT1, '', 'COUNT = 1')

    doubleFilteredStartsFsMin0 = Con(IsNull(doubleFilteredStartsFsMin), 0, doubleFilteredStartsFsMin)
    # check to see if the StreamLink/Region has a doubledFileteredStartsFsMin, if not give it one
    starting_basins_range0 = ZonalStatistics(rgForDoubleFilteredStartsWhole, 'Value', doubleFilteredStartsFsMin0, 'RANGE')
    rg_eq_range0 = Con(starting_basins_range0, rgForDoubleFilteredStartsWhole, '', 'VALUE = 0')
    flow_len_up =FlowLength(fd_for_cp, 'UPSTREAM')
    missing_rg_fl_up_min = ZonalStatistics(rg_eq_range0, 'Value', flow_len_up, 'MINIMUM')
    missing_starts = missing_rg_fl_up_min == flow_len_up
    missing_starts_fs_min = Con(missing_starts, fsMin)
    missing_starts_min = ZonalStatistics(rgForDoubleFilteredStartsWhole, 'Value', missing_starts_fs_min, 'MINIMUM')
    full_refiltered_starting_basins = CellStatistics([refiltered_starting_basins_min_replace_nulls, missing_starts_min], 'MINIMUM')

    return full_refiltered_starting_basins

    # return refiltered_starting_basins_min_replace_nulls



def startByFlowLengthDown(cost_path, flow_direction, basin):#, ProcSize):
    '''Determines the initiation point of a cost path by its flow length
    in the downstream direction.'''
    fd_on_cp = Con(cost_path, flow_direction)
    flowLenAlongCp = FlowLength(fd_on_cp)
    fsMaxFlDown = FocalStatistics(flowLenAlongCp, "RECTANGLE 3 3 CELL", 'MAXIMUM')
    fsMaxTF = flowLenAlongCp == fsMaxFlDown
    flowLenDnEqMax = Con(fsMaxTF, basin, where_clause = 'VALUE = 1')

    return flowLenDnEqMax

def startByFlowLengthUp(cost_path, flow_direction, basin):#, ProcSize):
    '''Determines the initiation point of a cost path by its flow length
    in the upstream direction.'''
    fd_on_cp = Con(cost_path, flow_direction)
    flowLenUpAlongCp = FlowLength(fd_on_cp, 'UPSTREAM')
    fsMaxFlUp = FocalStatistics(flowLenUpAlongCp, "RECTANGLE 3 3 CELL", 'MAXIMUM')
    flowLenUpEq0 = Con(flowLenUpAlongCp, basin, where_clause = 'VALUE = 0')
    flowLenUpEq0AtAllStart = Con(fsMaxFlUp, flowLenUpEq0, where_clause = 'VALUE <= ' + str(1.5*flow_direction.meanCellHeight))

    return flowLenUpEq0AtAllStart


def fixByInversionByStartingPath(ndPlus, fenceEl, invertTargetDEM, spot4Hole, ndRcls, bestestDEM, starting_path_min, log):

    ## To correct an inverted DEM, remember water must flow out the upstream ends (before inversion), not downstream!
    ## ndPlus is 2 at regions2Fix, else 1 in area to process
    fencedRegion2Fix = Pick(ndPlus, [fenceEl, invertTargetDEM])
    log.info('fencedRegion2Fix at time: ' + time.asctime())
    holeAtMin = Con(IsNull(spot4Hole) == 1, fencedRegion2Fix, '')
    
    fillNdBarrier = Fill(holeAtMin)
    log.info('fillNdBarrier at time: ' + time.asctime())

    ## Calculate the difference between the original inverted and filled inverted DEMs
    ## This should be the change to apply to the initial DEM make things flow
    fillNdDif = fillNdBarrier - fencedRegion2Fix

    # 2021.06.10 - identified issue with above approach - if hole is at a local minimum,
    # filling process will not lower to the hole elevation
    # 2021.07.14 - new corrective action - use elevation of hole to constrain correction elevation,
    # thus no elevations greater than hole elevation will be allowed to pass
    # requires grouping of all cost paths to unique, consistent zones
    correctionToApply_unfiltered = Con(ndRcls, bestestDEM - fillNdDif)

    # the flowpaths for streamlink and costpath don't always overlap
    # streamlink needs to use flow direction from cost distance to fix this
    starting_path_min_fs = FocalStatistics(starting_path_min, statistics_type = 'MINIMUM')
    
    unfiltered_correction_needs_adjustment = starting_path_min_fs < correctionToApply_unfiltered
    # above does not allow for corrections outside of starting stream link
    pathStartingNeedsCorrectionTrue = Con(unfiltered_correction_needs_adjustment, unfiltered_correction_needs_adjustment)
    log.info('correctionToApply at time: ' + time.asctime())
    isNullPathStartingCorrection = IsNull(pathStartingNeedsCorrectionTrue)
    filtered_corrections_everywhere = Con(isNullPathStartingCorrection == 0, starting_path_min_fs, correctionToApply_unfiltered)
    isnull_filtered_corrections = IsNull(filtered_corrections_everywhere)
    correctedDEM = Con(isnull_filtered_corrections, bestestDEM, filtered_corrections_everywhere)

    return correctedDEM


def doFlattener(fillTif, cntTif, cnt1rTif, surfaceElevFile, int1rMaxFile, buf_bnd, roadsFc, input_waterway, input_water, breakpolys, voidProc, voidFixTif, bigNoDataAreas, mediumNoDataAreas, cleanup, messages):
    try:
        arguments = [fillTif, cntTif, cnt1rTif, surfaceElevFile, int1rMaxFile, buf_bnd, roadsFc, input_waterway, input_water, breakpolys, voidProc, voidFixTif, bigNoDataAreas, mediumNoDataAreas, cleanup]

        for a in arguments:
            if a == arguments[0]:
                arg_str = str(a) + '\n'
            else:
                arg_str += str(a) + '\n'

        messages.addMessage("Tool: Executing with parameters:\n" + arg_str)

        huc12, huc8, proc_size = df.figureItOut(fillTif)

        if cleanup:
            # log to file only
            log, nowYmd, logName, startTime = df.setupLoggingNoCh(platform.node(), sys.argv[0], huc12)
        else:
            # log to file and console
            log, nowYmd, logName, startTime = df.setupLoggingNew(platform.node(), sys.argv[0], huc12)

            
        if not os.path.isdir(voidProc):
            os.makedirs(voidProc)

        for j in [voidFixTif, bigNoDataAreas, mediumNoDataAreas]:
            jdir = os.path.dirname(j)
            if not os.path.isdir(jdir):
                os.makedirs(jdir)

        startTime = time.time()
        log.info("Beginning execution: " + time.asctime())
        log.info("Tool: Executing with parameters:\n" + arg_str)
        messages.addMessage("Log file at " + logName)

    ## Set the environments
        arcpy.env.overwriteOutput = True
        arcpy.env.scratchWorkspace = voidProc

        sfldr = arcpy.env.scratchFolder
        sgdb = arcpy.env.scratchGDB
        arcpy.env.scratchWorkspace = sgdb#sfldr#
        arcpy.env.workspace = sgdb

        arcpy.env.snapRaster = fillTif#snapRaster

        arcpy.env.cellSize = proc_size

        gdb = sgdb + "\\"
        cp = voidProc + "\\"
        inm = 'in_memory\\'

    ####------------------------------------------------------------------------------

        log.info('log file is ' + logName)
        log.debug('starting raster copy')
    ## Copy the Pitfilled and Punched DEMs to the workspace (not needed, just better for display than TIFFs)
        pitFilledDEM = Raster(arcpy.CopyRaster_management(fillTif, opj(gdb, os.path.splitext(os.path.basename(fillTif))[0])))
        arcpy.env.extent = pitFilledDEM.extent

        bestestDEM = pitFilledDEM

        rect3x3Nbr = NbrRectangle(3, 3, 'CELL')
        circle3Nbr = NbrCircle(3, 'CELL')

        mspName = 'main_stem_poly'

    ####------------------------------------------------------------------------------
    ## Create a layer of WBD boundary to buffer and clip datasets with
        # huc12Lyr = arcpy.MakeFeatureLayer_management(huc12Fc, 'HUC12FCLayer', '"HUC12" = \'' + huc12 + "'")
        # bnd = arcpy.CopyFeatures_management(huc12Lyr, opj(gdb, "bnd_" + huc12))
        Clip = arcpy.CopyFeatures_management(buf_bnd, opj(gdb, 'buf_' + huc12 + '_1km'))#, '1000 METER')
        # buffer buffered boundary another 4km
        Clip5K = arcpy.Buffer_analysis(Clip, opj(gdb, 'buf_' + huc12 + '_5km'), '4000 METER')

        # waterways = arcpy.Clip_analysis(input_waterway, Clip, opj(gdb, 'waterways'))
        water = arcpy.Clip_analysis(input_water, Clip, opj(gdb, 'water'))
        water5K = arcpy.Clip_analysis(input_water, Clip5K, opj(gdb, 'water5K'))
        if arcpy.Exists(breakpolys):
            breaks = arcpy.Clip_analysis(breakpolys, Clip, opj(gdb, 'break_polys'))
        else:
            breaks = None

    ######------------------------------------------------------------------------------

    ## Create a DEM in meters to keep frrom using z-factor which calculates in script directory
        meterPfDEM = 0.01 * bestestDEM

    ## Calculate curvature to find areas where curvature is positive (channels)
        crv = Curvature(meterPfDEM, '', opj(cp, "pro_crv"), opj(cp, "pln_crv"))
        crv.save(opj(cp, "crv"))

        proCrv = Raster(opj(cp, "pro_crv"))
    ######------------------------------------------------------------------------------

        ndFixedList = []
    ## Define continuous areas where no ground returns were received
    # optional count raster
        if arcpy.Exists(cntTif) == True:
            log.debug('found count raster: ' + cntTif)
            class2CountRaw = Raster(cntTif)
            class2Count0Ws = ExtractByMask(Con(IsNull(class2CountRaw), 0, class2CountRaw), Clip)
            class2Count0Ws.save(opj(voidProc, "cntbe" + str(int(proc_size)) + "m0"))

            firstCountRaw = Raster(cnt1rTif)
            firstCount0Ws = ExtractByMask(Con(IsNull(firstCountRaw), 0, firstCountRaw), Clip)
            firstCount0Ws.save(opj(voidProc, "cnt1r" + str(int(proc_size)) + "m0"))

            cntRatio = 0.9
            log.debug('count ratio is ' + str(cntRatio))

            if cntRatio > 0.8:
                
                copyCount = arcpy.CopyRaster_management(class2Count0Ws, gdb + class2Count0Ws.name)

                arcpy.env.workspace = sgdb#voidProc

                if class2Count0Ws.minimum == 0.0:
                ## fix initial DEM in areas of No Data
                    ## first fix large ND areas that touch edge
                        ## make them run off edge
                    ## then small ND areas connected to the above
                        ## make them run off edge
                    ## then large, round ND areas that don't touch edge (lakes, ponds, lagoons)
                        ## make them flat
                    ## then small, thin ND areas
                        ## make them run to lowest point in ND

                    ndTF = Con(class2Count0Ws == 0, 1, 0)
                    count1or0 = Con(class2Count0Ws == 0, 0, 1)
                    fstCount1Sum = FocalStatistics(count1or0, rect3x3Nbr, 'SUM')
                    # remove single or double 0 count cells
                    filtered_count = Con(fstCount1Sum < 7, count1or0, 1)
                    # remove single or double 1 count cells
                    dbl_filtered_count = Con(fstCount1Sum > 2, filtered_count, 0)

                    ndAs1 = Con(dbl_filtered_count == 0, 1)
                    thin_no_data = Shrink(ndAs1, 1, 1)
                    # thin_expand = Expand(thin_no_data, 1, 1)
                    thin_expand2 = Expand(thin_no_data, 2, 1)
                    back_to_normal = Shrink(thin_expand2, 1, 1)
                    # core_no_data = Shrink(thin_no_data, 2, 1)
                    ndRegionsPre = RegionGroup(back_to_normal, 'EIGHT')
                    ndZst = ZonalStatistics(ndRegionsPre, 'VALUE', class2Count0Ws, 'MEAN')#firstCount0Ws, 'MEAN')
                    ndRegions = Con(ndZst < 0.25, ndRegionsPre)

    #             ## These are the no data areas we want to work on 
    #                 ndRegions = Con(nd2WorkTF > 0, RegionGroup(Con(nd2WorkTF, 1), 'EIGHT'))
                    ndRegions.save(cp + 'nd_rgns')

                    ndCleanRegions = ndRegions

                    arcpy.CopyRaster_management(ndCleanRegions, bigNoDataAreas)

                    
                    pitFilledDEM_fd = FlowDirection(pitFilledDEM)
                    pitFilledDEM_basins = Basin(pitFilledDEM_fd)
                    # pf_basins_rtp = arcpy.conversion.RasterToPolygon(pitFilledDEM_basins, opj(sgdb, 'pf_basins'), 'NO_SIMPLIFY')
                    # initialWs = pitFilledDEM_basins
                    # basin_min = ZonalStatistics(pitFilledDEM_basins, 'Value', pitFilledDEM, 'MINIMUM')
                    # deepest_cells_el = Con(basin_min == pitFilledDEM, pitFilledDEM)
                    
                    rDepth, MaxDepth, fillPrevLvl = df.findMaxDepth(pitFilledDEM)

                ## Find ND regions that intersect HUC 12 boundary and move in a couple cells to avoid edge issues
                    bndBufLine = arcpy.PolygonToLine_management(Clip, inm + 'bnd_buf_line')
                    # df.copyfc(verbose, bndBufLine, gdb)
                    bndBufRasterLinePre = Raster(arcpy.PolylineToRaster_conversion(bndBufLine, 'OBJECTID', opj(cp, 'buf_raster'), cellsize = str(proc_size)))
                    bndBufRasterLine = Expand(bndBufRasterLinePre, 1, bndBufRasterLinePre.maximum)

                ## see how thick the no data is and use maximum of 1 cell or nodata thickness to expand 
                    ndThickness = ZonalGeometry(ndRegions, 'VALUE', 'THICKNESS')
                    bndBufNdThick = Con(bndBufRasterLine, ndThickness)

                    if bndBufNdThick.maximum >= 2:
                        numCells = int(bndBufNdThick.maximum/2)
                    else:
                        numCells = 1
                    bndBufRasterLine3 = Expand(bndBufRasterLinePre, numCells, bndBufRasterLinePre.maximum)
                    
                    addInBnd = Con(pitFilledDEM, IsNull(bndBufRasterLine) + IsNull(bndBufRasterLine3))
                    nearBndNd = Con(addInBnd == 1, ndRegions)#goodCrvNdRegions)#ndRegions)

                ## Check to see if there are any potential rivers to enforce (from void analysis)
                    ptlRiversOnBnd = ZonalStatistics(ndCleanRegions, 'value', nearBndNd, 'MAXIMUM')

                    log.debug('here! at ' + time.asctime())
                    if ptlRiversOnBnd.maximum is not None:# > 0.0:
                        if ptlRiversOnBnd.maximum > 0.0:

                            # if True:
                            try:
                            ## Only enforce the 'big' rivers (more than 2 acres (8400 m2))
                                bigRivers = Con(ptlRiversOnBnd, ndCleanRegions, '', 'COUNT > ' + str(int(8400/(ndCleanRegions.meanCellHeight*ndCleanRegions.meanCellWidth))))#1000')

                                if bigRivers.maximum is not None:
                                    if bigRivers.maximum > 0.0:

                                        bigRiversArea = ZonalGeometry(bigRivers, 'Value', 'AREA')
                                        bigRiver = Con(bigRiversArea == bigRiversArea.maximum, bigRivers)

                                        # bring in cells that border the edge of the NoData river
                                        pf_fs5_min_el = FocalStatistics(bestestDEM, NbrRectangle(5,5), 'MINIMUM')
                                        ## Smooth the original DEM in areas of no data to create a focal minimum DEM
                                        bigRiverDEM = Con(bigRiver, pf_fs5_min_el)#bestestDEM)

                                        ndThRivers = Con(bigRiver, ndThickness)
                        
                                        ndSearchNbr = NbrCircle(min([10*proc_size, int(ndThRivers.maximum)]), 'CELL')
                                        
                                        fsNdWidthMin = FocalStatistics(bigRiverDEM, ndSearchNbr, 'MINIMUM')#"RECTANGLE " + ndThickCrit + " "  + ndThickCrit + " CELL", "MINIMUM")
                                        fs5_nd_maj = FocalStatistics(fsNdWidthMin, NbrCircle(5), 'MAJORITY')
                                        bc_fs5_nd_maj = BoundaryClean(fs5_nd_maj, 'DESCEND')
                                        
                                        ## Find the lowest edge (so we can force water flow to it)
                                        ndEndBuffers = RegionGroup(Con(nearBndNd, bigRiver), 'EIGHT')
                                        ndEndBuffersMin = ZonalStatistics(ndEndBuffers, 'value', bc_fs5_nd_maj, 'MINIMUM')#fs15CircleMin, 'MINIMUM')
                                        lowestBndBuf = Con(ndEndBuffersMin == ndEndBuffersMin.minimum, ndEndBuffers)
                                        
                                        ## Find the other end buffers
                                        otherNdEndBuffers = Con(IsNull(lowestBndBuf) == 1, ndEndBuffers)

                                        # make the big river main stem and branches flow (but not cleaned)
                                        big_river_nd_dem = Con(IsNull(bigRiver), pitFilledDEM, bc_fs5_nd_maj)
                                        big_river_nd_fd = FlowDirection(big_river_nd_dem)
                                        big_river_nd_basins = Basin(big_river_nd_fd)
                                        big_river_nd_flow = fence_and_flow_region(bigRiver, big_river_nd_basins, big_river_nd_dem)
                                        big_river_nd_flow_fd = FlowDirection(big_river_nd_flow)
                                        # big_river_nd_flow_fa = FlowAccumulation(big_river_nd_flow_fd)

                                        # run a downstream trace to separate mainstem and branches
                                        cp_other_nd_down = CostPath(otherNdEndBuffers, big_river_nd_flow, big_river_nd_flow_fd)
                                        distToMainStem = CostDistance(cp_other_nd_down, Con(bigRiver, 1))

                                        ## Define the main stem of the big river
                                        distance_threshold = (ndThRivers.mean + ndThRivers.maximum) / 2 * 2
        ##                                mainStem = Con(ndThRivers < distance_threshold, bigRiver)
                                        mainStem = Con(distToMainStem < distance_threshold, bigRiver)
                                        if mainStem.maximum is not None:#False
                                            mainStemPoly = arcpy.RasterToPolygon_conversion(mainStem, gdb + mspName)

                                            # now make the main stem flow
                                            main_stem_nd_dem = Con(IsNull(mainStem), pitFilledDEM, bc_fs5_nd_maj)
                                            just_main_stem_nd_dem = Con(mainStem, bc_fs5_nd_maj)
                                            main_stem_nd_fd = FlowDirection(main_stem_nd_dem)
                                            main_stem_nd_basins = Basin(main_stem_nd_fd)
                                            main_stem_nd_flow = fence_and_flow_region(mainStem, main_stem_nd_basins, main_stem_nd_dem)
                                            main_stem_nd_flow_fd = FlowDirection(main_stem_nd_flow)
            ####                                main_stem_nd_flow_basins = Basin(main_stem_nd_flow_fd)
            ####                                main_stem_nd_flow_basins_rtp = arcpy.conversion.RasterToPolygon(main_stem_nd_flow_basins, opj(sgdb, 'main_stem_nd_flow_basins'), 'NO_SIMPLIFY')
                                            main_stem_nd_basins_rtp = arcpy.conversion.RasterToPolygon(main_stem_nd_basins, opj(sgdb, 'main_stem_nd_basins'), 'NO_SIMPLIFY')

                                            basin0Deepest = ZonalStatistics(main_stem_nd_basins, 'VALUE', main_stem_nd_dem, 'MINIMUM')
                                            basin0DeepestCell = Con(basin0Deepest == main_stem_nd_dem, main_stem_nd_dem)

                                            all_other_nd_end_basins = Con(otherNdEndBuffers, main_stem_nd_basins)
                                            all_other_nd_area = ZonalGeometry(all_other_nd_end_basins, 'Value', 'AREA')
                                            max_basin_area_each_other_nd_end = ZonalStatistics(otherNdEndBuffers, 'VALUE', all_other_nd_area, 'MAXIMUM')
                                            biggest_basin_in_each_end = Con(all_other_nd_area == max_basin_area_each_other_nd_end, main_stem_nd_basins)
                                            recover_full_biggest = df.fullZoneByZs(biggest_basin_in_each_end, main_stem_nd_basins)
                                            deepest_in_full_biggest = Con(recover_full_biggest, basin0DeepestCell)

                                            fl_down = FlowLength(main_stem_nd_flow_fd)
                                            fl_down_in_deepest = Con(deepest_in_full_biggest, fl_down)
                                            basin_max_fld = ZonalStatistics(recover_full_biggest, 'VALUE', fl_down_in_deepest, 'MAXIMUM')
                                            eq_max_fld = Con(basin_max_fld == fl_down_in_deepest, 1)

                                            cp_deepest_max_down = CostPath(eq_max_fld, main_stem_nd_flow, main_stem_nd_flow_fd)

                                            # add nibble here to only use elevations on the cost path???
                                            # use nibble to replace mainstem elevation values with nearest value from cost path
                                            # but this still leaves a few dead ends we'll fix in a bit
                                            cp_deepest_nd_el = Con(cp_deepest_max_down, bc_fs5_nd_maj)
                                            nibble_nd = Nibble(bc_fs5_nd_maj, cp_deepest_nd_el)

                                            # now go back and fix some potential errors in the above DEM
                                            # the fixedDEM_flatter rivers sometimes dead end in a low spot since the ND focal stats fixing process
                                            # can create separate 'branches' of lower 'fixed' elevations and if the flow direction algorithm decides to choose the wrong one it gets stuck there
                                            # just_main_stem_nd_dem_regions = RegionGroup(just_main_stem_nd_dem, 'EIGHT')
                                            # just_main_stem_nibble_nd_dem = Con(mainStem, nibble_nd)
                                            # just_main_stem_nibble_nd_dem_regions = RegionGroup(just_main_stem_nibble_nd_dem, 'EIGHT')
                                            # bc_fs5_nd_maj_regions = RegionGroup(bc_fs5_nd_maj, 'EIGHT')
                                            
                                            # just_main_stem_nd_el_regions_on_cp = Con(cp_deepest_max_down, just_main_stem_nibble_nd_dem_regions)#double_down
                                            # cp_just_main_region_in_bc_region = ZonalStatistics(bc_fs5_nd_maj_regions, 'Value', just_main_stem_nd_el_regions_on_cp, 'MAJORITY')
                                            # full_ms_el_regions_not_eq_just_main_stem_regions = NotEqual(cp_just_main_region_in_bc_region, just_main_stem_nibble_nd_dem_regions)
                                            # flowing_el_in_not_eq_regions = Con(full_ms_el_regions_not_eq_just_main_stem_regions, main_stem_nd_flow)
                                            # fixed_main_stem_dem = Con(IsNull(flowing_el_in_not_eq_regions), main_stem_nd_dem, flowing_el_in_not_eq_regions)
                                            fixed_main_stem_dem = Con(IsNull(mainStem), main_stem_nd_dem, nibble_nd)
                                            # fixed_main_stem_dem2 = Con(IsNull(flowing_el_in_not_eq_regions), main_stem_nd_dem, nibble_nd)
                                            fixed_main_stem_fd = FlowDirection(fixed_main_stem_dem)
                                            fixed_main_stem_basins = Basin(fixed_main_stem_fd)
                                            

                                            fixed_main_stem_flow = fence_and_flow_region(mainStem, fixed_main_stem_basins, fixed_main_stem_dem)
                                            fixed_main_stem_flow_fd = FlowDirection(fixed_main_stem_flow)
                                            fixed_main_stem_flow_fa = FlowAccumulation(fixed_main_stem_flow_fd)


                                            # clean out the main stem
                                            filtered_down = startByFlowLengthDown(cp_deepest_max_down, fixed_main_stem_flow_fd, recover_full_biggest)
                                            filtered_up = startByFlowLengthUp(cp_deepest_max_down, fixed_main_stem_flow_fd, recover_full_biggest)
                                            # filtered_down = vf.startByFlowLengthDown(cp_deepest_max_down, fixed_main_stem_flow_fd, recover_full_biggest)
                                            # filtered_up = vf.startByFlowLengthUp(cp_deepest_max_down, fixed_main_stem_flow_fd, recover_full_biggest)
                                            double_filtered = Con(filtered_down, filtered_up)
                                            double_filtered_regions = RegionGroup(double_filtered, 'EIGHT')
                                            max_dbl_filter_regions = ZonalStatistics(double_filtered, 'VALUE', double_filtered_regions, 'MAXIMUM')
                                            eq_max_dbl_filter = Con(double_filtered_regions == max_dbl_filter_regions, double_filtered_regions)
                                            cp_max_down = CostPath(eq_max_dbl_filter, fixed_main_stem_flow, fixed_main_stem_flow_fd)


                                            sl_cp_for_main_stem = startingBasinsByStreamLink(eq_max_dbl_filter, fixed_main_stem_flow, fixed_main_stem_flow_fd, fixed_main_stem_basins, fixed_main_stem_dem)

                                            # can't use flowing main stem DEM above becuase it was made to flow and wipes out many low spots that need cleaning
                                            invertTargetDEM = df.createInvertDEM(fixed_main_stem_dem)
                                            step_ndRcls_flatter, step_expFix_flatter, step_ndPlus_flatter = df.setupInversionZones(cp_max_down)#cp_double_down)
                                            fence_el_flatter = df.getFenceEl(invertTargetDEM)
                                            # this DEM only is fixed on the cost path
                                            fixedDEM_flatter = fixByInversionByStartingPath(step_ndPlus_flatter, fence_el_flatter, invertTargetDEM, double_filtered, step_ndRcls_flatter, fixed_main_stem_dem, sl_cp_for_main_stem, log)
                                            # fixedDEM_flatter = vf.fixByInversionByStartingPath(step_ndPlus_flatter, fence_el_flatter, invertTargetDEM, double_filtered, step_ndRcls_flatter, fixed_main_stem_dem, sl_cp_for_main_stem, log)
                                            fixedDEM_flatter_fd = FlowDirection(fixedDEM_flatter)
                                            fixedDEM_flatter_fa = FlowAccumulation(fixedDEM_flatter_fd)
                                            # use nibble to replace mainstem elevation values with nearest value from cost path
                                            sl_cp_for_main_stem_dem = Con(sl_cp_for_main_stem, fixedDEM_flatter)
                                            nibble_main = Nibble(fixedDEM_flatter, sl_cp_for_main_stem_dem)
                                            isn_big_river = IsNull(bigRiver)
                                            flatter_mainstem = Con(isn_big_river, pitFilledDEM, nibble_main)

                                            # the above DEM can still go into some dead end low spots/disconnections due to the masking
                                            fill_fixedDEM_flatter = Fill(fixedDEM_flatter)
                                            flowing_mainstem = Con(IsNull(fixed_main_stem_flow), pitFilledDEM, fill_fixedDEM_flatter)
                                            # fill_fd = FlowDirection(flowing_mainstem)
                                            # fill_fa = FlowAccumulation(fill_fd)

                                            ## Find the branches of the big river, don't clean them yet
                                            distant_branches = Con(distToMainStem >= distance_threshold, bigRiver)
                                            rg_distant_branches = RegionGroup(distant_branches, 'EIGHT')
                                            big_distant_branches_pre = Con(rg_distant_branches, rg_distant_branches, '', 'COUNT >= 5')

                                            fstNdSum = FocalStatistics(ndTF, rect3x3Nbr, 'SUM')
                                            nd2WorkTF = Con(pitFilledDEM, Con(fstNdSum > 1, ndTF, 0))#class2Count0Ws == 0, 1, 0)
                                            ndExpand = Expand(Con(nd2WorkTF, 1), 1, 1)#nd2WorkTF, 1, 1)
                                            ndRegionsExp = RegionGroup(ndExpand, 'EIGHT')

                                            branches_starter = Con(big_distant_branches_pre, ndRegionsExp)
                                            full_starter = df.fullZoneByZs(branches_starter, ndRegionsExp)
                                            mainStemExpand = Expand(mainStem, 1, mainStem.maximum)
                                            full_starter_not_main = Con(IsNull(mainStemExpand), full_starter)
                                            rg_full_starter = RegionGroup(full_starter_not_main, 'EIGHT')
                                            full_starter_with_distant = ZonalStatistics(rg_full_starter, 'Value', big_distant_branches_pre, 'MAXIMUM')
                                            distToMainStem2 = CostDistance(cp_other_nd_down, Con(full_starter, 1))

                                            ndRegions_branches_partials = Con(full_starter_with_distant, ndRegions)
                                            expand_branches_partials = Expand(Con(ndRegions_branches_partials, 1), 1, 1)#rg_branches_partials.maximum)
                                            expand_branches = RegionGroup(expand_branches_partials, 'EIGHT')


            ####                                branches_table = arcpy.CopyRows_management(big_distant_branches, opj(gdb, 'branches_tbl'))
            ####                                branches_values = [row[0] for row in arcpy.da.SearchCursor(branches_table, ['VALUE'])]
            ####
            ####                                # find deepest cells at max distance from main stem
            ####                                expand_branches = Expand(big_distant_branches, 1, branches_values)
                                            fs_dist = FocalStatistics(distToMainStem2, statistics_type = 'MEAN')
                                            expand_deepest_dist = Con(expand_branches, Con(basin0DeepestCell, fs_dist))

                                            expand_branches_max_dist = ZonalStatistics(expand_branches, 'Value', expand_deepest_dist, 'MAXIMUM')
                                            eq_branches_max_dist = Con(expand_branches_max_dist == expand_deepest_dist, basin0DeepestCell)

                                            path_dem = Con(CellStatistics([expand_branches, mainStem]), flowing_mainstem)#pitFilledDEM)

                                            branches_thickness = ZonalGeometry(big_distant_branches_pre, 'Value', 'THICKNESS')
                                            branch_fs_dist = int(((branches_thickness.mean / pitFilledDEM.meanCellHeight) * 2 - 1)//2*2+1)
                                            fs_min_el_expand_area = FocalStatistics(path_dem, NbrRectangle(branch_fs_dist, branch_fs_dist), 'MINIMUM')
                                            # trim the focal statistics area to the no data plus a couple cells to constrain resulting cost paths
                                            trim_area = CellStatistics([mainStem, expand_branches])
                                            fs_min_el_expand_trim = Con(trim_area, fs_min_el_expand_area)
                                            dif_from_fs_min_trim = flowing_mainstem - fs_min_el_expand_trim
                                            difPlus1_branches = dif_from_fs_min_trim + 1

                                            cost_dist_from_deepest_neg = CostDistance(lowestBndBuf, difPlus1_branches, out_backlink_raster = opj(sgdb, 'back_lnk'))
                                            cost_dist_from_deepest = Con(cost_dist_from_deepest_neg < 1, 1, cost_dist_from_deepest_neg)
                                            cp_favor_deeper = CostPath(eq_branches_max_dist, cost_dist_from_deepest, opj(sgdb, 'back_lnk'))

                                            cp_favor_deeper_all = cp_favor_deeper

                                            cp_to_zonal = Con(eq_branches_max_dist, cp_favor_deeper_all)
                                            if cp_to_zonal.minimum > 2:
                                                cp_to_recode = df.fullZoneByZs(cp_to_zonal, cp_favor_deeper_all)
                                                cp_minimum = ZonalStatistics(cp_to_recode, 'Value', eq_branches_max_dist, 'MINIMUM')
                                            else:
                                                log.warning('UHOH - branches not unique, elevation recoding will not go well!')

                                            desired_cp = Con(eq_branches_max_dist, cp_favor_deeper_all)
                                            full_desired_cp = df.fullZoneByZs(desired_cp, cp_favor_deeper_all)
                                            full_desired_cp_el = ZonalStatistics(full_desired_cp, 'Value', eq_branches_max_dist, 'MINIMUM')

                                            invertTargetDEM_branches = df.createInvertDEM(flowing_mainstem)#RiverDEM)#filled_fence_no_hole)#
                                            step_ndRcls_branches, step_expFix_branches, step_ndPlus_branches = df.setupInversionZones(full_desired_cp)#cp_favor_deeper_all)#cp_hole)
                                            fence_el_branches = df.getFenceEl(invertTargetDEM_branches)

                                            fixedDEM_branches = fixByInversionByStartingPath(step_ndPlus_branches, fence_el_branches, invertTargetDEM_branches, eq_branches_max_dist, step_ndRcls_branches, flowing_mainstem, full_desired_cp_el, log)
                                            # fixedDEM_branches = vf.fixByInversionByStartingPath(step_ndPlus_branches, fence_el_branches, invertTargetDEM_branches, eq_branches_max_dist, step_ndRcls_branches, flowing_mainstem, full_desired_cp_el, log)
                                            branch_new_el = Con(full_desired_cp_el, fixedDEM_branches)
                                            flatter_just_mainstem = Con(isn_big_river == 0, flowing_mainstem)

                                            flowing_mainstem_with_branches = Con(IsNull(branch_new_el), flatter_just_mainstem, branch_new_el)
                                            flowing_mainstem_with_branches_fs =FocalStatistics(flowing_mainstem_with_branches, NbrRectangle(5,5), 'MAXIMUM')
                                            flatter_branches_max_1k = flowing_mainstem_with_branches_fs + 1000
                                            flatter_branches_fenced = Con(IsNull(flowing_mainstem_with_branches), flatter_branches_max_1k, flowing_mainstem_with_branches)
                                            flatter_just_mainstem_with_hole = Con(flatter_branches_fenced != flatter_branches_fenced.minimum, flatter_branches_fenced)
                                            filled_flatter = Fill(flatter_just_mainstem_with_hole)

                                            fixedDEM_branches_and_main = Con(IsNull(flowing_mainstem_with_branches), pitFilledDEM, filled_flatter)
                                            flatter_branches_fd = FlowDirection(fixedDEM_branches_and_main)
                                            flatter_branches_basins = Basin(flatter_branches_fd)

                                            flatter_fence_and_flow = fence_and_flow_region(flowing_mainstem_with_branches, flatter_branches_basins, fixedDEM_branches)
                                            flatterRiverDEM = Con(IsNull(flatter_fence_and_flow), fixedDEM_branches, flatter_fence_and_flow)
                                            fixedDEM_branches_flowing_fd = FlowDirection(flatterRiverDEM)
                                            fixedDEM_branches_flowing_basins = Basin(fixedDEM_branches_flowing_fd)

                                            flatter_branches_basins_rtp = arcpy.conversion.RasterToPolygon(fixedDEM_branches_flowing_basins, opj(sgdb, 'flatter_branches_basins'), 'NO_SIMPLIFY')
            ##                                flatter_brandhes_fa = FlowAccumulation(flatter_branches_fd)

                                    else:
                                        log.debug('no big rivers to clean - NoData must be big enough and on borders')

                                else:
                                    log.debug('no big rivers at all - NoData must be big enough and on borders')

                            except:
        ##                    except Exception as err:
        ##                        log.debug('fail on big river clean')
        ##                        log.debug(err.message)
        ##                        arcpy.AddError(err.message)

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

                if arcpy.Exists(gdb + mspName) and 'flatterRiverDEM' in locals():
                    ndCleanNotConnectedPre = Con(IsNull(mainStem), ndCleanRegions)
                
                    # flatterRiverDEM is not always created, e.g. failure after creation due to empty zones in 102702040206
                else:
                    log.warning('No rivers found, ndCleanNotConnected not successfully created, pfDEM becoming flatterRiverDEM for further No Data processing')
                    ndCleanNotConnectedPre = ndCleanRegions
                    flatterRiverDEM = pitFilledDEM

            ## Now fix initial DEM in large No Data areas not connected to main stem (i.e. lakes, ponds, lagoons)

                # test the OSM water polygons to see if they contained water
                water_only_polys = arcpy.analysis.Select(water, 'just_water', "fclass = 'water' OR fclass = 'reservoir'")
                water_only_desc = arcpy.da.Describe(water_only_polys)
                arcpy.management.DeleteIdentical(water_only_polys, "Shape")
                # eliminate duplicate OSM_ids
                if df.testForZero(water_only_polys):
                    water_to_flatten_list = []
                    
                    water_only_shrink = arcpy.analysis.Buffer(water_only_polys, buffer_distance_or_field = '-200 METERS')
                    # test to see if there are big water areas from OSM, use these
                    if df.testForZero(water_only_shrink):
                        f = []
                        with arcpy.da.SearchCursor(water_only_shrink, [water_only_desc['OIDFieldName'], 'osm_id'], sql_clause = ('DISTINCT osm_id', None)) as scur:
                            for srow in scur:
                                f.append(srow[1])
                        # c = df.getFrsAsList(water_only_polys, 'osm_id', '')
                        # big_water = df.selectByList(f, water_only_desc['OIDFieldName'], water5K, 'big_water_5k', '_0', inm)
                        big_water = df.selectByList(f, 'osm_id', water5K, 'big_water_5k', '_0', inm)

                        # find 'large' lakes/reservoirs that we will flatten
                        water_only_polys5K = arcpy.analysis.Select(big_water, 'just_water5K', "fclass = 'water' OR fclass = 'reservoir'")
                        
                        # if df.testForZero(water_only_shrink200):
                            # interior_water = arcpy.conversion.PolygonToRaster(water_only_shrink200_clip, water_only_desc['OIDFieldName'])
                            # b = [a[0] for a in arcpy.da.SearchCursor(water_only_polys, water_only_desc['OIDFieldName'])]
                            # water_only_shrink200_clip = arcpy.Clip_analysis(d, Clip5K)
                        water_only5K_shrink200 = arcpy.analysis.Buffer(water_only_polys5K, buffer_distance_or_field = '-200 METERS')
                            # int_dif = Raster(int1rMaxFile) - Raster(int1rMinFile)
                        #     edge_water = 
                        # use edge detection from int1rMaxFile to find waterline break
                        water_only_raster = arcpy.conversion.PolygonToRaster(water_only5K_shrink200, water_only_desc['OIDFieldName'])
                        # water_only_shrink_std = ZonalStatistics(water_only_raster, 'Value', pitFilledDEM, 'STD')
                        # water_only_shrink_rdepth = ZonalStatistics(water_only_raster, 'Value', rDepth, 'MEAN')
                        water_only_count = Con(water_only_raster, class2Count0Ws)
                        water_only_countTF = Con(water_only_count > 0, 1, 0)
                        try:
                            # water_only_int = Con(water_only_raster, int1rMaxFile)

                            # redo the NoData zones in the big water area, use focal stats to remove 'rogue' returns'
                            # count_width = 5
                            # fs5_cnt_sum = FocalStatistics(class2Count0Ws, NbrRectangle(count_width, count_width), 'SUM')
                            # fs5_cnt_voids = Con(fs5_cnt_sum < (0.2*count_width**2), 0, class2Count0Ws)
                            # fs5_cnt0 = Con(fs5_cnt_voids == 0, 1)

                            # shrinkAndExpand = 1
                            # shrinkBigWater = Shrink(fs5_cnt0, shrinkAndExpand, 1)
                            # expandBigWater = Expand(shrinkBigWater, shrinkAndExpand, 1)

                            # rg_fs5_cnt0 = RegionGroup(expandBigWater, 'EIGHT')
                            # partial_cnt0_rg = Con(water_only_raster, expandBigWater)
                            # full_cnt0_big_water = df.fullZoneByZs(partial_cnt0_rg, expandBigWater)#rg_fs5_cnt0)

                            # based on an analysis of OSM water and wetland polygons on HUC102500091002, this separates water from non-water well
                            water_to_flatten_osm = Con(water_only_countTF.mean > 0.4, water_only_raster)
                            # water_to_flatten_osm.save('flat_osm')

                            partial_nd_osm = Con(water_to_flatten_osm, ndCleanNotConnectedPre)
                            full_nd_osm = df.fullZoneByZs(partial_nd_osm, ndCleanNotConnectedPre)

                            water_to_flatten_max_value = int(ndCleanNotConnectedPre.maximum + 1)
                            cell_stat_nd = CellStatistics([full_nd_osm, water_to_flatten_osm], 'MAXIMUM')
                            water_to_flatten_osm_value = Con(cell_stat_nd, water_to_flatten_max_value)
                            water_to_flatten_osm_value.save('flat_osm')
                            water_to_flatten_list.append(water_to_flatten_osm_value)

                            intRaster = Raster(int1rMaxFile)
                            high_int = intRaster < 32
                            low_int = intRaster > 248
                            low_high_int = low_int | high_int

                            small_nd = Con(IsNull(water_to_flatten_osm_value), ndCleanNotConnectedPre)

                            # use these regions to aggregate the original basins. Analyze the edge elevations
                            # lagoons should have fairly even edge elevations all the way around, unlike other basins
                            # edge_shrunk = 

                        #     water_only_polys_not_large = df.antiSelectByList(f, water_only_desc['OIDFieldName'], water, 'little_water_5k', '_0', inm)
                        except:
                            small_nd = ndCleanNotConnectedPre

                    else:
                        small_nd = ndCleanNotConnectedPre
                    #     water_only_polys_not_large = arcpy.CopyFeatures_management(water_only_polys)

                    log.debug('analyzing for small voids/open water areas at ' + time.asctime())
                    # maxNdTh = ZonalStatistics(small_nd, 'VALUE', ndThickness, 'MAXIMUM')
                    #refine ndCleanNotConnected by curvature, intensity
                    ## if high max thickness and mostly flooded to flow, keep
                    pondFillCrit = 0.8
                    rDepthPos = ExtractByAttributes(rDepth, 'VALUE > 0')
                    rDepthPosAs1 = Con(rDepthPos, 1)
                    # rDepthPosExpand = Expand(rDepthPos, 7, 1)
                    # ptlPonds = Con(maxNdTh > proc_size*2.5, Con(rDepthPosExpand, small_nd))

                    # if df.testForZero(ptlPonds):
                        # pondCountRatio = Lookup(ptlPonds, 'COUNT')*1.0/Lookup(small_nd, 'COUNT')
                        # ponds2Flatten = Con(pondCountRatio > pondFillCrit, small_nd)#ndCleanNotConnectedPre)
                        # ponds2FlattenAs1 = Con(ponds2Flatten, 1)
                        # pondsShrink = smoothByShrink(ponds2FlattenAs1, 2, 1)
                        # pondsShrinkTF = IsNull(pondsShrink)
                        # rgPondsShrinkTF = RegionGroup(pondsShrinkTF, 'EIGHT')
                        # pondsShrink_no_small_voids = Con(rgPondsShrinkTF, 1, '', 'LINK = 0 or (LINK = 1 AND Count < 100)')
                        # pondsShrinkRg = RegionGroup(pondsShrink_no_small_voids, 'EIGHT')
    ####                    pondsShrinkRg = RegionGroup(pondsShrink, 'EIGHT')
                        # log.debug('pondsShrinkRg here! at ' + time.asctime())

                hag = surfaceElevFile - pitFilledDEM
                # ndCleanNotConnected = CellStatistics([ndCleanNotConnectedPre, water_to_flatten], 'MAXIMUM')
                if len(water_to_flatten_list) > 0:
                    water_to_flatten_rg = CellStatistics(water_to_flatten_list, 'MAXIMUM')

                    ponds2smoothRg = water_to_flatten_rg

                    if ponds2smoothRg is not False:
                        try:
                            # remove ponds that are actually tall buildings in a low spot
                            ponds_hag = ZonalStatistics(ponds2smoothRg, 'Value', hag, 'MEAN')
                            hag0 = Con(IsNull(hag), 0, hag)
                            ponds_hag0 = ZonalStatistics(ponds2smoothRg, 'Value', hag0, 'MEAN')
                            # threshold at 200 cm
                            hag0_threshold = 200
                            ponds2keep = Con(ponds_hag0 < hag0_threshold, ponds2smoothRg)

                            # int2keep = Con(hag < hag0_threshold, low_high_int)
                            # int1rPlus1 = intRaster + 1
                        except:
                            ponds2keep = ponds2smoothRg

                        if ponds2keep.maximum is not None:

                            ponds2FlattenPoly = arcpy.RasterToPolygon_conversion(ponds2keep, opj(gdb, 'ponds_2_flatten'))

                            # get a good estimate of the actual pond elevation 
                            flat_fs5_range_el = FocalStatistics(flatterRiverDEM, NbrRectangle(5,5), 'RANGE')
                            stability_threshold = 100
                            stable_flatterRiver = Con(flat_fs5_range_el < stability_threshold, flatterRiverDEM)
                            minPondEl = ZonalStatistics(ponds2keep, 'VALUE', stable_flatterRiver, 'MINIMUM')#pitFilledDEM, 'MINIMUM')
                            majPondEl = ZonalStatistics(ponds2keep, 'VALUE', stable_flatterRiver, 'MAJORITY')#pitFilledDEM, 'MAJORITY')
                            # let's go with minimum plus 50% of the way to the majority
                            bestPondEl = Int(minPondEl + 0.5*(majPondEl-minPondEl))

                            # we don't want the ponds to be completely flat or Flow Direction goes nuts
                            fstBestPondEl = FocalStatistics(bestPondEl, circle3Nbr, 'MAXIMUM')
                            fstBestPondElPlus1 = fstBestPondEl + 1
                            flatPondDemWithPits = Con(IsNull(bestPondEl), flatterRiverDEM, bestPondEl)

                            pondElLE_FstBest = flatPondDemWithPits <= fstBestPondEl
                            pondElAndPlus1 = Con(pondElLE_FstBest == 1, fstBestPondElPlus1, flatPondDemWithPits)
                            flatPondDEM = Con(IsNull(pondElAndPlus1), flatterRiverDEM, pondElAndPlus1)
                            flatPondDemFd = FlowDirection(flatPondDEM)
                            flatPondDemBasins = Basin(flatPondDemFd)
                            fp_basins_rtp = arcpy.conversion.RasterToPolygon(flatPondDemBasins, opj(sgdb, 'flat_ponds_basins'), 'NO_SIMPLIFY')

                            ndFixedList.append(ponds2keep)#Flatten)
                            log.debug(str(ndFixedList))
                            log.debug('nd here! at ' + time.asctime())

                        else:
                            flatPondDEM = flatterRiverDEM
                    else:
                        flatPondDEM = flatterRiverDEM
                else:
                    flatPondDEM = flatterRiverDEM


        ####----------------------------------------------------------------------------
            # define the smaller no data areas to process in next step
                if len(ndFixedList) > 0:
                    log.debug('tracking the no data areas')
                    cs = CellStatistics(ndFixedList, 'MAXIMUM')
                    fixedNdPartial = Con(cs, ndRegions, '', 'VALUE >=1')
                    fixedNdWhole = df.fullZoneByZs(fixedNdPartial, ndRegions)
                    moreNdRegions2Fix = Con(IsNull(fixedNdWhole), ndRegions)
                else:
                    moreNdRegions2Fix = ndRegions#goodCrvNdRegions
                moreNdRegions2Fix.save(mediumNoDataAreas)
    ##
    ##    ####----------------------------------------------------------------------------
    ##            ## Fix bridge elevations

                bridge_fix_list = []

                # Clip = opj(sgdb, 'buf_' + huc12 + '_1km')
                inm = 'in_memory'
                hucRoadsIncService = arcpy.Clip_analysis(roadsFc, Clip, opj(inm, "roads_all"))
                if 'fclass' in df.getfields(hucRoadsIncService):
                    hucRoadsWhole = arcpy.Select_analysis(hucRoadsIncService, opj(inm, 'roads'), 'fclass <> \'service\' AND fclass <> \'path\'AND fclass <> \'cycleway\'AND fclass <> \'footway\'')
                else:
                    hucRoadsWhole = arcpy.Select_analysis(hucRoadsIncService, opj(inm, 'roads'), 'type <> \'service\' AND type <> \'path\'')

                hucRoads = arcpy.FeatureToLine_management(hucRoadsWhole, opj(inm, 'roads_broken'))
                hucRoadsBuffer = arcpy.Buffer_analysis(hucRoads, opj(inm, 'roads_buffer'), str(proc_size) + ' METERS')

                hucRoadsRaster = arcpy.PolygonToRaster_conversion(hucRoadsBuffer, 'oneway', opj(sgdb, 'roads_rast'), cellsize = str(proc_size))

                partialNdRegions = Con(hucRoadsRaster, moreNdRegions2Fix)

                # expand rDepth by a cell to make sure it hits road
                rDepth_fs = FocalStatistics(rDepth)
                partial_nd_in_rdepth = Con(rDepth_fs, partialNdRegions)

                if df.testForZero(partial_nd_in_rdepth):#False

                    full_nd_2_test = df.fullZoneByZs(partial_nd_in_rdepth, moreNdRegions2Fix)

                    # make sure the fill region is deep enough to be worth correcting
                    rDepth_rg = RegionGroup(rDepthPosAs1, 'EIGHT')
                    rDepth_rg_min = ZonalStatistics(rDepth_rg, 'VALUE', pitFilledDEM, 'MINIMUM')
                    rDepth_rg_max = ZonalStatistics(rDepth_rg, 'VALUE', pitFilledDEM, 'MAXIMUM')
                    # use 1.5x IA 2010 RMSE as threshold
                    rDepth_test = (rDepth_rg_max - rDepth_rg_min) > 27#75
                    rDepth_filtered = Con(rDepth_test, rDepth_rg_min)
                    ####            rDepth_min_bridges = Con(bridge_voids, rDepth_rg_min)

                    # some no data areas on larger streams/rivers may be very big, trim to road buffer
                    full_nd_in_road_buffer = Con(hucRoadsRaster, full_nd_2_test)

                    full_regroup_all = RegionGroup(full_nd_in_road_buffer, 'EIGHT')
                    # remove regions with area less than 5 cells at 2m resolution?
                    size_threshold = 5 * pow(2,2) / pow(proc_size, 2)
                    full_regroup_big = Con(full_regroup_all, full_regroup_all, '', "COUNT > " + str(size_threshold))
    ##                rtp_possible_bridges = arcpy.conversion.RasterToPolygon(full_regroup_big, opj(sgdb, 'maybe_bridges'), 'NO_SIMPLIFY')

                    hag_gt50 = Con(hag > 50, hag)
                    hag_gt50_range = FocalStatistics(hag_gt50, '', 'RANGE')
    ##                hag_gt50_std = FocalStatistics(hag_gt50, '', 'STD')
                    hag0_gt50_range = Con(IsNull(hag_gt50_range), 0, hag_gt50_range)
                    hag_median = ZonalStatistics(full_regroup_big, 'VALUE', hag0_gt50_range, 'MEDIAN')
                    # if range > 150, likely vegetation due to excessive range from HUC102000010607 and 102500091002
                    range_threshold = 150
                    full_regroup = Con(hag_median < range_threshold, full_regroup_big)
                    
                    full_nd_table = arcpy.CopyRows_management(full_regroup, opj(sgdb, 'full_nd_tbl'))
                    full_nd_values = [row[0] for row in arcpy.da.SearchCursor(full_nd_table, ['VALUE'])]

                    for j in range(1, 3):
        ##            j = 1
                        if j == 1:
                            nd_2_test = full_regroup
                        else:
                            failed_first_time = Con(IsNull(bridge_test_passed), full_regroup)
                            nd_2_test = failed_first_time

                        if len(full_nd_values) > 4000:
                            full_nd_values1 = full_nd_values[:int(len(full_nd_values)/2)]
                            full_nd_values2 = full_nd_values[int(len(full_nd_values)/2):]
                            expand_full_nd1 = Expand(nd_2_test, j, full_nd_values1)
                            expand_full_nd2 = Expand(nd_2_test, j, full_nd_values2)
                            expand_full_nd = CellStatistics([expand_full_nd1, expand_full_nd2], 'MAXIMUM')
                        else:
                            expand_full_nd = Expand(nd_2_test, j, full_nd_values)

                        just_expand_nd = Con(IsNull(nd_2_test), expand_full_nd)

                        just_expand_min = ZonalStatistics(just_expand_nd, 'VALUE', pitFilledDEM, "MINIMUM")
                        just_expand_max = ZonalStatistics(just_expand_nd, 'VALUE', pitFilledDEM, "MAXIMUM")

                        if j == 1:
                            low_thresh = 0.2
                            high_thresh = 0.8

                        else:
                            low_thresh = 0.25
                            high_thresh = 0.75

                        min_plus_25 = just_expand_min + low_thresh * (just_expand_max - just_expand_min)
                        min_plus_75 = just_expand_min + high_thresh * (just_expand_max - just_expand_min)

                        le_twenty = Con(pitFilledDEM <= min_plus_25, 1)
                        ge_eighty = Con(pitFilledDEM >= min_plus_75, 1)

                        le_twenty_rg = RegionGroup(le_twenty, 'EIGHT')
                        ge_eighty_rg = RegionGroup(ge_eighty, 'EIGHT')

                        expand_le_rg_var = ZonalStatistics(just_expand_nd, 'VALUE', le_twenty_rg, 'VARIETY')
                        expand_ge_rg_var = ZonalStatistics(just_expand_nd, 'VALUE', ge_eighty_rg, 'VARIETY')

                        bridge_voids = (expand_le_rg_var >= 2) & (expand_ge_rg_var >= 2)

                        bridge_test_TF = ZonalStatistics(expand_full_nd, 'Value', bridge_voids, 'MAXIMUM')

                        bridge_test_passed = Con(bridge_test_TF > 0, expand_full_nd)

                        # make sure the elevations of the fill region and no data expansion are similar
                        expand_full_nd_rdepth_min = ZonalStatistics(expand_full_nd, 'Value', rDepth_filtered, 'MINIMUM')
                        
                        expand_full_nd_expand_min = ZonalStatistics(expand_full_nd, 'Value', just_expand_min, 'MINIMUM')
                        rdepth_expand_min_dif = expand_full_nd_expand_min - expand_full_nd_rdepth_min
        ####                # if not, filter out those with a large difference in elevation
        ####                rDepth_filtered2 = Con(rdepth_expand_min_dif < 50, rDepth_filtered)

                        bridge_voids_min_el = ZonalStatistics(expand_full_nd, 'Value', rDepth_filtered, 'MINIMUM')

                        if j == 1:
                            # recode the bridge and lowest elevations to fill region minimum
                            void_and_low = CellStatistics([nd_2_test, le_twenty], 'MAXIMUM')
                        else:
                            pass
                        # get a path between lowest elevations and highest elevations for recoding elevation to edge/outside of void
                        bridge_voids_le_rg = Con(bridge_test_passed, le_twenty_rg)
                        bridge_voids_le_rg_min = ZonalStatistics(bridge_test_passed, 'VALUE', bridge_voids_le_rg, 'MINIMUM')
                        bridge_voids_le_rg_max = ZonalStatistics(bridge_test_passed, 'VALUE', bridge_voids_le_rg, 'MAXIMUM')
                        bridge_voids_min = Con(bridge_voids_le_rg_min == bridge_voids_le_rg, bridge_voids_le_rg)
                        bridge_voids_max = Con(bridge_voids_le_rg_max == bridge_voids_le_rg, bridge_voids_le_rg)
                        rel_el = pitFilledDEM - bridge_voids_min_el
                        cd = CostDistance(bridge_voids_min, rel_el, out_backlink_raster = opj(sgdb, 'back_lnk_bridge'))
                        bridge_cp = CostPath(bridge_voids_max, cd, opj(sgdb, 'back_lnk_bridge'))

                        min_location_in_fr = rDepth_rg_min == pitFilledDEM
                        fr_min_eq_bridge_voids_min_el = rDepth_rg_min == bridge_voids_min_el
                        fr_min_eq_as_fr = Con(fr_min_eq_bridge_voids_min_el, rDepth_rg)
                        full_fr_min_eq = df.fullZoneByZs(fr_min_eq_as_fr, rDepth_rg)
                        fr_min_location_for_bridge = Con(full_fr_min_eq, min_location_in_fr)

                        rel_el2 = pitFilledDEM - rDepth_rg_min
                        rel_el_cs_min = CellStatistics([rel_el, rel_el2], 'MINIMUM')
                        # can't have 0 in cost distance in Pro
                        rel_el_cs_min_plus_1 = rel_el_cs_min + 1
                        fr_min_location_for_bridge_T = Con(fr_min_location_for_bridge, rDepth_rg)
                        cd_to_fr_min = CostDistance(fr_min_location_for_bridge_T, rel_el_cs_min_plus_1, out_backlink_raster = opj(sgdb, 'back_link_fr_min'))
                        void_and_low_pre = CellStatistics([nd_2_test, le_twenty], 'MAXIMUM')
                        cp_to_fr_min = CostPath(void_and_low_pre, cd_to_fr_min, opj(sgdb, 'back_link_fr_min'))
        ##                cp_to_fr_min = CostPath(bridge_test_passed, cd_to_fr_min, opj(sgdb, 'back_link_fr_min'))
        ##                cp_to_fr_min2 = CostPath(bridge_cp, cd_to_fr_min, opj(sgdb, 'back_link_fr_min'))
                        # we don't want to change the minimum elevation
                        cp_to_fr_min_shorter = Con(IsNull(fr_min_location_for_bridge_T), cp_to_fr_min)
                        # this is more than just the true bridges, filter to true bridges using new_bridge_el
        ##                void_and_low = CellStatistics([nd_2_test, le_twenty, cp], 'MAXIMUM')
                        void_and_low= CellStatistics([void_and_low_pre, cp_to_fr_min_shorter], 'MAXIMUM')


                        # increase min el by 1 so it will flow to deepest cell in fr
                        new_bridge_el_plus1 = Con(bridge_test_TF == 1, bridge_voids_min_el + 1)
                        fr_min_plus1 = rDepth_rg_min + 1
                        fr_bridge_el_plus1 = CellStatistics([new_bridge_el_plus1, fr_min_plus1], 'MINIMUM')

                        bridge_fix_el = Con(void_and_low, fr_bridge_el_plus1)
                        bridge_fix_el.save(opj(sgdb, 'b_fx_el' + str(j)))
                        bridge_fix_list.append(bridge_fix_el)

                    if len(bridge_fix_list) > 0:
                        bridge_fix_els = CellStatistics(bridge_fix_list, 'MINIMUM')

                        # replace bridge and nearbly low elevations in DEM
                        bridge_fixed_dem = Con(IsNull(bridge_fix_els), flatPondDEM, bridge_fix_els)

                        ndAllFixedDEM = bridge_fixed_dem
                    else:
                        ndAllFixedDEM = flatPondDEM

                else:
                    ndAllFixedDEM = flatPondDEM
                    
                ndAllFixedDEM_fd = FlowDirection(ndAllFixedDEM)
                ndAllFixedDEM_fa = FlowAccumulation(ndAllFixedDEM_fd)
                fa_gt_20ac = ExtractByAttributes(ndAllFixedDEM_fa, 'VALUE > ' + str(20*4047 / pow(proc_size, 2)))
                fa_gt_20ac_int = Int(fa_gt_20ac)
                stf = StreamToFeature(fa_gt_20ac_int, ndAllFixedDEM_fd, opj(sgdb, 'nd_streams2'))

                ## check to see if DEM is 'right' size compared to pitfill DEM area, before save
                vfArea = ZonalGeometry(Con(IsNull(flatPondDEM) == 0, 1), 'VALUE', 'AREA')
                pfArea = ZonalGeometry(Con(IsNull(fillTif) == 0, 1), 'VALUE', 'AREA')
                areaRatio = vfArea.maximum / pfArea.maximum
                try:
                    assert(areaRatio >= 0.97),'VoidFixed DEM failed area data check'
                    ndAllFixedDEM.save(voidFixTif)
                    gdbCopy = arcpy.CopyRaster_management(ndAllFixedDEM, opj(sgdb, os.path.splitext(ndAllFixedDEM.name)[0]))
                except AssertionError:
                    log.error('VoidFixed DEM failed area data check', exc_info = True)



            else:
                log.warning('WARNING - Insufficient Return DEM')

        else:
            log.warning('WARNING - No Return DEM')

    ##except AssertionError:

    except:
        winsound.Beep(500, 1000)
        winsound.Beep(1000, 1000)
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

if __name__ == "__main__":
##if True:

    if len(sys.argv) == 1:
        arcpy.AddMessage("Whoo, hoo! Running from Python Window!")
        cleanup = False

        parameters = ["C:/Program Files/ArcGIS/Pro/bin/Python/envs/arcgispro-py3/pythonw.exe",
	"C:/DEP/Scripts/basics/cmd_flattener_DEM.pyt",
	"C:/DEP/LiDAR_Current/elev_FLib_mean18/07080105/ef3m070801050901.tif",
	"C:/DEP/LiDAR_Current/count_Lib/07080105/cbe3m070801050901.tif",
    "C:/DEP/LiDAR_Current/count_Lib/07080105/cfr3m070801050901.tif",
	"C:/DEP/LiDAR_Current/surf_el_Lib/07080105/frmax3m070801050901.tif",
	"C:/DEP/LiDAR_Current/int_Lib/07080105/fr_int_max3m070801050901.tif",
	"C:/DEP/Man_Data_ACPF/dep_ACPF2022/07080105/idepACPF070801050901.gdb/buf_070801050901",
	"C:/DEP/Basedata_Summaries/Basedata_26915.gdb/roads_merge",
	"C:/DEP/Basedata_Summaries/Basedata_26915.gdb/waterways",
	"C:/DEP/Basedata_Summaries/Basedata_26915.gdb/water",
	"C:/DEP/LiDAR_Current/bl_Lib/07080105/breaks_07080105.gdb/break_polys_070801050901",
	"C:/DEP_Proc/DEMProc/Void_dem2013_3m_070801050901",
	"C:/DEP/LiDAR_Current/elev_VLib_mean18/07080105/ev3m070801050901.tif",
	"C:/DEP/LiDAR_Current/voids_Lib_mean18/07080105/bigvds3m070801050901.tif",
	"C:/DEP/LiDAR_Current/voids_Lib_mean18/07080105/medvds3m070801050901.tif"]

        for i in parameters[2:]:
            sys.argv.append(i)
    else:
        arcpy.AddMessage("Whoo, hoo! Command-line enabled!")
        # DO NOT clean up the folder after done processing - matcher needs this data
        cleanup = False

    messages = msgStub()

    # input_dem, output_dem, plib_metadata, depressions_fc, depth_threshold, area_threshold, procDir = [i for i in sys.argv[1:]]
    fillTif, cntTif, cnt1rTif, surfaceElevFile, int1rMaxFile, buf_bnd, roadsFc, input_waterway, input_water, breakpolys, voidProc, voidFixTif, bigNoDataAreas, mediumNoDataAreas = [i for i in sys.argv[1:]]

    doFlattener(fillTif, cntTif, cnt1rTif, surfaceElevFile, int1rMaxFile, buf_bnd, roadsFc, input_waterway, input_water, breakpolys, voidProc, voidFixTif, bigNoDataAreas, mediumNoDataAreas, cleanup, messages)
    arcpy.AddMessage("Back from doing!")
