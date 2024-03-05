# coding: utf-8
## 2022.06.01 - Identified issues with expLowAsp1, expLowAsp when this is not being created locally, these lines are commented out and replaced with revised code - bkgelder

import arcpy
from arcpy.sa import *
import os
import time
import traceback
import sys
import platform
import winsound
sys.path.append(os.getcwd())
# sys.path.append('O:\\DEP\\Scripts\\basics')
# fix for embedded code (using 'load code' in ArcGIS Pro) fails when testing len(sys.argv)
if not hasattr(sys, 'argv'):
    sys.argv = ['']
import dem_functions as df
from os.path import join as opj


def step_zonal_fill(step_ws, step_dem, base_dem, counter):
    # zonal fill the initial watersheds
    step_zf = ZonalFill(step_ws, step_dem)
    # calculate where zonal fill is an elevation increase (needs to change)
    # alway refer back to initDEM so thickness evaluation is always against initial DEM
    step_zf_dif = step_zf - base_dem#
##    step_zf_dif.save('zfd_' + str(counter))
    step_zf_dem = Con(step_zf_dif >= 0, step_zf, step_dem)
##    step_zf_dem.save('zfdem_' + str(counter))

    # calculate how thick the positive fill is
    stepPosFill = Con(step_zf_dif > 0, step_ws)
##    stepPosFill.save('pos_fill' + str(counter))
    stepZThick = ZonalGeometry(stepPosFill, 'Value', 'THICKNESS')
##    stepZThick.save('zthick_' + str(counter))
    #now get this data all over the basin for generating new DEM
    step_ws_zf_thick = ZonalStatistics(step_ws, 'Value', stepZThick, 'MAXIMUM')

    return step_zf_dem, step_ws_zf_thick, step_zf, step_zf_dif

def create_thin_update_dem(step_ws_zf_thickness, step_zf_dem, step_dem, counter, thickThreshold):
    # create new DEM and figure out where water will flow
    step_tu_DEM = Con(step_ws_zf_thickness, step_zf_dem, step_dem, 'VALUE < ' + str(thickThreshold))
##    step_tu_DEM.save('tudemz_' + str(counter))
    stepTuFd = FlowDirection(step_tu_DEM)
    step_tu_basins = Basin(stepTuFd)

    return step_tu_DEM, step_tu_basins

def create_zf_update_dem(step_ws_zf_tested, step_zf_dem, step_dem, counter):
    # create new DEM and figure out where water will flow
    step_tu_DEM = Con(step_ws_zf_tested, step_zf_dem, step_dem)#, 'VALUE = 1')
##    step_tu_DEM.save('tudemz_' + str(counter))
    stepTuFd = FlowDirection(step_tu_DEM)
    step_tu_basins = Basin(stepTuFd)

    return step_tu_DEM, step_tu_basins


def create_dem_filled_to_holes(step_tu_basins, step_deepest_cell, base_dem, counter):
    '''recover initial DEM around initial sinks of each new watershed
       need to keep each new sink (deep zonal fill area) separate, can't just route to deepest sink in tuws
       and must remember that deepest zonal filled area does not always have the deepest sink
    first find deepest cells in each tuws
    step_tu_basins = potential aggregation basins due to zonal filling
    step_deepest_cell = deepest cell in potential aggregation basins, used as hole
    base_dem = DEM to use when filling to holes'''
    stepDeepestCellInTuWs = ZonalStatistics(step_tu_basins, 'VALUE', step_deepest_cell, 'MINIMUM')
    stepDeepestCellAs1 = Con(stepDeepestCellInTuWs == step_deepest_cell, 1)
##    stepDeepestCellAs1.save('stpdpst1_' + str(counter))
    next_step_deepest_cells_el = Con(stepDeepestCellAs1 == 1, step_deepest_cell)
    stepHolesAtDeepestDEM = Con(IsNull(stepDeepestCellAs1), base_dem)

    counter += 1
    log.warning('counter: ' + str(counter))
##    stepHolesAtDeepestDEM.save('stp_dem' + str(counter))
    stepfill2Deepest = Fill(stepHolesAtDeepestDEM)
    next_step_dem_no_holes = Con(IsNull(stepfill2Deepest), base_dem, stepfill2Deepest)#fillInitDEM, stepfill2Deepest)

    return next_step_dem_no_holes, next_step_deepest_cells_el, stepHolesAtDeepestDEM, counter
    



def fixByInversionByPath(ndPlus, fenceEl, invertTargetDEM, spot4Hole, ndRcls, bestestDEM, unique_path, log):

    ## To correct an inverted DEM, remember water must flow out the upstream ends (before inversion), not downstream!
    ## ndPlus is 2 at regions2Fix, else 1 in area to process
    fencedRegion2Fix = Pick(ndPlus, [fenceEl, invertTargetDEM])#originalDEM])
    log.warning('fencedRegion2Fix at time: ' + time.asctime())
    holeAtMin = Con(IsNull(spot4Hole) == 1, fencedRegion2Fix, '')
    
    fillNdBarrier = Fill(holeAtMin)
    log.warning('fillNdBarrier at time: ' + time.asctime())
##                                    faFence = FlowAccumulation(FlowDirection(fillNdBarrier))

    ## Calculate the difference between the original inverted and filled inverted DEMs
    ## This should be the change to apply to the initial DEM make things flow
    fillNdDif = fillNdBarrier - fencedRegion2Fix

    # 2021.06.10 - identified issue with above approach - if hole is at a local minimum,
    # filling process will not lower to the hole elevation
    # corrective action - use elevation of hole to constrain correction elevation,
    # thus no elevations greater than hole elevation will be allowed to pass
    # requires grouping of all cost paths to unique, consistent zones
    correctionToApply_unfiltered = Con(ndRcls, bestestDEM - fillNdDif)#originalDEM - fillNdDif)
    spot4HoleElevation = Con(spot4Hole, bestestDEM)
    log.warning('spot4Hole at time: ' + time.asctime())
    pathHoleElevation = ZonalStatistics(unique_path, 'VALUE', spot4HoleElevation, 'MINIMUM')
    correctionToApply = Con(correctionToApply_unfiltered > pathHoleElevation, pathHoleElevation, correctionToApply_unfiltered)
    log.warning('correctionToApply at time: ' + time.asctime())
    correctedDEM = Int(Con(IsNull(correctionToApply), bestestDEM, correctionToApply))
    
    return correctedDEM

def fixByInversionByStartingPath(ndPlus, fenceEl, invertTargetDEM, spot4Hole, ndRcls, bestestDEM, starting_path_min, log):

    ## To correct an inverted DEM, remember water must flow out the upstream ends (before inversion), not downstream!
    ## ndPlus is 2 at regions2Fix, else 1 in area to process
    fencedRegion2Fix = Pick(ndPlus, [fenceEl, invertTargetDEM])
    log.warning('fencedRegion2Fix at time: ' + time.asctime())
    holeAtMin = Con(IsNull(spot4Hole) == 1, fencedRegion2Fix, '')
    
    fillNdBarrier = Fill(holeAtMin)
    log.warning('fillNdBarrier at time: ' + time.asctime())

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
    log.warning('correctionToApply at time: ' + time.asctime())
    isNullPathStartingCorrection = IsNull(pathStartingNeedsCorrectionTrue)
    filtered_corrections_everywhere = Con(isNullPathStartingCorrection == 0, starting_path_min_fs, correctionToApply_unfiltered)
    isnull_filtered_corrections = IsNull(filtered_corrections_everywhere)
    correctedDEM = Con(isnull_filtered_corrections, bestestDEM, filtered_corrections_everywhere)

    return correctedDEM

def startByFlowLengthDown(cost_path, flow_direction, basin):#, proc_size):
    '''Determines the initiation point of a cost path by its flow length
    in the downstream direction.'''
    fd_on_cp = Con(cost_path, flow_direction)
    flowLenAlongCp = FlowLength(fd_on_cp)
    fsMaxFlDown = FocalStatistics(flowLenAlongCp, "RECTANGLE 3 3 CELL", 'MAXIMUM')
    fsMaxTF = flowLenAlongCp == fsMaxFlDown
    flowLenDnEqMax = Con(fsMaxTF, basin, where_clause = 'VALUE = 1')

    return flowLenDnEqMax

def startByFlowLengthUp(cost_path, flow_direction, basin):#, proc_size):
    '''Determines the initiation point of a cost path by its flow length
    in the upstream direction.'''
    fd_on_cp = Con(cost_path, flow_direction)
    flowLenUpAlongCp = FlowLength(fd_on_cp, 'UPSTREAM')
    fsMaxFlUp = FocalStatistics(flowLenUpAlongCp, "RECTANGLE 3 3 CELL", 'MAXIMUM')
    flowLenUpEq0 = Con(flowLenUpAlongCp, basin, where_clause = 'VALUE = 0')
    flowLenUpEq0AtAllStart = Con(fsMaxFlUp, flowLenUpEq0, where_clause = 'VALUE <= ' + str(1.5*flow_direction.meanCellHeight))

    return flowLenUpEq0AtAllStart

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

    return refiltered_starting_basins_min_replace_nulls

def create_no_pass_through_basins(basin, starts, cp, dem):
    '''
    Function filters starting points for cost paths to find those that do not
    start in 'pass through' basins that are lower on the stream network. These
    may have initation points but are not a real local minimum for inversion
    basin - The basin zones to test
    starts - The initiating points for cost paths
    cp - The cost path from the starts
    '''
    # find which cost path starts do not have basins passing through them
    # (but are not multi-cell deepest)
    basin_cp_starts = ZonalStatistics(basin, 'Value', starts, 'VARIETY')
    basin_cp_min = ZonalStatistics(basin, 'Value', cp, 'MINIMUM')
    non_pass_through_basins_pre = Con(basin_cp_min > 2, basin_cp_starts)
    starts_no_pass_through = Con(non_pass_through_basins_pre, starts)

    # find those cost paths that start with multiple cells at deepest elevation
    # the cost path tool identifies each traversed cell at this deepest elevation as another
    # cost path start. Thus they are excluded in the previous selection because cost path is > 2 very quickly
    basin_cp_max = ZonalStatistics(basin, 'Value', cp, 'MAXIMUM')
    non_pass_through_basins_pre2 = Con(basin_cp_max > 2, basin_cp_starts)
####    starts_no_pass_through2 = Con(non_pass_through_basins_pre2, starts)

    basin_cp_starts_full = df.fullZoneByZs(basin_cp_starts, basin)
    basin_cp_starts_full_min = ZonalStatistics(basin_cp_starts_full, 'VALUE', dem, 'MINIMUM')
####    eq_full_min = Con(dem == basin_cp_starts_full_min, dem)
    eq_full_min_basin = Con(dem == basin_cp_starts_full_min, basin)
####    double_deep_eq_min = Con(eq_full_min_basin, eq_full_min_basin, '', 'COUNT > 1')
    ####full_double_deep = df.fullZoneByZs(double_deep_eq_min, basin)
    starts_no_pass_through2 = Con(eq_full_min_basin, eq_full_min_basin, '', 'COUNT > 1')

    starts_no_pass_through_all = CellStatistics([starts_no_pass_through, starts_no_pass_through2], 'MAXIMUM')
    starts_no_pass_through_full = df.fullZoneByZs(starts_no_pass_through_all, basin)

    return starts_no_pass_through_all, starts_no_pass_through_full

def get_cost_path_favor_deeper(pre_zf_deepest_cells_el, fd_for_cp, pre_zf_dem, step_dem_cp, stepWs, stepDEMwithHoles, step_deepest_cells_el, fixed_dif, counter):
    '''This function is designed to find the 'best' cost path from previous
    deepest cells to new deepest cells (which are a subset of the previous
    ones). It does this by creating potential pathways between these points
    and then favoring those that traverse the deepest portions of the pathway.
    pre_zf_deepest_cells_el = used to trace previous end points to new end points after step
    fd_for_cp = flow direction used for back link in cost path tracing
    step_dem_cp = DEM used for cost path tracing
    pre_zf_dem = used to calculate current fill depth for extended area traversal as well as weight new cost path
    stepWs = used to get minimum value for each new step watershed
    stepDEMwithHoles =
    step_deepest_cells_el =
    fixed_dif = to find areas where DEM has already been fixed
    counter = '''
    
    # first try - get the downstream path from the previous deepest cells to the new deepest cells
    current_cp = CostPath(pre_zf_deepest_cells_el, step_dem_cp, fd_for_cp)

    # then see where we filled since previous DEM and add 0 depth connected cells of same elevation
    # REMEMBER - flow direction will allow you to traverse a cell with 0 depth of flow
    current_fill_depth = step_dem_cp - pre_zf_dem#dem_no_holes_after_fill_after_crv_test - pre_zf_dem#fillInitDEM
    ##fill_depth_pos = Con(current_fill_depth > 0, current_fill_depth)
    pos_dif_step_el = Con(current_fill_depth > 0, step_dem_cp)
    # speed up region grouping of elevation by only looking out 3 cells in each direction
    fs7_pos_dif_max = FocalStatistics(pos_dif_step_el, NbrRectangle(7,7), 'MAXIMUM')
    reduced_step_DEM = Con(step_dem_cp <= fs7_pos_dif_max, step_dem_cp)
    reduced_step_DEM_rg = RegionGroup(reduced_step_DEM, 'EIGHT')
    # valid regions must have had some non 0 fill depth
    reduced_step_DEM_rg_in_fill = ZonalStatistics(reduced_step_DEM_rg, 'Value', pos_dif_step_el, 'MAXIMUM')

    # now add in areas where we initially fixed the DEM
    ss_fix_dif_neg = Con(fixed_dif < 0, step_dem_cp)

    # finally add in the cells at the minimum elevation of the DEM with holes
    # for some reason the current_cp above doesn't always flow to the hole, it ends on these flats
    step_dem_cp_holes_min = ZonalStatistics(stepWs, 'VALUE', stepDEMwithHoles, 'MINIMUM')
    eq_step_dem_cp_holes_min = Con(step_dem_cp == step_dem_cp_holes_min, step_dem_cp)

    # now join all four of these into a potential path that we evaluate for best path
    cp_and_pos_fill = CellStatistics([current_cp, reduced_step_DEM_rg_in_fill, ss_fix_dif_neg, eq_step_dem_cp_holes_min], 'MAXIMUM')

    ####cp_and_pos_fill = CellStatistics([current_cp, fill_depth_pos], 'MAXIMUM')
    cp_and_pos_fill_elev = Con(cp_and_pos_fill, pre_zf_dem)
    fs_min_el_combo_area = FocalStatistics(cp_and_pos_fill_elev, statistics_type = 'MINIMUM')
    difFromFSMinimumCp = Con(cp_and_pos_fill_elev, pre_zf_dem - fs_min_el_combo_area)
    difPlus1 = difFromFSMinimumCp + 1
    cost_dist_from_deepest = CostDistance(step_deepest_cells_el, difPlus1, out_backlink_raster = opj(sgdb, 'back_lnk' + str(counter)))
    cp_favor_deeper = CostPath(pre_zf_deepest_cells_el, cost_dist_from_deepest, opj(sgdb, 'back_lnk' + str(counter)))

    return cp_favor_deeper, cost_dist_from_deepest, difPlus1

def create_cp_dem_from_stops(all_cp_stops, pre_zf_dem):
    '''Create a dem that drains to a new set of cost path stops by using the
    stops to add new holes to a dem then fill to those holes. Also calculate flow
    direction.
    '''
    dem_holes_all_cp_stops = Con(IsNull(all_cp_stops), pre_zf_dem)
    fill_2_stops = Fill(dem_holes_all_cp_stops)
    fill_2_stops_fd = FlowDirection(fill_2_stops)
    cp_dem = Con(IsNull(fill_2_stops), pre_zf_dem, fill_2_stops)
    cp_fd = FlowDirection(cp_dem)#dem_no_holes_2_all_cp_stops)

    return cp_dem, cp_fd, fill_2_stops, fill_2_stops_fd

def find_last_deepest_go_cell(go_raster, cp, step_fill_fd, pre_zf_deepest):
    '''Find the last deepest cell in a cost path (aka cost path stopping point)
    that is in the cost path so we can stop the cp there.
    '''
    # find cost path areas outside allowed, find local max distance from #1
    # add new stops into this area here, like roads or levees
    cp_go = Con(go_raster, cp)#_favor_deeper)#expLowAsp
    cp_go1 = Con(cp_go, 1)
    cp_go1_rg = RegionGroup(cp_go1, 'EIGHT')

    # get max fa cell in upper cost path segments traversing in allowed regions
    stepfill2Fa = FlowAccumulation(step_fill_fd)#stepfill2Fd)
    stepfill2Fa_fs = FocalStatistics(stepfill2Fa, NbrRectangle(3,3), 'MAXIMUM')
    pre_zf_deepest_max_fa_fs = Con(pre_zf_deepest, stepfill2Fa_fs)

    # find the pre_zf_deepest cell that needs to be retained in the upper good cost path segments
    cp_go1_rg_max_fa = ZonalStatistics(cp_go1_rg, 'Value', pre_zf_deepest_max_fa_fs, 'MAXIMUM')
    eq_go1_max_fa = Con(pre_zf_deepest_max_fa_fs == cp_go1_rg_max_fa, 1)

    return eq_go1_max_fa



def get_cost_path_with_new_breaks(expLowAspTF, roadsRaster, cp_favor_deeper, cost_dist_from_deepest, dif_plus1, cp_counter, stepfill2Fd, pre_zf_deepest, step_deepest, pre_zf_dem, log):
    '''Revises a cost path to stop at various restrictions including profile
    curvature and roads. This requires adding additional stops to the cost
    path, thus slicing the original path into multiple segments where aggregation
    and DEM inversion is allowable.
    '''

    nogo_counter = 10*(cp_counter // 10)

    roadsRasterTF = IsNull(roadsRaster)

    ####    go =
    cp_nogo = Con(expLowAspTF == 0, cp_favor_deeper)
    cp_nogo_big = Con(cp_nogo, cp_nogo, '', 'COUNT > 2')
    arcpy.BuildRasterAttributeTable_management(cp_nogo_big)
    cp_nogo_big_count = int(arcpy.management.GetCount(cp_nogo_big).getOutput(0))
    ##    while cp_nogo_big_count > 2:

##    eq_go1_max_fa = find_last_deepest_go_cell(expLowAsp, cp_favor_deeper, stepfill2Fd, pre_zf_deepest)
    eq_go1_max_fa = find_last_deepest_go_cell(expLowAspTF, cp_favor_deeper, stepfill2Fd, pre_zf_deepest)
    eq_go1_max_fa.save(opj(sgdb, 'eq_go1_max_fa' + str(nogo_counter)))

    # aggregate stops into one raster
    all_cp_stops = CellStatistics([step_deepest, eq_go1_max_fa], 'MINIMUM')

    cp_dem_2_stops, cp_dem_2_stops_fd, fill_2_stops, fill_2_stops_fd = create_cp_dem_from_stops(all_cp_stops, pre_zf_dem)

    # calculate new cost path based on dem with additional stops and all previous
    # stops that were removed as starting points
    ##pre_zf_deepest_that_matter = Con(IsNull(all_cp_stops), pre_zf_deepest)
    ##cp_to_new_holes = CostPath(pre_zf_deepest_that_matter, fill_2_stops, fill_2_stops_fd)#cp_dem_2_stops, cp_dem_2_stops_fd)
    ####

    ws_for_new_stops = Watershed(fill_2_stops_fd, eq_go1_max_fa)
    pre_zf_deepest_in_nogo = Con(expLowAspTF == 0, pre_zf_deepest)
    pre_zf_deepest_in_nogo_in_new_ws = Con(ws_for_new_stops, pre_zf_deepest_in_nogo)
    pre_zf_deepest_out_new_ws = Con(IsNull(ws_for_new_stops), pre_zf_deepest)
    pre_zf_deepest_remaining = CellStatistics([pre_zf_deepest_in_nogo_in_new_ws, pre_zf_deepest_out_new_ws], 'MAXIMUM')
    cp_from_deepest_remaining = CostPath(pre_zf_deepest_remaining, cost_dist_from_deepest, opj(sgdb, 'back_lnk' + str(cp_counter)))



    # now try a different way to stop the bad cost paths from being created
    cp_nogo = Con(expLowAspTF == 0, cp_favor_deeper)
    cp_nogo_big = Con(cp_nogo, cp_nogo, '', 'COUNT > 2')
    arcpy.BuildRasterAttributeTable_management(cp_nogo_big)
    cp_nogo_big_count = int(arcpy.management.GetCount(cp_nogo_big).getOutput(0))
    previous_stops = CellStatistics([step_deepest, eq_go1_max_fa], 'MINIMUM')
    previous_cp_nogo_big_count = cp_nogo_big_count
    previous_nogo_delta = 50
    nogo_run_start = nogo_counter
    nogo_run_count = nogo_counter - nogo_run_start
    log.warning('cp_nogo_big count is: ' + str(cp_nogo_big_count) + ' for iteration: ' + str(nogo_counter))

    while cp_nogo_big_count > 5 and nogo_run_count < 10 and previous_nogo_delta > 0:

##        log.warning(f'on nogo iteration: {nogo_counter}')
        log.warning('on nogo iteration: ' + str(nogo_counter))

        # get the max fa cells from pre_zf deepest that are the furthest along in 'good' traversal areas
        # we'll make flow from other cells aggregate to these cells
        cp_nogo_cost_dist = Con(cp_nogo, cost_dist_from_deepest)

        fs_nogo_cost_dist = FocalStatistics(cp_nogo_cost_dist, statistics_type = 'MAXIMUM')
        eq_fs_max = Con(cp_nogo_cost_dist == fs_nogo_cost_dist, 1)
        # turn these local maxima into regions and calculate their watershed
        eq_fs_max_rg = RegionGroup(eq_fs_max, 'EIGHT')
        new_cp_from_long_nogo = CostPath(eq_fs_max_rg, cost_dist_from_deepest, opj(sgdb, 'back_lnk' + str(cp_counter)))

        # get new cost paths in the nogo region
        new_cp_from_long_nogo_still_in_nogo = Con(cp_nogo, new_cp_from_long_nogo)
        new_cp_from_long_nogo_still_in_nogo_rg = RegionGroup(new_cp_from_long_nogo_still_in_nogo, 'EIGHT')
        ws_new_cp_2xnogo_rg = Watershed(stepfill2Fd, new_cp_from_long_nogo_still_in_nogo_rg)

        # get max fa cell in upper cost path segments traversing in allowed regions
        stepfill2Fa = FlowAccumulation(stepfill2Fd)
        stepfill2Fa_fs = FocalStatistics(stepfill2Fa, NbrRectangle(3,3), 'MAXIMUM')
        pre_zf_deepest_max_fa_fs = Con(pre_zf_deepest, stepfill2Fa_fs)

        # find the pre_zf_deepest cell that needs to be retained in the upper good cost path segments
        # these cells will be added to the 
        ws_nogo_max_fa = ZonalStatistics(ws_new_cp_2xnogo_rg, 'Value', pre_zf_deepest_max_fa_fs, 'MAXIMUM')
        eq_nogo_max_fa = Con(pre_zf_deepest_max_fa_fs == ws_nogo_max_fa, pre_zf_deepest)#1)
        eq_nogo_max_fa.save(opj(sgdb, 'eq_nogo_max_fa' + str(nogo_counter)))

        all_cp_stops = CellStatistics([previous_stops, eq_nogo_max_fa], 'MINIMUM')
        ##all_cp_stops = CellStatistics([step_deepest, eq_go1_max_fa, eq_nogo_max_fa, eq_nogo3_max_fa, eq_nogo4_max_fa], 'MINIMUM')

        cp_dem_2_stops, cp_dem_2_stops_fd, fill_2_stops, fill_2_stops_fd = create_cp_dem_from_stops(all_cp_stops, pre_zf_dem)

            # calculate new cost path based on dem with additional stops and all previous
            # stops that were removed as starting points
        pre_zf_deepest_that_matter = Con(IsNull(all_cp_stops), pre_zf_deepest)
        nogo_counter += 1
        nogo_run_count = nogo_counter - nogo_run_start
        cp_to_new_holes = CostPath(pre_zf_deepest_that_matter, fill_2_stops, fill_2_stops_fd)
        cp_to_new_holes.save(opj(sgdb, 'cp_2_new_holes' + str(nogo_counter)))

        cp_nogo_next = Con(expLowAspTF == 0, cp_to_new_holes)#favor_deeper)
        cp_nogo_big_next = Con(cp_nogo_next, cp_nogo_next, '', 'COUNT > 2')
        arcpy.BuildRasterAttributeTable_management(cp_nogo_big_next)
        cp_nogo_big_count_next = int(arcpy.management.GetCount(cp_nogo_big_next).getOutput(0))
        log.warning('cp_nogo_big count is: ' + str(cp_nogo_big_count_next) + ' for iteration: ' + str(nogo_counter))

        stepfill2Fd = fill_2_stops_fd
        pre_zf_dem = fill_2_stops
        cp_favor_deeper = cp_to_new_holes
        previous_stops = all_cp_stops
        previous_nogo_delta = cp_nogo_big_count - cp_nogo_big_count_next
        cp_nogo_big_count = cp_nogo_big_count_next
        cp_nogo = cp_nogo_next

    #keep all the deepest cells in the cost paths that did not get analyzed yet
    skipped_analysis_stops = Con(cp_nogo_next, pre_zf_deepest)
    all_cp_stops_final = CellStatistics([previous_stops, skipped_analysis_stops], 'MINIMUM')
    cp_dem_2_stops, cp_dem_2_stops_fd, fill_2_stops, fill_2_stops_fd = create_cp_dem_from_stops(all_cp_stops_final, pre_zf_dem)

    return cp_to_new_holes, cp_dem_2_stops, cp_dem_2_stops_fd, all_cp_stops_final, fill_2_stops, fill_2_stops_fd
##    return cp_to_new_holes, fill_2_stops, fill_2_stops_fd, all_cp_stops



def cost_path_simplify_and_starting_basins(cp, cp_dem, cp_fd, basins, pre_zf_dem):
    '''Simplify a cost path by reducing the number of starting cells by Flowlength
    trimming down/up. This results in a cost path that has more unique paths, less
    merged paths. Then remove very short cost paths that do not traverse multiple
    watersheds as they will not result in changes to the DEM. (And these short paths
    sometimes cause corruptions during inversion).
    
    cp - cost path to further simplify
    cp_dem - dem flows to all cost path stops, has no holes
    cp_fd - flow direction for above dem#to all cost path stops from dem with no holes
    basins - basins to be used for simplification, pre zonal filling steps
    pre_zf_dem - dem before zonal filling steps
    '''
    
    flCrvEqMax = startByFlowLengthDown(cp, cp_fd, basins)
    flowLengthUpStarts = startByFlowLengthUp(cp, cp_fd, basins)
    doubleFilteredStarts = Con(flCrvEqMax, flowLengthUpStarts)
    cp_to_new_holes_double = CostPath(doubleFilteredStarts, cp_dem, cp_fd)
    cp_to_new_holes_double_rg = RegionGroup(Con(cp_to_new_holes_double, 1), 'EIGHT')
    rg_to_new_holes_double_basin_var = ZonalStatistics(cp_to_new_holes_double_rg, 'Value', basins, 'VARIETY')
    cp_to_new_holes_dbl_filtered = Con(rg_to_new_holes_double_basin_var, cp, '', 'VALUE > 1')

    no_pass_basins_cells, no_pass_basins_full = create_no_pass_through_basins(basins, doubleFilteredStarts, cp_to_new_holes_dbl_filtered, pre_zf_dem)

    return cp_to_new_holes_dbl_filtered, no_pass_basins_cells, no_pass_basins_full, doubleFilteredStarts



def new_cp_to_enforce_after_inversion(all_cp_stops, fixed_dem, fixed_basins):
    '''Takes a dem and creates shorter cost paths to force flow again'''
    
    # and one final attempt to get water to flow to this steps deepest cells
    dem_holes_at_deepest = Con(IsNull(all_cp_stops), fixed_dem)
    fill_dem_to_deepest = Fill(dem_holes_at_deepest)
    fd_to_deepest = FlowDirection(fill_dem_to_deepest)
    ##fa_fill_dem_to_deepest2 = FlowAccumulation(fd_fill_dem_to_deepest2)
    fixed_dem_basin_min = ZonalStatistics(fixed_basins, 'Value', fixed_dem, 'MINIMUM')
    eq_fixed_dem_basin_min = Con(fixed_dem == fixed_dem_basin_min, fixed_basins)
    cp_for_next_fix = CostPath(eq_fixed_dem_basin_min, fill_dem_to_deepest, fd_to_deepest)

    return cp_for_next_fix, fill_dem_to_deepest, fd_to_deepest


def doCleaner(fillTif, voidFixTif, roadsFc, voidProc, xElevFile, yElevFile):
    try:
        huc12, huc8, proc_size = df.figureItOut(fillTif)

        if cleanup:
            # log to file only
            log, nowYmd, logName, startTime = df.setupLoggingNoCh(platform.node(), sys.argv[0], huc12, '')#'_' + version)
            arcpy.SetLogHistory(False)
        else:
            # log to file and console
            log, nowYmd, logName, startTime = df.setupLoggingNew(platform.node(), sys.argv[0], huc12, '')#'_' + version)

        log.info(outputString)

        ##try:
        arcpy.CheckOutExtension('Spatial')
        arcpy.env.overwriteOutput = True

        ## Set the environments
        # control where scratchFolder and GDB are created
        arcpy.env.scratchWorkspace = voidProc
        sfldr = arcpy.env.scratchFolder
        sgdb = arcpy.env.scratchGDB
        arcpy.env.scratchWorkspace = sfldr
        ##    arcpy.env.workspace = sgdb
        arcpy.env.workspace = sfldr

        arcpy.env.snapRaster = fillTif#snapRaster

        arcpy.env.cellSize = proc_size
        arcpy.env.snapRaster = fillTif
        arcpy.env.extent = fillTif

        log.warning('sys.argv is: ' + str(sys.argv) + '\n')

        #voidFixTif is a hydro-flattened DEM generated using automated flattening processes
        if arcpy.Exists(voidFixTif):
            initDEM = Raster(voidFixTif)#bridgeFixTif)

        else:
            initDEM = Raster(fillTif)
            # expLowAspTF = Raster(opj(voidProc, 'crv_tf'))

        # update the crv TF raster to use bridgeFix curvature
        meterBfDEM = 0.01 * initDEM
        pro_crv_name = "pro_crv5"
        crv2 = Curvature(meterBfDEM, '', opj(sgdb, pro_crv_name), opj(sgdb, "pln_crv5"))
##        crv2.save(opj(sgdb, "crv5"))
        proCrv2 = Float(Raster(opj(sgdb, pro_crv_name)))

        aspect, tpi, aspLow, expLowAspTF = df.channelized_areas(meterBfDEM, sgdb, proCrv2)


        counter = 0
        sfx = '_' + str(counter)

        # generate initial watersheds from DEM
        fd = FlowDirection(initDEM)
        basin0 = Basin(fd)
        basin0_rtp = arcpy.conversion.RasterToPolygon(basin0, opj(sgdb, 'basin0'), 'NO_SIMPLIFY')

        basin0Deepest = ZonalStatistics(basin0, 'VALUE', initDEM, 'MINIMUM')
        basin0DeepestCell = Con(basin0Deepest == initDEM, initDEM)

        initDemWHoles = Con(IsNull(basin0DeepestCell) == 1, initDEM)
        fill2InitHoles = Fill(initDemWHoles)
        fillInitDEM = Con(IsNull(fill2InitHoles), initDEM, fill2InitHoles)

        # set a threshold to define how wide an overflow path is allowed to be
        # overall goal is to keep zonal filling watersheds until overflow path is thicker than threshold
        thickThreshold = proc_size / 1.3 #allow for one diagonal fill cell #3.0

        Clip = opj(sgdb, 'buf_' + huc12 + '_1km')
##        if not arcpy.Exists(Clip):
##            huc12Lyr = arcpy.MakeFeatureLayer_management(HUC12FC, 'HUC12FCLayer', '"HUC12" = \'' + huc12 + "'")
##            Clip = arcpy.Buffer_analysis(huc12Lyr, opj(gdb, 'buf_' + huc12 + '_1km'), '1000 METER')
        inm = 'in_memory'
        hucRoadsIncService = arcpy.Clip_analysis(roadsFc, Clip, opj(inm, "roads_all"))
        if 'fclass' in df.getfields(hucRoadsIncService):
            hucRoadsWhole = arcpy.Select_analysis(hucRoadsIncService, opj(inm, 'roads'), 'fclass <> \'service\' AND fclass <> \'path\'AND fclass <> \'cycleway\'AND fclass <> \'footway\'')
        else:
            hucRoadsWhole = arcpy.Select_analysis(hucRoadsIncService, opj(inm, 'roads'), 'type <> \'service\' AND type <> \'path\'')

        hucRoads = arcpy.FeatureToLine_management(hucRoadsWhole, opj(inm, 'roads_broken'))

        hucRoadsRaster = arcpy.PolylineToRaster_conversion(hucRoads, 'oneway', opj(sgdb, 'roads_rast'), cellsize = proc_size)


    ##-----------------------------------------------------------------------------
        # process fill regions that are very thin and shallow

        stepDEM_thin = fillInitDEM#initDEM
        stepDeepestCellEl_thin = basin0DeepestCell
        log.warning(arcpy.GetCount_management(basin0))
        stepWs_thin = basin0

        #set up a change threshold, less than 1% difference to stop
        currentCount = int(arcpy.GetCount_management(stepWs_thin).getOutput(0))
        previousCount = currentCount*2
        log.warning('previousCount: ' + str(previousCount) + ' and currentCount: ' + str(currentCount))
        while previousCount > currentCount * 1.1 and counter < 1:
    ####    while counter < :

            previousCount = int(arcpy.GetCount_management(stepWs_thin).getOutput(0))
            
            # zonal fill the initial watersheds
            # baseDEM is stepDEM for this run
            stepZfDem_thin, stepWsZThick_thin, stepZf_thin, stepZfDif_thin = step_zonal_fill(stepWs_thin, stepDEM_thin, fillInitDEM, counter)

            # create new temporary ZF DEM and figure out where water will flow
            step_tu_dem_thin, stepTuBasins_thin = create_thin_update_dem(stepWsZThick_thin, stepZfDem_thin, stepDEM_thin, counter, thickThreshold)

            # recover initial DEM around initial sinks of each new watershed
            #   need to keep each new sink (deep zonal fill area) separate, can't just route to deepest sink in tuws
            #   and must remember that deepest zonal filled area does not always have the deepest sink
            # first find deepest cells in each tuws
            stepNoHolesAtDeepestDEM_thin, stepDeepestCellEl_thin, stepDEM_thin, counter = create_dem_filled_to_holes(stepTuBasins_thin, stepDeepestCellEl_thin, fillInitDEM, counter)

            stepfill2Fd_thin = FlowDirection(stepNoHolesAtDeepestDEM_thin)#fill2Deepest)

            stepWs_thin = Basin(stepfill2Fd_thin)
            stepWs_thin.save('step_ws' + str(counter))
            currentCount = int(arcpy.GetCount_management(stepWs_thin).getOutput(0))
            log.warning('previousCount: ' + str(previousCount) + ' and currentCount: ' + str(currentCount))


        cp_after_thin = CostPath(basin0DeepestCell, stepDEM_thin, stepfill2Fd_thin)

        # now remove cost paths that cross the road
        raster_2erase_cp = hucRoadsRaster
        cp_crossing_stop_thin = Con(cp_after_thin, raster_2erase_cp)

        cp_crossing_stop_thin_rg = RegionGroup(cp_crossing_stop_thin, 'EIGHT')
        ws_cp_crossing_stop_thin = Watershed(stepfill2Fd_thin, cp_crossing_stop_thin_rg)
        zs_min_ws_cp = ZonalStatistics(ws_cp_crossing_stop_thin, 'Value', basin0DeepestCell, 'MINIMUM')
        extra_cp_stops = Con(zs_min_ws_cp == basin0DeepestCell, basin0DeepestCell)
        all_cp_stops_thin = CellStatistics([stepDeepestCellEl_thin, extra_cp_stops], 'MINIMUM')

        dem_holes_all_cp_stops_thin = Con(IsNull(all_cp_stops_thin), fillInitDEM)
        fill_2_all_cp_stops_thin = Fill(dem_holes_all_cp_stops_thin)
        dem_no_holes_2_all_cp_stops_thin = Con(IsNull(fill_2_all_cp_stops_thin), all_cp_stops_thin, fill_2_all_cp_stops_thin)
        dem_no_holes_2_all_cp_stops_thin_fd = FlowDirection(dem_no_holes_2_all_cp_stops_thin)
        cpAtEnd_thin = CostPath(basin0DeepestCell, dem_no_holes_2_all_cp_stops_thin, dem_no_holes_2_all_cp_stops_thin_fd)

        flEqMax_thin = startByFlowLengthDown(cpAtEnd_thin, stepfill2Fd_thin, basin0)#cost_path, flow_direction, basin):#, proc_size):
        flowLenUpEq0AtAllStart_thin = startByFlowLengthUp(cpAtEnd_thin, stepfill2Fd_thin, basin0)
        doubleFilteredStarts_thin = Con(flEqMax_thin, flowLenUpEq0AtAllStart_thin)

        cp_after_double_thin = CostPath(doubleFilteredStarts_thin, dem_no_holes_2_all_cp_stops_thin, dem_no_holes_2_all_cp_stops_thin_fd)

        starting_basins_thin, starting_full_basins_thin = create_no_pass_through_basins(basin0, doubleFilteredStarts_thin, cp_after_double_thin, fillInitDEM)

        startingBasinsMinSl_thin = startingBasinsByStreamLink(doubleFilteredStarts_thin, dem_no_holes_2_all_cp_stops_thin, dem_no_holes_2_all_cp_stops_thin_fd, starting_full_basins_thin, fillInitDEM)

        invertTargetDEM_thin = df.createInvertDEM(initDEM)
        fenceEl_thin = df.getFenceEl(invertTargetDEM_thin)
        ndRcls_thin, expFix_thin, ndPlus_thin = df.setupInversionZones(cp_after_double_thin)
        fixedThinDEM = fixByInversionByStartingPath(ndPlus_thin, fenceEl_thin, invertTargetDEM_thin, doubleFilteredStarts_thin, ndRcls_thin, initDEM, startingBasinsMinSl_thin, log)

        fixDif_thin = initDEM - fixedThinDEM
        fix_dif_pos = Con(fixDif_thin > 0, fixDif_thin)
        super_shallow_threshold = 9/2
        fix_dif_pos_as_1 = Con(fixDif_thin > 0, 1)
        fix_dif_pos_rg = RegionGroup(fix_dif_pos_as_1, 'EIGHT')
        max_fix_dif = ZonalStatistics(fix_dif_pos_rg, 'Value', fix_dif_pos, 'MAXIMUM')
        thin_shallow_fixes = Con(max_fix_dif < super_shallow_threshold, fix_dif_pos)
        thin_shallow_fixed_DEM = Con(IsNull(thin_shallow_fixes), initDEM, fixedThinDEM)
        thin_shallow_fixed_DEM.save(xElevFile)
        thin_shallow_fixed_fd = FlowDirection(thin_shallow_fixed_DEM)
        thin_shallow_fixed_basins = Basin(thin_shallow_fixed_fd)
        thin_shallow_fixed_basins_rtp = arcpy.RasterToPolygon_conversion(thin_shallow_fixed_basins, opj(sgdb, 'thin_shallow_fixed_basins'), 'NO_SIMPLIFY')
        thin_shallow_count = int(arcpy.GetCount_management(thin_shallow_fixed_basins).getOutput(0))
        log.warning('thin_shallow_count: ' + str(thin_shallow_count))
        # <Result '16860'>

        thin_shallow_fixed_dif = thin_shallow_fixed_DEM - initDEM

    ##-----------------------------------------------------------------------------
        # process fill regions that are one cell thicker and in high curvature

        log.warning('now moving on to high curvature')
        # now moving on to high curvature
        thin_shallow_fixed_basins_deepest_el = ZonalStatistics(thin_shallow_fixed_basins, 'Value', initDEM, 'MINIMUM')
        thin_shallow_fixed_deepest_cell_el = Con(thin_shallow_fixed_basins_deepest_el == initDEM, initDEM)
        tcounter = 10
        stepDEM_crv = thin_shallow_fixed_DEM
        stepWs_crv = thin_shallow_fixed_basins
        stepDeepestCellEl_crv = thin_shallow_fixed_deepest_cell_el
        log.warning(arcpy.GetCount_management(stepWs_crv))
        #set up a change threshold, less than 1% difference to stop
        currentCount = int(arcpy.GetCount_management(stepWs_crv).getOutput(0))
        previousCount = currentCount*2
        log.warning('previousCount: ' + str(previousCount) + ' and currentCount: ' + str(currentCount))
        thickThreshold += proc_size #increase by one cell in size, e.g. 3.0/1.3 + 3.0
        while previousCount > currentCount * 1.01 and tcounter < 20:
    ##    while tcounter < previousCount > currentCount * 1.01 and tcounter < 20:
    ##    while tcounter < 11:
                    
            previousCount = int(arcpy.GetCount_management(stepWs_crv).getOutput(0))

            stepZfDem_crv, stepWsZThick_crv, stepZf_crv, stepZfDif_crv = step_zonal_fill(stepWs_crv, stepDEM_crv, thin_shallow_fixed_DEM, tcounter)
            stepZfDem_crv.save(opj(sgdb, 'zfdem_' + str(tcounter)))
            stepWsZThick_crv.save(opj(sgdb, 'zfthick_' + str(tcounter)))
            stepZf_crv.save(opj(sgdb, 'zf_' + str(tcounter)))
    ##        stepZfDif_crv.save(opj(sgdb, 'zfdif_' + str(tcounter)))

            # calculate compactness for greater thickness fill regions
    ####        crv_zf_dif_pos = Con(stepZfDif_crv > 0, stepWs_crv)#fixedDEMCrvFix2_basin)
            crv_zf_dif_max = ZonalStatistics(stepWs_crv, 'Value', stepZfDif_crv, 'MAXIMUM')
            depth_threshold = 75
            crv_zf_dif_TF = Con(crv_zf_dif_max < depth_threshold, 1, 0)

            crv_zf_thick_TF = stepWsZThick_crv < thickThreshold
            crv_ws_zf_true = crv_zf_dif_TF & crv_zf_thick_TF
            crv_ws_zf_true.save(opj(sgdb, 'crv_zf_TF' + str(tcounter)))

    ####        zf_crv_TF = stepWsZThick_crv < thickThreshold
    ####        zf_crv_TF.save(opj(sgdb, 'zf_crv_TF' + str(tcounter)))
            stepThinUpdateDEM_crv, stepTuBasins_crv = create_zf_update_dem(crv_ws_zf_true, stepZfDem_crv, stepDEM_crv, tcounter)#, thickThreshold)

            stepDEM_crv, stepDeepestCellEl_crv, stepHolesAtDeepestDEM_crv, tcounter = create_dem_filled_to_holes(stepTuBasins_crv, stepDeepestCellEl_crv, thin_shallow_fixed_DEM, tcounter)#, stepDEM_crv, tcounter)
            stepDEM_crv.save(opj(sgdb, 'stp_dem_' + str(tcounter)))
            stepDeepestCellEl_crv.save(opj(sgdb, 'stp_dpst_' + str(tcounter)))
            stepfill2Fd_crv = FlowDirection(stepDEM_crv)
            # these not really needed

            stepWs_crv = Basin(stepfill2Fd_crv)
            stepWs_crv.save('step_ws' + str(tcounter))
            currentCount = int(arcpy.GetCount_management(stepWs_crv).getOutput(0))
            log.warning('previousCount: ' + str(previousCount) + ' and currentCount: ' + str(currentCount))



        stepTuBasins_min = ZonalStatistics(stepTuBasins_crv, 'Value', stepDEM_crv, 'MINIMUM')
        #find basins that have multiple minima 
        eq_step_tu_basin_min_ws = Con(stepTuBasins_min == stepDEM_crv, stepWs_crv)
        stepTuBasins_min_var = ZonalStatistics(stepTuBasins_crv, 'Value', eq_step_tu_basin_min_ws, 'VARIETY')
        log.warning("try a new way of identifying low cell doubles")
        # try a new way of identifying low cell doubles
        stepDEM_crv_go = Con(expLowAspTF, stepDEM_crv)
##        stepDEM_crv_go = Con(expLowAsp1, stepDEM_crv)
        stepTuBasins_go_min = ZonalStatistics(stepTuBasins_crv, 'Value', stepDEM_crv_go, 'MINIMUM')
        stepTuBasins_go_min_ws = Con(stepTuBasins_go_min == stepDEM_crv_go, stepWs_crv)
        stepTuBasins_min_var = ZonalStatistics(stepTuBasins_crv, 'Value', stepTuBasins_go_min_ws, 'VARIETY')
        stepTuBasins_min_var_gt1 = Con(stepTuBasins_min_var > 1, stepTuBasins_crv)
        var_gt1_go = Con(stepTuBasins_min_var_gt1, expLowAspTF)
##        var_gt1_go = Con(stepTuBasins_min_var_gt1, expLowAsp1)
        var_gt1_go_rg = RegionGroup(var_gt1_go, 'EIGHT')
        min_var_gt1_go_rg = Con(eq_step_tu_basin_min_ws, var_gt1_go_rg)
        min_go_rg_var = ZonalStatistics(stepTuBasins_min_var_gt1, 'VALUE', min_var_gt1_go_rg, 'VARIETY')
        min_go_rg_var_eq1 = Con(min_go_rg_var == 1, stepTuBasins_min_var_gt1)
        minDif = stepDEM_crv_go - stepTuBasins_min
        stepTuBasins_min_var_eq1 = Con(min_go_rg_var_eq1, stepTuBasins_min)
        minDif = stepDEM_crv_go - stepTuBasins_min_var_eq1
        minDifPlus1 = minDif + 1
        stepTuBasins_min_ws = ZonalStatistics(stepTuBasins_min_var_eq1, 'VALUE', eq_step_tu_basin_min_ws, 'MINIMUM')
        eq_stepTuBasins_min_ws = Con(eq_step_tu_basin_min_ws == stepTuBasins_min_ws, eq_step_tu_basin_min_ws)
        cost_dist_from_min = CostDistance(eq_stepTuBasins_min_ws, minDifPlus1, out_backlink_raster = opj(sgdb, 'min_back_link_crv'))
        not_eq_stepTuBasins_min_ws = Con(eq_step_tu_basin_min_ws != stepTuBasins_min_ws, eq_step_tu_basin_min_ws)
        cp_from_other_minima = CostPath(not_eq_stepTuBasins_min_ws, cost_dist_from_min, opj(sgdb, 'min_back_link_crv'))
        cp_new_elevation = Con(cp_from_other_minima, stepTuBasins_min)


        cp_favor_deeper_crv, cost_dist_from_deepest_crv, dif_plus1 = get_cost_path_favor_deeper(thin_shallow_fixed_deepest_cell_el, stepfill2Fd_crv, thin_shallow_fixed_DEM, stepDEM_crv, stepWs_crv, stepHolesAtDeepestDEM_crv, stepDeepestCellEl_crv, thin_shallow_fixed_dif, tcounter)

        cp_to_new_holes_crv, dem_cp_stops_crv, fd_cp_stops_crv, all_cp_stops_crv, fill_2_cp_stops_crv, fill_2_cp_stops_crv_fd = get_cost_path_with_new_breaks(
            expLowAspTF, hucRoadsRaster, cp_favor_deeper_crv, cost_dist_from_deepest_crv, dif_plus1, tcounter, stepfill2Fd_crv, thin_shallow_fixed_deepest_cell_el, stepDeepestCellEl_crv, thin_shallow_fixed_DEM, log)

        cp_simple_crv, no_pass_basins_crv, no_pass_basins_crv_full, doublefiltered_crv = cost_path_simplify_and_starting_basins(cp_to_new_holes_crv, dem_cp_stops_crv, fd_cp_stops_crv, thin_shallow_fixed_basins, thin_shallow_fixed_DEM)

        startingStreamLinksForCrvPossibleErrorFix = startingBasinsByStreamLink(doublefiltered_crv, dem_cp_stops_crv, fd_cp_stops_crv, no_pass_basins_crv_full, thin_shallow_fixed_DEM)

        log.warning('creating invert DEM at time: ' + time.asctime())
        invertTargetDEMCrv = df.createInvertDEM(thin_shallow_fixed_DEM)
        log.warning('fencing invert DEM at time: ' + time.asctime())
        fenceElCrv = df.getFenceEl(invertTargetDEMCrv)
        step_ndRclsCrv, step_expFixCrv, step_ndPlusCrv = df.setupInversionZones(cp_simple_crv)#cp_to_new_holes_dbl_filtered)
        fixedDEMCrv_pre_new_el = fixByInversionByStartingPath(step_ndPlusCrv, fenceElCrv, invertTargetDEMCrv, no_pass_basins_crv, step_ndRclsCrv, thin_shallow_fixed_DEM, startingStreamLinksForCrvPossibleErrorFix, log)

        # bring in the new elevations to join the neighboring fill regions with the same minimum elevation
        fixedDEMCrv = Con(IsNull(cp_new_elevation), fixedDEMCrv_pre_new_el, cp_new_elevation)

        log.warning('flow direction, basins, rtp at time: ' + time.asctime())
        fixedDEMCrv_fd = FlowDirection(fixedDEMCrv)
        fixedDEMCrv_basin = Basin(fixedDEMCrv_fd)
        fixedDEMCrv_basin_rtp = arcpy.RasterToPolygon_conversion(fixedDEMCrv_basin, opj(sgdb, 'Crv_basins'), 'NO_SIMPLIFY')



        # now do a fix that will try to eliminate very short cost path segements
        # this is typically where cost path missed the deepest cell by one or two cells

        cp_after_inversion_crv2, fill_dem_after_inversion_crv2, fd_after_inversion_crv2 = new_cp_to_enforce_after_inversion(all_cp_stops_crv, fixedDEMCrv, fixedDEMCrv_basin)

        # should repeat above 
        cp_simple_crv2, no_pass_basins_crv2, no_pass_basins_crv2_full, doublefiltered_crv2 = cost_path_simplify_and_starting_basins(cp_after_inversion_crv2, fill_dem_after_inversion_crv2, fd_after_inversion_crv2, fixedDEMCrv_basin, fixedDEMCrv)

        # find the very short flowpath segments where deeper cells are just off flowpath
        cpToDeepestGT2_short = Con(cp_simple_crv2, cp_simple_crv2, '', 'VALUE > 2 AND COUNT < 3')
        cpToDeepestGT2_short_range = ZonalStatistics(cpToDeepestGT2_short, 'Value', thin_shallow_fixed_DEM, 'RANGE')
        cpToDeepest_short_slope = Con(cpToDeepestGT2_short_range < 72, cpToDeepestGT2_short)
        # now get just those short cost paths
        shortCpToDeepest = CostPath(cpToDeepest_short_slope, cost_dist_from_deepest_crv, opj(sgdb, 'back_lnk' + str(tcounter)))
        # revise that from cpToDeepest
        cpToDeepest_el = Con(shortCpToDeepest, fixedDEMCrv)

        # find the deepest elevation on these short segments
        cpToDeepestGT2_short_min = ZonalStatistics(cpToDeepest_short_slope, 'Value', fixedDEMCrv, 'MINIMUM')
        cpToDeepest_fs_min_el = FocalStatistics(cpToDeepestGT2_short_min, "RECTANGLE 3 3 CELL", 'MAXIMUM')

        # replace the elevations on cost path that are not lower the above short segments
        cp_fs_min_elev = Con(cpToDeepest_fs_min_el < cpToDeepest_el, cpToDeepest_fs_min_el, cpToDeepest_el)
        fixedDEMCrv_short_cp_fixed_DEM = Con(IsNull(cp_fs_min_elev), fixedDEMCrv, cp_fs_min_elev)
        ####    crv_start2a_no_pass_basins_full = df.fullZoneByZs(crv_start2a_no_pass_basins, fixedDEMCrv_basin)

        log.warning('creating starting basins at time: ' + time.asctime())
        sl4_crv = startingBasinsByStreamLink(doublefiltered_crv2, fill_dem_after_inversion_crv2, fd_after_inversion_crv2, no_pass_basins_crv2_full, fixedDEMCrv_short_cp_fixed_DEM)
        log.warning('creating invert DEM at time: ' + time.asctime())
        invertTargetDEMCrvFix = df.createInvertDEM(fixedDEMCrv_short_cp_fixed_DEM)
        fenceElCrvFix = df.getFenceEl(invertTargetDEMCrvFix)
        log.warning('fencing invert DEM at time: ' + time.asctime())
        step_ndRclsCrvFix, step_expFixCrvFix, step_ndPlusCrvFix = df.setupInversionZones(cp_simple_crv2)#cpForFixedCrvFix_Starts)
        fixedDEMCrvFix = fixByInversionByStartingPath(step_ndPlusCrvFix, fenceElCrvFix, invertTargetDEMCrvFix, no_pass_basins_crv2, step_ndRclsCrvFix, fixedDEMCrv_short_cp_fixed_DEM, sl4_crv, log)
        log.warning('flow direction, basins, rtp at time: ' + time.asctime())
        fixedDEMCrvFix_fd = FlowDirection(fixedDEMCrvFix)
        fixedDEMCrvFix_basin = Basin(fixedDEMCrvFix_fd)
        fixedDEMCrvFix_basin_rtp = arcpy.RasterToPolygon_conversion(fixedDEMCrvFix_basin, opj(sgdb, 'CrvFix_basins'), 'NO_SIMPLIFY')


    ##    # and one final attempt to get water to flow to this steps deepest cells
        cp_after_inversion_crv3, fill_dem_after_inversion_crv3, fd_after_inversion_crv3 = new_cp_to_enforce_after_inversion(all_cp_stops_crv, fixedDEMCrvFix, fixedDEMCrvFix_basin)

        cp_simple_crv3, no_pass_basins_crv3, no_pass_basins_crv_full3, starts_crv3 = cost_path_simplify_and_starting_basins(cp_after_inversion_crv3, fill_dem_after_inversion_crv3, fd_after_inversion_crv3, fixedDEMCrvFix_basin, fixedDEMCrvFix)


        ####    crv_start3a_no_pass_basins_full = df.fullZoneByZs(crv_start3a_no_pass_basins, fixedDEMCrvFix_basin)
        log.warning('creating starting basins at time: ' + time.asctime())
        sl3_crv = startingBasinsByStreamLink(starts_crv3, fill_dem_after_inversion_crv3, fd_after_inversion_crv3, no_pass_basins_crv_full3, fixedDEMCrvFix)
        log.warning('creating invert DEM at time: ' + time.asctime())
        invertTargetDEMCrvFix2 = df.createInvertDEM(fixedDEMCrvFix)
        log.warning('fencing invert DEM at time: ' + time.asctime())
        fenceElCrvFix2 = df.getFenceEl(invertTargetDEMCrvFix2)
        step_ndRclsCrvFix2, step_expFixCrvFix2, step_ndPlusCrvFix2 = df.setupInversionZones(cp_simple_crv3)
        fixedDEMCrvFix2 = fixByInversionByStartingPath(step_ndPlusCrvFix2, fenceElCrvFix2, invertTargetDEMCrvFix2, no_pass_basins_crv3, step_ndRclsCrvFix2, fixedDEMCrvFix, sl3_crv, log)
        fixedDEMCrvFix2.save(yElevFile)
        log.warning('flow direction, basins, rtp at time: ' + time.asctime())
        # fixedDEMCrvFix2_fd = FlowDirection(fixedDEMCrvFix2)
        # fixedDEMCrvFix2_basin = Basin(fixedDEMCrvFix2_fd)
        # fixedDEMCrvFix2_basin_rtp = arcpy.RasterToPolygon_conversion(fixedDEMCrvFix2_basin, opj(sgdb, 'CrvFix2_basins'), 'NO_SIMPLIFY')

        # crv_fill_dif = fixedDEMCrvFix2 - initDEM


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
        log.warning(pymsg + "\n")

        if arcpy.GetMessages(2) not in pymsg:
            msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
            arcpy.AddError(msgs)
            log.warning(msgs)
        sys.exit(1)

    finally:
        
        if 'logName' in locals():
            log.warning("Ending script execution at " + time.asctime())
            log.warning("Script execution lasted " + str(time.time()-startTime) + " seconds or " + str((time.time()-startTime)/60) + " minutes\n")






if __name__ == "__main__":

    outputString = 'system arguments are ' + str(sys.argv) + '\n'

    if 0 <= len(sys.argv) <= 1:
        cleanup = False
        parameters = ["C:/Program Files/ArcGIS/Pro/bin/Python/envs/arcgispro-py3/pythonw.exe",
        "O:/DEP/Scripts/basics/cmd_clean_shallow_crv.py",
        "O:/DEP/LiDAR_2013/elev_FLib_mean18_26915/07020008/ef2m070200080303.tif",
        "O:/DEP/LiDAR_2013/elev_VLib_mean18_26915/07020008/ev2m070200080303.tif",
        "O:/DEP/Basedata_Summaries/Basedata_26915.gdb/roads_merge",
        "D:/DEP_Proc/DEMProc/Void_dem2013_2m_070200080303",
        "O:/DEP/LiDAR_2013/elev_VLib_mean18_26915/07020008/ex2m070200080303.tif",
        "O:/DEP/LiDAR_2013/elev_VLib_mean18_26915/07020008/ey2m070200080303.tif"]

        for i in parameters[2:]:
            sys.argv.append(i)
    else:
        arcpy.AddMessage("Whoo, hoo! Command-line enabled!")
        # DO NOT clean up the folder after done processing - matcher needs this data
        cleanup = False

    messages = msgStub()

    # input_dem, output_dem, plib_metadata, depressions_fc, depth_threshold, area_threshold, procDir = [i for i in sys.argv[1:]]
    fillTif, voidFixTif, roadsFc, voidProc, xElevFile, yElevFile = [i for i in sys.argv[1:]]

    doCleaner(fillTif, voidFixTif, roadsFc, voidProc, xElevFile, yElevFile, cleanup, messages)
    arcpy.AddMessage("Back from doing!")
