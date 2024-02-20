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
## 2020-02-26 bkgelder - additional error handling code to handle CostDistance errors with goodSlope and cutCost
# 2021.01.07 bkgelder - updating to handl 
# 2021.02.19 bkgelder - changed dfs_2_cut to dfs2cut to mirror name in matcher

# Import system modules
import arcpy
import pickle
import sys
import string
import os
import traceback
import datetime
import time
import platform
from arcpy.sa import *
import dem_functions as df
import winsound

outputString = 'system arguments are ' + str(sys.argv) + '\n'

# this code makes it so one can copy the whole program into an ArcGIS Python Window and run it
# It populates the system arguments just like it was called from the command line
if len(sys.argv) == 1:
    cleanup = False
    parameters = ["O:/DEP/LiDAR_2013/elev_FLib_mean18_26915/10230006/ef3m102300060602.tif",
	"O:/DEP/LiDAR_2013/elev_PLib_mean18_26915/10230006/ep3m102300060602.tif",
	"O:/DEP/LiDAR_2013/elev_CLib_mean18_26915/10230006/ec3m102300060602.tif",
	"D:/DEP_Proc/DEMProc/Cut_dem2013_102300060602",
	"D:/DEP_Proc/DEMProc/Cut_dem2013_102300060602/search_102300060602_mean18.pkl",
	"O:/DEP/Basedata_Summaries/Basedata_26915.gdb/roads_merge",
	"D:/DEP/Man_Data_ACPF/dep_ACPF2020/10230006/idepACPF102300060602.gdb/cuts_prelim_mean18_dem2013_3m_102300060602",
	"D:/DEP/Man_Data_ACPF/dep_ACPF2020/10230006/idepACPF102300060602.gdb/cuts_final_mean18_dem2013_3m_102300060602",
	"D:/DEP/Man_Data_ACPF/dep_ACPF2020/10230006/idepACPF102300060602.gdb/dprsns2cut_mean18_dem2013_3m_102300060602"]
    for i in parameters[2:]:
        sys.argv.append(i)

    outputString += 'running via shell'#; parameters were: ' + str(parameters)

else:
    cleanup = True
    outputString += 'parameters were passed in via command line'
    ## set verbose (True = all output, False = minimal output saved)
    verbose = True
    if not verbose:
        arcpy.SetLogHistory(False)
    
try:

    if cleanup:
        # log to file only
        log, nowYmd, logName, startTime = df.setupLoggingNoCh(platform.node(), sys.argv[0], huc12)
    else:
        # log to file and console
        log, nowYmd, logName, startTime = df.setupLoggingNew(platform.node(), sys.argv[0], huc12)

    # input and output elevation files
    fillTif = sys.argv[1]
    punchTif = sys.argv[2]
    cutTif = sys.argv[3]
    # local processing directory
    CutProc = sys.argv[4]
    # search distance file between matching and cutting
    searchPickleFile = sys.argv[5]
    # roads file
    hucRoads = sys.argv[6]
    # good and best cuts file for repeat hydro-enforcement
    goodcutsfc = sys.argv[7]
    bestcutsfc = sys.argv[8]
    # depressions to try and drain
    depressions2cutfc = sys.argv[9]


    huc12 = os.path.basename(fillTif)[-16:-4]
    huc8 = huc12[:8]
    # this one not necessary?
    srOutCode = os.path.split(os.path.dirname(os.path.dirname(fillTif)))[-1][-5:]
    interpType = os.path.split(os.path.dirname(os.path.dirname(fillTif)))[-1][-12:-6]
    ProcSize = int(os.path.basename(fillTif)[2:-17])

    startTime = time.time()
    log.warn("Beginning logging for script at " + str(time.asctime()))
    log.warn(outputString)
    ##log.debug('This message should go to just the log file')
    ##log.warn('This one goes to both')

    log.warn('sys.argv is: ' + str(sys.argv) + '\n')

# Create output directories
## Set the environments

    arcpy.CheckOutExtension("Spatial")
    arcpy.env.overwriteOutput = True
    arcpy.env.cellSize = ProcSize
    arcpy.env.snapRaster = fillTif
    arcpy.env.scratchWorkspace = CutProc

    sfldr = arcpy.env.scratchFolder
    sgdb = arcpy.env.scratchGDB
    arcpy.env.scratchWorkspace = sfldr#
    arcpy.env.workspace = sgdb

    gdb = sgdb + '\\'
    cp = CutProc + "\\"
    inm = 'in_memory\\'

## Set the environments
    arcpy.env.workspace = gdb

    arcpy.env.scratchWorkspace = CutProc

# set up output raster names
    ## Fill region curvature and slope criteria for selecting FRs to enforce, cut lines to keep
    frSlpCrit = 7.5
    frCrvCrit = 0.25
    RMSE = 18.0 #cm

    ofElList = ["FR_OF_EL", "LONG"]
    ofElFld = ofElList[0]
    frFld = 'FILL_RGN'
    fillLvlFld = 'FILL_LVL'

    lcpMeanCrvFld = 'LCP_FR_MEAN_CRV'
    lcpMaxSlpFld = 'LCP_FR_MAX_SLP'
    lcpMaxElFld = 'LCP_FR_MAX_EL'

    minElFld = 'FR_MIN_EL'
    cutElFld = 'FR_CUT_EL'
    filOfElFld = 'MIN_AFT_CUT'
    maxMaxBfrFld = 'MAX_MAX_bfr_dist'
    maxFrOfDistFld = 'MAX_FR_OF_DIST'
    maxFillFld = 'MAX_FILL'

    wsSearchDistFld = 'WS_SRCH_DST'

    minFrDistList = ['MIN_WS_DST', 'DOUBLE']
    minFrDistList3 = ['MIN_WS_DST3', 'DOUBLE']
    minFrDistFld = minFrDistList3[0]

    ofPassFld = 'fp_pass'

    frPctDropList = ['MAX_PCT_DROP', 'DOUBLE']
    frPctDropFld = frPctDropList[0]

    medianFrFld = "Median_FR"        

    frAllCrvFrac = ['ALL_FR_CRV_FRAC', 'DOUBLE']
    frAllCrvFracFld = frAllCrvFrac[0]
    frThknsFld = 'FR_THKNS'

    cmb_score_fld = "Combo_score"

    frDepthFld = 'FR_DEPTH'
    frAreaFld = 'FR_AREA'
    frVolFld = 'FR_VOLUME'

    minElDifFld = 'MIN_EL_DIF'
    minElDif = [minElDifFld, 'long']

    maxElDifFld = 'MAX_EL_DIF'
    maxElDif = [maxElDifFld, 'long']

    rect3x3Nbr = NbrRectangle(3, 3, 'CELL')
    bsNbr = NbrRectangle(5, 5, 'CELL')
    
    annulus1x15Nbr = NbrAnnulus(1, 3, "CELL")
    circle15Nbr = NbrCircle(15, 'CELL')

    minSummaryFld = 'MIN_MIN'
    revCutElFld = 'REV_CUT_EL'

    bsMinUpElFld = 'bs_min_el_up'
    bsMinDnElFld = 'bs_min_el_dn'

    allDnCellsFld = 'all_dn'
    wsLvlFld = 'ws_lvl'

    crestMin = 9.0 # Meters

    gridfield = 'gridcode'
    gridfield2 = 'grid_code'

    goodCutsList = []

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

    searchPickleFileObject = open(searchPickleFile, 'rb')
    maxSearchDistList = pickle.load(searchPickleFileObject)
    searchPickleFileObject.close()

    log.warn('log file is ' + logName)


    punchGdb = os.path.join(gdb, os.path.splitext(os.path.basename(punchTif))[0])
    if not arcpy.Exists(punchGdb):
        punchGdbResult = arcpy.CopyRaster_management(punchTif, punchGdb)
    punchedDEMNoHoles = Con(IsNull(punchGdb) == 1, fillTif, punchTif)
    bestestDEM = punchedDEMNoHoles#Raster(punchGdb)

    arcpy.env.extent = bestestDEM.extent

    hucRoadsRaster = arcpy.PolylineToRaster_conversion(hucRoads, 'oneway', cp + 'roads_rast', cellsize = ProcSize)

    meterDEM = 0.01 * bestestDEM
    slopePct = Slope(meterDEM, 'PERCENT_RISE')
####    flats2 = Con(slopePct == 0.0, 1, 0)
####    noInteriorFlatsDEM = Con(flats2 == 0, bestestDEM, '')

## Calculate curvature to find areas where curvature is positive (channels)
    crv = Curvature(meterDEM, '', cp + "pro_crv", cp + "pln_crv")
    proCrv = Raster(cp + 'pro_crv')

    frMaxSlopeFld = 'FR_MAX_SLP'
    frMeanSlopeFld = 'FR_MEAN_SLP'

    inmDfs2Cut = arcpy.CopyFeatures_management(gdb + 'dfs2cut_' + huc12, inm + 'dfs_2_cut_' + huc12)
##    inmDfs2Cut = arcpy.CopyFeatures_management(gdb + 'dfs2cut_' + huc12, inm + 'dfs_2_cut_' + huc12)
    goodDslvAll = arcpy.CopyFeatures_management(gdb + 'good_int_dslv_all', inm + 'good_int_dslv_all')
    goodDnDslvAll = arcpy.CopyFeatures_management(gdb + 'good_dn_pts_all', inm + 'good_dn_pts_all')
    goodUpDslvAll = arcpy.CopyFeatures_management(gdb + 'good_up_pts_all', inm + 'good_up_pts_all')



    sfx = ''
    ofSfx = ''

## Initially the cut DEM is the punched DEM
    DEMpreviousCuts = bestestDEM

    arcpy.JoinField_management(goodDslvAll, frFld, inmDfs2Cut, frFld, [cutElFld, frThknsFld, minElFld, ofElFld, frAllCrvFracFld])#, compactFld

    cutItrtr = df.intersectingFeaturesUniqueIteration4(goodDslvAll, str(ProcSize*3.0) + ' METERS', inm + 'gnt_fp', ofPassFld, frFld, minElFld)

    arcpy.JoinField_management(inmDfs2Cut, frFld, goodDslvAll, frFld, ofPassFld)
    gdbSSdfs = df.copyfc(verbose, inmDfs2Cut, gdb)
    arcpy.JoinField_management(goodUpDslvAll, frFld, goodDslvAll, frFld, ofPassFld)
    arcpy.JoinField_management(goodDnDslvAll, frFld, goodDslvAll, frFld, ofPassFld)



    # max of 18 iterations, just cause it might stop this one from crashing
    cutItrtrMax = min([cutItrtr, 9])
    for ofPass in range(0, cutItrtrMax):
        sfx = '_' + str(ofPass)
        log.debug('beginning sfx ' + sfx + ' and ofPass ' + ofSfx + ' at ' + time.asctime())
        goodDslvAllLayer = arcpy.MakeFeatureLayer_management(goodDslvAll, "rgn_Lyr" + sfx + ofSfx, '"' + ofPassFld + '" = ' + str(ofPass), gdb)
        if df.testForZero(goodDslvAllLayer):#int(arcpy.GetCount_management(goodDslvAllLayer).getOutput(0)) > 0:
            goodDslvAllLyrFc = arcpy.CopyFeatures_management(goodDslvAllLayer, inm + "dp_cls_bfr2_lyr" + sfx + ofSfx)
            df.copyfc(verbose, goodDslvAllLyrFc, gdb)

            goodDnLayer = arcpy.MakeFeatureLayer_management(goodDnDslvAll, "dn_Lyr" + sfx + ofSfx, '"' + ofPassFld + '" = ' + str(ofPass), gdb)
            goodDnLyrFc = arcpy.CopyFeatures_management(goodDnLayer, inm + "dn_lyr" + sfx + ofSfx)
            df.copyfc(verbose, goodDnLyrFc, gdb)

            goodUpLayer = arcpy.MakeFeatureLayer_management(goodUpDslvAll, "up_Lyr_prlm" + sfx + ofSfx, '"' + ofPassFld + '" = ' + str(ofPass), gdb)
            goodUpLyrFc = arcpy.CopyFeatures_management(goodUpLayer, inm + "up_lyr" + sfx + ofSfx)
            df.copyfc(verbose, goodUpLyrFc, gdb)
            log.debug('set up layers for sfx ' + sfx + ' at ' + time.asctime())

            inDEM = DEMpreviousCuts

            upCellsIter = arcpy.PointToRaster_conversion(goodUpLayer, frFld, 'up_prlm' + sfx)

            maxCostDist = max(maxSearchDistList) * 100
            distToUp = CostDistance(upCellsIter, Con(bestestDEM, 1))#, maxCostDist)

            deepCloseBfr2Fr = arcpy.PolygonToRaster_conversion(goodDslvAllLayer, frFld, cp + 'cut_fr' + sfx + ofSfx, '', '', ProcSize)

            psblNewEl3 = arcpy.PolygonToRaster_conversion(goodDslvAllLayer, minElFld, cp + 'cut_el3' + sfx + ofSfx, '', '', ProcSize)

            cells2CutTo = arcpy.PointToRaster_conversion(goodDnLayer, frFld, 'dn_cells' + sfx)
            if df.testForZero(cells2CutTo):

                goodSlope = Con(cells2CutTo, cells2CutTo)#upSlope <= slopeMean, cells2CutTo)
                goodSlope.save('gd_slp' + sfx)

##                try:
                minElOutside = ZonalStatistics(cells2CutTo, 'value', bestestDEM, 'minimum')
                minElOverall = ZonalStatistics(deepCloseBfr2Fr, 'value', minElOutside, 'minimum')
                if minElOverall.maximum is not None:
                    meanCutDif = DEMpreviousCuts - (psblNewEl3 + minElOverall)/2.0
                    
                    cutCost = Con(meanCutDif >= 0, meanCutDif, 0)#cutElDif >= 0, cutElDif, 0)

                    cutCostDist = CostDistance(goodSlope, cutCost, '', cp + 'bklink' + sfx + ofSfx)#maxCostDist

                    lcp3 = CostPath(upCellsIter, cutCostDist, cp + 'bklink' + sfx + ofSfx, path_type = 'each_zone')
        ##                                                lcp3.save(cp + 'lcp3' + sfx + ofSfx)
                    log.debug('did lcp3 for sfx ' + sfx + ' at ' + time.asctime())
                    log.debug('lcp3.max for sfx ' + sfx + ' is ' + str(lcp3.maximum))
                    log.debug('deepCloseBfr2Fr count for sfx ' + sfx + ' is ' + arcpy.GetCount_management(deepCloseBfr2Fr).getOutput(0))

                    lcpFr = Con(lcp3, deepCloseBfr2Fr)
                    log.debug('lcpFr.max for sfx ' + sfx + ' is ' + str(lcpFr.maximum))
                    # needed for subsequent line in some rasters to avoid error
                    arcpy.BuildRasterAttributeTable_management(lcpFr)
                    log.debug('lcpFr count for sfx ' + sfx + ' is ' + arcpy.GetCount_management(lcpFr).getOutput(0))
                    log.debug('lcpFr name for sfx ' + sfx + ' is ' + str(lcpFr))
        ##                                                lcpFr.save(cp + 'lcp_fr' + sfx + ofSfx)
                    lcpFrPoly = arcpy.RasterToPolyline_conversion(lcpFr, inm + 'lcp3_cuts' + sfx + ofSfx, simplify = 'NO_SIMPLIFY')
                    log.debug('did lcp3 FrPoly for sfx ' + sfx + ' at ' + time.asctime())

                    zstSrchLcpFrSlp = ZonalStatisticsAsTable(lcpFr, 'value', slopePct, inm + 'zst_srch_lcp_fr_slp1' + sfx)
                    df.addCalcJoin(lcpFrPoly, gridfield2, zstSrchLcpFrSlp, 'value', ['LCP_FR_MAX_SLP', 'DOUBLE'], '!MAX!')
                    log.debug('did lcp3 stats 1 for sfx ' + sfx + ' at ' + time.asctime())
                    zstSrchLcpFrCrv = ZonalStatisticsAsTable(lcpFr, 'value', proCrv, inm + 'zst_srch_lcp_fr_crv1' + sfx)
                    df.addCalcJoin(lcpFrPoly, gridfield2, zstSrchLcpFrCrv, 'value', ['LCP_FR_MEAN_CRV', 'DOUBLE'], '!MEAN!')
                    log.debug('did lcp3 stats 2 for sfx ' + sfx + ' at ' + time.asctime())
                    zstSrchLcpFrEl = ZonalStatisticsAsTable(lcpFr, 'value', bestestDEM, inm + 'zst_srch_lcp_fr_el1' + sfx)
                    df.addCalcJoin(lcpFrPoly, gridfield2, zstSrchLcpFrEl, 'value', [lcpMaxElFld, 'DOUBLE'], '!MAX!')
                    df.copyfc(verbose, lcpFrPoly, gdb)
                    log.debug('did lcp3 stats 3 for sfx ' + sfx + ' at ' + time.asctime())

            ## Figure out if not crossing a road is an option
                    noRoadsCutCost = Con(IsNull(hucRoadsRaster), cutCost)
                    noRoadsCutCostDist = CostDistance(cells2CutTo, noRoadsCutCost, '', cp + 'nr_bkl' + sfx)#maxCostDist
                    try:
                        lcpNr = CostPath(upCellsIter, noRoadsCutCostDist, cp + 'nr_bkl' + sfx, path_type = 'each_zone')
                        log.debug('did lcpNr for sfx ' + sfx + ' at ' + time.asctime())
                        lcpNrFr = Con(lcpNr, deepCloseBfr2Fr)
        ##                                                lcpFr.save(cp + 'lcp_fr' + sfx + ofSfx)
                        lcpNrFrPoly = arcpy.RasterToPolyline_conversion(lcpNrFr, inm + 'lcp5_cuts' + sfx + ofSfx, simplify = 'NO_SIMPLIFY')

                        zstSrchLcpNrFrSlp = ZonalStatisticsAsTable(lcpNrFr, 'value', slopePct, inm + 'zst_srch_lcp_fr_slp3' + sfx)
                        df.addCalcJoin(lcpNrFrPoly, gridfield2, zstSrchLcpNrFrSlp, 'value', ['LCP_FR_MAX_SLP', 'DOUBLE'], '!MAX!')
                        log.debug('did lcpNr stats 1 for sfx ' + sfx + ' at ' + time.asctime())
                        zstSrchLcpNrFrCrv = ZonalStatisticsAsTable(lcpNrFr, 'value', proCrv, inm + 'zst_srch_lcp_fr_crv3' + sfx)
                        df.addCalcJoin(lcpNrFrPoly, gridfield2, zstSrchLcpNrFrCrv, 'value', ['LCP_FR_MEAN_CRV', 'DOUBLE'], '!MEAN!')
                        log.debug('did lcpNr stats 2 for sfx ' + sfx + ' at ' + time.asctime())
                        zstSrchLcpNrFrEl = ZonalStatisticsAsTable(lcpNrFr, 'value', bestestDEM, inm + 'zst_srch_lcp_fr_el3' + sfx)
                        df.addCalcJoin(lcpNrFrPoly, gridfield2, zstSrchLcpNrFrEl, 'value', [lcpMaxElFld, 'DOUBLE'], '!MAX!')
                        log.debug('did lcpNr stats 3 for sfx ' + sfx + ' at ' + time.asctime())

                        df.copyfc(verbose, lcpNrFrPoly, gdb)

                ## If road cross is only option (no results otherwise), medianXFld = 1, or if dfs is non-linear use initial lcp
                        lcpNrStats = arcpy.Statistics_analysis(lcpNrFrPoly, inm + 'lcpNr_stats' + sfx, [[gridfield2, 'COUNT']], gridfield2)
                        df.copytbl(verbose, lcpNrStats, gdb)
                        df.addCalcJoin(goodDslvAllLyrFc, frFld, lcpNrStats, gridfield2, ['CNT_LCP_NR', 'LONG'], '!COUNT_' + gridfield2 + '!')
                        
                        lcpStats = arcpy.Statistics_analysis(lcpFrPoly, inm + 'lcp_stats' + sfx, [[gridfield2, 'COUNT']], gridfield2)
                        df.copytbl(verbose, lcpStats, gdb)
                        df.addCalcJoin(goodDslvAllLyrFc, frFld, lcpStats, gridfield2, ['CNT_LCP', 'LONG'], '!COUNT_' + gridfield2 + '!')

                        df.copyfc(verbose, goodDslvAllLyrFc, gdb)

                        goodLCP3CutsList = df.getFrsAsList(goodDslvAllLyrFc, frFld, '(CNT_LCP >= 0 AND CNT_LCP_NR IS NULL AND ((' + ofElFld + ' - ' + minElFld + ' > ' + str(3*RMSE) + ') OR ' + frThknsFld + ' > ' + str(ProcSize * 2) + '))')
                        goodLCP3CutsSel = df.buildSelection(goodLCP3CutsList, gridfield2)

                        goodLCP3Cuts = arcpy.Select_analysis(lcpFrPoly, inm + 'good_lcp3_cuts' + sfx + ofSfx, goodLCP3CutsSel)
                        df.copyfc(verbose, goodLCP3Cuts, gdb)

                        goodLCP12CutsList = df.getFrsAsList(goodDslvAllLyrFc, frFld, '(CNT_LCP >= 0 AND CNT_LCP_NR >=0)')
                        goodLCP12CutsSel = df.buildSelection(goodLCP12CutsList, gridfield2)

                        goodLCP12Cuts = arcpy.Select_analysis(lcpNrFrPoly, inm + 'good_lcp5_cuts' + sfx + ofSfx, goodLCP12CutsSel)
                        df.copyfc(verbose, goodLCP12Cuts, gdb)
                        log.debug('finished goodLCP12Cuts')

                        goodCuts4LvlPrlm = arcpy.Merge_management([goodLCP3Cuts, goodLCP12Cuts], inm + 'good_prlm_lcp_cuts' + sfx)
                        log.debug('finished goodCuts4LvlPrlm')

                    except arcpy.ExecuteError:#ExecuteError:
                        msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
                        arcpy.AddError(msgs)
                        log.info('Trying to handle COSTPATH FROM cells is 0 error for sfx ' + sfx)
                        log.info(msgs)
                        # intercept COSTPATH error, 'ERROR 010045: COSTPATH: The number of FROM cells is 0.'
                        if '010045' in msgs:
                            lcpStats = arcpy.Statistics_analysis(lcpFrPoly, inm + 'lcp_stats' + sfx, [[gridfield2, 'COUNT']], gridfield2)
                            df.copytbl(verbose, lcpStats, gdb)
                            df.addCalcJoin(goodDslvAllLyrFc, frFld, lcpStats, gridfield2, ['CNT_LCP', 'LONG'], '!COUNT_' + gridfield2 + '!')

                            goodLCP3CutsList = df.getFrsAsList(goodDslvAllLyrFc, frFld, '(CNT_LCP >= 0 AND ((' + ofElFld + ' - ' + minElFld + ' > ' + str(3*RMSE) + ') OR ' + frThknsFld + ' > ' + str(ProcSize * 2) + '))')
                            goodLCP3CutsSel = df.buildSelection(goodLCP3CutsList, gridfield2)

                            goodLCP3Cuts = arcpy.Select_analysis(lcpFrPoly, inm + 'good_lcp3_cuts' + sfx + ofSfx, goodLCP3CutsSel)
                            df.copyfc(verbose, goodLCP3Cuts, gdb)

                            goodCuts4LvlPrlm = arcpy.CopyFeatures_management(goodLCP3Cuts, inm + 'good_prlm_lcp_cuts' + sfx)
                            log.debug('No no-road cuts for sfx ' + sfx)

                        else:
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
                            log.warn(pymsg + "\n")

                            sys.exit(1)
                            
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
                        log.warn(pymsg + "\n")

                        if arcpy.GetMessages(2) not in pymsg:
                            msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
                            arcpy.AddError(msgs)
                            log.warn(msgs)

                        sys.exit(1)
                        
##                except:
                else:
                    log.warning('minElOverall.maximum is None, no valid cells/points to cut for sfx: ' + sfx)
                    log.warning('Typical error follows')
                    error = '''WARNING - PYTHON ERRORS:
Traceback info:
  File "O:\DEP\Scripts\basics\cmd_cutter.py", line 316, in <module>
    cutCostDist = CostDistance(goodSlope, cutCost, '', cp + 'bklink' + sfx + ofSfx)#maxCostDist

Error Info:
ERROR 999999: Error executing function.
Class not registered

ERROR 010004: All cells in Raster D:\DEP_Proc\CutProc\Cut_070101061007\ifthe_ras7 have the NODATA value. Stop execution.
ERROR 010067: Error in executing grid expression.
Failed to execute (CostDistance).'''
                    log.warning(error)

                    
####                    winsound.Beep(500, 1000)
####                    winsound.Beep(1000, 1000)
####                    # Get the traceback object
####                    tb = sys.exc_info()[2]
####                    tbinfo = traceback.format_tb(tb)[0]
####
####                    # Concatenate information together concerning the error into a message string
####                    pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
####                    # Return python error messages for use in script tool or Python Window
####                    arcpy.AddError(pymsg)
####                    # Print Python error messages for use in Python / Python Window
####                    log.warn(pymsg + "\n")
####
####                    if arcpy.GetMessages(2) not in pymsg:
####                        msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
####                        arcpy.AddError(msgs)
####                        log.warn(msgs)
####
########                    sys.exit(1)

                validCrossArea = lcpMaxSlpFld + ' >= ' + str(0.667*frSlpCrit) + ' OR ' + lcpMeanCrvFld + ' > ' + str(frCrvCrit)
                goodCuts4Lvl = arcpy.Select_analysis(goodCuts4LvlPrlm, inm + 'good_lcp_cuts' + sfx, validCrossArea)
                goodCutsList.append(goodCuts4Lvl)

                arcpy.JoinField_management(goodCuts4Lvl, gridfield2, inmDfs2Cut, frFld, cutElFld)
                DEMwGoodCutsLvl, cutRaster = df.createCLDEM(DEMpreviousCuts, gdb, goodCuts4Lvl, 'cl_lvl_dem', sfx, cutElFld, ProcSize)

                DEMpreviousCuts = DEMwGoodCutsLvl

    arcpy.env.workspace = cp
    fr0 = Raster('fr0_0')
    goodCutsAll = df.condenseDataLvls(goodCutsList, inm + 'good_cuts_all')

    arcpy.JoinField_management(goodCutsAll, gridfield2, inmDfs2Cut, frFld, cutElFld)

    DEMwGoodCutsAll, cutRaster = df.createCLDEM(bestestDEM, gdb, goodCutsAll, 'cl_dem', sfx, cutElFld, ProcSize)

    fillAfterGoodCuts = Fill(DEMwGoodCutsAll)

    frList = arcpy.ListRasters(fr0.name.split('_')[0] + '_*')
    for rast in frList:
        if rast[3] not in string.digits:
            frList.remove(rast)
    maxFr = CellStatistics(frList, 'MAXIMUM')
    frList.sort()
    for index, frl in enumerate(frList):
        sfx = '_' + str(index)
        zstFillAftGoodCuts = ZonalStatisticsAsTable(frl, 'value', fillAfterGoodCuts, inm + 'zst_fil_aft_gd' + sfx)
    zstFillAftGoodCutsAll = df.condenseTableLvls(zstFillAftGoodCuts, inm, inm + 'fil_aft_gd_all')
    df.addCalcJoin(inmDfs2Cut, frFld, zstFillAftGoodCutsAll, 'VALUE', [filOfElFld, 'LONG'], '!MIN!')
####                                df.addCalcJoin(inmDfs2Cut, frFld, drainZST, 'VALUE', [filOfElFld, 'LONG'], '!MIN!')
    arcpy.JoinField_management(goodCutsAll, gridfield2, zstFillAftGoodCutsAll, 'VALUE', filOfElFld)
    arcpy.JoinField_management(goodCutsAll, gridfield2, inmDfs2Cut, frFld, [ofElFld, minElFld, maxFillFld])
    df.copytbl(verbose, zstFillAftGoodCutsAll, gdb)
##    df.copyfc(verbose, goodCutsAll, gdb)

    goodcutsCopy = arcpy.CopyFeatures_management(goodCutsAll, goodcutsfc)

    notBetterCutFrs = df.getFrsAsList(goodCutsAll, gridfield2, filOfElFld + ' > 0.65*(' + ofElFld + ' - ' + minElFld + ') + ' + minElFld + ' AND ' + filOfElFld + ' > 0.65*(' + maxFillFld + ' - ' + minElFld + ') + ' + minElFld)

    if len(notBetterCutFrs) > 0:
        betterCutsAll = arcpy.Select_analysis(goodCutsAll, inm + 'better_cuts_all', df.buildAntiSelection(notBetterCutFrs, gridfield2))
        df.copyfc(verbose, betterCutsAll, gdb)
    else:
        betterCutsAll = arcpy.CopyFeatures_management(goodCutsAll, inm + 'better_cuts_all')
        df.copyfc(verbose, betterCutsAll, gdb)

    bestCutsAll = arcpy.CopyFeatures_management(betterCutsAll, bestcutsfc)#inm + 'best_cuts_all')

    df.addCalcJoin(inmDfs2Cut, frFld, bestCutsAll, gridfield2, ['cut_length', 'DOUBLE'], '!Shape.Length!')

    depressions2cutCopy = arcpy.CopyFeatures_management(inmDfs2Cut, depressions2cutfc)#inm + 'best_cuts_all')

    sfx = '_0'
    DEMwBestCuts, bestCutRaster = df.createCLDEM(bestestDEM, gdb, bestCutsAll, 'cl_dem_fnl', sfx, cutElFld, ProcSize)

    ## check to see if 'right' size compared to pitfill before save, use 0.95 instead of 0.97 due to voids at holes
    pfArea = ZonalGeometry(Con(IsNull(fillTif) == 0, 1), 'VALUE', 'AREA')
    cutArea = ZonalGeometry(Con(IsNull(DEMwBestCuts) == 0, 1), 'VALUE', 'AREA')
    areaRatio = cutArea.maximum / pfArea.maximum
    try:
        assert(areaRatio >= 0.9),'VoidFixed DEM failed area data check'
        DEMwBestCuts.save(cutTif)
        log.info("COMPLETION! Succesful save of cut DEM at " + time.asctime())
    except AssertionError:
        log.error('VoidFixed DEM failed area data check', exc_info = True)


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
    log.warn(pymsg + "\n")

    if arcpy.GetMessages(2) not in pymsg:
        msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
        arcpy.AddError(msgs)
        log.warn(msgs)
    sys.exit(1)

finally:
    
    if 'logName' in locals():
        log.warn("Ending script execution at " + time.asctime())
        log.warn("Script execution lasted " + str(time.time()-startTime) + " seconds or " + str((time.time()-startTime)/60) + " minutes\n")

##    log.removeHandler(fh)
##    del fh
##    log.removeHandler(ch)
##    del ch

    arcpy.ResetEnvironments()
