# 2024.02.19 bkgelder - converted to Python 3.9

# Import system modules
import arcpy
import pickle
import sys
import string
from os.path import join as opj
import traceback
import time
import platform
from arcpy.sa import *
##from dem_functions import *
sys.path.append("C:\\DEP\\Scripts\\basics")
import dem_functions as df
import winsound



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
        doCutter(parameters[0].valueAsText, parameters[1].valueAsText, parameters[2].valueAsText, parameters[3].valueAsText, parameters[4].valueAsText, parameters[5].valueAsText, parameters[6].valueAsText, parameters[7].valueAsText, parameters[8].valueAsText, parameters[9].valueAsText, cleanup, messages)
        return

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return



def doCutter(input_dem, huc_roads, dfs_2_cut_fc, good_dslv_fc, good_up_dslv_fc, good_dn_dslv_fc, search_distance_file, output_dem, good_cuts_fc, best_cuts_fc, depressions2cut_fc, proc_dir, match_depth, cleanup, messages):

    try:
        arguments = [input_dem, huc_roads, dfs_2_cut_fc, good_dslv_fc, good_up_dslv_fc, good_dn_dslv_fc, search_distance_file, output_dem, good_cuts_fc, best_cuts_fc, depressions2cut_fc, proc_dir, match_depth, cleanup]

        for a in arguments:
            if a == arguments[0]:
                arg_str = str(a) + '\n'
            else:
                arg_str += str(a) + '\n'

        messages.addMessage("Tool: Executing with parameters:\n" + arg_str)

        huc12, huc8, ProcSize = df.figureItOut(input_dem)

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

        log.info('log file is ' + logName)

        log.info("Tool: Executing with parameters:\n" + arg_str)

    # Create output directories
    ## Set the environments
        # control where scratchFolder and GDB are created
        ## Make sure output locations exist
        arcpy.env.scratchWorkspace = proc_dir

        sgdb = arcpy.env.scratchGDB
        arcpy.env.scratchWorkspace = sgdb
        arcpy.env.workspace = sgdb

        arcpy.env.snapRaster = input_dem

        arcpy.env.cellSize = input_dem

        arcpy.env.extent = input_dem

        inm = 'in_memory'

        arcpy.CheckOutExtension("Spatial")
        arcpy.env.overwriteOutput = True



    # set up output raster names
        ## Fill region curvature and slope criteria for selecting FRs to enforce, cut lines to keep
        frSlpCrit = 7.5
        frCrvCrit = 0.25

        ofElList = ["FR_OF_EL", "LONG"]
        ofElFld = ofElList[0]
        frFld = 'FILL_RGN'

        lcpMeanCrvFld = 'LCP_FR_MEAN_CRV'
        lcpMaxSlpFld = 'LCP_FR_MAX_SLP'
        lcpMaxElFld = 'LCP_FR_MAX_EL'

        minElFld = 'FR_MIN_EL'
        cutElFld = 'FR_CUT_EL'
        filOfElFld = 'MIN_AFT_CUT'
        maxFillFld = 'MAX_FILL'

        minFrDistList = ['MIN_WS_DST', 'DOUBLE']
        minFrDistList3 = ['MIN_WS_DST3', 'DOUBLE']
        minFrDistFld = minFrDistList3[0]

        ofPassFld = 'fp_pass'

        frPctDropList = ['MAX_PCT_DROP', 'DOUBLE']
        frPctDropFld = frPctDropList[0]

        frAllCrvFrac = ['ALL_FR_CRV_FRAC', 'DOUBLE']
        frAllCrvFracFld = frAllCrvFrac[0]
        frThknsFld = 'FR_THKNS'

        minElDifFld = 'MIN_EL_DIF'
        minElDif = [minElDifFld, 'long']

        maxElDifFld = 'MAX_EL_DIF'
        maxElDif = [maxElDifFld, 'long']

        gridfield = 'gridcode'
        gridfield2 = 'grid_code'

        goodCutsList = []

        searchPickleFileObject = open(search_distance_file, 'rb')
        maxSearchDistList = pickle.load(searchPickleFileObject)
        searchPickleFileObject.close()

        hucRoadsRaster = arcpy.PolylineToRaster_conversion(huc_roads, 'oneway', opj(proc_dir, 'roads_rast'), cellsize = ProcSize)

        meterDEM = 0.01 * Raster(input_dem)
        # slopePct = Slope(meterDEM, 'PERCENT_RISE')
        slopePct = Raster(opj(proc_dir, 'slope_pct'))

    ## Calculate curvature to find areas where curvature is positive (channels)
        crv = Curvature(meterDEM, '', opj(proc_dir, "pro_crv"), opj(proc_dir, "pln_crv"))
        proCrv = Raster(opj(proc_dir, 'pro_crv'))

    #     inmDfs2Cut = arcpy.CopyFeatures_management(opj(sgdb, 'dfs2cut_' + huc12), opj(inm, 'dfs_2_cut_' + huc12))
    # ##    inmDfs2Cut = arcpy.CopyFeatures_management(opj(sgdb, 'dfs2cut_' + huc12), opj(inm, 'dfs_2_cut_' + huc12))
    #     goodDslvAll = arcpy.CopyFeatures_management(opj(sgdb, 'good_int_dslv_all'), opj(inm, 'good_int_dslv_all'))
    #     goodDnDslvAll = arcpy.CopyFeatures_management(opj(sgdb, 'good_dn_pts_all'), opj(inm, 'good_dn_pts_all'))
    #     goodUpDslvAll = arcpy.CopyFeatures_management(opj(sgdb, 'good_up_pts_all'), opj(inm, 'good_up_pts_all'))
        inmDfs2Cut = arcpy.CopyFeatures_management(dfs_2_cut_fc, opj(inm, 'dfs_2_cut_' + huc12))
    ##    inmDfs2Cut = arcpy.CopyFeatures_management(opj(sgdb, 'dfs2cut_' + huc12), opj(inm, 'dfs_2_cut_' + huc12))
        goodDslvAll = arcpy.CopyFeatures_management(good_dslv_fc, opj(inm, 'good_int_dslv_all'))
        goodDnDslvAll = arcpy.CopyFeatures_management(good_dn_dslv_fc, opj(inm, 'good_dn_pts_all'))
        goodUpDslvAll = arcpy.CopyFeatures_management(good_up_dslv_fc, opj(inm, 'good_up_pts_all'))

        sfx = ''
        ofSfx = ''

    ## Initially the cut DEM is the punched DEM
        DEMpreviousCuts = input_dem

        arcpy.JoinField_management(goodDslvAll, frFld, inmDfs2Cut, frFld, [cutElFld, frThknsFld, minElFld, ofElFld, frAllCrvFracFld])#, compactFld

        cutItrtr = df.intersectingFeaturesUniqueIteration4(goodDslvAll, str(ProcSize*3.0) + ' METERS', opj(inm, 'gnt_fp'), ofPassFld, frFld, minElFld)

        arcpy.JoinField_management(inmDfs2Cut, frFld, goodDslvAll, frFld, ofPassFld)
        gdbSSdfs = df.copyfc(verbose, inmDfs2Cut, sgdb)
        arcpy.JoinField_management(goodUpDslvAll, frFld, goodDslvAll, frFld, ofPassFld)
        arcpy.JoinField_management(goodDnDslvAll, frFld, goodDslvAll, frFld, ofPassFld)


        # max of 18 iterations, just cause it might stop this one from crashing
        cutItrtrMax = min([cutItrtr, 9])
        for ofPass in range(0, cutItrtrMax):
            sfx = '_' + str(ofPass)
            log.debug('beginning sfx ' + sfx + ' and ofPass ' + ofSfx + ' at ' + time.asctime())
            goodDslvAllLayer = arcpy.MakeFeatureLayer_management(goodDslvAll, "rgn_Lyr" + sfx + ofSfx, '"' + ofPassFld + '" = ' + str(ofPass), sgdb)
            if df.testForZero(goodDslvAllLayer):#int(arcpy.GetCount_management(goodDslvAllLayer).getOutput(0)) > 0:
                goodDslvAllLyrFc = arcpy.CopyFeatures_management(goodDslvAllLayer, inm + "dp_cls_bfr2_lyr" + sfx + ofSfx)
                df.copyfc(verbose, goodDslvAllLyrFc, sgdb)

                goodDnLayer = arcpy.MakeFeatureLayer_management(goodDnDslvAll, "dn_Lyr" + sfx + ofSfx, '"' + ofPassFld + '" = ' + str(ofPass), sgdb)
                goodDnLyrFc = arcpy.CopyFeatures_management(goodDnLayer, inm + "dn_lyr" + sfx + ofSfx)
                df.copyfc(verbose, goodDnLyrFc, sgdb)

                goodUpLayer = arcpy.MakeFeatureLayer_management(goodUpDslvAll, "up_Lyr_prlm" + sfx + ofSfx, '"' + ofPassFld + '" = ' + str(ofPass), sgdb)
                goodUpLyrFc = arcpy.CopyFeatures_management(goodUpLayer, inm + "up_lyr" + sfx + ofSfx)
                df.copyfc(verbose, goodUpLyrFc, sgdb)
                log.debug('set up layers for sfx ' + sfx + ' at ' + time.asctime())

                upCellsIter = arcpy.PointToRaster_conversion(goodUpLayer, frFld, 'up_prlm' + sfx)

                # maxCostDist = max(maxSearchDistList) * 100
                # distToUp = CostDistance(upCellsIter, Con(input_dem, 1))#, maxCostDist)

                deepCloseBfr2Fr = arcpy.PolygonToRaster_conversion(goodDslvAllLayer, frFld, opj(proc_dir, 'cut_fr' + sfx + ofSfx), '', '', ProcSize)

                psblNewEl3 = arcpy.PolygonToRaster_conversion(goodDslvAllLayer, minElFld, opj(proc_dir, 'cut_el3' + sfx + ofSfx), '', '', ProcSize)

                cells2CutTo = arcpy.PointToRaster_conversion(goodDnLayer, frFld, 'dn_cells' + sfx)
                if df.testForZero(cells2CutTo):

                    goodSlope = Con(cells2CutTo, cells2CutTo)#upSlope <= slopeMean, cells2CutTo)
                    goodSlope.save('gd_slp' + sfx)

    ##                try:
                    minElOutside = ZonalStatistics(cells2CutTo, 'value', input_dem, 'minimum')
                    minElOverall = ZonalStatistics(deepCloseBfr2Fr, 'value', minElOutside, 'minimum')
                    if minElOverall.maximum is not None:
                        meanCutDif = DEMpreviousCuts - (psblNewEl3 + minElOverall)/2.0
                        
                        cutCost = Con(meanCutDif >= 0, meanCutDif, 0)#cutElDif >= 0, cutElDif, 0)

                        back_link = opj(proc_dir, 'bklink' + sfx + ofSfx)
                        cutCostDist = CostDistance(goodSlope, cutCost, '', back_link)#maxCostDist

                        lcp3 = CostPath(upCellsIter, cutCostDist, back_link, path_type = 'each_zone')
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
                        lcpFrPoly = arcpy.RasterToPolyline_conversion(lcpFr, opj(inm, 'lcp3_cuts' + sfx + ofSfx), simplify = 'NO_SIMPLIFY')
                        log.debug('did lcp3 FrPoly for sfx ' + sfx + ' at ' + time.asctime())

                        zstSrchLcpFrSlp = ZonalStatisticsAsTable(lcpFr, 'value', slopePct, opj(inm, 'zst_srch_lcp_fr_slp1' + sfx))
                        df.addCalcJoin(lcpFrPoly, gridfield2, zstSrchLcpFrSlp, 'value', ['LCP_FR_MAX_SLP', 'DOUBLE'], '!MAX!')
                        log.debug('did lcp3 stats 1 for sfx ' + sfx + ' at ' + time.asctime())
                        zstSrchLcpFrCrv = ZonalStatisticsAsTable(lcpFr, 'value', proCrv, opj(inm, 'zst_srch_lcp_fr_crv1' + sfx))
                        df.addCalcJoin(lcpFrPoly, gridfield2, zstSrchLcpFrCrv, 'value', ['LCP_FR_MEAN_CRV', 'DOUBLE'], '!MEAN!')
                        log.debug('did lcp3 stats 2 for sfx ' + sfx + ' at ' + time.asctime())
                        zstSrchLcpFrEl = ZonalStatisticsAsTable(lcpFr, 'value', input_dem, opj(inm, 'zst_srch_lcp_fr_el1' + sfx))
                        df.addCalcJoin(lcpFrPoly, gridfield2, zstSrchLcpFrEl, 'value', [lcpMaxElFld, 'DOUBLE'], '!MAX!')
                        df.copyfc(verbose, lcpFrPoly, sgdb)
                        log.debug('did lcp3 stats 3 for sfx ' + sfx + ' at ' + time.asctime())

                ## Figure out if not crossing a road is an option
                        noRoadsCutCost = Con(IsNull(hucRoadsRaster), cutCost)
                        near_back_link = opj(proc_dir, 'nr_bkl' + sfx)
                        noRoadsCutCostDist = CostDistance(cells2CutTo, noRoadsCutCost, '', near_back_link)#maxCostDist
                        try:
                            lcpNr = CostPath(upCellsIter, noRoadsCutCostDist, near_back_link, path_type = 'each_zone')
                            log.debug('did lcpNr for sfx ' + sfx + ' at ' + time.asctime())
                            lcpNrFr = Con(lcpNr, deepCloseBfr2Fr)
            ##                                                lcpFr.save(cp + 'lcp_fr' + sfx + ofSfx)
                            lcpNrFrPoly = arcpy.RasterToPolyline_conversion(lcpNrFr, opj(inm, 'lcp5_cuts' + sfx + ofSfx), simplify = 'NO_SIMPLIFY')

                            zstSrchLcpNrFrSlp = ZonalStatisticsAsTable(lcpNrFr, 'value', slopePct, opj(inm, 'zst_srch_lcp_fr_slp3' + sfx))
                            df.addCalcJoin(lcpNrFrPoly, gridfield2, zstSrchLcpNrFrSlp, 'value', ['LCP_FR_MAX_SLP', 'DOUBLE'], '!MAX!')
                            log.debug('did lcpNr stats 1 for sfx ' + sfx + ' at ' + time.asctime())
                            zstSrchLcpNrFrCrv = ZonalStatisticsAsTable(lcpNrFr, 'value', proCrv, opj(inm, 'zst_srch_lcp_fr_crv3' + sfx))
                            df.addCalcJoin(lcpNrFrPoly, gridfield2, zstSrchLcpNrFrCrv, 'value', ['LCP_FR_MEAN_CRV', 'DOUBLE'], '!MEAN!')
                            log.debug('did lcpNr stats 2 for sfx ' + sfx + ' at ' + time.asctime())
                            zstSrchLcpNrFrEl = ZonalStatisticsAsTable(lcpNrFr, 'value', input_dem, opj(inm, 'zst_srch_lcp_fr_el3' + sfx))
                            df.addCalcJoin(lcpNrFrPoly, gridfield2, zstSrchLcpNrFrEl, 'value', [lcpMaxElFld, 'DOUBLE'], '!MAX!')
                            log.debug('did lcpNr stats 3 for sfx ' + sfx + ' at ' + time.asctime())

                            df.copyfc(verbose, lcpNrFrPoly, sgdb)

                    ## If road cross is only option (no results otherwise), medianXFld = 1, or if dfs is non-linear use initial lcp
                            lcpNrStats = arcpy.Statistics_analysis(lcpNrFrPoly, opj(inm, 'lcpNr_stats' + sfx), [[gridfield2, 'COUNT']], gridfield2)
                            df.copytbl(verbose, lcpNrStats, sgdb)
                            df.addCalcJoin(goodDslvAllLyrFc, frFld, lcpNrStats, gridfield2, ['CNT_LCP_NR', 'LONG'], '!COUNT_' + gridfield2 + '!')
                            
                            lcpStats = arcpy.Statistics_analysis(lcpFrPoly, opj(inm, 'lcp_stats' + sfx), [[gridfield2, 'COUNT']], gridfield2)
                            df.copytbl(verbose, lcpStats, sgdb)
                            df.addCalcJoin(goodDslvAllLyrFc, frFld, lcpStats, gridfield2, ['CNT_LCP', 'LONG'], '!COUNT_' + gridfield2 + '!')

                            df.copyfc(verbose, goodDslvAllLyrFc, sgdb)

                            goodLCP3CutsList = df.getFrsAsList(goodDslvAllLyrFc, frFld, '(CNT_LCP >= 0 AND CNT_LCP_NR IS NULL AND ((' + ofElFld + ' - ' + minElFld + ' > ' + str(3*float(match_depth)) + ') OR ' + frThknsFld + ' > ' + str(ProcSize * 2) + '))')
                            goodLCP3CutsSel = df.buildSelection(goodLCP3CutsList, gridfield2)

                            goodLCP3Cuts = arcpy.Select_analysis(lcpFrPoly, opj(inm, 'good_lcp3_cuts' + sfx + ofSfx), goodLCP3CutsSel)
                            df.copyfc(verbose, goodLCP3Cuts, sgdb)

                            goodLCP12CutsList = df.getFrsAsList(goodDslvAllLyrFc, frFld, '(CNT_LCP >= 0 AND CNT_LCP_NR >=0)')
                            goodLCP12CutsSel = df.buildSelection(goodLCP12CutsList, gridfield2)

                            goodLCP12Cuts = arcpy.Select_analysis(lcpNrFrPoly, opj(inm, 'good_lcp5_cuts' + sfx + ofSfx), goodLCP12CutsSel)
                            df.copyfc(verbose, goodLCP12Cuts, sgdb)
                            log.debug('finished goodLCP12Cuts')

                            goodCuts4LvlPrlm = arcpy.Merge_management([goodLCP3Cuts, goodLCP12Cuts], opj(inm, 'good_prlm_lcp_cuts' + sfx))
                            log.debug('finished goodCuts4LvlPrlm')

                        except arcpy.ExecuteError:#ExecuteError:
                            msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
                            arcpy.AddError(msgs)
                            log.info('Trying to handle COSTPATH FROM cells is 0 error for sfx ' + sfx)
                            log.info(msgs)
                            # intercept COSTPATH error, 'ERROR 010045: COSTPATH: The number of FROM cells is 0.'
                            if '010045' in msgs:
                                lcpStats = arcpy.Statistics_analysis(lcpFrPoly, opj(inm, 'lcp_stats' + sfx), [[gridfield2, 'COUNT']], gridfield2)
                                df.copytbl(verbose, lcpStats, sgdb)
                                df.addCalcJoin(goodDslvAllLyrFc, frFld, lcpStats, gridfield2, ['CNT_LCP', 'LONG'], '!COUNT_' + gridfield2 + '!')

                                goodLCP3CutsList = df.getFrsAsList(goodDslvAllLyrFc, frFld, '(CNT_LCP >= 0 AND ((' + ofElFld + ' - ' + minElFld + ' > ' + str(3*float(match_depth)) + ') OR ' + frThknsFld + ' > ' + str(ProcSize * 2) + '))')
                                goodLCP3CutsSel = df.buildSelection(goodLCP3CutsList, gridfield2)

                                goodLCP3Cuts = arcpy.Select_analysis(lcpFrPoly, opj(inm, 'good_lcp3_cuts' + sfx + ofSfx), goodLCP3CutsSel)
                                df.copyfc(verbose, goodLCP3Cuts, sgdb)

                                goodCuts4LvlPrlm = arcpy.CopyFeatures_management(goodLCP3Cuts, opj(inm, 'good_prlm_lcp_cuts' + sfx))
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
                                log.warning(pymsg + "\n")

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
                            log.warning(pymsg + "\n")

                            if arcpy.GetMessages(2) not in pymsg:
                                msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
                                arcpy.AddError(msgs)
                                log.warning(msgs)

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
    ####                    log.warning(pymsg + "\n")
    ####
    ####                    if arcpy.GetMessages(2) not in pymsg:
    ####                        msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
    ####                        arcpy.AddError(msgs)
    ####                        log.warning(msgs)
    ####
    ########                    sys.exit(1)

                    validCrossArea = lcpMaxSlpFld + ' >= ' + str(0.667*frSlpCrit) + ' OR ' + lcpMeanCrvFld + ' > ' + str(frCrvCrit)
                    goodCuts4Lvl = arcpy.Select_analysis(goodCuts4LvlPrlm, opj(inm, 'good_lcp_cuts' + sfx), validCrossArea)
                    goodCutsList.append(goodCuts4Lvl)

                    arcpy.JoinField_management(goodCuts4Lvl, gridfield2, inmDfs2Cut, frFld, cutElFld)
                    DEMwGoodCutsLvl, cutRaster = df.createCLDEM(DEMpreviousCuts, sgdb, goodCuts4Lvl, 'cl_lvl_dem', sfx, cutElFld, ProcSize)

                    DEMpreviousCuts = DEMwGoodCutsLvl

        arcpy.env.workspace = proc_dir
        goodCutsAll = df.condenseDataLvls(goodCutsList, opj(inm, 'good_cuts_all'))

        arcpy.JoinField_management(goodCutsAll, gridfield2, inmDfs2Cut, frFld, cutElFld)

        DEMwGoodCutsAll, cutRaster = df.createCLDEM(input_dem, sgdb, goodCutsAll, 'cl_dem', sfx, cutElFld, ProcSize)

        fillAfterGoodCuts = Fill(DEMwGoodCutsAll)

        fr0 = Raster('fr0_0')
        frList = arcpy.ListRasters(fr0.name.split('_')[0] + '_*')
        for rast in frList:
            if rast[3] not in string.digits:
                frList.remove(rast)
        maxFr = CellStatistics(frList, 'MAXIMUM')
        frList.sort()
        for index, frl in enumerate(frList):
            sfx = '_' + str(index)
            zstFillAftGoodCuts = ZonalStatisticsAsTable(frl, 'value', fillAfterGoodCuts, opj(inm, 'zst_fil_aft_gd' + sfx))
        zstFillAftGoodCutsAll = df.condenseTableLvls(zstFillAftGoodCuts, inm, opj(inm, 'fil_aft_gd_all'))
        df.addCalcJoin(inmDfs2Cut, frFld, zstFillAftGoodCutsAll, 'VALUE', [filOfElFld, 'LONG'], '!MIN!')
    ####                                df.addCalcJoin(inmDfs2Cut, frFld, drainZST, 'VALUE', [filOfElFld, 'LONG'], '!MIN!')
        arcpy.JoinField_management(goodCutsAll, gridfield2, zstFillAftGoodCutsAll, 'VALUE', filOfElFld)
        arcpy.JoinField_management(goodCutsAll, gridfield2, inmDfs2Cut, frFld, [ofElFld, minElFld, maxFillFld])
        df.copytbl(verbose, zstFillAftGoodCutsAll, sgdb)
    ##    df.copyfc(verbose, goodCutsAll, sgdb)

        goodcutsCopy = arcpy.CopyFeatures_management(goodCutsAll, good_cuts_fc)

        notBetterCutFrs = df.getFrsAsList(goodCutsAll, gridfield2, filOfElFld + ' > 0.65*(' + ofElFld + ' - ' + minElFld + ') + ' + minElFld + ' AND ' + filOfElFld + ' > 0.65*(' + maxFillFld + ' - ' + minElFld + ') + ' + minElFld)

        if len(notBetterCutFrs) > 0:
            betterCutsAll = arcpy.Select_analysis(goodCutsAll, opj(inm, 'better_cuts_all'), df.buildAntiSelection(notBetterCutFrs, gridfield2))
            df.copyfc(verbose, betterCutsAll, sgdb)
        else:
            betterCutsAll = arcpy.CopyFeatures_management(goodCutsAll, opj(inm, 'better_cuts_all'))
            df.copyfc(verbose, betterCutsAll, sgdb)

        bestCutsAll = arcpy.CopyFeatures_management(betterCutsAll, best_cuts_fc)#inm + 'best_cuts_all')

        df.addCalcJoin(inmDfs2Cut, frFld, bestCutsAll, gridfield2, ['cut_length', 'DOUBLE'], '!Shape.Length!')

        depressions2cutCopy = arcpy.CopyFeatures_management(inmDfs2Cut, depressions2cut_fc)#inm + 'best_cuts_all')

        sfx = '_0'
        DEMwBestCuts, bestCutRaster = df.createCLDEM(input_dem, sgdb, bestCutsAll, 'cl_dem_fnl', sfx, cutElFld, ProcSize)

        ## check to see if 'right' size compared to pitfill before save, use 0.95 instead of 0.97 due to voids at holes
        pfArea = ZonalGeometry(Con(IsNull(input_dem) == 0, 1), 'VALUE', 'AREA')
        cutArea = ZonalGeometry(Con(IsNull(DEMwBestCuts) == 0, 1), 'VALUE', 'AREA')
        areaRatio = cutArea.maximum / pfArea.maximum
        try:
            assert(areaRatio >= 0.9),'VoidFixed DEM failed area data check'
            DEMwBestCuts.save(output_dem)
            log.info("COMPLETION! Succesful save of cut DEM at " + time.asctime())
            # paraDict = {
            #     '\n\nACPF: DEM Generation and Pit Fill Tool     ' : '\nRun Date: %s' % nowYmd,
            #     # '\nUnknown Vintage Lidar Data: ' : False,#tiles_t_or_f,
            #     '\nEarliest 3DEP Lidar Data: ' : collect_starts_min,
            #     '\nLatest 3DEP Lidar Data: ' : collect_ends_max,
            #     '\nLatest 3DEP Lidar Data: ' : collect_majority
            #     }
           # clib_metadata_template = df.getMetadata(['clib', 'deriv'], proc_dir)
            # addMetadata(output_dem, paraDict, clib_metadata_template, log)
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
	"C:/DEP/Scripts/basics/cmd_cutter_DEM.pyt",
	"//EL3354-02/M$/DEP_bkg_search_newtest/LiDAR_Current/elev_PLib_mean18/07080105/ep3m070801050901.tif",
	"//EL3354-02/D$/DEP/Basedata_Summaries/Basedata_26915.gdb/roads_merge",
	"D:/DEP_Proc_bkg_search_newtest/DEMProc/Cut_dem2013_3m_070801050901/search_070801050901_mean18.pkl",
	"D:/DEP_Proc_bkg_search_newtest/DEMProc/Cut_dem2013_3m_070801050901/scratch.gdb/dfs2cut_070801050901",
	"D:/DEP_Proc_bkg_search_newtest/DEMProc/Cut_dem2013_3m_070801050901/scratch.gdb/good_int_dslv_all",
	"D:/DEP_Proc_bkg_search_newtest/DEMProc/Cut_dem2013_3m_070801050901/scratch.gdb/good_up_pts_all",
	"D:/DEP_Proc_bkg_search_newtest/DEMProc/Cut_dem2013_3m_070801050901/scratch.gdb/good_dn_pts_all",
	"//EL3354-02/M$/DEP_bkg_search_newtest/LiDAR_Current/elev_CLib_mean18/07080105/ec3m070801050901.tif",
	"//EL3354-02/D$/DEP_bkg_search_newtest/Man_Data_ACPF/dep_ACPF2022/07080105/idepACPF070801050901.gdb/cuts_prelim_mean18_dem2013_3m_070801050901",
	"//EL3354-02/D$/DEP_bkg_search_newtest/Man_Data_ACPF/dep_ACPF2022/07080105/idepACPF070801050901.gdb/cuts_final_mean18_dem2013_3m_070801050901",
	"//EL3354-02/D$/DEP_bkg_search_newtest/Man_Data_ACPF/dep_ACPF2022/07080105/idepACPF070801050901.gdb/dprsns2cut_mean18_dem2013_3m_070801050901",
	"D:/DEP_Proc_bkg_search_newtest/DEMProc/Cut_dem2013_3m_070801050901",
	"18.0"]

        for i in parameters[2:]:
            sys.argv.append(i)
    else:
        arcpy.AddMessage("Whoo, hoo! Command-line enabled!")
        # clean up the folder after done processing
        cleanup = True

    input_dem, huc_roads, dfs_2_cut_fc, good_dslv_fc, good_up_dslv_fc, good_dn_dslv_fc, search_distance_file, output_dem, good_cuts_fc, best_cuts_fc, depressions2cut_fc, proc_dir, match_depth = [i for i in sys.argv[1:]]
    messages = msgStub()

    doCutter(input_dem, huc_roads, dfs_2_cut_fc, good_dslv_fc, good_up_dslv_fc, good_dn_dslv_fc, search_distance_file, output_dem, good_cuts_fc, best_cuts_fc, depressions2cut_fc, proc_dir, match_depth, cleanup, messages)
    arcpy.AddMessage("Back from doing!")
