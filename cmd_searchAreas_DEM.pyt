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
# changed mergedMdnsFC to mergedMdnsHuc8FC - 2021.02.02 bkgelder
# fixed copying error when no mergedMdns were found
## 2022.06.09 - fixed some railroad processing outputs so they saved correctly to gdb using os.path.join(), previously concatenating string - bkgelder
# 2024.02.19 - moved code to Python 3

# Import system modules
import arcpy
import sys
import os
import traceback
import time
import platform
from math import sqrt, atan2, pi
from arcpy.sa import *
import dem_functions as df
from os.path import join as opj
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
        self.label = "Hole_Puncher"
        self.description = "Punches holes in a DEM and fills remaining to remove depressions shallower than a criteria"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""

        param0 = arcpy.Parameter(
            displayName="Input Elevation Model",
            datatype="DERasterDataset",
            parameterType='Required',
            direction="Input")
        
        param1 = arcpy.Parameter(
            displayName="Output Punched Elevation Model",
            datatype="DERasterDataset",
            parameterType='Required',
            direction="Output")
        
        param2 = arcpy.Parameter(
            displayName="Punched DEM metadata template",
            datatype="GPDataFile",
            parameterType='Required',
            direction="Input")
        
        param3 = arcpy.Parameter(
            name="depressions_fc",
            displayName="Punched depressions feature class",
            datatype="DEFeatureClass",
            parameterType='Required',
            direction="Output")
        
        param4 = arcpy.Parameter(
            displayName="Depression Punch Depth Threshold",
            datatype="GPString",
            parameterType='Required',
            direction="Input")
        
        param5 = arcpy.Parameter(
            displayName="Depression Punch Area Threshold",
            datatype="GPString",
            parameterType='Optional',
            direction="Input")
        
        param6 = arcpy.Parameter(
            displayName="Depression Punch Area Threshold",
            datatype="GPString",
            parameterType='Optional',
            direction="Input")
        
        param7 = arcpy.Parameter(
            name = "procDir",
            displayName="Local Processing Directory",
            datatype="DEFolder",
            parameterType='Optional',
            direction="Input")
        
        parameters = [param0, param1, param2, param3, param4, param5, param6, param7]
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
        doSearcher(params[0].valueAsText, params[1].valueAsText, params[2].valueAsText, params[3].valueAsText, params[4].valueAsText, params[5].valueAsText, params[6].valueAsText, params[7].valueAsText, cleanup, messages)
        return

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return



def doSearcher(input_dem, output_dem, plib_metadata, depressions_fc, depth_threshold, area_threshold, procDir, cleanup, messages):

    try:
        arguments = [input_dem, output_dem, plib_metadata, depressions_fc, depth_threshold, area_threshold, procDir, cleanup, messages]

        for a in arguments:
            if a == arguments[0]:
                arg_str = a + '\n'
            else:
                arg_str += a + '\n'

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

        log.info("Tool: Executing with parameters:\n" + arg_str)

        arcpy.CheckOutExtension("Spatial")
        arcpy.env.overwriteOutput = True

        arcpy.env.snapRaster = input_dem
        arcpy.env.cellSize = input_dem
        ProcSize = int(arcpy.Raster(input_dem).meanCellHeight)

        startTime = time.time()
        log.info("Beginning logging for script at " + str(time.asctime()))

        log.info('sys.argv is: ' + str(sys.argv) + '\n')

        log.info('log file is ' + logName)

    ####------------------------------------------------------------------------------

        if not os.path.isdir(procDir):
            os.makedirs(procDir)
        arcpy.env.scratchWorkspace = procDir

    ## Set the environments
    ##            ## If you set a scratch workspace first you can control where the scratchGDB or scratchFolder are created
    ##            ## otherwise it defaults to a user's temp folder
    ##            ## if you don't set anything it will go to 'in_memory'
    ##
    ##    sfldr = arcpy.env.scratchFolder
        gdb = arcpy.env.scratchGDB
        arcpy.env.scratchWorkspace = gdb
        arcpy.env.workspace = gdb

        inm = 'in_memory\\'

        huc8gdb = os.path.dirname(huc8RoadsFC)
        if not arcpy.Exists(huc8gdb):
            if not os.path.isdir(os.path.dirname(huc8gdb)):
                os.makedirs(os.path.dirname(huc8gdb))
            huc8gdbResult = arcpy.CreateFileGDB_management(os.path.dirname(huc8gdb), os.path.basename(huc8gdb))

    ####------------------------------------------------------------------------------
    ####    ## Create a layer of WBD boundary to buffer and clip datasets with
        huc8Lyr = arcpy.MakeFeatureLayer_management(huc8fc, 'HUC12FCLayer', '"HUC8" = \'' + huc8 + "'")
        # bnd = arcpy.CopyFeatures_management(huc8Lyr, opj(gdb, "bnd_" + huc8))
        bndBuffer = arcpy.Buffer_analysis(huc8Lyr, opj(gdb, 'buf_' + huc8 + '_1km'), '1000 METER')

        log.debug("Processing road clip and buffer distance at " + time.asctime())

        try:
            mergedMdnList = []

        ## Calculate additional search distance needed to cross roads, railroads, runways 
            hucRoadsIncService = arcpy.Clip_analysis(roadsFC, bndBuffer, inm + "roads_all")
            if 'fclass' in df.getfields(hucRoadsIncService):
                hucRoadsWhole = arcpy.Select_analysis(hucRoadsIncService, inm + 'roads', 'fclass <> \'service\' AND fclass <> \'path\'AND fclass <> \'cycleway\'AND fclass <> \'footway\'')
            else:
                hucRoadsWhole = arcpy.Select_analysis(hucRoadsIncService, inm + 'roads', 'type <> \'service\' AND type <> \'path\'')

            hucRoads = arcpy.FeatureToLine_management(hucRoadsWhole, inm + 'roads_broken')
            hucRoadsCopy = arcpy.CopyFeatures_management(hucRoads, huc8RoadsFC)

            if 'fclass' in df.getfields(hucRoads):
                # mwRoadsOld = arcpy.Select_analysis(hucRoads, opj(gdb, 'huc_roads_old'), 'oneway = \'F\' AND fclass NOT LIKE \'%_link\'')#<> \'motorway_link\'')
                mwRoads = arcpy.Select_analysis(hucRoads, opj(gdb, 'huc_roads_mw'), 'oneway <> \'B\' AND (fclass = \'motorway\' OR fclass = \'primary\' OR fclass = \'trunk\')')
            else:
                mwRoads = arcpy.Select_analysis(hucRoads, opj(gdb, 'huc_roads_mw'), '"oneway" = 1 AND "TYPE" NOT LIKE \'%_link\'')
            descMwRoads = arcpy.Describe(mwRoads)

            mwRoadsGNT = arcpy.GenerateNearTable_analysis(mwRoads, mwRoads, 'mw_roads_gnt', '50 METERS', closest = 'ALL', closest_count = 10, location = 'LOCATION', angle = 'ANGLE')

            mwRoadsExtraSearch = arcpy.TableSelect_analysis(mwRoadsGNT, 'mw_roads_xtra_search', '"NEAR_DIST" > 1')
            arcpy.JoinField_management(mwRoads, descMwRoads.OIDFieldName, mwRoadsExtraSearch, 'IN_FID', ['NEAR_DIST', 'FROM_X', 'FROM_Y', 'NEAR_X', 'NEAR_Y'])#osm_id')
            mwRoadsMdn = arcpy.Select_analysis(mwRoads, inm + 'mw_roads_mdn', 'NEAR_DIST > 0')

            if int(arcpy.GetCount_management(mwRoadsMdn).getOutput(0)) > 0:

                rdMedianNearLine = arcpy.XYToLine_management(mwRoadsMdn, inm + 'mdn_roads_near_lines', 'FROM_X', 'FROM_Y', 'NEAR_X', 'NEAR_Y', spatial_reference = mwRoads)
                rdMedianNearLineCentroid = arcpy.CreateFeatureclass_management(inm[:-1], 'rd_mdn_poly_centroid', 'POINT', spatial_reference = mwRoads)
                df.tryAddField(rdMedianNearLineCentroid, 'NEAR_DIST', 'DOUBLE')
                iCurCentroid = arcpy.da.InsertCursor(rdMedianNearLineCentroid, ['SHAPE@', 'NEAR_DIST'])
                with arcpy.da.SearchCursor(rdMedianNearLine, ['SHAPE@XY', 'SHAPE@LENGTH']) as sCurCentroid:
                    for sRowCentroid in sCurCentroid:
                        if sRowCentroid[0] != None:
                            iCurCentroid.insertRow([arcpy.Point(sRowCentroid[0][0],sRowCentroid[0][1]), sRowCentroid[1]])
                del iCurCentroid

                ftpRoadsBuf = arcpy.FeatureToPolygon_management([hucRoads, bndBuffer], inm + 'rd_buf_ftp')
                ftpRoadsBufLayer = arcpy.MakeFeatureLayer_management(ftpRoadsBuf, 'ftp_rd_buf_layer')
                # rdParallelNear = arcpy.SelectLayerByLocation_management(ftpRoadsBufLayer, 'CONTAINS_CLEMENTINI', rdMedianNearLineCentroid)

                rdMdnPrelim = arcpy.CopyFeatures_management(ftpRoadsBufLayer, inm + 'rd_mdn_polys_prelim')

    ##            ## filter out the polygons not really close to the raods with medians
                df.tryAddField(rdMdnPrelim, 'ID_RD_MDN', 'LONG')
                with arcpy.da.UpdateCursor(rdMdnPrelim, ['ID_RD_MDN', 'OID@']) as ucur:
                    for urow in ucur:
                        urow[0] = urow[1]
                        ucur.updateRow(urow)

                mwRoadsBuffer = arcpy.Buffer_analysis(mwRoads, buffer_distance_or_field = '50 METERS')
                rdMdnClip = arcpy.Clip_analysis(rdMdnPrelim, mwRoadsBuffer)

                valueDict = {r[0]:(r[1:]) for r in arcpy.da.SearchCursor(rdMdnClip, ['ID_RD_MDN', 'SHAPE@AREA'])}
                df.tryAddField(rdMdnPrelim, 'AREA_CLIP', 'DOUBLE')
                df.tryAddField(rdMdnPrelim, 'AREA_TRUE', 'DOUBLE')
                fields_to_join = ['AREA_CLIP']#'ID_RD_MDN', 
                with arcpy.da.UpdateCursor(rdMdnPrelim, ['ID_RD_MDN'] + fields_to_join + ['AREA_TRUE', 'SHAPE@AREA']) as ucur:
                    for urow in ucur:  
                        # store the Join value of the row being updated in a keyValue variable  
                        keyValue = urow[0]  
                        # verify that the keyValue is in the Dictionary  
                        if keyValue in valueDict:
                            for i, field in enumerate(fields_to_join):
                                urow[i+1] = valueDict[keyValue][i]
                            urow[-2] = urow[-1]
                            ucur.updateRow(urow)
                del valueDict
                rdMdn = arcpy.Select_analysis(rdMdnPrelim, where_clause = 'AREA_CLIP/AREA_TRUE > 0.99')
                mergedMdnList.append(rdMdn)

    ## ---------------------------------------------------------------------
        ## Find parallel road/railroad sections
    ##      Find railroads that are near and parallel to roads and define area between as a median
        ## Buffer the railroads so they become a polygon that we can intersect with road and RR search buffers and determine if points cross one of these features
            hucRailroadsPrelim = arcpy.Clip_analysis(rrsFC, bndBuffer, inm + "rrs_prelim")
            ## Start railroad processing
            if df.testForZero(hucRailroadsPrelim):
                types = ['abandoned', 'light_rail', 'rail', 'preserved', 'yard']
                if 'fclass' in df.getfields(hucRailroadsPrelim):
                    sel = df.buildStringSelection(types, 'fclass')
                else:# assume 'type' field in OSM data
                    sel = df.buildStringSelection(types, 'type')
                hucRailroads = arcpy.Select_analysis(hucRailroadsPrelim, inm + 'rrs', sel)
                df.copyfc(verbose, hucRailroads, gdb)

                if df.testForZero(hucRailroads):
                    ftpRoadsRRsBuf = arcpy.FeatureToPolygon_management([hucRoads, hucRailroads, bndBuffer], opj(gdb, 'rd_rr_buf_ftp'))
                    ## Slice railroads into short features so we can compare railroad and road bearings
                    hucRailDissolve = arcpy.CreateFeatureclass_management(gdb, 'hucRail10', template = hucRailroads, spatial_reference = hucRailroads)
                    initFields = df.getfields(hucRailDissolve)[2:-1]
                    descRRLong = arcpy.Describe(hucRailDissolve)

                    iCur = arcpy.da.InsertCursor(hucRailDissolve, [u'Shape@'] + initFields)
                    splitEvery = 250 #METERS
                    with arcpy.da.SearchCursor(hucRailroads, [u'Shape@'] + initFields) as lines:
                        for line in lines:
                            if line[0].length > splitEvery:
                                out_count = int(line[0].length/splitEvery)
                                for i in range(0, out_count):
                                    part = line[0].segmentAlongLine(i/float(out_count), ((i+1)/float(out_count)), True)
                                    partList = [part]
                                    for j in line[1:]:
                                        partList.append(j)
                                    iCur.insertRow(partList)
                            else:
                                iCur.insertRow(line)
                    del iCur, line, lines


                    df.tryAddField(hucRailDissolve, "HUC_FID", "LONG")
                    arcpy.CalculateField_management(hucRailDissolve, "HUC_FID", '!OBJECTID!', "PYTHON")

                    df.tryAddField(hucRailDissolve, 'BEARING', 'DOUBLE')
                    with arcpy.da.UpdateCursor(hucRailDissolve, ['OID@', 'SHAPE@', 'BEARING']) as ucur:
                        for urow in ucur:
                            ## Calculate the distance (hypotenuse) and angle between previous and current point
                            initPnt = urow[1].firstPoint
                            endPnt = urow[1].lastPoint
                            pnt_dx = endPnt.X - initPnt.X
                            pnt_dy = endPnt.Y - initPnt.Y
                            bearingAngle = atan2(pnt_dx, pnt_dy)*(360.0/(2*pi))

                            urow[2] = bearingAngle

                            ucur.updateRow(urow)

    ##                            First find roads that are in this buffer area
                    rrRdBufferDist = '40 METERS'
                    hucRRsBuffer = arcpy.Buffer_analysis(hucRailDissolve, inm + 'rrs_bfr', rrRdBufferDist, line_end_type = 'FLAT')

                    roadsInRRsBuffer = arcpy.Clip_analysis(hucRoads, hucRRsBuffer, inm + 'rds_in_rr_bfr')
                    df.copyfc(verbose, roadsInRRsBuffer, gdb)
                    df.tryAddField(roadsInRRsBuffer, "HUC_RD_ID", "LONG")
                    arcpy.CalculateField_management(roadsInRRsBuffer, "HUC_RD_ID", '!OBJECTID!', "PYTHON")

                ## Calculate bearing of road section
                    df.tryAddField(roadsInRRsBuffer, 'BEARING', 'DOUBLE')
                    with arcpy.da.UpdateCursor(roadsInRRsBuffer, ['OID@', 'SHAPE@', 'BEARING']) as ucur:
                        for urow in ucur:
                            ## Calculate the distance (hypotenuse) and angle between previous and current point
                            initPnt = urow[1].firstPoint
                            endPnt = urow[1].lastPoint
                            pnt_dx = endPnt.X - initPnt.X
                            pnt_dy = endPnt.Y - initPnt.Y
                            bearingAngle = atan2(pnt_dx, pnt_dy)*(360.0/(2*pi))
                            urow[2] = bearingAngle
                            ucur.updateRow(urow)

                    df.copyfc(verbose, roadsInRRsBuffer, gdb)

                    rdsRRsGNT = arcpy.GenerateNearTable_analysis(hucRailDissolve, roadsInRRsBuffer, inm + 'rds_in_rrs_gnt', rrRdBufferDist, location = 'LOCATION', angle = "ANGLE", closest = "ALL")
                    df.addCalcJoin(rdsRRsGNT, 'IN_FID', hucRailDissolve, 'HUC_FID', ['RR_BEARING', 'DOUBLE'], '!BEARING!')
                    df.addCalcJoin(rdsRRsGNT, 'NEAR_FID', roadsInRRsBuffer, 'HUC_RD_ID', ['RD_BEARING', 'DOUBLE'], '!BEARING!')

                    df.tryAddField(rdsRRsGNT, 'BEAR_DIF', 'DOUBLE')
                    with arcpy.da.UpdateCursor(rdsRRsGNT, ['BEAR_DIF', 'RR_BEARING', 'RD_BEARING']) as ucur:
                        for urow in ucur:
                            if urow[1] != None and urow[2] != None:
                                urow[0] = df.angleDif(urow[1], urow[2])
                                ucur.updateRow(urow)

                    angleCrit = 20 #Degrees
                    rdsRRsPrll = arcpy.TableSelect_analysis(rdsRRsGNT, opj(gdb, 'rds_rrs_gnt_prll'), '((BEAR_DIF > '+str(0-angleCrit)+' AND BEAR_DIF < ' + str(0+angleCrit)+') OR BEAR_DIF > '+str(180-angleCrit)+ ' OR BEAR_DIF < '+str(-180+angleCrit)+ ')')# AND NEAR_DIST > 0')

                    rdsRRsPrllStats = arcpy.Statistics_analysis(rdsRRsPrll, opj(gdb, 'rds_rrs_prll_stats'), [['BEAR_DIF', 'RANGE']], 'IN_FID')
                    arcpy.JoinField_management(hucRailDissolve, 'OBJECTID', rdsRRsPrllStats, 'IN_FID', 'RANGE_BEAR_DIF')
                    # prllRdsOnBothSidesRrs = arcpy.Select_analysis(hucRailDissolve, opj(gdb, 'prll_rds_both_sides_rrs'), 'RANGE_BEAR_DIF > '+str(180-2*angleCrit))

                ## Now get just the roads that parallel the railroads and do the intersection again
                    arcpy.JoinField_management(roadsInRRsBuffer, 'OBJECTID', rdsRRsPrll, 'NEAR_FID', 'BEAR_DIF')
                    prllIntRoads = arcpy.Select_analysis(roadsInRRsBuffer, opj(gdb, 'rd_int_rr_prll'), 'BEAR_DIF IS NOT NULL')
                    if int(arcpy.GetCount_management(prllIntRoads).getOutput(0)) > 0:
                        prllroadsInRRsBuffer = arcpy.Intersect_analysis([prllIntRoads, hucRRsBuffer], opj(gdb, 'rds_prll_int_rrs'))

                    ## Now create some points that are at centroid of nearest line between road and rr 
                        rrNearLine = arcpy.XYToLine_management(rdsRRsPrll, opj(gdb, 'rr_near_rd_lines'), 'FROM_X', 'FROM_Y', 'NEAR_X', 'NEAR_Y', spatial_reference = mwRoads)
                        arcpy.JoinField_management(roadsInRRsBuffer, 'OBJECTID', rdsRRsPrll, 'NEAR_FID', 'FROM_X')

                        rrRdNearLineCentroid = arcpy.CreateFeatureclass_management(gdb, 'rr_rd_line_centroid', 'POINT', spatial_reference = mwRoads)
                        df.tryAddField(rrRdNearLineCentroid, 'NEAR_DIST', 'DOUBLE')
                        iCurCentroid = arcpy.da.InsertCursor(rrRdNearLineCentroid, ['SHAPE@', 'NEAR_DIST'])
                        with arcpy.da.SearchCursor(rrNearLine, ['SHAPE@XY', 'SHAPE@LENGTH']) as sCurCentroid:
                            for sRowCentroid in sCurCentroid:
                                if sRowCentroid[0] != None:
                                    iCurCentroid.insertRow([arcpy.Point(sRowCentroid[0][0],sRowCentroid[0][1]), sRowCentroid[1]])
                        del iCurCentroid

                    ## Now use clip to whittle everything down to just the area we want to use
                        # something around here is creating a 'scratch.shp' file in the MedianProc folder
                        rrMdnClip = arcpy.Clip_analysis(ftpRoadsRRsBuf, hucRRsBuffer, opj(gdb, 'rr_mdn_clip'))#rrNearLineCentroidBuffer
                        hucRoadsBuffer = arcpy.Buffer_analysis(prllroadsInRRsBuffer, inm + 'rd_line_bfr', rrRdBufferDist)
                        rdRrMdnClip = arcpy.Clip_analysis(rrMdnClip, hucRoadsBuffer, opj(gdb, 'rd_rr_mdn_clip'))
                        rdRrMdnSngl = arcpy.MultipartToSinglepart_management(rdRrMdnClip, opj(gdb, 'rd_rr_mdn_sngl'))

                        ftpLayer2 = arcpy.MakeFeatureLayer_management(rdRrMdnSngl, 'ftp_layer2')#ftpRoadsRRsBuf, 'ftp_layer')

                        # rrParallelNear = arcpy.SelectLayerByLocation_management(ftpLayer2, 'CONTAINS_CLEMENTINI', rrRdNearLineCentroid) 
                        rdRrMdnBig = arcpy.CopyFeatures_management(ftpLayer2, opj(gdb, 'rr_mdn_polys'))
                    ## Clip by railroad buffer then road buffer to limit median buffer area

                        try:
                            hucRoadsBufferFlat = arcpy.Buffer_analysis(prllIntRoads, inm + 'rd_line_bfr', rrRdBufferDist, line_end_type = 'FLAT')
                        except:
                            hucRoadsBufferFlat = arcpy.Buffer_analysis(prllIntRoads, inm + 'rd_line_bfr2', rrRdBufferDist)#, line_end_type = 'FLAT')
    ##                        log.warning('WARNING:flat buffer failure on ' + huc12)
                        rdRrMdnJoinClip = arcpy.Clip_analysis(rdRrMdnBig, hucRoadsBufferFlat, opj(gdb, 'rd_rr_mdn_join_clip'))
                        rdRrMJCftp = arcpy.FeatureToPolygon_management([rdRrMdnJoinClip, prllIntRoads], opj(gdb, 'rd_rr_ftp'))

                        ftpLayer3 = arcpy.MakeFeatureLayer_management(rdRrMJCftp, 'ftp_layer3')#ftpRoadsRRsBuf, 'ftp_layer')

                        ftpSel = arcpy.SelectLayerByLocation_management(ftpLayer3, 'INTERSECT', hucRailDissolve)#prllIntRoads))
                        ftpSel2 = arcpy.SelectLayerByLocation_management(ftpSel, 'INTERSECT', prllIntRoads, selection_type = 'SUBSET_SELECTION')
                        rdRrMdnFinal = arcpy.CopyFeatures_management(ftpSel2, opj(gdb, 'rr_mdn_polys3'))
                        mergedMdnList.append(rdRrMdnFinal)#MdnSeparate)

    ## ------------------------------------------------------------------------------

                ## Find railyards and other double tracks
                if df.testForZero(hucRailroads):#int(arcpy.GetCount_management(hucRailroadsPrelim).getOutput(0)) > 0:
                    rrSearchDist = 100 #METERS
                    slPts, slLines, slLinesL, slLinesR, slLinesSin, slPtsStart = df.stationLines10(gdb, hucRailDissolve, 'rr_sl_pts', 'rr_sl_lines', rrSearchDist/2.0, rrSearchDist, hucRailroads, 'rr_sl_pts_start')

                    rrSlInt = arcpy.Intersect_analysis([slLines, hucRailDissolve], inm + 'rr_sl_int', output_type = "POINT")
                    df.addCalcJoin(rrSlInt, 'arcid', hucRailDissolve, 'HUC_FID', ['sl_BEARING', 'DOUBLE'], '!BEARING!')
                    df.tryAddField(rrSlInt, 'BEAR_DIF', 'DOUBLE')
                    with arcpy.da.UpdateCursor(rrSlInt, ['BEAR_DIF', 'BEARING', 'sl_BEARING']) as ucur:
                        for urow in ucur:
                            if urow[1] != None and urow[2] != None:
                                urow[0] = df.angleDif(urow[1], urow[2])
                                ucur.updateRow(urow)
                    rrSelf = arcpy.Select_analysis(rrSlInt, inm + 'rr_self', 'HUC_FID = arcid')
                    rrNonSelf = arcpy.Select_analysis(rrSlInt, inm + 'rr_non_self_prll', 'HUC_FID <> arcid AND (BEAR_DIF < -170 OR (BEAR_DIF > -10 AND BEAR_DIF < 10) OR BEAR_DIF > 170)')
                    df.tryAddField(rrNonSelf, 'NEAR_DIST', 'DOUBLE')

                    with arcpy.da.SearchCursor(rrSelf, ['SHAPE@XY', 'HUC_FID', 'arcid', 'cs_id'], 'HUC_FID = arcid') as scur:
                        for srow in scur:
                            with arcpy.da.UpdateCursor(rrNonSelf, ['SHAPE@XY', 'HUC_FID', 'arcid', 'cs_id', 'NEAR_DIST'], 'cs_id = ' + str(srow[3])) as ucur:
                                for urow in ucur:
                                    urow[4] = sqrt(pow(srow[0][0]-urow[0][0], 2) + pow(srow[0][1]-urow[0][1], 2))
                                    ucur.updateRow(urow)
                    
                    df.tryAddField(rrNonSelf, 'NEAR_GAP', 'DOUBLE')
                    prev_in_fid = -1
                    with arcpy.da.UpdateCursor(rrNonSelf, ['cs_id', 'NEAR_DIST', 'NEAR_GAP'], sql_clause = (None, 'ORDER BY cs_id, NEAR_DIST')) as ucur:
                        for urow in ucur:
                            if prev_in_fid == urow[0]:
                                urow[2] = urow[1] - prev_near_dist
                                ucur.updateRow(urow)

                            prev_in_fid = urow[0]
                            prev_near_dist = urow[1]
                    rrSlStats = arcpy.Statistics_analysis(rrNonSelf, inm + 'rr_sl_stats', [['NEAR_DIST', 'MIN'], ['NEAR_DIST', 'MEAN'], ['NEAR_DIST', 'MAX'], ['NEAR_DIST', 'COUNT'], ['NEAR_GAP', 'MEAN'], ['NEAR_GAP', 'MAX']], 'cs_id')
                    parallelRRs = arcpy.TableSelect_analysis(rrSlStats, inm + 'rr_gnt_parallel', 'COUNT_NEAR_DIST > 4 AND MEAN_NEAR_GAP < 30')
                    arcpy.JoinField_management(parallelRRs, 'cs_id', slLines, 'cs_id', 'arcid')
                    parallelRRsSummary = arcpy.Statistics_analysis(parallelRRs, inm + 'rr_prll_smry', [['MAX_NEAR_DIST', 'MAX']], 'arcid')
                    arcpy.AlterField_management(parallelRRsSummary, 'MAX_MAX_NEAR_DIST', 'MAX_NEAR_DIST')
                    arcpy.JoinField_management(hucRailDissolve, 'HUC_FID', parallelRRsSummary, 'arcid', 'MAX_NEAR_DIST')
                    hucRRsParallelNear = arcpy.Select_analysis(hucRailDissolve, inm + 'rrs_prll_near', 'MAX_NEAR_DIST > 0')
                                
                    if df.testForZero(hucRRsParallelNear) > 0:
                        rrParallelNearBuffer = arcpy.Buffer_analysis(hucRRsParallelNear, inm + 'rr_prll_bfr', str(rrSearchDist/2.0) + ' METERS')
                        rrPrllNearBufferRaster = arcpy.PolygonToRaster_conversion(rrParallelNearBuffer, 'BUFF_DIST', opj(gdb, 'rr_prl_bfr'), cellsize = ProcSize)
                        rrPrllRG = RegionGroup(Int(rrPrllNearBufferRaster))
                        rrPrllRgPoly = arcpy.RasterToPolygon_conversion(rrPrllRG, opj(gdb, 'rr_prll_rg_poly'), raster_field = 'VALUE')
                        rrPrllInt = arcpy.Intersect_analysis([rrPrllRgPoly, hucRRsParallelNear], inm + 'rr_prll_rg_int')

                        gridfield = 'gridcode'
                        gridfield2 = 'grid_code'
                        rrPrllBg = arcpy.MinimumBoundingGeometry_management(rrPrllInt, inm + 'rr_prll_bg', 'CONVEX_HULL', group_option = 'LIST', group_field = gridfield)
                        mergedMdnList.append(rrPrllBg)

    ## ------------------------------------------------------------------------------
                                
        ## Buffer the runways so they become a polygon that we can intersect with road and RR search buffers and determine if points cross one of these features
            hucRunways = arcpy.Clip_analysis(apFC, bndBuffer, inm + "runway_so")
            if int(arcpy.GetCount_management(hucRunways).getOutput(0)) > 0:
                apSrchExp = '2*!Width!/3.28'
                apBfrFld = 'rw_bfr_dist'#bfrFldList[2]
                df.tryAddField(hucRunways, apBfrFld, "SHORT")
                arcpy.CalculateField_management(hucRunways, apBfrFld, apSrchExp, "PYTHON")
                rwsBfr = arcpy.Buffer_analysis(hucRunways, inm + "rw_srch_bfr", apBfrFld, "FULL", "ROUND", "LIST", apBfrFld)
                mergedMdnList.append(rwsBfr)

    ## ------------------------------------------------------------------------------

        ## Merge all the transportation buffers together
            if len(mergedMdnList) > 0:
                mergedMdns = arcpy.Merge_management(mergedMdnList, inm + 'rd_rr_rd_rw_mrg')
                if df.testForZero(mergedMdns):
            ##                            mergedMdns = arcpy.Merge_management([rdMdnSeparate, rrRdMergedMdnZoned, rrMdnFinal, rwsBfr], gdb + 'rd_rr_rd_rw_mrg')
                    df.tryAddField(mergedMdns, 'T', 'SHORT')
                    arcpy.CalculateField_management(mergedMdns, 'T', '1', 'PYTHON')
                    mergedMdnsRasterT = Raster(arcpy.PolygonToRaster_conversion(mergedMdns, 'T', opj(gdb, 'mrgd_mdn'), cellsize = ProcSize))
                    mergedMdnsRaster = Raster(arcpy.PolygonToRaster_conversion(mergedMdns, arcpy.Describe(mergedMdns).OIDFieldName, opj(gdb, 'mrgd_mdns'), cellsize = ProcSize))
                    mergedMdnsCost = CostDistance(Con(IsNull(mergedMdnsRasterT) == 1, 1), Con(IsNull(mergedMdnsRasterT) == 1, 0, 1))
    ####                mergedMdnsCost.save(cp + 'mdn_cost')

                    mergedMdnsCostZst = ZonalStatisticsAsTable(mergedMdnsRaster, 'VALUE', mergedMdnsCost, inm + 'mrgd_mdn_zst')
                ## populate new search distance field - 2x because cost distance maxes at center, ##1.5x to cut at 60 degree angle, 
                    df.addCalcJoin(mergedMdns, arcpy.Describe(mergedMdns).OIDFieldName, mergedMdnsCostZst, 'VALUE', ['bfr_dist', 'DOUBLE'], '2*!MAX!')

                    mergedMdnsCopy = arcpy.CopyFeatures_management(mergedMdns, mergedMdnsHuc8FC)
                else:
                    mergedMdns = None
    ##                mergedMdnsCost = None

            else:
                mergedMdns = None
    ##            mergedMdnsCost = None
            
            log.debug("Transport processing copying done at " + time.asctime())

        except Exception as e:
    ##        log.debug(e.message)
            arcpy.AddError(e.message)

            # Get the traceback object
            tb = sys.exc_info()[2]
            tbinfo = traceback.format_tb(tb)[0]

            # Concatenate information together concerning the error into a message string
            pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
            print(pymsg)

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
	"C:/DEP/Scripts/basics/cmd_searchAreas_DEM.pyt",
	"C:/DEP/LiDAR_Current/elev_FLib_mean18/07080105/ef3m070801050901.tif",
	"D:/DEP/Basedata_Summaries/elev_PLib_mean18/07080105/ep3m070801050901.tif",
	"D:/DEP/Basedata_Summaries/elev_PLib_mean18/07080105/ep3m070801050901.tif",
	"D:/DEP/Basedata_Summaries/elev_PLib_mean18/07080105/ep3m070801050901.tif",
	"C:/DEP/Man_Data_ACPF/dep_ACPF2022/07080105/idepACPF070801050901.gdb/dprsns_mean18_dem2013_3m_070801050901",
	"5.0",
	"500",
	"C:/DEP_Proc/DEMProc/Cut_dem2013_3m_070801050901"]

        for i in parameters[2:]:
            sys.argv.append(i)
    else:
        arcpy.AddMessage("Whoo, hoo! Command-line enabled!")
        # clean up the folder after done processing
        cleanup = True

    messages = msgStub()

    # inputs
    input_dem = sys.argv[1]
    huc8fc = sys.argv[2]
    roadsFC = sys.argv[3]
    rrsFC = sys.argv[4]
    apFC = sys.argv[5]

    # outputs
    mergedMdnsHuc8FC = sys.argv[6]
    huc8RoadsFC = sys.argv[7]

    # local processing directory
    procDir = sys.argv[8]

    doSearcher(input_dem, huc8fc, roadsFC, rrsFC, apFC, mergedMdnsHuc8FC, huc8RoadsFC, procDir, cleanup, messages)
    arcpy.AddMessage("Back from doing!")

