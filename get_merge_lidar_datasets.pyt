# -*- coding: utf-8 -*-
'''A Python program to get feature classes describing lidar datasets available
via the Entwine Point cloud format. This intersects the USGS WESM data on exact project
boundaries with the EPT JSON that shows generalized boundaries and has web addresses
for the data. This dataset is then used in later programs to figure out which EPT
datasets to request.'''
import arcpy
# coding: utf-8

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
        self.label = "EPT_WESM_download"
        self.description = "Creates a feature class to enable EPT downloads"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""

        output_ept_wesm_file = arcpy.Parameter(
            name="ept_wesm_features",
            displayName="Output EPT WESM Feature",
            datatype="DEFeatureClass",
            parameterType='Required',
            direction="Output")
        params = [output_ept_wesm_file]#None
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
        doEPT(parameters[0].valueAsText, cleanup, messages)
        return

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return

def download_file(url, local_filename, log):
    """downloading a file fast using requests"""
    ##https://stackoverflow.com/questions/16694907/download-large-file-in-python-with-requests
    import requests
    import shutil
    log.info('downloading ' + url)
    with requests.get(url, stream=True) as r:
        with open(local_filename, 'wb') as f:
            shutil.copyfileobj(r.raw, f)

def doEPT(ept_wesm_file, cleanup, messages):
    import datetime
    import os
    import sys
    import platform
    import traceback
    import time
    import urllib.request

    import dem_functions as df
    from os.path import join as opj

    # if True:
    try:
        arguments = [ept_wesm_file, cleanup]

        for a in arguments:
            if a == arguments[0]:
                arg_str = str(a) + '\n'
            else:
                arg_str += str(a) + '\n'

        messages.addMessage("Tool: Executing with parameters:\n" + arg_str)

        arcpy.env.overwriteOutput = True


        huc12 = 'XXXXXXXXXXXX'
        if cleanup:
            # log to file only
            log, nowYmd, logName, startTime = df.setupLoggingNoCh(platform.node(), sys.argv[0], huc12)
        else:
            # log to file and console
            log, nowYmd, logName, startTime = df.setupLoggingNew(platform.node(), sys.argv[0], huc12)

        log.info("Beginning execution: " + time.asctime())
        log.info("Log file at " + logName)
        messages.addMessage("Log file at " + logName)

        log.info("Tool: Executing with parameters:\n" + arg_str)

        #create names of outputs so we can see test if it's been run recently
        assert ept_wesm_file.find('.gdb') != -1, "Output must be located in '.gdb' (File Geodatabase)"

        eptDir = os.path.dirname(os.path.dirname(ept_wesm_file))
        log.debug(f"checking for dir: {eptDir}")

        if not os.path.isdir(eptDir):
            log.debug(f"making dir: {eptDir}")
            os.makedirs(eptDir)

        #get geoJSON from https://raw.githubusercontent.com/hobuinc/usgs-lidar/master/boundaries/resources.geojson
        now_ymd_string = nowYmd[:10]#ept_wesm_file[-10:]#datetime.datetime.today().replace(day=1)
        ept_first_of_month_name = "ept_resources_" + now_ymd_string
        ept_4269_first_of_month_name = "ept_resources_epsg4269_" + now_ymd_string
        wesm_first_of_month_name = "main_wesm_" + now_ymd_string
        requested_gdb = os.path.basename(os.path.dirname(ept_wesm_file))
        ept_gdb_path = opj(eptDir, requested_gdb)#'ept.gdb')

        if not arcpy.Exists(ept_gdb_path):
            log.debug(f"making ept gdb: {ept_gdb_path}")
            ept_gdb = arcpy.CreateFileGDB_management(os.path.dirname(ept_gdb_path), os.path.basename(ept_gdb_path))

        arcpy.env.workspace = ept_gdb_path#'in_memory'
        ept_features_path = opj(ept_gdb_path, ept_first_of_month_name)
        ept_4269_features_path = opj(ept_gdb_path, ept_4269_first_of_month_name)
        if not arcpy.Exists(ept_features_path) or not arcpy.Exists(ept_wesm_file):
            # pass # get the ept file
            ept_download_location = opj(eptDir, ept_first_of_month_name + '.geojson')
            if not os.path.exists(ept_download_location):
                log.info(f'downloading from: https://raw.githubusercontent.com/hobuinc/usgs-lidar/master/boundaries/resources.geojson')
                log.info('requesting ept to ' + ept_download_location)
                ept_response = urllib.request.urlretrieve('https://raw.githubusercontent.com/hobuinc/usgs-lidar/master/boundaries/resources.geojson', ept_download_location)
                # download_file('https://raw.githubusercontent.com/hobuinc/usgs-lidar/master/boundaries/resources.geojson', ept_download_location, log)
            # requests.request()
            wesm_download_location = opj(eptDir, wesm_first_of_month_name + '.gpkg')
            if not os.path.exists(wesm_download_location):
                log.info('requesting wesm to ' + wesm_download_location)
                download_file('https://rockyweb.usgs.gov/vdelivery/Datasets/Staged/Elevation/metadata/WESM.gpkg', wesm_download_location, log)
##                wesm_response = urllib.request.urlretrieve('https://rockyweb.usgs.gov/vdelivery/Datasets/Staged/Elevation/metadata/WESM.gpkg', wesm_download_location)
                assert os.path.getsize(wesm_download_location) > 1000000000, "Check WESM download address, should be larger than 1 GB"

            log.info('projecting WESM to EPSG 4269 (NAD83)')
            # do WESM first so map will be in epsg 4269 (NAD83)
            main_wesm_copy = arcpy.conversion.FeatureClassToFeatureClass(opj(wesm_download_location, 'main.wesm'), ept_gdb_path, 'wesm_from_gpkg')
            workunit_lower_field = 'workunit_lower'
            if not workunit_lower_field in df.getfields(main_wesm_copy):
                arcpy.AddField_management(main_wesm_copy, workunit_lower_field, 'TEXT')
            arcpy.CalculateField_management(main_wesm_copy, workunit_lower_field, '!workunit!.lower()', 'PYTHON3')

            log.info('figuring out name conversion between WESM and EPT')
            # if not arcpy.Exists(ept_features_path):
            arcpy.env.outputZFlag = "Disabled"
            ept_features_4326 = arcpy.conversion.JSONToFeatures(ept_download_location, ept_features_path, "POLYGON")
            ept_features = arcpy.management.Project(ept_features_4326, ept_4269_features_path, 4269)
            # else:
            #     ept_features = ept_features_path
            opr_year_field = 'opr_year'
            if not opr_year_field in df.getfields(ept_features):
                arcpy.AddField_management(ept_features, opr_year_field, 'TEXT')
            name_lower_field = 'name_lower'
            if not name_lower_field in df.getfields(ept_features):
                arcpy.AddField_management(ept_features, name_lower_field, 'TEXT')
            alt_name_field = 'alt_name'
            if not alt_name_field in df.getfields(ept_features):
                arcpy.AddField_management(ept_features, alt_name_field, 'TEXT')
            alt2_name_field = 'alt2_name'
            if not alt2_name_field in df.getfields(ept_features):
                rslt3 = arcpy.AddField_management(ept_features, alt2_name_field, 'TEXT')
                log.info(rslt3.getMessages())
            with arcpy.da.UpdateCursor(ept_features, [opr_year_field, 'url', name_lower_field, 'name', alt_name_field, alt2_name_field]) as ucur:
            # with arcpy.da.SearchCursor(ept_features, [opr_year_field, 'url', name_lower_field, 'name', alt_name_field, alt2_name_field], "name = 'USGS_LPC_IN_ET_B5_Allen_2012__LAS_2016' or name = 'USGS_LPC__IL_District7_Lawrence_2014_LAS_2017'") as ucur:
                for urow in ucur:
                    url_year = os.path.basename(os.path.dirname(urow[1]))[-4:]
                    # print(urow)
                    urow[0] = url_year
                    # code a few special ones
                    if urow[3] == 'IA_FullState':
                        urow[2] = 'IA_STATEWIDE_2008'.lower()
                    elif urow[3] == 'MN_FullState':
                        # just chose one of the IWI_REDRIVER collections to represent the state...
                        urow[2] = '	IWI_REDRIVER_I_2009'.lower()
                    else:
                        urow[2] = urow[3].lower().replace('-', '_')
                        # new = urow[3].lower().replace('-', '_')
                        if urow[3].startswith('USGS_LPC'):
                            # some names needed leading or trailing underscores stripped after this step to facilitate match with wesm
                            urow[4] = '_'.join(urow[2].split('_')[2:-2]).strip('_')
                            urow[5] = '_'.join(urow[2].split('_')[2:]).strip('_')
                            # new1 = '_'.join(new.split('_')[2:-2]).strip('_')
                            # new2 = '_'.join(new.split('_')[2:]).strip('_')
                        # print(f"new1: {new1}, new2: {new2}")
                    ucur.updateRow(urow)
            # #             #get WESM from https://rockyweb.usgs.gov/vdelivery/Datasets/Staged/Elevation/metadata/WESM.gpkg


            rslt1 = arcpy.management.JoinField(main_wesm_copy, workunit_lower_field, ept_features, name_lower_field, "id;count;url;opr_year")
            log.info(rslt1.getMessages())
            rslt2 = arcpy.management.JoinField(main_wesm_copy, workunit_lower_field, ept_features, alt_name_field, "id;count;url;opr_year")
            log.info(rslt2.getMessages())
            arcpy.management.JoinField(main_wesm_copy, workunit_lower_field, ept_features, alt2_name_field, "id;count;url;opr_year")

            select_not_null4 = arcpy.analysis.Select(main_wesm_copy, 'wesm_ept_valid_' + now_ymd_string, where_clause = "id_12 IS NOT NULL OR id_1 IS NOT NULL OR id IS NOT NULL")
            select_is_null4 = arcpy.analysis.Select(main_wesm_copy, where_clause = "id_12 IS NULL AND id_1 IS NULL AND id IS NULL")
            mn_ept_features = arcpy.analysis.Select(ept_features, 'mn_fullstate', "name = 'MN_FullState'")

            # sj_data = arcpy.analysis.Select(sj, 'sj_data', "url_12_13 IS NOT Null")
            # mn_sj_data = arcpy.analysis.Clip(sj_data, ept_mn_fullstate)
            mn_select_is_null4 = arcpy.analysis.Clip(select_is_null4, mn_ept_features)#sj_data, ept_mn_fullstate)
            mn_sj_data = arcpy.analysis.SpatialJoin(mn_select_is_null4, mn_ept_features, 'sj_null_wesm_ept')

            # synchronize fields between ept and wesm data
            arcpy.management.DeleteField(mn_sj_data, 'Join_Count; TARGET_FID')
            arcpy.management.DeleteField(mn_sj_data, 'Shape_Length_1; Shape_Area_1')
            arcpy.management.AddFields(select_not_null4, [['id_12_13', 'DOUBLE'], ['count_12_13', 'DOUBLE'], ['url_12_13', 'TEXT'], ['opr_year_12_13', 'TEXT'], ['name_lower', 'TEXT'], ['alt_name', 'TEXT'], ['alt2_name', 'TEXT'], ['name', 'TEXT']])

            select_not_null4_copy = arcpy.management.CopyFeatures(select_not_null4, ept_wesm_file)#ept_features_path)#'wesm_copy')

            log.info('figure out Minnesota statewide data vagaries')
            # srow_copy = arcpy.management.CopyFeatures(srow[1], 'srow_copy')
            mn_ept_neg_buffer = arcpy.analysis.Buffer(mn_ept_features, 'mn_neg_buffer', '-667 METERS')
            missing_mn = arcpy.analysis.Erase(mn_ept_neg_buffer, select_not_null4_copy)

            sp_missing = arcpy.MultipartToSinglepart_management(missing_mn)

            missing_mn_big = arcpy.Select_analysis(sp_missing, where_clause= 'SHAPE_AREA > 0.003')
            big_fields = df.getfields(missing_mn_big)

            fields = ['workunit', 'workunit_id', 'project', 'project_id', 'collect_start', 'collect_end', 'ql', 'spec', 'p_method', 'dem_gsd_meters', 'horiz_crs', 'vert_crs', 'geoid', 'lpc_pub_date', 'lpc_category', 'lpc_reason', 'sourcedem_pub_date', 'sourcedem_category', 'sourcedem_reason', 'onemeter_category', 'onemeter_reason', 'seamless_category', 'seamless_reason', 'lpc_link', 'sourcedem_link', 'metadata_link', 'workunit_lower', 'id', 'count', 'url', 'opr_year', 'id_1', 'count_1', 'url_1', 'opr_year_1', 'id_12', 'count_12', 'url_12', 'opr_year_12', 'name', 'id_12_13', 'count_12_13', 'url_12_13', 'opr_year_12_13', 'name_lower', 'alt_name', 'alt2_name']
            icur = arcpy.da.InsertCursor(select_not_null4_copy, ['OID@', 'SHAPE@'] + big_fields[2:-4] + ['workunit_id', 'collect_start'])

            # give the Minnesota data it's own fictitious workunit_id (and maybe date?)
            # with arcpy.da.SearchCursor(mn_sj_data, ['OID@', 'SHAPE@'] + fields, "OBJECTID = 15") as mn_scur:
            with arcpy.da.SearchCursor(missing_mn_big, ['OID@', 'SHAPE@'] + big_fields[2:-4]) as mn_scur:
                for srow in mn_scur:
                    print(srow[0:5])
                    icur.insertRow(list(srow) + [-10000, datetime.datetime(2008, 4, 21, 0, 0)])

            del icur

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

    return


if __name__ == "__main__":
    import sys

    if len(sys.argv) == 1:
        arcpy.AddMessage("Whoo, hoo! Running from Python Window!")
        cleanup = False

        parameters = ["C:/Program Files/ArcGIS/Pro/bin/Python/envs/arcgispro-py3/pythonw.exe",
    "C:/DEP/Scripts/basics/cmd_ept_wesm_processing.py",
    "C:/DEP/Elev_Base_Data/ept/ept.gdb/ept_resources_2023_05_22"]

        for i in parameters[2:]:
            sys.argv.append(i)
    else:
        arcpy.AddMessage("Whoo, hoo! Command-line enabled!")
        # clean up the folder after done processing
        cleanup = True

    # ept_wesm_file = "C:/DEP/Elev_Base_Data/ept/ept.gdb/ept_resources_2023_05_20"
    ept_wesm_file = sys.argv[1]
    doEPT(ept_wesm_file, cleanup, msgStub())

    # arcpy.AddMessage("Back from doEPT!")