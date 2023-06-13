# -*- coding: utf-8 -*-

import arcpy

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
        self.label = "Tool"
        self.description = ""
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""

        param0 = arcpy.Parameter(
            name="ept_wesm_features",
            displayName="Buffered HUC12 Feature",
            datatype="DEFeatureClass",
            parameterType='Required',
            direction="Input")
        
        param1 = arcpy.Parameter(
            name="snap_raster",
            displayName="Snap Raster",
            datatype="DEFeatureClass",
            parameterType='Required',
            direction="Input")
        
        param2 = arcpy.Parameter(
            name="snap_raster",
            displayName="EPT WESM Merged Features",
            datatype="DEFeatureClass",
            parameterType='Required',
            direction="Input")
        
        param3 = arcpy.Parameter(
            displayName="DEM metadata template",
            datatype="GPDataFile",
            parameterType='Required',
            direction="Input")
        
        param4 = arcpy.Parameter(
            displayName="LiDAR derivatives metadata template",
            datatype="GPDataFile",
            parameterType='Required',
            direction="Input")
        
        param5 = arcpy.Parameter(
            displayName="Local Processing Directory",
            datatype="DEFolder",
            parameterType='Optional',
            direction="Input")
        
        param6 = arcpy.Parameter(
            displayName="Elevation Data Storage Directory",
            datatype="DEFolder",
            parameterType='Optional',
            direction="Input")
        
        param7 = arcpy.Parameter(
            displayName="7za and LASTools Software Directory",
            datatype="DEFolder",
            parameterType='Optional',
            direction="Input")
        
        param8 = arcpy.Parameter(
            displayName="PDAL.exe Location",
            datatype="GPDataFile",
            parameterType='Required',
            direction="Input")
        
        param9 = arcpy.Parameter(
            displayName="Integer Resolution/Ground Sample Distance of output rasters, multiples joined by comma",
            datatype="GPString",
            parameterType='Required',
            direction="Input")
        
        param10 = arcpy.Parameter(
            displayName="Output Pit-Filled Elevation Model",
            datatype="DERasterDataset",
            parameterType='Required',
            direction="Output")
        
        param11 = arcpy.Parameter(
            displayName="Output Bare Earth Minimum Elevation Model",
            datatype="DERasterDataset",
            parameterType='Required',
            direction="Output")
        
        param12 = arcpy.Parameter(
            displayName="Output First Return Maximum Elevation/Canopy Height Model",
            datatype="DERasterDataset",
            parameterType='Required',
            direction="Output")
        
        param13 = arcpy.Parameter(
            displayName="Output Bare Earth Return Count Raster",
            datatype="DERasterDataset",
            parameterType='Required',
            direction="Output")
        
        param14 = arcpy.Parameter(
            displayName="Output First Return Count Raster",
            datatype="DERasterDataset",
            parameterType='Required',
            direction="Output")
        
        param15 = arcpy.Parameter(
            displayName="Output Intensity First Return Minimum Raster",
            datatype="DERasterDataset",
            parameterType='Required',
            direction="Output")
        
        param16 = arcpy.Parameter(
            displayName="Output Intensity First Return Maximum Raster",
            datatype="DERasterDataset",
            parameterType='Required',
            direction="Output")
        
        param17 = arcpy.Parameter(
            displayName="Output Intensity Bare Earth Maximum Raster",
            datatype="DERasterDataset",
            parameterType='Required',
            direction="Output")
        
        param18 = arcpy.Parameter(
            displayName="Output HUC12 Merged Breakline Polygon Features",
            datatype="DEFeatureClass",
            parameterType='Required',
            direction="Output")
        
        param19 = arcpy.Parameter(
            displayName="Output HUC12 Merged Breakline Polyline Features",
            datatype="DEFeatureClass",
            parameterType='Required',
            direction="Output")
        
        param20 = arcpy.Parameter(
            name="ept_wesm_features",
            displayName="EPT WESM Feature for AOI",
            datatype="DEFeatureClass",
            parameterType='Required',
            direction="Output")
        
        param9.values = "3,2,1"#default value to create 3, 2, and 1 meter rasters
        
        params = [param0, param1, param2, param3,
                  param4, param5, param6, param7,
                  param8, param9, param10, param11,
                  param12, param13, param14, param15,
                  param16, param17, param18, param19,
                  param20]
        # params = [huc12_buf_fc, snap, ept_wesm_file, flib_metadata_template, derivative_metadata,
        #  procDir, eleBaseDir, softwareDir, pdal_exe, gsds,
        #  fElevFile, bareEarthReturnMinFile, firstReturnMaxFile, cntFile,cnt1rFile,
        #  int1rMinFile, int1rMaxFile, intBeMaxFile, breakpolys, breaklines, wesm_project_file]
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
                                                                     