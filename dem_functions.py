import arcpy, os, sys, traceback, math, datetime, logging
from arcpy.sa import *
import subprocess
import time
from os.path import join as opj

# a compendium of functions to enable DEM processing
# written by Brian Gelder, bkgelder@iastate.edu
# 2019.03.23
# Python 2.7, with an eye towards Python 3


def channelized_areas(meter_dem, proc_dir, pro_crv):
    aspect = Aspect(meter_dem)

# Calcualte range of aspect in radians
    deg2rad = math.pi / 180.0

    aspRad = aspect * deg2rad

    # pihalf = 0.5*math.pi
    # pi3halfs = 1.5*math.pi
    pi2x = 2.0*math.pi

# Use a 3x3 kernel (window) size (no comparison to center cell)
##    window = [(-1,-1),(-1,0),(-1,1),(0,-1),(0,0),(0,1),(1,-1),(1,0),(1,1)]
    window = [(-1,0),(0,-1),(0,1),(1,0)]

    refAsp = aspRad
    refAspAlt = refAsp - pi2x

# hack to calculate aspect difference by comparing by nearest four cells (window) 
    shifts  = []
    for index, cells in enumerate(window):
        aspShift = arcpy.Shift_management(aspRad, opj(proc_dir, 'asp_' + str(index)), cells[0]*aspRad.meanCellHeight, cells[1]*aspRad.meanCellHeight)
        compAsp = Raster(aspShift)

        minCellAspDif = CellStatistics([Abs(refAsp - compAsp), Abs(refAsp - (compAsp - pi2x)), Abs(refAspAlt - compAsp)], 'MINIMUM')
        minCellAspDif.save(opj(proc_dir, 'asp_dif_' + str(window.index(cells))))
        shifts.append(minCellAspDif.name)
    del index
    
    arcpy.env.workspace = proc_dir
    aspRange = CellStatistics(shifts, 'MAXIMUM')
    # aspRange.save(opj(proc_dir, 'aspRng2'))

## Create terrain position index for display
    
    tpi = makeTPI(meter_dem)

    aspLow = Con(tpi < 0, Con(aspRange, 1, 0, '"VALUE" > 1.05'))
    aspLow.save(opj(proc_dir, 'asplow'))
    expLowAsp = Con(IsNull(aspLow), 0, Expand(Con(pro_crv > 1.0, aspLow), 1, 1))
##    expLowAsp.save(cp + 'explowasp')
    expLowAsp1 = Con(expLowAsp == 1, expLowAsp)
    expLowAspTF = Con(IsNull(expLowAsp1), 0, 1)
    expLowAspTF.save(opj(proc_dir, 'crv_tf'))

    return aspect, tpi, aspLow, expLowAspTF
    

def getHucSelection(pdcAcres, pdcSlope, thresh):

    selStr = pdcAcres + ' >= 15 AND ' + pdcSlope + ' > 0 AND (( AGarea.VALUE_1 * 0.000247) / ' + pdcAcres + ' > %s)' %(thresh)

    return selStr
            
def MakeCatchList(inCATCHFCls, isAG, pdChnl, pdCatch, nAG, log, inHUC):
    '''Makes a list of catchments for DEP random sampling. Sorts results by
    percent agriculture and threshold adjusts to minimize number of watersheds
    processed.
    This was copied from cmd_gen1Flowpath_v6.py on 2020-02-20 to minimize code copies'''
    try:
        log.info("Select Catchments...")
        # Create column names
        pdcAcres = os.path.basename(pdCatch) + ".Acres"#"pdCatch" + huc12 + "_" + interpType + ".Acres"
        pdcSlope = os.path.basename(pdChnl) + ".Slope"#"pdChnl" + huc12 + "_" + interpType + ".Slope"
        pdcWSNO = os.path.basename(pdCatch) + ".WSNO"#"pdCatch" + huc12 + "_" + interpType + ".WSNO"


        # Tabulate Ag in each catchment (WSNO)
        taIsAg = TabulateArea(inCATCHFCls, "WSNO", isAG, "VALUE", "AGarea.dbf")
                
        catchLayer = arcpy.MakeFeatureLayer_management(inCATCHFCls, "CTCH.lyr")
        channelLayer = arcpy.MakeFeatureLayer_management(pdChnl, "CHNL.lyr")
        arcpy.AddJoin_management(catchLayer, "WSNO", taIsAg, "WSNO")
        arcpy.AddJoin_management(catchLayer, "WSNO", channelLayer, "WSNO") 

        # Select Catchments GT 25% Ag and 15 Acres
        thresh = .25#75
        maxCatch = 150
        minCatch = 100
####        selStr = '"%s.Acres" >= 15 AND "%s.Slope" > 0 AND (( AGarea.VALUE_1 * 0.000247) / "%s.Acres" > %s)' %(os.path.basename(pdCatch), os.path.basename(pdChnl), os.path.basename(pdCatch), thresh)
        selStr = getHucSelection(pdcAcres, pdcSlope, thresh)

    ## switched to list comprehension to handle errors when no catchments returned (del row would fail)
    ## bkgelder - 2019/03/20
    ## revised again to use with statement (for handling 0 returns)
    ## bkgelder - 2019/09/06
##        CatchList = [str(int(row[0])) for row in arcpy.da.SearchCursor(catchLayer, ["pdCatch%s.WSNO" %(inHUC)], selStr)]
        CatchList = []
        with arcpy.da.SearchCursor(catchLayer, [pdcWSNO], selStr) as scur:
####        with arcpy.da.SearchCursor(catchLayer, ["%s.WSNO" %(os.path.basename(pdCatch))], selStr) as scur:
            for srow in scur:
                CatchList.append(srow[0])
                                                
        
        catchCount = len(CatchList)
        log.info("ag sub-catchment count: " + str(catchCount) + " at " + str(thresh))

        if catchCount > 0:
            # Set higher or lower Catchments % Ag criteria to focus on subcatchments (in high ag), area GT 15 Acres
            if (catchCount < minCatch) or (catchCount > maxCatch):
                if catchCount > maxCatch:
                    while catchCount > maxCatch and thresh < 0.975:
                        thresh += 0.025
                        selStr = getHucSelection(pdcAcres, pdcSlope, thresh)
####                        selStr = '"%s.Acres" >= 15 AND "%s.Slope" > 0 AND (( AGarea.VALUE_1 * 0.000247) / "%s.Acres" > %s)' %(os.path.basename(pdCatch), os.path.basename(pdChnl), os.path.basename(pdCatch), thresh)
                
                        CatchList = [str(int(row[0])) for row in arcpy.da.SearchCursor(catchLayer, [pdcWSNO], selStr)]
####                        CatchList = [str(int(row[0])) for row in arcpy.da.SearchCursor(catchLayer, ["%s.WSNO" %(os.path.basename(pdCatch))], selStr)]
                        catchCount = len(CatchList)
                        log.info("   loop sub-catchment count: " + str(catchCount) + " at " + str(thresh))

                elif catchCount < minCatch:
                    while catchCount < minCatch and thresh > 0.25:
                        thresh -= 0.025
                        selStr = getHucSelection(pdcAcres, pdcSlope, thresh)
####                        selStr = '"%s.Acres" >= 15 AND "%s.Slope" > 0 AND (( AGarea.VALUE_1 * 0.000247) / "%s.Acres" > %s)' %(os.path.basename(pdCatch), os.path.basename(pdChnl), os.path.basename(pdCatch), thresh)

                        CatchList = [str(int(row[0])) for row in arcpy.da.SearchCursor(catchLayer, [pdcWSNO], selStr)]
####                        CatchList = [str(int(row[0])) for row in arcpy.da.SearchCursor(catchLayer, ["%s.WSNO" %(os.path.basename(pdCatch))], selStr)]
                        catchCount = len(CatchList)
                        log.info("   loop sub-catchment count: " + str(catchCount) + " at " + str(thresh))
        else:
            log.warning("agricultural sub-catchment count is 0, quitting")
        
        arcpy.Delete_management(catchLayer)
        arcpy.Delete_management(channelLayer)
        
        return(CatchList)
    except:
        
        # Get the traceback object
        #
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]

        # Concatenate information together concerning the error into a message string
        #
        pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
        msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"

        # Return python error messages for use in script tool or Python Window
        #
        ##            arcpy.AddError(pymsg)
        ##            arcpy.AddError(msgs)

        # Print Python error messages for use in Python / Python Window
        #
        log.warning(pymsg)
        log.warning(msgs)

        log.warning('failure on: ' + inHUC)

    

def defineLocalProc(node, uversion = ''):
    #Alan Kuutilla's naming for laptops
    if '-M' in node.upper():
        localProc = 'C:\\DEP_Proc'
    elif 'EL3354-02' in node.upper() or 'EL3321-02' in node.upper():
        localProc = 'D:\\DEP_Proc'
    elif 'DA214B-12' in node.upper() or 'DA214B-11' in node.upper():
        localProc = 'D:\\DEP_Proc'
    elif 'DEP' in node.upper():
        localProc = 'D:\\DEP_Proc'
    else:
        localProc = 'C:\\DEP_Proc'

    if uversion != '':
        localProc += uversion

    return localProc
    

def createBasicDirectories(node, ACPFyear, uversion = ''):
    '''Creates the most basic DEP directories, ACPF directory and DEP base directory as strings'''
    # depBase - location of most DEP inputs and outputs
    # acpfBase - location of DEP/ACPF management data geodatabases
    if 'EL3354-02' in node.upper():# or 'DA214B-11' in node.upper():
        acpfStart = 'D:\\DEP'
        depBase = 'M:\\DEP'
    elif 'EL3321-M10' in node.upper():
        # run everything local on laptop
        acpfStart = 'C:\\DEP'
        depBase = 'C:\\DEP'
    else:
        acpfStart = '\\\\EL3354-02\\D$\\DEP'#\\Man_Data_ACPF\\dep_ACPF' + ACPFyear
        depBase = '\\\\EL3354-02\\M$\\DEP'

    # basedata never changes with version
    basedataDir = opj(acpfStart, 'Basedata_Summaries')

    # see above
    acpfBase = opj(acpfStart + uversion, 'Man_Data_ACPF\\dep_ACPF'+ ACPFyear)
    # otherBase is raster residue cover maps
    otherBase = opj(depBase, 'Man_Data_Other')
    depBase += uversion

    return acpfBase, depBase, basedataDir, otherBase, acpfStart


def loadInterpDict():
    '''Creates a dicationary based on the four main categories of Terrain pyramid
    options that are then used to create different DEMs.'''
    interpDict = {"ZMEAN": 'mean18', "ZMINMAX": 'mnmx18', "ZMIN": 'min18', "ZMAX": 'max18'}
    return interpDict

def loadBasicVariablesDict(node, ACPFyear, uversion = ''):
    '''Creates the locations of the HUC, County, and State data is stored as this
    depends on the required processing steps and year of the input data.'''
##    node = listy[0]
##    ACPFyear = listy[1]
##    j = ()

    acpfBase, depBase, basedataDir, otherBase, acpfStart = createBasicDirectories(node, ACPFyear, uversion)

    if int(ACPFyear) > 2021:
        MWHUC12_ACPF = 'MW_HUC12_v2022'
        MWHUC8_ACPF = 'MW_HUC8_v2022'
        MWHUC2_ACPF = 'MW_HUC2_v2022'
    elif int(ACPFyear) > 2018:
        MWHUC12_ACPF = 'MW_HUC12_v2019'
        MWHUC8_ACPF = 'MW_HUC8_v2019'
        MWHUC2_ACPF = 'MW_HUC2_v2019'
    else:
        MWHUC12_ACPF = 'MW_HUC12_v2013'
        MWHUC8_ACPF = 'MW_HUC8_v2013'
        MWHUC2_ACPF = 'MW_HUC2_v2013'

    basicDict = {
    "depBase" : depBase,
    "acpfBase" : acpfBase,
    "acpfStart" : acpfStart,
    "basedataDir" : basedataDir,
    "otherBase" : otherBase,
    "mwHuc12s5070" : opj(basedataDir, 'Basedata_5070.gdb', MWHUC12_ACPF),
    "mwHuc8s5070" : opj(basedataDir, 'Basedata_5070.gdb', MWHUC8_ACPF),
    "mwHuc2s5070" : opj(basedataDir, 'Basedata_5070.gdb', MWHUC2_ACPF),

    "mwCounties5070" : opj(basedataDir, 'Basedata_5070.gdb','MW_Counties_PCS'),
    "mwStates5070" : opj(basedataDir, 'Basedata_5070.gdb','MW_States'),
    "main_wesm" : opj(basedataDir, 'WESM.gpkg','main.wesm'),
        # huc_8_fields gdb
    "huc8_fields" : opj(acpfBase, 'fields_by_huc8.gdb')
    }

    return basicDict, MWHUC12_ACPF, MWHUC8_ACPF, MWHUC2_ACPF


##def loadVariablesDict(listy):
def loadVariablesDict(node, ACPFyear, huc12, outEPSG, interpType, cellSize, nowYmd, version = ''):
    '''This function creates the locations where the DEP data is stored as this
    depends on the required processing steps and year of the input data. The output
    list is subsetted based on the subset argument. Version adds a version identifier
    to all output and creates a copy ACPF directory.
    node - platform.node()
    ACPFyear - year of ACPF data to use, as string
    huc12 - 12 digit string
    outEPSG - coordinate reference system as string
    interpType - Terrain code abbreviation (e.g. mean18, mnmx18)
    cellSize - cell size as integer
    nowYmd - %Y_%m_%d_%H_%M_%S to append to filenames
    version - identification string to code things uniquely for testing
    '''

    assert len(huc12) == 12, 'HUC12 code too short'
    assert 4 <= len(outEPSG) <= 5, 'coordinate system EPSG should be 4 or 5 characters'
    assert type(cellSize) == type(1), 'cellSize should be integer'

    huc12_2022 = ['051201130904', '051402020705', '050800020907', '071200030304', '071200030302', '071200011308', '040400010510', '040400010103', '040400010302', '040400010303', '040400010401', '040400010403', '040400010501', '040400010502', '040400010503', '040400010504', '040400010505', '040400010508', '071200011004', '071200011008', '040400010204', '040400010104', '040400010201', '040400010205', '071200010201', '071200030303', '071200030203', '071200030407', '071100090107', '102300030501', '040400010402', '071401010401', '071100090403', '071100090401', '071401010206', '102100100206', '102001011107', '102002020103', '102702010101', '102002020211', '102200031003', '040400010506', '071300111006', '102400010106', '051402060703', '040801010504', '040700030401', '041000100503', '080101000201', '040400010604', '040301140201', '040301140101', '040601040401', '040400020101', '040400010102', '040400010605', '040400010602', '040400010601', '040400010105', '040101010105', '042600000101', '040400010507', '041800000101', '041800000102', '090300010301', '040101010101', '040101010208', '040101010207', '040101011504', '040101010203', '040101011403', '040101010108', '040101010107', '040101011404', '040101011503', '040101011502', '040101010103', '040101010106', '040101011401', '040101010104', '040101011402', '040101010102', '040101011501', '040900020102', '040900020901', '040900020702', '040900020304', '040900010504', '040900020203', '040900020301', '040900020602', '040900020401', '040900020101', '040900020701', '040900010501', '040900040602', '040900040603', '040900040502', '040900040601', '040900040702', '040900040703', '040900040701', '040900021008', '040900021007', '040900021006', '040900021005', '040900040503', '040900040501', '040900020503', '040900021004', '040900021001', '040900020902', '040900020801', '040900021009', '040900020303', '040900020803', '040900020302', '040900020703', '040900021002', '040900020402', '040900020501', '040900020502', '040900020201', '040900020802', '040900021003', '040900020603', '040900020204', '040900020202', '040900020504', '040900020601', '041000130106', '041000130405', '041000130104', '041000130101', '041000130102', '041000130103', '041000130105', '041000130107', '041000130111', '041000130108', '041000130109', '041000130110', '041000130112', '041000130201', '041000130301', '041000130202', '041000130203', '041000130204', '041000130302', '041000130303', '041000130304', '041000130404', '041000130401', '041000130402', '041000130403', '041000130305', '041000130306', '041000130307', '041000130308', '041000130309', '041000130406', '041000130407', '041900000100', '040900010502', '101000042706', '101102050609', '101101010101', '101101012903', '042000021003', '042000021005', '042000021105', '042000021104', '042000021006', '042000021004', '042000021002', '042000021007', '042000021001', '042000021103', '042000021101', '042000021102', '042000020301', '042000020302', '042000020303', '042000020903', '042000020603', '042000020606', '042000020602', '042000020601', '042000020902', '042000020706', '042000020901', '042000020604', '042000020801', '042000020201', '042000020401', '042000020202', '042000020704', '042000020702', '042000020802', '042000020609', '042000020701', '042000020402', '042000020501', '042000020502', '042000020608', '042000020703', '042000020607', '042000020605', '040301020104', '040900010503', '042600000102', '041800000200', '042000020203', '042000020705', '042600000200', '051401010304', '101702040602', '070200080601', '090201060401', '090201060301', '103001010202', '103001010203', '051202070601', '080101000103', '080102010506', '102001030507', '102002020105', '102100090605', '102200040406', '102200010801', '102200030805', '102200031002', '102200031006', '101800091707', '101800140102', '070900021401', '071200060402', '050902030204', '040400010206', '042400020200', '101701040109', '101701040202', '101701040304', '101701040103', '101701040107', '101701040101', '101701040102', '101701040104', '101701040105', '101701040106', '101701040108', '101701040110', '101701040201', '101701040203', '101701040204', '101701040207', '101701040205', '101701040206', '101701040301', '101701040305', '101701040302', '101701040303', '071401010403', '071401010502', '041000010201', '041000010206', '040400010606', '071200011002', '040400010509', '071200030306', '040400010301', '041900000200', '090201030805']

    huc12_2019 = ['080203020101', '071100090201', '071100090301', '071402040501', '071402040506', '051201120204', '070200090301', '101702040504', '040400010603', '051201090801', '102300020301', '071200021201', '071200030104', '071200030107', '051201130901', '051201130902', '051201130903', '051201130904', '071200011308', '071200020703', '071402010206', '071200030304', '071200030305', '051402020705', '071300020104', '051201141004', '071300120401', '071000030403', '102400010106', '070600030704', '070700051804', '051402040501', '051402060507', '051402060508', '051402060509', '102500100303', '102500100304', '102500110303', '102500110304', '102500160804', '102701020906', '102701030510', '110100070804', '070400080904', '102300030501', '102801011707', '102801011708', '102801021404', '070300040403', '071000010502', '071000010801', '102802020410', '102901030205', '102901050204', '070102060903', '070300050404', '070400010209', '070102070401', '070400030601', '070400060502', '070600010504', '070102070602', '070200080503', '103001011204', '070102010101', '070102020504', '080103000101', '080103000102', '080103000103', '080103000104', '080103000201', '080103000202', '080103000203', '080103000204', '080103000301', '080103000302', '080103000303', '080103000304', '080103000305', '080103000306', '080103000307', '080103000308', '080103000309', '101900180606', '101800140808', '102500090701', '102500090703', '102500090901', '101900180703', '101900180705', '101900180706', '102001010901', '102001010902', '102001010903', '102001011004', '102001011005', '102001011101', '102001011102', '102500160301', '102500160303', '102001011106', '102001011107', '102001030206', '102002020101', '102002020103', '102702010101', '102002020211', '102702030401', '102702060101', '102702060301', '102702060901', '102100090503', '102200031003', '102200031005', '070700031401', '040301020102', '040301020103', '040301020105', '040301020106', '040301020107', '040301020108', '040301020109', '040301020110', '040301020201', '040301020204', '040301020205', '040301020305', '040301020401', '040301020403', '040301020404', '040301020405', '040301020406', '040301020407', '070400070701', '070500051205', '070500010904', '070500020704', '090203020104', '041800000200', '040103020109', '041800000100', '040103010504', '040103010503', '040103010601', '040103010602', '040103010603', '040103010607', '040103010608', '040103010705', '040103010806', '040103010807', '040103010901', '040103010902', '040103010903', '040103010904', '040103010906', '040103011002', '040103011003', '040103011004', '040103011008', '040103011101', '040103011105', '040103020703', '040103020705', '090300050302', '040101010802', '040102011604', '090300010302', '090300070501', '040102010602', '040101010401', '040101011301', '090300091423', '071402040504', '071402040505', '071402040502', '090203140101', '090203140102', '090203140103', '090203140104', '090203140105', '090203140106', '090203140107', '090203120606', '090203120601', '090203120502', '090203140408', '090203110702', '090203120602', '090203140604', '090203110701', '090203140504', '090203140603', '090203120501', '090203140508', '090203140606', '090203140607', '090203140409', '090203110803', '090203140407', '090203110703', '090203140808', '090203140405', '090203110804', '090203140801', '090203140507', '090203140608', '090203140809', '090203140602', '090203140406', '090203140806', '040400020503', '040400020502', '040400020501', '040400020401', '040400020403', '040400020306', '040400020102', '040400030606', '040400030601', '040301011203', '040301011202', '040301011109', '040301010705', '040301010704', '040301010703', '040301010702', '040301010605', '040301010604', '040301010205', '040301010101', '040302040405', '040302040401', '040302040106', '040301030205', '040301030104', '040301030103', '040301040506', '040301050607', '040301050606', '040301050605', '040301080913', '040400020101', '041900000002', '041900000001', '040302040101', '041800000300', '040301020104', '040301020101', '040301020402', '040301020111']
    if int(ACPFyear) > 2021:
        if huc12 in huc12_2022:
            ACPFDEMyear = '2022'
        elif huc12 in huc12_2019:
            ACPFDEMyear = '2019'
        else:
            ACPFDEMyear = '2013'
##        print('ACPFDEMyear is: ' + ACPFDEMyear)

    elif int(ACPFyear) > 2018:
        if huc12 in huc12_2019:# and int(ACPFyear) > 2018:
            ACPFDEMyear = '2019'
##        print('ACPFDEMyear is: ' + ACPFDEMyear)
        else:
            ACPFDEMyear = '2013'
##        print('ACPFDEMyear is: ' + ACPFDEMyear)

    else:
        ACPFDEMyear = '2013'

    dem_year = '_dem' + ACPFDEMyear

    DEMyear = 'Current'

    if version != '':
        if version.startswith('_') or version.endswith('_'):
            version = version.replace('_', '')
        uversion = '_' + version
        version = version + '_'
    else:
        uversion = ''

    uinterpType = '_' + interpType
    ucellSize = '_' + str(cellSize) + 'm'
    uhuc12 = '_' + huc12


    try:
##        print('uinterpType: ' + uinterpType)
##        print('uversion: ' + uversion)
##        print('dem_year: ' + dem_year)
##        print('node: ' + node)
        # localProc - local processing directory (a local drive)
        # get the local processing directory (usually D:)
        localProc = defineLocalProc(node, uversion)
##        print('localProc: ' + localProc)
        # basedataDir - location of the PCS specific datasets

        locationsDict, MWHUC12_ACPF, MWHUC8_ACPF, MWHUC2_ACPF = loadBasicVariablesDict(node, ACPFyear, uversion)
        acpfBase = locationsDict['acpfBase']
        depBase = locationsDict['depBase']
        basedataDir = locationsDict['basedataDir']
        otherBase = locationsDict['otherBase']
        locationsDictNoVersion, MWHUC12_ACPF_no_version, MWHUC8_ACPF_no_version, MWHUC2_ACPF_no_version = loadBasicVariablesDict(node, ACPFyear, '')
        otherBaseNoVersion = locationsDictNoVersion['otherBase']
        acpfBaseNoVersion = locationsDictNoVersion['acpfBase']

####        acpfBase, depBase, basedataDir = createBasicDirectories(node, ACPFyear, uversion)
##        print(', '.join([acpfBase, depBase, basedataDir]))

        lidarOutputDir = opj(depBase, 'LiDAR_' + DEMyear)
        
        depFlowpaths = opj(depBase, 'DEP_Flowpaths')

        depDocumentation = opj(basedataDir, 'Documentation')

        depMetadata = opj(depBase, 'toolMetadata')

        gullyOutputs = opj(depBase, 'Gully_Outputs')

        depOutput = opj(locationsDict['acpfStart'], 'DEP_Outputs')

        huc8 = huc12[0:8]
        uhuc12 = '_' + huc12

        pElevDir = os.path.join(lidarOutputDir, '_'.join(['elev', 'PLib', interpType]))# + '_' + outEPSG)# + uversion)
        fElevDir = os.path.join(lidarOutputDir, '_'.join(['elev', 'FLib', interpType]))# + '_' + outEPSG)# + uversion)
        cElevDir = os.path.join(lidarOutputDir, '_'.join(['elev', 'CLib', interpType]))# + '_' + outEPSG)# + uversion)
        vElevDir = os.path.join(lidarOutputDir, '_'.join(['elev', 'VLib', interpType]))# + '_' + outEPSG)# + uversion)

        firstReturnDir = os.path.join(lidarOutputDir, '_'.join(['surf', 'el', 'Lib']))# + outEPSG)# + uversion)
        countDir = os.path.join(lidarOutputDir, '_'.join(['count', 'Lib']))# + outEPSG)# + uversion)
        intDir = os.path.join(lidarOutputDir, '_'.join(['int', 'Lib']))# + outEPSG)# + uversion)

        copyDir = os.path.join(lidarOutputDir, '_'.join(['copy', 'gdb', interpType]))# + '_' + outEPSG)# + uversion)

        holesDir = os.path.join(lidarOutputDir, '_'.join(['holes', 'Lib', interpType]))# + '_' + outEPSG)# + uversion)

        voidDir = os.path.join(lidarOutputDir, '_'.join(['voids', 'Lib', interpType]))# + '_' + outEPSG)# + uversion)

        huc8Dir = os.path.join(lidarOutputDir, '_'.join(['huc8', outEPSG]))# + uversion)

##        huc8gdb = opj(huc8Dir, 'huc_' + huc8 + uversion + '.gdb')
        huc8gdb = opj(huc8Dir, 'huc_' + huc8 + '.gdb')

##        copyGDB = os.path.join(copyDir, huc8, 'copies_' + huc12 + uversion + '.gdb')
        copyGDB = os.path.join(copyDir, huc8, 'copies_' + huc12 + '.gdb')
        
        ACPFDir = os.path.join(acpfBase, huc8, 'idepACPF' + huc12 + '.gdb')
##        if not arcpy.Exists(ACPFDir):
##            if not os.path.isdir(os.path.dirname(ACPFDir)):
##                os.makedirs(os.path.dirname(ACPFDir))
##            acpfOriginal = ACPFDir.replace(uversion, '')
##            acpfVersionCopy = arcpy.Copy_management(acpfOriginal, ACPFDir)
##        localACPFDir = ACPFDir

        first_of_month = datetime.datetime.today().replace(day=1)
        ept_first_of_month_name = "ept_resources_" + first_of_month.strftime("%Y_%m_%d")

        locationsDict.update({
        # pit filled elevation file
        "fElevFile" : os.path.join(fElevDir, huc8, 'ef' + str(cellSize) + 'm' + huc12 + '.tif'),
        # punched elevation file
        "pElevFile" : os.path.join(pElevDir, huc8, 'ep' + str(cellSize) + 'm' + huc12 + '.tif'),
        # cut elevaiton file
        "cElevFile" : os.path.join(cElevDir, huc8, 'ec' + str(cellSize) + 'm' + huc12 + '.tif'),
        # lake-fixed elevation file
        "vElevFile" : os.path.join(vElevDir, huc8, 'ev' + str(cellSize) + 'm' + huc12 + '.tif'),
        # bridge-fixed elevation file
        "wElevFile" : os.path.join(vElevDir, huc8, 'ew' + str(cellSize) + 'm' + huc12 + '.tif'),
        # channel-fixed elevation files
        "xElevFile" : os.path.join(vElevDir, huc8, 'ex' + str(cellSize) + 'm' + huc12 + '.tif'),
        # channel-fixed elevation files
        "yElevFile" : os.path.join(vElevDir, huc8, 'ey' + str(cellSize) + 'm' + huc12 + '.tif'),
        # channel-fixed elevation files
        "zElevFile" : os.path.join(vElevDir, huc8, 'ez' + str(cellSize) + 'm' + huc12 + '.tif'),
        # first return max elevation file
        "surfaceElevFile" : os.path.join(firstReturnDir, huc8, 'frmax' + str(cellSize) + 'm' + huc12 + '.tif'),
        # first return min elevation file
        "firstReturnMinFile" : os.path.join(firstReturnDir, huc8, 'frmin' + str(cellSize) + 'm' + huc12 + '.tif'),
        # last return min elevation file
        "lastReturnMinFile" : os.path.join(firstReturnDir, huc8, 'lrmin' + str(cellSize) + 'm' + huc12 + '.tif'),
        # bareearth return min elevation file, created from lasdataset to raster, not terrain to raster
        "bareEarthReturnMinFile" : os.path.join(firstReturnDir, huc8, 'bemin' + str(cellSize) + 'm' + huc12 + '.tif'),

        # holes that were punched in punched elevation file
        "holesFile" : os.path.join(holesDir, huc8, 'holes' + str(cellSize) + 'm' + huc12 + '.tif'),

        # bare earth return count file
        "cntFile" : os.path.join(countDir, huc8, 'cbe' + str(cellSize) + 'm' + huc12 + '.tif'),
        # first return count file
        "cnt1rFile" : os.path.join(countDir, huc8, 'cfr' + str(cellSize) + 'm' + huc12 + '.tif'),

        # first return min intensity file
        "int1rMinFile" : os.path.join(intDir, huc8, 'fr_int_min' + str(cellSize) + 'm' + huc12 + '.tif'),
        # first return max intensity file
        "int1rMaxFile" : os.path.join(intDir, huc8, 'fr_int_max' + str(cellSize) + 'm' + huc12 + '.tif'),
        # bare earth return min intensity file
        "intBeMinFile" : os.path.join(intDir, huc8, 'be_int_min' + str(cellSize) + 'm' + huc12 + '.tif'),
        # bare earth return max intensity file
        "intBeMaxFile" : os.path.join(intDir, huc8, 'be_int_max' + str(cellSize) + 'm' + huc12 + '.tif'),

        # documentation folder
        "docFolder" : depDocumentation,

        # metadata documents
        "flib_metadata" : opj(depMetadata, 'FLib_DEMs2022_mTemplate.xml'),
        "derivative_metadata" : opj(depMetadata, 'FLib_Derivatives2022_mTemplate.xml'),

        # no data areas (big voids) to fix first (likely rivers/lakes/lagoons)
        "bigNoData" : opj(voidDir, huc8, 'bigvds' + str(cellSize) + 'm' + huc12 + '.tif'),
        # no data areas (smaller voids) to fix second (likely bridge decks/culverts)
        "mediumNoData" : opj(voidDir, huc8, 'medvds' + str(cellSize) + 'm' + huc12 + '.tif'),
        # no data areas (smallest voids) to fix third (likely streams)
        "smallNoData" : opj(voidDir, huc8, 'smlvds' + str(cellSize) + 'm' + huc12 + '.tif'),

        # other datasets where feature, thus name, is based on DEM year vintage
##      new version of file names
        "pdCatch" : opj(ACPFDir, "pdCatch" + uinterpType + dem_year + ucellSize + uhuc12),
        "pdChnl" : opj(ACPFDir, "pdChnl" + uinterpType + dem_year + ucellSize + uhuc12),
        "wShed" : opj(ACPFDir, "WShed" + uinterpType + dem_year + ucellSize + uhuc12),

        # old version of feature class names
####        "pdCatch" : opj(ACPFDir, "pdCatch" + huc12 + uversion + uinterpType + dem_year),
####        "pdChnl" : opj(ACPFDir, "pdChnl" + huc12 + uversion + uinterpType + dem_year),
####        "wShed" : opj(ACPFDir, "WShed" + huc12 + uversion + uinterpType + dem_year),

        "fieldBoundaries" : opj(ACPFDir, "FB" + huc12),
        "LU6" : opj(ACPFDir, "LU6_" + huc12),
        "SSURGO" : opj(ACPFDir, 'gSSURGO'),
        # "rc_table" : opj(ACPFDir, "ResCover_" + huc12),
        "manField" : 'Management_CY_' + ACPFyear,
        "rotField" : 'CropRotatn_CY_' + ACPFyear,
        "tillField": 'Till_code_CY_' + ACPFyear,
        "resCoverField": 'Adj_RC_CY_' + ACPFyear,

            # for running cutter
        "watershedBoundaries" : opj(ACPFDir, "bnd" + huc12),
        "bufferedBoundaries" : opj(ACPFDir, "buf_" + huc12),
        "wesm_project_boundaries" : opj(ACPFDir, "_".join(["wesm", ept_first_of_month_name, huc12])),
        "snapraster" : opj(basedataDir, 'Basedata_' + outEPSG + '.gdb/Snap1m'),
        "roadsfc" : opj(basedataDir, 'Basedata_' + outEPSG + '.gdb/roads_merge'),
        "rrsfc" : opj(basedataDir, 'Basedata_' + outEPSG + '.gdb/railways_merge'),
        "rwsfc" : opj(basedataDir, 'Basedata_' + outEPSG + '.gdb/runways'),
        "waterwaysfc" : opj(basedataDir, 'Basedata_' + outEPSG + '.gdb/waterways'),
        "waterfc" : opj(basedataDir, 'Basedata_' + outEPSG + '.gdb/water'),

        "mn_dnr_forest_roads" : opj(basedataDir, 'Basedata_' + outEPSG + '.gdb/mn_dnr_roads'),
        "us_fs_forest_roads" : opj(basedataDir, 'Basedata_' + outEPSG + '.gdb/us_fs_roads'),
        "us_fs_forest_trails" : opj(basedataDir, 'Basedata_' + outEPSG + '.gdb/us_fs_trails'),

        # other datasets where feature, thus name, is based on DEM year vintage
        "goodcutsfc" : opj(ACPFDir, "cuts_prelim" + uinterpType + dem_year + ucellSize + uhuc12),
        "bestcutsfc" : opj(ACPFDir, "cuts_final" + uinterpType + dem_year + ucellSize + uhuc12),
        "depressionsfc" : opj(ACPFDir, "dprsns" + uinterpType + dem_year + ucellSize + uhuc12),
        "depressions2cutfc" : opj(ACPFDir, "dprsns2cut" + uinterpType + dem_year + ucellSize + uhuc12),
        "best_cut_stf" : opj(ACPFDir, "best_cut_stf" + uinterpType + dem_year + ucellSize + uhuc12),
        "void_pit_dif" : opj(ACPFDir, "void_pit_dif" + uinterpType + dem_year + ucellSize + uhuc12),
        "punch_pit_dif" : opj(ACPFDir, "punch_pit_dif" + uinterpType + dem_year + ucellSize + uhuc12),
        "fil_dif_cut" : opj(ACPFDir, "fil_dif_cut" + uinterpType + dem_year + ucellSize + uhuc12),
        "fil_dif_delta" : opj(ACPFDir, "fil_dif_delta" + uinterpType + dem_year + ucellSize + uhuc12),
        "fil_dif_pit" : opj(ACPFDir, "fil_dif_pit" + uinterpType + dem_year + ucellSize + uhuc12),
        "layerPath" : opj(depBase, 'Layers'),
        "overviewPath" : opj(depBase, 'Overview_Maps'),

##        "goodcutsfc" : opj(ACPFDir, "cuts_prelim" + uversion + uinterpType + dem_year + _huc12),
##        "bestcutsfc" : opj(ACPFDir, "cuts_final" + uversion + uinterpType + dem_year + _huc12),
##        "depressionsfc" : opj(ACPFDir, "dprsns" + uversion + uinterpType + dem_year + _huc12),
##        "depressions2cutfc" : opj(ACPFDir, "dprsns2cut" + uversion + uinterpType + dem_year + _huc12),
##        "layerPath" : opj(lidarOutputDir, 'Layers'),

        "mwHuc12s" : opj(basedataDir, 'Basedata_' + outEPSG + '.gdb', MWHUC12_ACPF),
        "mwHuc12s5070_stats" : opj(basedataDir, 'Basedata_5070.gdb', MWHUC12_ACPF + '_'.join(['', 'Status', str(cellSize) + 'm', interpType])),
        "mwHuc8s" : opj(basedataDir, 'Basedata_' + outEPSG + '.gdb', MWHUC8_ACPF),
        "mwHuc2s" : opj(basedataDir, 'Basedata_' + outEPSG + '.gdb', MWHUC2_ACPF),

        "medianProcDir" : opj(localProc,  'Median_Proc', '_'.join(['Medians', outEPSG, huc8])),
        "mergedMdnsHuc8FC" : opj(huc8gdb, 'rd_rr_rd_rw_mrg_' + huc8),
        "huc8Roads" : opj(huc8gdb, 'roads_' + huc8),

        "GordRaster" : opj(depFlowpaths, '_'.join(['HUC12', 'GridOrder', interpType]), huc8, 'gord_' + huc12 + '.tif'),
        "FPath_Mosaic_out" : opj(depFlowpaths, '_'.join(['HUC12', 'FlowPaths', interpType]), huc8, 'fp' + huc12 + '.tif'),
        "FPLen_Mosaic_out" : opj(depFlowpaths, '_'.join(['HUC12', 'FPLengths', interpType]), huc8, 'fpLen' + huc12 + '.tif'),

        # rasters to use in calculating gully erosion outputs and metrics
        "tau_slope" : opj(gullyOutputs, 'Tau_D8_Slope' + '_' + interpType, huc8, 'sd8_' + str(cellSize) + 'm_' + huc12 + '.tif'),
        "ag_fd" : opj(gullyOutputs, 'AG_FD' + '_' + interpType, huc8, 'fd_' + str(cellSize) + 'm_' + huc12 + '.tif'),
        "ag_fa" : opj(gullyOutputs, 'AG_FA' + '_' + interpType, huc8, 'fa_' + str(cellSize) + 'm_' + huc12 + '.tif'),
        "ag_plan_crv" : opj(gullyOutputs, 'AG_Plan_Crv' + '_' + interpType, huc8, 'pln_crv_' + str(cellSize) + 'm_' + huc12 + '.tif'),

        # gully erosion terrain index rasters
        "stream_power_raster" : opj(gullyOutputs, 'spi' + '_' + interpType, huc8, 'stream_power_' + str(cellSize) + 'm_' + huc12 + '.tif'),
        "modified_stream_power_raster" : opj(gullyOutputs, 'mod_spi' + '_' + interpType, huc8, 'mod_strm_pwr_' + str(cellSize) + 'm_' + huc12 + '.tif'),
        "compound_topo_index_raster" : opj(gullyOutputs, 'cti' + '_' + interpType, huc8, 'compound_topo_' + str(cellSize) + 'm_' + huc12 + '.tif'),
        "specific_contributing_area_raster" : opj(gullyOutputs, 'contrib_area' + '_' + interpType, huc8, 'spec_con_area_' + str(cellSize) + 'm_' + huc12 + '.tif'),
        "flup_raster" : opj(gullyOutputs, 'flow_length' + '_' + interpType, huc8, 'flup_' + str(cellSize) + 'm_' + huc12 + '.tif'),

        "samples" : opj(ACPFDir, 'smpl' + str(cellSize) + 'm_' + interpType + huc12),
        "samples_defined" : opj(ACPFDir, 'smpldef' + str(cellSize) + 'm_' + interpType + huc12),
        "nulls" : opj(ACPFDir, 'null' + str(cellSize) + 'm_' + interpType + huc12),
        "null_flowpaths" : opj(ACPFDir, 'null_flowpaths' + str(cellSize) + 'm_' + interpType + huc12),
        "areas" : opj(ACPFDir, 'areas' + str(cellSize) + 'm_' + interpType + huc12),
        "snaps" : opj(ACPFDir, 'snaps' + str(cellSize) + 'm_' + interpType + huc12),
        # "tillages" : opj(ACPFDir, 'tillage' + str(cellSize) + 'm' + huc12),

##        "samples" : opj(ACPFDir, 'smpl' + str(cellSize) + 'm_' + interpType + '_' + huc12),
##        "samples_defined" : opj(ACPFDir, 'smpldef' + str(cellSize) + 'm_' + interpType + '_' + huc12),
##        "nulls" : opj(ACPFDir, 'null' + str(cellSize) + 'm_' + interpType + '_' + huc12),
##        "areas" : opj(ACPFDir, 'areas' + str(cellSize) + 'm_' + interpType + '_' + huc12),
##        "snaps" : opj(ACPFDir, 'snaps' + str(cellSize) + 'm_' + interpType + '_' + huc12),
##        "tillages" : opj(ACPFDir, 'tillage' + str(cellSize) + 'm' + '_' + huc12),
##
        "blFile" : opj(lidarOutputDir, 'bl_Lib\\' + huc8 + '\\breaklines_' + huc12 + '.shp'),
        "breaklines" : opj(lidarOutputDir, 'bl_Lib', huc8, "breaks_" + huc8 + ".gdb", 'break_lines_' + huc12),
        "breakpolys" : opj(lidarOutputDir, 'bl_Lib', huc8, "breaks_" + huc8 + ".gdb", 'break_polys_' + huc12),
        "lasTilesDrive" : opj(basedataDir, 'Basedata_5070.gdb','IA_MN_NE_Tiles_Merge'),
        "kansasTilesPath" : opj(basedataDir, 'Basedata_5070.gdb','Kansas_tiles_available_actual_extents'),

        "stripsFlumes_old" : opj(depBase, 'temp\\STRIPS2_flume_points\\STRIPS2_flume_points.gdb\\flume_points'),
        "ephemeralGullies_old" : opj(depBase, 'temp\\Gord_cmwang\\ruraldata\\data' + huc12 + '.gdb\\Gullyhead'),

        "stripsFlumes" : opj(otherBaseNoVersion, 'Gully_data\\STRIPS2_flume_points\\STRIPS2_flume_points.gdb\\flume_points'),
        "ephemeralGullies" : opj(otherBaseNoVersion, 'Gully_data\\cmwang\\data' + huc12 + '.gdb\\Gullyhead'),
        "nrcs_gullies" : opj(otherBaseNoVersion, 'Gully_data\\NRCS\\NRCS_Gulley_' + huc12 + '.gdb\\GullyHeads_' + huc12),

        # processing directories
        "peukProcDir" : opj(localProc, 'DEMProc', 'Catch' + dem_year + '_' + str(cellSize) + 'm_' + huc12),
        "fpProcDir" : opj(localProc, 'DEMProc', 'Flowpath' + dem_year + '_' + str(cellSize) + 'm_' + huc12),
        "manProcDir" : opj(localProc, 'DEMProc', 'Manage' + dem_year + '_' + str(cellSize) + 'm_' + huc12),
        "samplerProcDir" : opj(localProc, 'DEMProc', 'Sample' + dem_year + '_' + str(cellSize) + 'm_' + huc12),
        "fProcDir" : opj(localProc, "DEMProc", "LAS" + dem_year + '_' + str(cellSize) + 'm_' + huc12),
        "cutProcDir" : opj(localProc, "DEMProc", "Cut" + dem_year + '_' + str(cellSize) + 'm_' + huc12),
        "voidProcDir" : opj(localProc, "DEMProc", "Void" + dem_year + '_' + str(cellSize) + 'm_' + huc12),
        "otherProcDir" : opj(localProc, "DEMProc", "Other" + dem_year + '_' + str(cellSize) + 'm_' + huc12),

        # focal statistics kernels
        "EKernelFile" : opj(depBase, 'kernels', 'FlowE_Kernel.txt'),
        "SEKernelFile" : opj(depBase, 'kernels', 'FlowSE_Kernel.txt'),
        "SKernelFile" : opj(depBase, 'kernels', 'FlowS_Kernel.txt'),
        "SWKernelFile" : opj(depBase, 'kernels', 'FlowSW_Kernel.txt'),
        "WKernelFile" : opj(depBase, 'kernels', 'FlowW_Kernel.txt'),
        "NWKernelFile" : opj(depBase, 'kernels', 'FlowNW_Kernel.txt'),
        "NKernelFile" : opj(depBase, 'kernels', 'FlowN_Kernel.txt'),
        "NEKernelFile" : opj(depBase, 'kernels', 'FlowNE_Kernel.txt'),

        # residue cover directory
        "rescoverDir" : opj(os.path.dirname(basedataDir), 'Man_Data_Other\\GEE_residue_cover'),

        # DEP SSURGO soils directory
        # altered SOL file to remove SOL = ACPF + 1 year assumption, on consultation with David James, 2021.10.28 bkgelder
##        "soilsDir" : opj(acpfBase, 'Man_Data_ACPF\\dep_WEPP_SOL'+ str(int(ACPFyear) + 1)),
        "soilsDir" : opj(os.path.dirname(acpfBaseNoVersion), 'dep_WEPP_SOL'+ str(int(ACPFyear) + 1)),
        # Dave switched it back, 2022.04.03 bkgelder
##        "soilsDir" : opj(os.path.dirname(acpfBase), 'dep_WEPP_SOL'+ str(int(ACPFyear))),

        "eleBaseDir" : opj(depBase, 'Elev_Base_Data'),

        "projRasterDir" : opj(depBase, 'proj_rasters')})


        if int(ACPFyear) <= 2015:
            geeResidueMap = opj(otherBaseNoVersion, 'Iowa_RC.gdb\\huc' + huc8 + '_ACPF2017')
        elif int(ACPFyear) <= 2025:
            geeResidueMap = opj(otherBaseNoVersion, 'GEE_residue_cover', 'residue_cover_' + ACPFyear + '.tif')
        if int(ACPFyear) <= 2017:
            mnResidueMap = opj(otherBaseNoVersion, 'Minnesota\\2017_combo_agonly_recount_combo_residue_times2.img')
        elif int(ACPFyear) <= 2020:
            mnResidueMap = opj(otherBaseNoVersion, 'Minnesota\\minnesota_2020_boax_rowcrop_noforage_residue_times2_rc200_new_combo2019_2020_model.img')
        elif int(ACPFyear) == 2021:
            mnResidueMap = opj(otherBaseNoVersion, 'Minnesota\\minnesota_2021_boax_rowcrop_noforage_residue_times2_rc200_new_combo2019_2020_model.img')
        elif int(ACPFyear) == 2022:
            mnResidueMap = opj(otherBaseNoVersion, 'Minnesota\\minnesota_s2_boa_wgs15x_residue_times2_rowcrops_noforage_final_w_flood_mask_out.tif')
        
        locationsDict.update({
        "irrigationMap" : opj(otherBaseNoVersion, 'lanid2011-2017', 'lanid2017.tif'),
        "canopyCoverMap" : opj(otherBaseNoVersion, 'Forest_Cover', 'LF2022_CC_220_CONUS', 'Tif', 'LC22_CC_220.tif'),
        "mnResidueMap" : mnResidueMap,
        "geeResidueMap" : geeResidueMap,
##        "mnTillageTable" : opj(locationsDictNoVersion["huc8_fields"], 'huc' + huc8 + '_mn_rc' + ACPFyear),
##        "geeTillageTable" : opj(locationsDictNoVersion["huc8_fields"], 'huc' + huc8 + '_gee_rc' + ACPFyear),
        "tillageTable" : opj(ACPFDir, 'huc' + huc12 + '_till' + ACPFyear),
        # "geeTillageTable" : opj(ACPFDir, 'huc' + huc12 + '_till' + ACPFyear),
        "mnResidueTable" : opj(ACPFDir, 'huc' + huc12 + '_mn_rc' + ACPFyear),
        "geeResidueTable" : opj(ACPFDir, 'huc' + huc12 + '_gee_rc' + ACPFyear),
        "eptBaseDir" : opj(locationsDict["eleBaseDir"], 'ept'),

        "ept_wesm_monthly_file" : opj(locationsDict["eleBaseDir"], 'ept', 'ept.gdb', ept_first_of_month_name),
        "sampleOutput" : opj(depOutput, '_'.join(['smpl', 'acpf' + ACPFyear, nowYmd]), os.path.basename(locationsDict['samples']) + '.json'),
        "sampleDefOutput" : opj(depOutput, '_'.join(['smpldef', 'acpf' + ACPFyear, nowYmd]), os.path.basename(locationsDict['samples']) + '.json'),
        # nulls are stored with same basename as samples (they are null samples)
        "nullOutput" : opj(depOutput, '_'.join(['smpl', 'acpf' + ACPFyear, nowYmd]), os.path.basename(locationsDict['nulls']) + '.json'),
####        "nullDefOutput" : opj(depOutput, '_'.join(['smpl', 'acpf' + ACPFyear, nowYmd]), os.path.basename(locationsDict['nulls']) + '.json'),
        "fieldBoundariesOutput": opj(depOutput, '_'.join(['fb', 'acpf' + ACPFyear, nowYmd]), os.path.basename(locationsDict['fieldBoundaries']) + '.json'),
        "LU6Output": opj(depOutput, '_'.join(['lu6', 'acpf' + ACPFyear, nowYmd]), os.path.basename(locationsDict['LU6']) + '.csv'),
        "areaOutput" : opj(depOutput, '_'.join(['areas', 'acpf' + ACPFyear, nowYmd]), os.path.basename(locationsDict['areas']) + '.json'),
        # "tillageOutput" : opj(depOutput, '_'.join(['tillage', 'acpf' + ACPFyear, nowYmd]), os.path.basename(locationsDict['tillages']) + '.json'),
        "snapOutput" : opj(depOutput, '_'.join(['snap', 'acpf' + ACPFyear, nowYmd]), os.path.basename(locationsDict['snaps']) + '.json'),


        "pickleDistanceFile" : opj(locationsDict['cutProcDir'], 'search_' + huc12 + '_' + interpType + '.pkl')})

        return locationsDict

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

    
    

def nukedir(dir):
    if dir[-1] == os.sep: dir = dir[:-1]
    files = os.listdir(dir)
    for file in files:
        if file == '.' or file == '..': continue
        path = dir + os.sep + file
        if os.path.isdir(path):
            nukedir(path)
        else:
            os.unlink(path)
    os.rmdir(dir)

def loadHucs(scriptName, fc, HUCList = []):
    '''create a tuple of HUC12 code, EPGS code, and state abbreviation as well as a name of
    a 'stopFile' for which the existence is checked at each iteration of a loop'''

    try:
        
        if len(HUCList) > 0:
            where = buildStringSelection(HUCList, 'HUC12')
            hucSearch = ''
            pickleFileSuffix = 'pickle'
        else:
            hucSearch = os.path.basename(scriptName).split('_')[2]
            stateSearch = os.path.basename(scriptName).split('_')[3][:2]

            pickleFileSuffix = hucSearch + '_' + stateSearch

            if 'all' in hucSearch:
                if 'MW' in stateSearch:
                    where = "States NOT LIKE \'%CN%\' AND HUC12 IS NOT NULL AND HUC12 <> \'\'"
                else:
                    where = "States LIKE '%" + stateSearch + "%'"
            elif 'MW' in stateSearch:
                where = "HUC12 LIKE '" + hucSearch + "%'"
            elif 'MW' not in stateSearch:
                if 'all' not in hucSearch:
                    where = "HUC12 LIKE '" + hucSearch + "%' AND States LIKE '%" + stateSearch + "%' AND HUC12 IS NOT NULL AND HUC12 <> \'\'"
            else:
                where = ""

        fields = ['HUC12', 'PCSCode_By_County', 'States']

        HUCandPCSList = []
        with arcpy.da.SearchCursor(fc, fields, where_clause = where, sql_clause = [None, 'ORDER BY HUC12']) as scur:
            for srow in scur:
                HUCandPCSList.append([srow[0], srow[1], srow[2]])

        if 'build' in scriptName.lower():
            stopFile = 'stop_bulk_pit_' + pickleFileSuffix + '.txt'
        elif 'cut' in scriptName.lower():
            stopFile = 'stop_bulk_cut_' + pickleFileSuffix + '.txt'
        elif 'chnl' in scriptName.lower():
            stopFile = 'stop_bulk_chnl_' + pickleFileSuffix + '.txt'
        elif 'copier' in scriptName.lower():
            stopFile = 'stop_bulk_copier_' + pickleFileSuffix + '.txt'
        elif 'flow' in scriptName.lower():
            stopFile = 'stop_bulk_flow_' + pickleFileSuffix + '.txt'
        elif 'random' in scriptName.lower():
            stopFile = 'stop_bulk_random_' + pickleFileSuffix + '.txt'
        elif 'sample' in scriptName.lower():
            stopFile = 'stop_bulk_sample_' + pickleFileSuffix + '.txt'
        else:
            stopFile = 'stop_bulk_stuff_' + pickleFileSuffix + '.txt'
        print('stopFile is: ' + stopFile)

        # handle start and stride (end could be added)
##        cellSize = 3
        basenameNoext = os.path.splitext(os.path.basename(scriptName))[0]
        if 'tart' in basenameNoext:
            tart = basenameNoext.find('tart')
            start = int(basenameNoext[tart + len('tart'):])
        else:
            start = 0
        shortlist = HUCandPCSList[start:]

        if 'tride' in scriptName:
            tride = scriptName.find('tride')
            afterTride = scriptName[tride:].split('_')[0]
            stride = int(afterTride[len('tride'):])
        else:
            stride = 1
        shortlist = shortlist[::stride]

        if scriptName.lower().find('negstride') > -1:
            shortlist.reverse()

        return shortlist, stopFile, hucSearch, where

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


def setupLoggingNoCh(node, scriptName, huc12 = '000000000000', version = ''):
    # create logger with name 'example'
    log = logging.getLogger('example')
    log.setLevel(logging.DEBUG)

    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(levelname)s - %(message)s')

    nowYmd = datetime.datetime.strftime(datetime.datetime.now(), '%Y_%m_%d_%H_%M_%S')
    logsDir = defineLocalProc(node)
    logName = os.path.join(logsDir, 'Logs', os.path.splitext(os.path.basename(scriptName))[0] + '_' + huc12 + '_' + nowYmd + '.txt')
    if not os.path.isdir(os.path.dirname(logName)):
        os.makedirs(os.path.dirname(logName))

    # create file handler to log debug messages, new log file each time
    fh = logging.FileHandler(logName, mode = 'w')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    log.addHandler(fh)

    startTime = time.time()
    log.info("Beginning logging for script at " + str(time.asctime()))
    log.info("Logging output to: " + logName)

    return log, nowYmd, logName, startTime

def setupLoggingNew(node, scriptName, huc12 = '000000000000', version = ''):
    # create logger with name 'example'
    log = logging.getLogger('example')
    log.setLevel(logging.DEBUG)

    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(levelname)s - %(message)s')
    ##formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    nowYmd = datetime.datetime.strftime(datetime.datetime.now(), '%Y_%m_%d_%H_%M_%S')
    logsDir = defineLocalProc(node)
    logName = os.path.join(logsDir, 'Logs', os.path.splitext(os.path.basename(scriptName))[0] + '_' + huc12 + '_' + nowYmd + '.txt')
    if not os.path.isdir(os.path.dirname(logName)):
        os.makedirs(os.path.dirname(logName))

    # create file handler to log debug messages, new log file each time
    fh = logging.FileHandler(logName, mode = 'w')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    log.addHandler(fh)

    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    log.addHandler(ch)

    startTime = time.time()
    log.info("Beginning logging for script at " + str(time.asctime()))
    log.info("Logging output to: " + logName)
####    log.warn(outputString)

    return log, nowYmd, logName, startTime

def setupLogging(logsDir, scriptName, huc12):
    # create logger with name 'example'
    log = logging.getLogger('example')
    log.setLevel(logging.DEBUG)

    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(levelname)s - %(message)s')
    ##formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    nowYmd = datetime.datetime.strftime(datetime.datetime.now(), '%Y_%m_%d_%H_%M_%S')
    logName = os.path.join(logsDir, 'Logs', os.path.splitext(os.path.basename(scriptName))[0] + '_' + huc12 + '_' + nowYmd + '.txt')
    if not os.path.isdir(os.path.dirname(logName)):
        os.makedirs(os.path.dirname(logName))

    # create file handler to log debug messages, new log file each time
    fh = logging.FileHandler(logName, mode = 'w')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    log.addHandler(fh)

    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    log.addHandler(ch)

    startTime = time.time()
    log.info("Beginning logging for script at " + str(time.asctime()))
    log.info("Logging output to: " + logName)
####    log.warn(outputString)

    return log, nowYmd, logName, startTime
    

def MakeHUClist(inHUCFCls):
	
	HUCList = []
	
	# Make a feature layer of the HUC12 submitted
	arcpy.AddMessage("make features...")
	arcpy.MakeFeatureLayer_management (inHUCFCls, "HUCselect")

	# Select TileIDs
	Rows = arcpy.SearchCursor("HUCselect")
	for row in Rows:
		rowstr = str(row.getValue("HUC12"))
		HUCList.append(rowstr)		

	arcpy.AddMessage("HUC count: " + str(len(HUCList)))
	return(HUCList)

def tryAddField(inFC, fld2Add, fldType):
    flds = getfields(inFC)
    if fld2Add not in flds:
        arcpy.AddField_management(inFC, fld2Add, fldType)

def figureItOut(inputRaster):
    '''figure out huc12 code and projected coordinate system from a formatted raster name'''
    basename = os.path.basename(inputRaster)
    fillbaseNoExt = os.path.splitext(basename)[0]
    huc12 = fillbaseNoExt[-12:]

    # assume name format like 'ef' or 'ep' or 'ec'
    end = fillbaseNoExt[2:]
    # end = fillbaseNoExt.split('ef')[1]
        
    proc_size = int(end.split('m')[0])

    ## pcs is projected cordinate system EPSG code
    # rasterDesc = arcpy.Describe(inputRaster)
    # proc_size = rasterDesc.meanCellHeight
####    pcs = rasterDesc.spatialReference.PCSCode
####    pcs = sys.argv[1].split(os.path.sep)[-3].split('_')[-1]
    huc8 = huc12[:8]

    return huc12, huc8, proc_size

def splitall(path):
##https://www.oreilly.com/library/view/python-cookbook/0596001673/ch04s16.html
    allparts = []
    while 1:
        parts = os.path.split(path)
        if parts[0] == path:  # sentinel for absolute paths
            allparts.insert(0, parts[0])
            break
        elif parts[1] == path: # sentinel for relative paths
            allparts.insert(0, parts[1])
            break
        else:
            path = parts[0]
            allparts.insert(0, parts[1])
    return allparts    
##-------------------------------------------------------------------------
## ---- Modules ----
## from fill_Pits.py
    
			
def condenseStats(inZone, inValue, outTable, statType, sfx):
##    start2 = time.time()
    if sfx[-2:] == "_0" or sfx[-2:] == "_1":
        if arcpy.Exists('in_memory\\' + outTable) == False:
            zst = ZonalStatisticsAsTable(inZone, "VALUE", inValue, 'in_memory\\' + outTable, "", statType)
        else:
            zst1 = ZonalStatisticsAsTable(inZone, "VALUE", inValue, 'in_memory\\' + outTable + sfx, "", statType)
            zst = arcpy.Append_management(zst1, 'in_memory\\' + outTable, "NO_TEST")
            arcpy.Delete_management(zst1)
    else:
        zst1 = ZonalStatisticsAsTable(inZone, "VALUE", inValue, 'in_memory\\' + outTable + sfx, "", statType)
        zst = arcpy.Append_management(zst1, 'in_memory\\' + outTable, "NO_TEST")
        arcpy.Delete_management(zst1)

##    middle = time.time()

    return zst
      
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------

##def ZonalStatistics(in_zone, zone_field, in_value, statistic):
##    if in_zone.maximum < 65535:#32767:
##        zs = arcpy.sa.ZonalStatistics(in_zone, zone_field, in_value, statistic)
##    else:
##        log.debug('using overridden zonal statistics')
##        if statistic in getfields(in_zone):
##            arcpy.DeleteField_management(in_zone, statistic)
##        zst = ZonalStatisticsAsTable(in_zone, zone_field, in_value, gdb + 'agaasdekfr', '', 'ALL')
##        arcpy.JoinField_management(in_zone, zone_field, zst, 'VALUE', statistic)
##        zs = Lookup(in_zone, statistic)
##    return zs


## Cutter processing functions

def testForZero(dataset):
    if isinstance(dataset, arcpy.Raster):
##    if type(dataset) == 'Raster':
        if dataset.maximum is None:
            fcount = 0
        else:
            if dataset.hasRAT != True:
                arcpy.BuildRasterAttributeTable_management(dataset)
            try:
                fcount = int(arcpy.GetCount_management(dataset).getOutput(0))
            except:
                fcount = 0
    else:
        try:
            fcount = int(arcpy.GetCount_management(dataset).getOutput(0))
        except:
            fcount = 0
    return fcount

def ZonalStatistics(in_zone_data, zone_field, in_value_raster, statistics_type='#', ignore_nodata = '#'):
####    if type(in_zone_data) is Raster:
####        if not in_zone_data.hasRAT:# is not True:
####            arcpy.BuildRasterAttributeTable_management(in_zone_data)
####    else:
####        desc = arcpy.Describe(in_zone_data)
####        if 'Raster' in desc.datasetType:
####            arcpy.BuildRasterAttributeTable_management(in_zone_data)

    if testForZero(in_zone_data) == 0:
        raise ZeroFeaturesError(in_zone_data, 'No valid zones for Zonal Statistics')
    else:
        result = arcpy.sa.ZonalStatistics(in_zone_data, zone_field, in_value_raster, statistics_type, ignore_nodata)
        return result

##class NoFeaturesError(Exception):
##    pass

class ZeroFeaturesError(Exception):
    """Raised when arcpy attempts zonal calculations on an empty raster/feature class.
        Often this will cause Python to hard crash (to restart)"""

    def __init__(self, zone, msg):
        self.zone = zone
        # a message to pass out on error
        self.msg = msg
    
    
def tno(fc):
    if 'gdb' in fc.getOutput(0):
        name = fc.getOutput(0).split('.')[1][4:]
    elif 'in_memory' in fc.getOutput(0):
        name = fc.getOutput(0).split('\\')[-1]
    else:
        nameList = fc.getOutput(0).split('.')
        name = nameList[len(nameList)-1]

    return name


def calcTop1Distinct(inTable, fill2cut, cmb_score_fld, fldsort, ws_tmp):
    a = set()

    rank = 1
    prevFr = 0
    prevMatchWs = 0

    matchWs = []

    with arcpy.da.UpdateCursor(inTable, (fill2cut, cmb_score_fld, 'Pair_Rank', ws_tmp.name), sql_clause = (None, 'ORDER BY ' + fill2cut + ' ASC, ' + cmb_score_fld + ' ' + fldsort)) as ucursor:
        for row in ucursor:
            a.add(row[0])
            # starting on new fill region, rank = 1 and matchWs = []
            if  row[0] != prevFr:
                rank = 1
                matchWs = []
                row[2] = rank
                rank += 1
                matchWs.append(row[3])
            # starting on new fill region/ws combo, rank is next
            if row[0] == prevFr and row[3] not in matchWs:
                row[2] = rank
                rank += 1
                matchWs.append(row[3])
            prevFr = row[0]
            ucursor.updateRow(row)

    return sorted(a)

def intersectingFeaturesUniqueIteration4(inputFC, bufferSpec, nearTable, passField, frFld, sortFld):

    try:
    ## frFld is unique identifier for each feature
    ## sortFld defines order in which to process features, such as by largest area or by elevation
    ## passField stores the unique iteration to which each field within bufferSpec is assigned
        tryAddField(inputFC, passField, 'SHORT')
    ## the resulting table ('nearTable') has every 'near feature' the two input FCs (the same in this example) that are within the 'bufferSpec' distance
    ## store 'nearTable' in in_memory location for best performace since we will open a search cursor on in it many times
        nearfcTable = arcpy.GenerateNearTable_analysis(inputFC, inputFC, nearTable, bufferSpec, "NO_LOCATION", "NO_ANGLE", "ALL")

    ## initialize some lists to store features that are being coded for processing in current iteration
    ## cutNowOrNearList stores all features being assigned to the current processing iteration or features within bufferSpec
        ## THIS PROCESS could be SIMPLIFIED if we didn't want to include a buffer distance between features
        cutNowOrNearList = []
    ## cutNowList stores all features being assigned ot the current processing iteration (becomes a list of lists, one list for each iteration)
        cutNowList = []
    ## cutItrtr is the number of times the inputFC will have to be processed to process all unique steps, increase for each updateCursor stepthrough
        cutItrtr = 0
    ## Use modflag to detect when the previous iteration of the update cursor returned no records (previous cursor was empty)
        modflag = 1
        while modflag > 0:
            cutNowList.append([])
            modflag = 0
        ## Find all features not yet assigned a pass number (starting at 0)
            ## pass number updates at end of each time through while loop
            with arcpy.da.UpdateCursor(inputFC, ('OID@', passField, frFld, sortFld), passField + ' IS NULL', sql_clause =(None, 'ORDER BY ' + sortFld + ' ASC')) as ucur1:

            ## Get full list of all inputFC objectids that are not assigned a pass number
                for urow in ucur1:
                    ## Check to see if the feature is in our list of items not to process, if not assign number and add to list for rest of iteration
                    if urow[0] not in cutNowOrNearList:
                        cutNowOrNearList.append(urow[0])
                        cutNowList[cutItrtr].append(urow[2])
                        urow[1] = cutItrtr
                        ucur1.updateRow(urow)
                        modflag += 1

                    ## Find any features that intersect or are near this inputFC and add them to the cutNowOrNearList so they won't get assigned to this pass
                        with arcpy.da.SearchCursor(nearfcTable, ('NEAR_FID', 'IN_FID'), 'IN_FID = ' + str(urow[0])) as scursor:
                            for srow in scursor:
                                cutNowOrNearList.append(srow[0])
            ## Clear for the next pass
            cutNowOrNearList = []
            cutItrtr +=1

    ## Roll back cut iterator by 1 to remove the increment from when the upadte cursor was empty, clear last entry from cutNoweList
        cutItrtr -= 1
        cutNowList.pop()

    except Exception as e:
##        log.debug(e.message)
        arcpy.AddError(e.message)

    finally:
        return cutItrtr#, cutNowList



def createCLDEM(DEM2Mod, gdb, cutFC, outDEMname, sfx, cutElFld, ProcSize):
    try:

##        cutFcGDB = arcpy.CopyFeatures_management(cutFC, opj(gdb, os.path.basename(cutFC)))# + cutFC.getOutput(0).split('\\')[-1])
        cutFcGDB = arcpy.CopyFeatures_management(cutFC, gdb + cutFC.getOutput(0).split('\\')[-1])
        arcpy.AddField_management(cutFcGDB, "Inv_Line_Ord", "LONG")
        arcpy.CalculateField_management(cutFcGDB, "Inv_Line_Ord", str(int(DEM2Mod.maximum))+'-1*!' + cutElFld + '!', "PYTHON")
        cutRaster = arcpy.PolylineToRaster_conversion(cutFcGDB, cutElFld, gdb + "CL_All" + sfx, "", "Inv_Line_Ord", str(ProcSize))

## If the cut line elevation is higher than the original elevation, don't modify the elevation value
        cutLineDif = Minus(cutRaster, DEM2Mod)

        cutLineDifAllCells = Con(IsNull(cutLineDif) == 1, 0, cutLineDif)

        cutLineDEM = Con(cutLineDifAllCells >= 0, DEM2Mod, cutRaster)
        cutLineDEM.save(outDEMname + sfx)

##        log.debug("Finished createCLDEM at " + time.asctime())

        return cutLineDEM, cutRaster

    except Exception as e:
##        log.debug(e.message)
        arcpy.AddError(e.message)

        # Get the traceback object
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]

        # Concatenate information together concerning the error into a message string
        pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
        print(pymsg)
##        log.warn(pymsg)


def createCLDEMnoCopy(DEM2Mod, gdb, cutFcGDB, outDEMname, sfx, cutElFld, ProcSize):
    try:

##        cutFcGDB = arcpy.CopyFeatures_management(cutFC, gdb + cutFC.getOutput(0).split('\\')[-1])
        arcpy.AddField_management(cutFcGDB, "Inv_Line_Ord", "LONG")
        arcpy.CalculateField_management(cutFcGDB, "Inv_Line_Ord", str(int(DEM2Mod.maximum))+'-1*!' + cutElFld + '!', "PYTHON")
        cutRaster = arcpy.PolylineToRaster_conversion(cutFcGDB, cutElFld, gdb + "CL_All" + sfx, "", "Inv_Line_Ord", str(ProcSize))

## If the cut line elevation is higher than the original elevation, don't modify the elevation value
        cutLineDif = Minus(cutRaster, DEM2Mod)

        cutLineDifAllCells = Con(IsNull(cutLineDif) == 1, 0, cutLineDif)

        cutLineDEM = Con(cutLineDifAllCells >= 0, DEM2Mod, cutRaster)
        cutLineDEM.save(outDEMname + sfx)

##        log.debug("Finished createCLDEM at " + time.asctime())

        return cutLineDEM, cutRaster

    except Exception as e:
##        log.debug(e.message)
        arcpy.AddError(e.message)

        # Get the traceback object
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]

        # Concatenate information together concerning the error into a message string
        pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
        print(pymsg)
##        log.warn(pymsg)




def addCalcJoin(targetTbl, targetField, joinTbl, joinField, newFieldList, oldField):
    tryAddField(joinTbl, newFieldList[0], newFieldList[1])
    arcpy.CalculateField_management(joinTbl, newFieldList[0], oldField, 'python')
    arcpy.JoinField_management(targetTbl, targetField, joinTbl, joinField, newFieldList[0])


def getfields(infc, fieldString = ''):
    fields = [fld.name for fld in arcpy.ListFields(infc, fieldString)]
##    del fld
    return fields

        
def buildSelection(inList, field):
    sel = ''
    for index, item in enumerate(inList):
        if index == 0:
            sel = field + ' = ' + str(item)
        else:
            sel += ' OR ' + field + ' = ' + str(item)
    return sel
        
def buildAndSelection(inList, field, andSel):
    sel = ''
    for index, item in enumerate(inList):
        if index == 0:
            sel = field + ' = ' + str(item) + andSel
        else:
            sel += ' OR ' + field + ' = ' + str(item) + andSel
    return sel
        
def buildStringSelection(inList, field):
	for index, item in enumerate(inList):
		if index == 0:
			sel = field + ' = \'' + str(item) + '\''
		else:
			sel += ' OR ' + field + ' = \'' + str(item) + '\''
	return sel

def buildAntiSelection(inList, field):
	for index, item in enumerate(inList):
		if index == 0:
			sel = field + ' <> ' + str(item)
		else:
			sel += ' AND ' + field + ' <> ' + str(item)
	return sel


def selectByList(in_list, in_field, in_tbl, out_tbl, sfx, inm):
    list_length = len(in_list)
    list_threshold = 14999
    if list_length > list_threshold:
        iter_needed = list_length//list_threshold + 1
        append_list = []
        slice_denom = int(len(in_list)/iter_needed)
        for i in range(iter_needed):
            ofSfx = '_' + str(i)
            sub_list = in_list[i*slice_denom:(i+1)*slice_denom]
            if i == iter_needed - 1:
                sub_list = in_list[i*slice_denom:]
            if len(sub_list) > 1:
                gc_sel4step = in_field + " IN " + str(tuple(sub_list))
            else:
                gc_sel4step = buildSelection(sub_list, in_field)
            if i == 1:
                sel_init = arcpy.Select_analysis(in_tbl, opj(inm, out_tbl + sfx), gc_sel4step)
            else:
                sel_later = arcpy.Select_analysis(in_tbl, opj(inm, out_tbl + sfx + ofSfx), gc_sel4step)
                append_list.append(sel_later)
        sel_out = arcpy.Append_management(append_list, sel_init)
    else:
        if len(in_list) > 1:
            gc_sel4step = in_field + " IN " + str(tuple(in_list))
##            gc_sel4step = in_field + " IN " + str(tuple(sub_list))
        else:
            gc_sel4step = buildSelection(in_list, in_field)
        sel_out = arcpy.Select_analysis(in_tbl, opj(inm, out_tbl + sfx), gc_sel4step)
    return sel_out
        
##    stringLength = len(in_list) * (len(in_field) + len(' = ') + len(str(in_list[-1])) + len(' OR '))
##    if stringLength > 150000:
##        iterNeeded = stringLength//150000 + 1
##        appendList = []
##        sliceDenom = int(len(in_list)/iterNeeded)
##        for i in range(iterNeeded):
##            ofSfx = '_' + str(i)
##            sublist = in_list[i*sliceDenom:(i+1)*sliceDenom]
##            if i == iterNeeded - 1:
##                sublist = in_list[i*sliceDenom:]
##            
##            gcSel4step = buildSelection(sublist, in_field)
##            if i == 1:
##                selInit = arcpy.Select_analysis(in_tbl, inm + out_tbl + sfx, gcSel4step)
##            else:
##                selLater = arcpy.Select_analysis(in_tbl, inm + out_tbl + sfx + ofSfx, gcSel4step)
##                appendList.append(selLater)
##        selOut = arcpy.Append_management(appendList, selInit)
##    else:
##        gcSel4step = buildSelection(in_list, in_field)
##        selOut = arcpy.Select_analysis(in_tbl, inm + out_tbl + sfx, gcSel4step)
##    return selOut


def antiSelectByList(in_list, in_field, in_tbl, out_tbl, sfx, inm):
    stringLength = len(in_list) * (len(in_field) + len(' <> ') + len(str(in_list[-1])) + len(' OR '))
    if stringLength > 150000:
        iterNeeded = stringLength//150000 + 1
        appendList = []
        sliceDenom = int(len(in_list)/iterNeeded)
        for i in range(iterNeeded):
            ofSfx = '_' + str(i)
            sublist = in_list[i*sliceDenom:(i+1)*sliceDenom]
            if i == iterNeeded - 1:
                sublist = in_list[i*sliceDenom:]
            
            gcSel4step = buildAntiSelection(sublist, in_field)
            if i == 1:
                selInit = arcpy.Select_analysis(in_tbl, inm + out_tbl + sfx, gcSel4step)
            else:
                selLater = arcpy.Select_analysis(in_tbl, inm + out_tbl + sfx + ofSfx, gcSel4step)
                appendList.append(selLater)
        antiSelOut = arcpy.Append_management(appendList, selInit)
    else:
        gcSel4step = buildAntiSelection(in_list, in_field)
        antiSelOut = arcpy.Select_analysis(in_tbl, inm + out_tbl + sfx, gcSel4step)
    return antiSelOut


def conByList(in_list, in_field, in_rast, cp):
    initWs = arcpy.env.workspace
    arcpy.env.workspace = cp

    list_length = len(in_list)
        
    if list_length == 1:
        gcSel4step = buildSelection(in_list, in_field)
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

##    stringLength = len(in_list) * (len(in_field) + len(' = ') + len(str(in_list[-1])) + len(' OR '))
##    if stringLength > 150000:
####        log.warn('using conByList with string Length = ' +str(stringLength))
##        iterNeeded = stringLength//150000 + 1
##        appendList = []
##        sliceDenom = int(len(in_list)/iterNeeded)
##        for i in range(iterNeeded):
##            ofSfx = '_' + str(i)
##            sublist = in_list[i*sliceDenom:(i+1)*sliceDenom]
##            if i == iterNeeded - 1:
##                sublist = in_list[i*sliceDenom:]
##            
##            gcSel4step = buildSelection(sublist, in_field)
##            con4Later = Con(in_rast, in_rast, '', gcSel4step)
##            con4Later.save('c4l_' + str(i))
##            appendList.append(con4Later.name)
####            appendList.append(str(con4Later).split('\\')[-1])
####        log.debug(appendList)
##        selOutFloat = CellStatistics(appendList)
##        selOut = Int(selOutFloat)
####        for i in appendList:
####            arcpy.Delete_management(i)
##    else:
##        gcSel4step = buildSelection(in_list, in_field)
##        selOut = Con(in_rast, in_rast, '', gcSel4step)
##
##    arcpy.env.workspace = initWs
##
##    return selOut



def copyfc(verbose, infc, gdb):
    if verbose == True:
##        out = arcpy.CopyFeatures_management(infc, os.path.join(gdb, os.path.basename(str(infc)) + infc.getOutput(0).split('\\')[-1])
##        out = arcpy.CopyFeatures_management(infc, opj(gdb, infc.getOutput(0).split('\\')[-1]))
        out = arcpy.CopyFeatures_management(infc, opj(gdb, os.path.basename(str(infc))))
    else:
        out = None
    return out

def copyfc2(infc, gdb):
    out = arcpy.CopyFeatures_management(infc, gdb + infc.getOutput(0).split('\\')[-1])
    return out

def copytbl(verbose, intbl, gdb):
    if verbose == True:
        out = arcpy.CopyRows_management(intbl, gdb + intbl.getOutput(0).split('\\')[-1])
    else:
        out = None
    return out

def copytbl2(intbl, gdb):
    out = arcpy.CopyRows_management(intbl, gdb + intbl.getOutput(0).split('\\')[-1])
    return out

def condDelete(verbose, item):
    if not verbose:
        arcpy.Delete_management(item)
        del item
   

def getFrsAsList(fc, idFld, sel):
    selList = set()
    with arcpy.da.SearchCursor(fc, idFld, sel) as scur:
        for srow in scur:
            selList.add(srow[0])
    return list(selList)


def CreateInitialWs(inDEM):
        fd_tmp = FlowDirection(inDEM)
        snk_tmp = Sink(fd_tmp)
        initialSnkNull0 = Con(IsNull(snk_tmp), 0, snk_tmp)
##        initialSnkNull0.save(cp + 'initial_sink')
        initialWS = Watershed(fd_tmp, snk_tmp)
        arcpy.BuildRasterAttributeTable_management(initialWS)

        return initialWS, initialSnkNull0, fd_tmp


## from djHolePuncher.py
def findMaxDepth(inDEM):
        ## Determine how deep the areas that need fill are
        rFill = Fill(inDEM)

        rDepth = Minus(rFill, inDEM)
     
        MaxDepth = rDepth.maximum
##        p2f('Current maximum depth: ' + str(MaxDepth))
##        arcpy.AddMessage("Current maximum depth: " + str(MaxDepth))

        return(rDepth, MaxDepth, rFill)


## This version still has difficulties with 'connecting' fill regions that aren't really connected (just connected by fill depth zero)
## It is possible that some depressions greater than threshhold, could not be fill regions initially but then appear if they
##  are connected by zero depth fill to a deeper depression but, once filled, actually drained to a different fill region
## The alternative, filling up initial watersheds, is more messy since we don't know which sink is the deepest in a fill region until you're done
def punchHolesByDepth4(inDEM, rDepth, minBsnDpth, minBsnArea, ProcSize, index, cp, initialPunchWs, filledDEM, initialPunchSnkNull0, gdb, verbose, version):
##        log.debug('beginning stack for index ' + str(index) + ' at ' + time.asctime())
        sfx = "_" + str(index)
    #### Punch holes in the bottom of sinks that are deeper or fill regions larger than the threshhold
    ####  - inDEM: the current verion of the surface
    ####  - rDepth: the current depth of fill
    ####  - minBsnDepth: the depth of allowed sinks (i.e. all fill regions deeper than this are punched)
    ####  - minBsnArea: the area of allowed sinks (i.e. all fill regions larger than this are punched)
##        log.debug("Punch holes for " + str(sfx))
## Define fill regions (and cost surface)

        wsMinFilledEl = ZonalStatistics(initialPunchWs, 'VALUE', filledDEM, 'MINIMUM')
        fillDif = filledDEM - wsMinFilledEl
        fillGE0 = Con(fillDif == 0, wsMinFilledEl)
        rgMinFilledEl = RegionGroup(fillGE0, 'EIGHT')
        rgMinFilledEl.save('rg2raw' + sfx)


## Calculate depth and thickness statistics for fill regions
        zstMaxDepth = ZonalStatisticsAsTable(rgMinFilledEl, "VALUE", rDepth, 'in_memory\\zst_max_dpth')
        copytbl(verbose, zstMaxDepth, gdb)
        joinMax = arcpy.JoinField_management(rgMinFilledEl, "VALUE", zstMaxDepth, "VALUE", "MAX")

## Select the appropriate fill regions for enforcement (deep enough or large enough)
        if version.find('10.5') > -1 or version.find('10.6') > -1:
            fldName = arcpy.ListFields(rgMinFilledEl, 'MAX*')[0].name
            rgToPunch = Con(rgMinFilledEl, rgMinFilledEl, "", '"' + fldName + '" > ' + str(minBsnDpth) + ' OR "COUNT" * ' + str(ProcSize**2) + ' > ' + str(minBsnArea))
        else:
            rgToPunch = Con(rgMinFilledEl, rgMinFilledEl, "", '"MAX" > ' + str(minBsnDpth) + ' OR "COUNT" * ' + str(ProcSize**2) + ' > ' + str(minBsnArea))
        arcpy.BuildRasterAttributeTable_management(rgToPunch)
        rgToPunch.save('rg2punch' + sfx)


    ## Define sink and watershed number by minimum elevation sink (minimum value  if tie)
        sinkEl = Con(initialPunchSnkNull0, inDEM)
        rgMinSinkEl = ZonalStatistics(rgToPunch, 'VALUE', sinkEl, 'MINIMUM')
        sinksAtMin = Con(sinkEl == rgMinSinkEl, Con(initialPunchSnkNull0 > 0, initialPunchSnkNull0))
        rgMinSinkVal = ZonalStatistics(rgToPunch, 'VALUE', sinksAtMin, 'MINIMUM')

        holes2Punch = Con(initialPunchSnkNull0 == rgMinSinkVal, initialPunchSnkNull0)
        holes2Punch.save(cp + 'snkunq' + sfx)


        demWithHoles = Con(IsNull(holes2Punch), inDEM, '')
        demWithHoles.save(cp + "nwdm" + sfx)

                            ## Define the watersheds for that region
        fillLvl = Fill(demWithHoles)

        if index > 0:
            snkPrev = Con(IsNull(inDEM), initialPunchSnkNull0)
            cumSinks = CellStatistics([holes2Punch, snkPrev], 'MINIMUM')#[snkUnique, snkPrev], 'MINIMUM')
        else:
            cumSinks = holes2Punch#snkUnique

        fdTemp = FlowDirection(fillLvl)

        wsLvl = Watershed(fdTemp, cumSinks)

    ## To remove minor fill regions (fill > 0, less than criteria, and connected by fill depth >= 0 to other fill regions)
    ## make sure values are the same between watershed and rg min sink value
        fullFr0 = Con(inDEM <= wsMinFilledEl, rgMinSinkVal)
        fullFr0.save(cp + 'fr01' + sfx)

    ## Code as fill region where fill region and watershed number agree, leave blank otherwise (should be shallowish FRs connected by zero fill depth, unlikely to be deep)
        fr0 = Con(fullFr0 == wsLvl, fullFr0)

        return demWithHoles, holes2Punch, fr0, fillLvl, fdTemp, wsLvl, cumSinks



def condenseDataLvls(listToAppend, wsAndOutFc):
##    log.debug('condenseDataLvls fcs list is ' + str(listToAppend))
    if len(listToAppend) > 0:
        first = listToAppend.pop()
##        log.debug('getfields for ' + str(first) + ': ' + str(getfields(first)))
        outFc = arcpy.CopyFeatures_management(first, wsAndOutFc)
    if len(listToAppend) > 0:
##        for fc in listToAppend:
##            log.debug('getfields for ' + str(fc) + ': ' + str(getfields(fc)))
        arcpy.Append_management(listToAppend, outFc)

    for fc in listToAppend:
        arcpy.Delete_management(fc)
        del fc

    return outFc

##def condenseDataLvls(inFc, inWs, wsAndOutFc):
##    prevWs = arcpy.env.workspace
##    arcpy.env.workspace = inWs
##    log.debug('condenseDataLvls searching in ' + str(inWs))
##    searchTerm = '_'.join(inFc.getOutput(0).split('\\')[-1].split('_')[:-1]) + '*'
##    log.debug('condenseDataLvls searchTerm is ' + str(searchTerm))
##    fcs = arcpy.ListFeatureClasses(searchTerm)
##    log.debug('condenseDataLvls fcs list is ' + str(fcs))
##    if len(fcs) > 0:
##        first = fcs.pop()
##        log.debug('getfields for ' + str(first) + ': ' + str(getfields(first)))
##        outFc = arcpy.CopyFeatures_management(first, wsAndOutFc)
##    if len(fcs) > 0:
##        for fc in fcs:
##            log.debug('getfields for ' + str(fc) + ': ' + str(getfields(fc)))
##        arcpy.Append_management(fcs, outFc)
##
##    arcpy.env.workspace = prevWs
##
##    return outFc

def condenseTableLvls(inFc, inWs, wsAndOutFc):
    prevWs = arcpy.env.workspace
    arcpy.env.workspace = inWs
    fcs = arcpy.ListTables('_'.join(inFc.getOutput(0).split('\\')[-1].split('_')[:-1]) + '*')
    if len(fcs) > 0:
        first = fcs.pop()
        outFc = arcpy.CopyRows_management(first, wsAndOutFc)
    if len(fcs) > 0:
        arcpy.Append_management(fcs, outFc)

    arcpy.env.workspace = prevWs

    return outFc



def angleDif(angle1, angle2):
    difference = angle2-angle1
    while difference < -180: difference += 360
    while difference > 180: difference -= 360
    return difference
    

# def addRasterToGroup(CutProc, in_raster, frame, group, order = 'AUTO_ARRANGE'):
#     rstrLayer = arcpy.MakeRasterLayer_management(in_raster, in_raster.name + '_lyr')
#     rstrLayerFile = arcpy.SaveToLayerFile_management(rstrLayer, os.path.join(CutProc, in_raster.name + '.lyr'))
#     rstrLFObj = mp.Layer(rstrLayerFile.getOutput(0))
#     mp.AddLayerToGroup(df, group, rstrLFObj, order)

# def addRasterToMap(CutProc, in_raster, frame, order = 'AUTO_ARRANGE'):
#     rstrLayer = arcpy.MakeRasterLayer_management(in_raster, in_raster.name + '_lyr')
#     rstrLayerFile = arcpy.SaveToLayerFile_management(rstrLayer, os.path.join(CutProc, in_raster.name + '.lyr'))
#     rstrLFObj = mp.Layer(rstrLayerFile.getOutput(0))
#     mp.AddLayer(df, rstrLFObj, order)

# def addFeatureToMap(CutProc, in_fc, frame, order = 'AUTO_ARRANGE'):
#     featLayer = arcpy.MakeFeatureLayer_management(in_fc, in_fc.getOutput(0).split('\\')[-1] + '_lyr')
#     featLayerFile = arcpy.SaveToLayerFile_management(featLayer, os.path.join(CutProc, in_fc.getOutput(0).split('\\')[-1] + '.lyr'))
#     featLFObj = mp.Layer(featLayerFile.getOutput(0))
#     mp.AddLayer(df, featLFObj, order)

def curvatureBasedStreams(crvThresh, crvName, newProCrv, cutRaster, moreThanHalfwayBigFa, highFa, cp, fillAll, fdAll, maxNoCrvLength, ProcSize, lengthCrit, fld, potentialCrvStretch, ndTF, verbose, gdb):
##    crvCostPrelim1 = Expand(Con(newProCrv > crvThresh, 1, 0) + Con(IsNull(cutRaster), 0, 1), 1, 1)
    crvCostPrelim0 = Expand(Con(RegionGroup(Con(Con(newProCrv > crvThresh, 1, 0) + Con(ndTF == 1, 1, 0) + Con(IsNull(cutRaster), 0, 1) >= 1, 1), 'EIGHT'), 1, '', "COUNT > 2"), 1, 1)
    crvCostPrelim1 = Con(IsNull(crvCostPrelim0), 0, crvCostPrelim0)

    crvCostPrelim2 = crvCostPrelim1 + moreThanHalfwayBigFa
    crvCostPre = Con(crvCostPrelim2 >= 1, 1, 2)#1, crvCostPrelim)
    if verbose:
        crvCostPre.save(cp + 'crv_cp' + crvName)

## remove No Crv gaps greater than criteria (arbitrary)
    rgCrvCost = RegionGroup(Con(potentialCrvStretch, crvCostPre), 'EIGHT')
    goodCrvStretch = Con(rgCrvCost, crvCostPre, '', 'LINK = 1 OR LINK = 2 AND COUNT < ' + str(int(maxNoCrvLength/ProcSize)))

## cost distance via high profile curvature from bottomw of WS
    distFromCrv = CostDistance(highFa, Con(goodCrvStretch, crvCostPre))#, 1000)
    if verbose:
        distFromCrv.save(cp + 'dst_' + crvName)

## extract area where flow length is more than 98% in high profile curvature areas
    mostlyCrv05 = Con(distFromCrv < fld*1.02, distFromCrv)
    if verbose:
        mostlyCrv05.save(cp + 'mst_' + crvName)

## stream has to start at long high curvature point (GT criteria)
    rgGoodStart = Con(rgCrvCost, rgCrvCost, '', 'LINK = 1 AND COUNT > ' + str(int(maxNoCrvLength/ProcSize)))
    goodStartInMostly = Con(mostlyCrv05, rgGoodStart)
    connectedDsBig = CostPath(goodStartInMostly, fillAll, fdAll)

    prunedFlowpaths = pruneStreamByLength(connectedDsBig, fdAll, lengthCrit, ProcSize, fillAll)
    stfCrv = StreamToFeature(prunedFlowpaths, fdAll, gdb + 'stream_from_crv_' + crvName)
    stfCrv2 = StreamToFeature(connectedDsBig, fdAll, gdb + 'stream_from_crv_orgnl_' + crvName)
    stfCrv3 = StreamToFeature(Int(distFromCrv), fdAll, gdb + 'stream_from_crv_orgnl3_' + crvName)

##    log.debug('curve ' + str(crvThresh) + ' streams done')

    return stfCrv, distFromCrv


def pruneStreamByLength(distFromCrv, fdAll, lengthCrit, ProcSize, fillAll, cp, verbose):
    ## Calculate distance to upstream end
    fluCrv = FlowLength(Con(distFromCrv >= 0, fdAll), 'UPSTREAM')
    ## Separate out those segments that are less than dangle criteria + cell size (these are upstream ends)
    fluCrvLow1s = Con(fluCrv < (lengthCrit + ProcSize), 1)
    fluCrvLow1sSl = StreamLink(fluCrvLow1s, fdAll)

    ## Determine if upstream ends are long enough to save at all
    ##  MAY HAVE PROBLEMS WITH PARALLEL FLOWPATH ENDS
    rgFluCrvLow1s = RegionGroup(fluCrvLow1s, 'EIGHT')
    maxFluRg = ZonalStatistics(rgFluCrvLow1s, 'value', fluCrv, 'maximum')

    longEnoughSl = Con(maxFluRg > (lengthCrit - ProcSize), fluCrvLow1sSl)

##    cmbRgSl = Combine([rgFluCrvLow1s, fluCrvLow1sSl])
##    maxFluCmb = ZonalStatistics(cmbRgSl, 'value', fluCrv, 'maximum')

    ## Determine if the upstream segments are branching or not by calculating stream link variety
    fsSlVariety = FocalStatistics(longEnoughSl, "RECTANGLE 3 3 CELL", 'VARIETY')
    fsSlVar3Exp = Expand(fsSlVariety, 1, 3)

    ## Trim out those that don't 3 stream links (no branching) and are less than dangle criteria
    zsMaxVarSl = ZonalStatistics(longEnoughSl, 'VALUE', fsSlVar3Exp, 'MAXIMUM')
##    zstMaxVarSl = ZonalStatisticsAsTable(fluCrvLow1sSl, 'VALUE', fsSlVar3Exp, gdb + 'zst_max_var_sl', '','ALL')
    var1Sls = Con(zsMaxVarSl == 1, longEnoughSl)
    var1MaxFlu = ZonalStatistics(var1Sls, 'VALUE', fluCrv, 'MAXIMUM')
    ## no confluence flowpaths longer than criteria
    highFluVar1 = Con(var1MaxFlu > (lengthCrit - 2* ProcSize), var1Sls)

    ## Select long ends with confluences and find lowest flow length upstream and remove that branch
    if fsSlVar3Exp.maximum >= 3.0:
        rgVar3 = RegionGroup(Con(fsSlVar3Exp == 3, 1), 'EIGHT')
        ## Select those with confluences
        var3Sls = Con(rgVar3, longEnoughSl)
        ## Find lowest FLU in each confluence
        minFluRgVar3 = ZonalStatistics(rgVar3, 'VALUE', fluCrv, 'MINIMUM')
        ## Find lowest FLU for each branch in each confluence
        minFluSlVar3 = ZonalStatistics(var3Sls, 'VALUE', fluCrv, 'MINIMUM')
        ## Select branches higher than minimum
        slFluHigherThanMin = Con(minFluSlVar3 > minFluRgVar3, var3Sls)
        var3SlFluPopulated = ZonalStatistics(fluCrvLow1sSl, 'VALUE', slFluHigherThanMin, 'MINIMUM')

        ## Join no confluence starting points and confluence starting points
        cellsPopulated = CellStatistics([var3SlFluPopulated, highFluVar1], 'MAXIMUM')
    else:
        cellsPopulated = highFluVar1

    ## Calculate downstream route
    if verbose:
        cellsPopulated.save(cp + 'cells_pop')
    prunedFlowpathsPre = CostPath(cellsPopulated, fillAll, fdAll)
    prunedFlowpaths = Con(prunedFlowpathsPre, 1)

    return prunedFlowpaths

def overflowFRs(wsStorCrit, DEMwCuts, prevSinks, zsWsStor, zsFilDif, index, verbose):

    try:
        while zsWsStor.minimum < wsStorCrit and zsWsStor.minimum is not None:
##            log.debug('starting next loop, index was ' + str(index) + ' and zsWsStor min was ' + str(zsWsStor.minimum))
            index += 1
            sfx = '_' + str(index)
            
            testflowSnks = Con(zsWsStor > wsStorCrit, prevSinks)
            if verbose:
                testflowSnks.save('tst_snk' + sfx)

        ## figure out what area drains to the sinks being tested                                    
            crvDEM = Con(Con(IsNull(testflowSnks), 0, 1) == 0, DEMwCuts, '')
            filCrvDEM = Fill(crvDEM)
            testflowWss = Watershed(FlowDirection(filCrvDEM), testflowSnks)
        ## What is this in mean storage depth?
            zsWsStor = ZonalStatistics(testflowWss, 'value', zsFilDif, 'mean')
            if verbose:
                testflowWss.save('tst_flw_ws' + sfx)
                zsWsStor.save('zs_ws' + sfx)

            prevSinks = testflowSnks

##        log.debug('after last loop index was ' + str(index) + ' and zsWsStor min was ' + str(zsWsStor.minimum))

        if zsWsStor.minimum == None:
            testflowWss = None

        return prevSinks, zsWsStor, testflowWss, index

    except Exception as err:
##        log.debug(err.message)
        arcpy.AddError(err.message)

        

def stationLines10(LC, infc, ptFC, lineFC, spacing, length, fc4sr, ptFCstart):
    ## LC is a fgdb workspace
    ## infc is the feature class to process (create cross section lines from)
    ## ptFC is the point feature class to be created
    ## lineFC is the line feature class to be created
    ## spacing is the linear distance between points (and lines) (in meters)
    ## length is the lenght of the cross section lines (in meters)
    ## fc4sr is the feature class to use to define spatial reference
## New Point Feature
    if arcpy.Exists(ptFC) == True:
        arcpy.Delete_management(ptFC)
##        log.debug("Deleted NewFeature " + ptFC)

    newFCPts = arcpy.CreateFeatureclass_management(LC, ptFC, "Point", "", "", "", fc4sr)

    newFCPtstart = arcpy.CreateFeatureclass_management(LC, ptFCstart, "Point", "", "", "", fc4sr)

## New Line Feature

    if arcpy.Exists(lineFC) == True:
        arcpy.Delete_management(lineFC)

##    newFCLine = arcpy.CreateFeatureclass_management('in_memory', lineFC, "Polyline", "", "", "", fc4sr)
    newFCLine = arcpy.CreateFeatureclass_management(LC, lineFC, "Polyline", "", "", "", fc4sr)

    newFCLineL = arcpy.CreateFeatureclass_management(LC, lineFC + '_l', "Polyline", "", "", "", fc4sr)

    newFCLineR = arcpy.CreateFeatureclass_management(LC, lineFC + '_r', "Polyline", "", "", "", fc4sr)

    newFCLineSin = arcpy.CreateFeatureclass_management(LC, lineFC + '_sin', "Polyline", "", "", "", fc4sr)

    arcpy.AddField_management(newFCPts, 'arcid', 'long')
    arcpy.AddField_management(newFCPts, 'cs_id', 'long')

    arcpy.AddField_management(newFCPtstart, 'arcid', 'long')
    arcpy.AddField_management(newFCPtstart, 'cs_id', 'long')

    arcpy.AddField_management(newFCLine, 'arcid', 'long')
    arcpy.AddField_management(newFCLine, 'cs_id', 'long')

    arcpy.AddField_management(newFCLineL, 'arcid', 'long')
    arcpy.AddField_management(newFCLineL, 'cs_id', 'long')
    arcpy.AddField_management(newFCLineL, 'cs_id_l', 'long')

    arcpy.AddField_management(newFCLineR, 'arcid', 'long')
    arcpy.AddField_management(newFCLineR, 'cs_id', 'long')
    arcpy.AddField_management(newFCLineR, 'cs_id_r', 'long')

    arcpy.AddField_management(newFCLineSin, 'arcid', 'long')

    try:
        edit = arcpy.da.Editor(LC)
    # Edit session is started without an undo/redo stack for versioned data
    #  (for second argument, use False for unversioned data)
        edit.startEditing(False, False)

    # Start an edit operation
        edit.startOperation()
        
        csid = 0

        iCurPts = arcpy.da.InsertCursor(newFCPts, ["SHAPE@", 'arcid', 'cs_id'])
        iCurPtstart = arcpy.da.InsertCursor(newFCPtstart, ["SHAPE@", 'arcid', 'cs_id'])
        iCurLine = arcpy.da.InsertCursor(newFCLine, ["SHAPE@", 'arcid', 'cs_id'])
        iCurLineL = arcpy.da.InsertCursor(newFCLineL, ["SHAPE@", 'arcid', 'cs_id', 'cs_id_l'])
        iCurLineR = arcpy.da.InsertCursor(newFCLineR, ["SHAPE@", 'arcid', 'cs_id', 'cs_id_r'])
        iCurLineS = arcpy.da.InsertCursor(newFCLineSin, ["SHAPE@", 'arcid'])

        with arcpy.da.SearchCursor(infc, ['SHAPE@', 'HUC_FID']) as scur:
##        with arcpy.da.SearchCursor(infc, ['SHAPE@', 'arcid']) as scur:
            for srow in scur:
                partnum = 0
                partcount = srow[0].partCount

                while partnum < partcount:
                    part = srow[0].getPart(partnum)
                    pnt = part.next()
                    pntcount = 0

                    ## statnPtsSpacing is the distance between Station Points along the line
                    statnPtsSpacing = spacing
                    ## statnSctnLeftovers is the unallocated distance between Station Points along the line when you move to a new part
                    statnSctnLeftovers = 0.00
                    ## statnLineLength is the width of the perpindicular line at the station point
                    statnLineLength = length

                    # Enter while loop for each vertex
                    while pnt:
                        # Print x,y coordinates of current point
                        if pntcount > 0:

                            ## Calculate the distance (hypotenuse) and angle between previous and current point
                            pnt_dx = pnt.X - prevPtX
                            pnt_dy = pnt.Y - prevPtY
                            hypotenuse = math.sqrt(pnt_dx**2 + pnt_dy**2)
                            angle = math.atan2(pnt_dx, pnt_dy)
                            perp = angle + math.pi

                            try:
                                eastness = pnt_dx/hypotenuse
                                northness = pnt_dy/hypotenuse
                            except:
                                eastness = 0.0
                                northness = 0.0

                            statnPtsInSctn = (hypotenuse-(statnSctnLeftovers * statnPtsSpacing))/statnPtsSpacing + 1

                            ## Calculate statnPtsSpacing distance between starting point and next point
                            statnPtDX = statnPtsSpacing*eastness
                            statnPtDY = statnPtsSpacing*northness

                            statnLineLeftDX = statnLineLength*-northness
                            statnLineLeftDY = statnLineLength*eastness
                                
        ## Calculate for all sections
                            for i in range(0, int(statnPtsInSctn)):
            ## Set point X & Y
                                newPt = arcpy.Point(prevPtX + (statnSctnLeftovers + i)*statnPtDX, prevPtY + (statnSctnLeftovers + i)*statnPtDY)
                                iCurPts.insertRow((newPt, srow[1], csid))

                                newPt4Line1 = arcpy.Point(newPt.X + statnLineLeftDX, newPt.Y + statnLineLeftDY)
                                iCurPtstart.insertRow((newPt4Line1, srow[1], csid))
                                
                                newPt4Line2 = arcpy.Point(newPt.X - statnLineLeftDX, newPt.Y - statnLineLeftDY)
                                array = arcpy.Array([newPt4Line1, newPt4Line2])
                                polyline = arcpy.Polyline(array)

                                iCurLine.insertRow((polyline, srow[1], csid))

                                arrayL = arcpy.Array([newPt, newPt4Line2])
                                polylineL = arcpy.Polyline(arrayL)
                                csidl = csid *10

                                iCurLineL.insertRow((polylineL, srow[1], csid, csidl))

                                arrayR = arcpy.Array([newPt, newPt4Line1])
                                polylineR = arcpy.Polyline(arrayR)
                                csidr = csid *10 + 1

                                iCurLineR.insertRow((polylineR, srow[1], csid, csidr))

                                if pntcount == 1 and i == 0:
                                    arrayS = arcpy.Array([prevNewPt, newPt])
                                else:
                                    arrayS.append(newPt)

                                csid += 1

                                prevNewPt = newPt

                            ## Calculate how much distance was left between last station point in line part and next point in line feature
                            statnSctnRemainingFraction = (hypotenuse-(statnSctnLeftovers * statnPtsSpacing))%statnPtsSpacing/statnPtsSpacing
                            ## Calculate how much distance beetween station points was left over
                            statnSctnLeftovers = 1-statnSctnRemainingFraction 

                        if pntcount == 0:
                            prevNewPt = pnt
##                        else:
##                            prevNewPt = newPt


                        prevPtX = pnt.X
                        prevPtY = pnt.Y
                           
                        pnt = part.next()

                        pntcount += 1

                    partnum += 1

                # Add the line segment to the last point
                lastPoint = arcpy.Point(prevPtX, prevPtY)
                arrayS.append(lastPoint)
        ######        arrayS = arcpy.Array([newPt, lastPoint])
                polylineS = arcpy.Polyline(arrayS)

                iCurLineS.insertRow((polylineS, srow[1]))

        # Stop the edit operation.
        edit.stopOperation()

        # Stop the edit session and save the changes
        edit.stopEditing(True)

    except Exception as err:
##        log.debug(err.message)
        arcpy.AddError(err.message)

        # Get the traceback object
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]

        # Concatenate information together concerning the error into a message string
        pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
        msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"

        # Return python error messages for use in script tool or Python Window
        arcpy.AddError(pymsg)
        arcpy.AddError(msgs)

##        # Print Python error messages for use in Python / Python Window
##        log.warn(pymsg + "\n")
##        log.warn(msgs)

    finally:
        del iCurLine
        del iCurLineL
        del iCurLineR
        del iCurLineS
        del iCurPts
        del iCurPtstart
        del edit

    return newFCPts, newFCLine, newFCLineL, newFCLineR, newFCLineSin, newFCPtstart


## creates a Terrain Position Index (TPI) raster using a 3 cell annulus mean
def makeTPI(DEM):
    annulus1x3Nbr = NbrAnnulus(1, 3, "CELL")
    focalMeanAnnulus = FocalStatistics(DEM, annulus1x3Nbr, 'MEAN')
    tpi = DEM - focalMeanAnnulus

    return tpi
    
def expandByFst(zones, distance):
	dist = distance*2 + 1
	exp1ThinMax = FocalStatistics(zones, 'RECTANGLE ' + str(dist) + ' ' + str(dist) + ' CELL', 'MAXIMUM')
	exp1ThinMaj = FocalStatistics(zones, 'RECTANGLE ' + str(dist) + ' ' + str(dist) + ' CELL', 'MAJORITY')
	exp1Thin = Con(IsNull(exp1ThinMaj), exp1ThinMax, exp1ThinMaj)
	return exp1Thin

def createInvertDEM(inDEM):
    invertDEM = int(inDEM.maximum) - inDEM
    return invertDEM

def getFenceEl(invDEM):
    ras = invDEM + int((invDEM.maximum - invDEM.minimum)/10)
    return ras

def setupInversionZones(regions2Fix):
##                                ## Regions2Fix are areas to clean channel, code to 1, else 0
    ndRcls = Con(IsNull(regions2Fix), 0, 1)
    ## Expand the region to fix by 2 to create an area to burn in
    expFix = Expand(ndRcls, 2, 1)

    ## Create a 'fence' around the area we want to fix,
    ## should be larger than elevation difference in area
    ndPlus = expFix + ndRcls

    return ndRcls, expFix, ndPlus

def fixByInversion(ndPlus, fenceEl, invertTargetDEM, spot4Hole, ndRcls, bestestDEM):

    ## To correct an inverted DEM, remember water must flow out the upstream ends (before inversion), not downstream!
    ## ndPlus is 2 at regions2Fix, else 1 in area to process
    fencedRegion2Fix = Pick(ndPlus, [fenceEl, invertTargetDEM])#originalDEM])

    holeAtMin = Con(IsNull(spot4Hole) == 1, fencedRegion2Fix, '')
    
    fillNdBarrier = Fill(holeAtMin)
##                                    faFence = FlowAccumulation(FlowDirection(fillNdBarrier))

    ## Calculate the difference between the original inverted and filled inverted DEMs
    ## This should be the change to apply to the initial DEM make things flow
    fillNdDif = fillNdBarrier - fencedRegion2Fix

    correctionToApply = Con(ndRcls, bestestDEM - fillNdDif)#originalDEM - fillNdDif)
    correctedDEM = Int(Con(IsNull(correctionToApply), bestestDEM, correctionToApply))#originalDEM, correctionToApply))

    return correctedDEM


def fixByInversionByPath(ndPlus, fenceEl, invertTargetDEM, spot4Hole, ndRcls, bestestDEM, unique_path):

    ## To correct an inverted DEM, remember water must flow out the upstream ends (before inversion), not downstream!
    ## ndPlus is 2 at regions2Fix, else 1 in area to process
    fencedRegion2Fix = Pick(ndPlus, [fenceEl, invertTargetDEM])#originalDEM])

    holeAtMin = Con(IsNull(spot4Hole) == 1, fencedRegion2Fix, '')
    
    fillNdBarrier = Fill(holeAtMin)
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
    pathHoleElevation = ZonalStatistics(unique_path, 'VALUE', spot4HoleElevation, 'MINIMUM')
    correctionToApply = Con(correctionToApply_unfiltered > pathHoleElevation, pathHoleElevation, correctionToApply_unfiltered)
    correctedDEM = Int(Con(IsNull(correctionToApply), bestestDEM, correctionToApply))
    
    return correctedDEM


def createInvertAndExpand(expFix, originalDEM, ndThRivers, nd2WorkTF, bsNbr):
    ## Smooth and expand the original DEM in areas of no data
    targetDEM = Con(expFix, originalDEM)

    ## Minimum thickness for an area to possible be a 'pond' or river
    ndThickCrit = str(int(ndThRivers.maximum))

## Create smoothed DEM (used in wide rivers) and invert DEM of area to fix
    fsNdWidthMin = FocalStatistics(targetDEM, "RECTANGLE " + ndThickCrit + " "  + ndThickCrit + " CELL", "MINIMUM")
    ndSmoothDEM = Con(nd2WorkTF, fsNdWidthMin, targetDEM)

    invertTargetDEM = createInvertDEM(ndSmoothDEM)

## Expand DEM to create valid values beyond edges
    filter5x5InvertDEM = FocalStatistics(invertTargetDEM, bsNbr, 'MEAN')
    expInvertDEM = Con(IsNull(invertTargetDEM), filter5x5InvertDEM, invertTargetDEM)

    return invertTargetDEM, expInvertDEM


def calcDsFd(cmbTable, wsLvl, goodFd, fdFld, inm, gdb, sfx, version, huc12, ):
    try:
        wsBearingTbl = arcpy.CreateTable_management(inm, 'ws_bear_tbl' + sfx)
        arcpy.AddField_management(wsBearingTbl, 'WS', 'LONG')
        arcpy.AddField_management(wsBearingTbl, 'DS_BEARING', 'DOUBLE')
        icur = arcpy.da.InsertCursor(wsBearingTbl, ['WS', 'DS_BEARING'])
        if version.find('10.5') > -1 or version.find('10.6') > -1:
            goodFdName = getfields(cmbTable)[-1]
        else:
            goodFdName = str(goodFd).split('\\')[-1]
            
        with arcpy.da.SearchCursor(cmbTable, ['VALUE', wsLvl.name, goodFdName], sql_clause = (None, 'ORDER BY ' + wsLvl.name + ', ' + goodFdName)) as scur:
            for i, srow in enumerate(scur):
                if i == 0:
                    prevWs = srow[1]
                    bearingList = []
                curWs = srow[1]
                if curWs != prevWs:
                    ## calculate bearing for ws
                    numberFd = len(bearingList)
                    if numberFd > 4:
##                        log.warn('WARNING: number of bearing values is ' + str(numberFd) + ' for ' + huc12 + ' and ws ' + str(curWs))
                        print('WARNING: number of bearing values is ' + str(numberFd) + ' for ' + huc12 + ' and ws ' + str(curWs))
                    dx = 0
                    dy = 0
                    for j in range(0, numberFd):
                        dx += bearingList[j][0]
                        dy += bearingList[j][1]
                    bearingAngle = round(math.atan2(dx, dy)*(360.0/(2*math.pi)), 2)
                    ## store previous result
                    icur.insertRow([prevWs, bearingAngle])
                    ## reset list for more bearing calculations
                    prevWs = srow[1]
                    bearingList = []

                fd = srow[2]
                if fd == 1:
                    easting = 100.0
                    northing = 0.0
                elif fd == 2:
                    easting = 70.71
                    northing = -70.71
                elif fd == 4:
                    easting = 0.0
                    northing = -100.0                
                elif fd == 8:
                    easting = -70.71
                    northing = -70.71
                elif fd == 16:
                    easting = -100.0
                    northing = 0.0
                elif fd == 32:
                    easting = -70.71
                    northing = 70.71
                elif fd == 64:
                    easting = 0.0
                    northing = 100.0                
                elif fd == 128:
                    easting = 70.71
                    northing = 70.71
                bearingList.append([easting, northing])

        del icur
    ##    copytbl(verbose, wsBearingTbl, gdb)

        return wsBearingTbl

    except Exception as e:
##        log.debug(e.message)
        arcpy.AddError(e.message)

        # Get the traceback object
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]

        # Concatenate information together concerning the error into a message string
        pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
##        log.warn(pymsg)


def fullZoneByZs(partialZone, fullZone):
    """Retrieve or recover the full extent of zone(s) from a partial extraction of those zone(s)
    """
    zonal = ZonalStatistics(fullZone, 'VALUE', partialZone, 'MAXIMUM')
    resurrectedFull = Con(zonal, fullZone, where_clause = 'VALUE > 0')
##    resurrectedFull = Con(ZonalStatistics(fullZone, 'VALUE', partialZone, 'MAXIMUM') > 0, fullZone)
    return resurrectedFull

def rescaleExtent(layer_extent, fraction = 0.10):
    width = layer_extent.width
    layer_extent.XMin = layer_extent.XMin - fraction/2.0*width
    layer_extent.XMax = layer_extent.XMax + fraction/2.0*width
    height = layer_extent.height
    layer_extent.YMin = layer_extent.YMin - fraction/2.0*height
    layer_extent.YMax = layer_extent.YMax + fraction/2.0*height

    return layer_extent

def wait_for_timeout(proc, seconds):
    """Wait for a process to finish or raise exception after timeout
    possible alternative method - https://codereview.stackexchange.com/questions/142828/python-executer-that-kills-processes-after-a-timeout
    """
    start = time.time()
    end = start + seconds
    interval = min(seconds / 1000.0, 0.5)

    while True:
        result = proc.poll()
        if result is not None:
            return result
        if time.time() >= end:
            raise RuntimeError("Process timed out")
        time.sleep(interval)

def printParameters(parameters):
    for o,i in enumerate(parameters):
        if o == 0:
            print('\t["' + i.replace('\\', '/') + '",')
        elif o == len(parameters) - 1:
            print('\t"' + i.replace('\\', '/') + '"]')
        else:
            print('\t"' + i.replace('\\', '/') + '",')

def runAndPrint(runList, timeout, env):
    subp = subprocess.Popen(runList, env=env, stdout = subprocess.PIPE, stderr=subprocess.PIPE)#, shell = True)
    print('\tRunning process pid ' + str(subp.pid) + ' at ' + str(time.asctime()))
    try:
        wait_for_timeout(subp, timeout)
    except RuntimeError:
        print('\tRuntime Error - killing process pid ' + str(subp.pid) + ' after ' + str(timeout/60.0) + ' minutes')
        #will take a bit for child's child processes to die (las2las, arcpy functions)
        subp.terminate()

    except:
        print('\tsubprocess exception returned ' + str(subp.returncode))
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

    if subp.returncode != 0:
        print('\tsubprocess failure returned exit code: ' + str(subp.returncode))
        comms = subp.communicate()
        for comm in comms:#stdout, stderr
            lines = comm.splitlines()
            if comm == comms[0]:
                print('\t' + 'STD OUT')
            elif comm == comms[1]:
                print('\n\t' + 'STD ERR')
            for line in lines:
                print('\t' + line)
        
    return subp

def runAndPrintAlways(runList, timeout, env):
    subp = subprocess.Popen(runList, env=env, stdout = subprocess.PIPE, stderr=subprocess.PIPE)#, shell = True)
    print('\tRunning process pid ' + str(subp.pid) + ' at ' + str(time.asctime()))
    try:
        wait_for_timeout(subp, timeout)
    except RuntimeError:
        print('\tRuntime Error - killing process pid ' + str(subp.pid) + ' after ' + str(timeout/60.0) + ' minutes')
        #will take a bit for child's child processes to die (las2las, arcpy functions)
        subp.terminate()

    except:
        print('\tsubprocess exception returned ' + str(subp.returncode))
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

##    if subp.returncode != 0:
    print('\tsubprocess failure returned exit code: ' + str(subp.returncode))
    comms = subp.communicate()
    for comm in comms:#stdout, stderr
        lines = comm.splitlines()
        if comm == comms[0]:
            print('\t' + 'STD OUT')
        elif comm == comms[1]:
            print('\n\t' + 'STD ERR')
        for line in lines:
            print('\t' + line)
    
    return subp


def joinDict(in_data, in_field, join_data, join_field, fields_to_join, in_fields_to_add = []):
    '''A function to use Python dictionaries to join data to a table using update cursors.
    Much faster than using arcpy.JoinField_management or setting up a Join.
    in_data - the table to which data will be added
    in_field - the field on which to join in_data as a string
    join_data - the table from which data will be added
    join_field - the field on which to join the data as a string
    fields_to_join - list of field names to join, all strings
    (optional arg) in_fields_to_add - new field names, in respective order,
    to add to the in_data table, field types come from fields_to_join,
    defualt is to use names from fields_to_join'''

    if in_fields_to_add != []:
        assert len(in_fields_to_add) == len(fields_to_join), "in_fields_to_add is not of same length as fields_to_join"

### Build a dictionary from a da SearchCursor with unique key values storing the attributes to join
    valueDict = {r[0]:(r[1:]) for r in arcpy.da.SearchCursor(join_data, [join_field] + fields_to_join)}

    if in_fields_to_add == []:
        in_fields_to_add = fields_to_join

    for i, field in enumerate(fields_to_join):
        ### Crosswalk the field object types to Add Field parameters
        fieldObj = arcpy.ListFields(join_data, field)[0]
        if fieldObj.type == 'SmallInteger':
            newType = 'SHORT'
        elif fieldObj.type == 'Integer':
            newType = 'LONG'
        elif fieldObj.type == 'Single':
            newType = 'FLOAT'
        elif fieldObj.type == 'Double':
            newType = 'DOUBLE'
        elif fieldObj.type == 'String':
            newType = 'TEXT'
        else:
            newType = fieldObj.type
        ##  Add the field with appropriate type
        if newType == 'TEXT':
            arcpy.AddField_management(in_data, in_fields_to_add[i], newType, field_length = fieldObj.length)
        else:
            arcpy.AddField_management(in_data, in_fields_to_add[i], newType)
      
    with arcpy.da.UpdateCursor(in_data, [in_field] +  in_fields_to_add) as ucur:  
        for urow in ucur:  
            # store the Join value of the row being updated in a keyValue variable  
            keyValue = urow[0]  
            # verify that the keyValue is in the Dictionary  
            if keyValue in valueDict:
                for i, field in enumerate(fields_to_join):
                    urow[i+1] = valueDict[keyValue][i]
                ucur.updateRow(urow)
      
    del valueDict


def visualizeExtent(exob, outfc, frameExtent):
    '''A function to take an extent object and convert it into a polygon for
    visualization during troubleshooting.
    exob - an extent object
    outfc - a path and feature class name'''
    XMAX = frameExtent.XMax  
    XMIN = frameExtent.XMin  
    YMAX = frameExtent.YMax  
    YMIN = frameExtent.YMin  
    pnt1 = arcpy.Point(XMIN, YMIN)  
    pnt2 = arcpy.Point(XMIN, YMAX)  
    pnt3 = arcpy.Point(XMAX, YMAX)  
    pnt4 = arcpy.Point(XMAX, YMIN)  
    array = arcpy.Array()  
    array.add(pnt1)  
    array.add(pnt2)  
    array.add(pnt3)  
    array.add(pnt4)  
    array.add(pnt1)  
    polygon = arcpy.Polygon(array)  
    feature = arcpy.CopyFeatures_management(polygon, outfc)
    return feature

def cleanupOther(procDir, log = None, sgdb = None, inm = None):
    '''use arcpy.da.walk to traverse a directory and delete anything found during the traversal.
    Then proceed to delete anything found in a RAM directory
    omit or pass None for log and logging will be avoided
    omit or pass None for sgdb and/or inm to forego deleting these locations'''
    for root, listDirs, listFiles in arcpy.da.Walk(procDir):
        if log is not None:
            log.info('deleting files in ' + root)
##        else:
##            print('deleting files in ' + root)
        for f in listFiles:
            if arcpy.env.workspace != root:
                arcpy.env.workspace = root
            try:
                arcpy.Delete_management(f)
            except:
                if log is not None:
                    log.warning("couldn't delete " + f)
                else:
                    print("couldn't delete " + f)

    if sgdb is not None:
        try:
            arcpy.Delete_management(sgdb)
        except:
            if log is not None:
                log.warning("couldn't delete " + sgdb)

    if inm is not None:
        arcpy.env.workspace = inm
        fcs = arcpy.ListFeatureClasses()
        for fc in fcs:
            arcpy.Delete_management(fc)
        del fcs
        tbls = arcpy.ListTables()
        for tbl in tbls:
            arcpy.Delete_management(tbl)
        del tbls

    del root, listDirs, listFiles
