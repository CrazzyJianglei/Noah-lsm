import f90nml
import netCDF4
import numpy
import sys
import math
import random
import time
import os
from netCDF4 import Dataset
import warnings

warnings.filterwarnings("ignore")

# lis_module
from numpy import tile


class lisdomain:
    def __init__(self):
        self.nch = 0
        self.glbnch = 0
        self.lsm = 0
        self.soil = 0
        self.elev = 0
        self.ngrid = 0
        self.glbngrid = 0
        self.domain = 0
        self.landcover = 0
        self.gnc = 0
        self.gnr = 0
        self.lnc = int(0)
        self.lnr = int(0)
        self.pnc = 0
        self.pnr = 0
        self.ic = 0
        self.ir = 0
        self.maxt = 0
        self.npesx = 0
        self.npesy = 0

        self.mina = 0
        self.udef = 0
        self.gridDesc = [0 for x in range(0, 50)]
        self.soil_gridDesc = [0 for x in range(0, 6)]
        self.elev_gridDesc = [0 for x in range(0, 6)]
        self.lc_gridDesc = [0 for x in range(0, 6)]


d = lisdomain()


class lisforcing:
    def __init__(self):
        self.force = 0
        self.ecor = 0
        self.nforce = 0
        self.nf = 0
        self.nmif = 0
        self.rstflag = 0
        self.gridchange = 0
        self.interp = 0
        self.latmax = 0
        self.shortflag = 0
        self.longflag = 0
        self.findtime1 = 0
        self.findtime2 = 0
        self.findagrtime1 = 0
        self.findagrtime2 = 0
        self.f00_flag = 0
        self.f06_flag = 0
        self.obsforc = 0


f = lisforcing()


class lisparameters:
    mfile: str

    def __init__(self):
        self.lai = 0
        self.nt = 0
        self.vclass = 0
        self.laiflag = 0
        self.saiflag = 0
        # 对应fortran integer型，将初始值赋值为0
        self.mfile = ''
        self.vfile = ''
        self.safile = ''
        self.clfile = ''
        self.po1file = ''
        self.po2file = ''
        self.po3file = ''
        self.slfile = ''
        self.aspfile = ''
        self.curvfile = ''
        self.sifile = ''
        self.avhrrdir = ''
        self.modisdir = ''
        self.iscfile = ''
        self.elevfile = ''
        # 对应fortran character*70型，将初始值赋值为''空的字符串
        self.lai_gridDesc = [0 for x in range(0, 6)]
        self.laitime = 0.0
        self.saitime = 0.0


p = lisparameters()


class lisoutput:
    def __init__(self):
        self.wfor = 0
        self.wtil = 0
        self.wout = 0
        self.wsingle = 0
        self.wparam = 0
        self.expcode = 0
        self.startcode = 0
        self.foropen = 0
        self.numoutf = 0
        self.fidgm = 0
        self.rcgm = 0
        self.fidtm = 0
        self.rctm = 0
        self.plevel = 0
        self.odir = ''
        self.dfile = ''


o = lisoutput()


class listime:
    def __init__(self):
        self.sss = 0
        self.sdoy = 0
        self.smn = 0
        self.shr = 0
        self.sda = 0
        self.smo = 0
        self.syr = 299
        self.endcode = 0
        self.ess = 0
        self.emn = 0
        self.edoy = 0
        self.ehr = 0
        self.eda = 0
        self.emo = 0
        self.eyr = 0
        self.tstype = 0
        self.ts = 0
        self.tscount = 0
        self.yyyymmdd = 0
        self.hhmmss = 0
        self.doy = 0
        self.yr = 0
        self.mo = 0
        self.da = 0
        self.hr = 0
        self.mn = 0
        self.ss = 0
        self.endtime = 0
        self.pda = 0

        self.time = 0.0
        self.etime = 0.0
        # real*8
        self.gmt = 0.0
        self.egmt = 0.0
        self.sgmt = 0.0
        # real


t = listime()
tt = listime()


class lisassimil:
    def __init__(self):
        self.daalg = 0
        self.daobs = 0
        self.davar = 0
        self.nstv = 0
        self.noutv = 0
        self.ndav = 0
        self.nens = 0
        self.rensem = 0
        self.ndva = 0

        self.obslo = 0.0
        self.obsla = 0.0
        self.obelo = 0.0
        self.obela = 0.0
        self.oblor = 0.0
        self.oblar = 0.0

        self.fvlfn = ''
        self.svlfn = ''

        self.pSLGA = True
        self.pFORCE = True
        # fortran logical 型
        self.win = 0


a = lisassimil()


class lisdec:
    def __init__(self):
        self.d = d
        self.f = f
        self.p = p
        self.t = t
        self.o = o
        self.a = a


lis = lisdec()


# grid_module
class griddec:
    def __init__(self):
        self.lat = 0.0
        self.lon = 0.0
        self.forcing = [0 for x in range(0, 10)]
        self.fgrd = [0 for x in range(0, 13)]
        self.elev = 0.0
        self.SLOPE = 0.0
        self.aspect = 0.0
        self.curv = 0.0


griddec = griddec()


# tile_module
class tiledec:
    def __init__(self):
        self.col = 0
        self.row = 0
        self.index = 0
        self.vegt = 0
        self.ensem = 0
        self.fgrd = 0.0


tiledec = tiledec()


# noah_module
class noahdec:
    def __init__(self):
        self.ts = 0  # Timestep (seconds)   时间步长（秒）
        self.maxt = 0  # Maximum tiles per grid  每个网格最大的tiles数
        self.SIBVEG = 0  # UMD to SIB Vegetation Class Index value UMD至SIB植被等级指数值
        self.NSLAY = 0  # Number of NOAH soil layers (4) NOAH的土地层数
        self.COUNT = 0
        self.ZOBSOIL = [0 for x in range(0, 1)]  # Zobler Soil Classes (LIS%NCH) Zobler土壤等级（LIS%NCH）
        self.VEGP = [0.0 for x in range(0, 7)]  # Static vegetation parameter values, dim(NOAH_NVEGP) 静态植被参数
        self.VEGIP = 0.0  # Interpolated Green Fraction from monthly parameters 基于月参数的插值绿度分数
        self.VEGMP1 = 0.0  # Month 1 Greenness Fraction Value 第一个月绿色分数值
        self.VEGMP2 = 0.0  # Month 2 Greenness Fraction Value 第二个月绿色分数值
        self.ALBSF1 = 0.0  # Date 1 Snow-Free Albedo Value 日期1无雪反照率值
        self.ALBSF2 = 0.0  # Date 2 Snow-Free Albedo Value 日期2无雪反照率值
        self.SOILP = [0.0 for x in range(0, 10)]  # Static soil parameter values, dim(NOAH_NSOILP) 静态土壤参数值
        self.ALBSF = 0.0  # Quarterly Snow-Free Albedo dataset 季度无雪反照率数据集
        self.MXSNALB = 0.0  # Maximum snow ALBEDO dataset 最大积雪反照率数据集
        self.TEMPBOT = 0.0  # Bottom boundary temperature 底部边界温度
        # -------------------------------------------------------------------------
        # NOAH-State Variables NAOH状态变量
        # -------------------------------------------------------------------------
        self.T1 = 0.0  # NOAH Skin Temperature (K)   NOAH 表层温度
        self.CMC = 0.0  # NOAH Canopy Water Content NOAH 冠层含水量
        self.SNOWH = 0.0  # NOAH Actual Snow depth (m) NOAH 实际积雪深度
        self.SNEQV = 0.0  # NOAH Water Equivalent Snow Depth (m) NOAH 水当量雪深
        self.STC = [0.0 for x in range(0, 4)]  # NOAH Soil Temperaure (4 layers) NOAH 土壤温度
        self.SMC = [0.0 for x in range(0, 4)]  # NOAH Soil Moisture (4 layers) NOAH 土壤水分
        self.SH2O = [0.0 for x in range(0, 4)]  # NOAH Liquid-only soil moisture (4 layers) NOAH 纯液体土壤水分
        self.CH = 0.0  # NOAH Heat/moisture exchange coef. NOAH 热湿交换系数
        self.CM = 0.0  # NOAH Momentum exchange coef. NOAH 动量交换系数
        self.forcing = [0.0 for x in range(0, 10)]  # TILE forcing. tile 强迫
        self.vegt = 0.0  # vegetation type of tile tile 的植被种类
        # -----------------------------------------------------------------------
        #  NOAH-Output variables NOAH 输出变量
        # -----------------------------------------------------------------------
        self.swnet = 0.0
        self.lwnet = 0.0
        self.qle = 0.0
        self.qh = 0.0
        self.qg = 0.0
        self.snowf = 0.0
        self.rainf = 0.0
        self.evap = 0.0
        self.qs = 0.0
        self.qsb = 0.0
        self.qsm = 0.0
        self.avgsurft = 0.0
        self.ALBEDO = 0.0
        self.swe = 0.0
        self.soilmoist1 = 0.0
        self.soilmoist2 = 0.0
        self.soilmoist3 = 0.0
        self.soilmoist4 = 0.0
        self.soilwet = 0.0
        self.ecanop = 0.0
        self.canopint = 0.0
        self.tveg = 0.0
        self.esoil = 0.0
        self.rootmoist = 0.0
        self.soilm_prev = 0.0
        self.swe_prev = 0.0
        self.vegmp1 = 0.0
        self.vegmp2 = 0.0


# noahdec = noahdec()


# noahdrv_module
class noahdrvdec:
    def __init__(self):
        self.noahopen = 0  # Keeps track of opening files 跟踪打开的文件
        self.numout = 0  # Counts number of output times for Noah 统计NOAH的输出次数
        self.noah_nvegp = 0  # Number of static vegetation parameter 静态植被参数数量
        self.noah_nsoilp = 0  # Number of static soil parameters 静态土壤参数数量
        self.noah_zst = 0  # Number of Zobler soil classes Zobler 土壤等级的数量
        self.noah_gflag = 0  # Time flag to update gfrac files 更新 gfrac 文件的时间标志
        self.noah_albtime = 0  # Time flag to update ALBEDO files 更新反照率文件的时间标志
        self.noah_aflag = 0  # Time flag to update ALBEDO files 更新反照率文件的时间标志
        self.noah_albdchk = 0  # Day check to interpolate ALB values 日间检查以插入 ALB 值
        self.noah_gfracdchk = 0  # Day check to interpolate gfrac value 日间检查以插入 gfrac 值
        self.varid = [0 for x in range(0, 21)]  # For netcdf output 对于 netcdf 输出
        self.noah_rfile = ''  # NOAH Active Restart File NOAH主动重启文件
        self.noah_vfile = ''  # NOAH Static Vegetation Parameter File NOAH 静态植被参数文件
        self.noah_sfile = ''  # NOAH Soil Parameter File  NOAH土壤参数文件
        self.noah_mgfile = ''  # NOAH Monthly Veg Green Frac NOAH 每月植被绿度
        self.noah_albfile = ''  # NOAH Quart Snow-free ALBEDO NOAH 季度无雪反照率
        self.noah_mxsnal = ''  # NOAH GLDAS max snow ALBEDO NOAH GLDAS 最大雪反射率
        self.noah_tbot = ''  # NOAH GLDAS Bottom Temp NOAH GLDAS 底部温度
        self.noah_gridDesc = [0.0 for x in range(0, 6)]  # grid definition for LAI dataset LAI 数据集的网格定义
        self.noah_gfractime = 0.0  # Time flag to update gfrac files 更新 gfrac 文件的时间标志
        self.noah_ism = 0.0  # NOAH Initial Soil Moisture (m3/m3) NOAH 初始土壤水分 (m3/m3)
        self.noah_it = 0.0  # NOAH Initial Soil Temperature (K) NOAH 初始土壤温度 (K)
        self.writeintn = 0.0  # NOAH Output Interval (hours) NOAH 输出间隔（小时）


noahdrv = noahdrvdec()


# noahderda_module
class noahderda:
    def __init__(self):
        self.delta = 0.0  # the delta value in calculating Ep #  计算Ep的 delta 值
        self.RR = 0.0  # the RR value in calculating Ep 计算Ep时的RR值
        self.Ep = 0.0  # the potential evaporation 潜在蒸发
        self.K = 0.0  # the K Value used to calculate soil heat flux 用于计算土壤热通量的 K 值
        self.soildrymoist = 0.0  # the soil dry moist for the first layer  第一层土壤干湿
        self.CMC = 0.0  # the water quantity captured by the vegetation cnnopy 植被cnnopy捕获的水量
        self.FXEXP = 0.0
        self.n = 0.0  # the index used to calculte the transpiration 用于计算蒸腾的指数
        self.Bc = 0.0  # vegetation resistence coefficient 植被抵抗系数


noahder = noahderda()

# spmdMod
zd_spmd = {
    'masterproc': True,
    'iam': 0,
    'npes': 1,
    'c_str': [0 for x in range(50)],
    'c_end': [0 for x in range(50)],
    'r_str': [0 for x in range(50)],
    'r_end': [0 for x in range(50)],
    'local_nc': [0 for x in range(50)],
    'local_nr': [0 for x in range(50)],
    'deltas': [0 for x in range(50)],
    'offsets': [0 for x in range(50)],
    'tdeltas': [0 for x in range(50)],
    'toffsets': [0 for x in range(50)],
    'gdeltas': [0 for x in range(50)],
    'goffsets': [0 for x in range(50)],
    'nchs': [0 for x in range(50)],
    'ngrids': [0 for x in range(50)],
}


def spmd_init():
    ier = 0  # 并行计算用到的变量，这里暂不用管。
    if zd_spmd['iam'] == 0:
        zd_spmd['masterproc'] = True
    else:
        zd_spmd['masterproc'] = False
    if zd_spmd['masterproc']:
        print('Starting GSFC-LIS')
        print('** Number of Procs **', zd_spmd['npes'])
    zd_spmd['c_str'] = [0 for x in range(0, zd_spmd['npes'])]
    zd_spmd['c_end'] = [0 for x in range(0, zd_spmd['npes'])]
    zd_spmd['r_str'] = [0 for x in range(0, zd_spmd['npes'])]
    zd_spmd['r_end'] = [0 for x in range(0, zd_spmd['npes'])]
    zd_spmd['local_nc'] = [0 for x in range(0, zd_spmd['npes'])]
    zd_spmd['local_nr'] = [0 for x in range(0, zd_spmd['npes'])]
    zd_spmd['deltas'] = [0 for x in range(0, zd_spmd['npes'])]
    zd_spmd['offsets'] = [0 for x in range(0, zd_spmd['npes'])]
    zd_spmd['tdeltas'] = [0 for x in range(0, zd_spmd['npes'])]
    zd_spmd['toffsets'] = [0 for x in range(0, zd_spmd['npes'])]
    zd_spmd['gdeltas'] = [0 for x in range(0, zd_spmd['npes'])]
    zd_spmd['goffsets'] = [0 for x in range(0, zd_spmd['npes'])]
    zd_spmd['nchs'] = [0 for x in range(0, zd_spmd['npes'])]
    zd_spmd['ngrids'] = [0 for x in range(0, zd_spmd['npes'])]


# readdomain_default
def readdomain_default():
    run_dd = [0.0 for i in range(7)]
    param_dd = [0.0 for i in range(6)]
    nml = f90nml.read('D:/NoahDate/NOAH/lis.crd')
    for i in range(0, 7):
        run_dd[i] = nml['run_domain']['run_dd'][i]
    for i in range(0, 6):
        param_dd[i] = nml['param_domain']['param_dd'][i]
    if zd_spmd['masterproc']:
        print('DOMAIN details:')
    if zd_spmd['npes'] == 1:
        lis.d.npesx = 1
        lis.d.npesy = 1
    if lis.d.npesx * lis.d.npesy != zd_spmd['npes']:
        print('428 Layout does not match the number of processors...')
        print('429 npex, npey, ', lis.d.npesx, 'x', lis.d.npesy, '!=', zd_spmd['npes'])
        print('430 Stopping program ..')
        # endrun()
    lat_str = [0.0 for i in range(lis.d.npesy)]
    lon_str = [0.0 for i in range(lis.d.npesx)]
    lat_end = [0.0 for i in range(lis.d.npesy)]
    lon_end = [0.0 for i in range(lis.d.npesx)]
    # if (run_dd[4] - run_dd[2]) % run_dd[6] == 0:
    #     lnc = (int((run_dd[4] - run_dd[2]) / run_dd[6]) + 1) / lis.d.npesx
    #     lnr = (int((run_dd[3] - run_dd[1]) / run_dd[5]) + 1) / lis.d.npesy
    # else:
    #     lnc = (int((run_dd[4] - run_dd[2]) / run_dd[6]) + 1) / lis.d.npesx
    #     lnr = (int((run_dd[3] - run_dd[1]) / run_dd[5]) + 1) / lis.d.npesy

    lnc = (int((run_dd[4] - run_dd[2]) / run_dd[6]) + 1) / lis.d.npesx
    lnr = (int((run_dd[3] - run_dd[1]) / run_dd[5]) + 1) / lis.d.npesy

    lis.d.lnc = int(lnc)
    lis.d.lnr = int(lnr)

    for i in range(1, lis.d.npesy + 1):
        lat_str[i - 1] = run_dd[1] + (i - 1) * lnr * run_dd[6]
    for i in range(1, lis.d.npesx + 1):
        lon_str[i - 1] = run_dd[2] + (i - 1) * lnc * run_dd[5]
    for i in range(1, lis.d.npesy + 1):
        if (i + 1) > lis.d.npesy:
            lat_end[i - 1] = run_dd[3]
        else:
            lat_end[i - 1] = lat_str[i] - run_dd[6]

    for i in range(1, lis.d.npesx + 1):
        if (i + 1) > lis.d.npesx:
            lon_end[i - 1] = run_dd[4]
        else:
            lon_end[i - 1] = lat_str[i] - run_dd[5]
    count = 0
    for i in range(1, lis.d.npesx + 1):
        for j in range(1, lis.d.npesy + 1):
            if zd_spmd['iam'] == count:
                lis.d.gridDesc[3] = lat_str[j - 1]
                lis.d.gridDesc[4] = lon_str[i - 1]
                lis.d.gridDesc[6] = lat_end[j - 1]
                lis.d.gridDesc[7] = lon_end[i - 1]
            count = count + 1
    print('468 local domain ', zd_spmd['iam'], lis.d.gridDesc[3], lis.d.gridDesc[6], lis.d.gridDesc[4],
          lis.d.gridDesc[7])
    zd_spmd['r_str'][zd_spmd['iam']] = round((lis.d.gridDesc[3] - run_dd[1]) / run_dd[5]) + 1
    zd_spmd['r_end'][zd_spmd['iam']] = round((lis.d.gridDesc[6] - run_dd[1]) / run_dd[5]) + 1
    zd_spmd['c_str'][zd_spmd['iam']] = round((lis.d.gridDesc[4] - run_dd[2]) / run_dd[6]) + 1
    zd_spmd['c_end'][zd_spmd['iam']] = round((lis.d.gridDesc[7] - run_dd[2]) / run_dd[6]) + 1
    print('473 global indices ', zd_spmd['iam'], ',', zd_spmd['c_str'][zd_spmd['iam']],
          zd_spmd['r_str'][zd_spmd['iam']],
          zd_spmd['c_end'][zd_spmd['iam']], zd_spmd['r_end'][zd_spmd['iam']])
    lis.d.gnc = round((run_dd[4] - run_dd[2]) / run_dd[6]) + 1
    lis.d.gnr = round((run_dd[3] - run_dd[1]) / run_dd[5]) + 1
    lis.d.gridDesc[0] = run_dd[0]
    lis.d.gridDesc[8] = run_dd[5]
    lis.d.gridDesc[43] = param_dd[0]
    lis.d.gridDesc[44] = param_dd[1]
    lis.d.gridDesc[46] = param_dd[2]
    if lis.d.gridDesc[0] == 0:
        lis.d.gridDesc[9] = run_dd[6]
        lis.d.gridDesc[49] = param_dd[4]
    lis.d.gridDesc[47] = param_dd[3]
    lis.d.gridDesc[48] = param_dd[5]
    if lis.d.gridDesc[0] == 0:
        lis.d.gridDesc[42] = round((lis.d.gridDesc[46] - lis.d.gridDesc[43]) / lis.d.gridDesc[48]) + 1
        lis.d.gridDesc[41] = round((lis.d.gridDesc[47] - lis.d.gridDesc[44]) / lis.d.gridDesc[49]) + 1
    if lis.d.gridDesc[0] == 0:
        lis.d.gridDesc[5] = 128
        lis.d.gridDesc[10] = 64
        lis.d.gridDesc[19] = 255
        if lis.d.gridDesc[6] < lis.d.gridDesc[3]:
            print('495 lat2 must be greater than lat1')
            print('496 Stopping run...')
            # endrun()
        if lis.d.gridDesc[7] < lis.d.gridDesc[4]:
            print('499 lon2 must be greater than lon1')
            print('500 Stopping run...')
            # endrun()
        if (lis.d.gridDesc[3] - lis.d.gridDesc[43]) % lis.d.gridDesc[8] == 0:
            tl = (round((lis.d.gridDesc[3] - lis.d.gridDesc[43]) / lis.d.gridDesc[8]))
            lis.d.gridDesc[3] = lis.d.gridDesc[43] + tl * lis.d.gridDesc[8]
            print('505 modified gridDesc(4) ', lis.d.gridDesc[3])
        if (lis.d.gridDesc[4] - lis.d.gridDesc[44]) % lis.d.gridDesc[9] == 0:
            tl = (round((lis.d.gridDesc[4] - lis.d.gridDesc[44]) / lis.d.gridDesc[9]))
            lis.d.gridDesc[4] = lis.d.gridDesc[44] + tl * lis.d.gridDesc[9]
            print('509 modified gridDesc(5) ', lis.d.gridDesc[4])
        if (lis.d.gridDesc[46] - lis.d.gridDesc[6]) % lis.d.gridDesc[8] == 0:
            tl = (int((lis.d.gridDesc[46] - lis.d.gridDesc[6]) / lis.d.gridDesc[8]))
            lis.d.gridDesc[6] = lis.d.gridDesc[46] - tl * lis.d.gridDesc[8]
            print('513 modified gridDesc(7) ', lis.d.gridDesc[6])
        if (lis.d.gridDesc[47] - lis.d.gridDesc[7]) % lis.d.gridDesc[9] == 0:
            tl = (int((lis.d.gridDesc[47] - lis.d.gridDesc[7]) / lis.d.gridDesc[9]))
            lis.d.gridDesc[7] = lis.d.gridDesc[47] - tl * lis.d.gridDesc[9]
            print('517 modified gridDesc(8) ', lis.d.gridDesc[7])
        lis.d.gridDesc[1] = round((lis.d.gridDesc[7] - lis.d.gridDesc[4]) / lis.d.gridDesc[9]) + 1
        lis.d.gridDesc[2] = round((lis.d.gridDesc[6] - lis.d.gridDesc[3]) / lis.d.gridDesc[8]) + 1
    for k in range(40, 50 + 1):
        print('(', k, ',', lis.d.gridDesc[k - 1], ')')
    lis.d.pnc = lis.d.gridDesc[41]
    lis.d.pnr = lis.d.gridDesc[42]
    print('524 local domain', '(', lis.d.lnc, lis.d.lnr, ')')
    print('525 parameter domain', '(', lis.d.pnc, lis.d.pnr, ')')
    print('526 running domain', '(', lis.d.gnc, lis.d.gnr, ')')
    zd_spmd['local_nc'][zd_spmd['iam']] = lis.d.lnc  # 下标不用-1
    zd_spmd['local_nr'][zd_spmd['iam']] = lis.d.lnr  # 下标不用-1
    zd_spmd['deltas'][zd_spmd['iam']] = lis.d.lnc * lis.d.lnr  # 下标不用-1
    if zd_spmd['masterproc']:
        zd_spmd['offsets'][0] = 0
        for i in range(1, zd_spmd['npes']):
            zd_spmd['offsets'][i] = zd_spmd['offsets'][i - 1] + zd_spmd['deltas'][i - 1]
    lat_str = []
    lon_str = []
    lat_end = []
    lon_end = []


# read_umdavhrr_mask
def read_umdavhrr_mask():  # 参数为二维列表。
    import numpy as np  # 加载数学函数包，同时也是对矩阵数组的操作包
    import struct
    nc, nr = 0, 0
    line1, line2 = 0, 0
    num_lon_pts = 0
    # 创建一个宽度为3，高度为4的数组 ,
    # myList = [([0] * 3) for i in range(4)]
    # mask = [[1 for j in range(0, lis.d.lnc)] for i in range(0, lis.d.lnr)]

    # mask = [([0.0] * int(lis.d.lnc)) for i in range(0, int(lis.d.lnr))]
    # mask = [([0]*lis.d.lnc) for i in range(0, lis.d.lnr)]
    # localmask=np.zeros([lis.d.lnc,lis.d.lnr], dtype='f4')

    if lis.d.gridDesc[8] != 0.01:
        print('556 MSG: maketiles -- Reading ', lis.p.mfile.rstrip(), ' (', zd_spmd['iam'], ')')
        # vf1list = [nc, nr, mask]
        localmask = vegmask_from_1km()
        print('559 MSG: maketiles -- Done Reading ', lis.p.mfile.rstrip(), ' (', zd_spmd['iam'], ')')
    else:
        print('561 MSG: maketiles -- Reading ', lis.p.mfile.rstrip(), ' (', zd_spmd['iam'], ')')
        binFile = open(lis.p.mfile, 'rb')
        # content = f.readlines()
        localmask = np.zeros([int(lis.d.lnc), int(lis.d.lnr)], dtype='f4')
        nc = int((lis.d.lc_gridDesc[3] - lis.d.lc_gridDesc[1]) / lis.d.lc_gridDesc[5]) + 1
        print('566 执行到了else语句')
        line1 = int((lis.d.gridDesc[3] - lis.d.lc_gridDesc[0]) / lis.d.lc_gridDesc[4]) + 1
        line2 = int((lis.d.gridDesc[4] - lis.d.lc_gridDesc[1]) / lis.d.lc_gridDesc[5]) + 1
        for c in range(0, int(lis.d.lnc)):
            for r in range(0, int(lis.d.lnr)):
                glnc = line2 + c
                glnr = line1 + r
                line = glnr * nc + glnc
                binFile.seek(line * 4)
                context = binFile.read(4)
                tt, = struct.unpack('f', context)
                # print('579', tt)
                localmask[c, r], = struct.unpack('f', context)
                # print(localmask[c, r])
        binFile.close()
        print('583 MSG: maketiles -- Done ', lis.p.mfile.rstrip())
    return localmask


# read_umdavhrr_lc
def read_umdavhrr_lc():  # 参数为三维列表。
    import numpy as np
    # fgrd = [[[0.0] * lis.d.lnc for i in range(int(lis.d.lnr))] for i in range(int(lis.p.nt))]
    isum = 0.0
    nc, nr, veg = int(lis.d.lnc), int(lis.d.lnr), int(lis.p.nt)
    cf1list = [nc, nr, veg]
    if lis.d.gridDesc[8] == 0.01:
        fgrd = np.zeros([int(lis.d.lnc), int(lis.d.lnr), int(lis.p.nt)], dtype='f4')
        # cf1list = [nc, nr, veg]
        vegcover_from_1km(cf1list)
        nc, nr, veg = cf1list[0], cf1list[1], cf1list[2]
        print('597 MSG: maketiles -- Reading ',
              '598 D:/NoahDate/NOAH/CN_LUCC_2005_Geographic_UMDcode_output.dat'.rstrip(), ' (', zd_spmd['iam'], ')')
        line1 = int((lis.d.gridDesc[3] - lis.d.lc_gridDesc[0]) / lis.d.lc_gridDesc[4]) + 1
        line2 = int((lis.d.gridDesc[4] - lis.d.lc_gridDesc[1]) / lis.d.lc_gridDesc[5]) + 1
        num_lon_pts = nc
        for i in range(0, int(lis.d.lnc)):
            for j in range(0, int(lis.d.lnr)):
                for k in range(0, lis.p.nt):
                    fgrd[k][i][j] = veg[k][i + line2 - 1][j + line1 - 1]
    else:
        fgrd = vegcover_from_1km(cf1list)
    return fgrd


# read_elev_gtopo30
def read_elev_gtopo30(reg):  # 对列表包装reg=[elev],目的给实参一个指定大小的列表
    line1, line2, nc_dom = 0, 0, 0
    if lis.f.ecor >= 1:
        reg[0] = [[0.0] * lis.d.lnc for i in range(lis.d.lnr)]
        line1 = round((lis.d.gridDesc[3] - lis.d.elev_gridDesc[0]) / lis.d.gridDesc[8]) + 1
        line2 = round((lis.d.gridDesc[4] - lis.d.elev_gridDesc[1]) / lis.d.gridDesc[9]) + 1
        nc_dom = round((lis.d.elev_gridDesc[3] - lis.d.elev_gridDesc[1]) / lis.d.elev_gridDesc[4]) + 1
        print('619 Reading elevation file ...', lis.p.elevfile)
        lis_open_file(21, file=lis.p.elevfile, form='unformatted', access='direct', recl=4, status='old',
                      script='getelev.pl')
        f = open(lis.p.elevfile, 'r')
        content = f.readlines()
        for r in range(1, lis.d.lnr + 1):
            for c in range(1, lis.d.lnc + 1):
                glnc = line2 + c - 1
                glnr = line1 + r - 1
                line = (glnr - 1) * nc_dom + glnc
                reg[0][c - 1][r - 1] = content[line - 1]
        f.close()
        lis_log_msg('631 MSG: read_elev_gtopo30 -- Done reading elevation difference file')


# calculate_domveg
def calculate_domveg(fgrd):  # 两个列表参数。
    import numpy as np
    lnc = int(lis.d.lnc)
    lnr = int(lis.d.lnr)
    npt = int(lis.p.nt)
    tsum = np.zeros([lnc, lnr], dtype='f4')
    fvt = np.zeros([npt], dtype='f4')

    for r in range(1, lnr + 1):
        for c in range(1, lnc + 1):  # c的取值范围为1~lnc
            rsum = 0.0
            for t in range(1, npt + 1):
                # print('577',c,r,lnc)
                if fgrd[c - 1, r - 1, t - 1] < lis.d.mina:
                    fgrd[c - 1, r - 1, t - 1] = 0.0
                rsum = rsum + fgrd[c - 1, r - 1, t - 1]
            if rsum > 0.0:
                for t in range(1, npt + 1):
                    fgrd[c - 1, r - 1, t - 1] = fgrd[c - 1, r - 1, t - 1] / rsum
                ksum = 0.0
                for t in range(1, npt + 1):
                    ksum = ksum + fgrd[c - 1, r - 1, t - 1]
                if ksum < 0.9999 or ksum > 1.0001:
                    print('666 Error1 in vegetation tiles', ksum, c, r)

    pveg = np.zeros([lnc, lnr, npt], dtype='f4')
    for r in range(1, lnr + 1):
        for c in range(1, lnc + 1):
            for t in range(1, npt + 1):
                fvt[t - 1] = fgrd[c - 1, r - 1, t - 1]
            for i in range(1, npt + 1):
                max = 0.0
                t = 0
                for j in range(1, npt + 1):
                    if fvt[j - 1] > max:
                        if fgrd[c - 1, r - 1, t - 1] > 0:
                            max = fvt[j - 1]
                            t = j
                if t > 0:
                    pveg[c - 1, r - 1, t - 1] = i
                    fvt[t - 1] = -999.0

    for r in range(1, lnr + 1):
        for c in range(1, lnc + 1):
            rsum = 0.0
            for t in range(1, npt + 1):
                if pveg[c - 1, r - 1, t - 1] < 1:
                    fgrd[c - 1, r - 1, t - 1] = 0.0
                    pveg[c - 1, r - 1, t - 1] = 0
                if pveg[c - 1, r - 1, t - 1] > lis.d.maxt:
                    fgrd[c - 1, r - 1, t - 1] = 0.0
                    fgrd[c - 1, r - 1, t - 1] = 0
                rsum = rsum + fgrd[c - 1, r - 1, t - 1]
            if rsum > 0.0:
                for t in range(1, npt + 1):
                    fgrd[c - 1, r - 1, t - 1] = fgrd[c - 1, r - 1, t - 1] / rsum
                ksum = 0.0
                for t in range(1, npt + 1):
                    ksum = ksum + fgrd[c - 1, r - 1, t - 1]
                tsum[c - 1, r - 1] = ksum
                if ksum < 0.9999 or ksum > 1.0001:
                    print('703 Error2 in vegetation tiles', ksum, c, r)
    return tsum


# create_vegtilespace
def create_vegtilespace(fgrd, tsum, localmask, elev):  # 参数为4个列表。
    import numpy as np

    n_ensem = 1
    lnc = int(lis.d.lnc)
    lnr = int(lis.d.lnr)
    npt = int(lis.p.nt)
    # tsum = np.zeros([lnc, lnr], dtype='f4')
    locallat, locallon = 0.0, 0.0
    domveg = [[0] * 50 for i in range(50)]
    landnveg = 5
    lis.d.nch = 0

    for t in range(1, npt + 1):
        for r in range(1, lnr + 1):
            for c in range(1, lnc + 1):
                #                print('649',localmask.shape)
                if 0.99 < localmask[c - 1][r - 1] < 3.01:  # localmask是列表，列表的下标是[][]
                    if fgrd[c - 1, r - 1, t - 1] > 0.0:
                        lis.d.nch = lis.d.nch + 1
                    if tsum[c - 1, r - 1] == 0.0 and t == landnveg:
                        lis.d.nch = lis.d.nch + 1
    print('725 MSG: maketiles -- nch', lis.d.nch, ' (', zd_spmd['iam'], ')')

    zd_lisdrv['tile'] = [tiledec for i in range(int(lis.d.nch))]
    lis.d.ngrid = 0
    for r in range(1, lnr + 1):
        for c in range(1, lnc + 1):
            if 0.99 < localmask[c - 1][r - 1] < 3.01:
                lis.d.ngrid = lis.d.ngrid + 1

    count = 0
    zd_lisdrv['grid'] = [griddec for i in range(int(lis.d.ngrid))]
    zd_lisdrv['gindex'] = [[0.0] * lnr for i in range(lnc)]

    for r in range(1, lnr + 1):
        for c in range(1, lnc + 1):
            zd_lisdrv['gindex'][c - 1][r - 1] = -1
            if 0.99 < localmask[c - 1][r - 1] < 3.01:
                locallat = lis.d.gridDesc[3] + (r - 1) * lis.d.gridDesc[8]
                locallon = lis.d.gridDesc[4] + (c - 1) * lis.d.gridDesc[9]
                zd_lisdrv['grid'][count].lat = locallat
                zd_lisdrv['grid'][count].lon = locallon
                zd_lisdrv['grid'][count].fgrd = fgrd[c - 1, r - 1, :]
                zd_lisdrv['gindex'][c - 1][r - 1] = count
                count = count + 1
    print('747 MSG: maketiles -- ngrid', lis.d.ngrid, ' (', zd_spmd['iam'], ')')
    if lis.o.wparam == 1:
        domveg = np.zeros([lnc, lnr], dtype='f4')
    count = 0
    for r in range(1, lnr + 1):
        for c in range(1, lnc + 1):
            for t in range(1, npt + 1):
                for m in range(1, n_ensem):
                    if localmask[c - 1][r - 1] > 0.99:
                        if fgrd[c - 1, r - 1, t - 1] > 0.0:
                            count = count + 1
                            zd_lisdrv['tile'][count - 1].ensem = m
                            zd_lisdrv['tile'][count - 1].row = r
                            zd_lisdrv['tile'][count - 1].col = c
                            zd_lisdrv['tile'][count - 1].index = zd_lisdrv['gindex'][c - 1, r - 1]
                            zd_lisdrv['tile'][count - 1].vegt = t
                            if lis.o.wparam == 1:
                                domveg[c - 1, r - 1] = t * 1.0
                            zd_lisdrv['tile'][count - 1].fgrd = fgrd[c - 1, r - 1, t - 1]
                            if lis.f.ecor >= 1:
                                if elev[c - 1][r - 1] == -9999:  # elev为列表[c][r]
                                    elev[c - 1][r - 1] = 0.0
                                zd_lisdrv['grid'][zd_lisdrv['tile'][count - 1].index - 1].elev = elev[c - 1][r - 1]
                        if tsum[r - 1][c - 1] == 0.0 and t == landnveg:
                            zd_lisdrv['tile'][count - 1].ensem = m
                            zd_lisdrv['tile'][count - 1].row = r
                            zd_lisdrv['tile'][count - 1].col = c
                            zd_lisdrv['tile'][count - 1].index = zd_lisdrv['gindex'][c - 1][r - 1]
                            zd_lisdrv['tile'][count - 1].vegt = t
                            if lis.o.wparam == 1:
                                domveg[c - 1, r - 1] = t * 1.0
                            zd_lisdrv['tile'][count - 1].fgrd = 1.0
                            if lis.f.ecor >= 1:
                                if elev[c - 1, r - 1] == -9999:
                                    elev[c - 1, r - 1] = 0.0
                                zd_lisdrv['grid'][zd_lisdrv['tile'][count - 1].index - 1].elev = elev[c - 1][r - 1]
    zd_lisdrv['ntiles_pergrid'] = np.zeros([int(lis.d.gnc * lis.d.gnr)], dtype='f4')
    zd_lisdrv['s_tid'] = np.zeros([int(lis.d.gnc * lis.d.gnr)], dtype='f4')
    ntpg = np.zeros([int(lis.d.gnc * lis.d.gnr)], dtype='f4')
    for t in range(1, int(lis.d.lnc * lis.d.lnr) + 1):
        zd_lisdrv['ntiles_pergrid'][t - 1] = 0
        ntpg[t - 1] = 0
        zd_lisdrv['s_tid'][t - 1] = 0
    for t in range(1, int(lis.d.nch) + 1):
        index = zd_lisdrv['gindex'][zd_lisdrv['tile'][t - 1].row - 1][zd_lisdrv['tile'][t - 1].col - 1]
        if index != -1:
            gid = int(zd_lisdrv['tile'][t - 1].col + zd_lisdrv['tile'][t - 1].row * lis.d.lnc)
            ntpg[gid] = ntpg[gid] + 1
    gtmp = np.zeros([int(lis.d.gnc), int(lis.d.gnr)], dtype='f4')
    gtmp1 = np.zeros([int(lis.d.gnc), int(lis.d.gnr)], dtype='f4')
    gtmp1 = ntpg
    count = 1
    for n in range(1, zd_spmd['npes'] + 1):
        for r in range(zd_spmd['r_str'][n - 1], zd_spmd['r_end'][n - 1]):
            for c in range(zd_spmd['c_str'][n - 1], zd_spmd['c_end'][n - 1]):
                gtmp[c - 1, r - 1] = gtmp1[count - 1]
                count = count + 1
    count = 1
    for r in range(1, int(lis.d.gnr) + 1):
        for c in range(1, int(lis.d.gnc) + 1):
            zd_lisdrv['ntiles_pergrid'][count - 1] = gtmp[c - 1, r - 1]
            count = count + 1
    for t in range(1, int(lis.d.gnc * lis.d.gnr) + 1):
        if zd_lisdrv['ntiles_pergrid'][t - 1] > 0:
            zd_lisdrv['s_tid'][t - 1] = 1
            c = t + 1
            break
        else:
            zd_lisdrv['s_tid'][t - 1] = 0
    for t in range(c, int(lis.d.gnc * lis.d.gnr) + 1):
        zd_lisdrv['s_tid'][t - 1] = zd_lisdrv['s_tid'][t - 2] + zd_lisdrv['ntiles_pergrid'][t - 2]

    del gtmp
    del gtmp1

    if lis.o.wparam == 1:
        del domveg
    print('824 create_vegtilespace done')


# createtiles_latlon
def createtiles_latlon():
    # liselev = [[0.0]*50 for i in range(50)]

    # localmask=[[0.0] * int(lis.d.lnc) for i in range(int(lis.d.lnr))]
    localmask = read_umdavhrr_mask()
    liselev = [[]]
    if lis.f.ecor > 0:
        liselev = [[0.0] * int(lis.d.lnc) for i in range(0, int(lis.d.lnr))]  # 循环变量在外面，表明建立的列表维数越靠前
        read_elev_gtopo30(liselev)

    # fgrd = np.zeros[int(lis.d.lnc), int(lis.d.lnr), int(lis.p.nt), dtype='f4']
    fgrd = read_umdavhrr_lc()
    # tsum = [[0.0]*lis.d.lnc for i in range(lis.d.lnr)]
    tsum = calculate_domveg(fgrd)
    create_vegtilespace(fgrd, tsum, localmask, liselev)

    del liselev
    del localmask
    print('846 createtiles_latlon')
    del fgrd
    del tsum


# time_manager
zd_timemag = {
    'uninit_int': -999999999,
    'rst_type': -999999999,  # cal# endar type
    'rst_nstep': -999999999,  # current step number
    'rst_step_days': -999999999,  # days component of timestep size
    'rst_step_sec': -999999999,  # seconds component of timestep size
    'rst_start_ymd': -999999999,  # start date
    'rst_start_tod': -999999999,  # start time of day
    'rst_stop_ymd': -999999999,  # stop date
    'rst_stop_tod': -999999999,  # stop time of day
    'rst_ref_ymd': -999999999,  # reference date
    'rst_ref_tod': -999999999,  # reference time of day
    'rst_curr_ymd': -999999999,  # current date
    'rst_curr_tod': -999999999  # current time of day
}


def date2time(a):
    global yrdays
    days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 30]
    if a[3] % 4 == 0 and a[3] % 100 != 0 or a[3] % 400 == 0:  # correct for y2k
        yrdays = 366
    else:
        yrdays = 365
    a[1] = 0
    for k in range(0, a[4] - 1):
        a[1] = a[1] + days[k]
    a[1] = a[1] + a[5]
    if yrdays == 366 and a[4] > 1:
        a[1] = a[1] + 1
    a[0] = (float(a[3]) + (
            (((((float(a[8]) / 60.0) + float(a[7])) / 60.0) + float(a[6])) / 24.0) + float(a[1] - 1)) / float(
        yrdays))
    a[2] = (((float(a[8]) / 60.0) + float(a[7])) / 60.0) + float(a[6])


def timemgr_init(lt):
    lt = t
    advacelist = [lt.etime, lt.edoy, lt.egmt, lt.eyr, lt.emo, lt.eda, lt.ehr, lt.emn, lt.ess]
    date2time(advacelist)
    lt.etime, lt.edoy, lt.egmt, lt.eyr, lt.emo, lt.eda, lt.ehr, lt.emn, lt.ess = advacelist[0], advacelist[1], \
                                                                                 advacelist[2], advacelist[3] \
        , advacelist[4], advacelist[5], advacelist[6], advacelist[7], advacelist[8]
    lt.tscount = 0
    timemgr_print(lt)


def timemgr_print(lt):
    lt = t
    if zd_spmd['masterproc'] == True:
        print(' 902 ************************************************')
        print("\033[1;36m%s\033[0m" % "Timestep size (seconds):", "\033[1;36m%s\033[0m" % lt.ts)
        print("\033[1;36m%s\033[0m" % "Start date (ymd tod):", "\033[1;36m%s\033[0m" % lt.syr,
              "\033[1;36m%s\033[0m" % lt.smo, "\033[1;36m%s\033[0m" % lt.sda)
        print("\033[1;36m%s\033[0m" % "Stop date (ymd tod):", "\033[1;36m%s\033[0m" % lt.eyr,
              "\033[1;36m%s\033[0m" % lt.emo, "\033[1;36m%s\033[0m" % lt.eda)
        print("\033[1;36m%s\033[0m" % "Current step number:", "\033[1;36m%s\033[0m" % lt.tscount)
        print("\033[1;36m%s\033[0m" % "Current date (ymd tod):", "\033[1;36m%s\033[0m" % lt.yr,
              "\033[1;36m%s\033[0m" % lt.mo, "\033[1;36m%s\033[0m" % lt.da)
        print(' 911 ************************************************')


def advance_timestep(lt):
    lt = t
    days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, lt.mo]
    lt.ss = lt.ss + lt.ts
    while lt.ss > 59:
        lt.ss = lt.ss - 60
        lt.mn = lt.mn + 1
    while lt.mn > 59:
        lt.mn = lt.mn - 60
        lt.hr = lt.hr + 1
    while lt.hr >= 24:
        lt.hr = lt.hr - 24
        lt.da = lt.da + 1
    if lt.yr % 4 == 0 and lt.yr % 100 != 0 or lt.yr % 400 == 0:
        days[2] = 29
    else:
        days[2] = 28
    tda = days[lt.mo]
    while lt.da > tda:
        lt.da = lt.da - days[lt.mo]
        lt.mo = lt.mo + 1
    while lt.mo > 12:
        lt.mo = lt.mo - 12
        lt.yr = lt.yr + 1
    if lt.tstype == 2:
        if lt.mo == 1 or lt.mo == 3 or lt.mo == 5 or lt.mo == 7 or lt.mo == 8 or lt.mo == 10 or lt.mo == 12:
            lt.ts = 2678400
            # noahdrv.WRITEINTN=744
        if lt.mo == 4 or lt.mo == 6 or lt.mo == 9 or lt.mo == 11:
            lt.ts = 2592000
        # noahdrv.WRITEINTN=720
        if lt.mo == 2:
            if lt.yr % 4 == 0 and lt.yr % 100 != 0 or lt.yr % 400 == 0:
                lt.ts = 2505600
            else:
                lt.ts = 2419200
    advacelist = [lt.time, lt.doy, lt.gmt, lt.yr, lt.mo, lt.da, lt.hr, lt.mn, lt.ss]
    date2time(advacelist)
    lt.time = advacelist[0]
    lt.doy = advacelist[1]
    lt.gmt = advacelist[2]
    lt.tscount = lt.tscount + 1
    if zd_spmd['masterproc'] == True:
        print('957 GSFC-LIS time: ', format(lt.mo, '>2.0f'), '/', format(lt.da, '>2.0f'), '/', format(lt.yr, '>4.0f'),
              ' ',
              format(lt.hr, '>2.0f'), ':', format(lt.mn, '>2.0f'), ':', format(lt.ss, '>2.0f'))
        # 24 format(a16,i2,a1,i2,a1,i4,1x,i2,a1,i2,a1,i2)#write 将后面内容以24格式写到显示器上。
    if lt.endcode == 0:  # end at #real-time date (tbd)
        print('961 warning: do not know how to stop in #real-time')
    if lt.endcode == 1:  # end on date specified in lis.crd file
        advacelist2 = [lt.etime, lt.edoy, lt.egmt, lt.eyr, lt.emo, lt.eda, lt.ehr, lt.emn, lt.ess]
        date2time(advacelist2)  # 两次调用该函数，两次改变下面三个变量的值。
        lt.etime = advacelist2[0]
        lt.edoy = advacelist2[1]
        lt.egmt = advacelist2[2]
        if lt.time >= lt.etime:
            lt.endtime = 1
            print('970 GSFC-LDAS run completed')


def get_curr_calday(lt, offset=""):
    lt = t
    days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    get_curr_calday = 0
    if lt.yr % 4 == 0 and lt.yr % 100 != 0 or lt.yr % 400 == 0:
        days[2] = 29
    else:
        days[2] = 28
    if lt.mo != 1:
        for i in range(0, lt.mo - 1):  # 从0开始，与fortran的下标1对应。
            get_curr_calday = get_curr_calday + days[i]
    get_curr_calday = get_curr_calday + lt.da + lt.hr / 24 + \
                      lt.mn / (24 * 60) + lt.ss / (24 * 60 * 60)
    if offset:
        if int(offset) > 0:
            get_curr_calday = get_curr_calday + int(offset) / (24 * 60 * 60)
        elif int(offset) < 0:
            get_curr_calday = get_curr_calday - int(offset) / (24 * 60 * 60)
    return get_curr_calday


def tick(tk):
    days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 31]
    tk[8] = tk[8] + tk[9]
    while tk[8] > 59:
        tk[8] = tk[8] - 60
        tk[7] = tk[7] + 1
    while tk[8] < 0:
        tk[8] = tk[8] + 60
        tk[7] = tk[7] - 1
    while tk[7] > 59:
        tk[7] = tk[7] - 60
        tk[6] = tk[6] + 1
    while tk[7] < 0:
        tk[7] = tk[7] + 60
        tk[6] = tk[6] - 1
    while tk[6] > 23:
        tk[6] = tk[6] - 24
        tk[5] = tk[5] + 1
    while tk[6] < 0:
        tk[6] = tk[6] + 24
        tk[5] = tk[5] - 1
    if tk[3] % 4 == 0 and tk[3] % 100 != 0 or tk[3] % 40 == 0:  # correct for y2k
        days[2] = 29
    else:
        days[2] = 28
    while tk[5] > days[tk[4]]:
        tk[5] = tk[5] - days[tk[4]]
        tk[4] = tk[4] + 1
    while tk[5] < 1:
        prvmo = tk[4] - 1
        if tk[4] == 1:
            prvmo = 12
        tk[5] = tk[5] + days[prvmo]
        if prvmo == 12:
            tk[4] = prvmo
            tk[3] = tk[3] - 1
        else:
            tk[4] = prvmo
    while tk[4] > 12:
        tk[4] = tk[4] - 12
        tk[3] = tk[3] + 1
    while tk[4] < 1:
        tk[4] = tk[4] + 12
        tk[3] = tk[3] - 1
    date2time(tk)  # 通过该函数进一步修改tk。


def time2date(t2d):
    days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 30]
    t2d[3] = int(t2d[0])
    tmp = int(t2d[0])
    if t2d[3] % 4 == 0 and t2d[3] % 100 != 0 or t2d[3] % 400 == 0:  # correct for y2k
        yrdays = 366
    else:
        yrdays = 365
    if yrdays == 366:
        days[2] = 29  #
    else:
        days[2] = 28
    t2d[1] = int((tmp - t2d[3]) * float(yrdays)) + 1
    tmp = int((tmp - t2d[3]) * float(yrdays)) + 1  # 默认以float替代Fortran中的双精度实型变量。
    t2d[6] = round((tmp - t2d[1]) * 24.0)
    tmp = round((tmp - t2d[1]) * 24.0)
    t2d[7] = int((tmp - t2d[6]) * 60.0)
    tmp = int((tmp - t2d[6]) * 60.0)
    ss = round((tmp - t2d[7]) * 60.0)
    t2d[4] = 1
    while t2d[1] > 0:
        t2d[1] = t2d[1] - days[t2d[4]]
        t2d[4] = t2d[4] + 1
    t2d[4] = t2d[4] - 1
    t2d[5] = t2d[1] + days[t2d[4]]
    t2d[2] = (((float(ss) / 60.0) + float(t2d[7])) / 60.0) + float(t2d[6])
    if t2d[2] == 24:
        t2d[2] = 0
        t2d[5] = t2d[5] + 1
        if t2d[5] > days[t2d[4]]:
            t2d[5] = 1
            t2d[4] = t2d[4] + 1
            if t2d[4] > 12:
                t2d[4] = 1
                t2d[3] = t2d[3] + 1


def get_step_size(lt):
    lt = t
    get_step_size = lt.ts
    return get_step_size


def get_nstep(lt):
    # lt = t
    get_nstep1 = lt.tscount
    return get_nstep1


def get_curr_date(gcd):
    gcd[0] = t
    gcd[1] = t.yr
    gcd[2] = t.mo
    gcd[3] = t.da
    gcd[4] = t.ss + t.mn * 60 + t.hr * 3600


def unpack_ymd_tod(uyt):  # ymd, tod, year, month, day, hour, min, sec
    uyt[2] = uyt[0] / 10000
    uyt[3] = (uyt[0] - uyt[2] * 10000) / 100
    uyt[4] = uyt[0] - uyt[2] * 10000 - uyt[3] * 100
    uyt[5] = uyt[1] / 10000
    uyt[6] = (uyt[1] - uyt[5] * 10000) / 100
    uyt[7] = uyt[1] - uyt[5] * 10000 - uyt[6] * 100


def pack_ymd_tod(pyt):  # ymd, tod, year, month, day, hour, min, sec
    pyt[0] = pyt[2] * 10000 + pyt[3] * 100 + pyt[4]
    pyt[1] = pyt[5] * 10000 + pyt[6] * 100 + pyt[7]


def convert_time_step(cts):  # days, sec, time_step
    if cts[2] >= 86400:
        cts[0] = cts[2] / 86400
        cts[1] = cts[2] - cts[0] * 86400
    else:
        cts[0] = 0
        cts[1] = cts[2]


def is_last_step(lt):
    lt = t
    is_last_step = False
    if lt.time >= lt.etime:
        is_last_step = True
    return is_last_step


def get_curr_time(days, seconds):
    return


def timemgr_write_restart(ftn_unit):
    # set_restart_values()
    RC = 0
    if RC != 0:
        print('1137 ERR: timemgr_write_restart -- WRITE iostat= ', RC, ' on i/o unit = ', ftn_unit)
        # endrun()


def timemgr_read_restart(ftn_unit):
    ftn_unit = input()
    RC = 0
    if RC != 0:
        print('1145 ERR: timemgr_read_restart -- READ iostat= ', RC, ' on i/o unit = ', ftn_unit)
        # endrun()


def chkrc(RC, mes):
    print(mes)
    # if ( RC == esmf_success ):
    # endrun()


def get_prev_date(yr, mon, day, tod):
    print('1156 get_prev_date not implemented, stopping..')
    # endrun()


def get_start_date(yr, mon, day, tod):
    print('1161 get_start_date not implemented, stopping..')
    # endrun()


def get_ref_date(yr, mon, day, tod):
    print('1166 get_ref_date not implemented, stopping..')
    # endrun()


# Readnoahcrd
def readnoahcrd(noahdrv):
    lsm = 0
    # noahdrv = noahdrvdec
    nml = f90nml.read('D:/NoahDate/NOAH/lis.crd')
    noahdrv.writeintn = nml['noah']['noahdrv']['WRITEINTN']
    noahdrv.noah_rfile = nml['noah']['noahdrv']['NOAH_RFILE']
    noahdrv.noah_mgfile = nml['noah']['noahdrv']['NOAH_MGFILE']
    noahdrv.noah_albfile = nml['noah']['noahdrv']['NOAH_ALBFILE']
    noahdrv.noah_tbot = nml['noah']['noahdrv']['NOAH_TBOT']
    noahdrv.noah_mxsnal = nml['noah']['noahdrv']['NOAH_MXSNAL']
    noahdrv.noah_vfile = nml['noah']['noahdrv']['NOAH_VFILE']
    noahdrv.noah_sfile = nml['noah']['noahdrv']['NOAH_SFILE']
    noahdrv.noah_ism = nml['noah']['noahdrv']['NOAH_ISM']
    noahdrv.noah_it = nml['noah']['noahdrv']['NOAH_IT']
    noahdrv.noah_nvegp = nml['noah']['noahdrv']['NOAH_NVEGP']
    noahdrv.noah_nsoilp = nml['noah']['noahdrv']['NOAH_NSOILP']
    for i in range(0, 6):
        noahdrv.noah_gridDesc[i] = nml['noah']['noahdrv']['noah_gridDesc'][i]
    print('1189 Running NOAH LSM:')
    print('1190 NOAH Active Restart File:', noahdrv.noah_rfile)
    noahdrv.noah_gfractime = 0.0
    noahdrv.noah_albtime = 0
    noahdrv.noah_albdchk = 0
    noahdrv.noah_gfracdchk = 0
    noahdrv.noahopen = 0
    noahdrv.noah_zst = 9


# readcard
def set_output_counters(soc):  # 列表soc=[numoutf]
    soc[0] = 0


def check_timestep(ct):  # ct=[timestep]
    CH = ' '
    if ct[0] > 3600:
        if zd_spmd['masterproc'] == True:
            print('--------------------------------------')
            print('1209 Warning, user timestep > 3600 seconds###')
            print('resetting lis.ts to 3600####')
            print('--------------------------------------')
    if ct[0] < 1:
        if zd_spmd['masterproc'] == True:
            print('1214 Timestep can not be less than 1 minute, reset to 15 min.')
        ct[0] = 15 * 60


def set_time(st):  # st=[time]
    # st = t
    st.yr = st.syr
    st.mo = st.smo
    st.da = st.sda
    st.hr = st.shr
    st.mn = st.smn
    st.ss = st.sss


def openfile(of):  # of = [name, odir, expcode, ffile]
    # global fbase,fcode,fname,fmkdir #声明全局变量，在if中使用。
    zfc1 = (of[1].rstrip()).lstrip()
    zfc2 = '/EXP' + str(of[2])  # + '/'.rstrip()
    listzfc2 = '/EXP' + str(of[2]) + '/'.rstrip()
    zfc3 = (of[3].rstrip()).lstrip()

    of[0] = zfc1 + zfc2  # + zfc3
    mkdir = 'mkdir -p ' + of[0]
    os.system(mkdir)
    return


def readcard():
    global name
    name = ''
    ttmp = ''
    # iargc = 0 #fortran 中用来获取命令行参数的数量，这里注释掉
    i = 0
    nml = f90nml.read('D:/NoahDate/NOAH/lis.crd')
    # 读取对应的结构体driver并赋值lis
    lis.d.domain = nml['driver']['LIS']['d']['DOMAIN']
    lis.d.lsm = nml['driver']['LIS']['d']['LSM']
    lis.f.force = nml['driver']['LIS']['f']['FORCE']
    lis.d.landcover = nml['driver']['LIS']['d']['LANDCOVER']
    lis.d.soil = nml['driver']['LIS']['d']['SOIL']
    lis.d.elev = nml['driver']['LIS']['d']['ELEV']
    lis.p.lai = nml['driver']['LIS']['p']['LAI']
    # 读取对应的结构体lis_run_inputs并赋值lis
    lis.o.expcode = nml['lis_run_inputs']['LIS']['o']['EXPCODE']
    lis.p.vclass = nml['lis_run_inputs']['LIS']['p']['VCLASS']
    lis.p.nt = nml['lis_run_inputs']['LIS']['p']['NT']
    lis.f.nf = nml['lis_run_inputs']['LIS']['f']['NF']
    lis.f.nmif = nml['lis_run_inputs']['LIS']['f']['NMIF']
    lis.f.ecor = nml['lis_run_inputs']['LIS']['f']['ECOR']
    lis.o.wfor = nml['lis_run_inputs']['LIS']['o']['WFOR']
    lis.f.interp = nml['lis_run_inputs']['LIS']['f']['INTERP']
    lis.o.wsingle = nml['lis_run_inputs']['LIS']['o']['WSINGLE']
    lis.o.wparam = nml['lis_run_inputs']['LIS']['o']['WPARAM']
    lis.o.wtil = nml['lis_run_inputs']['LIS']['o']['WTIL']
    lis.o.wout = nml['lis_run_inputs']['LIS']['o']['WOUT']
    lis.o.startcode = nml['lis_run_inputs']['LIS']['o']['STARTCODE']
    lis.t.sss = nml['lis_run_inputs']['LIS']['t']['SSS']
    lis.t.smn = nml['lis_run_inputs']['LIS']['t']['SMN']
    lis.t.shr = nml['lis_run_inputs']['LIS']['t']['SHR']
    lis.t.sda = nml['lis_run_inputs']['LIS']['t']['SDA']
    lis.t.smo = nml['lis_run_inputs']['LIS']['t']['SMO']
    lis.t.syr = nml['lis_run_inputs']['LIS']['t']['SYR']
    lis.t.endcode = nml['lis_run_inputs']['LIS']['t']['ENDCODE']
    lis.t.ess = nml['lis_run_inputs']['LIS']['t']['ESS']
    lis.t.emn = nml['lis_run_inputs']['LIS']['t']['EMN']
    lis.t.ehr = nml['lis_run_inputs']['LIS']['t']['EHR']
    lis.t.eda = nml['lis_run_inputs']['LIS']['t']['EDA']
    lis.t.emo = nml['lis_run_inputs']['LIS']['t']['EMO']
    lis.t.eyr = nml['lis_run_inputs']['LIS']['t']['EYR']
    lis.t.tstype = nml['lis_run_inputs']['LIS']['t']['TSTYPE']
    lis.t.ts = nml['lis_run_inputs']['LIS']['t']['TS']
    lis.d.udef = nml['lis_run_inputs']['LIS']['d']['UDEF']
    lis.o.odir = nml['lis_run_inputs']['LIS']['o']['ODIR']
    lis.o.dfile = nml['lis_run_inputs']['LIS']['o']['DFILE']
    lis.d.maxt = nml['lis_run_inputs']['LIS']['d']['MAXT']
    lis.d.mina = nml['lis_run_inputs']['LIS']['d']['MINA']
    # 读取对应的结构体proc_layout并赋值lis
    lis.d.npesx = nml['proc_layout']['lis']['d']['npesx']
    lis.d.npesy = nml['proc_layout']['lis']['d']['npesy']
    # 读取对应的结构体data_assimilation并赋值lis
    lis.a.daalg = nml['data_assimilation']['LIS']['a']['DAALG']
    lis.a.daobs = nml['data_assimilation']['LIS']['a']['DAOBS']
    lis.a.davar = nml['data_assimilation']['LIS']['a']['DAVAR']
    lis.a.nstv = nml['data_assimilation']['LIS']['a']['NSTV']
    lis.a.noutv = nml['data_assimilation']['LIS']['a']['NOUTV']
    lis.a.ndav = nml['data_assimilation']['LIS']['a']['NDAV']
    lis.a.nens = nml['data_assimilation']['LIS']['a']['NENS']
    lis.a.rensem = nml['data_assimilation']['LIS']['a']['RENSEM']
    lis.a.ndva = nml['data_assimilation']['LIS']['a']['NDVA']
    lis.a.obslo = nml['data_assimilation']['LIS']['a']['obslo']
    lis.a.obsla = nml['data_assimilation']['LIS']['a']['obsla']
    lis.a.obelo = nml['data_assimilation']['LIS']['a']['obelo']
    lis.a.obela = nml['data_assimilation']['LIS']['a']['obela']
    lis.a.oblor = nml['data_assimilation']['LIS']['a']['oblor']
    lis.a.oblar = nml['data_assimilation']['LIS']['a']['oblar']
    lis.a.fvlfn = nml['data_assimilation']['LIS']['a']['FVLFN']
    lis.a.svlfn = nml['data_assimilation']['LIS']['a']['SVLFN']
    lis.a.pSLGA = nml['data_assimilation']['LIS']['a']['pSLGA']
    lis.a.pFORCE = nml['data_assimilation']['LIS']['a']['pFORCE']
    # lis.a.win = nml['data_assimilation']['LIS']['a']['win']
    lis.d.ic = 1
    lis.d.ir = 1
    for i in range(1, len(sys.argv)):
        pre = sys.argv[i]  # 获取命令行特定位置的参数内容。
        if pre.rstrip() == '-ic':
            pre1 = sys.argv[i + 1]
            lis.d.ic = format(pre1, '>3.0f')
            ttmp = pre1.rstrip() + '-'  # 字符串拼接。
        if pre.rstrip() == '-ir':
            pre1 = sys.argv[i + 1]
            lis.d.ir = format(pre1, '>3.0f')
            ttmp = ttmp.rstrip() + pre1.rstrip()  # 字符串拼接。
    oplist = [name, lis.o.odir, lis.o.expcode, lis.o.dfile]
    # openfile(oplist)
    # name = oplist[0]
    # lis.o.odir=oplist[1]
    # lis.o.expcode=oplist[2]
    # lis.o.dfile=oplist[3]
    # 88 format(a4,25x,a3,5x,16a)  #Fortran中格式的定义，这里予以注释
    # 89 format(20x,a49)
    # ------------------------------------------------------------------------
    # Impose limit on time step
    # ------------------------------------------------------------------------
    ctlist = [lis.t.ts]
    check_timestep(ctlist)
    lis.t.ts = ctlist[0]

    set_time(lis.t)

    soclist = [lis.o.numoutf]
    set_output_counters(soclist)
    lis.o.numoutf = soclist[0]

    dtlist = [lis.t.time, lis.t.doy, lis.t.gmt, lis.t.yr, lis.t.mo, lis.t.da, lis.t.hr, lis.t.mn, lis.t.ss]
    date2time(dtlist)
    lis.t.time = dtlist[0]
    # lis.t.doy = dtlist[1]
    lis.t.gmt = dtlist[2]
    # lis.t.yr = dtlist[3]
    # lis.t.mo = dtlist[4]
    # lis.t.da = dtlist[5]
    # lis.t.hr = dtlist[6]
    # lis.t.mn = dtlist[7]
    # lis.t.ss = dtlist[8]

    if zd_spmd['masterproc'] == True:
        print('********** 1360 GSFC-LIS driver **********')
        print('experiment code: ', '-', lis.o.expcode, '-')
        print('starting time: ', lis.t.smo, '/', lis.t.sda, '/', lis.t.syr)
        print('# ending time: ', lis.t.emo, '/', lis.t.eda, '/', lis.t.eyr)
        print('  ')
        print('1365 forcing details:')
    # 读取对应的结构体landcover并赋值lis
    lis.p.mfile = nml['landcover']['LIS']['p']['MFILE']
    lis.p.vfile = nml['landcover']['LIS']['p']['VFILE']
    for i in range(0, 6):
        lis.d.lc_gridDesc[i] = nml['landcover']['LIS']['d']['lc_gridDesc'][i]
    # 读取对应的结构体elevation并赋值lis
    lis.p.elevfile = nml['elevation']['LIS']['p']['ELEVFILE']
    lis.p.slfile = nml['elevation']['LIS']['p']['SLFILE']
    for i in range(0, 6):
        lis.d.elev_gridDesc[i] = nml['elevation']['LIS']['d']['elev_gridDesc'][i]
    # 读取对应的结构体soils并赋值lis
    lis.p.safile = nml['soils']['LIS']['p']['SAFILE']
    lis.p.clfile = nml['soils']['LIS']['p']['CLFILE']
    lis.p.iscfile = nml['soils']['LIS']['p']['ISCFILE']
    lis.p.po1file = nml['soils']['LIS']['p']['PO1FILE']
    lis.p.po2file = nml['soils']['LIS']['p']['PO2FILE']
    lis.p.po3file = nml['soils']['LIS']['p']['PO3FILE']
    lis.p.sifile = nml['soils']['LIS']['p']['SIFILE']
    for i in range(0, 6):
        lis.d.soil_gridDesc[i] = nml['soils']['LIS']['d']['soil_gridDesc'][i]
    # 读取对应的结构体lai并赋值lis
    lis.p.avhrrdir = nml['lai']['LIS']['p']['AVHRRDIR']
    lis.p.modisdir = nml['lai']['LIS']['p']['MODISDIR']
    for i in range(0, 6):
        lis.p.lai_gridDesc[i] = nml['lai']['LIS']['p']['lai_gridDesc'][i]

    lis.p.laitime = 0.0
    lis.p.saitime = 0.0
    if lis.p.lai == 2:
        if zd_spmd['masterproc']:
            print("1396 Using AVHRR Satellite LAI")
    if lis.p.lai == 3:
        if zd_spmd['masterproc']:
            print("1399 Using MODIS Satellite LAI")
    if lis.a.daalg == 0:
        if zd_spmd['masterproc']:
            print("1402 No Data assimilation")
        lis.a.nens = 1
    lis.f.rstflag = 1
    lis.f.gridchange = 1
    lis.o.foropen = 0
    if lis.d.maxt > lis.p.nt:
        lis.d.maxt = lis.p.nt
    if lis.d.maxt < 1:
        lis.d.maxt = 1
    if zd_spmd['masterproc']:
        print('1412 miscellaneous details:')
        if lis.p.vclass == 1:
            print('1414 avhrr umd vegetation', lis.p.vclass)
        if lis.p.vclass == 2:
            print('1416 modis umd vegetation', lis.p.vclass)
        print('1417 mask file:', lis.p.mfile)
        print('1418 vegetation file: ', lis.p.vfile)
        if lis.d.soil == 1:
            print('1420 original vegetation-based soil scheme')
        if lis.d.soil == 2:
            print('1422 reynolds soils')
        if lis.d.soil == 1:
            print('1424 original vegetation-based soil scheme')
        if lis.d.soil == 2:
            print('1426 reynolds soils')


# Lisdrv_module
zd_lisdrv = {
    'grid': [],
    'tile': [],
    'gindex': [[]],
    'ntiles_pergrid': [],
    's_tid': []
}


def setup_timeMgr():
    timemgr_init(lis.t)
    if zd_spmd['masterproc'] == True:
        print('1442 time manager initialized..')


def LIS_config():
    spmd_init()
    readcard()
    setup_timeMgr()
    lis.t.endtime = 0


def LIS_ticktime():
    curSec = 0
    ierr = 0
    advance_timestep(lis.t)


def LIS_domain_init():
    readdomain_default()
    createtiles_latlon()
    print('1461 lisdrv_module')
    zd_spmd['nchs'][zd_spmd['iam']] = lis.d.nch
    zd_spmd['tdeltas'][zd_spmd['iam']] = lis.d.nch
    zd_spmd['ngrids'][zd_spmd['iam']] = lis.d.ngrid
    zd_spmd['gdeltas'][zd_spmd['iam']] = lis.d.ngrid
    lis.d.glbnch = 0
    for i in range(0, zd_spmd['npes']):
        lis.d.glbnch = lis.d.glbnch + zd_spmd['nchs'][i]
    lis.d.glbngrid = 0
    for i in range(0, zd_spmd['npes']):
        lis.d.glbngrid = lis.d.glbngrid + zd_spmd['ngrids'][i]
    if zd_spmd['masterproc']:
        zd_spmd['toffsets'][0] = 0
        for i in range(1, zd_spmd['npes']):
            zd_spmd['toffsets'][i] = zd_spmd['toffsets'][i - 1] + zd_spmd['tdeltas'][i - 1]
        zd_spmd['goffsets'][0] = 0
        for i in range(1, zd_spmd['npes']):
            zd_spmd['goffsets'][i] = zd_spmd['goffsets'][i - 1] + zd_spmd['gdeltas'][i - 1]


def setnch(nc, nr, maxt):
    lis.d.glbnch = nc * nr * maxt  # may be too big


def getdomain():
    d = lis.d.domain
    return d


def getlsm():
    lsmno = lis.d.lsm
    return lsmno


def getnch():  # result(n)  给返回值换一个名字
    n = lis.d.glbnch
    return n


def getnc():
    ncol = lis.d.lnc
    return ncol


def getnr():  # result(nrow)
    nrow = lis.d.lnr
    return nrow


def gettileindex(lat, lon):  # result(k)
    k = -1
    for t in range(0, lis.d.nch):
        if zd_lisdrv['grid'][zd_lisdrv['tile'][t].index - 1].lat == lat and zd_lisdrv['grid'][
            zd_lisdrv['tile'][t].index - 1].lon == lon:
            k = t
    return k


def getmaxt():  # result(maxt)
    maxt = lis.d.maxt
    return maxt


def getforcing():  # result(f)
    f = lis.f.force
    return f


def getnmif():  # result(f)
    f = lis.f.nmif
    return f


def LIS_endofrun():  # result(finish)
    if zd_spmd['masterproc']:
        finish = is_last_step(lis.t)
    return finish


def LIS_backuptime(on):
    if on == 0:  # backup
        lis.tt = lis.t
    else:  # resume
        lis.t = lis.tt
    lis.f.findtime1 = 1
    lis.f.findtime2 = 1


def LISdrv_clear():
    zd_lisdrv['grid'] = []
    zd_lisdrv['gindex'] = []
    zd_lisdrv['tile'] = []
    zd_lisdrv['s_tid'] = []
    zd_lisdrv['ntiles_pergrid'] = []


# Noah_varder
zd_varder = {
    # 'noah':[noahdec for x in range(100)],
    # 'noahbk':[noahdec for x in range(100)]
    'noah': [],
    'noahbk': []
}


def noah_varder_ini(nch):
    readnoahcrd(noahdrv)
    zd_varder['noah'] = [noahdec() for x in range(0, nch)]
    if lis.a.daalg == 2:
        zd_varder['noahbk'] = [0 for x in range(0, nch)]
    if lis.a.daalg == 3:
        zd_varder['noahbk'] = [0 for x in range(0, nch)]
    if lis.a.daalg == 4:
        zd_varder['noahbk'] = [0 for x in range(0, nch)]
    if lis.a.daalg == 5:
        zd_varder['noahbk'] = [0 for x in range(0, nch)]
    if lis.o.wout == 4 or lis.o.wout == 5:
        zd_varder['noahbk'] = [0 for x in range(0, nch)]


def noah_clear():
    zd_varder['noah'] = []
    if lis.a.daalg == 2 or lis.a.daalg == 3 or lis.a.daalg == 4 or lis.a.daalg == 5:
        zd_varder['noahbk'] = []


# Noah_reset
def noah_reset(flag):
    if flag == 0:
        zd_varder['noahbk'] = zd_varder['noah']
    else:
        zd_varder['noah'] = zd_varder['noahbk']


# Mapvegc
def mapvegc(vegt):
    # global sibveg
    sibveg = int(0)
    if vegt[0] == 1:
        sibveg = 4
    if vegt[0] == 2:
        sibveg = 1
    if vegt[0] == 3:
        sibveg = 5
    if vegt[0] == 4:
        sibveg = 2
    if vegt[0] == 5:
        sibveg = 3
    if vegt[0] == 6:
        sibveg = 3
    if vegt[0] == 7:
        sibveg = 6
    if vegt[0] == 8:
        sibveg = 8
    if vegt[0] == 9:
        sibveg = 9
    if vegt[0] == 10:
        sibveg = 7
    if vegt[0] == 11:
        sibveg = 12
    if vegt[0] == 12:
        sibveg = 11
    if vegt[0] == 13:
        sibveg = 11
    if vegt[0] == 13:
        sibveg = 7
    vegt[0] = sibveg


# Noah_setvegparms
def noah_setvegparms():
    print('1635 MSG: setnoahp -- Calling MAPVEGC to convert UMD to SIB', zd_spmd['iam'])
    print('1636 DBG: setnoahp -- nch', lis.d.nch, zd_spmd['iam'])
    for n in range(0, lis.d.nch):
        maplist = [zd_lisdrv['tile'][n].vegt]
        mapvegc(maplist)  # Fortran中将结构体赋值给动态数组后，该调用方式在Python中是否正确？。L
        zd_lisdrv['tile'][n].vegt = maplist[0]
        zd_varder['noah'][n].vegt = zd_lisdrv['tile'][n].vegt  # Fortran中将结构体赋值给动态数组后，该调用方式在Python中是否正确？。L
    print('DBG: setnoahp -- left MAPVEGC', zd_spmd['iam'])
    f = open(noahdrv.noah_vfile, 'r')
    content = f.readlines()
    value = []
    for i in range(0, noahdrv.noah_nvegp):
        column = content[i].replace(" ", ",").split(',')[0:lis.p.nt]
        value.append(column)
    f.close()
    for i in range(0, lis.d.nch):
        for j in range(0, noahdrv.noah_nvegp):
            zd_varder['noah'][i].VEGP[j] = value[j][zd_lisdrv['tile'][i].vegt]


# Noah_getlsmstate
def noah_getlsmstate(StVarMN):
    # StVarMN = [[0.0 for x in range(lis.d.nch)] for x in range(lis.a.ndav+lis.a.ndva)]
    for i in range(0, lis.d.nch):  # python lis.d.nch为列，后面为行。跟fortran数组排法相同。
        StVarMN[0][i] = zd_varder['noah'][i].qle  # 将第一行的所有元素赋值。
        StVarMN[1][i] = zd_varder['noah'][i].T1
        StVarMN[2][i] = zd_varder['noah'][i].STC[0]
        StVarMN[3][i] = zd_varder['noah'][i].SMC[0]
        StVarMN[4][i] = zd_varder['noah'][i].CMC


# Noah_setlsmstate
def noah_setlsmstate(StVarMN):
    for i in range(0, lis.d.nch):  # python lis.d.nch为列，后面为行。跟fortran数组排法相同。
        zd_varder['noah'][i].T1 = StVarMN[1][i]
        zd_varder['noah'][i].STC[0] = StVarMN[2][i]
        zd_varder['noah'][i].SMC[0] = StVarMN[3][i]
        zd_varder['noah'][i].CMC = StVarMN[4][i]


# create_gdd
zd_save = {
    'flag': 0  # 用字典的形式来替代fortran中的save属性。在下次调用函数是保留上次函数的返回值
}  # 且该字典仅供randou函数调用。


def randou(rand_array, number):
    for i in range(1, number + 1):
        if zd_save['flag'] == 0:
            zd_save['flag'] = 1
        ran1 = random.random()
        rand_array[i - 1] = ran1


def GAMMA(a):
    pi = 3.14159265358979324
    cg = 7
    p = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
         771.32342877765313, -176.61502916214059, 12.507343278686905,
         -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7]
    x = a
    if x < 0.5:
        g = pi / (math.sin(pi * x) * GAMMA(1.0 - x))
    else:
        x = x - 1.0
        t = p[0]
        for i in range(1, cg + 2):
            t = t + p[i] / (x + i)
        w = x + cg + 0.5
        g = math.sqrt(2.0 * pi) * pow(w, x + 0.5) * math.exp(-w) * t
    return g


def ggd_ensemble(ggdlist):  # ggdlist = [xx_array,numberGGD,aa,maxstd] 仅第一个元素发生变化
    import numpy as np
    index = 1
    rs = 0.0
    num1, sum = 0, 0.0
    group = np.array([0.0, 0.0])
    group1 = []
    numberGGD = int(ggdlist[1])
    aa = float(ggdlist[2])
    maxstd = float(ggdlist[3])
    aa_array = [0.0 for i in range(numberGGD)]
    yy_array = [0.0 for i in range(numberGGD)]
    s = 1 / aa - math.floor(1 / aa)
    v = math.exp(1.0) / (math.exp(1.0) + s)
    gama = GAMMA(3 / aa) / GAMMA(1 / aa)
    while index < numberGGD:
        rs = v
        while rs == v:
            rs = random.random()
        group[0] = rs
        rs1 = random.random()
        group[1] = rs1
        if group[0] < v:
            x1 = pow(group[0] / v, 1 / s)
            g1 = group[1] * pow(x1, s - 1)
        else:
            x1 = 1 - math.log((group[0] - v) / (1 - v))
            g1 = group[1] * math.exp(0 - x1)
        if g1 > pow(x1, s - 1) * math.exp(0 - x1):
            group = [0 for i in range(2)]
            continue
        else:
            num1 = math.floor(1 / aa)
            sum = 0.0
            if num1 >= 1:
                group1 = [0 for i in range(num1)]
                randou(group1, num1)
                for kk in range(1, num1 + 1):
                    sum = sum + math.log(group1[kk - 1])
            x1 = x1 - sum
            x1 = x1 / pow(pow(gama, 0.5), aa)
            x1 = pow(x1, 1 / aa)
            if x1 <= maxstd:
                aa_array[index - 1] = x1
                aa_array[index] = 0 - x1
                index = index + 2
            if group1 != []:
                del group1
    ggdlist[0] = aa_array
    yy_array = []
    aa_array = []


def brinv(a, n, l):  # l在调用时直接给实参赋值1，参数不再包装，a为列表传值，n不变。
    IS = [0.0 for i in range(n)]
    JS = [0.0 for i in range(n)]
    l = 1
    t, d = 0.0, 0.0
    for k in range(1, n + 1):
        d = 0.0
        for i in range(k, n + 1):
            for j in range(k, n + 1):
                if abs(a[j - 1][i - 1]) > d:
                    d = abs(a[j - 1][i - 1])
                    IS[k - 1] = i
                    JS[k - 1] = j
        if d + 1.0 == 1.0:
            l = 0
            print(' ' + '1776 ERR* * NOT INV')  # format 20的格式。
        for j in range(1, n + 1):
            t = a[j - 1][k - 1]
            a[j - 1][k - 1] = a[j - 1][IS[k - 1] - 1]
            a[j - 1][IS[k - 1] - 1] = t
        for i in range(1, n + 1):
            t = a[k - 1][i - 1]
            a[k - 1][i - 1] = a[JS[k - 1] - 1][i - 1]
            a[JS[k - 1] - 1][i - 1] = t
        a[k - 1][k - 1] = 1 / a[k - 1][k - 1]
        for j in range(1, n + 1):
            if j != k:
                a[j - 1][k - 1] = a[j - 1][k - 1] * a[k - 1][k - 1]
        for i in range(1, n + 1):
            if i != k:
                for j in range(1, n + 1):
                    if j != k:
                        a[j - 1][i - 1] = a[j - 1][i - 1] - a[k - 1][i - 1] * a[j - 1][k - 1]
        for i in range(1, n + 1):
            if i != k:
                a[k - 1][i - 1] = -a[k - 1][i - 1] * a[k - 1][k - 1]
    for k in range(n, 0, -1):
        for j in range(1, n + 1):
            t = a[j - 1][k - 1]
            a[j - 1][k - 1] = a[j - 1][JS[k - 1] - 1]
            a[j - 1][JS[k - 1] - 1] = t
        for i in range(1, n + 1):
            t = a[k - 1][i - 1]
            a[k - 1][i - 1] = a[IS[k - 1] - 1][i - 1]
            a[IS[k - 1] - 1][i - 1] = t
    IS = []
    JS = []


def brmul(a, b, m, n, k, c):  # 对参数c列表进行更新。
    for i in range(1, m + 1):
        for j in range(1, k + 1):
            c[j - 1][i - 1] = 0.0
            for l in range(1, n + 1):
                c[j - 1][i - 1] = c[j - 1][i - 1] + a[l - 1][i - 1] * b[j - 1][l - 1]


# evap_limits_module
zd_evap = {
    'aa': [],
    'std': [],
    'minstd': [],
    'maxstd': [],
    'minval': [],
    'maxval': [],
    'irs': 0
}


# Noah_qa_sm
def noah_qa_sm(nqm):  # nqm = [sm,tt]
    if nqm[0] <= float(zd_evap['minstd'][nqm[1]]):
        nqm[0] = float(zd_evap['minstd'][nqm[1]])
    if nqm[0] >= float(zd_evap['maxstd'][nqm[1]]):
        nqm[0] = float(zd_evap['maxstd'][nqm[1]])


# Lis_indices_module
zd_indices = {
    'lis_nc_working': 0,
    'lis_nr_working': 0,
    'lis_nc_data': 0,
    'lis_nr_data': 0,
    'lis_tnroffset': 0,
    'lis_grid_offset': 0,
    'lis_g2l_row_offset': 0,
    'lis_g2l_col_offset': 0
}


def lis_global_to_local_row_offset(offset):
    lis_get_run_slat = lis.d.gridDesc[3]
    lis_get_data_slat = lis.d.gridDesc[43]
    lis_get_run_lat_res = lis.d.gridDesc[8]  # 列表元素，初始化均为0，要注意该列表是否有赋值或者写入情况。L
    lis_global_to_local_row_offset = round((lis_get_run_slat - lis_get_data_slat) / lis_get_run_lat_res) - offset
    return lis_global_to_local_row_offset


def lis_global_to_local_col_offset():
    lis_get_run_wlon = lis.d.gridDesc[4]
    lis_get_data_wlon = lis.d.gridDesc[44]
    lis_get_run_lon_res = lis.d.gridDesc[9]  # 列表元素，初始化均为0，要注意该列表是否有赋值或者写入情况。L
    lis_global_to_local_col_offset = round((lis_get_run_wlon - lis_get_data_wlon) / lis_get_run_lon_res)
    return lis_global_to_local_col_offset


def lis_set_indices():
    zd_indices['lis_nc_working'] = lis.d.lnc
    zd_indices['lis_nr_working'] = lis.d.lnr
    zd_indices['lis_tnroffset'] = 0
    zd_indices['lis_grid_offset'] = 0
    zd_indices['lis_nc_data'] = lis.d.lnc
    zd_indices['lis_nr_data'] = lis.d.lnr
    zd_indices['lis_g2l_row_offset'] = lis_global_to_local_row_offset(zd_indices['lis_tnroffset'])
    zd_indices['lis_g2l_col_offset'] = lis_global_to_local_col_offset()
    if lis.o.plevel == 2:
        print('1877 DBG: lis_set_indices -- lis_nc_working', zd_indices['lis_nc_working'])
        print('1878 DBG: lis_set_indices -- lis_nr_working', zd_indices['lis_nr_working'])
        print('1879 DBG: lis_set_indices -- lis%d%lnc', lis.d.lnc)
        print('1880 DBG: lis_set_indices -- lis%d%lnr', lis.d.lnr)
        print('1881 DBG: lis_set_indices -- lis_nc_data', zd_indices['lis_nc_data'])
        print('1882 DBG: lis_set_indices -- lis_nr_data', zd_indices['lis_nr_data'])
        print('1883 DBG: lis_set_indices -- lis%d%gnc', lis.d.gnc)
        print('1884 DBG: lis_set_indices -- lis%d%gnr', lis.d.gnr)
        print('1885 DBG: lis_set_indices -- lis_tnroffset', zd_indices['lis_tnroffset'])
        print('1886 DBG: lis_set_indices -- lis_grid_offset', zd_indices['lis_grid_offset'])


# Baseforcing_module
zd_basefor_mod = {
    'glbdata1': [[]],
    'glbdata2': [[]],
    'modelelev': []
}


def LIS_baseforcing_init():
    gridDesci = [0 for x in range(0, 50)]
    allocate_forcing_mem()


def get():
    return


def time_interp():
    return


def allocate_forcing_mem():
    nmif = getnmif()
    zd_basefor_mod['glbdata1'] = [[0.0] * nmif for i in range(lis.d.ngrid)]
    zd_basefor_mod['glbdata2'] = [[0.0] * nmif for i in range(lis.d.ngrid)]


def baseforcing_clear():
    zd_basefor_mod['glbdata1'] = []
    zd_basefor_mod['glbdata2'] = []


# evapda_module
class evapdadec:
    def __init__(self):
        self.sig = [[0] * 10 for i in range(0, 10)]
        self.stv = [0 for x in range(0, 10)]


evapdadec = evapdadec()


# evapobs_module
class evapobsdec:
    def __init__(self):
        self.evap = [0 for x in range(0, 1)]
        self.mcv = [0 for x in range(0, 1)]


evapobsdec = evapobsdec()


# noah_SMdaOUT
def noah_SMdaOUT(SMStVar):
    for i in range(0, lis.d.nch):
        zd_varder['noah'][i].lwnet = zd_varder['noah'][i].lwnet + (5.67e-8) * (
                pow(SMStVar[1][i], 4.0) - pow(zd_varder['noah'][i].T1, 4.0))
        zd_varder['noah'][i].T1 = SMStVar[1][i]
        zd_varder['noah'][i].STC[0] = SMStVar[2][i]
        zd_varder['noah'][i].rootmoist = zd_varder['noah'][i].rootmoist + (
                SMStVar[3][i] - zd_varder['noah'][i].SMC[0]) * 0.1 * 1000.0
        zd_varder['noah'][i].SMC[0] = SMStVar[3][i]
        zd_varder['noah'][i].soilmoist1 = zd_varder['noah'][i].SMC[0] * 1000.0 * 0.1
        zd_varder['noah'][i].CMC = SMStVar[4][i]
        zd_varder['noah'][i].canopint = zd_varder['noah'][i].CMC * 1000.0


# Evapobsdrv_module
zd_evapob = {
    'evapobsdir': '',
    'nob': 0,
    'readfl': '',
    'oer': 0.0,
    'mer': 0.0,
    'merl': 0.0,
    'evapobs': [[evapobsdec] * 10 for i in range(10)]  # evapobs为二维数组
}


def evapobs_setup():
    print('1969 reading namelist soil_moisture_da')
    nml = f90nml.read('D:/NoahDate/NOAH/lis.crd')
    # 打开该文件。读取namelist中soil_moisture_da中的变量。L
    zd_evapob['readfl'] = 'readfl'
    zd_evapob['oer'] = nml['soil_moisture_da']['oer']
    zd_evapob['evapobsdir'] = nml['soil_moisture_da']['evapobsdir']
    zd_evapob['nob'] = nml['soil_moisture_da']['nob']
    zd_evapob['mer'] = nml['soil_moisture_da']['mer']
    zd_evapob['merl'] = nml['soil_moisture_da']['merl']


# Evapenkfdrv_module
zd_evapenk = {
    'nst': 0,
    'nou': 0,
    'nda': 0,
    'nen': 0,
    'ntl': 0,
    'ndv': 0,
    'win': 0,
    'rs': [-99 for x in range(0, 1)],
    'StVarEvap': [[[0 for i in range(10)] for i in range(10)] for i in range(10)],
    'StVarMN': [[0] * 10 for i in range(10)],
    'Ensem4D': [[[[0 for i in range(10)] for i in range(10)] for i in range(10)] for i in range(10)],
    'smda': [evapdadec for i in range(10)]  # 指针数组。
}


def evap_enkf_init():
    import numpy as np
    zd_evapenk['nst'] = lis.a.nstv  # number of lsm state variables
    zd_evapenk['nou'] = lis.a.noutv  # number of lsm output variables
    zd_evapenk['nda'] = lis.a.ndav  # number of lsm state variables to be assimilated
    zd_evapenk['ndv'] = lis.a.ndva  # number of model diagnostic variable to be assimilated
    zd_evapenk['nen'] = lis.a.nens  # number of ensemble members for EnKF
    zd_evapenk['ntl'] = lis.d.nch  # number of tiles of the domain
    zd_evapenk['win'] = lis.a.win
    zd_evapenk['StVarEvap'] = [
        [[0 for i in range(zd_evapenk['ntl'])] for i in range(zd_evapenk['nda'] + zd_evapenk['ndv'])] for i in
        range(zd_evapenk['nen'])]
    print('2009 evapken', zd_evapenk['ntl'], zd_evapenk['ndv'], zd_evapenk['ndv'])
    zd_evapenk['StVarMN'] = [[0] * zd_evapenk['ntl'] for i in range(zd_evapenk['nda'] + zd_evapenk['ndv'])]
    zd_evapenk['Ensem4D'] = [[[[0 for i in range(lis.a.win)] for i in range(zd_evapenk['ntl'])] for i in
                              range(zd_evapenk['nda'] + zd_evapenk['ndv'])] for i in range(lis.a.nens)]
    print('2013 Setups for state variable limits in ', 'E:\data\SVLFN')  # 注意该位置要修改。
    f = open(lis.a.svlfn, 'r')  # 注意该位置要修改。

    for i in range(0, zd_evapenk['nda'] + zd_evapenk['ndv']):
        data = f.readline()
        data_current = data.strip('\n').split(' ')
        zd_evap['aa'].append(data_current[0])
        zd_evap['std'].append(data_current[1])
        zd_evap['minstd'].append(data_current[2])
        zd_evap['maxstd'].append(data_current[3])
        zd_evap['minval'].append(data_current[4])
        zd_evap['maxval'].append(data_current[5])
        print('1835', i, zd_evap['aa'][i], zd_evap['std'][i], zd_evap['minstd'][i], zd_evap['maxstd'][i],
              zd_evap['minval'][i], zd_evap['maxval'][i])
    f.close()
    print('2028', np.array(zd_evapenk['StVarEvap']).shape)
    print('2029', zd_evap)
    zd_evapenk['StVarEvap'][0][0][0] = -9999


def evap_enkf_setst(nenth):
    noah_reset(1)
    for j in range(0, zd_evapenk['nda'] + zd_evapenk['ndv']):
        for i in range(0, zd_evapenk['ntl']):
            zd_evapenk['StVarMN'][j][i] = zd_evapenk['StVarEvap'][nenth - 1][j][i]
    noah_setlsmstate(zd_evapenk['StVarMN'])


def evap_enkf_update():
    ir = 0
    obsst = [0.0 for i in range(zd_evapenk['ndv'])]
    RR = [[[0 for i in range(zd_evapenk['win'])] for i in range(zd_evapenk['ndv'])] for i in range(zd_evapenk['ndv'])]
    Mhxb = [0 for i in range(zd_evapenk['ndv'])]
    Hxb = [[0] * zd_evapenk['ndv'] for i in range(zd_evapenk['nen'])]
    tmpStVar = [[0] * zd_evapenk['nen'] for i in range(zd_evapenk['ndv'])]
    Tyb = [[0] * zd_evapenk['nen'] for i in range(zd_evapenk['ndv'])]
    obssm = [[0] * zd_evapenk['ndv'] for i in range(zd_evapenk['nen'])]
    Yhxb = [[0] * zd_evapenk['ndv'] for i in range(zd_evapenk['nen'])]
    Pao = [[0] * zd_evapenk['nen'] for i in range(zd_evapenk['nen'])]
    ObsR = [[0] * zd_evapenk['ndv'] for i in range(zd_evapenk['ndv'])]
    PhO = [[0] * zd_evapenk['nen'] for i in range(zd_evapenk['nen'])]
    Xb = [[0] * zd_evapenk['nda'] for i in range(zd_evapenk['nen'])]
    Ph = [[0] * zd_evapenk['nen'] for i in range(zd_evapenk['nen'])]
    Pa = [[0] * zd_evapenk['nen'] for i in range(zd_evapenk['nen'])]
    Mxb = [[0] * zd_evapenk['nda'] for i in range(zd_evapenk['nen'])]
    CYM = [0 for i in range(zd_evapenk['ndv'])]
    CYO = [0 for i in range(zd_evapenk['nen'])]
    maindir1 = '/mnt/data1/chensh/00000000/00000000000000000000'
    zd_evapenk['StVarMN'] = [[0] * 10 for i in range(10)]
    for k in range(0, zd_evapenk['nen']):
        for j in range(0, zd_evapenk['nda'] + zd_evapenk['ndv']):
            for i in range(0, zd_evapenk['ntl']):
                zd_evapenk['StVarMN'][j][i] = zd_evapenk['StVarMN'][j][i] + zd_evapenk['StVarEvap'][k][j][i]
    for j in range(0, zd_evapenk['nda'] + zd_evapenk['ndv']):
        for i in range(0, zd_evapenk['ntl']):
            zd_evapenk['StVarMN'][j][i] = zd_evapenk['StVarMN'][j][i] / zd_evapenk['nen']
    for t in range(0, zd_evapenk['ntl']):
        daflag = 1
        for winn in range(0, zd_evapenk['win']):
            for j in range(0, zd_evapenk['ndv']):
                ggd_list = [CYO, zd_evapenk['nen'], zd_evap['aa'][j], zd_evap['maxstd'][j]]
                ggd_ensemble(ggd_list)
                CYO = ggd_list[0]
                obsst[j] = zd_evapob['evapobs'][t][winn].evap[j]
                for i in range(0, zd_evapenk['nen']):
                    obssm[i][j] = obsst[j] + CYO[i]  # sqrt(abs(CV))*gasdev(rs)
                    nqmlist = [obssm[i][j], j]
                    noah_qa_sm(nqmlist)
                    obssm[i][j] = nqmlist[0]
                    CYM[j] = CYM[j] + obssm[i][j]
            for x in range(0, zd_evapenk['ndv']):
                CYM[x] = CYM[x] / zd_evapenk['nen']
            for j in range(0, zd_evapenk['ndv']):
                for k in range(0, zd_evapenk['ndv']):
                    for i in range(0, zd_evapenk['nen']):
                        ObsR[k][j] = ObsR[k][j] + (obssm[i][j] - CYM[j]) * (obssm[i][k] - CYM[k])
                    RR[k][j][winn] = ObsR[k][j] / (zd_evapenk['nen'])  # /(zd_evapenk['nen']-1)
        for winn in range(0, zd_evapenk['win']):
            for k in range(0, zd_evapenk['nen']):
                for i in range(0, zd_evapenk['ndv']):
                    Hxb[k][i] = zd_evapenk['Ensem4D'][k][0:zd_evapenk['ndv']][t][winn]
                    Tyb[k][i] = zd_evapenk['Ensem4D'][k][0:zd_evapenk['ndv']][t][winn]
                    Mhxb[i] = Mhxb[i] + Hxb[k][i]
            # Mhxb[:]=Mhxb[:]/zd_evapenk['nen']
            for i in range(0, zd_evapenk['ndv']):
                Mhxb[i] = Mhxb[i] / zd_evapenk['nen']
            for k in range(1, zd_evapenk['nen']):
                for i in range(0, zd_evapenk['ndv']):
                    Yhxb[k][i] = obsst[i] - Hxb[k][i]
                    Hxb[k][i] = Hxb[k][i] - Mhxb[i]
                    Tyb[i][k] = Tyb[i][k] - Mhxb[i]
            for i in range(0, zd_evapenk['ndv']):
                for j in range(0, zd_evapenk['ndv']):
                    ObsR[i][j] = RR[i][j][winn]
            brinv(ObsR, zd_evapenk['ndv'], ir)  # invert matrix
            ir = 1
            brmul(Tyb, ObsR, zd_evapenk['nen'], zd_evapenk['ndv'], zd_evapenk['ndv'], tmpStVar)
            brmul(tmpStVar, Hxb, zd_evapenk['nen'], zd_evapenk['ndv'], zd_evapenk['nen'], Pao)
            brmul(tmpStVar, Yhxb, zd_evapenk['nen'], zd_evapenk['ndv'], zd_evapenk['nen'], PhO)
            for i in range(0, zd_evapenk['nen']):
                for j in range(0, zd_evapenk['nen']):
                    Pa[i][j] = Pao[i][j] + Pa[i][j]
                    Ph[i][j] = PhO[i][j] + Ph[i][j]
        for k in range(1, zd_evapenk['nen']):
            Pa[k - 1][k - 1] = Pa[k - 1][k - 1] + zd_evapenk['nen'] - 1
            for i in range(0, zd_evapenk['nda']):
                Xb[k - 1][i] = \
                    zd_evapenk['StVarEvap'][k - 1][zd_evapenk['ndv']:zd_evapenk['ndv'] + zd_evapenk['nda'] - 1][t] - \
                    zd_evapenk['StVarMN'][zd_evapenk['ndv']:zd_evapenk['ndv'] + zd_evapenk['nda'] - 1][t]
        brinv(Pa, zd_evapenk['nen'], ir)
        brmul(Pa, Ph, zd_evapenk['nen'], zd_evapenk['nen'], zd_evapenk['nen'], PhO)
        brmul(Xb, PhO, zd_evapenk['nda'], zd_evapenk['nen'], zd_evapenk['nen'], Mxb)
        zd_evapenk['StVarMN'] = [[0] * zd_evapenk['ntl'] for i in range(zd_evapenk['ndv'] + zd_evapenk['nda'])]
        for k in range(1, zd_evapenk['nen'] + 1):
            for i in range(zd_evapenk['ndv'], zd_evapenk['ndv'] + zd_evapenk['nda']):
                for j in range(0, zd_evapenk['nda']):
                    zd_evapenk['StVarMN'][i][t] = zd_evapenk['StVarMN'][i][t] + Mxb[k - 1][j]
        for k in range(zd_evapenk['ndv'], zd_evapenk['ndv'] + zd_evapenk['nda']):
            zd_evapenk['StVarMN'][k][t] = zd_evapenk['StVarMN'][k][t] / zd_evapenk['nen']
        # for k in range(zd_evapenk['ndv'],zd_evapenk['ndv']+zd_evapenk['nda']):
        # noah_qa_sm(obssm[i][j],j) #该列表i和j从哪里来？弄清之后在调用。
    noah_reset(1)
    noah_SMdaOUT(zd_evapenk['StVarMN'])
    Mhxb = []
    Hxb = []
    Tyb = []
    obssm = []
    Yhxb = []
    Pao = []
    ObsR = []
    PhO = []
    Xb = []
    Ph = []
    Pa = []
    Mxb = []
    CYM = []
    CYO = []
    RR = []
    print('2151 evap_enkf_update complete the data assimilation')


def evap_en4d_bkpe():
    tmpVar = [0 for i in range(zd_evapenk['nen'])]  # 数组初始值为0
    noah_reset(0)
    noah_getlsmstate(zd_evapenk['StVarMN'])
    for t in range(0, zd_evapenk['nen']):
        for i in range(0, zd_evapenk['ntl']):
            for j in range(0, zd_evapenk['nda'] + zd_evapenk['ndv']):
                zd_evapenk['StVarEvap'][t][j][i] = zd_evapenk['StVarMN'][j][i]
    for t in range(0, zd_evapenk['ntl']):
        for j in range(zd_evapenk['ndv'], zd_evapenk['nda'] + zd_evapenk['ndv']):
            ggdlist = [tmpVar, zd_evapenk['nen'], zd_evap['aa'][j], zd_evap['maxstd'][j]]
            ggd_ensemble(ggdlist)
            tmpVar = ggdlist[0]
            for i in range(0, zd_evapenk['nen']):
                if j == 4:
                    tmpVar[i] = tmpVar[i] / 100
                if j == 5:
                    tmpVar[i] = tmpVar[i] / 1000
                tmpVar[i] = zd_evapenk['StVarEvap'][i][j][t] + tmpVar[i]
                qalist = [tmpVar[i], j]
                noah_qa_sm(qalist)
                # tmpVar=qalist[0]
                zd_evapenk['StVarEvap'][i][j][t] = qalist[0]


def evap_en4d_getst(winth, nenth):
    noah_getlsmstate(zd_evapenk['StVarMN'])
    for i in range(0, zd_evapenk['ntl']):
        for j in range(0, zd_evapenk['nda'] + zd_evapenk['ndv']):
            zd_evapenk['Ensem4D'][nenth][j][i][winth] = zd_evapenk['StVarMN'][j][i]  # use StVarMN as tmpStVar


def evap_en4d_clear():
    if lis.a.daalg == 4 or lis.a.daalg == 3:
        zd_evapenk['StVarEvap'] = []
        zd_evapenk['StVarMN'] = []
        zd_evapenk['Ensem4D'] = []


# dataassim_module
def LIS_dataassim_init():
    if (lis.a.daalg > 0):
        evap_enkf_init()
        evapobs_setup()


def LIS_dataassim_clear():
    evap_en4d_clear()


# Lis_log_msg
def lis_log_msg(msg):
    print(time.strftime('2214 %Y-%m-%dT%H:%M:%S', time.localtime()), '', msg.rstrip(), '(', zd_spmd['iam'], ')')


def lis_log_blocked_msg(msg):
    ierr = 0
    lis_log_msg(msg)
    if zd_spmd['npes'] > 1:
        lis_log_msg('DBG: lis_log_blocked_msg -- waiting at barrier')
        lis_log_msg('DBG: lis_log_blocked_msg -- passed barrier')


# lis_openfileMod
zd_opfiMod = {
    'use_opendap_server': True
}


def lis_set_filename(file, time_offset=""):  # file为单独列表

    # if time_offset:
    #   # file[0] = opendap_data_prefix.rstrip()+'/'+ciam.strip()+'/'+"var_"+time_offset+".bin"
    # else:
    #   # file[0] = opendap_data_prefix.rstrip()+'/'+ciam.strip()+'/'+"var.bin"
    return


def lis_open_file(unit, file, form="", status="", access="", recl="", script="", time_offset=""):
    ios = 0
    form_use = ''  # 定义字符串变量
    status_use = ''
    access_use = ''
    chscript_use = ''
    cunit = ''
    use_opendap_server = ''
    if not form:
        form_use = 'unformatted'
    else:
        if form.strip() == 'unformatted' or form.strip() == 'formatted':
            form_use = form.strip()
    if not 'status':
        status_use = 'old'
    else:
        if status.strip() == 'old' or status.strip() == 'new' or status.strip() == 'replace' or status.strip() == 'unknown':
            status_use = status.strip()
    if not 'access':
        if lis.d.gridDesc[8] == 0.01:
            access_use = 'direct'
        else:
            access_use = 'sequential'
    else:
        if access.strip() == 'sequential' or access.strip() == 'direct':
            access_use = access.strip()
    if not 'script':
        script_use = 'none'
    else:
        script_use = script.strip()
    if use_opendap_server:
        if script_use != 'none':
            if not 'time_offset':
                retrieve_data(file, script_use)
            else:
                lis_log_msg('MSG: lis_open_file -- Opening ' + file.rstrip())
    if access_use == 'sequential':  # 两份打开文件有什么不同？
        f = open(file, 'r')
    else:
        f = open(file, 'r')
    cunit = str(format(unit, '>4.0f'))
    if ios != 0:
        lis_log_msg('ERR: lis_open_file -- Cannot open file ' + file.rstrip() + ' on unit ' + cunit.lstrip())
    else:
        lis_log_msg('MSG: lis_open_file -- Successfully opened ' + file.rstrip() + ' on unit ' + cunit.lstrip())


def lis_read_file_new(lrf):  # lrf = [unit,array,gridDesc]
    line1 = 0
    line2 = 0
    line = 0
    c = 0
    r = 0
    glnc = 0
    glnr = 0
    nc = 0
    line1 = round((lis.d.gridDesc[3] - lrf[2][0]) / lrf[2][4]) + 1
    line2 = round((lis.d.gridDesc[4] - lrf[2][1]) / lrf[2][5]) + 1
    nc = round((lrf[2][3] - lrf[2][1]) / lrf[2][5]) + 1
    f = open(lrf[0], 'rb')
    content = f.readlines()
    for r in range(0, int(lis.d.lnr)):
        for c in range(0, lis.d.lnc):
            glnc = line2 + c - 1
            glnr = line1 + r - 1
            line = (glnr - 1) * nc + glnc
            lrf[1][c][r] = content[c]  # 该赋值方式是否正确？
    f.close()


def retrieve_data(file, script, time_offset=""):
    exists = False
    try1 = 1
    if not exists and try1 < 4:  # keep trying to retrieve file
        if not time_offset:
            retrieve_script(script)
        else:
            retrieve_script(script)
        try:
            f = open(file)
            f.close()
        except IOError:
            print("2322 文件不存在")
            exists = True
        try1 = try1 + 1
    else:
        if not exists:  # error, could not retrieve the file
            lis_log_msg('ERR: lis_open_file -- ' + 'Could not retrieve data file ' + file.rstrip())
            # endrun()
        else:
            sys.exit()


def retrieve_script(script):
    lis_log_msg('MSG: lis_open_file -- Retrieving data via ' + script.rstrip() + ' script')


# Resample_module
def vegmask_from_1km():  # vf1 = [nc,nr,mask]
    # import lisdrv_module
    # from lisdrv_module import lis
    # global In,out
    import struct
    import numpy as np
    # the number of the column of the lis.p.VFILE file   pycharm数组下标从0开始，fortran从1开始
    ic = int((lis.d.lc_gridDesc[3] - lis.d.lc_gridDesc[1]) / lis.d.lc_gridDesc[5]) + 1
    # the number of the row of the lis.p.VFILE file
    ir = int((lis.d.lc_gridDesc[2] - lis.d.lc_gridDesc[0]) / lis.d.lc_gridDesc[4]) + 1
    # the number of the vegetation types
    nt = int(lis.p.nt)

    # 模拟区域在mask及vfile文件中的范围，以mask及vfile的原始空间分辨率计算，即读取范围
    cs = int((lis.d.gridDesc[4] - lis.d.lc_gridDesc[1]) / lis.d.lc_gridDesc[5])  # 经度开始位置，列开始位置
    ce = int((lis.d.gridDesc[7] - lis.d.lc_gridDesc[1]) / lis.d.lc_gridDesc[5])  # 经度结束位置，列结束位置
    rs = int((lis.d.gridDesc[3] - lis.d.lc_gridDesc[0]) / lis.d.lc_gridDesc[4])  # 纬度开始位置，行开始位置
    re = int((lis.d.gridDesc[6] - lis.d.lc_gridDesc[0]) / lis.d.lc_gridDesc[4])  # 纬度结束位置，行结束位置

    iresc = int(lis.d.gridDesc[9] / lis.d.lc_gridDesc[5])
    iresr = int(lis.d.gridDesc[8] / lis.d.lc_gridDesc[4])
    if (ce - cs + 1) % iresc == 0:
        nc = int((ce - cs) / iresc * 1.0)  # vf1[0]
    else:
        nc = int((ce - cs) / iresc * 1.0) + 1

    if (re - rs + 1) % iresr == 0:
        nr = int((re - rs) / iresr * 1.0)  # vf1[1]
    else:
        nr = int((re - rs) / iresr * 1.0) + 1

    if nc * iresc < (ce - cs):
        nc = nc + 1
    if nr * iresr < (re - rs):
        nr = nr + 1
    # In = [[0.0] * ic for i in range(0, iresr)]
    In = np.zeros([ce - cs + 1, re - rs + 1], dtype='f4')
    # out = [[[0.0] * nc for j in range(0, 1)] for i in range(0, nt+1)] #对应fortran中的0：nt
    out = np.zeros([nc, nr, nt + 1], dtype='f4')
    mask = [[1 for j in range(0, nr)] for i in range(0, nc)]
    # print('2147',np.array(mask).shape)
    binFile = open(lis.p.vfile, 'rb')  # 文件打开代号为51
    # content = binFile.readlines()

    for j in range(rs, re + 1):
        # for i in range(cs, ce+1):
        grec = (j * ic + cs) * 4
        binFile.seek(grec, 0)
        for i in range(cs, ce + 1):
            context = binFile.read(4)
            # print('2158',context)
            In[i - cs, j - rs], = struct.unpack('>f', context)  # 小端('little- endian')用 ！，大端则不带！
    # print('2161', In)
    for j in range(0, nr):
        minj = j * iresr
        maxj = min(ir, (j + 1) * iresr, re - rs + 1)  # 对j进行修改。

        for i in range(0, nc):
            mini = i * iresc
            maxi = min(ic, (i + 1) * iresc, ce - cs + 1)
            for j1 in range(minj, maxj):
                for i1 in range(mini, maxi):
                    # print(i,i1)
                    # print(In.shape)
                    iveg = int(In[i1, j1])  # veg type
                    out[i, j, iveg] = out[i, j, iveg] + 1.0
            if out[i, j, 0] > 0.5 * (maxi - mini + 1) * maxj:
                mask[i][j] = 0
    # print('2176',mask)           #
    binFile.close()
    del In
    del out
    print('2401 vegmask_from_1km done')
    return mask


def vegcover_from_1km(cf1list):  # cf1=[nc,nr,out]
    import struct
    import numpy as np
    # the number of the column of the lis.p.VFILE file   pycharm数组下标从0开始，fortran从1开始
    ic = int((lis.d.lc_gridDesc[3] - lis.d.lc_gridDesc[1]) / lis.d.lc_gridDesc[5]) + 1
    # the number of the row of the lis.p.VFILE file
    ir = int((lis.d.lc_gridDesc[2] - lis.d.lc_gridDesc[0]) / lis.d.lc_gridDesc[4]) + 1
    # the number of the vegetation types
    nt = int(lis.p.nt)

    # 模拟区域在mask及vfile文件中的范围，以mask及vfile的原始空间分辨率计算，即读取范围
    cs = int((lis.d.gridDesc[4] - lis.d.lc_gridDesc[1]) / lis.d.lc_gridDesc[5])  # 经度开始位置，列开始位置
    ce = int((lis.d.gridDesc[7] - lis.d.lc_gridDesc[1]) / lis.d.lc_gridDesc[5])  # 经度结束位置，列结束位置
    rs = int((lis.d.gridDesc[3] - lis.d.lc_gridDesc[0]) / lis.d.lc_gridDesc[4])  # 纬度开始位置，行开始位置
    re = int((lis.d.gridDesc[6] - lis.d.lc_gridDesc[0]) / lis.d.lc_gridDesc[4])  # 纬度结束位置，行结束位置

    iresc = int((lis.d.gridDesc[9] / lis.d.lc_gridDesc[5]))
    iresr = int((lis.d.gridDesc[8] / lis.d.lc_gridDesc[4]))
    if (ce - cs + 1) % iresc == 0:
        nc = int((ce - cs) / iresc * 1.0)  # vf1[0]
    else:
        nc = int((ce - cs) / iresc * 1.0) + 1

    if (re - rs + 1) % iresr == 0:
        nr = int((re - rs) / iresr * 1.0)  # vf1[1]
    else:
        nr = int((re - rs) / iresr * 1.0) + 1

    if nc * iresc < (ce - cs):
        nc = nc + 1
    if nr * iresr < (re - rs):
        nr = nr + 1

    In = np.zeros([ce - cs + 1, re - rs + 1], dtype='f4')
    out = np.zeros([nc, nr, nt + 1], dtype='f4')
    binFile = open(lis.p.vfile, 'rb')  # 文件打开代号为51

    for j in range(rs, re + 1):
        # for i in range(cs, ce+1):
        grec = (j * ic + cs) * 4
        binFile.seek(grec, 0)
        for i in range(cs, ce + 1):
            context = binFile.read(4)
            # print('2158',context)
            # data = np.fromfile(context, dtype='>f4', count=-1)
            In[i - cs, j - rs], = struct.unpack('>f', context)  # 小端('little- endian')用 ！，大端则不带！

    for j in range(0, nr):
        minj = j * iresr
        maxj = min(ir, (j + 1) * iresr, re - rs + 1)  # 对j进行修改。
        for i in range(0, nc):
            mini = i * iresc
            maxi = min(ic, (i + 1) * iresc, ce - cs + 1)
            for j1 in range(minj, maxj):
                for i1 in range(mini, maxi):
                    iveg = int(In[i1, j1])  # veg type
                    out[i, j, iveg] = out[i, j, iveg] + 1.0  # out列表的第三维不需要-1.
    binFile.close()
    cf1list[0] = nc
    cf1list[1] = nr
    cf1list[2] = out
    print("2466 vegcover_from_1km done")
    del In
    return out


def fillup(fil):  # fil = [out, mask, outout, paraDesc,maskDesc, nc,nr, undef, rad, fill]
    # import numpy as np
    # import  lisdrv_module
    # from lisdrv_module import lis
    for j in range(0, int(lis.d.lnr)):
        for i in range(0, int(lis.d.lnc)):
            tmp = 0.0
            # print('2244', np.array(fil[1]).shape)
            if fil[1][i][j] > 0.5 and fil[0][i][j] == fil[7]:
                ifound = 0
                for irad in range(0, fil[8]):
                    for j1 in range(j - 1 - irad, j - 1 + irad + 1):
                        jx = j1
                        if jx <= 0:
                            jx = 0
                        if jx >= fil[6]:
                            jx = fil[6] - 1
                        for i1 in range(i - 1 - irad, i - 1 + irad + 1):
                            ix = i1
                            if ix <= 0:
                                ix = 0
                            if ix >= fil[5]:
                                ix = fil[5] - 1
                            if fil[0][ix][jx] != fil[7]:
                                ifound = ifound + 1
                                tmp = tmp + fil[0][ix][jx]
                    if ifound >= 1:
                        exit()
                if ifound >= 1:
                    fil[2][i][j] = tmp / ifound
                else:
                    fil[2][i][j] = fil[9]
            if fil[1][i][j] > 0.5 and fil[0][i][j] != fil[7]:
                fil[2][i][j] = fil[0][i][j]
    print('2506 fillup done')


def avgdata_from_1km(af1):  # af1=[inputf,outout,soil_gridDesc,mask_gridDesc]
    import struct
    import numpy as np
    # global out
    undef = -9999.0
    # out = [[0] * 50 for i in range(50)]
    # In = [[0] * 50 for i in range(50)]
    # mask = [[0] * 50 for i in range(50)]
    ios1 = 0
    # the number of the column of the lis.p.VFILE file   pycharm数组下标从0开始，fortran从1开始
    ic = int((af1[2][3] - af1[2][1]) / af1[2][5]) + 1
    # the number of the row of the lis.p.VFILE file
    ir = int((af1[2][2] - af1[2][0]) / af1[2][4]) + 1

    # 模拟区域在mask及vfile文件中的范围，以mask及vfile的原始空间分辨率计算，即读取范围
    cs = int((lis.d.gridDesc[4] - lis.d.lc_gridDesc[1]) / lis.d.lc_gridDesc[5])  # 经度开始位置，列开始位置
    ce = int((lis.d.gridDesc[7] - lis.d.lc_gridDesc[1]) / lis.d.lc_gridDesc[5])  # 经度结束位置，列结束位置
    rs = int((lis.d.gridDesc[3] - lis.d.lc_gridDesc[0]) / lis.d.lc_gridDesc[4])  # 纬度开始位置，行开始位置
    re = int((lis.d.gridDesc[6] - lis.d.lc_gridDesc[0]) / lis.d.lc_gridDesc[4])  # 纬度结束位置，行结束位置

    # the zoom scale of the parameter file from original resolution to current simulation resolution
    # 参数文件从原始分辨率到当前模拟分辨率的缩放比例
    iresc = int(lis.d.gridDesc[9] / af1[2][5])
    iresr = int(lis.d.gridDesc[8] / af1[2][4])
    # the number of the column of the simulating required parameter file 模拟所需参数文件的列号
    nc = int(ic / iresc * 1.0)
    # the number of the row of the simulating required parameter file 模拟所需参数文件的行数
    nr = int(ir / iresr * 1.0)
    if (ce - cs + 1) % iresc == 0:
        nc = int((ce - cs) / iresc * 1.0)  # vf1[0]
    else:
        nc = int((ce - cs) / iresc * 1.0) + 1

    if (re - rs + 1) % iresr == 0:
        nr = int((re - rs) / iresr * 1.0)  # vf1[1]
    else:
        nr = int((re - rs) / iresr * 1.0) + 1

    if nc * iresc < (ce - cs):
        nc = nc + 1
    if nr * iresr < (re - rs):
        nr = nr + 1

    In = np.zeros([ce - cs + 1, re - rs + 1], dtype='f4')
    out = np.zeros([nc, nr], dtype='f4')

    binFile = open(str(af1[0]), 'rb')

    cnt = 0

    for j in range(rs, re + 1):
        # for i in range(cs, ce+1):
        grec = (j * ic + cs) * 4
        binFile.seek(grec, 0)
        for i in range(cs, ce + 1):
            context = binFile.read(4)
            # print('2334', context)
            In[i - cs, j - rs], = struct.unpack('>f', context)  # 小端('little- endian')用！，大端则不带！

    for j in range(0, nr):
        minj = j * iresr
        maxj = min(ir, (j + 1) * iresr, re - rs + 1)
        for i in range(0, nc):
            mini = i * iresc
            maxi = min(ic, (i + 1) * iresc, ce - cs + 1)
            cnt = 0
            for j1 in range(minj, maxj):
                for i1 in range(mini, maxi):
                    if In[i1, j1] != undef:  # In为np数组，np.zeos生成的数组下标为[,],列表下标为[][]
                        out[i, j] = out[i, j] + In[i1, j1]
                        cnt = cnt + 1
            if cnt > 0:
                out[i, j] = out[i, j] / cnt
            else:
                out[i, j] = undef
    # vf1list = [nc,nr,mask]
    mask = vegmask_from_1km()
    # mask1=[[mask[j][i] for i in range(nr)] for j in range(nc)] #行列互换
    outout = out
    filist = [out, mask, outout, af1[1], af1[2], af1[3], nc, nr, undef, 2, 0.0]
    fillup(filist)
    af1[1] = filist[2]
    del In
    del out
    del mask
    # del mask1
    print('2595', 'avgdata_from_1km done')


# noah_settbot
def noah_settbot():
    placetbot = [[0.0] * int(zd_indices['lis_nc_data']) for i in range(int(zd_indices['lis_nr_data']))]
    # out = [[0.0] * 50 for i in range(50)] #注意在Fortran中out是否有特殊意思。
    out = [[]]
    lsflist = [noahdrv.noah_tbot]
    lis_set_filename(lsflist)
    noahdrv.noah_tbot = lsflist[0]
    # out = [[0.0] * 50 for i in range(50)]
    print('2607 Opening GLDAS TBOT File', noahdrv.noah_tbot)
    if lis.d.gridDesc[8] != 0.01:
        aflist = [noahdrv.noah_tbot, out, noahdrv.noah_gridDesc, lis.d.lc_gridDesc]
        avgdata_from_1km(aflist)
        noahdrv.noah_tbot = aflist[0]
        out = aflist[1]
        noahdrv.noah_gridDesc = aflist[2]
        lis.d.lc_gridDesc = aflist[3]
        placetbot = out
        out = []  # 二维列表释放。
    else:
        lis_open_file(12, noahdrv.noah_tbot, "unformatted", 'old', 'direct', 1, 'gettbot.pl')
        lrflist = [noahdrv.noah_tbot, placetbot, noahdrv.noah_gridDesc]
        lis_read_file_new(lrflist)  # 这里将代号12参数换成对应文件
        placetbot = lrflist[1]
        noahdrv.noah_gridDesc = lrflist[2]
    print('2623 MSG: setnoahp -- Read  TBOT file', ' (', zd_spmd['iam'], ')')
    for i in range(0, lis.d.nch):
        if placetbot[zd_lisdrv['tile'][i].col][zd_lisdrv['tile'][i].row - zd_indices['lis_tnroffset']] != -9999.00:
            zd_varder['noah'][i].tempbot = placetbot[zd_lisdrv['tile'][i].row - zd_indices['lis_tnroffset'] - 1][
                zd_lisdrv['tile'][i].col - 1]
    placetbot = []


# Noah_setmxalb
def noah_setmxalb():
    tmpalb = [[0] * int(zd_indices['lis_nc_data']) for i in range(int(zd_indices['lis_nr_data']))]
    lsflist = [noahdrv.noah_mxsnal]
    lis_set_filename(lsflist)
    noahdrv.noah_mxsnal = lsflist[0]
    out = []  # [[0] * 100 for i in range(100)]
    print('2638 MSG: setnoahp -- Opening GLDAS MAXSNALB File', noahdrv.noah_mxsnal.rstrip(), zd_spmd['iam'])
    if lis.d.gridDesc[8] != 0.01:  # 该行有错误？
        aflist = [noahdrv.noah_mxsnal, out, noahdrv.noah_gridDesc, lis.d.lc_gridDesc]
        avgdata_from_1km(aflist)
        noahdrv.noah_mxsnal = aflist[0]
        out = aflist[1]
        noahdrv.noah_gridDesc = aflist[2]
        lis.d.lc_gridDesc = aflist[3]
        for i in range(0, int(lis.d.lnc)):
            for j in range(0, int(lis.d.lnr)):
                tmpalb[j][i] = out[i][j]
        # 源程序中该赋值语句是否正确？
        out = []
    else:
        lis_open_file(12, file=noahdrv.noah_mxsnal, access='direct', form='unformatted', recl=1,
                      script='getmaxsnalb.pl')
        lrflist = [noahdrv.noah_mxsnal, tmpalb, noahdrv.noah_gridDesc]
        lis_read_file_new(lrflist)
        tmpalb = lrflist[1]
        noahdrv.noah_gridDesc = lrflist[2]
    for i in range(0, lis.d.nch):
        if tmpalb[zd_lisdrv['tile'][i].row - zd_indices['lis_tnroffset'] - 1][zd_lisdrv['tile'][i].col - 1] != -9999.00:
            zd_varder['noah'][i].mxsnalb = tmpalb[zd_lisdrv['tile'][i].row - zd_indices['lis_tnroffset'] - 1][
                zd_lisdrv['tile'][i].col - 1]
    tmpalb = []


# Read_faosand
def read_faosand(array):
    out = [[]]  # 注意out是否是Fortran内置数组?.L
    print('2668 MSG: Reading FAO sand file')
    lsflist = [lis.p.safile]  # 配置文件中 lis.p.safile=sand60_1KM.1gd4r
    lis_set_filename(lsflist)
    lis.p.safile = lsflist[0]
    if lis.d.gridDesc[8] != 0.01:
        aflist = [lis.p.safile, out, lis.d.soil_gridDesc, lis.d.lc_gridDesc]
        avgdata_from_1km(aflist)
        lis.p.safile = aflist[0]
        out = aflist[1]
        lis.d.soil_gridDesc = aflist[2]
        lis.d.lc_gridDesc = aflist[3]
        for i in range(0, int(lis.d.lnc)):
            for j in range(0, int(lis.d.lnr)):
                array[j][i] = out[i][j]
        # if out != []:
        out = []
    else:
        lis_open_file(15, lis.p.safile, 'unformatted', 'old', 'direct', 1, 'getsand.pl')  # 读取文件代号是怎么回事？。L
        line1 = int((lis.d.gridDesc[3] - lis.d.soil_gridDesc[0]) / lis.d.gridDesc[8]) + 1
        line2 = int((lis.d.gridDesc[4] - lis.d.soil_gridDesc[1]) / lis.d.gridDesc[9]) + 1
        nc_dom = int((lis.d.soil_gridDesc[3] - lis.d.soil_gridDesc[1]) / lis.d.soil_gridDesc[5]) + 1
        f = open(lis.p.safile, 'rb')
        content = f.readlines()
        for r in range(0, int(lis.d.lnr)):
            for c in range(0, lis.d.lnc):
                glnc = line2 + c
                glnr = line1 + r
                line = (glnr - 1) * nc_dom + glnc - 1
                array[c][r] = content[line]


# Read_faoclay
def read_faoclay(array):
    out = [[0.0] * 100 for i in range(100)]
    lsflist = [lis.p.clfile]  # 配置文件中 lis.p.clfile = clay60_10KM.1gd4r
    lis_set_filename(lsflist)
    lis.p.clfile = lsflist[0]
    if lis.d.gridDesc[8] != 0.01:
        afllist = [lis.p.clfile, out, lis.d.soil_gridDesc, lis.d.lc_gridDesc]
        avgdata_from_1km(afllist)
        lis.p.clfile = afllist[0]
        out = afllist[1]
        lis.d.soil_gridDesc = afllist[2]
        lis.d.lc_gridDesc = afllist[3]
        for i in range(0, int(lis.d.lnc)):
            for j in range(0, int(lis.d.lnr)):
                array[j][i] = out[i][j]
        # array[1:lis_nc_data,1:lis.d.lnr]=out[1:lis.d.lnc,1:lis.d.lnr]
        out = []
    else:
        lis_open_file(15, lis.p.safile, 'unformatted', 'old', 'direct', 1, 'getsand.pl')  # 读取文件代号是怎么回事？。L
        line1 = int((lis.d.gridDesc[3] - lis.d.soil_gridDesc[0]) / lis.d.gridDesc[8]) + 1
        line2 = int((lis.d.gridDesc[4] - lis.d.soil_gridDesc[1]) / lis.d.gridDesc[9]) + 1
        nc_dom = int((lis.d.soil_gridDesc[3] - lis.d.soil_gridDesc[1]) / lis.d.soil_gridDesc[5]) + 1
        f = open(lis.p.safile, 'rb')
        content = f.readlines()
        for r in range(0, int(lis.d.lnr)):
            for c in range(0, lis.d.lnc):
                glnc = line2 + c - 1
                glnr = line1 + r - 1
                line = (glnr - 1) * nc_dom + glnc
                array[c][r] = content[line]
                # read(15,rec=line) array(c,r) 读取文件，将文件内容复制给数组，后期注意改写。L
        f.close()


# Soiltype
def soiltype(nc, nr, sand, clay, silt, soiltyp):  # 参数silt没有起作用。
    # class = 0 # 该变量定义涉及Python常量，是否有特殊意义？。L
    cl, sa = 0.0, 0.0
    print('2738 Using Zobler Texture classification..')
    for j in range(0, nr):
        for i in range(0, nc):
            if clay[j][i] < 0.00:
                soiltyp[j][i] = -99
            else:
                cl = clay[j][i]
                sa = sand[j][i]
            if cl < 0.23:
                if sa < 0.50:
                    soiltyp[j][i] = 8
                else:
                    if sa < 0.75:
                        soiltyp[j][i] = 4
                    else:
                        soiltyp[j][i] = 1
            else:
                if cl < 0.28:
                    if sa < 0.45:
                        soiltyp[j][i] = 8
                    else:
                        soiltyp[j][i] = 7
                else:
                    if cl < 0.37:
                        if sa < 0.2:
                            soiltyp[j][i] = 2  # ! silty clay loam
                        else:
                            if sa < 0.43:
                                soiltyp[j][i] = 6  # ! clay loam
                            else:
                                soiltyp[j][i] = 7  # ! sandy clay
                    else:
                        if cl < 0.41:
                            if sa < 0.2:
                                soiltyp[j][i] = 2  # ! light clay
                            else:
                                if sa < 0.43:
                                    soiltyp[j][i] = 6
                                else:
                                    soiltyp[j][i] = 5
                        else:
                            if sa < 0.43:
                                soiltyp[j][i] = 3
                            else:
                                soiltyp[j][i] = 5


# Noah_setsoils
def noah_setsoils():  # 注意定义全局变量
    print('2787 MSG: setnoahp -- reading soil and clay files', '(', zd_spmd['iam'], ')')
    lis_data_nc = int(zd_indices['lis_nc_data'])
    lis_data_nr = int(zd_indices['lis_nr_data'])

    sand1 = [[0.0] * lis_data_nc for i in range(lis_data_nr)]
    clay1 = [[0.0] * lis_data_nc for i in range(lis_data_nr)]
    silt1 = [[0.0] * lis_data_nc for i in range(lis_data_nr)]
    read_faosand(sand1)
    read_faoclay(clay1)
    print('2796 MSG: setnoahp -- read sand and clay files', ' (', zd_spmd['iam'], ')')
    soiltyp = [[0.0] * lis_data_nc for i in range(lis_data_nr)]
    soiltype(lis_data_nc, lis_data_nr, sand1, clay1, silt1, soiltyp)
    f = open(noahdrv.noah_sfile, 'r')
    content = f.readlines()
    basicset = []
    for i in range(0, noahdrv.noah_nsoilp):
        column = content[i].replace(" ", ",").split(',')[0:noahdrv.noah_zst]
        basicset.append(column)
    f.close()
    print('2806 msg: setnoahp -- read sfile ', noahdrv.noah_sfile.replace("", ""), ' (', zd_spmd['iam'], ')')
    placesltyp = [0 for x in range(0, lis.d.nch)]
    for i in range(0, lis.d.nch):
        placesltyp[i] = soiltyp[zd_lisdrv['tile'][i].row - zd_indices['lis_tnroffset'] - 1][
            zd_lisdrv['tile'][i].col - 1]
        zd_varder['noah'][i].zobsoil = placesltyp[i]
    soiltyp = []
    # absoft_release_cache()  该函数源程序为注释状态。L
    for i in range(0, lis.d.nch):
        k = placesltyp[i]
        for j in range(0, noahdrv.noah_nsoilp):
            zd_varder['noah'][i].SOILP[j] = basicset[j][k - 1]
    placesltyp = []
    sand1 = []
    clay1 = []
    # absoft_release_cache()  该函数源程序为注释状态。L


# Noah_gfrac
def noah_gfrac(f1list=None):
    value1 = [[0.0] * int(zd_indices['lis_nc_data']) for i in
              range(int(zd_indices['lis_nr_data']))]  # temporary value holder for mo1 mo1的临时值
    value2 = [[0.0] * int(zd_indices['lis_nc_data']) for i in
              range(int(zd_indices['lis_nr_data']))]  # temporary value holder for mo2 mo2的临时值
    out = []
    time1, doy1, gmt1 = 0.0, 0, 0.0
    time2, doy2, gmt2 = 0.0, 0, 0.0
    noahdrv.noah_gflag = 0
    zeroi = 0
    numi = 16
    if lis.t.da < 16:
        mo1 = lis.t.mo - 1
        yr1 = lis.t.yr
        if mo1 == 0:
            mo1 = 12
            yr1 = lis.t.yr - 1
        mo2 = lis.t.mo
        yr2 = lis.t.yr
    else:
        mo1 = lis.t.mo
        yr1 = lis.t.yr
        mo2 = lis.t.mo + 1
        yr2 = lis.t.yr
        if mo2 == 13:
            mo2 = 1
            yr2 = lis.t.yr + 1
    d2tlist = [time1, doy1, gmt1, yr1, mo1, numi, zeroi, zeroi, zeroi]
    date2time(d2tlist)
    time1 = d2tlist[0]
    doy1 = d2tlist[1]
    gmt = d2tlist[2]
    d2tlist2 = [time2, doy2, gmt2, yr2, mo2, numi, zeroi, zeroi, zeroi]
    date2time(d2tlist2)
    time2 = d2tlist2[0]
    doy2 = d2tlist2[1]
    gmt2 = d2tlist2[2]
    # ------------------------------------------------------------------------
    #  Weights to be used to interpolate greenness fraction values. 用于插值绿色分数值的权重
    # ------------------------------------------------------------------------
    wt1 = (time2 - lis.t.time) / (time2 - time1)
    wt2 = (lis.t.time - time1) / (time2 - time1)
    # ------------------------------------------------------------------------
    #  Determine if GFRAC files need to be updated 判断是否需要更新GFRAC文件
    # ------------------------------------------------------------------------
    if time2 > noahdrv.noah_gfractime:
        gfrac_flag = 1
    else:
        gfrac_flag = 0
    if gfrac_flag == 1:
        noahdrv.noah_gfractime = time2
        noahdrv.noah_gflag = 1
        # ------------------------------------------------------------------------
        # Open greenness fraction dataset of months corresponding to
        # time1 and time2 for selected LDAS domain and read data.
        # 打开与所选LDAS域的时间1和时间2对应的月份的绿色分数数据集并读取数据。
        # ------------------------------------------------------------------------
        # mm1 = str(format(mo1,'>2.2f'))
        # mm2 = str(format(mo2, '>2.2f'))
        if mo1 > 10:
            mm1 = str(mo1)
        else:
            mm1 = '0' + str(mo1)

        if mo2 > 10:
            mm2 = str(mo2)
        else:
            mm2 = '0' + str(mo2)

        if lis.d.gridDesc[8] != 0.01:
            file1 = noahdrv.noah_mgfile.rstrip() + 'gfrac_' + mm1 + '.1gd4r'
            lsflist = [file1]
            lis_set_filename(lsflist, time_offset=mm1)
            file1 = lsflist[0]
            af1list = [file1, out, noahdrv.noah_gridDesc, lis.d.lc_gridDesc]
            avgdata_from_1km(af1list)
            file1 = af1list[0]
            out = af1list[1]
            noahdrv.noah_gridDesc = af1list[2]
            lis.d.lc_gridDesc = af1list[3]
            for i in range(0, int(lis.d.lnc)):
                for j in range(0, int(lis.d.lnr)):
                    value1[j][i] = out[i][j]
            # value1[1:lis_nc_data,1:lis.d.lnr]=out[1:lis.d.lnc,1:lis.d.lnr]
            out = []
            file2 = noahdrv.noah_mgfile.rstrip() + 'gfrac_' + mm2 + '.1gd4r'  # 1gd4r'
            lsflist2 = [file2]
            lis_set_filename(lsflist2, time_offset=mm2)
            file2 = lsflist2[0]
            af1list2 = [file2, out, noahdrv.noah_gridDesc, lis.d.lc_gridDesc]
            avgdata_from_1km(af1list2)
            file2 = af1list2[0]
            out = af1list2[1]
            noahdrv.noah_gridDesc = af1list2[2]
            lis.d.lc_gridDesc = af1list2[3]
            value2 = out
            out = []
        else:
            file1 = noahdrv.noah_mgfile.rstrip() + 'gfrac_' + mm1 + '.1gd4r'
            lsflist = [file1]
            lis_set_filename(lsflist, time_offset=mm1)
            file1 = lsflist[0]
            lis_open_file(10, file1, status='old', form='unformatted', access='direct', recl=4, script='getgfrac.pl',
                          time_offset=mm1)
            lrfn = [file1, value1, noahdrv.noah_gridDesc]
            lis_read_file_new(lrfn)
            value1 = lrfn[1]
            file2 = noahdrv.noah_mgfile.rstrip() + 'gfrac_' + mm2 + '.1gd4r'
            lsflist2 = [file2]
            lis_set_filename(lsflist2, time_offset=mm2)
            file2 = lsflist2[0]
            lis_open_file(11, file2, status='old', form='unformatted', access='direct', recl=4, script='getgfrac.pl',
                          time_offset=mm2)
            lrfn2 = [file2, value2, noahdrv.noah_gridDesc]
            lis_read_file_new(lrfn2)
            value2 = lrfn2[1]
        # ------------------------------------------------------------------------
        # Assign MONTHLY vegetation greenness fractions to each tile. 为每个方块指定每月植被绿色度分数
        # ------------------------------------------------------------------------
        for i in range(0, lis.d.nch):
            if (value1[zd_lisdrv['tile'][i].row - zd_indices['lis_tnroffset']][
                    zd_lisdrv['tile'][i].col] != -9999.000) and (
                    value2[zd_lisdrv['tile'][i].row - zd_indices['lis_tnroffset']][
                        zd_lisdrv['tile'][i].col] != -9999.000):
                zd_varder['noah'][i].vegmp1 = value1[zd_lisdrv['tile'][i].row - zd_indices['lis_tnroffset']][
                    zd_lisdrv['tile'][i].col]
                zd_varder['noah'][i].vegmp2 = value2[zd_lisdrv['tile'][i].row - zd_indices['lis_tnroffset']][
                    zd_lisdrv['tile'][i].col]
    # ------------------------------------------------------------------------
    #  Interpolate greenness fraction values once daily 每天插值一次绿色分数值
    # ------------------------------------------------------------------------
    print('2956 noah_gfrac')
    if noahdrv.noah_gfracdchk != lis.t.da:
        noahdrv.noah_gflag = 1
        if lis.a.daalg == 0 and lis.a.pSLGA:
            for i in range(0, lis.d.nch):
                ran1 = random.random()
                zd_varder['noah'][i].VEGIP = (
                        (wt1 * zd_varder['noah'][i].vegmp1) + (wt2 * zd_varder['noah'][i].vegmp2))  # +ran1/2-0.25
                if zd_varder['noah'][i].VEGIP < 0:
                    zd_varder['noah'][i].VEGIP = 0 - zd_varder['noah'][i].VEGIP
                if zd_varder['noah'][i].VEGIP > 1:
                    zd_varder['noah'][i].VEGIP = 1.0
        else:
            for i in range(0, lis.d.nch):
                zd_varder['noah'][i].VEGIP = (
                        (wt1 * zd_varder['noah'][i].vegmp1) + (wt2 * zd_varder['noah'][i].vegmp2))
        noahdrv.noah_gfracdchk = lis.t.da
        print('2973 Done noah_gfrac', ' (', zd_spmd['iam'], ')')
        if lis.o.wparam == 1:  # 0-do not,1-write model parameters in output
            # allocate(gfracout(lisimpo.dimpo.lnc,lisimpo.dimpo.lnr))
            gfracout = -9999.0
            if lis.dimpo.gridDesc[8] != 0.01:
                for i in range(0, lis.dimpo.nch):
                    if (lis.d.gridDesc[6] >= zd_lisdrv['grid'][i].lat >= lis.d.gridDesc[4] <=
                            zd_lisdrv['grid'][i].lon <= lis.d.gridDesc[7]):
                        rindex = zd_lisdrv['tile'][i].row - round((lis.dimpo.gridDesc[3] - lis.dimpo.gridDesc[43]) \
                                                                  / lis.d.gridDesc[8])
                        cindex = zd_lisdrv['tile'][i].col - round((lis.dimpo.gridDesc[4] - lis.dimpo.gridDesc[44]) \
                                                                  / lis.d.gridDesc[9])
                        gfracout[rindex][cindex] = zd_varder['noah'][i].VEGIP * 1.0

            else:
                if lis.d.gridDesc[8] == 0.01:
                    for i in range(0, lis.d.nch):
                        gfracout[zd_lisdrv['tile'][i].row][zd_lisdrv['tile'][i].col] = zd_varder['noah'][i].VEGIP * 1.0
            with open("gfracout.bin", 'w') as f:
                f.write(gfracout)
            f.close()
            gfracout = []


# noah_alb
def noah_alb():
    # === Local Variables ======局部变量============
    albout = [[0.0] * 50 for i in range(50)]
    out = [[0.0] * 50 for i in range(50)]
    cindex, rindex = 0, 0
    line, line1, line2, glnc, glnr, nc, nr = 0, 0, 0, 0, 0, 0, 0
    ios1 = 0
    i, j, c, r = 0, 0, 0, 0  # loop counters  循环计数器
    janda, janmo = 0, 0  # january 31  1月31日
    aprda, aprmo = 0, 0  # april 30    4月30日
    julda, julmo = 0, 0  # july 31     7月31日
    octda, octmo = 0, 0  # october 31  10月31日
    yr = 0  # year of run   运行年份
    doy1 = 0  # temporary time variables  临时时间变量
    zeroi = 0  # integer number holders  整数持有
    albflag = 0  # flag to update ALBEDO values  更新反照率的标志
    time = 0.0  # current model time variable 当前模型时间变量
    jan31, apr30, jul31, oct31 = 0.0, 0.0, 0.0, 0.0  # dates of quarterly ALBEDO files  集度反照率文件日期
    qdif = 0.0  # difference between q1 and q2 times  q1和q2之间的不同
    timdif = 0.0  # difference between time and q1 time time 和 q1time之间的不同
    gmt1, gmt2 = 0.0, 0.0  # gmt values
    value1 = [[0.0] * int(zd_indices['lis_nc_data']) for i in
              range(int(zd_indices['lis_nr_data']))]  # Temporary value holder for QQ1  qq1的临时值
    value2 = [[0.0] * int(zd_indices['lis_nc_data']) for i in
              range(int(zd_indices['lis_nr_data']))]  # Temporary value holder for QQ2  qq1的临时值
    ran1 = 0

    valdif = [0.0 for i in range(lis.d.nch)]  # Difference of QQ2 and QQ1 ALBEDO    QQ2和QQ1的反照率不同
    qq1, qq2 = '', ''  # Filename places for quarter values  四分之一值的文件名位置
    file1, file2 = '', ''
    # === End Variable Definition ==========结束变量定义================

    zeroi = 0
    noahdrv.noah_aflag = 0
    # -------------------------------------------------------------------------
    # Determine Dates of the quarters in terms of Year (e.g., 1999.3)  根据年份确定日期
    # -------------------------------------------------------------------------
    time = lis.t.time
    yr = lis.t.yr
    # -------------------------------------------------------------------------
    #  January 31         1月31日
    # -------------------------------------------------------------------------
    janda = 31
    janmo = 1
    d2tlist = [jan31, doy1, gmt1, yr, janmo, janda, zeroi, zeroi, zeroi]
    date2time(d2tlist)
    jan31 = d2tlist[0]
    doy1 = d2tlist[1]
    gmt1 = d2tlist[2]
    # -------------------------------------------------------------------------
    #  April 30    4月30日
    # -------------------------------------------------------------------------
    aprda = 30
    aprmo = 4
    d2tlist2 = [apr30, doy1, gmt1, yr, aprmo, aprda, zeroi, zeroi, zeroi]
    date2time(d2tlist2)
    apr30 = d2tlist2[0]
    doy1 = d2tlist2[1]
    gmt1 = d2tlist2[2]
    # -------------------------------------------------------------------------
    #  July 31    7月31日
    # -------------------------------------------------------------------------
    julda = 31
    julmo = 7
    d2tlist3 = [jul31, doy1, gmt1, yr, julmo, julda, zeroi, zeroi, zeroi]
    date2time(d2tlist3)
    jul31 = d2tlist3[0]
    doy1 = d2tlist3[1]
    gmt1 = d2tlist3[2]
    # -------------------------------------------------------------------------
    #  October 31   10月31日
    # -------------------------------------------------------------------------
    octda = 31
    octmo = 10
    d2tlist4 = [oct31, doy1, gmt1, yr, octmo, octda, zeroi, zeroi, zeroi]
    date2time(d2tlist4)
    oct31 = d2tlist4[0]
    doy1 = d2tlist4[1]
    gmt1 = d2tlist4[2]
    # -------------------------------------------------------------------------
    # Determine which two quarterly ALBEDO files book-end model time.  确定哪两个季度反照率文件记录结束模型时间
    # -------------------------------------------------------------------------
    if jan31 <= time <= apr30:
        qq1 = "01"
        qq2 = "02"
        qdif = apr30 - jan31
        timdif = time - jan31
        albflag = 1
    if apr30 <= time <= jul31:
        qq1 = "02"
        qq2 = "03"
        qdif = jul31 - apr30
        timdif = time - apr30
        albflag = 2
    if jul31 <= time <= oct31:
        qq1 = "03"
        qq2 = "04"
        qdif = oct31 - jul31
        timdif = time - jul31
        albflag = 3
    if time >= oct31:
        qq1 = "04"
        qq2 = "01"
        qdif = (jan31 + 1.0) - oct31
        timdif = time - oct31
        albflag = 4
    if time < jan31:
        qq1 = "04"
        qq2 = "01"
        oct31 = oct31 - 1.0
        qdif = jan31 - oct31
        timdif = time - oct31
        albflag = 5
    if noahdrv.noah_albtime != albflag:
        noahdrv.noah_albtime = albflag
        noahdrv.noah_aflag = 1
        # -------------------------------------------------------------------------
        #  Open the needed two quarterly snow-free ALBEDO files   打开所需的两个季度无雪反射率文件
        # -------------------------------------------------------------------------
        if lis.d.gridDesc[8] != 0.01:
            file1 = noahdrv.noah_albfile.rstrip() + 'alb_' + qq1 + '.1gd4r'
            lsflist = [file1]
            lis_set_filename(lsflist, time_offset=qq1)
            file1 = lsflist[0]
            af1 = [file1, out, noahdrv.noah_gridDesc, lis.d.lc_gridDesc]
            avgdata_from_1km(af1)
            file1 = af1[0]
            out = af1[1]
            noahdrv.noah_gridDesc = af1[2]
            lis.d.lc_gridDesc = af1[3]
            for i in range(0, int(lis.d.lnc)):
                for j in range(0, int(lis.d.lnr)):
                    value2[j][i] = out[i][j]
            file2 = noahdrv.noah_albfile.rstrip() + 'alb_' + qq2 + '.1gd4r'
            lsflist2 = [file1]
            lis_set_filename(lsflist2, time_offset=qq1)
            file1 = lsflist2[0]
            af2 = [file2, out, noahdrv.noah_gridDesc, lis.d.lc_gridDesc]
            avgdata_from_1km(af2)
            file2 = af2[0]
            out = af2[1]
            noahdrv.noah_gridDesc = af2[2]
            lis.d.lc_gridDesc = af2[3]
            out = []
        else:
            file1 = noahdrv.noah_albfile.rstrip() + 'alb_' + qq2 + '.1gd4r'
            lsflist = [file1]
            lis_set_filename(lsflist, time_offset=qq1)
            file1 = lsflist[0]
            lis_open_file(10, file=file1, status='old', form='unformatted',
                          access='direct', recl=4, script='getalbedo.pl', time_offset=qq1)
            lrfnlist = [file1, value1, noahdrv.noah_gridDesc]
            lis_read_file_new(lrfnlist)
            value1 = lrfnlist[1]
            file2 = noahdrv.noah_albfile.rstrip() + 'alb_' + qq2 + '.1gd4r'
            lsflist2 = [file2]
            lis_set_filename(lsflist2, time_offset=qq2)
            file2 = lsflist2[0]
            lis_open_file(11, file=file2, status='old', form='unformatted',
                          access='direct', recl=4, script='getalbedo.pl', time_offset=qq2)
            lrfnlist2 = [file2, value2, noahdrv.noah_gridDesc]
            lis_read_file_new(lrfnlist2)
            value2 = lrfnlist2[1]
        # -------------------------------------------------------------------------
        # Assign quarterly ALBEDO fractions to each tile.  为每个区块指定季度反照率分数
        # -------------------------------------------------------------------------
        for i in range(0, lis.d.nch):
            if ((value1[zd_lisdrv['tile'][i].row - zd_indices['lis_tnroffset']][
                     zd_lisdrv['tile'][i].col] != -9999.000) and (
                    value2[zd_lisdrv['tile'][i].row - zd_indices['lis_tnroffset']][
                        zd_lisdrv['tile'][i].col] != -9999.000)):
                zd_varder['noah'][i].albsf1 = value1[zd_lisdrv['tile'][i].row - zd_indices['lis_tnroffset']][
                    zd_lisdrv['tile'][i].col]
                zd_varder['noah'][i].albsf2 = value2[zd_lisdrv['tile'][i].row - zd_indices['lis_tnroffset']][
                    zd_lisdrv['tile'][i].col]
    # -------------------------------------------------------------------------
    # Assign ALBEDO fractions to each zd_lisdrv['tile'] and interpolate daily.将反照率分数分配给zd_lisdrv['tile']并且每天插值
    # -------------------------------------------------------------------------
    if noahdrv.noah_albdchk != lis.t.da:
        noahdrv.noah_aflag = 1
        if lis.a.daalg == 0 and lis.a.pSLGA:
            for i in range(0, lis.d.nch):
                if zd_varder['noah'][i].albsf1 != -9999.000:
                    random.randint(0, 1)
                    valdif[i] = zd_varder['noah'][i].albsf2 - zd_varder['noah'][i].albsf1
                    zd_varder['noah'][i].albsf = ((timdif * valdif[i] / qdif) + zd_varder['noah'][
                        i].albsf1) * 5  # +ran1/2-0.25
                    if zd_varder['noah'][i].albsf < 0:
                        zd_varder['noah'][i].albsf = 0 - zd_varder['noah'][i].albsf
                    if zd_varder['noah'][i].albsf > 1:
                        zd_varder['noah'][i].albsf = 1.0
        else:
            for i in range(0, lis.d.nch):
                if zd_varder['noah'][i].albsf1 != -9999.000:
                    if 1168 <= i <= 1169:
                        print('274', i, line2 - 1 + zd_lisdrv['tile'][i].col,
                              line1 - 1 + zd_lisdrv['tile'][i].row - zd_indices['lis_tnroffset'],
                              zd_varder['noah'][i].albsf1, timdif, qdif, zd_varder['noah'][i].albsf2)
                    valdif[i] = zd_varder['noah'][i].albsf2 - zd_varder['noah'][i].albsf1
                    zd_varder['noah'][i].albsf = (timdif * valdif[i] / qdif) + zd_varder['noah'][i].albsf1
        noahdrv.noah_albdchk = lis.t.da
        if lis.o.wparam == 1:
            albout = [[0.0] * lis.d.lnc for i in range(lis.d.lnr)]
            albout = -9999.0
            if lis.d.gridDesc[8] != 0.01:
                for i in range(0, lis.d.nch):
                    if (lis.d.gridDesc[3] <= zd_lisdrv['grid'][i].lat <= lis.d.gridDesc[6]
                            and lis.d.gridDesc[4] <= zd_lisdrv['grid'][i].lon <= lis.d.gridDesc[7]):
                        rindex = zd_lisdrv['tile'][i].row - round(
                            (lis.d.gridDesc[3] - lis.d.gridDesc[43]) / lis.d.gridDesc[8])
                        cindex = zd_lisdrv['tile'][i].col - round(
                            (lis.d.gridDesc[4] - lis.d.gridDesc[44]) / lis.d.gridDesc[9])
                        albout[rindex][cindex] = zd_varder['noah'][i].albsf
            else:
                if lis.d.gridDesc[8] == 0.01:
                    for i in range(0, lis.d.nch):
                        albout[zd_lisdrv['tile'][i].row][zd_lisdrv['tile'][i].col] = zd_varder['noah'][i].albsf
            with open("albout.bin", 'w') as f:
                f.write(albout)
            albout = []


# Noah_dynsetup
def noah_dynsetup():
    t, n, ier = 0, 0, 0
    noah_gfrac()
    noah_alb()


# noah_coldstart
def noah_coldstart():
    if lis.o.startcode == 2:
        print('3230 MSG: noah_coldstart -- cold-starting noah', '...using ics from card file', ' (', zd_spmd['iam'],
              ')')
        print('3231 DBG: noah_coldstart -- nch', lis.d.nch, ' (', zd_spmd['iam'], ')')
        for t in range(0, lis.d.nch):
            zd_varder['noah'][t].T1 = noahdrv.noah_it
            # !noah(t)%T1=280.0
            zd_varder['noah'][t].CMC = 0.0004  # !0.04
            zd_varder['noah'][t].SNOWH = 0.0
            zd_varder['noah'][t].SNEQV = 0.0
            zd_varder['noah'][t].CH = 0.0150022404  # !0.0150022404
            zd_varder['noah'][t].CM = 0.0205970779  # !0.0205970779
            for l in range(0, 4):
                zd_varder['noah'][t].STC[l] = noahdrv.noah_it
                zd_varder['noah'][t].SMC[l] = noahdrv.noah_ism
                zd_varder['noah'][t].SH2O[l] = noahdrv.noah_ism

        lis.t.yr = lis.t.syr
        lis.t.mo = lis.t.smo
        lis.t.da = lis.t.sda
        lis.t.hr = lis.t.shr
        lis.t.mn = lis.t.smn
        lis.t.ss = lis.t.sss
        d2tlist = [lis.t.time, lis.t.doy, lis.t.gmt, lis.t.yr, lis.t.mo, lis.t.da, lis.t.hr, lis.t.mn, lis.t.ss]
        date2time(d2tlist)
        lis.t.time = d2tlist[0]
        lis.t.doy = d2tlist[1]
        lis.t.gmt = d2tlist[2]
        print('3256 MSG: noah_coldstart -- Using lis.crd start time ', lis.t.time, ' (', zd_spmd['iam'], ')')
    if lis.o.wout == 4 or lis.o.wout == 5:
        zd_varder['noahbk'] = zd_varder['noah']


# Noah_setup
def noah_setup():
    noah_setvegparms()
    noah_settbot()
    noah_setmxalb()
    noah_setsoils()
    noah_gfrac()
    noah_alb()
    noah_coldstart()
    for t in range(0, lis.d.nch):
        zd_varder['noah'][t].swnet = 0
        zd_varder['noah'][t].lwnet = 0
        zd_varder['noah'][t].qle = 0
        zd_varder['noah'][t].qh = 0
        zd_varder['noah'][t].qg = 0
        zd_varder['noah'][t].snowf = 0
        zd_varder['noah'][t].rainf = 0
        zd_varder['noah'][t].evap = 0
        zd_varder['noah'][t].qs = 0
        zd_varder['noah'][t].qsb = 0
        zd_varder['noah'][t].qsm = 0
        zd_varder['noah'][t].swe = 0
        zd_varder['noah'][t].soilmoist1 = 0
        zd_varder['noah'][t].soilmoist2 = 0
        zd_varder['noah'][t].soilmoist3 = 0
        zd_varder['noah'][t].soilmoist4 = 0
        zd_varder['noah'][t].soilwet = 0
        zd_varder['noah'][t].tveg = 0
        zd_varder['noah'][t].esoil = 0
        zd_varder['noah'][t].rootmoist = 0
        zd_varder['noah'][t].soilm_prev = 0
        zd_varder['noah'][t].swe_prev = 0
        zd_varder['noah'][t].ecanop = 0
        zd_varder['noah'][t].canopint = 0
        zd_varder['noah'][t].COUNT = 0


# makepdsn
def makepdsn(yesterday, beforeyester, kpds, hour, writeint):  # kpda为列表，直接传值。
    if kpds[15] != 0:
        kpds[10] = hour - writeint
        if kpds[10] < 0:
            kpds[10] = 24 - writeint
            kpds[20] = int(beforeyester[0:2])
            kpds[7] = int(beforeyester[2:4])
            kpds[8] = int(beforeyester[4:6])
            kpds[9] = int(beforeyester[6:8])
        else:
            kpds[20] = int(yesterday[0:2])
            kpds[7] = int(yesterday[2:4])
            kpds[8] = int(yesterday[4:6])
            kpds[9] = int(yesterday[6:8])
    else:
        kpds[20] = int(yesterday[0:2])
        kpds[7] = int(yesterday[2:4])
        kpds[8] = int(yesterday[4:6])
        kpds[9] = int(yesterday[6:8])
        kpds[10] = hour
    if kpds[7] == 0:
        kpds[7] = 100
    else:
        kpds[20] = kpds[20] + 1


# stats
def stats(sta):  # var,udef,nch,mean,stdev,min,max
    vsum = 0.0
    sta[3] = 0.0
    dev = 0.0
    sta[4] = 0.0
    sta[5] = 100000.0
    sta[6] = -100000.0
    count = 0
    for t in range(0, sta[2]):
        if sta[0][t] != sta[1]:
            count = count + 1
            vsum = vsum + sta[0][t]
            if sta[0][t] >= sta[6]:
                sta[6] = sta[0][t]
            if sta[0][t] <= sta[5]:
                sta[5] = sta[0][t]
    if vsum == 0:
        sta[6] = 0.0
        sta[5] = 0.0
    if count > 0:
        sta[3] = vsum / float(count)
    else:
        sta[3] = 0
    count = 0
    for t in range(0, sta[2]):
        if sta[0][t] != sta[1]:
            count = count + 1
            dev = dev + pow(sta[0][t] - sta[3], 2)
    if count >= 1:
        sta[3] = pow(dev * pow(float(count) - 1, -1), 0.5)
    else:
        sta[3] = 0


# Noah_binout
def noah_binout(ftn):  # 整数参数（文件代号）更改为文件名参数。
    # ftn, ftn_stats参数为两个整数。
    import math
    LVH2O = 2.501E+6  # Latent heat of vaporization  汽化潜热
    rainf = [0.0 for i in range(lis.d.nch)]
    snowf = [0.0 for i in range(lis.d.nch)]
    for t in range(0, lis.d.nch):
        if float(zd_varder['noah'][t].forcing[0]) < 273.15:
            rainf[t] = 0.0
            snowf[t] = zd_varder['noah'][t].forcing[7]
        else:
            rainf[t] = 0.0
            snowf[t] = 0.0
        zd_varder['noah'][t].swnet = zd_varder['noah'][t].swnet / float(zd_varder['noah'][t].COUNT)
        zd_varder['noah'][t].lwnet = (-1) * zd_varder['noah'][t].lwnet / float(zd_varder['noah'][t].COUNT)
        zd_varder['noah'][t].qle = zd_varder['noah'][t].qle / float(zd_varder['noah'][t].COUNT)
        zd_varder['noah'][t].qh = zd_varder['noah'][t].qh / float(zd_varder['noah'][t].COUNT)
        zd_varder['noah'][t].qg = zd_varder['noah'][t].qg / float(zd_varder['noah'][t].COUNT)
        zd_varder['noah'][t].snowf = zd_varder['noah'][t].snowf / float(zd_varder['noah'][t].COUNT)
        zd_varder['noah'][t].rainf = zd_varder['noah'][t].rainf / float(zd_varder['noah'][t].COUNT)
        zd_varder['noah'][t].evap = zd_varder['noah'][t].evap / float(zd_varder['noah'][t].COUNT)
        zd_varder['noah'][t].qs = zd_varder['noah'][t].qs / float(zd_varder['noah'][t].COUNT)
        zd_varder['noah'][t].qsb = zd_varder['noah'][t].qsb / float(zd_varder['noah'][t].COUNT)
        zd_varder['noah'][t].qsm = zd_varder['noah'][t].qsm / float(zd_varder['noah'][t].COUNT)

        zd_varder['noah'][t].ecanop = zd_varder['noah'][t].ecanop / (LVH2O * float(zd_varder['noah'][t].COUNT))
        zd_varder['noah'][t].tveg = zd_varder['noah'][t].tveg / (LVH2O * float(zd_varder['noah'][t].COUNT))
        zd_varder['noah'][t].esoil = zd_varder['noah'][t].esoil / (LVH2O * float(zd_varder['noah'][t].COUNT))

        zd_varder['noah'][t].soilm_prev = zd_varder['noah'][t].soilmoist1 + zd_varder['noah'][t].soilmoist2 + zd_varder[
            'noah'][t].soilmoist3 + zd_varder['noah'][t].soilmoist4
        zd_varder['noah'][t].swe_prev = zd_varder['noah'][t].swe

    writevar_bin_real(ftn, 'swnet')
    writevar_bin_real(ftn, 'lwnet')
    writevar_bin_real(ftn, 'qle')
    writevar_bin_real(ftn, 'qh')
    writevar_bin_real(ftn, 'qg')
    writevar_bin_real(ftn, 'snowf')
    writevar_bin_real(ftn, 'rainf')
    writevar_bin_real(ftn, 'evap')
    writevar_bin_real(ftn, 'qs')
    writevar_bin_real(ftn, 'qsb')
    writevar_bin_real(ftn, 'qsm')
    writevar_bin_real(ftn, 'SMC')
    writevar_bin_real(ftn, 'SNEQV')
    writevar_bin_real(ftn, 'avgsurft')
    writevar_bin_real(ftn, 'ALBEDO')
    # noah % swe = noah % swe / float(noah % count)
    writevar_bin_real(ftn, 'swe')
    # Subsurface State Variables 地下状态变量
    writevar_bin_real(ftn, 'stc1')
    writevar_bin_real(ftn, 'stc2')
    writevar_bin_real(ftn, 'stc3')
    writevar_bin_real(ftn, 'stc4')
    writevar_bin_real(ftn, 'soilmoist1')
    writevar_bin_real(ftn, 'soilmoist2')
    writevar_bin_real(ftn, 'soilmoist3')
    writevar_bin_real(ftn, 'soilmoist4')
    writevar_bin_real(ftn, 'soilwet')

    writevar_bin_real(ftn, 'ecanop')
    writevar_bin_real(ftn, 'tveg')
    writevar_bin_real(ftn, 'esoil')
    writevar_bin_real(ftn, 'rootmoist')
    writevar_bin_real(ftn, 'canopint')

    #
    # if lis.o.wfor == 1:
    #     if lis.f.force == 1:
    #
    #         writevar_bin_withstats_int(ftn, file2, math.sqrt(
    #             zd_varder['noah'].forcing[4] * zd_varder['noah'].forcing[4] + zd_varder['noah'].forcing[5] *
    #             zd_varder['noah'].forcing[5]), "Wind(m/s)", 1)
    #         writevar_bin_withstats_int(ftn, file2, rainf, "Rainfforc(kg/m2s)", 2)
    #         writevar_bin_withstats_int(ftn, file2, snowf, "Snowfforc(kg/m2s)", 2)
    #         writevar_bin_withstats_int(ftn, file2, zd_varder['noah'].forcing[0], "Tair(K)", 1)
    #         writevar_bin_withstats_int(ftn, file2, zd_varder['noah'].forcing[1], "Qair(kg/kg)", 2)
    #         writevar_bin_withstats_int(ftn, file2, zd_varder['noah'].forcing[6], "PSurf(Pa)", 1)
    #         writevar_bin_withstats_int(ftn, file2, zd_varder['noah'].forcing[2], "SWdown(W/m2)", 1)
    #         writevar_bin_withstats_int(ftn, file2, zd_varder['noah'].forcing[3], "LWdown(W/m2)", 1)
    #     else:
    #         if lis.f.force == 5:
    #             writevar_bin_withstats_int(file, file2, zd_varder['noah'].forcing[0], "Tair(K)", 1)
    #             writevar_bin_withstats_int(file, file2, zd_varder['noah'].forcing[1], "Qair(kg/kg)", 2)
    #             writevar_bin_withstats_int(file, file2, zd_varder['noah'].forcing[2], "SWdown(W/m2)", 1)
    #             writevar_bin_withstats_int(file, file2, zd_varder['noah'].forcing[3], "LWdown(W/m2)", 1)
    #             writevar_bin_withstats_int(file, file2, zd_varder['noah'].forcing[4], "Wind(m/s)", 1)
    #             writevar_bin_withstats_int(file, file2, zd_varder['noah'].forcing[5], "PSurf(Pa)", 1)
    #             writevar_bin_withstats_int(file, file2, zd_varder['noah'].forcing[6], "Rainfforc(kg/m2s)", 2)


def writevar_bin_real(ftn, var):
    import numpy as np

    LSUBS = 2.83E+6  # Latent heat of sublimation  升华潜热

    if zd_spmd['masterproc']:
        gtmp = np.zeros([lis.d.gnc, lis.d.gnr], dtype=np.float32)
        gtmp1 = np.zeros(lis.d.glbngrid, dtype=np.float32)
        var1 = np.zeros(lis.d.nch, dtype=np.float32)
        var0 = np.zeros(lis.d.nch, dtype=np.float32)
    # allocate(gtmp(lis%d%gnc,lis%d%gnr))
    # allocate(gtmp1(lis%d%glbngrid))
    # endif
    if var == 'swnet':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].swnet * zd_lisdrv['tile'][i].fgrd
    if var == 'lwnet':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].lwnet * zd_lisdrv['tile'][i].fgrd
    if var == 'qle':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].qle * zd_lisdrv['tile'][i].fgrd
    if var == 'qh':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].qh * zd_lisdrv['tile'][i].fgrd
    if var == 'qg':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].qg * zd_lisdrv['tile'][i].fgrd
    if var == 'snowf':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].snowf * zd_lisdrv['tile'][i].fgrd
    if var == 'rainf':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].rainf * zd_lisdrv['tile'][i].fgrd
    if var == 'evap':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].evap * zd_lisdrv['tile'][i].fgrd
    if var == 'qs':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].qs * zd_lisdrv['tile'][i].fgrd
    if var == 'qsb':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].qsb * zd_lisdrv['tile'][i].fgrd
    if var == 'qsm':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].qsm * zd_lisdrv['tile'][i].fgrd
    if var == 'SMC':
        for i in range(0, lis.d.nch):
            var0[i] = (zd_varder['noah'][i].SMC[0] * 1000.0 * 0.1 +
                       zd_varder['noah'][i].SMC[1] * 1000.0 * 0.3 + zd_varder['noah'][i].SMC[2] * 1000.0 * 0.6 +
                       zd_varder['noah'][i].SMC[3] * 1000.0 * 1.0 - zd_varder['noah'][i].soilm_prev) / float(
                zd_varder['noah'][i].COUNT)
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + var0[c] * zd_lisdrv['tile'][i].fgrd
    if var == 'SNEQV':
        for i in range(0, lis.d.nch):
            var0[i] = (zd_varder['noah'][i].SNEQV * 1000.0 - zd_varder['noah'][i].swe_prev) / float(
                zd_varder['noah'][i].COUNT)
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + var0[c] * zd_lisdrv['tile'][i].fgrd

    if var == 'avgsurft':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].avgsurft * zd_lisdrv['tile'][i].fgrd

    if var == 'ALBEDO':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].ALBEDO * zd_lisdrv['tile'][i].fgrd
    # noah % swe = noah % swe / float(noah % count)

    if var == 'swe':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].swe * zd_lisdrv['tile'][i].fgrd

    if var == 'stc0':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].STC[0] * zd_lisdrv['tile'][i].fgrd
    if var == 'stc1':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].STC[1] * zd_lisdrv['tile'][i].fgrd
    if var == 'stc2':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].STC[2] * zd_lisdrv['tile'][i].fgrd
    if var == 'stc3':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].STC[3] * zd_lisdrv['tile'][i].fgrd

    if var == 'soilmoist0':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].soilmoist1 * zd_lisdrv['tile'][i].fgrd
    if var == 'soilmoist1':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].soilmoist2 * zd_lisdrv['tile'][i].fgrd
    if var == 'soilmoist2':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].soilmoist3 * zd_lisdrv['tile'][i].fgrd
    if var == 'soilmoist3':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].soilmoist4 * zd_lisdrv['tile'][i].fgrd

    if var == 'soilwet':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].soilwet * zd_lisdrv['tile'][i].fgrd

    if var == 'ecanop':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].ecanop * zd_lisdrv['tile'][i].fgrd
    if var == 'tveg':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].tveg * zd_lisdrv['tile'][i].fgrd
    if var == 'esoil':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].esoil * zd_lisdrv['tile'][i].fgrd
    if var == 'rootmoist':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].rootmoist * zd_lisdrv['tile'][i].fgrd
    if var == 'canopint':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].canopint * zd_lisdrv['tile'][i].fgrd

    if var == 'SMC0':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].SMC[0] * zd_lisdrv['tile'][i].fgrd
    if var == 'SMC1':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].SMC[1] * zd_lisdrv['tile'][i].fgrd
    if var == 'SMC2':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].SMC[2] * zd_lisdrv['tile'][i].fgrd
    if var == 'SMC3':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].SMC[3] * zd_lisdrv['tile'][i].fgrd

    if var == 'SNOWH':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].SNOWH * zd_lisdrv['tile'][i].fgrd

    if var == 'SH2O1':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].SH2O[0] * zd_lisdrv['tile'][i].fgrd
    if var == 'SH2O2':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].SH2O[1] * zd_lisdrv['tile'][i].fgrd
    if var == 'SH2O3':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].SH2O[2] * zd_lisdrv['tile'][i].fgrd
    if var == 'SH2O4':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].SH2O[3] * zd_lisdrv['tile'][i].fgrd
    if var == 'CH':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].CH * zd_lisdrv['tile'][i].fgrd
    if var == 'CM':
        for i in range(0, lis.d.nch):
            c = zd_lisdrv['tile'][i].index
            var1[c] = var1[c] + zd_varder['noah'][c].CM * zd_lisdrv['tile'][i].fgrd

    # enddo
    # if (defined SPMD)
    # call MPI_GATHERV[var1,gdeltas[iam],\
    #         MPI_#real,gtmp1,gdeltas,goffsets,MPI_#real,0,MPI_COMM_WORLD,ierr)
    # else
    gtmp1 = var1
    # endif
    if zd_spmd['masterproc']:
        gtmp[0:int(lis.d.gnc), 0:int(lis.d.gnr)] = lis.d.udef
        count = 1
        for n in range(0, zd_spmd['npes']):
            for r in range(zd_spmd['r_str'][n], zd_spmd['r_end'][n]):
                for c in range(zd_spmd['c_str'][n], zd_spmd['c_end'][n]):
                    gid = c + (r - 1) * lis.d.gnc
                    ntiles = zd_lisdrv['ntiles_pergrid'][gid]
                    if ntiles != 0:
                        gtmp[c, r] = gtmp1[count]
                        count = count + 1

    ftn.write(gtmp)
    del gtmp
    del gtmp1
    del var1
    # endif


# end def writevar_bin_#real


# Drv_output_mod
def readvar_restart_int(ftn, var):
    gtmp = [ftn for i in range(lis.d.glbnch)]
    count = 0
    for r in range(zd_spmd['r_str'][zd_spmd['iam']], zd_spmd['r_end'][zd_spmd['iam']] + 1):
        for c in range(zd_spmd['c_str'][zd_spmd['iam']], zd_spmd['c_end'][zd_spmd['iam']] + 1):
            gid = c + (r - 1) * lis.d.gnc
            ntiles = zd_lisdrv['ntiles_pergrid'][gid]
            stid = zd_lisdrv['s_tid'][gid - 1]
            for t in range(0, ntiles):
                tid = stid + t - 1
                var[count] = gtmp[tid - 1]
                count = count + 1
    gtmp = []


def writevar_restart_int(file, var):  # ftn参数变为文件名参数，调用时注意
    if zd_spmd['masterproc']:
        gtmp = [0 for i in range(lis.d.glbnch)]
        gtmp1 = [var for i in range(lis.d.glbnch)]
    if zd_spmd['masterproc']:
        count = 1
        for n in range(0, zd_spmd['npes']):
            for r in range(zd_spmd['r_str'][n - 1], zd_spmd['r_end'][n - 1] + 1):
                for c in range(zd_spmd['c_str'][n - 1], zd_spmd['c_end'][n - 1] + 1):
                    gid = c + (r - 1) * lis.d.gnc
                    ntiles = zd_lisdrv['ntiles_pergrid'][gid - 1]
                    stid = zd_lisdrv['s_tid'][gid - 1]
                    for t in range(0, ntiles):
                        tid = stid + t - 1
                        gtmp[tid - 1] = gtmp1[count - 1]
                        count = count + 1
        with open(file, 'r') as f:
            f.write(gtmp)
        f.close()
        gtmp = []
        gtmp1 = []


def tile2grid(t, g, nch, nc, nr, udef, tile):  # g为列表。直接传值。
    zd_lisdrv['tile'] = [tiledec for x in range(nch)]  # type型列表.
    # g = [[0] * nc for i in range(nr)]
    for r in range(0, nr):
        for c in range(0, nc):
            if zd_lisdrv['gindex'][r][c] != -1:
                g[r][c] = udef
    for i in range(0, nch):
        c = tile[i].col
        r = tile[i].row
        g[r][c] = g[r][c] + t[i] * tile[i].fgrd


def drv_writevar_grib(ftn, var, kpds, lismask, writeint, today, yesterday, iret=None):
    global ttmp, gtmp
    kpds_t = [0 for x in range(200)]
    gridDesc = [0 for x in range(200)]
    ngrid = 0
    c = 0
    r = 0
    count = 0
    print('3739 WARNING: This routine is not supported')
    # endrun() #请确认该函数是否有用。
    gridDesc = [0 for x in range(200)]
    gridDesc[0] = round(lis.d.gridDesc[0])
    gridDesc[1] = round(lis.d.gridDesc[1])
    gridDesc[2] = round(lis.d.gridDesc[2])
    gridDesc[3] = round(lis.d.gridDesc[3] * 1000)
    gridDesc[4] = round(lis.d.gridDesc[4] * 1000)
    gridDesc[5] = round(lis.d.gridDesc[5])
    gridDesc[6] = round(lis.d.gridDesc[6] * 1000)
    gridDesc[7] = round(lis.d.gridDesc[7] * 1000)
    gridDesc[8] = round(lis.d.gridDesc[8] * 1000)
    gridDesc[9] = round(lis.d.gridDesc[9] * 1000)
    gridDesc[10] = round(lis.d.gridDesc[10])
    gridDesc[19] = round(lis.d.gridDesc[19])
    for i in range(0, 25):
        kpds_t[i] = kpds[i]  # 该赋值方法是否正确？
    ngrid = lis.d.lnc * lis.d.lnr
    ttmp = [0 for x in range(ngrid)]
    gtmp = [[0] * lis.d.lnc for i in range(lis.d.lnr)]
    count = 0
    tile2grid(var, gtmp, lis.d.glbnch, lis.d.lnc, lis.d.lnr, lis.d.udef, tile)
    lismask1 = [False for i in range(lis.d.lnc * lis.d.lnr)]
    for r in range(0, lis.d.lnr):
        for c in range(0, lis.d.lnc):
            ttmp[count] = gtmp[c][r]  # 赋值方式不对？
            if zd_lisdrv['gindex'][r][c] != -1:
                lismask1[count] = True
            count = count + 1
    makepdsn(today, yesterday, kpds_t, lis.t.hr, writeint)
    if iret != 0:
        print('3770 putgb failed for hour = ', lis.t.hr, ',', 'field =', ttmp)
        # endrun()  #列表内存释放问题，该如何释放？
    ttmp = []
    gtmp = []


def writevar_bin_withstats_int(file, file2, var, mvar, form):
    vmean, vstdev, vmin, vmax = 0.0, 0.0, 0.0, 0.0,
    var1 = [0 for i in range(lis.d.ngrid)]
    if zd_spmd['masterproc']:
        gtmp = [[0.0] * lis.d.gnc for i in range(lis.d.gnr)]
        gtmp1 = [0.0 for i in range(lis.d.glbngrid)]
    for i in range(0, lis.d.nch):
        c = zd_lisdrv['tile'][i].index
        var1[c] = var1[c] + var[i] * zd_lisdrv['tile'][i].fgrd
    gtmp1 = var1
    if zd_spmd['masterproc']:
        gtmp = [[lis.d.udef] * lis.d.gnc for i in range(lis.d.gnr)]
        count = 1
        for n in range(0, zd_spmd['npes']):
            for r in range(zd_spmd['r_str'][n], zd_spmd['r_end'][n] + 1):
                for c in range(zd_spmd['c_str'], zd_spmd['c_end'][n] + 1):
                    gid = c + (r - 1) * lis.d.gnc
                    ntiles = zd_lisdrv['ntiles_pergrid'][gid - 1]
                    if ntiles != 0:
                        gtmp[r][c] = gtmp1[count - 1]
                        count = count + 1
        with open(file, 'r') as f:
            f.write(gtmp)
        stalist = [gtmp, lis.d.udef, lis.d.gnc * lis.d.gnr, vmean, vstdev, vmin, vmax]
        stats(stalist)
        gtmp = stalist[0]
        vmean = stalist[3]
        vstdev = stalist[4]
        vmin = stalist[5]
        vmax = stalist[6]

        if form == 1:
            zfc999 = ' ' + mvar.rjust(18, ' ') + format(vmean, '>14.3f') + format(vstdev, '>14.3f') + format(vmin,
                                                                                                             '>14.3f') + format(
                vmax, '>14.3f')
            with open(file2, 'w') as f:
                f.write(zfc999)
            f.close()
        else:
            if form == 2:
                zfc998 = ' ' + mvar.rjust(18, ' ') + ('%14.3e' % vmean) + ('%14.3e' % vstdev) + ('%14.3e' % vmin) + (
                        '%14.3e' % vmax)
                with open(file2, 'w') as f:
                    f.write(zfc998)
        gtmp = []
        gtmp1 = []


# noahrst
def noahrst():
    RW = 0
    C = 0
    R = 0
    T = 0
    I = 0
    J = 0
    L = 0
    N = 0
    F = 0
    FOUND = 0
    YR = 0
    MO = 0
    DA = 0
    HR = 0
    MN = 0
    SS = 0
    VCLASS = 0
    NC = 0
    NR = 0
    NCH = 0
    RHOICE = 917.0  # ! Density of ICE
    WT1 = 0.0
    WT2 = 0.0
    FILEN = ''
    MKFYRMO = ''
    FNAME = ['' for i in range(80)]
    FBASE = ['' for i in range(80)]
    FSUBS = ['' for i in range(80)]
    FMKDIR = ['' for i in range(80)]
    FTIME = ['' for i in range(10)]
    FYRMODIR = ['' for i in range(80)]
    K = 0
    NOAAIC = 0
    curSec = 0
    NOAAIC = 0  # 常量
    print('3861 DBG: noahrst -- in noahrst', '(', zd_spmd['iam'], ')')
    if lis.o.startcode == 1:
        print('3863 noah restart file used:', noahdrv.noah_rfile)
        f = open(noahdrv.noah_rfile, 'r')
        content = f.readlines()
        vclass = content[0].replace(" ", ",").split(',')[0]
        nc = content[0].replace(" ", ",").split(',')[1]
        nr = content[0].replace(" ", ",").split(',')[2]
        nch = content[0].replace(" ", ",").split(',')[3]
        # vclass,nc,nr,nch读取noahdrv.noah_rfile所表示的文件并赋值给该四个变量。
        if vclass != lis.p.vclass:
            print(noahdrv.noah_rfile, '3881 vegetation class conflict')
            # endrun() 该函数注意弄清是否有作用？.L
        if nc != lis.d.gnc or nr != lis.d.gnr:
            print(noahdrv.noah_rfile, '3884 grid space mismatch - noah halted')
        if nch != lis.d.glbnch:
            print('3877 restart tile space mismatch, halting..')
        readvar_restart_int(40, zd_varder['noah'].T1)
        readvar_restart_int(40, zd_varder['noah'].CMC)
        readvar_restart_int(40, zd_varder['noah'].SNOWH)
        readvar_restart_int(40, zd_varder['noah'].SNEQV)
        for l in range(0, 4):
            readvar_restart_int(40, zd_varder['noah'].STC[1])
        for l in range(0, 4):
            readvar_restart_int(40, zd_varder['noah'].SMC[1])
        for l in range(0, 4):
            readvar_restart_int(40, zd_varder['noah'].SH2O[1])
        readvar_restart_int(40, zd_varder['noah'].CH)
        readvar_restart_int(40, zd_varder['noah'].CM)


# Noah_writerst
def noah_dump_restart(file):  # 函数参数由文件代号数字换成文件名（字符串）
    t = 0
    ierr = 0
    alloc_size = 0
    local_size = 0
    tile_loop = 0
    if zd_spmd['masterproc']:
        file.write(str(lis.p.vclass))
        file.write(str(lis.d.gnc))
        file.write(str(lis.d.gnr))
        file.write(str(lis.d.glbnch))
    writevar_bin_real(file, 'avgsurft')
    writevar_bin_real(file, 'canopint')
    writevar_bin_real(file, 'SNOWH')
    writevar_bin_real(file, 'SNEQV')
    writevar_bin_real(file, 'stc1')
    writevar_bin_real(file, 'stc2')
    writevar_bin_real(file, 'stc3')
    writevar_bin_real(file, 'stc4')
    writevar_bin_real(file, 'SMC0')
    writevar_bin_real(file, 'SMC1')
    writevar_bin_real(file, 'SMC2')
    writevar_bin_real(file, 'SMC3')
    writevar_bin_real(file, 'SH2O1')
    writevar_bin_real(file, 'SH2O2')
    writevar_bin_real(file, 'SH2O3')
    writevar_bin_real(file, 'SH2O4')
    writevar_bin_real(file, 'CH')
    writevar_bin_real(file, 'CM')


def noah_writerst():
    global f, filen
    c = 0
    r = 0
    t = 0
    i = 0
    j = 0
    l = 0
    n = 0
    if lis.t.gmt == 24 - noahdrv.writeintn or lis.t.endtime == 1:
        if zd_spmd['masterproc']:
            print('3935 Noah Active restart written: ', noahdrv.noah_rfile)
            if lis.t.mo < 10:
                strmo = '0' + str(lis.t.mo)
            else:
                strmo = str(lis.t.mo)

            if lis.t.da < 10:
                strda = '0' + str(lis.t.da)
            else:
                strda = str(lis.t.da)

            if lis.t.hr < 10:
                strhr = '0' + str(lis.t.hr)
            else:
                strhr = str(lis.t.hr)

            ftime = str(lis.t.yr) + strmo + strda + strhr
            fname = '/EXP' + str(lis.o.expcode) + '/NOAH' + str(lis.t.yr) + '/' + str(lis.t.yr) + strmo + strda
            fsubs = '.Noahrst'

            fbase = lis.o.odir

            filen = fbase + fname + ftime + fsubs
            fileID = open(filen, 'wb')
        # noah_dump_restart(fileID)
        if zd_spmd['masterproc']:
            fileID.close()
            print('3962 noah archive restart written: ', filen)


# readkpds
def readkpds(file, kpds):
    f = open(file, 'r')
    f.seek(29, 0)  # 将光标移动29个字符，并从该位置往后读。依次读6个字符大小
    kpds[4] = int(f.read(6).replace(" ", ""))
    kpds[5] = int(f.read(6).replace(" ", ""))
    kpds[6] = int(f.read(6).replace(" ", ""))
    kpds[13] = int(f.read(6).replace(" ", ""))
    kpds[14] = int(f.read(6).replace(" ", ""))
    kpds[15] = int(f.read(6).replace(" ", ""))
    kpds[21] = int(f.read(6).replace(" ", ""))
    f.close()
    if kpds[15] != 0:
        kpds[14] = noahdrv.writeintn


# Noah_gribout
def noah_gribout(ftn):
    vmean = 0.0
    vstdev = 0.0
    vmin = 0.0
    vmax = 0.0
    rainf = [0.0 for i in range(lis.d.glbnch)]
    snowf = [0.0 for i in range(lis.d.glbnch)]
    t = 0
    c = 0
    r = 0
    i = 0
    k = 0
    lismask = [[True] * lis.d.lnc for i in range(lis.d.lnr)]
    today = ''
    yesterday = ''
    tod = ['' for i in range(8)]
    yes = ['' for i in range(8)]
    ts = 0
    ts1 = 0
    doy1 = 0
    print('4002 WARNING: This routine is not supported in this version')
    # endrun() 请查明该函数是否有作用？
    interval = noahdrv.writeintn
    hr1 = lis.t.hr
    da1 = lis.t.da
    mo1 = lis.t.mo
    yr1 = lis.t.yr
    mn1 = lis.t.mn
    ss1 = 0
    ts1 = -3600 * 24
    dummygmt = 1.0
    dummytime = 1.0
    tod = list(format(yr1, '>4.0f') + format(mo1, '>2.0f') + format(da1, '>2.0f'))
    for i in range(0, 8):
        if tod[i] == ' ':
            tod[i] = '0'
    today = ''.join(tod)
    ticklist = [dummytime, doy1, dummygmt, yr1, mo1, da1, hr1, mn1, ss1, ts1]
    tick(ticklist)
    dummytime = ticklist[0]
    doy1 = ticklist[1]
    dummygmt = ticklist[2]
    yr1 = ticklist[3]
    mo1 = ticklist[4]
    da1 = ticklist[5]
    hr1 = ticklist[6]
    mn1 = ticklist[7]
    ss1 = ticklist[8]
    ts1 = ticklist[9]

    yes = list(format(yr1, '>4.0f') + format(mo1, '>2.0f') + format(da1, '>2.0f'))
    for i in range(0, 8):
        if yes[i] == ' ':
            yes[i] = '0'
    yesterday = ''.join(yes)
    kpds = [0 for i in range(25)]
    kpds[0] = 221
    kpds[1] = 221
    kpds[3] = 192
    kpds[11] = 0
    kpds[12] = 1
    kpds[16] = int((noahdrv.writeintn * 3600.0) / lis.t.ts)
    kpds[17] = 0
    kpds[18] = 1
    kpds[29] = 0
    kpds[22] = 221
    kpds[23] = 0
    kpds[24] = 0
    file = '/content/drive/My Drive/Noah/KPDS_completenoah.tbl'  # 已设置成云盘路径。
    f = open(file, 'r')
    for k in range(0, 42):
        f.read()
    for c in range(0, lis.d.lnc):
        for r in range(0, lis.d.lnr):
            if zd_lisdrv['gindex'][c][r] > 0:
                lismask[r][c] = True
            else:
                lismask[r][c] = False
    for t in range(0, lis.d.glbnch):
        if zd_varder['noah'][t].forcing[0] < 273.15:
            rainf[t] = 0.0
            snowf[t] = zd_varder['noah'][t].forcing[8]
        else:
            rainf[t] = zd_varder['noah'][t].forcing[8]
            snowf[t] = 0.0
    readkpds(file, kpds)
    zd_varder['noah'].swnet = zd_varder['noah'].swnet / float(zd_varder['noah'].COUNT)
    drv_writevar_grib(ftn, zd_varder['noah'].swnet, kpds, lismask, interval, today, yesterday)
    readkpds(file, kpds)
    zd_varder['noah'].lwnet = (-1) * zd_varder['noah'].lwnet / float(zd_varder['noah'].COUNT)
    drv_writevar_grib(ftn, zd_varder['noah'].lwnet, kpds, lismask, interval, today, yesterday)
    readkpds(file, kpds)
    zd_varder['noah'].qle = zd_varder['noah'].qle / float(zd_varder['noah'].COUNT)
    drv_writevar_grib(ftn, zd_varder['noah'].qle, kpds, lismask, interval, today, yesterday)
    readkpds(file, kpds)
    zd_varder['noah'].qh = zd_varder['noah'].qh / float(zd_varder['noah'].COUNT)
    drv_writevar_grib(ftn, zd_varder['noah'].qh, kpds, lismask, interval, today, yesterday)
    readkpds(file, kpds)
    zd_varder['noah'].qg = zd_varder['noah'].qg / float(zd_varder['noah'].COUNT)
    drv_writevar_grib(ftn, zd_varder['noah'].qg, kpds, lismask, interval, today, yesterday)
    readkpds(file, kpds)
    zd_varder['noah'].snowf = zd_varder['noah'].snowf / float(zd_varder['noah'].COUNT)
    drv_writevar_grib(ftn, zd_varder['noah'].snowf, kpds, lismask, interval, today, yesterday)
    readkpds(file, kpds)
    zd_varder['noah'].rainf = zd_varder['noah'].rainf / float(zd_varder['noah'].COUNT)
    drv_writevar_grib(ftn, zd_varder['noah'].rainf, kpds, lismask, interval, today, yesterday)
    readkpds(file, kpds)
    zd_varder['noah'].evap = zd_varder['noah'].evap / float(zd_varder['noah'].COUNT)
    drv_writevar_grib(ftn, zd_varder['noah'].evap, kpds, lismask, interval, today, yesterday)
    readkpds(file, kpds)
    zd_varder['noah'].qs = zd_varder['noah'].qs / float(zd_varder['noah'].COUNT)
    drv_writevar_grib(ftn, zd_varder['noah'].qs, kpds, lismask, interval, today, yesterday)
    readkpds(file, kpds)
    zd_varder['noah'].qsb = zd_varder['noah'].qsb / float(zd_varder['noah'].COUNT)
    drv_writevar_grib(ftn, zd_varder['noah'].qsb, kpds, lismask, interval, today, yesterday)
    readkpds(file, kpds)
    zd_varder['noah'].qsm = zd_varder['noah'].qsm / float(zd_varder['noah'].COUNT)
    drv_writevar_grib(ftn, zd_varder['noah'].qsm, kpds, lismask, interval, today, yesterday)
    readkpds(file, kpds)
    drv_writevar_grib(ftn, (zd_varder['noah'].SMC[0] * 1000.0 * 0.1 + zd_varder['noah'].SMC[1] * 1000.0 * 0.3 +
                            zd_varder['noah'].SMC[2] * 1000.0 * 0.6 + zd_varder['noah'].SMC[3] * 1000.0 * 1.0 -
                            zd_varder['noah'].soilm_prev) / float(zd_varder['noah'].COUNT), kpds, lismask, interval,
                      today, yesterday)
    readkpds(file, kpds)
    drv_writevar_grib(ftn,
                      (zd_varder['noah'].SNEQV * 1000.0 - zd_varder['noah'].swe_prev) / float(zd_varder['noah'].COUNT),
                      kpds, lismask, interval, today, yesterday)
    readkpds(file, kpds)
    drv_writevar_grib(ftn, zd_varder['noah'].avgsurft, kpds, lismask, interval, today, yesterday)
    readkpds(file, kpds)
    drv_writevar_grib(ftn, zd_varder['noah'].ALBEDO, kpds, lismask, interval, today, yesterday)
    readkpds(file, kpds)
    drv_writevar_grib(ftn, zd_varder['noah'].swe, kpds, lismask, interval, today, yesterday)
    readkpds(file, kpds)
    drv_writevar_grib(ftn, zd_varder['noah'].STC[1], kpds, lismask, interval, today, yesterday)
    readkpds(file, kpds)
    drv_writevar_grib(ftn, zd_varder['noah'].STC[2], kpds, lismask, interval, today, yesterday)
    readkpds(file, kpds)
    drv_writevar_grib(ftn, zd_varder['noah'].STC[3], kpds, lismask, interval, today, yesterday)
    readkpds(file, kpds)
    drv_writevar_grib(ftn, zd_varder['noah'].STC[4], kpds, lismask, interval, today, yesterday)
    readkpds(file, kpds)
    drv_writevar_grib(ftn, zd_varder['noah'].soilmoist1, kpds, lismask, interval, today, yesterday)
    readkpds(file, kpds)
    drv_writevar_grib(ftn, zd_varder['noah'].soilmoist2, kpds, lismask, interval, today, yesterday)
    readkpds(file, kpds)
    drv_writevar_grib(ftn, zd_varder['noah'].soilmoist3, kpds, lismask, interval, today, yesterday)
    readkpds(file, kpds)
    drv_writevar_grib(ftn, zd_varder['noah'].soilmoist4, kpds, lismask, interval, today, yesterday)
    readkpds(file, kpds)
    drv_writevar_grib(ftn, zd_varder['noah'].soilwet, kpds, lismask, interval, today, yesterday)

    readkpds(file, kpds)
    zd_varder['noah'].ecanop = zd_varder['noah'].ecanop / float(zd_varder['noah'].COUNT)
    drv_writevar_grib(ftn, zd_varder['noah'].ecanop, kpds, lismask, interval, today, yesterday)
    readkpds(file, kpds)
    zd_varder['noah'].tveg = zd_varder['noah'].tveg / float(zd_varder['noah'].COUNT)
    drv_writevar_grib(ftn, zd_varder['noah'].tveg, kpds, lismask, interval, today, yesterday)
    readkpds(file, kpds)
    zd_varder['noah'].esoil = zd_varder['noah'].esoil / float(zd_varder['noah'].COUNT)
    drv_writevar_grib(ftn, zd_varder['noah'].esoil, kpds, lismask, interval, today, yesterday)
    readkpds(file, kpds)
    drv_writevar_grib(ftn, zd_varder['noah'].rootmoist, kpds, lismask, interval, today, yesterday)
    readkpds(file, kpds)
    drv_writevar_grib(ftn, zd_varder['noah'].canopint, kpds, lismask, interval, today, yesterday)
    if lis.o.wfor == 1:
        readkpds(file, kpds)
        drv_writevar_grib(ftn, math.sqrt(
            zd_varder['noah'].forcing[4] * zd_varder['noah'].forcing[4] + zd_varder['noah'].forcing[5] *
            zd_varder['noah'].forcing[5]), kpds, lismask, interval, today, yesterday)
        readkpds(file, kpds)
        drv_writevar_grib(ftn, rainf, kpds, lismask, interval, today, yesterday)
        readkpds(file, kpds)
        drv_writevar_grib(ftn, snowf, kpds, lismask, interval, today, yesterday)
        readkpds(file, kpds)
        drv_writevar_grib(ftn, zd_varder['noah'].forcing[0], kpds, lismask, interval, today, yesterday)
        readkpds(file, kpds)
        drv_writevar_grib(ftn, zd_varder['noah'].forcing[1], kpds, lismask, interval, today, yesterday)
        readkpds(file, kpds)
        drv_writevar_grib(ftn, zd_varder['noah'].forcing[6], kpds, lismask, interval, today, yesterday)
        readkpds(file, kpds)
        drv_writevar_grib(ftn, zd_varder['noah'].forcing[2], kpds, lismask, interval, today, yesterday)
        readkpds(file, kpds)
        drv_writevar_grib(ftn, zd_varder['noah'].forcing[3], kpds, lismask, interval, today, yesterday)
    f.close()


# Noah_baciof
zd_baciof = {
    'FD': [0 for x in range(999)],
    'BAOPIS': [0 for i in range(20)]
}


def baopen(lu, cfn, iret):  # 将iret包装为列表
    if lu < 1 or lu > 999:
        iret[0] = 6


def baclose(lu, iret):  # 将iret包装为列表
    if lu < 1 or lu > 999:
        iret[0] = 6
    if iret[0] == 0:
        zd_baciof['FD'][lu - 1] = 0


# noah_totinit
def noah_totinit():
    for t in range(0, lis.d.nch):
        if lis.t.gmt % noahdrv.writeintn == 0:
            zd_varder['noah'][t].soilm_prev = zd_varder['noah'][t].SMC[0] * 1000.0 * 0.1 + zd_varder['noah'][t].SMC[
                1] * 1000.0 * 0.3 + zd_varder['noah'][t].SMC[2] * 1000.0 * 0.6 + zd_varder['noah'][t].SMC[3] * 1000.0
            zd_varder['noah'][t].swe_prev = zd_varder['noah'][t].SNEQV * 1000.0
    for t in range(0, lis.d.nch):
        zd_varder['noah'][t].swnet = 0
        zd_varder['noah'][t].lwnet = 0
        zd_varder['noah'][t].qle = 0
        zd_varder['noah'][t].qh = 0
        zd_varder['noah'][t].qg = 0
        zd_varder['noah'][t].snowf = 0
        zd_varder['noah'][t].rainf = 0
        zd_varder['noah'][t].evap = 0
        zd_varder['noah'][t].qs = 0
        zd_varder['noah'][t].qsb = 0
        zd_varder['noah'][t].qsm = 0
        zd_varder['noah'][t].ecanop = 0
        zd_varder['noah'][t].tveg = 0
        zd_varder['noah'][t].esoil = 0
        zd_varder['noah'][t].COUNT = 0


# Noah_almaout
def noah_almaout():
    global fstats
    t = 0
    c = 0
    r = 0
    m = 0
    i = 0
    n = 0
    iret = 0
    ftn = 0
    ftn_stats = 0
    res = 0
    varids = [0 for x in range(10)]
    kpds = [[0] * 25 for i in range(32)]
    FEXIST = 0
    STATUS = 0
    folder_name = ''

    today = ''
    yesterday = ''
    mkfyrmo = ''
    filengb = ''
    mkfyrmo1 = ''
    fname = ['' for x in range(80)]
    fbase = ['' for x in range(40)]
    #    fmkdir = ['' for x in range(80)]
    ftime = ['' for x in range(8)]
    ftimec = ['' for x in range(4)]
    #   fyrmodir = ['' for x in range(26)]
    ftimeb = ['' for x in range(10)]
    fsubgb = ['' for x in range(9)]
    fbinname = ''
    temp1 = ''
    file = ''
    name = ''
    if lis.t.gmt % noahdrv.writeintn == 0:
        noahdrv.numout = noahdrv.numout + 1
        if lis.t.mo < 10:
            str1 = '0' + str(lis.t.mo)
        else:
            str1 = str(lis.t.mo)
        if lis.t.da < 10:
            str2 = '0' + str(lis.t.da)
        else:
            str2 = str(lis.t.da)

        ftime = str(lis.t.yr) + str1 + str2

        ftimec = str(lis.t.yr)

        fbase = lis.o.odir

        fyrmodir = '/EXP' + str(lis.o.expcode) + '/NOAH' + str(lis.t.yr) + '/' + str(lis.t.yr) + str1 + str2

        fmkdir = 'mkdir '
        mkfyrmo = fmkdir + fbase[0:c - 1] + fyrmodir

        os.system(mkfyrmo.rstrip())  # 批处理程序，系统函数

        if lis.t.hr < 10:
            str3 = '0' + str(lis.t.hr)
        else:
            str3 = str(lis.t.hr)

        ftimeb = '/' + str(lis.t.yr) + str1 + str2 + str3

        if lis.o.wout == 1:
            fsubgb = '.gs4r'  # 列表1-5赋值
        if lis.o.wout == 2:
            fsubgb = '.grb '  # 列表1-5赋值
        if lis.o.wout == 3:
            fsubgb = '.nc  '  # 列表1-5赋值

        os.system('mkdir -p ' + fbase + fyrmodir + ftimeb)

        filengb = fbase + fyrmodir + ftimeb + fsubgb

        if noahdrv.noahopen == 0:
            file = 'Noahstats.dat'  # 是否不需要文个件路径就可以了？
            olist = [name, lis.o.odir, lis.o.expcode, file]
            openfile(olist)
            name = olist[0]
            # lis.o.odir = olist[1]
            # lis.o.expcode = olist[2]
            # file = olist[3]
            # if lis.o.startcode == 1:
            #     f = open(name,'wb') #源程序这里两个打开文件的方式具体哪里不一样？
            # else:

            noahdrv.noahopen = 1
            name1 = name + '/'.rstrip() + file
            fstats = open(name1, 'w+')

        zfc1 = '       Statistical Summary of Noah output for:  ' + str(lis.t.mo) + '/' + str(lis.t.da) + '/' + str(
            lis.t.yr) + ' ' + str(lis.t.hr) + ':' + str(lis.t.mn) + ':' + str(lis.t.ss)

        fstats.write(zfc1)  # 为open(name, 'wb')则对zfc1写入报错
        zfc2 = 'Mean          Stdev          Min           Max'
        fstats.write('\n')
        fstats.write('\n')
        fstats.write(zfc2)  # 隔两行写入。

        if lis.o.wout == 1 or lis.o.wout == 4:
            # if zd_spmd['masterproc'] == True:
            fileid = open(filengb, 'wb')
            noah_binout(fileid)  # lis.o.wout == 1 所以ftn=0
            fileid.close()
        if lis.o.wout == 2:  # 暂时未调试
            baolist = [iret]
            ftn = 56
            baopen(ftn, filengb, baolist)
            iret = baolist[0]

            noah_gribout(ftn)

            balist = [iret]
            baclose(ftn, balist)
            iret = balist[0]
            ftn.close()

        noah_writestats(fstats)
        for t in range(0, int(lis.d.nch)):
            zd_varder['noah'][t].count = 0
        print(name)


def noah_writestats(ftn_stats):
    import numpy as np
    rainf = np.zeros(lis.d.nch, dtype=float)
    snowf = np.zeros(lis.d.nch, dtype=float)
    backup = np.zeros(lis.d.nch, dtype=float)

    for t in range(0, lis.d.nch):
        if zd_varder['noah'][t].forcing[0] < 273.15:
            rainf[t] = 0.0
            snowf[t] = zd_varder['noah'][t].forcing[6]
        else:
            rainf[t] = zd_varder['noah'][t].forcing[6]
            snowf[t] = 0.0

        # ---------------------------------------------------------------------------
        # General Energy Balance Components  一般能量平衡组件
        # ---------------------------------------------------------------------------
    vmean, vstdev, vmin, vmax = 0.0, 0.0, 0.0, 0.0
    for t in range(0, lis.d.nch):
        backup[t] = zd_varder['noah'][t].swnet
    staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
    stats(staresult)
    output = 'SWnet(W/m2)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
    ftn_stats.write(output)

    for t in range(0, lis.d.nch):
        backup[t] = zd_varder['noah'][t].lwnet
    staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
    stats(staresult)
    output = 'lwnet(W/m2)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
    ftn_stats.write(output)

    for t in range(0, lis.d.nch):
        backup[t] = zd_varder['noah'][t].qle
    staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
    stats(staresult)
    output = 'qle(W/m2)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
    ftn_stats.write(output)

    for t in range(0, lis.d.nch):
        backup[t] = zd_varder['noah'][t].qh
    staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
    stats(staresult)
    output = 'qh(W/m2)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
    ftn_stats.write(output)

    for t in range(0, lis.d.nch):
        backup[t] = zd_varder['noah'][t].qg
    staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
    stats(staresult)
    output = 'qg(W/m2)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
    ftn_stats.write(output)
    # ---------------------------------------------------------------------------
    # General Water Balance Components  一般水平衡组件
    # ---------------------------------------------------------------------------
    for t in range(0, lis.d.nch):
        backup[t] = zd_varder['noah'][t].snowf
    staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
    stats(staresult)
    output = 'Snowf(kg/m2s)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
    ftn_stats.write(output)

    for t in range(0, lis.d.nch):
        backup[t] = zd_varder['noah'][t].rainf
    staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
    stats(staresult)
    output = 'Rainf(kg/m2s)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
    ftn_stats.write(output)

    for t in range(0, lis.d.nch):
        backup[t] = zd_varder['noah'][t].evap
    staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
    stats(staresult)
    output = 'Evap(kg/m2s)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
    ftn_stats.write(output)

    for t in range(0, lis.d.nch):
        backup[t] = zd_varder['noah'][t].qs
    staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
    stats(staresult)
    output = 'qs(kg/m2s)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
    ftn_stats.write(output)

    for t in range(0, lis.d.nch):
        backup[t] = zd_varder['noah'][t].qsb
    staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
    stats(staresult)
    output = 'qsb(kg/m2s)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
    ftn_stats.write(output)

    for t in range(0, lis.d.nch):
        backup[t] = (zd_varder['noah'][t].SMC[0] * 1000.0 * 0.1 +
                     zd_varder['noah'][t].SMC[1] * 1000.0 * 0.3 +
                     zd_varder['noah'][t].SMC[2] * 1000.0 * 0.6 +
                     zd_varder['noah'][t].SMC[3] * 1000.0 - zd_varder['noah'][t].soilm_prev) / float(
            zd_varder['noah'][t].COUNT)
    staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
    stats(staresult)
    output = 'DelSoilMoist(kg/m2s)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
    ftn_stats.write(output)

    for t in range(0, lis.d.nch):
        backup[t] = (zd_varder['noah'][t].SNEQV * 1000.0 - zd_varder['noah'][t].swe_prev) / float(
            zd_varder['noah'][t].COUNT)
    staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
    stats(staresult)
    output = 'DelSWE(kg/m2s)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
    ftn_stats.write(output)
    # ---------------------------------------------------------------------------
    # Surface State Variables  表面状态变量
    # ---------------------------------------------------------------------------
    for t in range(0, lis.d.nch):
        backup[t] = zd_varder['noah'][t].avgsurft
    staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
    stats(staresult)
    output = 'AvgSurfT(K)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
    ftn_stats.write(output)

    for t in range(0, lis.d.nch):
        backup[t] = zd_varder['noah'][t].ALBEDO
    staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
    stats(staresult)
    output = 'Albedo(-)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
    ftn_stats.write(output)

    for t in range(0, lis.d.nch):
        backup[t] = zd_varder['noah'][t].swe
    staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
    stats(staresult)
    output = 'SWE(kg/m2)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
    ftn_stats.write(output)
    # ---------------------------------------------------------------------------
    # Subsurface State Variables 地下状态变量
    # ---------------------------------------------------------------------------
    for t in range(0, lis.d.nch):
        backup[t] = zd_varder['noah'][t].soilmoist1
    staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
    stats(staresult)
    output = 'SoilMoist1(kg/m2)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
    ftn_stats.write(output)

    for t in range(0, lis.d.nch):
        backup[t] = zd_varder['noah'][t].soilmoist2
    staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
    stats(staresult)
    output = 'SoilMoist2(kg/m2)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
    ftn_stats.write(output)

    for t in range(0, lis.d.nch):
        backup[t] = zd_varder['noah'][t].soilmoist3
    staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
    stats(staresult)
    output = 'SoilMoist3(kg/m2)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
    ftn_stats.write(output)

    for t in range(0, lis.d.nch):
        backup[t] = zd_varder['noah'][t].soilmoist4
    staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
    stats(staresult)
    output = 'SoilMoist4(kg/m2)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
    ftn_stats.write(output)

    for t in range(0, lis.d.nch):
        backup[t] = zd_varder['noah'][t].STC[0]
    staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
    stats(staresult)
    output = 'SoilTemp1 (K)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
    ftn_stats.write(output)

    for t in range(0, lis.d.nch):
        backup[t] = zd_varder['noah'][t].STC[1]
    staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
    stats(staresult)
    output = 'SoilTemp2 (K)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
    ftn_stats.write(output)

    for t in range(0, lis.d.nch):
        backup[t] = zd_varder['noah'][t].STC[2]
    staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
    stats(staresult)
    output = 'SoilTemp3 (K)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
    ftn_stats.write(output)

    for t in range(0, lis.d.nch):
        backup[t] = zd_varder['noah'][t].STC[3]
    staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
    stats(staresult)
    output = 'SoilTemp4 (K)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
    ftn_stats.write(output)

    for t in range(0, lis.d.nch):
        backup[t] = zd_varder['noah'][t].soilwet
    staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
    stats(staresult)
    output = 'SoilWet(-)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
    ftn_stats.write(output)

    # ---------------------------------------------------------------------------
    # Evaporation Components  蒸发成分
    # ---------------------------------------------------------------------------
    for t in range(0, lis.d.nch):
        backup[t] = zd_varder['noah'][t].ecanop
    staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
    stats(staresult)
    output = 'ECanop(kg/m2s)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
    ftn_stats.write(output)

    for t in range(0, lis.d.nch):
        backup[t] = zd_varder['noah'][t].tveg
    staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
    stats(staresult)
    output = 'tveg(kg/m2s)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
    ftn_stats.write(output)

    for t in range(0, lis.d.nch):
        backup[t] = zd_varder['noah'][t].esoil
    staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
    stats(staresult)
    output = 'esoil(kg/m2s)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
    ftn_stats.write(output)

    for t in range(0, lis.d.nch):
        backup[t] = zd_varder['noah'][t].rootmoist
    staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
    stats(staresult)
    output = 'RootMoist(kg/m2)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
    ftn_stats.write(output)

    for t in range(0, lis.d.nch):
        backup[t] = zd_varder['noah'][t].canopint
    staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
    stats(staresult)
    output = 'CanopInt(kg/m2s)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
    ftn_stats.write(output)

    if lis.o.wfor == 1:
        for t in range(0, lis.d.nch):
            backup[t] = zd_varder['noah'][t].forcing[4]
        staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
        stats(staresult)
        output = 'Wind(m/s)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
        ftn_stats.write(output)

        for t in range(0, lis.d.nch):
            backup[t] = rainf[t]
        staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
        stats(staresult)
        output = 'Rainfforc(kg/m2s)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
        ftn_stats.write(output)

        for t in range(0, lis.d.nch):
            backup[t] = snowf[t]
        staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
        stats(staresult)
        output = 'Snowfforc(kg/m2s)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
        ftn_stats.write(output)

        for t in range(0, lis.d.nch):
            backup[t] = zd_varder['noah'][t].forcing[0]
        staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
        stats(staresult)
        output = 'Tair(K)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
        ftn_stats.write(output)

        for t in range(0, lis.d.nch):
            backup[t] = zd_varder['noah'][t].forcing[1]
        staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
        stats(staresult)
        output = 'Qair(kg/kg)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
        ftn_stats.write(output)

        for t in range(0, lis.d.nch):
            backup[t] = zd_varder['noah'][t].forcing[6]
        staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
        stats(staresult)
        output = 'PSurf(Pa)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
        ftn_stats.write(output)

        for t in range(0, lis.d.nch):
            backup[t] = zd_varder['noah'][t].forcing[2]
        staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
        stats(staresult)
        output = 'SWdown (W/m2)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
        ftn_stats.write(output)

        for t in range(0, lis.d.nch):
            backup[t] = zd_varder['noah'][t].forcing[3]
        staresult = [backup, lis.d.udef, lis.d.glbnch, vmean, vstdev, vmin, vmax]
        stats(staresult)
        output = 'LWdown (W/m2)' + str(staresult[3]) + str(staresult[4]) + str(staresult[5]) + str(staresult[6])
        ftn_stats.write(output)
    # endif


# EOC
# end def noah_writestats

# Noah_output
def noah_singlegather(i, var):
    pass


def noah_singleout(lis, tile, var, i):
    pass


def noah_gather():
    pass


def noah_singleout(lis, tile, var, i):
    pass


def noah_output():
    i = 0
    var = [0.0 for i in range(lis.d.glbnch)]
    if lis.o.wsingle == 1:
        if lis.o.wfor == 0:
            num_vars = 30
        else:
            num_vars = 38
        if lis.t.gmt % noahdrv.writeintn == 0:
            if lis.d.gridDesc(9) == 0.01:
                os.system('/data1/pool/control/traffic-wait.pl')
            for i in range(0, num_vars):
                noah_singlegather(i, var)
                if zd_spmd['masterproc']:
                    noah_singleout(lis, tile, var, i)
            if lis.d.gridDesc(9) == 0.01:
                os.system('/data1/pool/control/traffic-go.pl')
            noah_totinit()
    if lis.t.gmt % noahdrv.writeintn == 0:
        if zd_spmd['npes'] > 1:
            noah_gather()
        if zd_spmd['masterproc']:
            noah_almaout()
        noah_totinit()


# Ncbinarydomain_module
class ncbinarydec:
    def __init__(self):
        self.ncold = 0
        self.nrold = 0
        self.nmif = 0
        self.month = 0
        self.nctime1 = 0.0
        self.nctime2 = 0.0
        self.ncdir = ''
        self.tempreso = 0.0
        self.fillvalue = 0.0
        self.startlat = 0.0
        self.startlon = 0.0
        self.latres = 0.0
        self.lonres = 0.0


ncbinarydec = ncbinarydec()
ncbinarydrv = ncbinarydec
zd_ncdomain_mod = {'mi': int(0), 'intowinsizelat': int(0), 'intowinsizelon': int(0),
                   'formermonth': int(0), 'formerstart': int(0), 'currentmonth': int(0), 'currentstart': int(0),
                   'startlon': int(0), 'countlon': int(0), 'startlat': int(0), 'countlat': int(0), 'gswp_nlon': int(0),
                   'gswp_nlat': int(0), 'time_step': int(0),
                   'n111': [[]],
                   'n121': [[]],
                   'n211': [[]],
                   'n221': [[]],
                   'w111': [[]],
                   'w121': [[]],
                   'w211': [[]],
                   'w221': [[]],
                   'months': [],  # [ 0 for i in range(12)],
                   'pcol': [[]],
                   'prow': [[]],
                   'formerforcing': [[[[]]]],
                   'Qair': [[[]]],
                   'Tair': [[[]]],
                   'Wind': [[[]]],
                   'PSurf': [[[]]],
                   'SWdown': [[[]]],
                   'LWdown': [[[]]],
                   'Prec': [[[]]]}


def removeinvalid(verifydata, countlon, countlat):  # 参数情况：数组(发生变化)、整数、整数
    fillvalue = ncbinarydrv.fillvalue
    # slati, elati, slon, elon,count,sum=0,0,0,0,0,0
    for k in range(0, lis.f.nf):
        for j in range(0, countlat):
            for i in range(0, countlon):
                if verifydata[i][j][k] == fillvalue:
                    slati = j
                    elati = j
                    slon = i
                    elon = i
                    while verifydata[j][i][k] == fillvalue:
                        slati = slati - 1
                        elati = elati + 1
                        slon = slon - 1
                        elon = elon + 1
                        if slati < 1:
                            slati = 1
                        if elati > countlat:
                            elati = countlat
                        if slon < 1:
                            slon = 1
                        if elon > countlon:
                            elon = countlon
                        sum = 0.0
                        count = 0
                        for m in range(slati - 1, elati):
                            for n in range(slon - 1, elon):
                                if verifydata([m][n][k] != fillvalue):
                                    sum = sum + verifydata[m][n][k]
                                    count = count + 1
                        if count > 0:
                            verifydata[j][i][k] = sum / count
    return verifydata


def defnatncbinary(gridDesci):  # 该参数未起作用,作用函数已注释.
    import f90nml
    nml = f90nml.read('D:/NoahDate/NOAH/lis.crd')
    ncbinarydrv.ncdir = nml['ncbinary']['ncbinarydrv']['NCDIR']
    ncbinarydrv.nrold = nml['ncbinary']['ncbinarydrv']['NROLD']
    ncbinarydrv.ncold = nml['ncbinary']['ncbinarydrv']['NCOLD']
    ncbinarydrv.month = nml['ncbinary']['ncbinarydrv']['month']
    ncbinarydrv.tempreso = nml['ncbinary']['ncbinarydrv']['TEMPRESO']
    ncbinarydrv.fillvalue = nml['ncbinary']['ncbinarydrv']['FILLVALUE']
    ncbinarydrv.startlat = nml['ncbinary']['ncbinarydrv']['startlat']
    ncbinarydrv.startlon = nml['ncbinary']['ncbinarydrv']['startlon']
    ncbinarydrv.latres = nml['ncbinary']['ncbinarydrv']['latres']
    ncbinarydrv.lonres = nml['ncbinary']['ncbinarydrv']['lonres']
    print('4774 Using ncbinary forcing')
    print('4775 ncbinary forcing directory :', ncbinarydrv.ncdir)
    ncbinarydrv.nctime1 = 3000.0
    ncbinarydrv.nctime2 = 0.0
    zd_ncdomain_mod['mi'] = ncbinarydrv.ncold * ncbinarydrv.nrold  # 修改字典中的值.


def rminvallatlon(datalon, minlon, lonres, countlon, datalat, minlat, latres,
                  countlat):  # 参数只有两个列表datalat,datalon发生变化，其余未变。
    fillvalue = ncbinarydrv.fillvalue
    if (minlat - int(minlat)) < 0.5:
        minlat1 = int(minlat) - 0.5
    else:
        minlat1 = int(minlat) + 0.5
    if (minlon - int(minlon)) < 0.5:
        minlon1 = int(minlon) - 0.5
    else:
        minlon1 = int(minlon) + 0.5
    for j in range(0, countlat):
        if datalat[j][0] == fillvalue:
            datalat[j][0] = minlat1 + j * latres
    for i in range(0, countlon):
        if datalon[0][i] == fillvalue:
            datalon[0][i] = minlon1 + i * lonres
    return datalat, datalon


def ncbinaryfile(nc):  # nc=[name,gdasdir, yr, mo, da, hr,ncmo]name，ncmo改变。
    # fmon = ['' for i in range(2)]
    yro = nc[2]
    moo = nc[3]
    dao = nc[4]
    hro = nc[5]
    fbase = nc[1]
    if nc[5] >= 8:  # hr>=8
        nc[5] = nc[5] - 8
    else:  # hr<8
        nc[5] = nc[5] - 8 + 24
        if nc[4] > 1:  # da>=1
            nc[4] = nc[4] - 1
        else:  # da=1
            if nc[3] > 1:  # mo>1
                nc[3] = nc[3] - 1
                if nc[3] == 3 or nc[3] == 5 or nc[3] == 7 or nc[3] == 8 or nc[3] == 10 or nc[3] == 12:
                    nc[4] = nc[4] - 1 + 31
                if nc[3] == 4 or nc[3] == 6 or nc[3] == 9 or nc[3] == 11:
                    nc[4] = nc[4] - 1 + 30
                if nc[3] == 2:
                    if (nc[2] % 4) == 0 and (nc[2] % 100) != 0 or (nc[2] % 400) == 0:
                        nc[4] = nc[4] - 1 + 29
                    else:
                        nc[4] = nc[4] - 1 + 28
            else:
                nc[4] = 31
                nc[3] = 12
                nc[2] = nc[2] - 1
    nc[6] = nc[3]
    fyear = str(nc[2])
    if nc[3] < 10:
        # fmon[1] = nc[3]  #将mo转换为字符串。
        # fmon[0] = 0 #fmon列表大小为2.
        fmon1 = str(0) + str(nc[3])
    else:
        fmon1 = str(nc[3])  # fmon为字符型列表。
    fname = fyear + fmon1  # fname 为字符串

    temp = fbase + '/' + fname
    nc[0] = temp  # 如果temp长度大于80，name则截取前80个字符作为新的字符串。
    nc[2] = yro
    nc[3] = moo
    nc[4] = dao
    nc[5] = hro


def getncutctime(yr, mo, da, hr, gcc):  # gcc=[uyr,umo,uda,uhr]前面4个参数未变化。
    isutc = 1  # 如果是utc时间，就要将当前模拟时间减去8，以获得强迫数据的正确时间，设isutc=1
    gcc[0] = yr  # 如果是中国时间，则不用处理，设isutc=0
    gcc[1] = mo
    gcc[2] = da
    gcc[3] = hr
    if isutc == 1:  # utc时间
        if gcc[3] >= 8:  # uhr>=8
            gcc[3] = gcc[3] - 8
        else:  # uhr<8
            gcc[3] = gcc[3] - 8 + 24
            if gcc[2] > 1:  # uda>=1
                gcc[2] = gcc[2] - 1
            else:  # uda=1
                if gcc[1] > 1:  # gcc[1]>1
                    gcc[1] = gcc[1] - 1
                    if (gcc[1] == 1 or gcc[1] == 3 or gcc[1] == 5 or gcc[1] == 7 or gcc[1] == 8 or gcc[1] == 10 or gcc[
                        1] == 12):
                        gcc[2] = gcc[2] - 1 + 31
                    if gcc[1] == 4 or gcc[1] == 6 or gcc[1] == 9 or gcc[1] == 11:
                        gcc[2] = gcc[2] - 1 + 30
                    if gcc[1] == 2:
                        if gcc[0] % 4 == 0 and gcc[0] % 100 != 0 or gcc[0] % 400 == 0:
                            gcc[2] = gcc[2] - 1 + 29
                        else:
                            gcc[2] = gcc[2] - 1 + 28
                else:  # gcc[1]=1
                    gcc[2] = 31
                    gcc[1] = 12
                    gcc[0] = gcc[0] - 1


def getcurrentnc(order, yr, mo, da, hr):  # 该程序中startintolon是否是调用了其他子程序的变量？
    import numpy as np
    uyr = 0
    umo = 0
    uda = 0
    uhr = 0
    k, c = 0, 0
    rminvalid = 1
    startintolat, endintolat = 0, 0
    gcclist = [uyr, umo, uda, uhr]
    getncutctime(yr, mo, da, hr, gcclist)
    uyr = gcclist[0]
    umo = gcclist[1]
    uda = gcclist[2]
    uhr = gcclist[3]
    localforcing = [[[0.0] * lis.f.nf for i in range(zd_ncdomain_mod['countlon'])] for i in
                    range(zd_ncdomain_mod['countlat'])]
    if (uyr % 4) == 0 and (uyr % 100) != 0 or (uyr % 400) == 0:
        zd_ncdomain_mod['months'] = [248, 232, 248, 240, 248, 240, 248, 248, 240, 248, 240, 248]
    else:
        zd_ncdomain_mod['months'] = [248, 224, 248, 240, 248, 240, 248, 248, 240, 248, 240, 248]
    if zd_ncdomain_mod['formermonth'] == umo:
        if lis.t.tstype == 0:
            k = (int(((uda - 1) * 24 + uhr) / ncbinarydrv.tempreso) + 1) - zd_ncdomain_mod['months'][umo - 1] + int(
                (lis.a.win * lis.t.ts) / (ncbinarydrv.tempreso * 3600)) + 1
            c = k + 0
        if lis.t.tstype == 1:
            k = (int(((uda - 1) * 24) / ncbinarydrv.tempreso) + 1) - zd_ncdomain_mod['months'][umo - 1] + int(
                (lis.a.win * lis.t.ts) / (ncbinarydrv.tempreso * 3600)) + 1
            c = k + 7
        print('4910', k, c)
        if k >= 0 and c < lis.a.win:
            for i in range(0, lis.f.nf):
                for r in range(k, c + 1):
                    for x in range(0, zd_ncdomain_mod['countlat']):
                        for y in range(0, zd_ncdomain_mod['countlon']):
                            localforcing[y][x][i] = localforcing[y][x][i] + zd_ncdomain_mod['formerforcing'][i][r][y][x]
                            # zd_ncdomain_mod['formerforcing']的维数lis.f.nf, lis.a.win,countlat,countlon

    if zd_ncdomain_mod['currentmonth'] == umo:
        if lis.t.tstype == 0:
            k = int(((uda - 1) * 24 + uhr) / ncbinarydrv.tempreso) + 1
            c = k + 0
        if lis.t.tstype == 1:
            k = int((uda - 1) * 24 / ncbinarydrv.tempreso) + 1
            c = k + 7
        for r in range(k, c):
            for i in range(0, zd_ncdomain_mod['countlon']):
                for j in range(0, zd_ncdomain_mod['countlat']):
                    localforcing[i][j][0] = localforcing[i][j][0] + zd_ncdomain_mod['Tair'][r][j][i]
                    localforcing[i][j][1] = localforcing[i][j][1] + zd_ncdomain_mod['Qair'][r][j][i]
                    localforcing[i][j][2] = localforcing[i][j][2] + zd_ncdomain_mod['SWdown'][r][j][i]
                    localforcing[i][j][3] = localforcing[i][j][3] + zd_ncdomain_mod['LWdown'][r][j][i]
                    localforcing[i][j][4] = localforcing[i][j][4] + zd_ncdomain_mod['Wind'][r][j][i]
                    localforcing[i][j][5] = localforcing[i][j][5] + zd_ncdomain_mod['PSurf'][r][j][i]
                    localforcing[i][j][6] = localforcing[i][j][6] + zd_ncdomain_mod['Prec'][r][j][i]
    for i in range(0, zd_ncdomain_mod['countlon']):
        for j in range(0, zd_ncdomain_mod['countlat']):
            for x in range(0, lis.f.nf):
                localforcing[i][j][x] = localforcing[i][j][x] / (c - k + 1)
    if rminvalid:
        localforcing = removeinvalid(localforcing, zd_ncdomain_mod['countlon'], zd_ncdomain_mod['countlat'])
    for k in range(0, int(lis.f.nf)):
        for r in range(0, int(lis.d.lnr)):
            for c in range(0, int(lis.d.lnc)):
                if zd_lisdrv['gindex'][c][r] != -1:
                    w1 = 0.0
                    startintolat = int(r / ncbinarydrv.latres * lis.d.gridDesc[8]) + 1 - (
                            zd_ncdomain_mod['intowinsizelat'] / 2 - 1)
                    endintolat = int(r / ncbinarydrv.latres * lis.d.gridDesc[8]) + 2 + zd_ncdomain_mod[
                        'intowinsizelat'] / 2 - 1
                    if startintolat < 0:
                        startintolat = 0
                    if endintolat >= zd_ncdomain_mod['countlat']:
                        endintolat = int(zd_ncdomain_mod['countlat']) - 1
                    startintolon = int(c / ncbinarydrv.lonres * lis.d.gridDesc[9]) + 1 - (
                            zd_ncdomain_mod['intowinsizelon'] / 2 - 1)
                    endintolon = int(c / ncbinarydrv.lonres * lis.d.gridDesc[9]) + 2 + zd_ncdomain_mod[
                        'intowinsizelon'] / 2 - 1
                    if startintolon < 0:
                        startintolon = 0
                    if endintolon >= zd_ncdomain_mod['countlon']:
                        endintolon = int(zd_ncdomain_mod['countlon']) - 1
                    for j in range(int(endintolat), int(startintolat) - 1, -1):  # endintolat, startintolat, -1
                        for i in range(int(startintolon), int(endintolon) + 1):  # endintolon
                            w1 = w1 + zd_ncdomain_mod['prow'][j - int(startintolat)][r] * \
                                 localforcing[j - int(startintolat)][i - int(startintolon)][k] * \
                                 zd_ncdomain_mod['pcol'][i - int(startintolon)][c]
                    kk = zd_lisdrv['gindex'][c][r]
                    # print('4031', type(kk))
                    if order == 1:
                        print('4971', np.array(zd_basefor_mod['glbdata1']).shape, kk, k, w1)
                        zd_basefor_mod['glbdata1'][kk][k] = w1
                    else:
                        zd_basefor_mod['glbdata2'][kk][k] = w1
    print('4981 ncbinarydomain_module')
    del localforcing
    # for i in range(0,zd_ncdomain_mod['countlon']):
    #     for j in range(0,zd_ncdomain_mod['countlat']):
    #         print(zd_ncdomain_mod['Tair'][1][j][i],zd_ncdomain_mod['SWdown'][1][j][i])


def currenttoformer():  # 实现对字典formerforcing列表的更新
    global value
    for j in range(1, lis.f.nf + 1):
        for i in range(0, lis.a.win):
            for r in range(0, zd_ncdomain_mod['countlat']):
                for c in range(0, zd_ncdomain_mod['countlon']):
                    if j == 1:
                        value = zd_ncdomain_mod['Tair'][zd_ncdomain_mod['time_step'] - lis.a.win + i][r][c]
                    if j == 2:
                        value = zd_ncdomain_mod['Qair'][zd_ncdomain_mod['time_step'] - lis.a.win + i][r][c]
                    if j == 3:
                        value = zd_ncdomain_mod['SWdown'][zd_ncdomain_mod['time_step'] - lis.a.win + i][r][c]
                    if j == 4:
                        value = zd_ncdomain_mod['LWdown'][zd_ncdomain_mod['time_step'] - lis.a.win + i][r][c]
                    if j == 5:
                        value = zd_ncdomain_mod['Wind'][zd_ncdomain_mod['time_step'] - lis.a.win + i][r][c]
                    if j == 6:
                        value = zd_ncdomain_mod['PSurf'][zd_ncdomain_mod['time_step'] - lis.a.win + i][r][c]
                    if j == 7:
                        value = zd_ncdomain_mod['Prec'][zd_ncdomain_mod['time_step'] - lis.a.win + i][r][c]
                    zd_ncdomain_mod['formerforcing'][j - 1][i][r][c] = value


def readcyclencbinary(yr, mo, da, hr):
    import struct
    import numpy as np
    ncmo = 0
    ffile = ''
    temp = []
    minlon = 0.0
    minlat = 0.0
    startintolon = 0
    endintolon = 0
    Lati = [[]]
    Long = [[]]
    rminvalid = 1
    outfile1 = "D:/NoahDate/noahout/parameter/weightnew"
    outfile2 = "D:/NoahDate/noahout/parameter/rcnew"
    nclist = [ffile, ncbinarydrv.ncdir, yr, mo, da, hr, ncmo]
    ncbinaryfile(nclist)
    ffile = nclist[0]  # 只对变化的量更新。
    ncmo = nclist[6]
    # fileid=12
    filename = ffile
    # pathDir = '/content/drive/My Drive/Noah'
    # 判断文件是否存在
    # if os.path.exists(pathDir + '/' + filename): #文件存在
    if os.path.exists(filename):  # 文件存在
        binFile = open(ffile, 'rb')  # 该文件打开代号为fileid 12。切记
        # content = f.readlines()
        if ncbinarydrv.month == -1:
            context = binFile.read(4 * 6)
            # print('2158',context)

            for i in range(0, 6):
                temp.append(struct.unpack("6f", context)[i])
            context = binFile.read(4)
            # print('4092',struct.unpack(context[6]))
            zd_ncdomain_mod['time_step'] = int(struct.unpack("i", context)[0])
            lstartlat = 1
            lstartlon = 1
            zd_ncdomain_mod['gswp_nlon'] = int((temp[1] - temp[0]) / temp[4]) + 1  # 弄清列表中的元素是否是整数。
            zd_ncdomain_mod['gswp_nlat'] = int((temp[3] - temp[2]) / temp[5])
            nc_startlat = temp[2]  # ncbinarydrv.startlat
            nc_startlon = temp[0]  # ncbinarydrv.startlon
            # the range of the simulation
            maxlat = lis.d.gridDesc[6]
            minlat = lis.d.gridDesc[3]
            maxlon = lis.d.gridDesc[7]
            minlon = lis.d.gridDesc[4]
            if minlat - int(minlat) > 0.5:
                zd_ncdomain_mod['startlat'] = int((minlat - nc_startlat) / temp[4]) + 1
            else:
                zd_ncdomain_mod['startlat'] = int((minlat - nc_startlat) / temp[4])
            if maxlat - int(maxlat) > 0.5:
                zd_ncdomain_mod['countlat'] = int((maxlat - nc_startlat) / temp[4]) + 2 - zd_ncdomain_mod[
                    'startlat'] + 1
            else:
                zd_ncdomain_mod['countlat'] = int((maxlat - nc_startlat) / temp[4]) + 1 - zd_ncdomain_mod[
                    'startlat'] + 1
            if zd_ncdomain_mod['countlat'] > zd_ncdomain_mod['gswp_nlat']:
                zd_ncdomain_mod['countlat'] = zd_ncdomain_mod['gswp_nlat']
            zd_ncdomain_mod['startlon'] = int((minlon - nc_startlon) / temp[5]) + 1
            if maxlon - int(maxlon) > 0.5:
                maxlon = int(maxlon) + 1.5
            else:
                maxlon = int(maxlon) + 0.5
            zd_ncdomain_mod['countlon'] = int((maxlon - nc_startlon) / temp[5]) + 1 - zd_ncdomain_mod['startlon'] + 1
            if zd_ncdomain_mod['countlon'] > zd_ncdomain_mod['gswp_nlon']:
                zd_ncdomain_mod['countlon'] = zd_ncdomain_mod['gswp_nlon']
            print('5072 readncbinary startlon,startlat,countlon,countlat,lonreso,latreso')
            print(zd_ncdomain_mod['startlon'], zd_ncdomain_mod['startlat'], zd_ncdomain_mod['countlon'],
                  zd_ncdomain_mod['countlat'], temp[5], temp[4])
            Lati = [[0.0] * zd_ncdomain_mod['countlon'] for i in range(zd_ncdomain_mod['countlat'])]
            Long = [[0.0] * zd_ncdomain_mod['countlon'] for i in range(zd_ncdomain_mod['countlat'])]
            fileP = 7
            for j in range(lstartlat, zd_ncdomain_mod['countlat'] + 1):
                for i in range(lstartlon, zd_ncdomain_mod['countlon'] + 1):
                    fileP = 7 + (zd_ncdomain_mod['startlat'] - 1 + j - 1) * zd_ncdomain_mod['gswp_nlon'] + \
                            zd_ncdomain_mod['startlon'] - 1 + i
                    binFile.seek(4 * fileP - 4, 0)
                    context = binFile.read(4)
                    Long[j - lstartlat][i - lstartlon] = struct.unpack("f", context)[0]
                    fileP = 7 + (zd_ncdomain_mod['startlat'] - 1 + zd_ncdomain_mod['gswp_nlat'] + j - 1) * \
                            zd_ncdomain_mod['gswp_nlon'] + +zd_ncdomain_mod['startlon'] - 1 + i
                    binFile.seek(4 * fileP - 4, 0)
                    context = binFile.read(4)
                    Lati[j - lstartlat][i - lstartlon] = struct.unpack("f", context)[0]
            zd_ncdomain_mod['formerforcing'] = [
                [[[0.0] * zd_ncdomain_mod['countlon'] for p in range(zd_ncdomain_mod['countlat'])] for i in
                 range(lis.a.win)] for j in range(lis.f.nf)]
        if zd_ncdomain_mod['SWdown']:
            # 如果分配内存，则释放。 zd_ncdomain_mod['formerforcing']的维数lis.f.nf, lis.a.win,countlat,countlon
            zd_ncdomain_mod['SWdown'] = []
            zd_ncdomain_mod['LWdown'] = []
            zd_ncdomain_mod['Wind'] = []
            zd_ncdomain_mod['PSurf'] = []
            zd_ncdomain_mod['Qair'] = []
            zd_ncdomain_mod['Tair'] = []
            zd_ncdomain_mod['Prec'] = []
        binFile.seek(24, 0)
        context = binFile.read(4)
        zd_ncdomain_mod['time_step'] = int(struct.unpack("i", context)[0])
        zd_ncdomain_mod['Qair'] = [[[0.0] * zd_ncdomain_mod['countlon'] for i in range(zd_ncdomain_mod['countlat'])] for
                                   i in range(zd_ncdomain_mod['time_step'])]
        zd_ncdomain_mod['Tair'] = [[[0.0] * zd_ncdomain_mod['countlon'] for i in range(zd_ncdomain_mod['countlat'])] for
                                   i in range(zd_ncdomain_mod['time_step'])]
        zd_ncdomain_mod['Wind'] = [[[0.0] * zd_ncdomain_mod['countlon'] for i in range(zd_ncdomain_mod['countlat'])] for
                                   i in range(zd_ncdomain_mod['time_step'])]
        zd_ncdomain_mod['PSurf'] = [[[0.0] * zd_ncdomain_mod['countlon'] for i in range(zd_ncdomain_mod['countlat'])]
                                    for i in range(zd_ncdomain_mod['time_step'])]
        zd_ncdomain_mod['SWdown'] = [
            [[float(0.0)] * zd_ncdomain_mod['countlon'] for i in range(zd_ncdomain_mod['countlat'])] for i in
            range(zd_ncdomain_mod['time_step'])]
        zd_ncdomain_mod['LWdown'] = [[[0.0] * zd_ncdomain_mod['countlon'] for i in range(zd_ncdomain_mod['countlat'])]
                                     for i in range(zd_ncdomain_mod['time_step'])]
        zd_ncdomain_mod['Prec'] = [[[0.0] * zd_ncdomain_mod['countlon'] for i in range(zd_ncdomain_mod['countlat'])] for
                                   i in range(zd_ncdomain_mod['time_step'])]
        for k in range(1, zd_ncdomain_mod['time_step'] + 1):
            for j in range(1, zd_ncdomain_mod['countlat'] + 1):
                for i in range(1, zd_ncdomain_mod['countlon'] + 1):
                    fileP = 7 + ((k - 1 + 2) * zd_ncdomain_mod['gswp_nlat'] + j - 1 + zd_ncdomain_mod['startlat'] - 1) * \
                            zd_ncdomain_mod['gswp_nlon'] + zd_ncdomain_mod['startlon'] - 1 + i
                    binFile.seek(4 * fileP - 4, 0)
                    context = binFile.read(4)
                    zd_ncdomain_mod['SWdown'][k - 1][j - 1][i - 1] = struct.unpack("f", context)[0]

                    fileP = 7 + ((k - 1 + 2 + zd_ncdomain_mod['time_step']) * zd_ncdomain_mod['gswp_nlat'] + j - 1 +
                                 zd_ncdomain_mod['startlat'] - 1) * zd_ncdomain_mod['gswp_nlon'] + zd_ncdomain_mod[
                                'startlon'] - 1 + i
                    binFile.seek(4 * fileP - 4, 0)
                    context = binFile.read(4)
                    zd_ncdomain_mod['LWdown'][k - 1][j - 1][i - 1] = struct.unpack("f", context)[0]

                    fileP = 7 + ((k - 1 + 2 + 2 * zd_ncdomain_mod['time_step']) * zd_ncdomain_mod['gswp_nlat'] + j - 1 +
                                 zd_ncdomain_mod['startlat'] - 1) * zd_ncdomain_mod['gswp_nlon'] + zd_ncdomain_mod[
                                'startlon'] - 1 + i
                    binFile.seek(4 * fileP - 4, 0)
                    context = binFile.read(4)
                    zd_ncdomain_mod['Wind'][k - 1][j - 1][i - 1] = struct.unpack("f", context)[0]

                    fileP = 7 + ((k - 1 + 2 + 3 * zd_ncdomain_mod['time_step']) * zd_ncdomain_mod['gswp_nlat'] + j - 1 +
                                 zd_ncdomain_mod['startlat'] - 1) * zd_ncdomain_mod['gswp_nlon'] + zd_ncdomain_mod[
                                'startlon'] - 1 + i
                    binFile.seek(4 * fileP - 4, 0)
                    context = binFile.read(4)
                    zd_ncdomain_mod['PSurf'][k - 1][j - 1][i - 1] = struct.unpack("f", context)[0]

                    fileP = 7 + ((k - 1 + 2 + 4 * zd_ncdomain_mod['time_step']) * zd_ncdomain_mod['gswp_nlat'] + j - 1 +
                                 zd_ncdomain_mod['startlat'] - 1) * zd_ncdomain_mod['gswp_nlon'] + zd_ncdomain_mod[
                                'startlon'] - 1 + i
                    binFile.seek(4 * fileP - 4, 0)
                    context = binFile.read(4)
                    zd_ncdomain_mod['Qair'][k - 1][j - 1][i - 1] = struct.unpack("f", context)[0]

                    fileP = 7 + ((k - 1 + 2 + 5 * zd_ncdomain_mod['time_step']) * zd_ncdomain_mod['gswp_nlat'] + j - 1 +
                                 zd_ncdomain_mod['startlat'] - 1) * zd_ncdomain_mod['gswp_nlon'] + zd_ncdomain_mod[
                                'startlon'] - 1 + i
                    binFile.seek(4 * fileP - 4, 0)
                    context = binFile.read(4)
                    zd_ncdomain_mod['Tair'][k - 1][j - 1][i - 1] = struct.unpack("f", context)[0]

                    fileP = 7 + ((k - 1 + 2 + 6 * zd_ncdomain_mod['time_step']) * zd_ncdomain_mod['gswp_nlat'] + j - 1 +
                                 zd_ncdomain_mod['startlat'] - 1) * zd_ncdomain_mod['gswp_nlon'] + zd_ncdomain_mod[
                                'startlon'] - 1 + i
                    binFile.seek(4 * fileP - 4, 0)
                    context = binFile.read(4)
                    zd_ncdomain_mod['Prec'][k - 1][j - 1][i - 1] = struct.unpack("f", context)[0]
        binFile.close()
    else:
        print('5172 readcyclencbinary ', ffile, 'do not exist')
    if ncbinarydrv.month == -1:  # 2
        zd_ncdomain_mod['intowinsizelat'] = 4  # the window size from observation to simulation in lat
        zd_ncdomain_mod['intowinsizelon'] = 4  # the window size from observation to simulation in lon
        Long1 = [[0.0] * zd_ncdomain_mod['countlon'] for i in range(1)]
        Lati1 = [[0.0] * 1 for i in range(zd_ncdomain_mod['countlat'])]
        # print('ncbinarydomain_module 1102', zd_ncdomain_mod['countlon'])
        # print('ncbinarydomain_module 1103', Long)
        for i in range(0, zd_ncdomain_mod['countlat']):
            Lati1[i][0] = Lati[i][0]
        for j in range(0, zd_ncdomain_mod['countlon']):
            Long1[0][j] = Long[0][j]

        if rminvalid:
            (Lati1, Long1) = rminvallatlon(Long1, minlon, temp[5], zd_ncdomain_mod['countlon'], Lati1, minlat, temp[4],
                                           zd_ncdomain_mod['countlat'])

        if zd_ncdomain_mod['intowinsizelat'] > lis.d.lnr:
            zd_ncdomain_mod['intowinsizelat'] = int(lis.d.lnr)
        if zd_ncdomain_mod['intowinsizelon'] > lis.d.lnc:
            zd_ncdomain_mod['intowinsizelon'] = int(lis.d.lnc)

        zd_ncdomain_mod['pcol'] = [[1] * int(lis.d.lnc) for i in range(zd_ncdomain_mod['intowinsizelon'])]
        zd_ncdomain_mod['prow'] = [[1] * int(lis.d.lnr) for i in range(zd_ncdomain_mod['intowinsizelon'])]
        reslat = lis.d.gridDesc[8]
        reslon = lis.d.gridDesc[9]
        for c in range(1, int(lis.d.lnc) + 1):
            loncurr = minlon + (c - 1) * reslon
            startintolon = int(c / temp[4] * reslon) + 1 - zd_ncdomain_mod['intowinsizelon'] / 2
            endintolon = int(c / temp[4] * reslon) + zd_ncdomain_mod['intowinsizelon'] / 2
            if startintolon < 1:
                startintolon = 1
                endintolon = startintolon + zd_ncdomain_mod['intowinsizelon'] - 1
            if endintolon >= zd_ncdomain_mod['countlon']:
                endintolon = zd_ncdomain_mod['countlon'] - 1
            #   startintolon = endintolon - zd_ncdomain_mod['intowinsizelon'] + 1

            for k in range(startintolon, endintolon + 1):
                for j in range(startintolon, k - 1 + 1):
                    zd_ncdomain_mod['pcol'][k - startintolon + 1 - 1][c - 1] = \
                        zd_ncdomain_mod['pcol'][k - startintolon + 1 - 1][c - 1] * (
                                loncurr - Long1[0][j - startintolon]) / (
                                Long1[0][k - startintolon] - Long1[0][j - startintolon])
                for j in range(k + 1, endintolon + 1):
                    zd_ncdomain_mod['pcol'][k - startintolon + 1 - 1][c - 1] = \
                        zd_ncdomain_mod['pcol'][k - startintolon + 1 - 1][c - 1] * (
                                loncurr - Long1[0][j - startintolon]) / (
                                Long1[0][k - startintolon] - Long1[0][j - startintolon])

            w2 = 0.0
            for k in range(startintolon, endintolon + 1):
                w2 = w2 + abs(zd_ncdomain_mod['pcol'][k - startintolon + 1 - 1][c - 1])
            for k in range(startintolon, endintolon + 1):
                zd_ncdomain_mod['pcol'][k - startintolon + 1 - 1][c - 1] = \
                    zd_ncdomain_mod['pcol'][k - startintolon + 1 - 1][c - 1] / w2

        for r in range(1, int(lis.d.lnr) + 1):
            latcurr = minlat + (r - 1) * reslat
            startintolat = int(r / temp[5] * reslat) + 1 - zd_ncdomain_mod['intowinsizelat'] / 2
            endintolat = int(r / temp[5] * reslat) + zd_ncdomain_mod['intowinsizelat'] / 2
            if startintolat < 1:
                startintolat = 1
                endintolat = startintolat + zd_ncdomain_mod['intowinsizelat'] - 1
            if endintolat >= zd_ncdomain_mod['countlat']:
                endintolat = zd_ncdomain_mod['countlat'] - 1
            # startintolat=endintolat-zd_ncdomain_mod['intowinsizelat']+1
            sumlagalat = 0.0
            for k in range(startintolat, endintolat + 1):
                for j in range(startintolat, k - 1 + 1):
                    zd_ncdomain_mod['prow'][k - startintolat + 1 - 1][r - 1] = \
                        zd_ncdomain_mod['prow'][k - startintolat + 1 - 1][r - 1] * (
                                latcurr - Lati1[j - startintolat][0]) / (
                                Lati1[k - startintolat][0] - Lati1[j - startintolat][0])
                for j in range(k + 1, endintolat + 1):
                    zd_ncdomain_mod['prow'][k - startintolat + 1 - 1][r - 1] = \
                        zd_ncdomain_mod['prow'][k - startintolat + 1 - 1][r - 1] * (
                                latcurr - Lati1[j - startintolat][0]) / (
                                Lati1[k - startintolat][0] - Lati1[j - startintolat][0])

            w3 = 0.0
            for k in range(startintolat, endintolat + 1):
                w3 = w3 + abs(zd_ncdomain_mod['prow'][k - startintolat + 1 - 1][r - 1])
            for k in range(startintolat, endintolat + 1):
                zd_ncdomain_mod['prow'][k - startintolat + 1 - 1][r - 1] = (zd_ncdomain_mod['prow'][
                    k - startintolat + 1 - 1][r - 1]) / w3

        Lati = []
        Long = []
        Long1 = []
        Lati1 = []
    ncbinarydrv.month = ncmo
    print('5269 readcyclencbinary is done')


def getcyclencbinary():
    time = 0.0
    doy = 0
    gmt = 0.0
    uyr, umo, uda, uhr = 0, 0, 0, 0
    uyr1, umo1, uda1, uhr1 = 0, 0, 0, 0

    yr = lis.tt.yr
    mo = lis.tt.mo
    da = lis.tt.da
    hr = lis.tt.hr
    mn = 0
    ss = 0
    ts = 0
    # tclist=[time,doy,gmt,yr,mo,da,hr,mn,ss,ts]
    # tick(tclist)
    # time=tclist[0]
    # doy=tclist[1]
    # gmt=tclist[2]
    # yr=tclist[3]
    # mo=tclist[4]
    # da=tclist[5]
    # hr=tclist[6]
    # mn=tclist[7]
    # ss=tclist[8]
    # ts=tclist[9]

    gctlist = [uyr, umo, uda, uhr]
    getncutctime(yr, mo, da, hr, gctlist)
    uyr = gctlist[0]
    umo = gctlist[1]
    uda = gctlist[2]
    uhr = gctlist[3]

    if ncbinarydrv.month == -1:
        readcyclencbinary(yr, mo, da, hr)
        currenttoformer()
        ncbinarydrv.month = mo
        zd_ncdomain_mod['currentmonth'] = umo
        zd_ncdomain_mod['formermonth'] = umo
    yr = lis.tt.yr
    mo = lis.tt.mo
    da = lis.tt.da
    hr = lis.tt.hr
    mn = 0
    ss = 0
    ts = (lis.a.win + 1) * lis.t.ts
    tclist2 = [time, doy, gmt, yr, mo, da, hr, mn, ss, ts]
    tick(tclist2)
    time, doy, gmt, yr, mo = tclist2[0], tclist2[1], tclist2[2], tclist2[3], tclist2[4]
    da, hr, mn, ss, ts = tclist2[5], tclist2[6], tclist2[7], tclist2[8], tclist2[9]
    gcclist = [uyr1, umo1, uda1, uhr1]
    getncutctime(yr, mo, da, hr, gcclist)
    uyr1, umo1, uda1, uhr1 = gcclist[0], gcclist[1], gcclist[2], gcclist[3]
    if umo != umo1:
        currenttoformer()
        readcyclencbinary(yr, mo, da, hr)
        zd_ncdomain_mod['formermonth'] = umo
        zd_ncdomain_mod['currentmonth'] = umo1


def getncbinary_timeindex(yr, mo, da, hr, dexlist):  # dexlist=[index] 包装单个参数.index
    uyr, umo, uda, uhr = 0, 0, 0, 0
    gcclist = [uyr, umo, uda, uhr]
    getncutctime(yr, mo, da, hr, gcclist)
    uyr, umo, uda, uhr = gcclist[0], gcclist[1], gcclist[2], gcclist[3]
    dexlist = int(((uda - 1) * 24 + uhr) / ncbinarydrv.tempreso) + 1
    return dexlist


def gindex(i, j):
    pass


def gindex(n, m):
    pass


def smoothweight(smoothdata, countlon, countlat):
    soiltyp = smoothdata
    count = 0
    sum = 0.0
    window = 2
    for j in range(0, countlat):
        for i in range(0, countlon):
            if gindex(i, j) != -1 and soiltyp(i, j) > 0:
                slati = j - window
                elati = j + window
                slon = i - window
                elon = i + window
                if slati < 1:
                    slati = 1
                if elati > countlat:
                    elati = countlat

                if slon < 1:
                    slon = 1
                if elon > countlon:
                    elon = countlon
            sum = 0.0
            count = 0
            for m in range(slati, elati):
                for n in range(slon, elon):
                    if gindex(n, m) != -1 and soiltyp(n, m) > 0:
                        sum = sum + soiltyp(n, m)
                        count = count + 1
            if count > 0:
                sum = sum / count
        smoothdata[i, j] = sum
    del soiltyp


def ncfilename(namelist, gdasdir, yr, mo, da, hr, intlist):
    # namelist=[name1,name2,name3,name4,name5,name6,name7],intlist=[ncmo,yro,moo,dao,hro]
    # 参数包装两个列表
    # 7个变量为常量
    base1 = 'temp_ITPCAS-CMFD_V0106_B-01_03hr_010deg_'
    base2 = 'shum_ITPCAS-CMFD_V0106_B-01_03hr_010deg_'
    base3 = 'srad_ITPCAS-CMFD_V0106_B-01_03hr_010deg_'
    base4 = 'lrad_ITPCAS-CMFD_V0106_B-01_03hr_010deg_'
    base5 = 'wind_ITPCAS-CMFD_V0106_B-01_03hr_010deg_'
    base6 = 'pres_ITPCAS-CMFD_V0106_B-01_03hr_010deg_'
    base7 = 'prec_ITPCAS-CMFD_V0106_B-01_03hr_010deg_'
    fbase = list(gdasdir.rjust(80, ' '))
    intlist[1] = yr
    intlist[2] = mo
    intlist[3] = da
    intlist[4] = hr
    fmon = ['' for x in range(2)]
    if intlist[1] >= 8:
        intlist[1] = intlist[1] - 8
    else:
        intlist[1] = intlist[1] - 8 + 24
        if intlist[3] > 1:
            intlist[3] = intlist[3] - 1
        else:
            if intlist[2] > 1:
                intlist[2] = intlist[2] - 1
                if intlist[2] == 3 or intlist[2] == 5 or intlist[2] == 7 or intlist[2] == 8 or intlist[2] == 10 or \
                        intlist[2] == 12:
                    intlist[3] = intlist[3] - 1 + 31
                if intlist[2] == 4 or intlist[2] == 6 or intlist[2] == 9 or intlist[2] == 11:
                    intlist[3] = intlist[3] - 1 + 30
                if intlist[2] == 2:
                    if intlist[1] / 4 == 0 and intlist[1] / 100 != 0 or intlist[1] / 400 == 0:
                        intlist[3] = intlist[3] - 1 + 29
                    else:
                        intlist[3] = intlist[3] - 1 + 28
            else:
                intlist[3] = 31
                intlist[2] = 12
                intlist[1] = intlist[1] - 1
    intlist[0] = intlist[2]
    fyear = format(intlist[1], '>4.0f')
    if intlist[2] < 10:
        fmon[0] = '0'
        fmon[1] = format(intlist[2], '>1.0f')
    else:
        fmon = list(format(intlist[2], '>2.0f'))
    fnameb = list(fyear + ''.join(fmon) + '.nc')
    c = 0
    for i in range(0, 79):
        if fbase[i] == ' ' and c == 0:
            c = i - 1
    namelist[0] = ''.join(fbase[0:c + 1]) + ''.rstrip() + base1 + ''.join(fnameb[0:10])
    namelist[1] = ''.join(fbase[0:c + 1]) + ''.rstrip() + base2 + ''.join(fnameb[0:10])
    namelist[2] = ''.join(fbase[0:c + 1]) + ''.rstrip() + base3 + ''.join(fnameb[0:10])
    namelist[3] = ''.join(fbase[0:c + 1]) + ''.rstrip() + base4 + ''.join(fnameb[0:10])
    namelist[4] = ''.join(fbase[0:c + 1]) + ''.rstrip() + base5 + ''.join(fnameb[0:10])
    namelist[5] = ''.join(fbase[0:c + 1]) + ''.rstrip() + base6 + ''.join(fnameb[0:10])
    namelist[6] = ''.join(fbase[0:c + 1]) + ''.rstrip() + base7 + ''.join(fnameb[0:10])


def freencbinary():
    if zd_ncdomain_mod['n111']:
        zd_ncdomain_mod['n111'] = []
        zd_ncdomain_mod['n121'] = []
        zd_ncdomain_mod['n211'] = []
        zd_ncdomain_mod['n221'] = []
        zd_ncdomain_mod['w111'] = []
        zd_ncdomain_mod['w121'] = []
        zd_ncdomain_mod['w211'] = []
        zd_ncdomain_mod['w221'] = []
    if zd_ncdomain_mod['SWdown']:
        zd_ncdomain_mod['SWdown'] = []
        zd_ncdomain_mod['LWdown'] = []
        zd_ncdomain_mod['Wind'] = []
        zd_ncdomain_mod['PSurf'] = []
        zd_ncdomain_mod['Qair'] = []
        zd_ncdomain_mod['Tair'] = []
        zd_ncdomain_mod['Prec'] = []
    if zd_ncdomain_mod['pcol']:
        zd_ncdomain_mod['pcol'] = []
        zd_ncdomain_mod['prow'] = []
    if zd_ncdomain_mod['formerforcing']:
        zd_ncdomain_mod['formerforcing'] = []


# getncbinary
def getncbinary():
    if lis.a.daalg > 0:
        getncbinaryda()
    else:
        getncbinarynonda()


def getncbinaryda():
    # EOP
    timenow = 0.0
    if zd_spmd['masterproc']:
        nstep = get_nstep(lis.t)
    # endif

    # -------------------------------------------------------------------
    # Determine the correct number of forcing variables 确定强制变量的正确数量
    # -------------------------------------------------------------------
    if nstep == 0:
        nforce = lis.f.nf
    else:
        nforce = lis.f.nf
    # endif
    lis.f.findtime1 = 0
    lis.f.findtime2 = 0
    movetime = 0
    # -------------------------------------------------------------------
    # Determine Required GSWP Data Times 确定所需的GSWP数据时间
    # (The previous hour & the future hour) 之前的时间和之后的时间
    # -------------------------------------------------------------------
    if lis.t.tstype == 0:  # simulate at a month resolution 以月分辨率模拟
        yr1 = lis.t.yr  # Previous Hour 之前的时间
        mo1 = lis.t.mo
        da1 = lis.t.da
        hr1 = int(lis.t.hr / 3) * 3 + 1  # hr1=3*((lis%t%hr)/3)
        mn1 = 0
        ss1 = 0
        ts1 = 0

        time1, doy1, gmt1 = 0.0, 0.0, 0.0
        tklist = [time1, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1]
        tick(tklist)
        time1, doy1, gmt1 = tklist[0], tklist[1], tklist[2]

        order = 1
        hr1 = lis.t.hr
        getcurrentnc(order, yr1, mo1, da1, hr1)  # line 158 in ncbinarydomain_module ncbinarydomain_模块中的第158行
        ncbinarydrv.nctime1 = time1

        order = 2
        nstep = int(hr1 / 3) * 3 + 1
        if nstep <= hr1:
            ts1 = lis.t.ts
            tklist = [timenow, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1]
            tick(tklist)
            timenow, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1 = tklist[0], tklist[1], tklist[2], tklist[3], tklist[
                4], tklist[5], tklist[6], tklist[7], tklist[8], tklist[9]
            getcurrentnc(order, yr1, mo1, da1, hr1)
        else:
            if hr1 >= 1:  # hr1>=1
                hr1 = hr1 - 1
            else:  # uhr<1
                hr1 = hr1 - 1 + 24
                if da1 > 1:  # da1>=1
                    da1 = da1 - 1
                else:  # uda=1
                    if mo1 > 1:  # mo1>1
                        mo1 = mo1 - 1
                        if (mo1 == 1 or mo1 == 3 or mo1 == 5 or mo1 == 7
                                or mo1 == 8 or mo1 == 10 or mo1 == 12):
                            da1 = da1 - 1 + 31
                        if mo1 == 4 or mo1 == 6 or mo1 == 9 or mo1 == 11:
                            da1 = da1 - 1 + 30
                        if mo1 == 2:
                            if (yr1 % 4 == 0 and yr1 % 100 != 0) or yr1 % 400 == 0:
                                da1 = da1 - 1 + 29
                            else:
                                da1 = da1 - 1 + 28
                    else:  # mo1=1
                        da1 = 31
                        mo1 = 12
                        yr1 = yr1 - 1
            getcurrentnc(order, yr1, mo1, da1, hr1)
        # end if #if(nstep  <=  hr1)
        hr1 = int(hr1 / 3) * 3 + 1
        ts1 = 0
        tklist = [time1, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1]
        tick(tklist)
        ncbinarydrv.nctime2 = tklist[0]
        # end if #if(LIS%t%TSTYPE  ==  0)

    if lis.t.tstype != 0:
        order = 1
        getcurrentnc(order, yr1, mo1, da1, hr1)


# end subroutine getncbinaryda

def getncbinarynonda():
    if zd_spmd['masterproc']:
        nstep = get_nstep(lis.t)
    # endif

    # -------------------------------------------------------------------
    # Determine the correct number of forcing variables 确定强制变量的正确数量
    # -------------------------------------------------------------------
    if nstep == 0:
        nforce = lis.f.nf
    else:
        nforce = lis.f.nf
    # endif
    lis.f.findtime1 = 0
    lis.f.findtime2 = 0
    movetime = 0
    # -------------------------------------------------------------------
    # Determine Required GSWP Data Times 确定所需的GSWP数据时间
    # (The previous hour & the future hour) 之前的小时 之后的小时
    # -------------------------------------------------------------------
    yr1 = lis.t.yr  # Time now 现在的时间
    mo1 = lis.t.mo
    da1 = lis.t.da
    hr1 = lis.t.hr
    mn1 = lis.t.mn
    ss1 = 0
    ts1 = 0
    timenow, doy1, gmt1 = 0.0, 0.0, 0.0
    if lis.t.tstype == 0:  # simulate at a month resolution 以月分辨率模拟

        tklist = [timenow, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1]
        tick(tklist)
        timenow, doy1, gmt1 = tklist[0], tklist[1], tklist[2]

        yr1 = lis.t.yr  # Previous Hour 之前的小时
        mo1 = lis.t.mo
        da1 = lis.t.da
        hr1 = int(lis.t.hr / 3) * 3  # hr1=3*((lis.t.hr)/3)
        mn1 = 0
        ss1 = 0
        ts1 = 0

        time1 = 0.0
        tklist2 = [time1, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1]
        tick(tklist2)
        time1, doy1, gmt1 = tklist2[0], tklist2[1], tklist2[2]

        yr2 = lis.t.yr  # Next Hour 下一个小时
        mo2 = lis.t.mo
        da2 = lis.t.da
        hr2 = int(lis.t.hr / 3) * 3
        mn2 = 0
        ss2 = 0
        ts2 = 3 * 60 * 60

        time2, doy2, gmt2 = 0.0, 0.0, 0.0
        tklist3 = [time2, doy2, gmt2, yr2, mo2, da2, hr2, mn2, ss2, ts2]
        tick(tklist3)
        time2 = tklist3[0]
        doy2 = tklist3[1]
        gmt2 = tklist3[2]
        yr2 = tklist3[3]
        mo2 = tklist3[4]
        da2 = tklist3[5]
        hr2 = tklist3[6]
        mn2 = tklist3[7]
        ss2 = tklist3[8]
        ts2 = tklist3[9]

        if timenow > ncbinarydrv.nctime2:
            movetime = 1
            lis.f.findtime2 = 1
        # endif

        if nstep == 0 or nstep == 1 or lis.f.rstflag == 1:
            lis.f.findtime1 = 1
            lis.f.findtime2 = 1
            # glbdata1 = 0
            # glbdata2 = 0
            movetime = 0
            lis.f.rstflag = 0
        # endif
        # print*,'140',lis%f%findtime1, movetime
        if lis.f.findtime1 == 1:
            #    print*, 'reading time1 data..'
            order = 1
            readncbinary(order, yr1, mo1, da1, hr1)
            ncbinarydrv.nctime1 = time1
        # endif

        if movetime == 1:
            ncbinarydrv.nctime1 = ncbinarydrv.nctime2
            lis.f.findtime2 = 1
            for f in range(0, nforce):
                for c in range(0, lis.d.ngrid):
                    zd_basefor_mod['glbdata1'][c][f] = zd_basefor_mod['glbdata2'][c][f]
                # enddo
            # enddo
        # endif
        if lis.f.findtime2 == 1:
            #    print*, 'reading time2 data..'
            order = 2
            readncbinary(order, yr2, mo2, da2, hr2)
            ncbinarydrv.nctime2 = time2
        # endif

    # endif #if(LIS%t%TSTYPE  ==  0)
    if lis.t.tstype != 0:
        order = 1
        readncbinary(order, yr1, mo1, da1, hr1)
    # endif


# end subroutine getncbinarynonda 结束子例程getncbinarynonda

# readncfile
def deallocate(Long1):
    pass


def deallocate(Lati1):
    pass


def readncbinary(order, yr, mo, da, hr):
    global Long1
    import struct
    import numpy as np
    ncmo = 0
    ffile = ''
    temp = []
    minlon = 0.0
    minlat = 0.0
    startintolon = 0
    endintolon = 0
    Lati = [[]]
    Long = [[]]
    rminvalid = 1
    outfile1 = "D:/NoahDate/noahout/parameter/weightnew"
    outfile2 = "D:/NoahDate/noahout/parameter/rcnew"
    nclist = [ffile, ncbinarydrv.ncdir, yr, mo, da, hr, ncmo]
    ncbinaryfile(nclist)
    ffile = nclist[0]  # 只对变化的量更新。
    ncmo = nclist[6]
    # fileid=12
    filename = ffile
    #  pathDir = '/content/drive/My Drive/Noah'
    # 判断文件是否存在
    #  if os.path.exists(pathDir + '/' + filename): #文件存在
    if ncbinarydrv.month != ncmo:
        starttime = 1
        fileid = 12
        if os.path.exists(filename):  # 文件存在
            binFile = open(ffile, 'rb')  # 该文件打开代号为fileid 12。切记
            # content = f.readlines()
            if ncbinarydrv.month == -1:
                context = binFile.read(4 * 6)
                # print('2158',context)

                for i in range(0, 6):
                    temp.append(struct.unpack("6f", context)[i])
                context = binFile.read(4)
                # print('4092',struct.unpack(context[6]))
                zd_ncdomain_mod['time_step'] = int(struct.unpack("i", context)[0])
                lstartlat = 1
                lstartlon = 1
                zd_ncdomain_mod['gswp_nlon'] = int((temp[1] - temp[0]) / temp[4]) + 1  # 弄清列表中的元素是否是整数。
                zd_ncdomain_mod['gswp_nlat'] = int((temp[3] - temp[2]) / temp[5])
                nc_startlat = temp[2]  # ncbinarydrv.startlat #
                nc_startlon = temp[0]  # ncbinarydrv.startlon #
                # the range of the simulation 模拟的范围
                maxlat = lis.d.gridDesc[6]
                minlat = lis.d.gridDesc[3]
                maxlon = lis.d.gridDesc[7]
                minlon = lis.d.gridDesc[4]
                if minlat - int(minlat) > 0.5:
                    zd_ncdomain_mod['startlat'] = int((minlat - nc_startlat) / temp[4]) + 1
                else:
                    zd_ncdomain_mod['startlat'] = int((minlat - nc_startlat) / temp[4])
                if maxlat - int(maxlat) > 0.5:
                    zd_ncdomain_mod['countlat'] = int((maxlat - nc_startlat) / temp[4]) + 2 - zd_ncdomain_mod[
                        'startlat'] + 1
                else:
                    zd_ncdomain_mod['countlat'] = int((maxlat - nc_startlat) / temp[4]) + 1 - zd_ncdomain_mod[
                        'startlat'] + 1
                if zd_ncdomain_mod['countlat'] > zd_ncdomain_mod['gswp_nlat']:
                    zd_ncdomain_mod['countlat'] = zd_ncdomain_mod['gswp_nlat']
                zd_ncdomain_mod['startlon'] = int((minlon - nc_startlon) / temp[5]) + 1
                if maxlon - int(maxlon) > 0.5:
                    maxlon = int(maxlon) + 1.5
                else:
                    maxlon = int(maxlon) + 0.5
                zd_ncdomain_mod['countlon'] = int((maxlon - nc_startlon) / temp[5]) + 1 - zd_ncdomain_mod[
                    'startlon'] + 1
                if zd_ncdomain_mod['countlon'] > zd_ncdomain_mod['gswp_nlon']:
                    zd_ncdomain_mod['countlon'] = zd_ncdomain_mod['gswp_nlon']
                # print('4757 readncbinary startlon,startlat,countlon,countlat,lonreso,latreso')
                # print(zd_ncdomain_mod['startlon'], zd_ncdomain_mod['startlat'], zd_ncdomain_mod['countlon'],\
                #       zd_ncdomain_mod['countlat'], temp[5], temp[4])
                Lati = [[0.0] * zd_ncdomain_mod['countlon'] for i in range(zd_ncdomain_mod['countlat'])]
                Long = [[0.0] * zd_ncdomain_mod['countlon'] for i in range(zd_ncdomain_mod['countlat'])]
                fileP = 7
                for j in range(lstartlat, zd_ncdomain_mod['countlat'] + 1):
                    for i in range(lstartlon, zd_ncdomain_mod['countlon'] + 1):
                        fileP = 7 + (zd_ncdomain_mod['startlat'] - 1 + j - 1) * zd_ncdomain_mod['gswp_nlon'] + \
                                zd_ncdomain_mod['startlon'] - 1 + i
                        binFile.seek(4 * fileP - 4, 0)
                        context = binFile.read(4)
                        Long[j - lstartlat][i - lstartlon] = struct.unpack("f", context)[0]
                        fileP = 7 + (zd_ncdomain_mod['startlat'] - 1 + zd_ncdomain_mod['gswp_nlat'] + j - 1) * \
                                zd_ncdomain_mod['gswp_nlon'] + +zd_ncdomain_mod['startlon'] - 1 + i
                        binFile.seek(4 * fileP - 4, 0)
                        context = binFile.read(4)
                        Lati[j - lstartlat][i - lstartlon] = struct.unpack("f", context)[0]
                # zd_ncdomain_mod['formerforcing'] = [[[[0.0] * zd_ncdomain_mod['countlon']
                # for p in range(zd_ncdomain_mod['countlat'])] for i in range(lis.a.win)] for j in range(lis.f.nf)]
            if ncbinarydrv.month != -1 and ncbinarydrv.month != ncmo:
                # 如果分配内存，则释放。 zd_ncdomain_mod['formerforcing']的维数lis.f.nf, lis.a.win,countlat,countlon
                zd_ncdomain_mod['SWdown'] = []
                zd_ncdomain_mod['LWdown'] = []
                zd_ncdomain_mod['Wind'] = []
                zd_ncdomain_mod['PSurf'] = []
                zd_ncdomain_mod['Qair'] = []
                zd_ncdomain_mod['Tair'] = []
                zd_ncdomain_mod['Prec'] = []
            binFile.seek(24, 0)
            context = binFile.read(4)
            zd_ncdomain_mod['time_step'] = int(struct.unpack("i", context)[0])
            zd_ncdomain_mod['Qair'] = [[[0.0] * zd_ncdomain_mod['countlon'] for i in range(zd_ncdomain_mod['countlat'])]
                                       for i in range(zd_ncdomain_mod['time_step'])]
            zd_ncdomain_mod['Tair'] = [[[0.0] * zd_ncdomain_mod['countlon'] for i in range(zd_ncdomain_mod['countlat'])]
                                       for i in range(zd_ncdomain_mod['time_step'])]
            zd_ncdomain_mod['Wind'] = [[[0.0] * zd_ncdomain_mod['countlon'] for i in range(zd_ncdomain_mod['countlat'])]
                                       for i in range(zd_ncdomain_mod['time_step'])]
            zd_ncdomain_mod['PSurf'] = [
                [[0.0] * zd_ncdomain_mod['countlon'] for i in range(zd_ncdomain_mod['countlat'])] for i in
                range(zd_ncdomain_mod['time_step'])]
            zd_ncdomain_mod['SWdown'] = [
                [[float(0.0)] * zd_ncdomain_mod['countlon'] for i in range(zd_ncdomain_mod['countlat'])] for i in
                range(zd_ncdomain_mod['time_step'])]
            zd_ncdomain_mod['LWdown'] = [
                [[0.0] * zd_ncdomain_mod['countlon'] for i in range(zd_ncdomain_mod['countlat'])] for i in
                range(zd_ncdomain_mod['time_step'])]
            zd_ncdomain_mod['Prec'] = [[[0.0] * zd_ncdomain_mod['countlon'] for i in range(zd_ncdomain_mod['countlat'])]
                                       for i in range(zd_ncdomain_mod['time_step'])]
            for k in range(1, zd_ncdomain_mod['time_step'] + 1):
                for j in range(1, zd_ncdomain_mod['countlat'] + 1):
                    for i in range(1, zd_ncdomain_mod['countlon'] + 1):
                        fileP = 7 + ((k - 1 + 2) * zd_ncdomain_mod['gswp_nlat'] + j - 1 + zd_ncdomain_mod[
                            'startlat'] - 1) * \
                                zd_ncdomain_mod['gswp_nlon'] + zd_ncdomain_mod['startlon'] - 1 + i
                        binFile.seek(4 * fileP - 4, 0)
                        context = binFile.read(4)
                        zd_ncdomain_mod['SWdown'][k - 1][j - 1][i - 1] = struct.unpack("f", context)[0]

                        fileP = 7 + ((k - 1 + 2 + zd_ncdomain_mod['time_step']) * zd_ncdomain_mod['gswp_nlat'] + j - 1 +
                                     zd_ncdomain_mod['startlat'] - 1) * zd_ncdomain_mod['gswp_nlon'] + zd_ncdomain_mod[
                                    'startlon'] - 1 + i
                        binFile.seek(4 * fileP - 4, 0)
                        context = binFile.read(4)
                        zd_ncdomain_mod['LWdown'][k - 1][j - 1][i - 1] = struct.unpack("f", context)[0]

                        fileP = 7 + ((k - 1 + 2 + 2 * zd_ncdomain_mod['time_step']) * zd_ncdomain_mod[
                            'gswp_nlat'] + j - 1 +
                                     zd_ncdomain_mod['startlat'] - 1) * zd_ncdomain_mod['gswp_nlon'] + zd_ncdomain_mod[
                                    'startlon'] - 1 + i
                        binFile.seek(4 * fileP - 4, 0)
                        context = binFile.read(4)
                        zd_ncdomain_mod['Wind'][k - 1][j - 1][i - 1] = struct.unpack("f", context)[0]

                        fileP = 7 + ((k - 1 + 2 + 3 * zd_ncdomain_mod['time_step']) * zd_ncdomain_mod[
                            'gswp_nlat'] + j - 1 +
                                     zd_ncdomain_mod['startlat'] - 1) * zd_ncdomain_mod['gswp_nlon'] + zd_ncdomain_mod[
                                    'startlon'] - 1 + i
                        binFile.seek(4 * fileP - 4, 0)
                        context = binFile.read(4)
                        zd_ncdomain_mod['PSurf'][k - 1][j - 1][i - 1] = struct.unpack("f", context)[0]

                        fileP = 7 + ((k - 1 + 2 + 4 * zd_ncdomain_mod['time_step']) * zd_ncdomain_mod[
                            'gswp_nlat'] + j - 1 +
                                     zd_ncdomain_mod['startlat'] - 1) * zd_ncdomain_mod['gswp_nlon'] + zd_ncdomain_mod[
                                    'startlon'] - 1 + i
                        binFile.seek(4 * fileP - 4, 0)
                        context = binFile.read(4)
                        zd_ncdomain_mod['Qair'][k - 1][j - 1][i - 1] = struct.unpack("f", context)[0]

                        fileP = 7 + ((k - 1 + 2 + 5 * zd_ncdomain_mod['time_step']) * zd_ncdomain_mod[
                            'gswp_nlat'] + j - 1 +
                                     zd_ncdomain_mod['startlat'] - 1) * zd_ncdomain_mod['gswp_nlon'] + zd_ncdomain_mod[
                                    'startlon'] - 1 + i
                        binFile.seek(4 * fileP - 4, 0)
                        context = binFile.read(4)
                        zd_ncdomain_mod['Tair'][k - 1][j - 1][i - 1] = struct.unpack("f", context)[0]

                        fileP = 7 + ((k - 1 + 2 + 6 * zd_ncdomain_mod['time_step']) * zd_ncdomain_mod[
                            'gswp_nlat'] + j - 1 +
                                     zd_ncdomain_mod['startlat'] - 1) * zd_ncdomain_mod['gswp_nlon'] + zd_ncdomain_mod[
                                    'startlon'] - 1 + i
                        binFile.seek(4 * fileP - 4, 0)
                        context = binFile.read(4)
                        zd_ncdomain_mod['Prec'][k - 1][j - 1][i - 1] = struct.unpack("f", context)[0]
            binFile.close()
        else:
            print('5870 readncbinary ', ffile, 'do not exist')
            exit()
        if ncbinarydrv.month == -1:  # 2
            zd_ncdomain_mod['intowinsizelat'] = 4  # the window size from observation to simulation in lat
            # lat中从观察到模拟的窗口大小
            zd_ncdomain_mod['intowinsizelon'] = 4  # the window size from observation to simulation in lon
            # lon中从观察到模拟的窗口大小
            Long1 = [[0.0] * zd_ncdomain_mod['countlon'] for i in range(1)]
            Lati1 = [[0.0] * 1 for i in range(zd_ncdomain_mod['countlat'])]
            # print('ncbinarydomain_module 1102', zd_ncdomain_mod['countlon'])
            # print('ncbinarydomain_module 1103', Long)
            for i in range(0, zd_ncdomain_mod['countlat']):
                Lati1[i][0] = Lati[i][0]
            for j in range(0, zd_ncdomain_mod['countlon']):
                Long1[0][j] = Long[0][j]

            if rminvalid:
                (Lati1, Long1) = rminvallatlon(Long1, minlon, temp[5], zd_ncdomain_mod['countlon'], Lati1, minlat,
                                               temp[4], zd_ncdomain_mod['countlat'])

            if zd_ncdomain_mod['intowinsizelat'] > lis.d.lnr:
                zd_ncdomain_mod['intowinsizelat'] = int(lis.d.lnr)
            if zd_ncdomain_mod['intowinsizelon'] > lis.d.lnc:
                zd_ncdomain_mod['intowinsizelon'] = int(lis.d.lnc)

            zd_ncdomain_mod['pcol'] = [[[1] * int(lis.d.lnc) for i in range(int(lis.d.lnr))] for j in
                                       range(zd_ncdomain_mod['intowinsizelon'] + 2)]
            zd_ncdomain_mod['prow'] = [[[1] * int(lis.d.lnc) for i in range(int(lis.d.lnr))] for j in
                                       range(zd_ncdomain_mod['intowinsizelat'] + 2)]
            reslat = lis.d.gridDesc[8]
            reslon = lis.d.gridDesc[9]
            for c in range(0, int(lis.d.lnc)):
                for r in range(0, int(lis.d.lnr)):
                    latcurr = minlat + (r - 1) * reslat
                    loncurr = minlon + (c - 1) * reslon

                    startintolat = int(r / temp[5] * reslat) + 1 - zd_ncdomain_mod['intowinsizelat'] / 2
                    endintolat = int(r / temp[5] * reslat) + 2 + zd_ncdomain_mod['intowinsizelat'] / 2
                    if startintolat < 0:
                        startintolat = 0
                        endintolat = startintolat + zd_ncdomain_mod['intowinsizelat'] - 1
                    if endintolat >= zd_ncdomain_mod['countlat']:
                        endintolat = zd_ncdomain_mod['countlat'] - 1
                    # startintolat=endintolat-zd_ncdomain_mod['intowinsizelat']+1

                    startintolon = int(c / temp[4] * reslon) + 1 - zd_ncdomain_mod['intowinsizelon'] / 2
                    endintolon = int(c / temp[4] * reslon) + 2 + zd_ncdomain_mod['intowinsizelon'] / 2
                    if startintolon < 0:
                        startintolon = 0
                        endintolon = startintolon + zd_ncdomain_mod['intowinsizelon'] - 1
                    if endintolon >= zd_ncdomain_mod['countlon']:
                        endintolon = zd_ncdomain_mod['countlon'] - 1
                    #   startintolon = endintolon - zd_ncdomain_mod['intowinsizelon'] + 1

                    sumlagalat = 0.0
                    for k in range(int(startintolat), int(endintolat + 1)):
                        for j in range(int(startintolat), k):
                            zd_ncdomain_mod['prow'][int(k - startintolat)][r][c] = \
                                zd_ncdomain_mod['prow'][int(k - startintolat)][r][c] * (
                                        latcurr - Lati1[int(j - startintolat)][0]) / (
                                        Lati1[int(k - startintolat)][0] - Lati1[int(j - startintolat)][0])
                        for j in range(k + 1, int(endintolat + 1)):
                            zd_ncdomain_mod['prow'][int(k - startintolat)][r][c] = \
                                zd_ncdomain_mod['prow'][int(k - startintolat)][r][c] * (
                                        latcurr - Lati1[int(j - startintolat)][0]) / (
                                        Lati1[int(k - startintolat)][0] - Lati1[int(j - startintolat)][0])

                    w3 = 0.0
                    for k in range(int(startintolat), int(endintolat) + 1):
                        w3 = w3 + abs(zd_ncdomain_mod['prow'][int(k - startintolat)][r][c])
                    for k in range(int(startintolat), int(endintolat) + 1):
                        zd_ncdomain_mod['prow'][int(k - startintolat)][r][c] = abs((zd_ncdomain_mod['prow'][
                            int(k - startintolat)][r][c])) / w3

                    for k in range(int(startintolon), int(endintolon) + 1):
                        for j in range(int(startintolon), k):
                            zd_ncdomain_mod['pcol'][int(k - startintolon)][r][c] = \
                                zd_ncdomain_mod['pcol'][int(k - startintolon)][r][c] * (
                                        loncurr - Long1[0][int(j - startintolon)]) / (
                                        Long1[0][int(k - startintolon)] - Long1[0][int(j - startintolon)])
                        for j in range(k + 1, int(endintolon + 1)):
                            zd_ncdomain_mod['pcol'][int(k - startintolon)][r][c] = \
                                zd_ncdomain_mod['pcol'][int(k - startintolon)][r][c] * (
                                        loncurr - Long1[0][int(j - startintolon)]) / (
                                        Long1[0][int(k - startintolon)] - Long1[0][int(j - startintolon)])

                    w2 = 0.0
                    for k in range(int(startintolon), int(endintolon) + 1):
                        w2 = w2 + abs(zd_ncdomain_mod['pcol'][int(k - startintolon)][r][c])
                    for k in range(int(startintolon), int(endintolon) + 1):
                        zd_ncdomain_mod['pcol'][int(k - startintolon)][r][c] = \
                            abs(zd_ncdomain_mod['pcol'][int(k - startintolon)][r][c]) / w2

            Lati = []
            Long = []
            Long1 = []
            Lati1 = []
        ncbinarydrv.month = ncmo
        print('5964 readcyclencbinary is done')
    # endif   !if (ncdrv % month.ne.mo) then
    localforcing = [[[0.0] * lis.f.nf for i in range(zd_ncdomain_mod['countlat'])] for i in
                    range(zd_ncdomain_mod['countlon'])]
    if lis.t.tstype == 0 or lis.t.tstype == 1 or lis.t.tstype == 2:
        if lis.t.tstype == 0:
            index = 0
            index = getncbinary_timeindex(yr, mo, da, hr, index)
            k = index
            c = k + 0
        if lis.t.tstype == 1:
            k = int((da - 1) * 24 / 3) + 1
            c = k + 7
        if lis.t.tstype == 2:
            k = 1
            c = zd_ncdomain_mod['time_step']

        for r in range(k - 1, c):
            for i in range(0, zd_ncdomain_mod['countlon']):
                for j in range(0, zd_ncdomain_mod['countlat']):
                    localforcing[i][j][0] = localforcing[i][j][0] + zd_ncdomain_mod['Tair'][r][j][i]
                    localforcing[i][j][1] = localforcing[i][j][1] + zd_ncdomain_mod['Qair'][r][j][i]
                    localforcing[i][j][2] = localforcing[i][j][2] + zd_ncdomain_mod['SWdown'][r][j][i]
                    localforcing[i][j][3] = localforcing[i][j][3] + zd_ncdomain_mod['LWdown'][r][j][i]
                    localforcing[i][j][4] = localforcing[i][j][4] + zd_ncdomain_mod['Wind'][r][j][i]
                    localforcing[i][j][5] = localforcing[i][j][5] + zd_ncdomain_mod['PSurf'][r][j][i]
                    localforcing[i][j][6] = localforcing[i][j][6] + zd_ncdomain_mod['Prec'][r][j][i]
        for i in range(0, zd_ncdomain_mod['countlon']):
            for j in range(0, zd_ncdomain_mod['countlat']):
                for x in range(0, lis.f.nf):
                    localforcing[i][j][x] = localforcing[i][j][x] / (c - k + 1)
    # else:
    # if (lis.t.tstype == 3):
    #
    if rminvalid:
        localforcing = removeinvalid(localforcing, zd_ncdomain_mod['countlon'], zd_ncdomain_mod['countlat'])
    reslat = lis.d.gridDesc[8]
    reslon = lis.d.gridDesc[9]
    for k in range(0, int(lis.f.nf)):
        for r in range(0, int(lis.d.lnr)):
            for c in range(0, int(lis.d.lnc)):
                if zd_lisdrv['gindex'][c][r] != -1:
                    w1 = 0.0
                    startintolat = int(r / ncbinarydrv.latres * reslat) + 1 - zd_ncdomain_mod['intowinsizelat'] / 2
                    endintolat = int(r / ncbinarydrv.latres * reslat) + 2 + zd_ncdomain_mod['intowinsizelat'] / 2
                    if startintolat < 0:
                        startintolat = 0
                        endintolat = startintolat + zd_ncdomain_mod['intowinsizelat'] - 1
                    if endintolat >= zd_ncdomain_mod['countlat']:
                        endintolat = zd_ncdomain_mod['countlat'] - 1
                        # startintolat=endintolat-zd_ncdomain_mod['intowinsizelat']+1

                    startintolon = int(c / ncbinarydrv.lonres * reslon) + 1 - zd_ncdomain_mod['intowinsizelon'] / 2
                    endintolon = int(c / ncbinarydrv.lonres * reslon) + 2 + zd_ncdomain_mod['intowinsizelon'] / 2
                    if startintolon < 0:
                        startintolon = 0
                        endintolon = startintolon + zd_ncdomain_mod['intowinsizelon'] - 1
                    if endintolon >= zd_ncdomain_mod['countlon']:
                        endintolon = zd_ncdomain_mod['countlon'] - 1
                    for j in range(int(endintolat), int(startintolat) - 1, -1):  # endintolat, startintolat, -1
                        for i in range(int(startintolon), int(endintolon) + 1):  # endintolon
                            w1 = w1 + zd_ncdomain_mod['prow'][j - int(startintolat)][r][c] * \
                                 localforcing[j - int(startintolat)][i - int(startintolon)][k] * \
                                 zd_ncdomain_mod['pcol'][i - int(startintolon)][r][c]
                    kk = zd_lisdrv['gindex'][c][r]
                    # print('4031', type(kk))
                    if order == 1:
                        # print('4998', np.array(zd_basefor_mod['glbdata1']).shape, kk, k, w1)
                        zd_basefor_mod['glbdata1'][kk][k] = w1
                    else:
                        zd_basefor_mod['glbdata2'][kk][k] = w1

    del localforcing
    if ncbinarydrv.month == -1:
        deallocate(Long1)
        deallocate(Lati1)

    print('6033 readncbinary')


# time_interp_ncbinary
def time_interp_ncbinary():
    # ==== Local Variables========局部变量=========
    ier = 0
    c, f, zdoy = 0, 0, 0
    bdoy, byr, bmo = 0, 0, 0
    bda, bhr, bmn = 0, 0, 0
    btime = 0.0
    wt1, wt2, czb, cze, czm, gmt1, gmt2 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    zw1, zw2 = 0.0, 0.0
    if zd_spmd['masterproc']:
        if get_nstep(lis.t) == 0:
            lis.f.nforce = lis.f.nf
        else:
            lis.f.nforce = lis.f.nf
    if lis.t.tstype == 0:  # simulate at a month resolution 以月分辨率模拟
        btime = ncbinarydrv.nctime1
        t2dlist = [btime, bdoy, gmt1, byr, bmo, bda, bhr, bmn]
        time2date(t2dlist)
        btime, bdoy, gmt1, byr, bmo, bda, bhr, bmn = t2dlist[0], t2dlist[1], t2dlist[2], t2dlist[3], t2dlist[4], \
                                                     t2dlist[5], t2dlist[6], t2dlist[7]
        btime = ncbinarydrv.nctime2
        t2dlist2 = [btime, bdoy, gmt2, byr, bmo, bda, bhr, bmn]
        time2date(t2dlist2)
        btime, bdoy, gmt2, byr, bmo, bda, bhr, bmn = t2dlist2[0], t2dlist2[1], t2dlist2[2], t2dlist2[3], t2dlist2[4], \
                                                     t2dlist2[5], t2dlist2[6], t2dlist2[7]
        wt1 = ((ncbinarydrv.nctime2 - lis.t.time) / (ncbinarydrv.nctime2 - ncbinarydrv.nctime1))
        wt2 = 1.0 - wt1
        for f in range(0, lis.f.nforce):
            for c in range(0, lis.d.ngrid):  # gdi(iam)
                zd_lisdrv['grid'][c].forcing[f] = zd_basefor_mod['glbdata1'][c][f] * wt1 + \
                                                  zd_basefor_mod['glbdata2'][c][f] * wt2
    if lis.t.tstype != 0:  # simulate at a other resolution 以其他分辨率模拟
        for f in range(0, lis.f.nforce):
            for c in range(0, lis.d.ngrid):
                zd_lisdrv['grid'][c].forcing[f] = zd_basefor_mod['glbdata1'][c][f]


# Obsfile_module
def inquirefile(existence):  # existence=[existence]列表包装
    # import lisdrv_module
    # from lisdrv_module import lis
    name = ['' for i in range(lis.a.ndva)]
    levles_filename(name)
    for i in range(0, lis.a.ndva):
        if os.path.exists(name[i]):
            existence[0] = True
        else:
            existence[0] = False


def read_evap():
    name = ['' for i in range(lis.a.ndva)]
    t = 0
    tt = 0
    i = 0
    io = 0
    c = 0
    r = 0
    cindex = 0
    count = 0
    obs = [[[]]]
    obslis = [[[]]]
    obsr, obsc = 0.0, 0.0
    obscol = 0
    obsrow = 0
    runsrow = 0
    runscol = 0
    runerow = 0
    runecol = 0
    win = 0
    lt = t
    evapobs = [[0] * lis.a.win for i in range(lis.d.nch)]  # 改变字典中该type类型数组的大小。
    # the column and row of the observation data 观察数据的列和行
    obscol = math.floor((lis.a.obelo - lis.a.obslo) / lis.a.oblor) + 1
    obsrow = math.floor((lis.a.obela - lis.a.obsla) / lis.a.oblar) + 1

    # the start row and column of simulation in the observation 观察中模拟的起始行和起始列
    runsrow = math.floor((lis.d.gridDesc[3] - lis.a.obsla) / lis.a.oblar) + 1
    runscol = math.floor((lis.d.gridDesc[4] - lis.a.obslo) / lis.a.oblor) + 1

    # the end row and column of simulation in the observation 观察中模拟的结束行和结束列观察中模拟的结束行和结

    runerow = math.floor((lis.d.gridDesc[6] - lis.a.obsla) / lis.a.oblar) + 1
    runecol = math.floor((lis.d.gridDesc[7] - lis.a.obslo) / lis.a.oblor) + 1

    # allocate temporal memory 分配时间内存
    obs = [[[0.0] * lis.a.ndav for i in range(runecol - runscol + 1)] for i in range(runerow - runsrow + 1)]
    obslis = [[[0.0] * int(lis.a.ndav) for i in range(int(lis.d.lnc))] for i in range(int(lis.d.lnr))]
    obs = 0
    obslis = 0
    print('6127 obsfile_module', runecol, runscol, runerow, runsrow)
    win = 1
    lis.t = lis.tt
    while win <= lis.a.win:  # win从1开始还是从0开始，需商榷
        levles_filename(name)  # 生成观测文件名
    for i in range(0, lis.a.ndva):
        f = open(name[i], 'r')
        content = f.readlines()
        for r in range(runsrow, runerow + 1):  # do r=1,lis.d.gnr
            for c in range(runscol, runecol + 1):  # do c=1, lis.d.gnc
                t = 2 * obscol * obsrow + (r - runsrow) * obscol + c - runscol + 1
                obs[r - runsrow + 1 - 1][c - runscol + 1 - 1][i] = content[t - 1]
        f.close()
        print('6140 Reading retrieval data ', name[i].rstrip())
        for r in range(0, lis.d.lnr):
            for c in range(0, lis.d.lnc):
                obsc = math.floor(c * lis.d.gridDesc[9] / lis.a.oblor)
                obsr = math.floor(r * lis.d.gridDesc[8] / lis.a.oblar)
                obslis[c][r][i] = obs[obsr][obsc][i]
        for r in range(0, lis.d.lnr):
            for c in range(0, lis.d.lnc):
                if zd_lisdrv['gindex'][c][r] != -1:
                    count = zd_spmd['c_str'][zd_spmd['iam'] + 1 - 1] + (c - 1) + (
                            zd_spmd['r_str'][zd_spmd['iam'] + 1 - 1] + (lis.d.lnr - r) - 1) * lis.d.gnc
                    evapobs[zd_lisdrv['gindex'][r][c]][win - 1].evap[i] = obslis[r][c][i]
        for t in range(0, lis.d.nch):
            if evapobs[t][win - 1].evap[i] != -999:
                if evapobs[t][win - 1].evap[i] < 0.0:
                    evapobs[t][win - 1].evap[i] = 0.0
                if evapobs[t][win - 1].evap[i] > 700:
                    evapobs[t][win - 1].evap[i] = 700
            evapobs[t][win - 1].mcv = pow(zd_evapob['oer'], 2)  # observation error covariance
    advance_timestep(lis.t)
    win = win + 1


def levles_filename(name):
    # import lisdrv_module
    # import evapobsdrv_module
    # from lisdrv_module import lis
    # from evapobsdrv_module import evapobsdir
    c = 0
    ndir = zd_evapob['evapobsdir']
    zfc = ndir.rjust(70, ' ')
    restzfc = ['' for i in range(80 - len(zfc))]
    fbase = list(zfc)
    fbase.extend(restzfc)
    for i in range(0, 80):
        if fbase[i] == ' ' and c == 0:
            c = i - 1
    ftimey = list(format(lis.t.yr, '>4.0f'))
    for i in range(0, 4):
        if ftimey[i] == ' ':
            ftimey[i] = '0'
    zfc2 = format(lis.t.yr, '>4.0f') + format(lis.t.mo, '>2.0f') + format(lis.t.da, '>2.0f')
    ftime = list(zfc2)
    for i in range(0, 8):
        if ftime[i] == ' ':
            ftime[i] = '0'
    # get the file name
    zfc3 = format(lis.t.yr, '>4.0f') + format(lis.t.mo, '>2.0f') + format(lis.t.da, '>2.0f') + format(lis.t.hr, '>2.0f')
    ftimeb = list(zfc3)
    for i in range(0, 10):
        if ftimeb[i] == ' ':
            ftimeb[i] = '0'
    # get all the observation files
    for k in range(1, lis.a.ndva + 1):
        if k == 1:
            zfc4 = ''.join(fbase[0:c]) + '/' + ''.join(ftimey[0:4]) + '/' + ''.join(ftime[0:8]) + '/' + ''.join(
                ftimeb[0:10]) + '.gs4r'
            name[k - 1] = zfc4[0:80]  # 如果zfc4的长度小于80，这里默认拿出全部，格式是正确的。
        if k == 2:
            zfc5 = ''.join(fbase[0:c]) + '/' + ''.join(ftime[0:4]) + '/' + ''.join(ftime[0:8]) + '/' + ''.join(
                ftimeb[0:10]) + '.les'
            name[k - 1] = zfc5[0:80]  # 如果zfc4的长度小于80，这里默认拿出全部，格式是正确的。
        if k == 3:
            zfc6 = ''.join(fbase[0:c]) + '/' + ''.join(ftime[0:4]) + '/' + ''.join(ftime[0:8]) + '/' + ''.join(
                ftimeb[0:10]) + '.lev'
            name[k - 1] = zfc6[0:80]  # 如果zfc4的长度小于80，这里默认拿出全部，格式是正确的。


# lsm_module
def LIS_lsm_init():
    lis_set_indices()


def LIS_setuplsm():
    noah_setup()


def LIS_readrestart():
    noahrst()


def LIS_force2tile():  # 其中noah_f2t函数未找到。
    # import lisdrv_module
    # from lisdrv_module import lis,grid,tile
    for t in range(0, lis.d.nch):
        index = zd_lisdrv['tile'][t].index
        # lsmf2t(lis.d.lsm, t, grid(index).forcing)
        noah_f2t(t, zd_lisdrv['grid'][index].forcing)


def noah_f2t(t, forcing):
    if lis.f.force == 8:
        zd_varder['noah'][t].forcing[0] = forcing[2]
        zd_varder['noah'][t].forcing[1] = forcing[3]
        zd_varder['noah'][t].forcing[2] = forcing[5]
        zd_varder['noah'][t].forcing[3] = forcing[6]
        zd_varder['noah'][t].forcing[4] = forcing[0]
        zd_varder['noah'][t].forcing[5] = forcing[0]
        zd_varder['noah'][t].forcing[6] = forcing[4]
        zd_varder['noah'][t].forcing[7] = forcing[1]
        zd_varder['noah'][t].forcing[8] = 0
    else:
        for f in range(0, lis.f.nforce):
            zd_varder['noah'][t].forcing[f] = forcing[f]
    return


def LIS_lsm_main():
    existence = True
    if lis.a.daalg == 2:
        evapobs_setup()
        existence = 0
        existence = [existence]
        inquirefile(existence)
        existence = existence[0]  # read obs card
        if existence:
            read_evap()
            evap_enkf_init()
        else:
            noah_main()
    else:
        noah_main()


def LIS_lsm_output():
    # import lisdrv_module
    # from lisdrv_module import lis
    noah_output()


def LIS_writerestart():
    noah_writerst()


def LIS_writerange():
    fbase = lis.o.odir
    noahdir1 = '/EXP' + str(lis.o.expcode) + '/'
    expcode1 = 'Range.txt'
    temp = fbase + noahdir1 + expcode1

    with open(temp, 'w') as f:
        f.write(str(lis.d.gridDesc[4]))
        f.write('\n' + str(lis.d.gridDesc[7]))
        f.write('\n' + str(lis.d.gridDesc[3]))
        f.write('\n' + str(lis.d.gridDesc[6]))
        f.write('\n' + str(lis.d.gridDesc[9]))
        f.write('\n' + str(lis.d.gridDesc[8]))
        if (lis.d.gridDesc[7] - lis.d.gridDesc[4]) % lis.d.gridDesc[9] == 0:
            colnumber = int((lis.d.gridDesc[7] - lis.d.gridDesc[4]) / lis.d.gridDesc[9])
        else:
            colnumber = int((lis.d.gridDesc[7] - lis.d.gridDesc[4]) / lis.d.gridDesc[9]) + 1
        if (lis.d.gridDesc[6] - lis.d.gridDesc[3]) % lis.d.gridDesc[8] == 0:
            rownumber = int((lis.d.gridDesc[6] - lis.d.gridDesc[3]) / lis.d.gridDesc[8])
        else:
            rownumber = int((lis.d.gridDesc[6] - lis.d.gridDesc[3]) / lis.d.gridDesc[8]) + 1
        f.write('\n' + str(colnumber))
        f.write('\n' + str(rownumber))
    f.close()


# noah_physics
# #INTERFACE:
def SFLX(ntl,
         ICE, DT, ZLVL, NSOIL, SLDPTH,
         LWDN, SOLDN, SFCPRS, PRCP, SFCTMP, Q2, SFCSPD,
         TH2, Q2SAT, DQSDT2,
         SLOPE, SHDFAC, SHDMIN, PTU, ALB, SNOALB,
         RSMIN, RGL, HS, SNUP, Z0, XLAI, NROOT,
         PSISAT, BEXP, DKSAT, SMCMAX, QUARTZ, DWSAT,
         SMCWLT, SMCREF, SMCDRY, F1, KDT, FRZX, FRZFACT, TBOT,
         CMC, T1, STC, SMC, SH2O, SNOWH, SNEQV, ALBEDO, CH, CM,
         EVP, ETA, SHEAT,
         EC, EDIR, ET, ETT, ESNOW, DRIP, DEW,
         BETA, ETP, SSOIL,
         FLX1, FLX2, FLX3,
         SNOMLT, SNCOVR,
         RUNOFF1, RUNOFF2, RUNOFF3,
         RC, PC, RCS, RCT, RCQ, RCSOIL,
         SOILW, SOILT, SOILM):
    # ----------------------------------------------------------------------
    # SUBROUTINE SFLX - VERSION 2.7 - June 2nd 2003 子程序SFLX-版本2.7-2003年6月2日
    # ----------------------------------------------------------------------
    # SUB-DRIVER FOR "NOAH/OSU LSM" FAMILY OF PHYSICS SUBROUTINES FOR A “NOAH/OSU LSM”系列物理子程序的子驱动程序
    # SOIL/VEG/SNOWPACK LAND-SURFACE MODEL TO UPDATE SOIL MOISTURE, SOIL 土壤/植被/积雪地表模型，用于更新土壤水分、土壤
    # ICE, SOIL TEMPERATURE, SKIN TEMPERATURE, SNOWPACK WATER CONTENT, 冰、土壤温度、皮肤温度、积雪含水量
    # SNOWDEPTH, AND ALL TERMS OF THE SURFACE ENERGY BALANCE AND SURFACE 雪深，表面能量平衡和表面
    # WATER BALANCE (EXCLUDING INPUT ATMOSPHERIC forcingS OF DOWNWARD 水平衡（不包括向下的输入大气作用力）
    # RADIATION AND PRECIP) 辐射和精度
    # ----------------------------------------------------------------------
    # SFLX ARGUMENT LIST KEY: SFLX参数列表键：
    # ----------------------------------------------------------------------
    #  C  CONFIGURATION INFORMATION 配置信息
    #  F  forcing DATA 强制数据
    #  I  OTHER (INPUT) forcing DATA 其他（输入）强制数据
    #  S  SURFACE CHARACTERISTICS 表面特征
    #  H  HISTORY (STATE) VARIABLES 历史（状态）变量
    #  O  OUTPUT VARIABLES 输出变量
    #  D  DIAGNOSTIC OUTPUT 诊断输出
    # ----------------------------------------------------------------------
    # 1. CONFIGURATION INFORMATION (C):  配置信息（C）：
    # ----------------------------------------------------------------------
    #   FFROZP     FRACTION OF FROZEN PRECIPITATION  冻结降水的FFROZP分数
    #   ICE	       SEA-ICE FLAG  (=1: SEA-ICE, =0: LAND) 海冰旗（=1：海冰，=0：陆地）
    #   DT	       TIMESTEP (SEC) (DT SHOULD NOT EXCEED 3600 SECS, RECOMMEND DT时间步长（秒）（建议DT不超过3600秒
    #                1800 SECS OR LESS) 1800秒或更少）
    #   ZLVL       HEIGHT (M) ABOVE GROUND OF ATMOSPHERIC forcing VARIABLES 大气强迫变量的ZLVL地面高度（M）
    #   NSOIL      NUMBER OF SOIL LAYERS (AT LEAST 2, AND NOT GREATER THAN  NSOIL土层数（至少2层，且不大于
    #                PARAMETER NSOLD SET BELOW)  参数NSOLD设置如下）
    #   SLDPTH     THE THICKNESS OF EACH SOIL LAYER (M) SLDPTH各土层的厚度（M）
    # ----------------------------------------------------------------------
    # 2. forcing DATA (F): 2、强制数据（F）：
    # ----------------------------------------------------------------------
    #   LWDN       LW DOWNWARD RADIATION (W M-2; POSITIVE, NOT NET LONGWAVE) LWDN LW向下辐射（W M-2；正，非净长波）
    #   SOLDN      SOLAR DOWNWARD RADIATION (W M-2; POSITIVE, NOT NET SOLAR) SOLDN太阳向下辐射（W M-2；正，非净太阳）
    #   SFCPRS     PRESSURE AT HEIGHT ZLVL ABOVE GROUND (PASCALS) 地面以上ZLVL高度的SFCPRS压力（帕斯卡）
    #   PRCP       PRECIP RATE (KG M-2 S-1) (NOTE, THIS IS A RATE) PRCP精确速率（KG M-2 S-1）（注意，这是一个速率）
    #   SFCTMP     AIR TEMPERATURE (K) AT HEIGHT ZLVL ABOVE GROUND 地面以上ZLVL高度的SFCTMP空气温度（K）
    #   TH2        AIR POTENTIAL TEMPERATURE (K) AT HEIGHT ZLVL ABOVE GROUND 地面以上ZLVL高度处的TH2空气位温度（K）
    #   Q2         MIXING RATIO AT HEIGHT ZLVL ABOVE GROUND (KG KG-1) 地面以上ZLVL高度处的Q2混合比（KG KG-1）
    # ----------------------------------------------------------------------
    # 3. OTHER forcing (INPUT) DATA (I): 其他强制（输入）数据（I）：
    # ----------------------------------------------------------------------
    #   SFCSPD     WIND SPEED (M S-1) AT HEIGHT ZLVL ABOVE GROUND 地面以上ZLVL高度处的风速（M S-1）
    #   Q2SAT      SAT MIXING RATIO AT HEIGHT ZLVL ABOVE GROUND (KG KG-1) 地面以上ZLVL高度的SAT混合比（KG-1）
    #   DQSDT2     SLOPE OF SAT SPECIFIC HUMIDITY CURVE AT T=SFCTMP T=SFCTMP时SAT比湿度曲线的斜率
    #                (KG KG-1 K-1)
    # ----------------------------------------------------------------------
    # 4. CANOPY/SOIL CHARACTERISTICS (S): 树冠/土壤特征：
    # ----------------------------------------------------------------------
    #   VEGTYP     VEGETATION TYPE (INTEGER INDEX) 植被类型（整数指数）
    #   SOILTYP    SOIL TYPE (INTEGER INDEX) 土壤类型（整数指数）
    #   SLOPETYP   CLASS OF SFC SLOPE (INTEGER INDEX) SFC斜率等级（整数指数）
    #   SHDFAC     AREAL FRACTIONAL COVERAGE OF GREEN VEGETATION 绿色植被的面积分数覆盖率
    #                (FRACTION= 0.0-1.0)
    #   SHDMIN     MINIMUM AREAL FRACTIONAL COVERAGE OF GREEN VEGETATION 绿色植被的最小面积分数覆盖率
    #                (FRACTION= 0.0-1.0) <= SHDFAC
    #   PTU        PHOTO THERMAL UNIT (PLANT PHENOLOGY FOR ANNUALS/CROPS) 光热单位（一年生/作物的植物物候）
    #                (NOT YET USED, BUT PASSED TO REDPRM FOR FUTURE USE IN
    #                VEG PARMS)
    #   ALB        BACKROUND SNOW-FREE SURFACE ALBEDO (FRACTION), FOR JULIAN  JULIAN的无雪地表反照率（分数）
    #                DAY OF YEAR (USUALLY FROM TEMPORAL INTERPOLATION OF 一年中的某一天（通常从
    #                MONTHLY MEAN VALUES' CALLING PROG MAY OR MAY NOT 月平均值的调用程序可能会也可能不会
    #                INCLUDE DIURNAL SUN ANGLE EFFECT) 包括日间太阳角度效应）
    #   SNOALB     UPPER BOUND ON MAXIMUM ALBEDO OVER DEEP SNOW (E.G. FROM 深雪（例如从
    #                ROBINSON AND KUKLA, 1985, J. CLIM. & APPL. METEOR.) ROBINSON和KUKLA，1985年，J.CLIM应用。流星。）
    #   TBOT       BOTTOM SOIL TEMPERATURE (LOCAL YEARLY-MEAN SFC AIR TBOT底部土壤温度（当地年平均SFC空气
    #                TEMPERATURE) 温度）
    # ----------------------------------------------------------------------
    # 5. HISTORY (STATE) VARIABLES (H): 历史（状态）变量（H）：
    # ----------------------------------------------------------------------
    #  CMC         CANOPY MOISTURE CONTENT (M)  冠层含水量（M）
    #  T1          GROUND/CANOPY/SNOWPACK) EFFECTIVE SKIN TEMPERATURE (K) T1地面/树冠/积雪）有效皮肤温度（K）
    #  STC(NSOIL)  SOIL TEMP (K) STC（NSOIL）土壤温度（K）
    #  SMC(NSOIL)  TOTAL SOIL MOISTURE CONTENT (VOLUMETRIC FRACTION) SMC（NSOIL）土壤总含水量（体积分数）
    #  SH2O(NSOIL) UNFROZEN SOIL MOISTURE CONTENT (VOLUMETRIC FRACTION) SH2O（NSOIL）未冻土壤含水量（体积分数）
    #                NOTE: FROZEN SOIL MOISTURE = SMC - SH2O 注：冻土水分=SMC-SH2O
    #  SNOWH       ACTUAL SNOW DEPTH (M) 雪实际雪深（M）
    #  SNEQV       LIQUID WATER-EQUIVALENT SNOW DEPTH (M) SNEQV液态水等效雪深（M）
    #                NOTE: SNOW DENSITY = SNEQV/SNOWH 注：雪密度=SNEQV/SNOWH
    #  ALBEDO      SURFACE ALBEDO INCLUDING SNOW EFFECT (UNITLESS FRACTION)  包括雪效应的反照率表面反照率（无单位分数）
    #                =SNOW-FREE ALBEDO (ALB) WHEN SNEQV=0, OR   =SNEQV=0时的无雪反照率（ALB），或
    #                =FCT(MSNOALB,ALB,VEGTYP,SHDFAC,SHDMIN) WHEN   SNEQV>0 =当SNEQV>0时的FCT（MSNOALB、ALB、VEGTYP、
    #  CH          SURFACE EXCHANGE COEFFICIENT FOR HEAT AND MOISTURE SHDFAC、SHDMIN） CH热湿表面交换系数
    #                (M S-1); NOTE: CH IS TECHNICALLY A CONDUCTANCE SINCE （M S-1）；注：CH在技术上是电导，因为
    #                IT HAS BEEN MULTIPLIED BY WIND SPEED. 它已乘以风速。
    #  CM          SURFACE EXCHANGE COEFFICIENT FOR MOMENTUM (M S-1); NOTE: CM动量表面交换系数（M S-1）；注：
    #                CM IS TECHNICALLY A CONDUCTANCE SINCE IT HAS BEEN CM在技术上是一种电导，因为它已经
    #                MULTIPLIED BY WIND SPEED.  CM IS NOT NEEDED IN SFLX 乘以风速。SFLX中不需要CM
    # ----------------------------------------------------------------------
    # 6. OUTPUT (O): 6、输出（O）：
    # ----------------------------------------------------------------------
    # OUTPUT VARIABLES NECESSARY FOR A COUPLED NUMERICAL WEATHER PREDICTION   耦合数值天气预报所需的输出变量
    # MODEL, E.G. NOAA/NWS/NCEP MESOSCALE ETA MODEL.  FOR THIS APPLICATION, #模式，例如NOAA/NWS/NCEP中尺度ETA模式。对于此应用程序，
    # THE REMAINING OUTPUT/DIAGNOSTIC/PARAMETER BLOCKS BELOW ARE NOT 以下剩余的输出/诊断/参数块不是
    # NECESSARY.  OTHER APPLICATIONS MAY REQUIRE DIFFERENT OUTPUT VARIABLES. 必要的。其他应用程序可能需要不同的输出变量。
    #   ETA        ACTUAL LATENT HEAT FLUX (W M-2: NEGATIVE, IF UP FROM SURFACE) ETA实际潜热通量（W M-2：负，如果从表面向上）
    #   SHEAT      SENSIBLE HEAT FLUX (W M-2: NEGATIVE, IF UPWARD FROM SHEAT显热通量（W M-2：负，如果从
    # ----------------------------------------------------------------------
    #   EC         CANOPY WATER EVAPORATION (W M-2)  EC冠层水分蒸发（W M-2）
    #   EDIR       DIRECT SOIL EVAPORATION (W M-2) EDIR直接土壤蒸发（W M-2）
    #   ET(NSOIL)  PLANT TRANSPIRATION FROM A PARTICULAR ROOT (SOIL) LAYER ET（NSOIL）植物从特定根（土壤）层的蒸腾作用
    #                 (W M-2)  （W M-2）
    #   ETT        TOTAL PLANT TRANSPIRATION (W M-2) ETT植物蒸腾总量（W M-2）
    #   ESNOW      SUBLIMATION FROM SNOWPACK (W M-2)  来自积雪的ESNOW升华（W M-2）
    #   DRIP       THROUGH-FALL OF PRECIP AND/OR DEW IN EXCESS OF CANOPY 雨水和/或露水滴落超过树冠
    #                WATER-HOLDING CAPACITY (M) 持水量（M）
    #   DEW        DEWFALL (OR FROSTFALL FOR T<273.15) (M) 露水降（或T<273.15时的霜降）（M）
    # ----------------------------------------------------------------------
    #   BETA       RATIO OF ACTUAL/POTENTIAL EVAP (DIMENSIONLESS) 实际/潜在蒸发排放比（无量纲）
    #   ETP        POTENTIAL EVAPORATION (W M-2) ETP潜在蒸发（W M-2）
    #   SSOIL      SOIL HEAT FLUX (W M-2: NEGATIVE IF DOWNWARD FROM SURFACE) SSOIL土壤热通量（W M-2：如果从表面向下，则为负值）
    # ----------------------------------------------------------------------
    #   FLX1       PRECIP-SNOW SFC (W M-2)   PRECIP-SNOW SFC（W M-2）
    #   FLX2       FREEZING RAIN LATENT HEAT FLUX (W M-2)  FLX2冻雨潜热通量（W M-2）
    #   FLX3       PHASE-CHANGE HEAT FLUX FROM SNOWMELT (W M-2)  来自融雪的FLX3相变热通量（W M-2）
    # ----------------------------------------------------------------------
    #   SNOMLT     SNOW MELT (M) (WATER EQUIVALENT)  融雪（M）（水当量）
    #   SNCOVR     FRACTIONAL SNOW COVER (UNITLESS FRACTION, 0-1)  SNCOVR分数积雪（无单位分数，0-1）
    # ----------------------------------------------------------------------
    #   RUNOFF1    SURFACE RUNOFF (M S-1), NOT INFILTRATING THE SURFACE  地表径流（M S-1），未渗入地表
    #   RUNOFF2    SUBSURFACE RUNOFF (M S-1), DRAINAGE OUT BOTTOM OF LAST 径流2地下径流（M S-1），最后一个底部排水
    #                SOIL LAYER 土层
    #   RUNOFF3    NUMERICAL TRUNCTATION IN EXCESS OF POROSITY (SMCMAX) 径流3数值扩展超过孔隙度（SMCMAX）
    #                FOR A GIVEN SOIL LAYER AT THE END OF A TIME STEP 对于时间步长结束时的给定土层
    # ----------------------------------------------------------------------
    #   RC         CANOPY RESISTANCE (S M-1)  冠层阻力（S M-1）
    #   PC         PLANT COEFFICIENT (UNITLESS FRACTION, 0-1) WHERE PC*ETP PC工厂系数（无单位分数，0-1），其中PC*ETP
    #                = ACTUAL TRANSPIRATION  =实际蒸腾
    #   XLAI       LEAF AREA INDEX (DIMENSIONLESS) XLAI叶面积指数（无量纲）
    #   RSMIN      MINIMUM CANOPY RESISTANCE (S M-1)  RSMIN最小冠层阻力（S M-1）
    #   RCS        INCOMING SOLAR RC FACTOR (DIMENSIONLESS) RCS入射太阳RC系数（无量纲）
    #   RCT        AIR TEMPERATURE RC FACTOR (DIMENSIONLESS) RCT空气温度RC系数（无量纲）
    #   RCQ        ATMOS VAPOR PRESSURE DEFICIT RC FACTOR (DIMENSIONLESS) RCQ大气蒸汽压差RC因子（无量纲）
    #   RCSOIL     SOIL MOISTURE RC FACTOR (DIMENSIONLESS) RC土壤水分RC系数（无量纲）
    # ----------------------------------------------------------------------
    # 7. DIAGNOSTIC OUTPUT (D): 7、诊断输出（D）：
    # ----------------------------------------------------------------------
    #   SOILW      AVAILABLE SOIL MOISTURE IN ROOT ZONE (UNITLESS FRACTION BETWEEN SMCWLT AND SMCMAX)
    # 根区土壤有效水分（SMCWLT和SMCMAX之间的无单位分数）
    #   SOILM      TOTAL SOIL COLUMN MOISTURE CONTENT (FROZEN+UNFROZEN) (M)
    # 土柱总含水量（冻结+未冻结）（M）
    # ----------------------------------------------------------------------
    # 8. PARAMETERS (P):  8、参数（P）：
    # ----------------------------------------------------------------------
    #   SMCWLT     WILTING POINT (VOLUMETRIC)   SMCWLT萎蔫点（体积）
    #   SMCDRY     DRY SOIL MOISTURE THRESHOLD WHERE DIRECT EVAP FRM TOP SMC干燥土壤水分阈值，其中直接蒸发FRM顶部
    #                LAYER ENDS (VOLUMETRIC) 层终点（体积）
    #   SMCREF     SOIL MOISTURE THRESHOLD WHERE TRANSPIRATION BEGINS TO SMCREF土壤水分阈值，其中蒸腾开始
    #                STRESS (VOLUMETRIC) 应力（体积）
    #   SMCMAX     POROSITY, I.E. SATURATED VALUE OF SOIL MOISTURE SMCMAX孔隙度，即土壤水分的饱和值
    #                (VOLUMETRIC) （体积）
    #   NROOT      NUMBER OF ROOT LAYERS, A FUNCTION OF VEG TYPE 根层的根数，是植物类型的函数
    # ----------------------------------------------------------------------
    # integer, PARAMETER:: NSOLD=20  整数，参数：：NSOLD=20
    # integer NTL  整数NTL

    # ----------------------------------------------------------------------
    # DECLARATIONS - PARAMETERS  声明-参数
    # ----------------------------------------------------------------------
    import numpy as np
    NSOLD = 20
    ZSOIL = np.zeros(NSOLD, dtype=float)
    RTDIS = np.zeros(NSOLD, dtype=float)
    TFREEZ = 273.15
    LVH2O = 2.501E+6
    LSUBS = 2.83E+6
    R = 287.04
    CP = 1004.5
    T24 = 0.0
    RCH, EPSCA, RR = 0.0, 0.0, 0.0
    TSNOW = 0.0

    CFACTR = 0.5  # CANOPY WATER PARAMETER  冠层水分参数
    CMCMAX = 0.5E-3  # CANOPY WATER PARAMETER 冠层水分参数
    RSMAX = 5000.0  # MAX. STOMATAL RESISTANCE 最大气孔阻力
    TOPT = 298.0  # OPTIMUM TRANSPIR AIR TEMP 最佳运输空气温度
    SBETA = -2.0  # TO CALC VEG EFFECT ON SHFLX 计算VEG对SHFLX的影响

    # real, PARAMETER:: FRZK=0.15
    # ICE content threshold in soil  土壤含冰量阈值
    # real, PARAMETER:: REFDK = 2.0E-6
    # Reference
    REFKDT = 3.0

    FXEXP = 2.0  # BARE SOIL EVAP EXP USED IN DEVAP   DEVAP中使用的裸土EVAP EXP
    CSOIL = 2.00E+6  # SOIL HEAT CAPACITY [J M-3 K-1]   土壤热容[J M-3 K-1]

    # -- SPECIFY DEPTH[M] OF LOWER BOUNDARY SOIL TEMPERATURE.   指定下限土壤温度的深度[M]
    #     PARAMETER(ZBOT = -3.0)
    ZBOT = -8.0

    # -- SPECIFY SNOW DISTRIBUTION SHAPE PARAMETER SALP - SHAPE PARAMETER   指定雪分布形状参数SALP-形状参数
    #    OF DISTRIBUTION FUNCTION OF SNOW COVER. (from ANDERSON, HYDRO-17)   积雪的分布函数。（安德森，HYDRO-17）
    #    BEST FIT IS WHEN SALP = 2.6
    #      PARAMETER(SALP = 2.6)
    # - Changed for version 2.6 June 2nd 2003 *   更改为2.6版2003年6月2日*
    SALP = 4.0

    # --  PARAMETER USED TO CALCULATE ROUGHNESS LENGTH OF HEAT.   用于计算热粗糙度长度的参数。
    #      PARAMETER(CZIL = 0.2)
    # - Changed for version 2.6 June 2nd 2003 *   更改为2.6版2003年6月2日*
    CZIL = 0.075

    # ----------------------------------------------------------------------
    #   INITIALIZATION   初始化
    # ----------------------------------------------------------------------
    RUNOFF1 = 0.0
    RUNOFF2 = 0.0
    RUNOFF3 = 0.0
    SNOMLT = 0.0
    DF1 = 0.0
    # ----------------------------------------------------------------------
    #  THE VARIABLE "ICE" IS A FLAG DENOTING SEA-ICE CASE   变量“ICE”是表示海冰情况的标志
    # ----------------------------------------------------------------------
    if ICE == 1:
        # ----------------------------------------------------------------------
        # SEA-ICE LAYERS ARE EQUAL THICKNESS AND SUM TO 3 METERS   海冰层厚度相等，总和为3米
        # ----------------------------------------------------------------------
        for KZ in range(0, NSOIL):
            ZSOIL[KZ] = -3. * float(KZ) / float(NSOIL - 1)
    else:
        # ----------------------------------------------------------------------
        # CALCULATE DEPTH (NEGATIVE) BELOW GROUND FROM TOP SKIN SFC TO BOTTOM OF  计算从顶部皮肤SFC到底部的地下深度（负）
        #   EACH SOIL LAYER.  NOTE:  SIGN OF ZSOIL IS NEGATIVE (DENOTING BELOW  各土层。注：ZSOIL的符号为负数（表示如下
        #   GROUND) 接地）
        # ----------------------------------------------------------------------
        ZSOIL[0] = -SLDPTH[0]
        for KZ in range(1, NSOIL):
            ZSOIL[KZ] = -SLDPTH[KZ] + ZSOIL[KZ - 1]
    # ----------------------------------------------------------------------
    # NEXT IS CRUCIAL CALL TO SET THE LAND-SURFACE PARAMETERS, INCLUDING   接下来是设置地表参数的关键调用，包括
    # SOIL-TYPE AND VEG-TYPE DEPENDENT PARAMETERS.  土壤类型和植被类型相关参数。
    # ----------------------------------------------------------------------
    #      CALL REDPRM (VEGTYP,SOILTYP,SLOPETYP,
    #     +      	   CFACTR,CMCMAX,RSMAX,TOPT,REFKDT,KDT,SBETA,
    #     O      	   SHDFAC,RSMIN,RGL,HS,ZBOT,FRZX,PSISAT,SLOPE,
    #     +      	   SNUP,SALP,BEXP,DKSAT,DWSAT,SMCMAX,SMCWLT,SMCREF,
    #     O      	   SMCDRY,F1,QUARTZ,FXEXP,RTDIS,SLDPTH,ZSOIL,
    #     +      	   NROOT,NSOIL,Z0,CZIL,XLAI,CSOIL,PTU)

    # ----------------------------------------------------------------------
    #  INITIALIZE PRECIPITATION LOGICALS.    初始化沉淀逻辑。
    # CALCULATE ROOT DISTRIBUTION.  PRESENT VERSION ASSUMES UNIFORM     计算根分布。当前版本假设一致
    # DISTRIBUTION BASED ON SOIL LAYER DEPTHS.    基于土层深度的分布。

    for I in range(0, int(NROOT)):
        RTDIS[I] = -SLDPTH[I] / ZSOIL[int(NROOT) - 1]

    # ----------------------------------------------------------------------
    SNOWNG = False
    FRZGRA = False

    # ----------------------------------------------------------------------
    # IF SEA-ICE CASE, ASSIGN DEFAULT WATER-EQUIV SNOW ON TOP  如果是海冰情况，则在顶部指定默认的WATER-EQUIV SNOW
    # ----------------------------------------------------------------------
    if ICE == 1:
        SNEQV = 0.01
        SNOWH = 0.05
    # endIF

    # ----------------------------------------------------------------------
    # IF INPUT SNOWPACK IS NONZERO, THEN COMPUTE SNOW DENSITY "SNDENS" AND  如果输入积雪非零，则计算雪密度“SNDENS”和
    #   SNOW THERMAL CONDUCTIVITY "SNCOND" (NOTE THAT CSNOW IS A FUNCTION  雪热导率“SNCOND”（注意，CSNOW是一个函数
    #   SUBROUTINE)  子程序）
    # ----------------------------------------------------------------------
    if SNEQV == 0.0:  # SNEQV是液态水等效雪深
        SNDENS = 0.0  # SNDENS是雪密度
        SNOWH = 0.0  # SNOWH是雪实际深度
        SNCOND = 1.0  # SNCOND是雪热导率
    else:
        if SNOWH != 0.0:
            SNDENS = SNEQV / SNOWH  # SNOWH是雪实际深度
            SNCOND = CSNOW(SNDENS)  # SNCOND是雪热导率 CSNOW是一个函数
        else:
            SNDENS = 0.0
            SNCOND = CSNOW(SNDENS)
        # if (SNOWH != 0.0):
        #     SNDENS = SNEQV / SNOWH
        #     SNCOND = CSNOW(SNDENS)
        # else:
        #     SNDENS = 0.0
        #     SNCOND = CSNOW(SNDENS)
    # endIF

    # ----------------------------------------------------------------------
    # DETERMINE IF IT'S PRECIPITATING AND WHAT KIND OF PRECIP IT IS.   确定是否有沉淀以及沉淀的类型。
    # IF IT'S PRCPING AND THE AIR TEMP IS COLDER THAN 0 C, IT'S SNOWING   如果是下雨，气温低于0摄氏度，那就是下雪了
    # IF IT'S PRCPING AND THE AIR TEMP IS WARMER THAN 0 C, BUT THE GRND   如果天气很冷，气温高于0摄氏度，但地面温度
    # TEMP IS COLDER THAN 0 C, FREEZING RAIN IS PRESUMED TO BE FALLING.   温度低于0摄氏度，预计将有冻雨
    # ----------------------------------------------------------------------
    if PRCP > 0.0:
        if SFCTMP <= TFREEZ:
            SNOWNG = True
        else:
            if T1 <= TFREEZ:
                FRZGRA = True
    # endIF
    # endIF
    # ----------------------------------------------------------------------
    # IF EITHER PRCP FLAG IS SET, DETERMINE NEW SNOWFALL (CONVERTING PRCP  如果设置了任一PRCP标志，则确定新的降雪量（转换PRCP
    # RATE FROM KG M-2 S-1 TO A LIQUID EQUIV SNOW DEPTH IN METERS) AND ADD  速率从KG M-2 S-1到液体当量雪深（米），并添加
    # IT TO THE EXISTING SNOWPACK.  它是现有积雪的一部分。
    # NOTE THAT SINCE ALL PRECIP IS ADDED TO SNOWPACK, NO PRECIP INFILTRATES 注意，由于所有PRECIP都添加到积雪中，因此没有PRECIP渗入
    # INTO THE SOIL SO THAT PRCP1 IS SET TO ZERO. 使PRCP1设置为零。
    # ----------------------------------------------------------------------
    if SNOWNG or FRZGRA:
        SN_NEW = PRCP * DT * 0.001
        SNEQV = SNEQV + SN_NEW
        PRCP1 = 0.0
        # ----------------------------------------------------------------------
        # UPDATE SNOW DENSITY BASED ON NEW SNOWFALL, USING OLD AND NEW SNOW. 使用旧雪和新雪，根据新降雪更新雪密度。
        # UPDATE SNOW THERMAL CONDUCTIVITY 更新雪导热系数
        # ----------------------------------------------------------------------
        SNDENS = SNOW_NEW(SFCTMP, SN_NEW, SNOWH, SNDENS)
        SNCOND = CSNOW(SNDENS)
    else:
        # ----------------------------------------------------------------------
        # PRECIP IS LIQUID (RAIN), HENCE SAVE IN THE PRECIP VARIABLE THAT  PRECIP是液体（雨），因此保存在PRECIP变量中
        # LATER CAN WHOLELY OR PARTIALLY INFILTRATE THE SOIL (ALONG WITH   之后可能会全部或部分渗入土壤（以及
        # ANY CANOPY "DRIP" ADDED TO THIS LATER)  任何雨篷“滴水”后来添加到此）
        # ----------------------------------------------------------------------
        PRCP1 = PRCP

    # endIF

    # ----------------------------------------------------------------------
    # DETERMINE SNOWCOVER AND ALBEDO OVER LAND.  确定陆地上的积雪和反照率。
    # ----------------------------------------------------------------------
    if ICE == 0:
        # ----------------------------------------------------------------------
        # IF SNOW DEPTH=0, SET SNOW FRACTION=0, ALBEDO=SNOW FREE ALBEDO. 如果雪深=0，则将雪分数设置为0，反照率设置为无雪反照率。
        # ----------------------------------------------------------------------
        if SNEQV == 0.0:
            SNCOVR = 0.0
            ALBEDO = ALB
        else:
            # ----------------------------------------------------------------------
            # DETERMINE SNOW FRACTIONAL COVERAGE.  确定雪覆盖率。
            # DETERMINE SURFACE ALBEDO MODIFICATION DUE TO SNOWDEPTH STATE. 确定由于积雪深度状态引起的地表反照率变化。
            # ----------------------------------------------------------------------
            SNCOVR = SNFRAC(SNEQV, SNUP, SALP, SNOWH, SNCOVR)
            ALBEDO = ALCALC(ALB, SNOALB, SHDFAC, SHDMIN, SNCOVR, TSNOW, ALBEDO)
        # endIF

    else:
        # ----------------------------------------------------------------------
        # SNOW COVER, ALBEDO OVER SEA-ICE 冰雪覆盖、海冰反照率
        # ----------------------------------------------------------------------
        SNCOVR = 1.0
        # changed in version 2.6 on June 2nd 2003  于2003年6月2日在2.6版中更改
        ALBEDO = 0.60
        ALBEDO = 0.65
    # endIF

    # ----------------------------------------------------------------------
    # THERMAL CONDUCTIVITY FOR SEA-ICE CASE  海冰情况下的导热系数
    # ----------------------------------------------------------------------
    if ICE == 1:
        DF1 = 2.2
    else:
        # ----------------------------------------------------------------------
        # NEXT CALCULATE THE SUBSURFACE HEAT FLUX, WHICH FIRST REQUIRES 接下来计算地下热流，首先需要
        # CALCULATION OF THE THERMAL DIFFUSIVITY.  TREATMENT OF THE 热扩散率的计算。治疗
        # LATTER FOLLOWS THAT ON PAGES 148-149 FROM "HEAT TRANSFER IN 后者遵循第148-149页的“热传递”
        # COLD CLIMATES", BY V. J. LUNARDINI (PUBLISHED IN 1981 《寒冷气候》，V.J.LUNARDINI著（1981年出版）
        # BY VAN NOSTRAND REINHOLD CO.) I.E. TREATMENT OF TWO CONTIGUOUS 由VAN NOSTRAND REINHOLD CO。），即处理两个相邻的
        # "PLANE PARALLEL" MEDIUMS (NAMELY HERE THE FIRST SOIL LAYER “平面平行”介质（即第一土层
        # AND THE SNOWPACK LAYER, IF ANY). THIS DIFFUSIVITY TREATMENT 以及积雪层（如果有的话）。这种扩散率处理
        # BEHAVES WELL FOR BOTH ZERO AND NONZERO SNOWPACK, INCLUDING THE 对于零积雪和非零积雪，包括
        # LIMIT OF VERY THIN SNOWPACK.  THIS TREATMENT ALSO ELIMINATES 极薄积雪的限制。这种治疗也消除了
        # THE NEED TO IMPOSE AN ARBITRARY UPPER BOUND ON SUBSURFACE 需要对次表面施加任意上限
        # HEAT FLUX WHEN THE SNOWPACK BECOMES EXTREMELY THIN. 积雪变得非常薄时的热通量。
        # ----------------------------------------------------------------------
        # FIRST CALCULATE THERMAL DIFFUSIVITY OF TOP SOIL LAYER, USING  首先计算表层土壤的热扩散率，使用
        # BOTH THE FROZEN AND LIQUID SOIL MOISTURE, FOLLOWING THE 冻土和液态土壤水分
        # SOIL THERMAL DIFFUSIVITY FUNCTION OF PETERS-LIDARD ET AL. PETERS-LIDARD等人的土壤热扩散函数。
        # (1998,JAS, VOL 55, 1209-1224), WHICH REQUIRES THE SPECIFYING （1998，JAS，第55卷，1209-1224），其中要求指定
        # THE QUARTZ CONTENT OF THE GIVEN SOIL CLASS (SEE ROUTINE REDPRM) 给定土壤类别的石英含量（见常规REDPRM）
        # ----------------------------------------------------------------------
        DF1 = TDFCND(SMC[0], QUARTZ, SMCMAX, SH2O[0])

        # ----------------------------------------------------------------------
        # NEXT ADD SUBSURFACE HEAT FLUX REDUCTION EFFECT FROM THE  接下来，添加来自
        # OVERLYING GREEN CANOPY, ADAPTED FROM SECTION 2.1.2 OF 上覆绿色树冠，改编自第2.1.2节
        # PETERS-LIDARD ET AL. (1997, JGR, VOL 102(D4)) PETERS-LIDARD等人（1997年，JGR，第102卷（D4））
        # ----------------------------------------------------------------------
        DF1 = DF1 * math.exp(SBETA * SHDFAC)
    # endIF

    # ----------------------------------------------------------------------
    # FINALLY "PLANE PARALLEL" SNOWPACK EFFECT FOLLOWING 最后“平面平行”积雪效应如下
    # V.J. LINARDINI REFERENCE CITED ABOVE. NOTE THAT DTOT IS 上文引用的V.J.LINARDINI参考文献。注意，DTOT是
    # COMBINED DEPTH OF SNOWDEPTH AND THICKNESS OF FIRST SOIL LAYER 积雪深度和第一土层厚度的组合深度
    # ----------------------------------------------------------------------
    DSOIL = -(0.5 * ZSOIL[0])

    if SNEQV == 0.:
        SSOIL = DF1 * (T1 - STC[0]) / DSOIL
    else:
        DTOT = SNOWH + DSOIL
        FRCSNO = SNOWH / DTOT
        FRCSOI = DSOIL / DTOT
        #
        # 1. HARMONIC MEAN (SERIES FLOW)  1、谐波平均值（串联流）
        #        DF1 = (SNCOND*DF1)/(FRCSOI*SNCOND+FRCSNO*DF1)
        DF1H = (SNCOND * DF1) / (FRCSOI * SNCOND + FRCSNO * DF1)
        # 2. ARITHMETIC MEAN (PARALLEL FLOW) 2、算术平均值（平行流）
        #        DF1 = FRCSNO*SNCOND + FRCSOI*DF1
        DF1A = FRCSNO * SNCOND + FRCSOI * DF1
        #
        # 3. GEOMETRIC MEAN (INTERMEDIATE BETWEEN HARMONIC AND ARITHMETIC MEAN) 3、几何平均值（介于调和平均值和算术平均值之间）
        #        DF1 = (SNCOND**FRCSNO)*(DF1**FRCSOI)
        # TEST - MBEK, 10 Jan 2002
        # weigh DF by snow fraction
        #        DF1 = DF1H*SNCOVR + DF1A*(1.0-SNCOVR)
        #        DF1 = DF1H*SNCOVR + DF1*(1.0-SNCOVR)
        DF1 = DF1A * SNCOVR + DF1 * (1.0 - SNCOVR)

        # ----------------------------------------------------------------------
        # CALCULATE SUBSURFACE HEAT FLUX, SSOIL, FROM FINAL THERMAL DIFFUSIVITY 根据最终热扩散率计算地下热通量SSOIL
        # OF SURFACE MEDIUMS, DF1 ABOVE, AND SKIN TEMPERATURE AND TOP 表面介质，DF1以上，皮肤温度和顶部
        # MID-LAYER SOIL TEMPERATURE 中层土壤温度
        # ----------------------------------------------------------------------
        SSOIL = DF1 * (T1 - STC[0]) / DTOT
    # endIF

    # ----------------------------------------------------------------------
    # DETERMINE SURFACE ROUGHNESS OVER SNOWPACK USING SNOW CONDITION FROM 根据以下条件确定积雪表面粗糙度：
    # THE PREVIOUS TIMESTEP. 上一个时间步长
    # ----------------------------------------------------------------------
    if SNCOVR > 0.:
        SNOWZ0(SNCOVR, Z0)
    # endIF

    # ----------------------------------------------------------------------
    # NEXT CALL ROUTINE SFCDIF TO CALCULATE THE SFC EXCHANGE COEF (CH) FOR 下一个调用例程SFCDIF计算的SFC交换COEF（CH）
    # HEAT AND MOISTURE. 热量和湿度。
    #
    # NOTE ###
    # COMMENT OUT CALL SFCDIF, IF SFCDIF ALREADY CALLED IN CALLING PROGRAM 如果调用程序中已经调用了SFCDIF，则注释掉CALL SFCDIF
    # (SUCH AS IN COUPLED ATMOSPHERIC MODEL). （例如在耦合大气模型中）
    #
    # NOTE ###
    # DO NOT CALL SFCDIF UNTIL AFTER ABOVE CALL TO REDPRM, IN CASE 在上述呼叫REDPRM之前，不要呼叫SFCDIF，以防
    # ALTERNATIVE VALUES OF ROUGHNESS LENGTH (Z0) AND ZILINTINKEVICH COEF 粗糙度长度（Z0）和齐林丁克维奇系数的替代值
    # (CZIL) ARE SET THERE VIA NAMELIST I/O. （CZIL）通过名称列表输入/输出设置在那里。
    #
    # NOTE ###
    # ROUTINE SFCDIF RETURNS A CH THAT REPRESENTS THE WIND SPD TIMES THE 例程SFCDIF返回一个CH，该CH表示风速SPD乘以
    # "ORIGINAL" NONDIMENSIONAL "CH" TYPICAL IN LITERATURE.  HENCE THE CH 文学中典型的“原始”无量纲“CH”。因此，CH
    # RETURNED FROM SFCDIF HAS UNITS OF M/S.  THE IMPORTANT COMPANION 从SFCDIF返回的单位为M/S。重要的伴星
    # COEFFICIENT OF CH, CARRIED HERE AS "RCH", IS THE CH FROM SFCDIF TIMES CH的系数，此处以“RCH”表示，是SFCDIF次的CH
    # AIR DENSITY AND PARAMETER "CP".  "RCH" IS COMPUTED IN "CALL PENMAN". 空气密度和参数“CP”。“RCH”在“CALL PENMAN”中计算。
    # RCH RATHER THAN CH IS THE COEFF USUALLY INVOKED LATER IN EQNS. RCH而不是CH是通常稍后在等式中调用的系数。
    #
    # NOTE ###
    # SFCDIF ALSO RETURNS THE SURFACE EXCHANGE COEFFICIENT FOR MOMENTUM, CM, SFCDIF还返回动量的表面交换系数CM，
    # ALSO KNOWN AS THE SURFACE DRAGE COEFFICIENT, BUT CM IS NOT USED HERE. 也称为表面阻力系数，但此处不使用CM。
    # ----------------------------------------------------------------------
    # CALC VIRTUAL TEMPS AND VIRTUAL POTENTIAL TEMPS NEEDED BY SUBROUTINES 计算子程序所需的虚拟温度和虚拟电位温度
    # SFCDIF AND PENMAN. SFCDIF和PENMAN
    # ----------------------------------------------------------------------
    T2V = SFCTMP * (1.0 + 0.61 * Q2)
    # ----------------------------------------------------------------------
    # COMMENT OUT BELOW 2 LINES IF CALL SFCDIF IS COMMENTED OUT, I.E. IN THE 如果CALL SFCDIF被注释掉，即在
    # COUPLED MODEL. 耦合模型。
    # ----------------------------------------------------------------------
    T1V = T1 * (1.0 + 0.61 * Q2)
    TH2V = TH2 * (1.0 + 0.61 * Q2)

    SFCDIF(ZLVL, Z0, T1V, TH2V, SFCSPD, CZIL, CM, CH)

    # ----------------------------------------------------------------------
    # CALCULATE TOTAL DOWNWARD RADIATION (SOLAR PLUS LONGWAVE) NEEDED IN  计算所需的总向下辐射（太阳加长波）
    # PENMAN EP SUBROUTINE THAT FOLLOWS PENMAN EP子程序
    # ----------------------------------------------------------------------
    FDOWN = SOLDN * (1.0 - ALBEDO) + LWDN
    # FDOWN = SOLNET + LWDN

    # ----------------------------------------------------------------------
    # CALL PENMAN SUBROUTINE TO CALCULATE POTENTIAL EVAPORATION (ETP), AND 调用PENMAN子程序计算潜在蒸发量（ETP），以及
    # OTHER PARTIAL PRODUCTS AND SUMS SAVE IN COMMON/RITE FOR LATER 其他部分乘积和总和保存在公共/仪式中供以后使用
    # CALCULATIONS. 计算
    # ----------------------------------------------------------------------
    (RCH, RR) = PENMAN(SFCTMP, SFCPRS, CH, T2V, TH2, PRCP, FDOWN, T24, SSOIL,
                       Q2, Q2SAT, ETP, RCH, EPSCA, RR, SNOWNG, FRZGRA,
                       DQSDT2, FLX2)

    # ----------------------------------------------------------------------
    # CALL CANRES TO CALCULATE THE CANOPY RESISTANCE AND CONVERT IT INTO PC 调用CANRES计算冠层阻力并将其转换为PC
    # IF NONZERO GREENNESS FRACTION 如果绿度分数非零
    # ----------------------------------------------------------------------
    if SHDFAC > 0.:
        # ----------------------------------------------------------------------
        #  FROZEN GROUND EXTENSION: TOTAL SOIL WATER "SMC" WAS REPLACED 冻土扩展：替换土壤总水分“SMC”
        #  BY UNFROZEN SOIL WATER "SH2O" IN CALL TO CANRES BELOW 通过调用下面的CANRES中的未冻结土壤水“SH2O”
        # ----------------------------------------------------------------------
        CANRES(SOLDN, CH, SFCTMP, Q2, SFCPRS, SH2O, ZSOIL, NSOIL,
               SMCWLT, SMCREF, RSMIN, RC, PC, NROOT, Q2SAT, DQSDT2,
               TOPT, RSMAX, RGL, HS, XLAI,
               RCS, RCT, RCQ, RCSOIL)

    # endIF

    # ----------------------------------------------------------------------
    # NOW DECIDE MAJOR PATHWAY BRANCH TO TAKE DEPENDING ON WHETHER SNOWPACK 现在，根据是否积雪，决定要走的主要道路
    # EXISTS OR NOT: 是否存在：
    # ----------------------------------------------------------------------
    ESNOW = 0.0
    if SNEQV == 0.0:
        NOPAC(ETP, ETA, PRCP, SMC, SMCMAX, SMCWLT,
              SMCREF, SMCDRY, CMC, CMCMAX, NSOIL, DT, SHDFAC,
              SBETA, Q2, T1, SFCTMP, T24, TH2, FDOWN, F1, SSOIL,
              STC, EPSCA, BEXP, PC, RCH, RR, CFACTR,
              SH2O, SLOPE, KDT, FRZX, PSISAT, ZSOIL,
              DKSAT, DWSAT, TBOT, ZBOT, RUNOFF1, RUNOFF2,
              RUNOFF3, EDIR, EC, ET, ETT, NROOT, ICE, RTDIS,
              QUARTZ, FXEXP, CSOIL,
              BETA, DRIP, DEW, FLX1, FLX2, FLX3)
    else:
        ETP1 = 0.0
        SNOPAC(ETP, ETA, PRCP, PRCP1, SNOWNG, SMC, SMCMAX, SMCWLT,
               SMCREF, SMCDRY, CMC, CMCMAX, NSOIL, DT,
               SBETA, DF1,
               Q2, T1, SFCTMP, T24, TH2, FDOWN, F1, SSOIL, STC, EPSCA,
               SFCPRS, BEXP, PC, RCH, RR, CFACTR, SNCOVR, SNEQV, SNDENS,
               SNOWH, SH2O, SLOPE, KDT, FRZX, PSISAT, SNUP,
               ZSOIL, DWSAT, DKSAT, TBOT, ZBOT, SHDFAC, RUNOFF1,
               RUNOFF2, RUNOFF3, EDIR, EC, ET, ETT, NROOT, SNOMLT,
               ICE, RTDIS, QUARTZ, FXEXP, CSOIL,
               BETA, DRIP, DEW, FLX1, FLX2, FLX3, ESNOW)
        ESNOW = ETA
    # endIF

    # ----------------------------------------------------------------------
    #   PREPARE SENSIBLE HEAT (H) FOR RETURN TO PARENT MODEL 准备显热（H）以返回到父模型
    # ----------------------------------------------------------------------
    SHEAT = -(CH * CP * SFCPRS) / (R * T2V) * (TH2 - T1)

    # ----------------------------------------------------------------------
    #  CONVERT UNITS AND/OR SIGN OF TOTAL EVAP (ETA), POTENTIAL EVAP (ETP), 转换总蒸发量（ETA）、潜在蒸发量（ETP）的单位和/或符号，
    #  SUBSURFACE HEAT FLUX (S), AND RUNOFFS FOR WHAT PARENT MODEL EXPECTS 地下热通量，以及母模型预期的径流
    #  CONVERT ETA FROM KG M-2 S-1 TO W M-2 将ETA从KG M-2 S-1转换为W M-2
    # ----------------------------------------------------------------------
    EVP = ETA
    ETA = ETA * LVH2O
    ETP = ETP * LVH2O

    # ----------------------------------------------------------------------
    # CONVERT THE SIGN OF SOIL HEAT FLUX SO THAT: 转换土壤热通量的符号，以便：
    #   SSOIL>0: WARM THE SURFACE  (NIGHT TIME) SSOIL>0：加热表面（夜间）
    #   SSOIL<0: COOL THE SURFACE  (DAY TIME) SSOIL<0：冷却表面（白天）
    # ----------------------------------------------------------------------
    SSOIL = -1.0 * SSOIL

    # ----------------------------------------------------------------------
    #  CONVERT RUNOFF3 (INTERNAL LAYER RUNOFF FROM SUPERSAT) FROM M TO M S-1 将径流3（来自SUPERSAT的内层径流）从M转换为M S-1
    #  AND ADD TO SUBSURFACE RUNOFF/DRAINAGE/BASEFLOW 并添加到地下径流/排水/基流中
    # ----------------------------------------------------------------------
    RUNOFF3 = RUNOFF3 / DT
    RUNOFF2 = RUNOFF2 + RUNOFF3

    # ----------------------------------------------------------------------
    # TOTAL COLUMN SOIL MOISTURE IN METERS (SOILM) AND ROOT-ZONE 柱状土壤总水分（单位：米）和根区
    # SOIL MOISTURE AVAILABILITY (FRACTION) RELATIVE TO POROSITY/SATURATION 相对于孔隙度/饱和度的土壤水分有效性（分数）
    # ----------------------------------------------------------------------
    SOILM = -1.0 * SMC[0] * ZSOIL[0]
    for K in range(1, NSOIL):
        SOILM = SOILM + SMC[K] * (ZSOIL[K - 1] - ZSOIL[K])
    # ----------------------------------------------------------------------
    # ROOT-ZONE SOIL MOISTURE AVAILABILITY (FRACTION) RELATIVE 根区土壤水分有效性（分数）相对
    # TO POROSITY/SATURATION (SOILW; aka, MSTAVRZ) 孔隙度/饱和度（SOILW；aka，MSTAVRZ）
    # ----------------------------------------------------------------------
    SOILWM = -1.0 * (SMCMAX - SMCWLT) * ZSOIL[0]
    SOILWW = -1.0 * (SMC[0] - SMCWLT) * ZSOIL[0]
    for K in range(1, int(NROOT)):
        SOILWM = SOILWM + (SMCMAX - SMCWLT) * (ZSOIL[K - 1] - ZSOIL[K])
        SOILWW = SOILWW + (SMC[K] - SMCWLT) * (ZSOIL[K - 1] - ZSOIL[K])
    SOILW = SOILWW / SOILWM

    # ----------------------------------------------------------------------
    # TOTAL COL SOIL MOISTURE AVAIL RELATIVE TO POROSITY/SATURATION (SOILT) 与孔隙度/饱和度（土壤）相关的土壤总水分利用率
    #  (aka, MSTAVTOT)
    # ----------------------------------------------------------------------
    SOILTM = -1.0 * (SMCMAX - SMCWLT) * ZSOIL[0]
    SOILTW = -1.0 * (SMC[0] - SMCWLT) * ZSOIL[0]
    for K in range(1, NSOIL):
        SOILTM = SOILTM + (SMCMAX - SMCWLT) * (ZSOIL[K - 1] - ZSOIL[K])
        SOILTW = SOILTW + (SMC[K] - SMCWLT) * (ZSOIL[K - 1] - ZSOIL[K])
    SOILT = SOILTW / SOILTM

    # ----------------------------------------------------------------------
    # END SUBROUTINE SFLX
    # ----------------------------------------------------------------------
    return (CMC, T1, STC, SMC,
            SH2O, SNOWH, SNEQV, ALBEDO, CH, CM,
            EVP, ETA, SHEAT,
            EC, EDIR, ET, ETT, ESNOW, DRIP, DEW,
            BETA, ETP, SSOIL,
            FLX1, FLX2, FLX3,
            SNOMLT, SNCOVR,
            RUNOFF1, RUNOFF2, RUNOFF3,
            RC, PC, RCS, RCT, RCQ, RCSOIL,
            SOILW, SOILT, SOILM)


# end

# -------  THE FOLLOWING SUBROUTINES ARE IN ALPHABETICAL ORDER  -------- 以下子程序按字母顺序排列：

# -- 2. PHYSICS SUBROUTINE ==>  SUBROUTINE ALCALC ---------------------- 物理子程序==>子程序 ALCALC

def ALCALC(ALB, SNOALB, SHDFAC, SHDMIN, SNCOVR, TSNOW, ALBEDO):
    # ----------------------------------------------------------------------
    # CALCULATE ALBEDO INCLUDING SNOW EFFECT (0 -> 1) 计算ALBEDO，包括雪效应（0->1）
    #   ALB     SNOWFREE ALBEDO                  ALB  无雪ALBEDO
    #   SNOALB  MAXIMUM (DEEP) SNOW ALBEDO       SNOALB 最大（深）雪反照率
    #   SHDFAC    AREAL FRACTIONAL COVERAGE OF GREEN VEGETATION     SHDFAC 绿色植被的面积分数覆盖率
    #   SHDMIN    MINIMUM AREAL FRACTIONAL COVERAGE OF GREEN VEGETATION    SHDMIN 绿色植被的最小面积分数覆盖率
    #   SNCOVR  FRACTIONAL SNOW COVER     SNCOVR 部分积雪
    #   ALBEDO  SURFACE ALBEDO INCLUDING SNOW EFFECT   ALBEDO 表面反照率，包括雪效应
    #   TSNOW   SNOW SURFACE TEMPERATURE (K)     TSNOW  雪表面温度（K）
    # ----------------------------------------------------------------------
    # real ALB, SNOALB, SHDFAC, SHDMIN, SNCOVR, ALBEDO, TSNOW  真实ALB, SNOALB, SHDFAC, SHDMIN, SNCOVR, ALBEDO,TSNOW

    # ----------------------------------------------------------------------
    # SNOALB IS ARGUMENT REPRESENTING MAXIMUM ALBEDO OVER DEEP SNOW, SNOALB是表示深雪上最大反照率的参数，
    # AS PASSED INTO SFLX, AND ADAPTED FROM THE SATELLITE-BASED MAXIMUM 传递到SFLX，并根据基于卫星的最大值进行调整
    # SNOW ALBEDO FIELDS PROVIDED BY D. ROBINSON AND G. KUKLA D.ROBINSON和G.KUKLA提供的雪地反照率
    # (1985, JCAM, VOL 24, 402-411) （1985年，JCAM，VOL 24402-411）
    # ----------------------------------------------------------------------
    # changed in version 2.6 on June 2nd 2003  2003年6月2日在版本2.6中更改
    # ALBEDO = ALB + (1.0-(SHDFAC-SHDMIN))*SNCOVR*(SNOALB-ALB)
    ALBEDO = ALB + SNCOVR * (SNOALB - ALB)
    if ALBEDO > SNOALB:
        ALBEDO = SNOALB

    #     BASE FORMULATION (DICKINSON ET AL., 1986, COGLEY ET AL, 1990) 基础方法（Dickinson等人，1986年；Cogley等人，1990年）
    #          IF (TSNOW.LE.263.16) THEN
    #            ALBEDO=SNOALB
    #          ELSE
    #            IF (TSNOW.LT.273.16) THEN
    #              TM=0.1*(TSNOW-263.16)
    #              ALBEDO=0.5*((0.9-0.2*(TM**3))+(0.8-0.16*(TM**3)))
    #            ELSE
    #              ALBEDO=0.67
    #            ENDIF
    #          ENDIF

    #     ISBA FORMULATION (VERSEGHY, 1991; BAKER ET AL., 1990) ISBA方法（Verseghy，1991；Baker等人，1990）
    #          IF (TSNOW.LT.273.16) THEN
    #            ALBEDO=SNOALB-0.008*DT/86400
    #          ELSE
    #            ALBEDO=(SNOALB-0.5)*EXP(-0.24*DT/86400)+0.5
    #          ENDIF

    # ----------------------------------------------------------------------
    # END SUBROUTINE ALCALC
    # ----------------------------------------------------------------------
    return ALBEDO


# end

# -- 3. PHYSICS SUBROUTINE ==>  SUBROUTINE CANRES  --------------------物理子程序==>子程序 CANRES
def CANRES(SOLAR, CH, SFCTMP, Q2, SFCPRS, SMC, ZSOIL, NSOIL,
           SMCWLT, SMCREF, RSMIN, RC, PC, NROOT, Q2SAT, DQSDT2,
           TOPT, RSMAX, RGL, HS, XLAI,
           RCS, RCT, RCQ, RCSOIL):
    # ----------------------------------------------------------------------
    # SUBROUTINE CANRES 子程序CANRES
    # ----------------------------------------------------------------------
    # CALCULATE CANOPY RESISTANCE WHICH DEPENDS ON INCOMING SOLAR RADIATION, 计算取决于入射太阳辐射的雨棚阻力
    # AIR TEMPERATURE, ATMOSPHERIC WATER VAPOR PRESSURE DEFICIT AT THE 空气温度、大气水蒸气压力赤字
    # LOWEST MODEL LEVEL, AND SOIL MOISTURE (PREFERABLY UNFROZEN SOIL 最低模型水位和土壤湿度（最好是未冻结的土壤
    # MOISTURE RATHER THAN TOTAL) 水分（而非总量）
    # ----------------------------------------------------------------------
    # SOURCE:  JARVIS (1976), NOILHAN AND PLANTON (1989, MWR), JACQUEMIN AND
    # NOILHAN (1990, BLM)
    # SEE ALSO:  CHEN ET AL (1996, JGR, VOL 101(D3), 7251-7268), EQNS 12-14
    # AND TABLE 2 OF SEC. 3.1.2
    # ----------------------------------------------------------------------
    # INPUT:
    #   SOLAR   INCOMING SOLAR RADIATION 太阳入射太阳辐射
    #   CH      SURFACE EXCHANGE COEFFICIENT FOR HEAT AND MOISTURE 热湿表面交换系数
    #   SFCTMP  AIR TEMPERATURE AT 1ST LEVEL ABOVE GROUND 地面以上第一层的空气温度
    #   Q2      AIR HUMIDITY AT 1ST LEVEL ABOVE GROUND 地面以上第一层的空气湿度
    #   Q2SAT   SATURATION AIR HUMIDITY AT 1ST LEVEL ABOVE GROUND 地面以上第一层的饱和空气湿度
    #   DQSDT2  SLOPE OF SATURATION HUMIDITY FUNCTION WRT TEMP 饱和湿度函数WRT温度的斜率
    #   SFCPRS  SURFACE PRESSURE 表面压力
    #   SMC     VOLUMETRIC SOIL MOISTURE 土壤体积含水量
    #   ZSOIL   SOIL DEPTH (NEGATIVE SIGN, AS IT IS BELOW GROUND) 土壤深度（负号，因为它低于地面）
    #   NSOIL   NO. OF SOIL LAYERS 土层数量
    #   NROOT   NO. OF SOIL LAYERS IN ROOT ZONE (1.LE.NROOT.LE.NSOIL) 根区土层数量（1.LE.NROOT.LE.NSOIL）
    #   XLAI    LEAF AREA INDEX 叶面积指数
    #   SMCWLT  WILTING POINT 萎蔫点
    #   SMCREF  REFERENCE SOIL MOISTURE (WHERE SOIL WATER DEFICIT STRESS 参考土壤湿度（其中土壤水分亏缺应力
    #             SETS IN) 开始）
    # RSMIN, RSMAX, TOPT, RGL, HS ARE CANOPY STRESS PARAMETERS SET IN 是树冠应力参数设置在
    #   SURBOUTINE REDPRM
    # OUTPUT:
    #   PC  PLANT COEFFICIENT 植物系数
    #   RC  CANOPY RESISTANCE 顶篷阻力
    # ---------------------------------------------------------------------

    CP = 1004.5
    RD = 287.04
    SIGMA = 5.67E-8
    SLV = 2.501000E6

    # ----------------------------------------------------------------------
    # INITIALIZE CANOPY RESISTANCE MULTIPLIER TERMS. 初始化冠层阻力乘数项。
    # ----------------------------------------------------------------------
    RCS = 0.0
    RCT = 0.0
    RCQ = 0.0
    RCSOIL = 0.0
    RC = 0.0

    # ----------------------------------------------------------------------
    # CONTRIBUTION DUE TO INCOMING SOLAR RADIATION 入射太阳辐射的贡献
    # ----------------------------------------------------------------------
    FF = 0.55 * 2.0 * SOLAR / (RGL * XLAI)
    RCS = (FF + RSMIN / RSMAX) / (1.0 + FF)
    RCS = max(RCS, 0.0001)

    # ----------------------------------------------------------------------
    # CONTRIBUTION DUE TO AIR TEMPERATURE AT FIRST MODEL LEVEL ABOVE GROUND 地面以上第一模型层的空气温度的贡献
    # RCT EXPRESSION FROM NOILHAN AND PLANTON (1989, MWR). 来自NOILHAN和PLANTON的RCT表达（1989，MWR）。
    # ----------------------------------------------------------------------
    RCT = 1.0 - 0.0016 * ((TOPT - SFCTMP) ** 2.0)
    RCT = max(RCT, 0.0001)

    # ----------------------------------------------------------------------
    # CONTRIBUTION DUE TO VAPOR PRESSURE DEFICIT AT FIRST MODEL LEVEL. 由于第一模型水平上的蒸汽压赤字造成的贡献。
    # RCQ EXPRESSION FROM SSIB 来自SSIB的RCQ表达式
    # ----------------------------------------------------------------------
    RCQ = 1.0 / (1.0 + HS * (Q2SAT - Q2))
    RCQ = max(RCQ, 0.01)

    # ----------------------------------------------------------------------
    # CONTRIBUTION DUE TO SOIL MOISTURE AVAILABILITY. 土壤水分有效性的贡献。
    # DETERMINE CONTRIBUTION FROM EACH SOIL LAYER, THEN ADD THEM UP. 确定各土层的贡献，然后将它们相加。
    # ----------------------------------------------------------------------
    GX = (SMC[1] - SMCWLT) / (SMCREF - SMCWLT)
    if GX > 1.:
        GX = 1.
    if GX < 0.:
        GX = 0.
    # ----------------------------------------------------------------------
    # USE SOIL DEPTH AS WEIGHTING FACTOR 使用土壤深度作为加权因子
    # ----------------------------------------------------------------------
    NSOLD = 20
    import numpy as np
    PART = np.zeros(NSOLD, dtype=float)
    PART[0] = (ZSOIL[0] / (ZSOIL[NROOT] + 1)) * GX
    # ----------------------------------------------------------------------
    # USE ROOT DISTRIBUTION AS WEIGHTING FACTOR 使用根分布作为权重因子
    #      PART(1) = RTDIS(1) * GX
    # ----------------------------------------------------------------------
    for K in range(1, NROOT):
        GX = (SMC[K] - SMCWLT) / (SMCREF - SMCWLT)
        if GX > 1.:
            GX = 1.
        if GX < 0.:
            GX = 0.
        # ----------------------------------------------------------------------
        # USE SOIL DEPTH AS WEIGHTING FACTOR 使用土壤深度作为加权因子
        # ----------------------------------------------------------------------
        PART[K] = ((ZSOIL[K] - ZSOIL[K - 1]) / ZSOIL[NROOT - 1]) * GX
    # ----------------------------------------------------------------------
    # USE ROOT DISTRIBUTION AS WEIGHTING FACTOR 使用根分布作为权重因子
    #        PART(K) = RTDIS(K) * GX
    # ----------------------------------------------------------------------

    for K in range(0, NROOT):
        RCSOIL = RCSOIL + PART[K]

    RCSOIL = max(RCSOIL, 0.0001)

    # ----------------------------------------------------------------------
    # DETERMINE CANOPY RESISTANCE DUE TO ALL FACTORS.  CONVERT CANOPY 根据所有因素确定冠层阻力。转换树冠
    # RESISTANCE (RC) TO PLANT COEFFICIENT (PC) TO BE USED WITH POTENTIAL 阻力（RC）与潜在使用的植物系数（PC）
    # EVAP IN DETERMINING ACTUAL EVAP.  PC IS DETERMINED BY: 确定实际蒸发量时的蒸发量。PC由下式确定：
    #   PC * LINERIZED PENMAN POTENTIAL EVAP =
    #   PENMAN-MONTEITH ACTUAL EVAPORATION (CONTAINING RC TERM).
    # ----------------------------------------------------------------------
    RC = RSMIN / (XLAI * RCS * RCT * RCQ * RCSOIL)

    #      TAIR4 = SFCTMP**4.
    #      ST1 = (4.*SIGMA*RD)/CP
    #      SLVCP = SLV/CP
    #      RR = ST1*TAIR4/(SFCPRS*CH) + 1.0
    RR = (4. * SIGMA * RD / CP) * (SFCTMP ** 4.) / (SFCPRS * CH) + 1.0
    DELTA = (SLV / CP) * DQSDT2

    PC = (RR + DELTA) / (RR * (1. + RC * CH) + DELTA)
    return


# ----------------------------------------------------------------------
# END SUBROUTINE CANRES 结束子程序CANRES
# ----------------------------------------------------------------------

def CSNOW(DSNOW):
    UNIT = 0.11631

    # ----------------------------------------------------------------------
    # CSNOW IN UNITS OF CAL/(CM*HR*C), RETURNED IN W/(M*C) CSNOW单位为CAL/（CM*HR*C），返回单位为W/（M*C）
    # BASIC VERSION IS DYACHKOVA EQUATION (1960), FOR RANGE 0.1-0.4  基本版本是DYACHKOVA方程（1960），范围为0.1-0.4
    # ----------------------------------------------------------------------
    C = 0.328 * 10 ** (2.25 * DSNOW)
    CSNOW1 = UNIT * C

    # ----------------------------------------------------------------------
    # DE VAUX EQUATION (1933), IN RANGE 0.1-0.6 DE VAUX方程（1933），范围为0.1-0.6
    # ----------------------------------------------------------------------
    #      CSNOW=0.0293*(1.+100.*DSNOW**2)

    # ----------------------------------------------------------------------
    # E. ANDERSEN FROM FLERCHINGER E 来自弗莱辛格的安徒生
    # ----------------------------------------------------------------------
    #      CSNOW=0.021+2.51*DSNOW**2

    # ----------------------------------------------------------------------
    # END FUNCTION CSNOW 结束函数CSNOW
    # ----------------------------------------------------------------------
    return CSNOW1


# end

# -- 4. PHYSICS SUBROUTINE ==>  SUBROUTINE DEVAP  --------------------- 物理子程序==>子程序DEVAP

def MIN(FX, param):
    pass


def DEVAP(EDIR1, ETP1, SMC, ZSOIL, SHDFAC, SMCMAX, BEXP,
          #      FUNCTION DEVAP (ETP1,SMC,ZSOIL,SHDFAC,SMCMAX,BEXP,
          DKSAT, DWSAT, SMCDRY, SMCREF, SMCWLT, FXEXP):
    # ----------------------------------------------------------------------
    # SUBROUTINE DEVAP DEVAP子程序
    # FUNCTION DEVAP 函数DEVAP
    # ----------------------------------------------------------------------
    # CALCULATE DIRECT SOIL EVAPORATION 计算土壤直接蒸发量
    # ----------------------------------------------------------------------
    # ----------------------------------------------------------------------
    # DIRECT EVAP A FUNCTION OF RELATIVE SOIL MOISTURE AVAILABILITY, LINEAR 直接蒸发量是相对土壤水分有效性的函数，线性
    # WHEN FXEXP=1. 当FXEXP=1时。
    # FX > 1 REPRESENTS DEMAND CONTROL 表示需求控制
    # FX < 1 REPRESENTS FLUX CONTROL 表示流量控制
    # ----------------------------------------------------------------------
    SRATIO = (SMC - SMCDRY) / (SMCMAX - SMCDRY)
    if SRATIO > 0.:
        fx = SRATIO ** FXEXP
        fx = max(MIN(fx, 1.), 0.)
    else:
        fx = 0.
    # endIF

    # ----------------------------------------------------------------------
    # ALLOW FOR THE DIRECT-EVAP-REDUCING EFFECT OF SHADE 考虑到阴影的直接蒸发减少效果
    # ----------------------------------------------------------------------
    #      DEVAP = FX * ( 1.0 - SHDFAC ) * ETP1
    EDIR1 = fx * (1.0 - SHDFAC) * ETP1

    # ----------------------------------------------------------------------
    # END SUBROUTINE DEVAP 结束子程序DEVAP
    # END FUNCTION DEVAP 端函数DEVAP
    # ----------------------------------------------------------------------
    return


# end

# -- 5. PHYSICS SUBROUTINE ==>  SUBROUTINE EVAPO  --------------------- 5.物理子例程==>子例程EVAPO

def EVAPO(ETA1, SMC, NSOIL, CMC, ETP1, DT, ZSOIL,
          SH2O,
          SMCMAX, BEXP, PC, SMCWLT, DKSAT, DWSAT,
          SMCREF, SHDFAC, CMCMAX,
          SMCDRY, CFACTR,
          EDIR1, EC1, ET1, ETT1, SFCTMP, Q2, NROOT, RTDIS, FXEXP, RETURN=None):
    # ----------------------------------------------------------------------
    # SUBROUTINE EVAPO 子程序EVAPO
    # ----------------------------------------------------------------------
    # CALCULATE SOIL MOISTURE FLUX.  THE SOIL MOISTURE CONTENT (SMC - A PER 计算土壤水分通量。土壤含水量（SMC-A每
    # UNIT VOLUME MEASUREMENT) IS A DEPENDENT VARIABLE THAT IS UPDATED WITH 单位体积测量值）是一个因变量，使用
    # PROGNOSTIC EQNS. THE CANOPY MOISTURE CONTENT (CMC) IS ALSO UPDATED. 预测方程。树冠含水量（CMC）也会更新。
    # FROZEN GROUND VERSION:  NEW STATES ADDED: SH2O, AND FROZEN GROUND 冻土版本：新增状态：SH2O和冻土
    # CORRECTION FACTOR, FRZFACT AND PARAMETER SLOPE. 校正系数、FRZFACT和参数斜率。
    # ----------------------------------------------------------------------
    # integer, PARAMETER:: NSOLD=20
    # ----------------------------------------------------------------------
    # EXECUTABLE CODE BEGINS HERE IF THE POTENTIAL EVAPOTRANSPIRATION IS 如果潜在蒸散量为：
    # GREATER THAN ZERO. 大于零。
    # ----------------------------------------------------------------------
    EDIR1 = 0.
    EC1 = 0.
    for K in range(0, NSOIL):
        ET1[K] = 0.

    ETT1 = 0.

    if ETP1 > 0.0:

        # ----------------------------------------------------------------------
        # RETRIEVE DIRECT EVAPORATION FROM SOIL SURFACE.  CALL THIS FUNCTION 恢复土壤表面的直接蒸发。调用此函数
        # ONLY IF VEG COVER NOT COMPLETE. 只有在植被覆盖不完整的情况下。
        # FROZEN GROUND VERSION:  SH2O STATES REPLACE SMC STATES. 冻土版本：SH2O状态取代SMC状态。
        # ----------------------------------------------------------------------
        if SHDFAC < 1.:
            DEVAP(EDIR1, ETP1, SH2O(1), ZSOIL(1), SHDFAC, SMCMAX,
                  # EDIR = DEVAP(ETP1,SH2O(1),ZSOIL(1),SHDFAC,SMCMAX,
                  BEXP, DKSAT, DWSAT, SMCDRY, SMCREF, SMCWLT, FXEXP)
        # endIF

        # ----------------------------------------------------------------------
        # INITIALIZE PLANT TOTAL TRANSPIRATION, RETRIEVE PLANT TRANSPIRATION, 初始化植物总蒸腾，检索植物蒸腾，
        # AND ACCUMULATE IT FOR ALL SOIL LAYERS. 并在所有土层中积累。
        # ----------------------------------------------------------------------
        if SHDFAC > 0.0:
            TRANSP(ET1, NSOIL, ETP1, SH2O, CMC, ZSOIL, SHDFAC, SMCWLT,
                   CMCMAX, PC, CFACTR, SMCREF, SFCTMP, Q2, NROOT, RTDIS)

    for K in range(0, NSOIL):
        ETT1 = ETT1 + ET1[K]

    # ----------------------------------------------------------------------
    # CALCULATE CANOPY EVAPORATION. 计算顶篷蒸发量。
    # IF STATEMENTS TO AVOID TANGENT LINEAR PROBLEMS NEAR CMC=0.0. 避免CMC=0.0附近切线线性问题的IF语句。
    # ----------------------------------------------------------------------
    if CMC > 0.0:
        EC1 = SHDFAC * ((CMC / CMCMAX) ** CFACTR) * ETP1
    else:
        EC1 = 0.0
    # endIF

    # ----------------------------------------------------------------------
    # EC SHOULD BE LIMITED BY THE TOTAL AMOUNT OF AVAILABLE WATER ON THE  EC应受到上可用水总量的限制
    # CANOPY.  -F.CHEN, 18-OCT-1994
    # ----------------------------------------------------------------------
    CMC2MS = CMC / DT
    EC1 = min(CMC2MS, EC1)
    # endIF
    # endIF

    # ----------------------------------------------------------------------
    # TOTAL UP EVAP AND TRANSP TYPES TO OBTAIN ACTUAL EVAPOTRANSP 合计蒸发量和运输类型，以获得实际蒸发量
    # ----------------------------------------------------------------------
    ETA1 = EDIR1 + ETT1 + EC1

    # ----------------------------------------------------------------------
    # END SUBROUTINE EVAPO 结束子例程EVAPO
    # ----------------------------------------------------------------------
    RETURN


# end


# -- 6. PHYSICS SUBROUTINE ==>  SUBROUTINE HRT ------------------------ 物理子程序==>子程序HRT

def HRT(RHSTS, STC, SMC, SMCMAX, NSOIL, ZSOIL, YY, ZZ1,
        TBOT, ZBOT, PSISAT, SH2O, DT, BEXP,
        F1, DF1, QUARTZ, CSOIL, AI, BI, CI):
    # ----------------------------------------------------------------------
    # SUBROUTINE HRT  HRT子程序
    # ----------------------------------------------------------------------
    # CALCULATE THE RIGHT HAND SIDE OF THE TIME TENDENCY TERM OF THE SOIL 计算土壤时间趋势项的右侧
    # THERMAL DIFFUSION EQUATION.  ALSO TO COMPUTE ( PREPARE ) THE MATRIX 热扩散方程。同时计算（准备）矩阵
    # COEFFICIENTS FOR THE TRI-DIAGONAL MATRIX OF THE IMPLICIT TIME SCHEME. 隐式时间格式的三对角矩阵的系数。
    # ----------------------------------------------------------------------
    # integer, PARAMETER:: NSOLD=20

    # ----------------------------------------------------------------------
    # DECLARE WORK ARRAYS NEEDED IN TRI-DIAGONAL IMPLICIT SOLVER 声明三对角隐式解算器中所需的工作数组
    # ----------------------------------------------------------------------
    # real AI(NSOLD)
    # real BI(NSOLD)
    # real CI(NSOLD)

    T0 = 273.15
    TBK, TAVG = 0.0, 0.0
    # ----------------------------------------------------------------------
    # SET SPECIFIC HEAT CAPACITIES OF AIR, WATER, ICE, SOIL MINERAL 设定空气、水、冰、土壤矿物的比热容
    # ----------------------------------------------------------------------
    CAIR = 1004.0
    CH2O = 4.2E6
    CICE = 2.106E6
    # NOTE: CSOIL NOW SET IN ROUTINE REDPRM AND PASSED IN 注：CSOIL现在在常规REDPRM中设置并通过
    #      PARAMETER(CSOIL = 1.26E6)

    # ----------------------------------------------------------------------
    # INITIALIZE LOGICAL FOR SOIL LAYER TEMPERATURE AVERAGING. 初始化土层温度平均的逻辑。
    # ----------------------------------------------------------------------
    ITAVG = True
    #      ITAVG = .FALSE.

    # ----------------------------------------------------------------------
    # BEGIN SECTION FOR TOP SOIL LAYER 表层土壤开始剖面图
    # ----------------------------------------------------------------------
    # CALC THE HEAT CAPACITY OF THE TOP SOIL LAYER 计算表层土壤的热容
    # ----------------------------------------------------------------------
    HCPCT = SH2O[0] * CH2O + (1.0 - SMCMAX) * CSOIL + (SMCMAX - SMC[0]) * CAIR \
            + (SMC[0] - SH2O[0]) * CICE

    # ----------------------------------------------------------------------
    # CALC THE MATRIX COEFFICIENTS AI, BI, AND CI FOR THE TOP LAYER 计算顶层的矩阵系数AI、BI和CI
    # ----------------------------------------------------------------------
    DDZ = 1.0 / (-0.5 * ZSOIL[1])
    AI[0] = 0.0
    CI[0] = (DF1 * DDZ) / (ZSOIL[0] * HCPCT)
    BI[0] = -CI[0] + DF1 / (0.5 * ZSOIL[0] * ZSOIL[0] * HCPCT * ZZ1)

    # ----------------------------------------------------------------------
    # CALCULATE THE VERTICAL SOIL TEMP GRADIENT BTWN THE 1ST AND 2ND SOIL 计算第一层和第二层土壤的垂直土壤温度梯度
    # LAYERS.  THEN CALCULATE THE SUBSURFACE HEAT FLUX. USE THE TEMP 层。然后计算地下热通量。使用温度
    # GRADIENT AND SUBSFC HEAT FLUX TO CALC "RIGHT-HAND SIDE TENDENCY 梯度和SUBSFC热通量计算“右侧趋势”
    # TERMS", OR "RHSTS", FOR TOP SOIL LAYER. “术语”或“RHST”，用于表层土壤。
    # ----------------------------------------------------------------------
    DTSDZ = (STC[0] - STC[1]) / (-0.5 * ZSOIL[1])
    SSOIL = DF1 * (STC[0] - YY) / (0.5 * ZSOIL[0] * ZZ1)
    RHSTS[0] = (DF1 * DTSDZ - SSOIL) / (ZSOIL[0] * HCPCT)

    # ----------------------------------------------------------------------
    # NEXT CAPTURE THE VERTICAL DIFFERENCE OF THE H[E]AT FLUX AT TOP AND 接下来，获取顶部通量的H[E]的垂直差，并
    # BOTTOM OF FIRST SOIL LAYER FOR USE IN HEAT FLUX CONSTRAINT APPLIED TO 用于热通量约束的第一土层底部
    # POTENTIAL SOIL FREEZING/THAWING IN ROUTINE SNKSRC. 常规SNKSRC中的潜在土壤冻结/解冻。
    # ----------------------------------------------------------------------
    QTOT = SSOIL - DF1 * DTSDZ

    # ----------------------------------------------------------------------
    # IF TEMPERATURE AVERAGING INVOKED (ITAVG=TRUE; ELSE SKIP): 如果调用了温度平均（ITAVG=TRUE；否则跳过）：
    # SET TEMP "TSURF" AT TOP OF SOIL COLUMN (FOR USE IN FREEZING SOIL 在土柱顶部设置温度“TSURF”（用于冻土
    # PHYSICS LATER IN FUNCTION SUBROUTINE SNKSRC).  IF SNOWPACK CONTENT IS 物理稍后在函数子程序SNKSRC中）。如果积雪含量为
    # ZERO, THEN TSURF EXPRESSION BELOW GIVES TSURF = SKIN TEMP.  IF 零，则下面的TSURF表达式给出TSURF=皮肤温度IF
    # SNOWPACK IS NONZERO (HENCE ARGUMENT ZZ1=1), THEN TSURF EXPRESSION 如果积雪为非零（因此参数ZZ1=1），则TSURF表达式
    # BELOW YIELDS SOIL COLUMN TOP TEMPERATURE UNDER SNOWPACK.  THEN 下面是积雪下的土壤柱顶部温度。然后
    # CALCULATE TEMPERATURE AT BOTTOM INTERFACE OF 1ST SOIL LAYER FOR USE 计算第一土层底部界面的温度以供使用
    # LATER IN FUNCTION SUBROUTINE SNKSRC 稍后在函数子例程SNKSRC中
    # ----------------------------------------------------------------------
    if ITAVG:
        TSURF = (YY + (ZZ1 - 1) * STC[0]) / ZZ1
        TBK = TBND(STC[0], STC[1], ZSOIL, ZBOT, 1, NSOIL)
    # endIF

    # ----------------------------------------------------------------------
    # CALCULATE FROZEN WATER CONTENT IN 1ST SOIL LAYER. 计算第一土层的冻结水含量。
    # ----------------------------------------------------------------------
    SICE = SMC[0] - SH2O[0]

    # ----------------------------------------------------------------------
    # IF FROZEN WATER PRESENT OR ANY OF LAYER-1 MID-POINT OR BOUNDING 如果存在冻结水，或第1层中点或边界
    # INTERFACE TEMPERATURES BELOW FREEZING, THEN CALL SNKSRC TO 界面温度低于冰点，然后调用SNKSRC
    # COMPUTE HEAT SOURCE/SINK (AND CHANGE IN FROZEN WATER CONTENT) 计算热源/散热器（以及冷冻水含量的变化）
    # DUE TO POSSIBLE SOIL WATER PHASE CHANGE 由于可能的土壤-水相变
    # ----------------------------------------------------------------------
    if ((SICE > 0.) or (TSURF < T0) or
            (STC[0] < T0) or (TBK < T0)):

        if ITAVG:
            TAVG = TMPAVG(TAVG, TSURF, STC[0], TBK, ZSOIL, NSOIL, 1)
        else:
            TAVG = STC[0]
        # endIF
        TSNSR = SNKSRC(TAVG, float(SMC[0]), float(SH2O[0]),
                       ZSOIL, NSOIL, SMCMAX, PSISAT, BEXP, DT, 1, QTOT)
        RHSTS[0] = RHSTS[0] - TSNSR / (ZSOIL[0] * HCPCT)
    # endIF

    # ----------------------------------------------------------------------
    # THIS ENDS SECTION FOR TOP SOIL LAYER. 这是表层土壤的结束部分。
    # ----------------------------------------------------------------------
    # INITIALIZE DDZ2
    # ----------------------------------------------------------------------
    DDZ2 = 0.0

    # ----------------------------------------------------------------------
    # LOOP THRU THE REMAINING SOIL LAYERS, REPEATING THE ABOVE PROCESS 循环穿过剩余土层，重复上述过程
    # (EXCEPT SUBSFC OR "GROUND" HEAT FLUX NOT REPEATED IN LOWER LAYERS) （下层未重复的SUBSFC或“地面”热通量除外）
    # ----------------------------------------------------------------------
    DF1K = DF1
    for k in range(1, NSOIL):
        # ----------------------------------------------------------------------
        # CALCULATE HEAT CAPACITY FOR THIS SOIL LAYER. 计算该土层的热容。
        # ----------------------------------------------------------------------
        HCPCT = SH2O[k] * CH2O + (1.0 - SMCMAX) * CSOIL + (SMCMAX - SMC[k]) * CAIR \
                + (SMC[k] - SH2O[k]) * CICE

        if k != NSOIL - 1:
            # ----------------------------------------------------------------------
            # THIS SECTION FOR LAYER 2 OR GREATER, BUT NOT LAST LAYER. 本节适用于第2层或更高层，但不适用于最后一层。
            # ----------------------------------------------------------------------
            # CALCULATE THERMAL DIFFUSIVITY FOR THIS LAYER. 计算该层的热扩散率。
            # ----------------------------------------------------------------------
            DF1N = TDFCND(SMC[k], QUARTZ, SMCMAX, SH2O[k])

            # ----------------------------------------------------------------------
            # CALC THE VERTICAL SOIL TEMP GRADIENT THRU THIS LAYER 计算穿过该层的垂直土壤温度梯度
            # ----------------------------------------------------------------------
            DENOM = 0.5 * (ZSOIL[k - 1] - ZSOIL[k + 1])
            DTSDZ2 = (STC[k] - STC[k + 1]) / DENOM

            # ----------------------------------------------------------------------
            # CALC THE MATRIX COEF, CI, AFTER CALC'NG ITS PARTIAL PRODUCT 计算矩阵系数CI，在计算其部分乘积之后
            # ----------------------------------------------------------------------
            DDZ2 = 2. / (ZSOIL[k - 1] - ZSOIL[k + 1])
            CI[k] = -DF1N * DDZ2 / ((ZSOIL[k - 1] - ZSOIL[k]) * HCPCT)

            # ----------------------------------------------------------------------
            # IF TEMPERATURE AVERAGING INVOKED (ITAVG=TRUE; ELSE SKIP):  CALCULATE 如果触发温度平均（ITAVG=真；否则跳过）计算
            # TEMP AT BOTTOM OF LAYER. 层底部温度。
            # ----------------------------------------------------------------------
            if ITAVG:
                TBK1 = TBND(STC[k], STC[k + 1], ZSOIL, ZBOT, k, NSOIL)
            # endIF
        else:
            # ----------------------------------------------------------------------
            # SPECIAL CASE OF BOTTOM SOIL LAYER:  CALCULATE THERMAL DIFFUSIVITY FOR 底层土壤的特殊情况：计算热扩散系数
            # BOTTOM LAYER. 底层。
            # ----------------------------------------------------------------------
            DF1N = TDFCND(SMC[k], QUARTZ, SMCMAX, SH2O[k])

            # ----------------------------------------------------------------------
            # CALC THE VERTICAL SOIL TEMP GRADIENT THRU BOTTOM LAYER. 计算穿过底层的垂直土壤温度梯度。
            # ----------------------------------------------------------------------
            DENOM = .5 * (ZSOIL[k - 1] + ZSOIL[k]) - ZBOT
            DTSDZ2 = (STC[k] - TBOT) / DENOM

            # ----------------------------------------------------------------------
            # SET MATRIX COEF, CI TO ZERO IF BOTTOM LAYER. 如果是底层，则将矩阵系数CI设置为零。
            # ----------------------------------------------------------------------
            CI[k] = 0.

            # ----------------------------------------------------------------------
            # IF TEMPERATURE AVERAGING INVOKED (ITAVG=TRUE; ELSE SKIP):  CALCULATE 如果触发温度平均（ITAVG=真；否则跳过）：计算
            # TEMP AT BOTTOM OF LAST LAYER. 最后一层底部的温度。
            # ----------------------------------------------------------------------
            if ITAVG:
                TBK1 = TBND(STC[k], TBOT, ZSOIL, ZBOT, k, NSOIL)
        # endIF
        # ----------------------------------------------------------------------
        # THIS ENDS SPECIAL LOOP FOR BOTTOM LAYER. 这结束了底层的特殊循环。
        # ----------------------------------------------------------------------
        # CALCULATE RHSTS FOR THIS LAYER AFTER CALC'NG A PARTIAL PRODUCT. 计算部分乘积后，计算该层的RHST。
        # ----------------------------------------------------------------------
        DENOM = (ZSOIL[k] - ZSOIL[k - 1]) * HCPCT
        RHSTS[k] = (DF1N * DTSDZ2 - DF1K * DTSDZ) / DENOM
        QTOT = -1.0 * DENOM * RHSTS[k]
        SICE = SMC[k] - SH2O[k]

        if ((SICE > 0.) or (TBK < T0) or
                (STC[k] < T0) or (TBK1 < T0)):
            if ITAVG:
                TAVG = TMPAVG(TAVG, TBK, STC[k], TBK1, ZSOIL, NSOIL, k)
            else:
                TAVG = STC[k]
            # endIF
            TSNSR = SNKSRC(TAVG, SMC[k], SH2O[k], ZSOIL, NSOIL,
                           SMCMAX, PSISAT, BEXP, DT, k, QTOT)
            RHSTS[k] = RHSTS[k] - TSNSR / DENOM
        # endIF

        # ----------------------------------------------------------------------
        # CALC MATRIX COEFS, AI, AND BI FOR THIS LAYER. 计算该层的矩阵系数、AI和BI。
        # ----------------------------------------------------------------------
        AI[k] = - DF1 * DDZ / ((ZSOIL[k - 1] - ZSOIL[k]) * HCPCT)
        BI[k] = -(AI[k] + CI[k])

        # ----------------------------------------------------------------------
        # RESET VALUES OF DF1, DTSDZ, DDZ, AND TBK FOR LOOP TO NEXT SOIL LAYER. 重置回路至下一土层的DF1、DTSDZ、DDZ和TBK值。
        # ----------------------------------------------------------------------
        TBK = TBK1
        DF1K = DF1N
        DTSDZ = DTSDZ2
        DDZ = DDZ2
    # ----------------------------------------------------------------------
    # END SUBROUTINE HRT 结束子例程HRT
    # ----------------------------------------------------------------------
    return


# end

# -- 7. PHYSICS SUBROUTINE ==>  SUBROUTINE HRTICE ---------------------- 物理子程序==>子程序HRTICE

def HRTICE(RHSTS, STC, NSOIL, ZSOIL, YY, ZZ1, DF1, AI, BI, CI):
    # ----------------------------------------------------------------------
    # SUBROUTINE HRTICE 子程序
    # ----------------------------------------------------------------------
    # CALCULATE THE RIGHT HAND SIDE OF THE TIME TENDENCY TERM OF THE SOIL 计算土壤时间趋势项的右侧
    # THERMAL DIFFUSION EQUATION IN THE CASE OF SEA-ICE PACK.  ALSO TO 海冰包情况下的热扩散方程。还
    # COMPUTE (PREPARE) THE MATRIX COEFFICIENTS FOR THE TRI-DIAGONAL MATRIX 计算（准备）三对角矩阵的矩阵系数
    # OF THE IMPLICIT TIME SCHEME. 隐式时间方案。
    # ----------------------------------------------------------------------
    # real, PARAMETER:: TBOT=271.16
    # ----------------------------------------------------------------------
    # SET A NOMINAL UNIVERSAL VALUE OF THE SEA-ICE SPECIFIC HEAT CAPACITY, 设定海水比热容的通用标称值，
    # HCPCT = 1880.0*917.0.
    # ----------------------------------------------------------------------
    TBOT = 271.16
    HCPCT = 1.72396E+6

    # ----------------------------------------------------------------------
    # THE INPUT ARGUMENT DF1 IS A UNIVERSALLY CONSTANT VALUE OF SEA-ICE 输入参数DF1是SEA-ICE的通用常量值
    # THERMAL DIFFUSIVITY, SET IN ROUTINE SNOPAC AS DF1 = 2.2. 热扩散率，设置为常规SNOPAC，DF1=2.2。
    # ----------------------------------------------------------------------
    # SET ICE PACK DEPTH.  USE TBOT AS ICE PACK LOWER BOUNDARY TEMPERATURE 设置冰袋深度。使用TBOT作为冰块下限温度
    # (THAT OF UNFROZEN SEA WATER AT BOTTOM OF SEA ICE PACK).  ASSUME ICE （海冰包底部的未冻结海水）。假设ICE
    # PACK IS OF N=NSOIL LAYERS SPANNING A UNIFORM CONSTANT ICE PACK 包装为N=覆盖均匀恒冰包装的NSOIL层
    # THICKNESS AS DEFINED BY ZSOIL(NSOIL) IN ROUTINE SFLX. 常规SFLX中ZSOIL（NSOIL）定义的厚度 。
    # ----------------------------------------------------------------------
    ZBOT = ZSOIL[NSOIL]

    # ----------------------------------------------------------------------
    # CALC THE MATRIX COEFFICIENTS AI, BI, AND CI FOR THE TOP LAYER 计算顶层的矩阵系数AI、BI和CI
    # ----------------------------------------------------------------------
    DDZ = 1.0 / (-0.5 * ZSOIL[1])
    AI[0] = 0.0
    CI[0] = (DF1 * DDZ) / (ZSOIL[0] * HCPCT)
    BI[0] = -CI[0] + DF1 / (0.5 * ZSOIL[0] * ZSOIL[0] * HCPCT * ZZ1)

    # ----------------------------------------------------------------------
    # CALC THE VERTICAL SOIL TEMP GRADIENT BTWN THE TOP AND 2ND SOIL LAYERS. 计算顶层和第二层土壤的垂直土壤温度梯度。
    # RECALC/ADJUST THE SOIL HEAT FLUX.  USE THE GRADIENT AND FLUX TO CALC 重新计算/调整土壤热流。使用梯度和焊剂计算
    # RHSTS FOR THE TOP SOIL LAYER. 顶部土层右侧。
    # ----------------------------------------------------------------------
    DTSDZ = (STC[0] - STC[1]) / (-0.5 * ZSOIL[1])
    SSOIL = DF1 * (STC[0] - YY) / (0.5 * ZSOIL[0] * ZZ1)
    RHSTS[0] = (DF1 * DTSDZ - SSOIL) / (ZSOIL[0] * HCPCT)

    # ----------------------------------------------------------------------
    # INITIALIZE DDZ2 初始化DDZ2
    # ----------------------------------------------------------------------
    DDZ2 = 0.0

    # ----------------------------------------------------------------------
    # LOOP THRU THE REMAINING SOIL LAYERS, REPEATING THE ABOVE PROCESS 循环穿过剩余土层，重复上述过程
    # ----------------------------------------------------------------------
    for K in range(1, NSOIL):
        if K != NSOIL - 1:
            # ----------------------------------------------------------------------
            # CALC THE VERTICAL SOIL TEMP GRADIENT THRU THIS LAYER. 计算穿过该层的垂直土壤温度梯度。
            # ----------------------------------------------------------------------
            DENOM = 0.5 * (ZSOIL[K - 1] - ZSOIL[K + 1])
            DTSDZ2 = (STC[K] - STC[K + 1]) / DENOM
            # ----------------------------------------------------------------------
            # CALC THE MATRIX COEF, CI, AFTER CALC'NG ITS PARTIAL PRODUCT. 计算其部分乘积后，计算矩阵COEF，CI。
            # ----------------------------------------------------------------------
            DDZ2 = 2. / (ZSOIL[K - 1] - ZSOIL[K + 1])
            CI[K] = -DF1 * DDZ2 / ((ZSOIL[K - 1] - ZSOIL[K]) * HCPCT)
        else:
            # ----------------------------------------------------------------------
            # CALC THE VERTICAL SOIL TEMP GRADIENT THRU THE LOWEST LAYER. 通过最低层计算垂直土壤温度梯度。
            # ----------------------------------------------------------------------
            DTSDZ2 = (STC[K] - TBOT) / (.5 * (ZSOIL[K - 1] + ZSOIL[K]) - ZBOT)
            # ----------------------------------------------------------------------
            # SET MATRIX COEF, CI TO ZERO. 将矩阵COEF，CI设置为零。
            # ----------------------------------------------------------------------
            CI[K] = 0.
        # endIF

        # ----------------------------------------------------------------------
        # CALC RHSTS FOR THIS LAYER AFTER CALC'NG A PARTIAL PRODUCT. 在计算部分产品后，计算该层的RHSTS。
        # ----------------------------------------------------------------------
        DENOM = (ZSOIL[K] - ZSOIL[K - 1]) * HCPCT
        RHSTS[K] = (DF1 * DTSDZ2 - DF1 * DTSDZ) / DENOM

        # ----------------------------------------------------------------------
        # CALC MATRIX COEFS, AI, AND BI FOR THIS LAYER. 计算该层的矩阵系数、AI和BI。
        # ----------------------------------------------------------------------
        AI[K] = - DF1 * DDZ / ((ZSOIL[K - 1] - ZSOIL[K]) * HCPCT)
        BI[K] = -(AI[K] + CI[K])

        # ----------------------------------------------------------------------
        # RESET VALUES OF DTSDZ AND DDZ FOR LOOP TO NEXT SOIL LYR. 重置回路至下一土壤LYR的DTSDZ和DDZ值。
        # ----------------------------------------------------------------------
        DTSDZ = DTSDZ2
        DDZ = DDZ2
    return


# ----------------------------------------------------------------------
# END SUBROUTINE HRTICE 结束子例程HRTICE
# ----------------------------------------------------------------------

# -- 8. PHYSICS SUBROUTINE ==>  SUBROUTINE HSTEP ---------------------- 物理子程序==>子程序HSTEP
def HSTEP(STCOUT, STCIN, RHSTS, DT, NSOIL, AI, BI, CI):
    # ----------------------------------------------------------------------
    # SUBROUTINE HSTEP 子例程HSTEP
    # ----------------------------------------------------------------------
    # CALCULATE/UPDATE THE SOIL TEMPERATURE FIELD. 计算/更新土壤温度场。
    # ----------------------------------------------------------------------
    NSOLD = 20
    import numpy as np
    RHSTSin = np.zeros(NSOIL, dtype=float)
    CIin = np.zeros(NSOLD, dtype=float)
    for K in range(0, NSOIL):
        RHSTS[K] = RHSTS[K] * DT
        AI[K] = AI[K] * DT
        BI[K] = 1. + BI[K] * DT
        CI[K] = CI[K] * DT
    # ----------------------------------------------------------------------
    # COPY VALUES FOR INPUT VARIABLES BEFORE CALL TO ROSR12 调用ROSR12之前复制输入变量的值
    # ----------------------------------------------------------------------
    for K in range(0, NSOIL):
        RHSTSin[K] = RHSTS[K]

    for K in range(0, NSOLD):
        CIin[K] = CI[K]
    # ----------------------------------------------------------------------
    # SOLVE THE TRI-DIAGONAL MATRIX EQUATION 求解三对角矩阵方程
    # ----------------------------------------------------------------------
    ROSR12(CI, AI, BI, CIin, RHSTSin, RHSTS, NSOIL)
    # ----------------------------------------------------------------------
    # CALC/UPDATE THE SOIL TEMPS USING MATRIX SOLUTION 使用矩阵溶液计算/更新土壤温度
    # ----------------------------------------------------------------------
    for K in range(0, NSOIL):
        STCOUT[K] = STCIN[K] + CI[K]

    return STCOUT


# ----------------------------------------------------------------------
# END SUBROUTINE HSTEP 结束子例程HSTEP
# ----------------------------------------------------------------------

# -- 9. PHYSICS SUBROUTINE ==>  SUBROUTINE NOPAC  ---------------------- 物理子程序==>子程序NOPAC

def NOPAC(ETP, ETA, PRCP, SMC, SMCMAX, SMCWLT,
          SMCREF, SMCDRY, CMC, CMCMAX, NSOIL, DT, SHDFAC,
          SBETA, Q2, T1, SFCTMP, T24, TH2, FDOWN, F1, SSOIL,
          STC, EPSCA, BEXP, PC, RCH, RR, CFACTR,
          SH2O, SLOPE, KDT, FRZFACT, PSISAT, ZSOIL,
          DKSAT, DWSAT, TBOT, ZBOT, RUNOFF1, RUNOFF2,
          RUNOFF3, EDIR, EC, ET, ETT, NROOT, ICE, RTDIS,
          QUARTZ, FXEXP, CSOIL,
          BETA, DRIP, DEW, FLX1, FLX2, FLX3):
    # ----------------------------------------------------------------------
    # SUBROUTINE NOPAC 次例行NOPAC
    # ----------------------------------------------------------------------
    # CALCULATE SOIL MOISTURE AND HEAT FLUX VALUES AND UPDATE SOIL MOISTURE 计算土壤水热通量值并更新土壤水分
    # CONTENT AND SOIL HEAT CONTENT VALUES FOR THE CASE WHEN NO SNOW PACK IS 无积雪包情况下的含量和土壤热含量值
    # PRESENT. 现在。
    # ----------------------------------------------------------------------
    import math
    CP = 1004.5
    SIGMA = 5.67E-8

    # ----------------------------------------------------------------------
    # EXECUTABLE CODE BEGINS HERE: 可执行代码从这里开始：
    # CONVERT ETP FROM KG M-2 S-1 TO MS-1 AND INITIALIZE DEW. 将ETP从KG M-2 S-1转换为MS-1并初始化DEW。
    # CONVERT PRCP FROM 'KG M-2 S-1' TO 'M S-1'. 将PRCP从“KG M-2 S-1”转换为“M S-1”
    # ----------------------------------------------------------------------
    PRCP1 = PRCP * 0.001
    ETP1 = ETP * 0.001
    DEW = 0.0
    ETA1 = 0.0
    DF1 = 0.0

    if ETP > 0.0:

        SMFLX(ETA1, SMC, NSOIL, CMC, ETP1, DT, PRCP1, ZSOIL,
              SH2O, SLOPE, KDT, FRZFACT,
              SMCMAX, BEXP, PC, SMCWLT, DKSAT, DWSAT,
              SMCREF, SHDFAC, CMCMAX,
              SMCDRY, CFACTR, RUNOFF1, RUNOFF2, RUNOFF3,
              EDIR, EC, ET, ETT, SFCTMP, Q2, NROOT, RTDIS, FXEXP,
              DRIP)
        # ----------------------------------------------------------------------
        #       CONVERT MODELED EVAPOTRANSPIRATION FM  M S-1  TO  KG M-2 S-1 转换模式蒸发输送FM M S-1至KG M-2 S-1
        # ----------------------------------------------------------------------
        ETA = ETA1 * 1000.0

    # ----------------------------------------------------------------------
    #        EDIR = EDIR1 * 1000.0
    #        EC = EC1 * 1000.0
    #        ETT = ETT1 * 1000.0
    #        ET(1) = ET1(1) * 1000.0
    #        ET(2) = ET1(2) * 1000.0
    #        ET(3) = ET1(3) * 1000.0
    #        ET(4) = ET1(4) * 1000.0
    # ----------------------------------------------------------------------

    else:
        # ----------------------------------------------------------------------
        # IF ETP < 0, ASSUME DEW FORMS (TRANSFORM ETP1 INTO DEW AND REINITIALIZE
        # 如果ETP<0，则假设DEW表（将ETP1转换为DEW并重新初始化
        # ETP1 TO ZERO). ETP1为零）。
        # ----------------------------------------------------------------------
        DEW = -ETP1
        ETP1 = 0.0

        # ----------------------------------------------------------------------
        # CONVERT PRCP FROM 'KG M-2 S-1' TO 'M S-1' AND ADD DEW AMOUNT. 将PRCP从“KG M-2 S-1”转换为“M S-1”，并添加DEW金额。
        # ----------------------------------------------------------------------
        PRCP1 = PRCP1 + DEW

        SMFLX(ETA1, SMC, NSOIL, CMC, ETP1, DT, PRCP1, ZSOIL,
              SH2O, SLOPE, KDT, FRZFACT,
              SMCMAX, BEXP, PC, SMCWLT, DKSAT, DWSAT,
              SMCREF, SHDFAC, CMCMAX,
              SMCDRY, CFACTR, RUNOFF1, RUNOFF2, RUNOFF3,
              EDIR, EC, ET, ETT, SFCTMP, Q2, NROOT, RTDIS, FXEXP,
              DRIP)

        # ----------------------------------------------------------------------
        # CONVERT MODELED EVAPOTRANSPIRATION FROM 'M S-1' TO 'KG M-2 S-1'. 将模式蒸发传输从“M S-1”转换为“KG M-2 S-1”。
        # ----------------------------------------------------------------------
        ETA = ETA1 * 1000.0
        # endIF
    # ----------------------------------------------------------------------
    # BASED ON ETP AND E VALUES, DETERMINE BETA 基于ETP和E值，确定BETA
    # ----------------------------------------------------------------------
    if ETP <= 0.0:
        BETA = 0.0
        if ETP < 0.0:
            BETA = 1.0
            ETA = ETP
    # endIF
    else:
        BETA = ETA / ETP
    # endIF

    # ----------------------------------------------------------------------
    # GET SOIL THERMAL DIFFUXIVITY/CONDUCTIVITY FOR TOP SOIL LYR, 获取表层土壤LYR的土壤热扩散率/电导率，
    # CALC. ADJUSTED TOP LYR SOIL TEMP AND ADJUSTED SOIL FLUX, THEN 计算调整的顶层土壤温度和调整的土壤通量，然后
    # CALL SHFLX TO COMPUTE/UPDATE SOIL HEAT FLUX AND SOIL TEMPS. 致电SHFLX计算/更新土壤热流和土壤温度。
    # ----------------------------------------------------------------------
    DF1 = TDFCND(SMC[0], QUARTZ, SMCMAX, SH2O[0])

    # ----------------------------------------------------------------------
    # VEGETATION GREENNESS FRACTION REDUCTION IN SUBSURFACE HEAT FLUX 地下热流中植被绿度分数的降低
    # VIA REDUCTION FACTOR, WHICH IS CONVENIENT TO APPLY HERE TO THERMAL 通过还原系数，此系数适用于热
    # DIFFUSIVITY THAT IS LATER USED IN HRT TO COMPUTE SUB SFC HEAT FLUX HRT中稍后用于计算亚SFC热流的扩散率
    # (SEE ADDITIONAL COMMENTS ON VEG EFFECT SUB-SFC HEAT FLX IN （请参阅关于VEG EFFECT SUB-SFC热FLX的附加注释
    # ROUTINE SFLX) 常规SFLX）
    # ----------------------------------------------------------------------
    DF1 = DF1 * math.exp(SBETA * SHDFAC)

    # ----------------------------------------------------------------------
    # COMPUTE INTERMEDIATE TERMS PASSED TO ROUTINE HRT (VIA ROUTINE 传递给常规HRT的计算机中间条款（通过常规
    # SHFLX BELOW) FOR USE IN COMPUTING SUBSURFACE HEAT FLUX IN HRT SHFLX下图）用于计算HRT中的地下热流
    # ----------------------------------------------------------------------
    YYNUM = FDOWN - SIGMA * T24
    YY = SFCTMP + (YYNUM / RCH + TH2 - SFCTMP - BETA * EPSCA) / RR
    ZZ1 = DF1 / (-0.5 * ZSOIL[0] * RCH * RR) + 1.0

    SHFLX(SSOIL, STC, SMC, SMCMAX, NSOIL, T1, DT, YY, ZZ1, ZSOIL,
          TBOT, ZBOT, SMCWLT, PSISAT, SH2O, BEXP, F1, DF1, ICE,
          QUARTZ, CSOIL)

    # ----------------------------------------------------------------------
    # SET FLX1 AND FLX3 (SNOPACK PHASE CHANGE HEAT FLUXES) TO ZERO SINCE 将FLX1和FLX3（SNOPACK相变热流）自设置为零
    # THEY ARE NOT USED HERE IN SNOPAC.  FLX2 (FREEZING RAIN HEAT FLUX) WAS SNOPAC中未使用。FLX2（冻雨热流）是
    # SIMILARLY INITIALIZED IN THE PENMAN ROUTINE. 在彭曼常规赛中同样首发。
    # ----------------------------------------------------------------------
    FLX1 = 0.0
    FLX3 = 0.0
    return


# ----------------------------------------------------------------------
# END SUBROUTINE NOPAC  结束子程序NOPAC
# ----------------------------------------------------------------------

# -- 10. PHYSICS SUBROUTINE ==>  SUBROUTINE PENMAN  -------------------- 物理子程序==>子程序彭曼
def PENMAN(SFCTMP, SFCPRS, CH, T2V, TH2, PRCP, FDOWN, T24,
           SSOIL, Q2, Q2SAT, ETP, RCH, EPSCA, RR, SNOWNG, FRZGRA,
           DQSDT2, FLX2):
    # ----------------------------------------------------------------------
    # SUBROUTINE PENMAN  彭曼支线
    # ----------------------------------------------------------------------
    # CALCULATE POTENTIAL EVAPORATION FOR THE CURRENT POINT.  VARIOUS 计算当前点的潜在蒸发量。各种
    # PARTIAL SUMS/PRODUCTS ARE ALSO CALCULATED AND PASSED BACK TO THE 还计算了部分金额/产品，并返回给
    # CALLING ROUTINE FOR LATER USE. 通话常规，供日后使用。
    # ----------------------------------------------------------------------
    CP = 1004.6
    CPH2O = 4.218E+3
    CPICE = 2.106E+3
    R = 287.04
    ELCP = 2.4888E+3
    LSUBF = 3.335E+5
    LSUBC = 2.501000E+6
    SIGMA = 5.67E-8

    # ----------------------------------------------------------------------
    # EXECUTABLE CODE BEGINS HERE: 可执行代码从这里开始：
    # ----------------------------------------------------------------------
    FLX2 = 0.0

    # ----------------------------------------------------------------------
    # PREPARE PARTIAL QUANTITIES FOR PENMAN EQUATION. 为彭曼方程准备部分数量。
    # ----------------------------------------------------------------------
    DELTA = ELCP * DQSDT2
    T24 = SFCTMP * SFCTMP * SFCTMP * SFCTMP
    RR = T24 * 6.48E-8 / (SFCPRS * CH) + 1.0
    RHO = SFCPRS / (R * T2V)
    RCH = RHO * CP * CH

    # ----------------------------------------------------------------------
    # ADJUST THE PARTIAL SUMS / PRODUCTS WITH THE LATENT HEAT 用潜热调整部分金额/产品
    # EFFECTS CAUSED BY FALLING PRECIPITATION. 下降降水引起的影响。
    # ----------------------------------------------------------------------
    if not SNOWNG:
        if PRCP > 0.0:
            RR = RR + CPH2O * PRCP / RCH
    else:
        RR = RR + CPICE * PRCP / RCH
    # endIF

    FNET = FDOWN - SIGMA * T24 - SSOIL

    # ----------------------------------------------------------------------
    # INCLUDE THE LATENT HEAT EFFECTS OF FRZNG RAIN CONVERTING TO ICE ON 包括冻雨转冰的潜热效应
    # IMPACT IN THE CALCULATION OF FLX2 AND FNET. FLX2和FNET计算中的影响。
    # ----------------------------------------------------------------------
    if FRZGRA:
        FLX2 = -LSUBF * PRCP
        FNET = FNET - FLX2
    # endIF

    # ----------------------------------------------------------------------
    # FINISH PENMAN EQUATION CALCULATIONS. 完成彭曼方程计算。
    # ----------------------------------------------------------------------
    RAD = FNET / RCH + TH2 - SFCTMP
    A = ELCP * (Q2SAT - Q2)
    EPSCA = (A * RR + RAD * DELTA) / (DELTA + RR)
    ETP = EPSCA * RCH / LSUBC
    noahder.Ep = ETP
    return (RCH, RR)  # (SNOWNG, FRZGRA)


# ----------------------------------------------------------------------
# END SUBROUTINE PENMAN 结束子程序彭曼
# ----------------------------------------------------------------------

# -- 11. PHYSICS SUBROUTINE ==>  SUBROUTINE ROSR12  --------------------
def ROSR12(P, A, B, C, D, DELTA, NSOIL):
    # ----------------------------------------------------------------------
    # SUBROUTINE ROSR12 ROSR12次例程
    # ----------------------------------------------------------------------
    # INVERT (SOLVE) THE TRI-DIAGONAL MATRIX PROBLEM SHOWN BELOW: 反演（求解）如下所示的三对角矩阵问题：
    # ###                                            ### ###  ###   ###  ###
    # #B(1), C(1),  0  ,  0  ,  0  ,   . . .  ,    0   # #      #   #      #
    # #A(2), B(2), C(2),  0  ,  0  ,   . . .  ,    0   # #      #   #      #
    # # 0  , A(3), B(3), C(3),  0  ,   . . .  ,    0   # #      #   # D(3) #
    # # 0  ,  0  , A(4), B(4), C(4),   . . .  ,    0   # # P(4) #   # D(4) #
    # # 0  ,  0  ,  0  , A(5), B(5),   . . .  ,    0   # # P(5) #   # D(5) #
    # # .                                          .   # #  .   # = #   .  #
    # # .                                          .   # #  .   #   #   .  #
    # # .                                          .   # #  .   #   #   .  #
    # # 0  , . . . , 0 , A(M-2), B(M-2), C(M-2),   0   # #P(M-2)#   #D(M-2)#
    # # 0  , . . . , 0 ,   0   , A(M-1), B(M-1), C(M-1)# #P(M-1)#   #D(M-1)#
    # # 0  , . . . , 0 ,   0   ,   0   ,  A(M) ,  B(M) # # P(M) #   # D(M) #
    # ###                                            ### ###  ###   ###  ###
    # ----------------------------------------------------------------------
    # ----------------------------------------------------------------------
    # INITIALIZE EQN COEF C FOR THE LOWEST SOIL LAYER 初始化最低土层的EQN COEF C
    # ----------------------------------------------------------------------
    C[NSOIL] = 0.0

    # ----------------------------------------------------------------------
    # SOLVE THE COEFS FOR THE 1ST SOIL LAYER 解决第一层土壤的覆土问题
    # ----------------------------------------------------------------------
    P[0] = -C[0] / B[0]
    DELTA[0] = D[0] / B[0]

    # ----------------------------------------------------------------------
    # SOLVE THE COEFS FOR SOIL LAYERS 2 THRU NSOIL 通过NSOIL解决土层2的COEFS
    # ----------------------------------------------------------------------
    for K in range(1, NSOIL):
        P[K] = -C[K] * (1.0 / (B[K] + A[K] * P[K - 1]))
        DELTA[K] = (D[K] - A[K] * DELTA[K - 1]) * (1.0 / (B[K] + A[K] * P[K - 1]))

    # ----------------------------------------------------------------------
    # SET P TO DELTA FOR LOWEST SOIL LAYER 将P设置为最低土层的三角洲
    # ----------------------------------------------------------------------
    P[NSOIL] = DELTA[NSOIL]

    # ----------------------------------------------------------------------
    # ADJUST P FOR SOIL LAYERS 2 THRU NSOIL 通过NSOIL调整土层2的P
    # ----------------------------------------------------------------------
    for K in range(1, NSOIL):
        KK = NSOIL - K + 1
        P[KK] = P[KK] * P[KK + 1] + DELTA[KK]

    return


# ----------------------------------------------------------------------
# END SUBROUTINE ROSR12 结束子例程ROSR12
# ----------------------------------------------------------------------

# -- 12. PHYSICS SUBROUTINE ==>  SUBROUTINE SFCDIF --------------------- 物理子程序==>子程序SFCDIF

def SFCDIF(ZLM, Z0, THZ0, THLM, SFCSPD, CZIL, AKMS, AKHS):
    # ----------------------------------------------------------------------
    # SUBROUTINE SFCDIF 子例程SFCDIF
    # ----------------------------------------------------------------------
    # CALCULATE SURFACE LAYER EXCHANGE COEFFICIENTS VIA ITERATIVE PROCESS. 通过迭代过程计算表层交换系数。
    # SEE CHEN ET AL (1997, BLM) 见CHEN等人（1997，BLM）
    # ----------------------------------------------------------------------

    # real WWST, WWST2, G, VKRM, EXCM, BETA, BTG, ELFC, WOLD, WNEW 真实值
    # real PIHF, EPSU2, EPSUST, EPSIT, EPSA, ZTMIN, ZTMAX, HPBL, SQVISC  真实值
    # real RIC, RRIC, FHNEU, RFC, RFAC, ZZ, PSLMU, PSLMS, PSLHU, PSLHS  真实值
    # real XX, PSPMU, YY, PSPMS, PSPHU, PSPHS, ZLM, Z0, THZ0, THLM  真实值
    # real SFCSPD, CZIL, AKMS, AKHS, ZILFC, ZU, ZT, RDZ, CXCH  真实值
    # real DTHV, DU2, BTGH, WSTAR2, USTAR, ZSLU, ZSLT, RLOGU, RLOGT  真实值
    # real RLMO, ZETALT, ZETALU, ZETAU, ZETAT, XLU4, XLT4, XU4, XT4  真实值
    # real XLU, XLT, XU, XT, PSMZ, SIMM, PSHZ, SIMH, USTARK, RLMN, RLMA  真实值
    # CC   ......REAL ZTFC

    # integer ITRMX, ILECH, ITR
    import math
    WWST = 1.2
    WWST2 = WWST * WWST
    G = 9.8
    VKRM = 0.40
    EXCM = 0.001
    BETA = 1. / 270.
    BTG = BETA * G
    ELFC = VKRM * BTG
    WOLD = .15
    WNEW = 1.0 - WOLD
    ITRMX = 5
    PIHF = 3.14159265 / 2.0
    # ----------------------------------------------------------------------
    EPSU2 = 1.E-4
    EPSUST = 0.07
    EPSIT = 1.E-4
    EPSA = 1.E-8
    ZTMIN = -5.
    ZTMAX = 1.
    HPBL = 1000.0
    SQVISC = 258.2
    # ----------------------------------------------------------------------
    RIC = 0.183
    RRIC = 1.0 / RIC
    FHNEU = 0.8
    RFC = 0.191
    RFAC = RIC / (FHNEU * RFC * RFC)

    # ----------------------------------------------------------------------
    # NOTE: THE TWO CODE BLOCKS BELOW DEFINE FUNCTIONS 注：以下两个代码块定义了功能
    # ----------------------------------------------------------------------
    # LECH'S SURFACE FUNCTIONS 莱赫的表面函数
    # ----------------------------------------------------------------------
    def PSLMU(ZZ):
        result = -0.96 * math.log(1.0 - 4.5 * ZZ)
        return result

    def PSLMS(ZZ):
        result = ZZ * RRIC - 2.076 * (1. - 1. / (ZZ + 1.))
        return result

    def PSLHU(ZZ):
        result = -0.96 * math.log(1.0 - 4.5 * ZZ)
        return result

    def PSLHS(ZZ):
        result = ZZ * RFAC - 2.076 * (1. - 1. / (ZZ + 1.))
        return result

    # ----------------------------------------------------------------------
    # PAULSON'S SURFACE FUNCTIONS 鲍尔森曲面函数
    # ----------------------------------------------------------------------
    def PSPMU(XX):
        result = -2. * math.log((XX + 1.) * 0.5) - math.log((XX * XX + 1.) * 0.5) + 2. * math.atan(XX) \
                 - PIHF
        return result

    def PSPMS(YY):
        result = 5. * YY
        return result

    def PSPHU(XX):
        result = -2. * math.log((XX * XX + 1.) * 0.5)
        return result

    def PSPHS(YY):
        result = 5. * YY
        return result

    # ----------------------------------------------------------------------
    # THIS ROUTINE SFCDIF CAN HANDLE BOTH OVER OPEN WATER (SEA, OCEAN) AND 该常规SFCDIF可在开阔水域（海洋）和
    # OVER SOLID SURFACE (LAND, SEA-ICE). 固体表面上（陆地、海洋）。
    # ----------------------------------------------------------------------
    ILECH = 0

    # ----------------------------------------------------------------------
    #     ZTFC: RATIO OF ZOH/ZOM  LESS OR EQUAL THAN 1 ZTFC：ZOH/ZOM的比率小于或等于1
    #     C......ZTFC=0.1
    #     CZIL: CONSTANT C IN Zilitinkevich, S. S.1995,:NOTE ABOUT ZT
    # ----------------------------------------------------------------------
    ZILFC = -CZIL * VKRM * SQVISC

    # ----------------------------------------------------------------------
    ZU = float(Z0)
    #     C.......ZT=Z0*ZTFC
    RDZ = 1. / ZLM
    CXCH = EXCM * RDZ
    DTHV = THLM - THZ0
    DU2 = max(SFCSPD * SFCSPD, EPSU2)

    # ----------------------------------------------------------------------
    # BELJARS CORRECTION OF USTAR 贝加尔斯修正案
    # ----------------------------------------------------------------------
    BTGH = BTG * HPBL
    # cc   If statements to avoid TANGENT LINEAR problems near zero
    if BTGH * AKHS * DTHV != 0.0:
        WSTAR2 = WWST2 * abs(BTGH * AKHS * DTHV) ** (2. / 3.)
    else:
        WSTAR2 = 0.0
    # endIF
    USTAR = max(math.sqrt(AKMS * math.sqrt(DU2 + WSTAR2)), EPSUST)

    # ----------------------------------------------------------------------
    # ZILITINKEVITCH APPROACH FOR ZT ZT的ZILITINKEVITCH方法
    # ----------------------------------------------------------------------
    ZT = math.exp(ZILFC * math.sqrt(USTAR * float(Z0))) * float(Z0)

    # ----------------------------------------------------------------------
    ZSLU = ZLM + ZU
    ZSLT = ZLM + ZT
    #     PRINT*,'ZSLT=',ZSLT
    #     PRINT*,'ZLM=',ZLM
    #     PRINT*,'ZT=',ZT
    #
    RLOGU = math.log(ZSLU / ZU)
    RLOGT = math.log(ZSLT / ZT)
    #
    RLMO = ELFC * AKHS * DTHV / math.pow(USTAR, 3)
    #     PRINT*,'RLMO=',RLMO
    #     PRINT*,'ELFC=',ELFC
    #     PRINT*,'AKHS=',AKHS
    #     PRINT*,'DTHV=',DTHV
    #     PRINT*,'USTAR=',USTAR

    for ITR in range(0, ITRMX):
        # ----------------------------------------------------------------------
        # 1./MONIN-OBUKKHOV LENGTH-SCALE 1./MONIN-OBUKHOV长度-比例
        # ----------------------------------------------------------------------
        ZETALT = max(ZSLT * RLMO, ZTMIN)
        RLMO = ZETALT / ZSLT
        ZETALU = ZSLU * RLMO
        ZETAU = ZU * RLMO
        ZETAT = ZT * RLMO

        if ILECH == 0:
            if RLMO < 0.:
                XLU4 = 1. - 16. * ZETALU
                XLT4 = 1. - 16. * ZETALT
                XU4 = 1. - 16. * ZETAU
                XT4 = 1. - 16. * ZETAT

                XLU = math.sqrt(math.sqrt(XLU4))
                XLT = math.sqrt(math.sqrt(XLT4))
                XU = math.sqrt(math.sqrt(XU4))
                XT = math.sqrt(math.sqrt(XT4))

                PSMZ = PSPMU(XU)
                #     PRINT*,'-----------1------------'
                #     PRINT*,'PSMZ=',PSMZ
                #     PRINT*,'PSPMU(ZETAU)=',PSPMU(ZETAU)
                #     PRINT*,'XU=',XU
                #     PRINT*,'------------------------'
                SIMM = PSPMU(XLU) - PSMZ + RLOGU
                PSHZ = PSPHU(XT)
                SIMH = PSPHU(XLT) - PSHZ + RLOGT
            else:
                ZETALU = min(ZETALU, ZTMAX)
                ZETALT = min(ZETALT, ZTMAX)
                PSMZ = PSPMS(ZETAU)
                #     PRINT*,'-----------2------------'
                #     PRINT*,'PSMZ=',PSMZ
                #     PRINT*,'PSPMS(ZETAU)=',PSPMS(ZETAU)
                #     PRINT*,'ZETAU=',ZETAU
                #     PRINT*,'------------------------'
                SIMM = PSPMS(ZETALU) - PSMZ + RLOGU
                PSHZ = PSPHS(ZETAT)
                SIMH = PSPHS(ZETALT) - PSHZ + RLOGT
                # endIF
        else:
            # ----------------------------------------------------------------------
            # LECH'S FUNCTIONS LECH的功能
            # ----------------------------------------------------------------------
            if RLMO < 0.:
                PSMZ = PSLMU(ZETAU)
                #     PRINT*,'-----------3------------'
                #     PRINT*,'PSMZ=',PSMZ
                #     PRINT*,'PSLMU(ZETAU)=',PSLMU(ZETAU)
                #     PRINT*,'ZETAU=',ZETAU
                #     PRINT*,'------------------------'
                SIMM = PSLMU(ZETALU) - PSMZ + RLOGU
                PSHZ = PSLHU(ZETAT)
                SIMH = PSLHU(ZETALT) - PSHZ + RLOGT
            else:
                ZETALU = MIN(ZETALU, ZTMAX)
                ZETALT = MIN(ZETALT, ZTMAX)
                #
                PSMZ = PSLMS(ZETAU)
                #     PRINT*,'-----------4------------'
                #     PRINT*,'PSMZ=',PSMZ
                #     PRINT*,'PSLMS(ZETAU)=',PSLMS(ZETAU)
                #     PRINT*,'ZETAU=',ZETAU
                #     PRINT*,'------------------------'
                SIMM = PSLMS(ZETALU) - PSMZ + RLOGU
                PSHZ = PSLHS(ZETAT)
                SIMH = PSLHS(ZETALT) - PSHZ + RLOGT
        # endIF
        # endIF
        # ----------------------------------------------------------------------
        # BELJAARS CORRECTION FOR USTAR BELJAARS USTAR修正
        # ----------------------------------------------------------------------
        USTAR = max(math.sqrt(AKMS * math.sqrt(DU2 + WSTAR2)), EPSUST)

        # ----------------------------------------------------------------------
        # ZILITINKEVITCH FIX FOR ZT ZT的ZILITINKEVITCH固定
        # ----------------------------------------------------------------------
        ZT = math.exp(ZILFC * math.sqrt(USTAR * float(Z0))) * float(Z0)

        ZSLT = ZLM + ZT
        RLOGT = math.log(ZSLT / ZT)
        # -----------------------------------------------------------------------
        USTARK = USTAR * VKRM
        AKMS = max(USTARK / SIMM, CXCH)
        AKHS = max(USTARK / SIMH, CXCH)
        # -----------------------------------------------------------------------
        # IF STATEMENTS TO AVOID TANGENT LINEAR PROBLEMS NEAR ZERO 避免接近零的切线线性问题的IF语句
        # -----------------------------------------------------------------------
        if BTGH * AKHS * DTHV != 0.0:
            WSTAR2 = WWST2 * math.pow(abs(BTGH * AKHS * DTHV), (2. / 3.))
        else:
            WSTAR2 = 0.0
        # endIF
        RLMN = ELFC * AKHS * DTHV / math.pow(USTAR, 3)
        # -----------------------------------------------------------------------
        RLMA = RLMO * WOLD + RLMN * WNEW
        # -----------------------------------------------------------------------
        #     IF(ABS((RLMN-RLMO)/RLMA).LT.EPSIT)    GO TO 110 如果（ABS（（RLMN-RLMO）/RLMA）.LT.EPSIT）转到110
        # -----------------------------------------------------------------------
        RLMO = RLMA

    return


# ----------------------------------------------------------------------
# END SUBROUTINE SFCDIF 结束子例程SFCDIF
# ----------------------------------------------------------------------

# -- 13. PHYSICS SUBROUTINE ==>  SUBROUTINE SHFLX  --------------------- 物理子例程==>子例程SHFLX

def SHFLX(SSOIL, STC, SMC, SMCMAX, NSOIL, T1, DT, YY, ZZ1, ZSOIL,
          TBOT, ZBOT, SMCWLT, PSISAT, SH2O, BEXP, F1, DF1, ICE,
          QUARTZ, CSOIL):
    # ----------------------------------------------------------------------
    # SUBROUTINE SHFLX 子例程SHFLX
    # ----------------------------------------------------------------------
    # UPDATE THE TEMPERATURE STATE OF THE SOIL COLUMN BASED ON THE THERMAL 基于热学原理更新土柱的温度状态
    # DIFFUSION EQUATION AND UPDATE THE FROZEN SOIL MOISTURE CONTENT BASED 基于扩散方程和冻土含水量的更新
    # ON THE TEMPERATURE. 基于温度。
    # ----------------------------------------------------------------------
    # integer, PARAMETER:: NSOLD=20
    import numpy as np
    T0 = 273.15
    NSOLD = 20
    RHSTS = np.zeros(NSOLD, dtype=float)
    AI = np.zeros(NSOLD, dtype=float)
    BI = np.zeros(NSOLD, dtype=float)
    CI = np.zeros(NSOLD, dtype=float)
    STCF = np.zeros(NSOLD, dtype=float)

    # ----------------------------------------------------------------------
    # HRT ROUTINE CALCS THE RIGHT HAND SIDE OF THE SOIL TEMP DIF EQN HRT常规计算土壤温度差右侧
    # ----------------------------------------------------------------------
    if ICE == 1:
        # ----------------------------------------------------------------------
        # SEA-ICE CASE SEA-ICE外壳
        # ----------------------------------------------------------------------
        HRTICE(RHSTS, STC, NSOIL, ZSOIL, YY, ZZ1, DF1, AI, BI, CI)
        HSTEP(STCF, STC, RHSTS, DT, NSOIL, AI, BI, CI)
    else:
        # ----------------------------------------------------------------------
        # LAND-MASS CASE 土地质量状况
        # ----------------------------------------------------------------------
        HRT(RHSTS, STC, SMC, SMCMAX, NSOIL, ZSOIL, YY, ZZ1, TBOT,
            ZBOT, PSISAT, SH2O, DT,
            BEXP, F1, DF1, QUARTZ, CSOIL, AI, BI, CI)

        STCF = HSTEP(STCF, STC, RHSTS, DT, NSOIL, AI, BI, CI)

    # endIF

    for I in range(0, NSOIL):
        STC[I] = STCF[I]

    # ----------------------------------------------------------------------
    # IN THE NO SNOWPACK CASE (VIA ROUTINE NOPAC BRANCH,) UPDATE THE GRND 在无雪包的情况下（通过NOPAC常规分行）更新地面
    # (SKIN) TEMPERATURE HERE IN RESPONSE TO THE UPDATED SOIL TEMPERATURE 此处的（皮肤）温度响应更新的土壤温度
    # PROFILE ABOVE.  (NOTE: INSPECTION OF ROUTINE SNOPAC SHOWS THAT T1 上述简介。（注：常规SNOPAC检查显示T1
    # BELOW IS A DUMMY VARIABLE ONLY, AS SKIN TEMPERATURE IS UPDATED 以下仅为假人变量，因为皮肤温度已更新
    # DIFFERENTLY IN ROUTINE SNOPAC) 不同于常规SNOPAC）
    # ----------------------------------------------------------------------
    T1 = (YY + (ZZ1 - 1.0) * STC[0]) / ZZ1

    # ----------------------------------------------------------------------
    # CALCULATE SURFACE SOIL HEAT FLUX 计算地表土壤热流
    # ----------------------------------------------------------------------
    SSOIL = DF1 * (STC[0] - T1) / (0.5 * ZSOIL[0])
    del RHSTS
    return


# ----------------------------------------------------------------------
# END SUBROUTINE SHFLX 结束子例程SHFLX
# ----------------------------------------------------------------------

# -- 14. PHYSICS SUBROUTINE ==>  SUBROUTINE SMFLX ---------------------- 物理子例程==>子例程SMFLX
def SMFLX(ETA1, SMC, NSOIL, CMC, ETP1, DT, PRCP1, ZSOIL,
          SH2O, SLOPE, KDT, FRZFACT,
          SMCMAX, BEXP, PC, SMCWLT, DKSAT, DWSAT,
          SMCREF, SHDFAC, CMCMAX,
          SMCDRY, CFACTR, RUNOFF1, RUNOFF2, RUNOFF3,
          EDIR, EC, ET, ETT, SFCTMP, Q2, NROOT, RTDIS, FXEXP,
          DRIP):
    # ----------------------------------------------------------------------
    # SUBROUTINE SMFLX 亚常规SMFLX
    # ----------------------------------------------------------------------
    # CALCULATE SOIL MOISTURE FLUX.  THE SOIL MOISTURE CONTENT (SMC - A PER 计算土壤水分通量。土壤含水量（SMC-A PER
    # UNIT VOLUME MEASUREMENT) IS A DEPENDENT VARIABLE THAT IS UPDATED WITH 单位体积测量）是一个相关变量，用更新
    # PROGNOSTIC EQNS. THE CANOPY MOISTURE CONTENT (CMC) IS ALSO UPDATED. 预测情商。顶篷含水量（CMC）也已更新。
    # FROZEN GROUND VERSION:  NEW STATES ADDED: SH2O, AND FROZEN GROUND 冻结地面版本：新增州：SH2O和冻结地面
    # CORRECTION FACTOR, FRZFACT AND PARAMETER SLOPE. 修正系数、摩擦系数和参数斜率。
    # ----------------------------------------------------------------------
    # integer, PARAMETER:: NSOLD=20
    TFREEZ = 273.15
    # ----------------------------------------------------------------------
    # EXECUTABLE CODE BEGINS HERE. 可执行代码从这里开始。
    # ----------------------------------------------------------------------
    import numpy as np
    DUMMY = 0.
    EDIR = 0.
    NSOLD = 20

    EC = 0.
    ETT = 0.
    ET = np.zeros(NSOIL, dtype=float)
    AI = np.zeros(NSOLD, dtype=float)
    BI = np.zeros(NSOLD, dtype=float)
    CI = np.zeros(NSOLD, dtype=float)
    SICE = np.zeros(NSOLD, dtype=float)
    RHSTT = np.zeros(NSOLD, dtype=float)
    SH2OFG = np.zeros(NSOLD, dtype=float)
    SH2OA = np.zeros(NSOLD, dtype=float)

    for K in range(0, NSOIL):
        ET[K] = 0.
    for K in range(0, NSOLD):
        AI[K] = 0.
        BI[K] = 0.
        CI[K] = 0.
        SICE[K] = 0.
        RHSTT[K] = 0.

    if ETP1 > 0.0:
        # ----------------------------------------------------------------------
        # RETRIEVE DIRECT EVAPORATION FROM SOIL SURFACE.CALL THIS FUNCTION 从土壤表面获取直接蒸发。调用此函数
        # ONLY IF VEG COVER NOT COMPLETE. 仅当VEG盖未完成时。
        # FROZEN GROUND VERSION: SH2O STATES REPLACE SMC STATES. 冰冻地面版本：SH2O状态取代SMC状态。
        # ----------------------------------------------------------------------
        if SHDFAC < 1.:
            EDIR = DEVAP(ETP1, SH2O[0], ZSOIL[0], SHDFAC, SMCMAX, BEXP, DKSAT, DWSAT, SMCDRY, SMCREF, SMCWLT, FXEXP)

        # ----------------------------------------------------------------------
        # INITIALIZE PLANT TOTAL TRANSPIRATION, RETRIEVE PLANT TRANSPIRATION, 初始化植物总蒸腾量，检索植物蒸腾量，
        # AND ACCUMULATE IT FOR ALL SOIL LAYERS. 并对所有土层进行累积。
        # ----------------------------------------------------------------------
        #     PRINT *, '2514 SMFLX', SHDFAC
        if SHDFAC > 0.0:
            TRANSP(ET, NSOIL, ETP1, SH2O, CMC, ZSOIL, SHDFAC, SMCWLT,
                   CMCMAX, PC, CFACTR, SMCREF, SFCTMP, Q2, NROOT, RTDIS)

            for K in range(0, NSOIL):
                ETT = ETT + ET
            # ----------------------------------------------------------------------
            # ! CALCULATE CANOPY EVAPORATION. 计算顶篷蒸发量。
            # ! IF   STATEMENTS    TO   AVOID  TANGENT   LINEAR   PROBLEMS   NEAR   CMC = 0.0.
            # 避免接近CMC的切线线性问题的IF语句=0.0。
            # ! ----------------------------------------------------------------------
            if CMC > 0.0:
                EC = SHDFAC * ((CMC / CMCMAX) ** CFACTR) * ETP1
            else:
                EC = 0.0

            noahder.CMC = CMC
            noahder.n = CFACTR
            noahder.Bc = ETP1
            # ! ----------------------------------------------------------------------
            # ! EC   SHOULD   BE    LIMITED    BY    THE      TOTAL      AMOUNT       OF       AVAILABLE
            # EC应受到可用总金额的限制
            # WATER        ON        THE 水在上面
            # ! CANOPY. - F.CHEN, 18 - OCT - 1994
            CMC2MS = CMC / DT
            EC = min(CMC2MS, EC)

    # ----------------------------------------------------------------------
    # TOTALUPEVAPAND TRANSp TYPES TO OBTAIN ACTUAL EVAPOTRANSP 获得实际蒸发量的总量和运输类型
    # ----------------------------------------------------------------------
    ETA1 = EDIR + ETT + EC
    # ----------------------------------------------------------------------
    # COMPUTE THE RIGHT HAND SIDE OF THE CANOPY EQN TERM ( RHSCT ) 计算顶篷右侧EQN术语（RHSCT）
    # ----------------------------------------------------------------------
    RHSCT = SHDFAC * PRCP1 - EC

    # ----------------------------------------------------------------------
    # CONVERT RHSCT (A RATE) TO TRHSCT (AN AMOUNT) AND ADD IT TO EXISTING
    # 将RHSCT（费率）转换为TRHSCT（金额），并将其添加到现有金额中
    # CMC.  IF RESULTING AMT EXCEEDS MAX CAPACITY, IT BECOMES DRIP AND WILL
    # CMC。如果结果AMT超过最大容量，它将滴水并将
    # FALL TO THE GRND. 摔倒在地上。
    # ----------------------------------------------------------------------
    DRIP = 0
    TRHSCT = DT * RHSCT
    EXCESS = CMC + TRHSCT
    if EXCESS > CMCMAX:
        DRIP = EXCESS - CMCMAX

    # ----------------------------------------------------------------------
    # PCPDRP IS THE COMBINED PRCP1 AND DRIP (FROM CMC) THAT GOES INTO THE PCPDRP是进入
    # SOIL 土壤
    # ----------------------------------------------------------------------
    PCPDRP = (1. - SHDFAC) * PRCP1 + DRIP / DT

    # ----------------------------------------------------------------------
    # STORE ICE CONTENT AT EACH SOIL LAYER BEFORE CALLING SRT & SSTEP 在致电SRT和SSTEP之前，在每层土壤中存储冰含量
    # ----------------------------------------------------------------------
    for I in range(0, NSOIL):
        SICE[I] = SMC[I] - SH2O[I]

    # ----------------------------------------------------------------------
    # CALL SUBROUTINES SRT AND SSTEP TO SOLVE THE SOIL MOISTURE  SRT和SSTEP解决土壤水分问题
    # TENDENCY EQUATIONS. 趋势方程。
    # IF THE INFILTRATING PRECIP RATE IS NONTRIVIAL,  如果渗入降水率是非零的，
    # (WE CONSIDER NONTRIVIAL TO BE A PRECIP TOTAL OVER THE TIME STEP （我们认为非民事行为是时间步的一个预兆
    # EXCEEDING ONE ONE-THOUSANDTH OF THE WATER HOLDING CAPACITY OF 超过
    # THE FIRST SOIL LAYER) 第一层土壤）
    # THEN CALL THE SRT/SSTEP SUBROUTINE PAIR TWICE IN THE MANNER OF 然后以以下方式调用SRT/SSTEP子例程对两次
    # TIME SCHEME "F" (IMPLICIT STATE, AVERAGED COEFFICIENT) 时间方案“F”（隐式状态，平均系数）
    # OF SECTION 2 OF KALNAY AND KANAMITSU (1988, MWR, VOL 116, 《卡尔奈和卡纳米特苏第2部分》（1988年，MWR，第116卷）
    # PAGES 1945-1958)TO MINIMIZE 2-DELTA-T OSCILLATIONS IN THE 第1945-1958页），以最小化
    # SOIL MOISTURE VALUE OF THE TOP SOIL LAYER THAT CAN ARISE BECAUSE 可能产生的表层土壤水分值
    # OF THE EXTREME NONLINEAR DEPENDENCE OF THE SOIL HYDRAULIC 土壤水力的极端非线性依赖
    # DIFFUSIVITY COEFFICIENT AND THE HYDRAULIC CONDUCTIVITY ON THE 扩散系数和导水率
    # SOIL MOISTURE STATE 土壤水分状况
    # OTHERWISE CALL THE SRT/SSTEP SUBROUTINE PAIR ONCE IN THE MANNER OF 否则，以以下方式调用SRT/SSTEP子例程对一次
    # TIME SCHEME "D" (IMPLICIT STATE, EXPLICIT COEFFICIENT) 时间方案“D”（隐式状态，显式系数）
    # OF SECTION 2 OF KALNAY AND KANAMITSU KALNAY和KANAMITSU第2节
    # PCPDRP IS UNITS OF KG/M**2/S OR MM/S, ZSOIL IS NEGATIVE DEPTH IN M PCPDRP单位为KG/M**2/S或MM/S，ZSOIL为负深度，单位为M
    # ----------------------------------------------------------------------
    #      IF ( PCPDRP .GT. 0.0 ) THEN
    # IF((PCPDRP * DT) > (0.001 * 1000.0 * (-ZSOIL(1)) * SMCMAX)):
    if 0 > 1:
        # ----------------------------------------------------------------------
        # FROZEN GROUND VERSION: 冰冻地面版本：
        # SMC STATES REPLACED BY SH2O STATES IN SRT SUBR.  SH2O & SICE STATES
        # SRT SUBR中的SMC状态替换为SH2O状态。SH2O和SICE州
        # INCLUDED IN SSTEP SUBR.  FROZEN GROUND CORRECTION FACTOR, FRZFACT 包含在SSTEP子合同中。冻结地面修正系数，冻结
        # ADDED.  ALL WATER BALANCE CALCULATIONS USING UNFROZEN WATER 补充。使用未冷冻水的所有水平衡计算
        # ----------------------------------------------------------------------
        SRT(RHSTT, EDIR1, ET1, SH2O, SH2O, NSOIL, PCPDRP, ZSOIL,
            DWSAT, DKSAT, SMCMAX, BEXP, RUNOFF1,
            RUNOFF2, DT, SMCWLT, SLOPE, KDT, FRZFACT, SICE, AI, BI, CI)
        SSTEP(SH2OFG, SH2O, DUMMY, RHSTT, RHSCT, DT, NSOIL, SMCMAX,
              CMCMAX, RUNOFF3, ZSOIL, SMC, SICE, AI, BI, CI)

        for K in range(0, NSOIL):
            SH2OA[K] = (SH2O(K) + SH2OFG(K)) * 0.5

        SRT(RHSTT, EDIR1, ET1, SH2O, SH2OA, NSOIL, PCPDRP, ZSOIL,
            DWSAT, DKSAT, SMCMAX, BEXP, RUNOFF1,
            RUNOFF2, DT, SMCWLT, SLOPE, KDT, FRZFACT, SICE, AI, BI, CI)

        SSTEP(SH2O, SH2O, CMC, RHSTT, RHSCT, DT, NSOIL, SMCMAX,
              CMCMAX, RUNOFF3, ZSOIL, SMC, SICE, AI, BI, CI)
    else:
        SRT(RHSTT, EDIR, ET, SH2O, SH2O, NSOIL, PCPDRP, ZSOIL,
            DWSAT, DKSAT, SMCMAX, BEXP, RUNOFF1,
            RUNOFF2, DT, SMCWLT, SLOPE, KDT, FRZFACT, SICE, AI, BI, CI)
        SSTEP(SH2O, SH2O, CMC, RHSTT, RHSCT, DT, NSOIL, SMCMAX,
              CMCMAX, RUNOFF3, ZSOIL, SMC, SICE, AI, BI, CI)

        for K in range(0, NSOIL):
            SH2OA[K] = (SH2O[K] + SH2OFG[K]) * 0.5

        SRT(RHSTT, EDIR, ET, SH2O, SH2OA, NSOIL, PCPDRP, ZSOIL,
            DWSAT, DKSAT, SMCMAX, BEXP, RUNOFF1,
            RUNOFF2, DT, SMCWLT, SLOPE, KDT, FRZFACT, SICE, AI, BI, CI)

        SSTEP(SH2O, SH2O, CMC, RHSTT, RHSCT, DT, NSOIL, SMCMAX,
              CMCMAX, RUNOFF3, ZSOIL, SMC, SICE, AI, BI, CI)

    del ET
    del AI
    del BI
    del CI
    return


# ----------------------------------------------------------------------
# END SUBROUTINE SMFLX 结束子例程SMFLX
# ----------------------------------------------------------------------

# -- 15. PHYSICS SUBROUTINE ==>  SUBROUTINE SNFRAC  ------------------- 物理子程序==>子程序SNFRAC
def SNFRAC(SNEQV, SNUP, SALP, SNOWH, SNCOVR):
    # ----------------------------------------------------------------------
    # SUBROUTINE SNFRAC 子例程 SNFRAC
    # ----------------------------------------------------------------------
    import math
    # CALCULATE SNOW FRACTION (0 -> 1) 计算雪分率（0->1）
    # SNEQV   SNOW WATER EQUIVALENT (M) SNEQV 雪水当量（M）
    # SNUP    THRESHOLD SNEQV DEPTH ABOVE WHICH SNCOVR=1 SNCOVR=1以上SNEQV深度的SNUP阈值
    # SALP    TUNING PARAMETER SALP调谐参数
    # SNCOVR  FRACTIONAL SNOW COVER SNCOVR部分积雪
    # ----------------------------------------------------------------------
    # real SNEQV, SNUP, SALP, SNCOVR, RSNOW, Z0N, SNOWH
    # ----------------------------------------------------------------------
    # SNUP IS VEG-CLASS DEPENDENT SNOWDEPTH THRESHHOLD (SET IN ROUTINE  SNUP是VEG-CLASS相关的雪深阈值（按常规设置
    # REDPRM) ABOVE WHICH SNOCVR=1. 红色），SNOCVR=1。
    # ----------------------------------------------------------------------
    if SNEQV < SNUP:
        RSNOW = SNEQV / SNUP
        SNCOVR = 1. - (math.exp(-SALP * RSNOW) - RSNOW * math.exp(-SALP))
    else:
        SNCOVR = 1.0
    # endIF

    Z0N = 0.035
    # FORMULATION OF DICKINSON ET AL. 1986  DICKINSON等人，1986年

    # SNCOVR=SNOWH/(SNOWH + 5*Z0N)  SNCOVR=雪/（雪+5*Z0N）

    # FORMULATION OF MARSHALL ET AL. 1994  MARSHALL等人1994年的公式
    # SNCOVR=SNEQV/(SNEQV + 2*Z0N)
    return SNCOVR


# ----------------------------------------------------------------------
# END SUBROUTINE SNFRAC 结束子程序SNFRAC
# ----------------------------------------------------------------------
def FRH2O(tkelv1, SMC, SH2O, SMCMAX, BEXP, PSIS):
    #
    # INPUT:
    #   TKELV.........Temperature (Kelvin)  温度（开尔文）
    #   SMC...........Total soil moisture content (volumetric) 土壤总含水量（体积）
    #   SH2O..........Liquid soil moisture content (volumetric) 液态土壤含水量（体积）
    #   SMCMAX........Saturation soil moisture content (from REDPRM) 饱和土壤含水量（来自REDPRM）
    #   B.............Soil type "B" parameter (from REDPRM) 土壤类型“B”参数（来自REDPRM）
    #   PSIS..........Saturated soil matric potential (from REDPRM) 饱和土壤基质电位（来自REDPRM）
    #
    # OUTPUT:
    #   FRH2O.........supercooled liquid water content. 过冷液态水含量。
    # --------------------------------------------------------------------
    import math
    CK = 8.0
    #      PARAMETER (CK=0.0) 参数 (CK=0.0)
    BLIM = 5.5
    #      PARAMETER (BLIM=7.0) 参数 (BLIM=7.0)
    ERROR = 0.005
    HLICE = 3.335E5
    GS = 9.81
    #       PARAMETER (DICE=920.0) 参数 (DICE=920.0)
    #       PARAMETER (DH2O=1000.0) 参数 (DH2O=1000.0)
    T0 = 273.15

    BX = BEXP
    if BEXP > BLIM:
        BX = BLIM

    # INITIALIZING ITERATIONS COUNTER AND ITERATIVE SOLUTION FLAG. 初始化迭代计数器和迭代解决方案标记。
    NLOG = 0
    KCOUNT = 0

    # ----------------------------------------------------------------------
    #  IF TEMPERATURE NOT SIGNIFICANTLY BELOW FREEZING (T0), SH2O = SMC 如果温度不明显低于冻结温度（T0），则SH2O=SMC
    # ----------------------------------------------------------------------

    if tkelv1 > T0 - 0.001:
        FRH2O = SMC
    else:

        if CK != 0.0:

            # -------------------------------------------------------------
            # OPTION 1: ITERATED SOLUTION FOR NONZERO CK 选项1：非零CK的迭代解
            # IN KOREN ET AL, JGR, 1999, EQN 17 IN KOREN等人，JGR，1999，EQN 17
            # -------------------------------------------------------------
            # INITIAL GUESS FOR SWL (frozen content) SWL的初始提示（冻结内容）
            SWL = SMC - SH2O

            # KEEP WITHIN BOUNDS.
            if SWL > (SMC - 0.02):
                SWL = SMC - 0.02
            if SWL < 0.:
                SWL = 0.
            # --------------------------------------------------------------
            #  START OF ITERATIONS 迭代的开始
            # --------------------------------------------------------------

            while NLOG < 10 and KCOUNT == 0:
                NLOG = NLOG + 1

                # if(SMC-SWL  ==  0  or  TKELV .eq. 0 .or. TKELV-T0  >  0) :
                # call printerror() 调用printerror（）
                # endif
                # print*,'8535',SMCMAX, SMC, SWL,BX, TKELV,T0,TKELV
                # DF = ALOG((PSIS * GS / HLICE) * ((1. + CK * SWL) ** 2.) * &
                #          (SMCMAX / (SMC - SWL)) ** BX) - ALOG(-(TKELV - T0) / TKELV)
                DF = math.log((PSIS * GS / HLICE) * (math.pow((1. + CK * SWL), 2.)) * \
                              math.pow(SMCMAX / (SMC - SWL), BX)) - math.log(-(tkelv1 - T0) / tkelv1)
                DENOM = 2. * CK / (1. + CK * SWL) + BX / (SMC - SWL)
                SWLK = SWL - DF / DENOM
                # BOUNDS USEFUL FOR MATHEMATICAL SOLUTION. 对数学解有用的边界。
                if SWLK > (SMC - 0.02):
                    SWLK = SMC - 0.02
                if SWLK < 0.:
                    SWLK = 0.
                # MATHEMATICAL SOLUTION BOUNDS APPLIED. 应用的数学解边界。
                DSWL = abs(SWLK - SWL)
                SWL = SWLK

                # ---------------------------------------------------------------
                # IF MORE THAN 10 ITERATIONS, USE EXPLICIT METHOD (CK=0 APPROX.)
                # 如果迭代次数超过10次，则使用显式方法（CK=0近似值）
                # WHEN DSWL LESS OR EQ. ERROR, NO MORE ITERATIONS REQUIRED.
                # 当DSWL小于或等于错误时，不再需要重复。
                # ---------------------------------------------------------------
                if DSWL <= ERROR:
                    KCOUNT = KCOUNT + 1
            # end if

            # print*,'9561'
            # ---------------------------------------------------------------
            #  END OF ITERATIONS 迭代结束
        # ---------------------------------------------------------------
        # BOUNDS APPLIED WITHIN DO-BLOCK ARE VALID FOR PHYSICAL SOLUTION. DO-BLOCK中应用的边界对物理解有效。
        FRH2O = SMC - SWL

    # ---------------------- END OPTION 1 -------------------

    # endIF

    # -------------------------------------------------------------
    # OPTION 2: EXPLICIT SOLUTION FOR FLERCHINGER EQ. i.e. CK=0 选项2：FLERCHINGER EQ的显式解决方案，即CK=0
    # IN KOREN ET AL., JGR, 1999, EQN 17 N KOREN等人，JGR，1999，EQN 17
    # ----------------------------------------------------------------
    if KCOUNT == 0:
        #  Print*,'Flerchinger used in NEW version. Iterations=',NLOG
        # FK = (((HLICE / (GS * (-PSIS))) * ((TKELV - T0) / TKELV)) **(-1 / BX)) * SMCMAX
        FK = (math.pow(((HLICE / (GS * (-PSIS))) * ((tkelv1 - T0) / tkelv1)), (-1 / BX))) * SMCMAX
        #  APPLY PHYSICAL BOUNDS TO FLERCHINGER SOLUTION 将物理界应用于FLERCHINGER溶液
        if FK < 0.02:
            FK = 0.02
        FRH2O = min(FK, SMC)
    return FRH2O


# ------------------- END  FRH2O ---------------------

def SNKSRC(TAVG, SMC, SH2O, ZSOIL, NSOIL,
           SMCMAX, PSISAT, BEXP, DT, K, QTOT):
    # ----------------------------------------------------------------------
    # FUNCTION SNKSRC SNKSRC功能
    # ----------------------------------------------------------------------
    # CALCULATE SINK/SOURCE TERM OF THE TERMAL DIFFUSION EQUATION. (SH2O) IS 计算扩散方程项的汇/源项。（SH2O）信息系统
    # AVAILABLE LIQUED WATER. 可用的液态水。
    # ----------------------------------------------------------------------

    DH2O = 1.0000E3
    HLICE = 3.3350E5
    T0 = 2.7315E2

    if K == 1:
        DZ = -ZSOIL[0]
    else:
        DZ = ZSOIL[K - 1] - ZSOIL[K]
    # endIF

    # ----------------------------------------------------------------------
    # VIA FUNCTION FRH2O, COMPUTE POTENTIAL OR 'EQUILIBRIUM' UNFROZEN  通过功能FRH2O、计算机电位或“平衡”解冻
    # SUPERCOOLED FREE WATER FOR GIVEN SOIL TYPE AND SOIL LAYER TEMPERATURE. 给定土壤类型和土壤层温度的过冷自由水。
    # FUNCTION FRH20 INVOKES EQN (17) FROM V. KOREN ET AL (1999, JGR, VOL. 功能FRH20从V.KOREN等人（1999年，JGR，卷。
    # 104, PG 19573).  (ASIDE:  LATTER EQN IN JOURNAL IN CENTIGRADE UNITS.
    # 104页，PG 19573页）。（旁白：以百分位为单位的日记中的后半部分。
    # ROUTINE FRH2O USE FORM OF EQN IN KELVIN UNITS.) KELVIN装置中EQN的常规FRH2O使用形式。）
    # ----------------------------------------------------------------------
    FREE = FRH2O(TAVG, SMC, SH2O, SMCMAX, BEXP, PSISAT)

    # ----------------------------------------------------------------------
    # IN NEXT BLOCK OF CODE, INVOKE EQN 18 OF V. KOREN ET AL (1999, JGR, 在下一个代码块中，请发票V.KOREN等人（1999，JGR，
    # VOL. 104, PG 19573.)  THAT IS, FIRST ESTIMATE THE NEW AMOUNTOF LIQUID 第104卷，PG 19573），即首先估计新的液体量
    # WATER, 'XH2O', IMPLIED BY THE SUM OF (1) THE LIQUID WATER AT THE BEGIN 水，“XH2O”，由（1）初始液态水的总和表示
    # OF CURRENT TIME STEP, AND (2) THE FREEZE OF THAW CHANGE IN LIQUID 当前时间步，以及（2）液体中THAW变化冻结
    # WATER IMPLIED BY THE HEAT FLUX 'QTOT' PASSED IN FROM ROUTINE HRT. 常规HRT传入的热流“QTOT”所暗示的水。
    # SECOND, DETERMINE IF XH2O NEEDS TO BE BOUNDED BY 'FREE' (EQUIL AMT) OR 第二，确定XH2O是否需要“自由”（等于AMT）或
    # IF 'FREE' NEEDS TO BE BOUNDED BY XH2O. 如果“免费”，则需要受到XH2O的约束。
    # ----------------------------------------------------------------------
    XH2O = SH2O + QTOT * DT / (DH2O * HLICE * DZ)

    # ----------------------------------------------------------------------
    # FIRST, IF FREEZING AND REMAINING LIQUID LESS THAN LOWER BOUND, THEN 首先，如果冷冻和剩余液体小于下限，则
    # REDUCE EXTENT OF FREEZING, THEREBY LETTING SOME OR ALL OF HEAT FLUX 减少冻结程度，从而释放部分或全部热流
    # QTOT COOL THE SOIL TEMP LATER IN ROUTINE HRT. 用常规HRT冷却土壤温度。
    # ----------------------------------------------------------------------
    if XH2O < SH2O and XH2O < FREE:
        if FREE > SH2O:
            XH2O = SH2O
        else:
            XH2O = FREE
    # endIF
    # endIF

    # ----------------------------------------------------------------------
    # SECOND, IF THAWING AND THE INCREASE IN LIQUID WATER GREATER THAN UPPER 第二，如果冷却和液态水中的增加大于上限
    # BOUND, THEN REDUCE EXTENT OF THAW, THEREBY LETTING SOME OR ALL OF HEAT 绑定，然后减少THAW的范围，从而释放部分或全部热量
    # FLUX QTOT WARM THE SOIL TEMP LATER IN ROUTINE HRT. FLUX QTOT在常规HRT中加热土壤温度。
    # ----------------------------------------------------------------------
    if XH2O > SH2O and XH2O > FREE:
        if FREE < SH2O:
            XH2O = SH2O
        else:
            XH2O = FREE
    # endIF
    # endIF

    if XH2O < 0.:
        XH2O = 0.
    if XH2O > SMC:
        XH2O = SMC

    # ----------------------------------------------------------------------
    # CALCULATE PHASE-CHANGE HEAT SOURCE/SINK TERM FOR USE IN ROUTINE HRT 计算常规HRT中使用的相变热源/散热器术语
    # AND UPDATE LIQUID WATER TO REFLCET FINAL FREEZE/THAW INCREMENT. 并更新液态水，以重新设置最终冷冻/冷却增量。
    # ----------------------------------------------------------------------
    SNKSRC1 = -DH2O * HLICE * DZ * (XH2O - SH2O) / DT
    SH2O = XH2O

    # ----------------------------------------------------------------------
    # END FUNCTION SNKSRC 结束功能SNKSRC
    # ----------------------------------------------------------------------
    return SNKSRC1


# end

# -- 16. PHYSICS SUBROUTINE ==>  SUBROUTINE SNOPAC --------------------- 物理子程序==>子程序SNOPAC

def SNOPAC(ETP, ETA, PRCP, PRCP1, SNOWNG, SMC, SMCMAX, SMCWLT,
           SMCREF, SMCDRY, CMC, CMCMAX, NSOIL, DT,
           SBETA, DF1, Q2, T1, SFCTMP,
           T24, TH2, FDOWN, F1, SSOIL, STC, EPSCA, SFCPRS,
           BEXP, PC, RCH, RR, CFACTR, SNCOVR, ESD, SNDENS,
           SNOWH, SH2O, SLOPE, KDT, FRZFACT, PSISAT, SNUP,
           ZSOIL, DWSAT, DKSAT, TBOT, ZBOT, SHDFAC, RUNOFF1,
           RUNOFF2, RUNOFF3, EDIR, EC, ET, ETT, NROOT, SNOMLT,
           ICE, RTDIS, QUARTZ, FXEXP, CSOIL,
           BETA, DRIP, DEW, FLX1, FLX2, FLX3, ESNOW):
    # ----------------------------------------------------------------------
    # SUBROUTINE SNOPAC 亚常规SNOPAC
    # ----------------------------------------------------------------------
    # CALCULATE SOIL MOISTURE AND HEAT FLUX VALUES & UPDATE SOIL MOISTURE 亚常规计算土壤水分和热流值并更新土壤水分
    # CONTENT AND SOIL HEAT CONTENT VALUES FOR THE CASE WHEN A SNOW PACK IS  箱子的含量和土壤热含量值当积雪包
    # PRESENT. 出席。
    # ----------------------------------------------------------------------
    # LOGICAL * 1 SNOWNG

    CP = 1004.5
    CPH2O = 4.218E+3
    CPICE = 2.106E+3
    ESDMIN = 1.E-6
    LSUBF = 3.335E+5
    LSUBC = 2.501000E+6
    LSUBS = 2.83E+6
    SIGMA = 5.67E-8
    TFREEZ = 273.15
    ETA1 = 0.0
    SSOIL1 = 0.0
    # PARAMETER(SNOEXP = 1.0)                                      <
    SNOEXP = 2.0

    # ----------------------------------------------------------------------
    # EXECUTABLE CODE BEGINS HERE: 可执行代码从这里开始：
    # CONVERT POTENTIAL EVAP (ETP) FROM KG M-2 S-1 TO M S-1 AND THEN TO AN
    # 将潜在EVAP（ETP）从KG M-2 S-1转换为M S-1，然后转换为
    # AMOUNT (M) GIVEN TIMESTEP (DT) AND CALL IT AN EFFECTIVE SNOWPACK 金额（M）给出时间（DT）并称之为有效雪包
    # REDUCTION AMOUNT, ETP2 (M).  THIS IS THE AMOUNT THE SNOWPACK WOULD BE 减少金额，ETP2（M）。这是积雪包的金额
    # REDUCED DUE TO EVAPORATION FROM THE SNOW SFC DURING THE TIMESTEP. 由于雪地SFC在最短时间内蒸发而减少。
    # EVAPORATION WILL PROCEED AT THE POTENTIAL RATE UNLESS THE SNOW DEPTH 除非雪深，否则蒸发将以潜在速率进行
    # IS LESS THAN THE EXPECTED SNOWPACK REDUCTION. 小于预期雪包减少量。
    # IF SEAICE (ICE=1), BETA REMAINS=1. 如果SEAICE（ICE=1），BETA REMAINS=1。
    # ----------------------------------------------------------------------
    PRCP1 = PRCP1 * 0.001

    ETP2 = ETP * 0.001 * DT
    BETA = 1.0
    if ICE != 1:
        if ESD < ETP2:
            BETA = ESD / ETP2
    # ----------------------------------------------------------------------
    # IF ETP<0 (DOWNWARD) THEN DEWFALL (=FROSTFALL IN THIS CASE).
    # ----------------------------------------------------------------------
    DEW = 0.0
    if ETP < 0.0:
        DEW = -ETP * 0.001

    # ----------------------------------------------------------------------
    # IF PRECIP IS FALLING, CALCULATE HEAT FLUX FROM SNOW SFC TO NEWLY
    # ACCUMULATING PRECIP.  NOTE THAT THIS REFLECTS THE FLUX APPROPRIATE FOR
    # THE NOT-YET-UPDATED SKIN TEMPERATURE (T1).  ASSUMES TEMPERATURE OF THE
    # SNOWFALL STRIKING THE GOUND IS =SFCTMP (LOWEST MODEL LEVEL AIR TEMP).
    # ----------------------------------------------------------------------
    FLX1 = 0.0
    if SNOWNG:
        FLX1 = CPICE * PRCP * (T1 - SFCTMP)
    else:
        if PRCP > 0.0:
            FLX1 = CPH2O * PRCP * (T1 - SFCTMP)
    # endIF
    # ----------------------------------------------------------------------
    # CALCULATE AN 'EFFECTIVE SNOW-GRND SFC TEMP' (T12) BASED ON HEAT FLUXES
    # BETWEEN THE SNOW PACK AND THE SOIL AND ON NET RADIATION.
    # INCLUDE FLX1 (PRECIP-SNOW SFC) AND FLX2 (FREEZING RAIN LATENT HEAT)
    # FLUXES.
    # FLX2 REFLECTS FREEZING RAIN LATENT HEAT FLUX USING T1 CALCULATED IN
    # PENMAN.
    # ----------------------------------------------------------------------
    DSOIL = -(0.5 * ZSOIL[0])
    DTOT = SNOWH + DSOIL
    DENOM = 1.0 + DF1 / (DTOT * RR * RCH)
    #      T12A = ( (FDOWN-FLX1-FLX2-SIGMA*T24)/RCH &
    #            + TH2 - SFCTMP - BETA*EPSCA ) / RR
    T12A = ((FDOWN - FLX1 - FLX2 - SIGMA * T24) / RCH
            + TH2 - SFCTMP - BETA * EPSCA) / RR
    T12B = DF1 * STC[1] / (DTOT * RR * RCH)
    T12 = (SFCTMP + T12A + T12B) / DENOM

    # ----------------------------------------------------------------------
    # IF THE 'EFFECTIVE SNOW-GRND SFC TEMP' IS AT OR BELOW FREEZING, NO SNOW
    # MELT WILL OCCUR.  SET THE SKIN TEMP TO THIS EFFECTIVE TEMP.  REDUCE
    # (BY SUBLIMINATION ) OR INCREASE (BY FROST) THE DEPTH OF THE SNOWPACK,
    # DEPENDING ON SIGN OF ETP.
    # UPDATE SOIL HEAT FLUX (SSOIL) USING NEW SKIN TEMPERATURE (T1)
    # SINCE NO SNOWMELT, SET ACCUMULATED SNOWMELT TO ZERO, SET 'EFFECTIVE'
    # PRECIP FROM SNOWMELT TO ZERO, SET PHASE-CHANGE HEAT FLUX FROM SNOWMELT
    # TO ZERO.
    # ----------------------------------------------------------------------
    if T12 <= TFREEZ:
        T1 = T12
        SSOIL = DF1 * (T1 - STC[0]) / DTOT
        # ESD = MAX(0.0, ESD-ETP2)
        ESD = max(0.0, ESD - ETP2)
        FLX3 = 0.0
        EX = 0.0
        SNOMLT = 0.0
    else:
        # ----------------------------------------------------------------------
        # IF THE 'EFFECTIVE SNOW-GRND SFC TEMP' IS ABOVE FREEZING, SNOW MELT
        # WILL OCCUR.  CALL THE SNOW MELT RATE,EX AND AMT, SNOMLT.  REVISE THE
        # EFFECTIVE SNOW DEPTH.  REVISE THE SKIN TEMP BECAUSE IT WOULD HAVE CHGD
        # DUE TO THE LATENT HEAT RELEASED BY THE MELTING. CALC THE LATENT HEAT
        # RELEASED, FLX3. SET THE EFFECTIVE PRECIP, PRCP1 TO THE SNOW MELT RATE,
        # EX FOR USE IN SMFLX.  ADJUSTMENT TO T1 TO ACCOUNT FOR SNOW PATCHES.
        # CALCULATE QSAT VALID AT FREEZING POINT.  NOTE THAT ESAT (SATURATION
        # VAPOR PRESSURE) VALUE OF 6.11E+2 USED HERE IS THAT VALID AT FRZZING
        # POINT.  NOTE THAT ETP FROM CALL PENMAN IN SFLX IS IGNORED HERE IN
        # FAVOR OF BULK ETP OVER 'OPEN WATER' AT FREEZING TEMP.
        # UPDATE SOIL HEAT FLUX (S) USING NEW SKIN TEMPERATURE (T1)
        # ----------------------------------------------------------------------
        #        T1 = TFREEZ * SNCOVR + T12 * (1.0 - SNCOVR)
        # Change made by MEK - Feb 2004
        #  Non-linear weighting of snow vs non-snow covered portions of gridbox
        #  so with SNOEXP = 2.0 (>1), surface skin temperature is higher than for
        #   the linear case (SNOEXP = 1).
        T1 = TFREEZ * SNCOVR ** SNOEXP + T12 * (1.0 - SNCOVR ** SNOEXP)

        QSAT = (0.622 * 6.11E2) / (SFCPRS - 0.378 * 6.11E2)
        ETP = RCH * (QSAT - Q2) / CP
        ETP2 = ETP * 0.001 * DT
        BETA = 1.0
        SSOIL = DF1 * (T1 - STC[0]) / DTOT

        # ----------------------------------------------------------------------
        # IF POTENTIAL EVAP (SUBLIMATION) GREATER THAN DEPTH OF SNOWPACK.
        # BETA<1
        # SNOWPACK HAS SUBLIMATED AWAY, SET DEPTH TO ZERO.
        # ----------------------------------------------------------------------
        if ESD < ETP2:
            #        IF (ESD .LE. ESNOW2) THEN
            # IF(ESD - ESNOW2 <= ESDMIN):
            BETA = ESD / ETP2
            ESD = 0.0
            EX = 0.0
            SNOMLT = 0.0
        else:
            # ----------------------------------------------------------------------
            # POTENTIAL EVAP (SUBLIMATION) LESS THAN DEPTH OF SNOWPACK, RETAIN
            #   BETA=1.
            # SNOWPACK (ESD) REDUCED BY POTENTIAL EVAP RATE
            # ETP3 (CONVERT TO FLUX)
            # ----------------------------------------------------------------------
            ESD = ESD - ETP2
            # ESD = ESD - ESNOW2
            ETP3 = ETP * LSUBC
            SEH = RCH * (T1 - TH2)
            T14 = T1 * T1
            T14 = T14 * T14
            FLX3 = FDOWN - FLX1 - FLX2 - SIGMA * T14 - SSOIL - SEH - ETP3
            # FLX3 = FDOWN - FLX1 - FLX2 - SIGMA * T14 - SSOIL - SEH - ETANRG
            if FLX3 <= 0.0:
                FLX3 = 0.0
            EX = FLX3 * 0.001 / LSUBF

            # ----------------------------------------------------------------------
            # SNOWMELT REDUCTION DEPENDING ON SNOW COVER
            # IF SNOW COVER LESS THAN 5% NO SNOWMELT REDUCTION
            # ***NOTE:  DOES 'IF' BELOW FAIL TO MATCH THE MELT WATER WITH THE MELT
            #           ENERGY?
            # ----------------------------------------------------------------------
            if SNCOVR > 0.05:
                EX = EX * SNCOVR
            SNOMLT = EX * DT

            # ----------------------------------------------------------------------
            # ESDMIN REPRESENTS A SNOWPACK DEPTH THRESHOLD VALUE BELOW WHICH WE
            # CHOOSE NOT TO RETAIN ANY SNOWPACK, AND INSTEAD INCLUDE IT IN SNOWMELT.
            # ----------------------------------------------------------------------
            if ESD - SNOMLT >= ESDMIN:
                ESD = ESD - SNOMLT
            else:
                # ----------------------------------------------------------------------
                # SNOWMELT EXCEEDS SNOW DEPTH
                # ----------------------------------------------------------------------
                EX = ESD / DT
                FLX3 = EX * 1000.0 * LSUBF
                SNOMLT = ESD
                ESD = 0.0
        # ----------------------------------------------------------------------
        # END OF 'ESD .LE. ETP2' IF-BLOCK
        # ----------------------------------------------------------------------
        PRCP1 = PRCP1 + EX
    # ----------------------------------------------------------------------
    # END OF 'T12 .LE. TFREEZ' IF-BLOCK
    # ----------------------------------------------------------------------
    # endIF
    # ----------------------------------------------------------------------
    # FINAL BETA NOW IN HAND, SO COMPUTE EVAPORATION.  EVAP EQUALS ETP
    # UNLESS BETA<1.
    # ----------------------------------------------------------------------
    ETA = BETA * ETP
    # ----------------------------------------------------------------------
    # SET THE EFFECTIVE POTNL EVAPOTRANSP (ETP1) TO ZERO SINCE THIS IS SNOW
    # CASE, SO SURFACE EVAP NOT CALCULATED FROM EDIR, EC, OR ETT IN SMFLX
    # (BELOW).
    # IF SEAICE (ICE=1) SKIP CALL TO SMFLX.
    # SMFLX RETURNS UPDATED SOIL MOISTURE VALUES.  IN THIS, THE SNOW PACK
    # CASE, ETA1 IS NOT USED IN CALCULATION OF EVAP.
    # ----------------------------------------------------------------------
    ETP1 = 0.0

    if ICE != 1:
        SMFLX(ETA1, SMC, NSOIL, CMC, ETP1, DT, PRCP1, ZSOIL,
              SH2O, SLOPE, KDT, FRZFACT,
              SMCMAX, BEXP, PC, SMCWLT, DKSAT, DWSAT,
              SMCREF, SHDFAC, CMCMAX,
              SMCDRY, CFACTR, RUNOFF1, RUNOFF2, RUNOFF3,
              EDIR, EC, ET, ETT, SFCTMP, Q2, NROOT, RTDIS, FXEXP,
              DRIP)
    # ----------------------------------------------------------------------
    # BEFORE CALL SHFLX IN THIS SNOWPACK CASE, SET ZZ1 AND YY ARGUMENTS TO
    # SPECIAL VALUES THAT ENSURE THAT GROUND HEAT FLUX CALCULATED IN SHFLX
    # MATCHES THAT ALREADY COMPUTER FOR BELOW THE SNOWPACK, THUS THE SFC
    # HEAT FLUX TO BE COMPUTED IN SHFLX WILL EFFECTIVELY BE THE FLUX AT THE
    # SNOW TOP SURFACE.  T11 IS A DUMMY ARGUEMENT SO WE WILL NOT USE THE
    # SKIN TEMP VALUE AS REVISED BY SHFLX.
    # ----------------------------------------------------------------------
    ZZ1 = 1.0
    YY = STC[0] - 0.5 * SSOIL * ZSOIL[0] * ZZ1 / DF1
    T11 = T1

    # ----------------------------------------------------------------------
    # SHFLX WILL CALC/UPDATE THE SOIL TEMPS.  NOTE:  THE SUB-SFC HEAT FLUX
    # (SSOIL1) AND THE SKIN TEMP (T11) OUTPUT FROM THIS SHFLX CALL ARE NOT
    # USED  IN ANY SUBSEQUENT CALCULATIONS. RATHER, THEY ARE DUMMY VARIABLES
    # HERE IN THE SNOPAC CASE, SINCE THE SKIN TEMP AND SUB-SFC HEAT FLUX ARE
    # UPDATED INSTEAD NEAR THE BEGINNING OF THE CALL TO SNOPAC.
    # ----------------------------------------------------------------------
    SHFLX(SSOIL1, STC, SMC, SMCMAX, NSOIL, T11, DT, YY, ZZ1, ZSOIL,
          TBOT, ZBOT, SMCWLT, PSISAT, SH2O, BEXP, F1, DF1, ICE,
          QUARTZ, CSOIL)

    # ----------------------------------------------------------------------
    # SNOW DEPTH AND DENSITY ADJUSTMENT BASED ON SNOW COMPACTION.  YY IS
    # ASSUMED TO BE THE SOIL TEMPERTURE AT THE TOP OF THE SOIL COLUMN.
    # ----------------------------------------------------------------------
    if ESD > 0.:
        SNOWPACK(ESD, DT, SNOWH, SNDENS, T1, YY)
    else:
        ESD = 0.
        SNOWH = 0.
        SNDENS = 0.
        SNCOND = 1.
        SNCOVR = 0.
    # endIF
    return


# ----------------------------------------------------------------------
# END SUBROUTINE SNOPAC
# ----------------------------------------------------------------------

# -- 17. PHYSICS SUBROUTINE ==>  SUBROUTINE SNOWPACK -------------------
def SNOWPACK(ESD, DTSEC, SNOWH, SNDENS, TSNOW, TSOIL):
    # ----------------------------------------------------------------------
    # SUBROUTINE SNOWPACK
    # ----------------------------------------------------------------------
    import math
    # CALCULATE COMPACTION OF SNOWPACK UNDER CONDITIONS OF INCREASING SNOW
    # DENSITY, AS OBTAINED FROM AN APPROXIMATE SOLUTION OF E. ANDERSON'S
    # DIFFERENTIAL EQUATION (3.29), NOAA TECHNICAL REPORT NWS 19, BY VICTOR
    # KOREN, 03/25/95.
    # ----------------------------------------------------------------------
    # ESD     WATER EQUIVALENT OF SNOW (M)
    # DTSEC   TIME STEP (SEC)
    # SNOWH   SNOW DEPTH (M)
    # SNDENS  SNOW DENSITY (G/CM3=DIMENSIONLESS FRACTION OF H2O DENSITY)
    # TSNOW   SNOW SURFACE TEMPERATURE (K)
    # TSOIL   SOIL SURFACE TEMPERATURE (K)
    #
    # SUBROUTINE WILL RETURN NEW VALUES OF SNOWH AND SNDENS
    # ----------------------------------------------------------------------

    C1 = 0.01
    C2 = 21.0
    G = 9.81
    KN = 4000.0

    # ----------------------------------------------------------------------
    # CONVERSION INTO SIMULATION UNITS
    # ----------------------------------------------------------------------
    SNOWHC = SNOWH * 100.
    ESDC = ESD * 100.
    DTHR = DTSEC / 3600.
    TSNOWC = TSNOW - 273.15
    TSOILC = TSOIL - 273.15

    # ----------------------------------------------------------------------
    # CALCULATING OF AVERAGE TEMPERATURE OF SNOW PACK
    # ----------------------------------------------------------------------
    TAVGC = 0.5 * (TSNOWC + TSOILC)

    # ----------------------------------------------------------------------
    # CALCULATING OF SNOW DEPTH AND DENSITY AS A RESULT OF COMPACTION
    #  SNDENS=DS0*(EXP(BFAC*ESD)-1.)/(BFAC*ESD)
    #  BFAC=DTHR*C1*EXP(0.08*TAVGC-C2*DS0)
    # NOTE: BFAC*ESD IN SNDENS EQN ABOVE HAS TO BE CAREFULLY TREATED
    # NUMERICALLY BELOW:
    #   C1 IS THE FRACTIONAL INCREASE IN DENSITY (1/(CM*HR))
    #   C2 IS A CONSTANT (CM3/G) KOJIMA ESTIMATED AS 21 CMS/G
    # ----------------------------------------------------------------------
    if ESDC > 1.E-2:
        ESDCX = ESDC
    else:
        ESDCX = 1.E-2
    # endIF
    BFAC = DTHR * C1 * math.exp(0.08 * TAVGC - C2 * SNDENS)

    #      DSX = SNDENS*((DEXP(BFAC*ESDC)-1.)/(BFAC*ESDC))
    # ----------------------------------------------------------------------
    # THE FUNCTION OF THE FORM (e**x-1)/x IMBEDDED IN ABOVE EXPRESSION
    # FOR DSX WAS CAUSING NUMERICAL DIFFICULTIES WHEN THE DENOMINATOR "x"
    # (I.E. BFAC*ESDC) BECAME ZERO OR APPROACHED ZERO (DESPITE THE FACT THAT
    # THE ANALYTICAL FUNCTION (e**x-1)/x HAS A WELL DEFINED LIMIT AS
    # "x" APPROACHES ZERO), HENCE BELOW WE REPLACE THE (e**x-1)/x
    # EXPRESSION WITH AN EQUIVALENT, NUMERICALLY WELL-BEHAVED
    # POLYNOMIAL EXPANSION.
    #
    # NUMBER OF TERMS OF POLYNOMIAL EXPANSION, AND HENCE ITS ACCURACY,
    # IS GOVERNED BY ITERATION LIMIT "IPOL".
    #      IPOL GREATER THAN 9 ONLY MAKES A DIFFERENCE ON DOUBLE
    #            PRECISION (RELATIVE ERRORS GIVEN IN PERCENT %).
    #       IPOL=9, FOR REL.ERROR <~ 1.6 E-6 % (8 SIGNIFICANT DIGITS)
    #       IPOL=8, FOR REL.ERROR <~ 1.8 E-5 % (7 SIGNIFICANT DIGITS)
    #       IPOL=7, FOR REL.ERROR <~ 1.8 E-4 % ...
    # ----------------------------------------------------------------------
    IPOL = 4
    PEXP = 0.
    for J in range(IPOL, 0, -1):
        #        PEXP = (1. + PEXP)*BFAC*ESDC/REAL(J+1)
        PEXP = (1. + PEXP) * BFAC * ESDCX / float(J + 1)

    PEXP = PEXP + 1.
    DSX = SNDENS * PEXP
    # ----------------------------------------------------------------------
    # ABOVE LINE ENDS POLYNOMIAL SUBSTITUTION
    # ----------------------------------------------------------------------
    #     END OF KOREAN FORMULATION

    #     BASE FORMULATION (COGLEY ET AL., 1990)
    #     CONVERT DENSITY FROM G/CM3 TO KG/M3
    #       DSM=SNDENS*1000.0

    #       DSX=DSM+DTSEC*0.5*DSM*G*ESD/ &
    #          (1E7*EXP(-0.02*DSM+KN/(TAVGC+273.16)-14.643))

    #     CONVERT DENSITY FROM KG/M3 TO G/CM3
    #       DSX=DSX/1000.0

    #     END OF COGLEY ET AL. FORMULATION

    # ----------------------------------------------------------------------
    # SET UPPER/LOWER LIMIT ON SNOW DENSITY
    # ----------------------------------------------------------------------
    if DSX > 0.40:
        DSX = 0.40
    if DSX < 0.05:
        DSX = 0.05
    SNDENS = DSX
    # ----------------------------------------------------------------------
    # UPDATE OF SNOW DEPTH AND DENSITY DEPENDING ON LIQUID WATER DURING
    # SNOWMELT.  ASSUMED THAT 13% OF LIQUID WATER CAN BE STORED IN SNOW PER
    # DAY DURING SNOWMELT TILL SNOW DENSITY 0.40.
    # ----------------------------------------------------------------------
    if TSNOWC >= 0.:
        DW = 0.13 * DTHR / 24.
        SNDENS = SNDENS * (1. - DW) + DW
        if SNDENS > 0.40:
            SNDENS = 0.40
    # endIF

    # ----------------------------------------------------------------------
    # CALCULATE SNOW DEPTH (CM) FROM SNOW WATER EQUIVALENT AND SNOW DENSITY.
    # CHANGE SNOW DEPTH UNITS TO METERS
    # ----------------------------------------------------------------------
    SNOWHC = ESDC / SNDENS
    SNOWH = SNOWHC * 0.01

    # ----------------------------------------------------------------------
    # END SUBROUTINE SNOWPACK
    # ----------------------------------------------------------------------
    return


# -- 18. PHYSICS SUBROUTINE ==>  SUBROUTINE SNOWZ0 ---------------------

def SNOWZ0(SNCOVR, Z0):
    # ----------------------------------------------------------------------
    # SUBROUTINE SNOWZ0
    # ----------------------------------------------------------------------
    # CALCULATE TOTAL ROUGHNESS LENGTH OVER SNOW
    # SNCOVR  FRACTIONAL SNOW COVER
    # Z0      ROUGHNESS LENGTH (m)
    # Z0S     SNOW ROUGHNESS LENGTH:=0.001 (m)
    # ----------------------------------------------------------------------
    # real SNCOVR, Z0, Z0S
    #      PARAMETER (Z0S=0.001)

    # CURRENT NOAH LSM CONDITION - MBEK, 09-OCT-2001
    Z0S = Z0
    Z0 = (1 - SNCOVR) * Z0 + SNCOVR * Z0S
    return


# ----------------------------------------------------------------------
# END SUBROUTINE SNOWZ0
# ----------------------------------------------------------------------

# -- 19. PHYSICS SUBROUTINE ==>  SUBROUTINE SNOW_NEW -------------------
def SNOW_NEW(TEMP, NEWSN, SNOWH, SNDENS):
    # SUBROUTINE SNOW_NEW
    # ----------------------------------------------------------------------
    # CALCULATE SNOW DEPTH AND DENSITITY TO ACCOUNT FOR THE NEW SNOWFALL.
    # NEW VALUES OF SNOW DEPTH & DENSITY RETURNED.
    #
    # TEMP    AIR TEMPERATURE (K)
    # NEWSN   NEW SNOWFALL (M)
    # SNOWH   SNOW DEPTH (M)
    # SNDENS  SNOW DENSITY (G/CM3=DIMENSIONLESS FRACTION OF H2O DENSITY)
    # CONVERSION INTO SIMULATION UNITS
    # ----------------------------------------------------------------------
    SNOWHC = SNOWH * 100.
    NEWSNC = NEWSN * 100.
    TEMPC = TEMP - 273.15

    # ----------------------------------------------------------------------
    # CALCULATING NEW SNOWFALL DENSITY DEPENDING ON TEMPERATURE
    # EQUATION FROM GOTTLIB L. 'A GENERAL RUNOFF MODEL FOR SNOWCOVERED
    # AND GLACIERIZED BASIN', 6TH NORDIC HYDROLOGICAL CONFERENCE,
    # VEMADOLEN, SWEDEN, 1980, 172-177PP.
    # -----------------------------------------------------------------------
    if TEMPC <= -15.:
        DSNEW = 0.05
    else:
        DSNEW = 0.05 + 0.0017 * (TEMPC + 15.) ** 1.5
    # endIF

    # ----------------------------------------------------------------------
    # ADJUSTMENT OF SNOW DENSITY DEPENDING ON NEW SNOWFALL
    # ----------------------------------------------------------------------
    HNEWC = NEWSNC / DSNEW
    SNDENS = (SNOWHC * SNDENS + HNEWC * DSNEW) / (SNOWHC + HNEWC)
    SNOWHC = SNOWHC + HNEWC
    SNOWH = SNOWHC * 0.01

    return SNDENS


# ----------------------------------------------------------------------
# END SUBROUTINE SNOW_NEW
# ----------------------------------------------------------------------

# -- 20. PHYSICS SUBROUTINE ==>  SUBROUTINE SRT - ----------------------

def SRT(RHSTT, EDIR, ET, SH2O, SH2OA, NSOIL, PCPDRP,
        ZSOIL, DWSAT, DKSAT, SMCMAX, BEXP, RUNOFF1,
        RUNOFF2, DT, SMCWLT, SLOPE, KDT, FRZX, SICE, AI, BI, CI):
    # ----------------------------------------------------------------------
    # SUBROUTINE SRT
    # ----------------------------------------------------------------------
    # CALCULATE THE RIGHT HAND SIDE OF THE TIME TENDENCY TERM OF THE SOIL
    # WATER DIFFUSION EQUATION.  ALSO TO COMPUTE ( PREPARE ) THE MATRIX
    # COEFFICIENTS FOR THE TRI-DIAGONAL MATRIX OF THE IMPLICIT TIME SCHEME.
    # ----------------------------------------------------------------------
    # ----------------------------------------------------------------------
    # FROZEN GROUND VERSION:
    # REFERENCE FROZEN GROUND PARAMETER, CVFRZ, IS A SHAPE PARAMETER OF
    # AREAL DISTRIBUTION FUNCTION OF SOIL ICE CONTENT WHICH EQUALS 1/CV.
    # CV IS A COEFFICIENT OF SPATIAL VARIATION OF SOIL ICE CONTENT.  BASED
    # ON FIELD DATA CV DEPENDS ON AREAL MEAN OF FROZEN DEPTH, AND IT CLOSE
    # TO CONSTANT = 0.6 IF AREAL MEAN FROZEN DEPTH IS ABOVE 20 CM.  THAT IS
    # WHY PARAMETER CVFRZ = 3 (INT{1/0.6*0.6}).
    # CURRENT LOGIC DOESN'T ALLOW CVFRZ BE BIGGER THAN 3
    # ----------------------------------------------------------------------
    import numpy as np
    import math
    CVFRZ = 3

    # ----------------------------------------------------------------------
    # DETERMINE RAINFALL INFILTRATION RATE AND RUNOFF.  INCLUDE THE
    # INFILTRATION FORMULE FROM SCHAAKE AND KOREN MODEL.
    # MODIFIED BY Q DUAN
    # ----------------------------------------------------------------------
    IOHINF = 1
    NSOLD = 20
    WDF = 0.
    WCND, WCND2, WDF2 = 0.0, 0.0, 0.0
    DMAX = np.zeros(NSOLD, dtype=float)
    # ----------------------------------------------------------------------
    # LET SICEMAX BE THE GREATEST, IF ANY, FROZEN WATER CONTENT WITHIN SOIL
    # LAYERS.
    # ----------------------------------------------------------------------
    SICEMAX = 0.0
    for KS in range(0, NSOIL):
        if SICE[KS] > SICEMAX:
            SICEMAX = SICE[KS]
    # ----------------------------------------------------------------------
    # DETERMINE RAINFALL INFILTRATION RATE AND RUNOFF
    # ----------------------------------------------------------------------
    PDDUM = PCPDRP
    RUNOFF1 = 0.0
    if PCPDRP != 0.0:
        # ----------------------------------------------------------------------
        # MODIFIED BY Q. DUAN, 5/16/94
        # ----------------------------------------------------------------------
        #        IF (IOHINF .EQ. 1) THEN
        DT1 = DT / 86400.
        SMCAV = SMCMAX - SMCWLT
        DMAX[0] = -ZSOIL[0] * SMCAV
        # ----------------------------------------------------------------------
        # FROZEN GROUND VERSION:
        # ----------------------------------------------------------------------
        DICE = -ZSOIL[0] * SICE[0]
        DMAX[0] = DMAX[0] * (1.0 - (SH2OA[0] + SICE[0] - SMCWLT) / SMCAV)
        DD = DMAX[0]

        for KS in range(1, NSOIL):
            # ----------------------------------------------------------------------
            # FROZEN GROUND VERSION:
            # ----------------------------------------------------------------------
            DICE = DICE + (ZSOIL[KS - 1] - ZSOIL[KS]) * SICE[KS]
            DMAX[KS] = (ZSOIL[KS - 1] - ZSOIL[KS]) * SMCAV
            DMAX[KS] = DMAX[KS] * (1.0 - (SH2OA[KS] + SICE[KS] - SMCWLT) / SMCAV)
            DD = DD + DMAX[KS]

        # ----------------------------------------------------------------------
        # VAL = (1.-EXP(-KDT*SQRT(DT1)))
        # IN BELOW, REMOVE THE SQRT IN ABOVE
        # ----------------------------------------------------------------------
        VAL = (1. - math.exp(-KDT * DT1))
        DDT = DD * VAL
        PX = PCPDRP * DT
        if PX < 0.0:
            PX = 0.0
        INFMAX = (PX * (DDT / (PX + DDT))) / DT

        # ----------------------------------------------------------------------
        # FROZEN GROUND VERSION:
        # REDUCTION OF INFILTRATION BASED ON FROZEN GROUND PARAMETERS
        # ----------------------------------------------------------------------
        FCR = 1.
        if DICE > 1.E-2:
            ACRT = CVFRZ * FRZX / DICE
            SUM = 1.
            IALP1 = CVFRZ - 1
            for J in range(0, IALP1):
                K = 1
                for JJ in range(J, IALP1):
                    K = K * JJ

                SUM = SUM + (ACRT ** (CVFRZ - J)) / float[K]

            FCR = 1. - math.exp(-ACRT) * SUM
        # endIF
        INFMAX = INFMAX * FCR

        # ----------------------------------------------------------------------
        # CORRECTION OF INFILTRATION LIMITATION:
        # IF INFMAX .LE. HYDROLIC CONDUCTIVITY ASSIGN INFMAX THE VALUE OF
        # HYDROLIC CONDUCTIVITY
        # ----------------------------------------------------------------------
        #         MXSMC = MAX ( SH2OA(1), SH2OA(2) )
        MXSMC = SH2OA[0]

        WDFCND(WDF, WCND, MXSMC, SMCMAX, BEXP, DKSAT, DWSAT,
               SICEMAX)

        INFMAX = max(INFMAX, WCND)
        INFMAX = min(INFMAX, PX)

        if PCPDRP > INFMAX:
            RUNOFF1 = PCPDRP - INFMAX
            PDDUM = INFMAX
    # ----------------------------------------------------------------------
    # TO AVOID SPURIOUS DRAINAGE BEHAVIOR, 'UPSTREAM DIFFERENCING' IN LINE
    # BELOW REPLACED WITH NEW APPROACH IN 2ND LINE:
    # 'MXSMC = MAX(SH2OA(1), SH2OA(2))'
    # ----------------------------------------------------------------------
    MXSMC = SH2OA[0]
    WDFCND(WDF, WCND, MXSMC, SMCMAX, BEXP, DKSAT, DWSAT,
           SICEMAX)

    # ----------------------------------------------------------------------
    # CALC THE MATRIX COEFFICIENTS AI, BI, AND CI FOR THE TOP LAYER
    # ----------------------------------------------------------------------
    DDZ = 1. / (-.5 * ZSOIL[1])
    AI[0] = 0.0
    BI[0] = WDF * DDZ / (-ZSOIL[0])
    CI[0] = -BI[0]

    # ----------------------------------------------------------------------
    # CALC RHSTT FOR THE TOP LAYER AFTER CALC'NG THE VERTICAL SOIL MOISTURE
    # GRADIENT BTWN THE TOP AND NEXT TO TOP LAYERS.
    # ----------------------------------------------------------------------
    DSMDZ = (SH2O[0] - SH2O[1]) / (-.5 * ZSOIL[1])
    RHSTT[0] = (WDF * DSMDZ + WCND - PDDUM + EDIR + ET[0]) / ZSOIL[0]
    SSTT = WDF * DSMDZ + WCND + EDIR + ET[0]

    # ----------------------------------------------------------------------
    # INITIALIZE DDZ2
    # ----------------------------------------------------------------------
    DDZ2 = 0.0

    # ----------------------------------------------------------------------
    # LOOP THRU THE REMAINING SOIL LAYERS, REPEATING THE ABV PROCESS
    # ----------------------------------------------------------------------
    for K in range(1, NSOIL):
        DENOM2 = (ZSOIL[K - 1] - ZSOIL[K])
        if K != NSOIL - 1:
            SLOPX = 1.

            # ----------------------------------------------------------------------
            # AGAIN, TO AVOID SPURIOUS DRAINAGE BEHAVIOR, 'UPSTREAM DIFFERENCING' IN
            # LINE BELOW REPLACED WITH NEW APPROACH IN 2ND LINE:
            # 'MXSMC2 = MAX (SH2OA(K), SH2OA(K+1))'
            # ----------------------------------------------------------------------
            MXSMC2 = SH2OA[K]
            WDFCND(WDF2, WCND2, MXSMC2, SMCMAX, BEXP, DKSAT, DWSAT,
                   SICEMAX)

            # ----------------------------------------------------------------------
            # CALC SOME PARTIAL PRODUCTS FOR LATER USE IN CALC'NG RHSTT
            # ----------------------------------------------------------------------
            DENOM = (ZSOIL[K - 1] - ZSOIL[K + 1])
            DSMDZ2 = (SH2O[K] - SH2O[K + 1]) / (DENOM * 0.5)

            # ----------------------------------------------------------------------
            # CALC THE MATRIX COEF, CI, AFTER CALC'NG ITS PARTIAL PRODUCT
            # ----------------------------------------------------------------------
            DDZ2 = 2.0 / DENOM
            CI[K] = -WDF2 * DDZ2 / DENOM2
        else:
            # ----------------------------------------------------------------------
            # SLOPE OF BOTTOM LAYER IS INTRODUCED
            # ----------------------------------------------------------------------
            SLOPX = SLOPE
            # ----------------------------------------------------------------------
            # RETRIEVE THE SOIL WATER DIFFUSIVITY AND HYDRAULIC CONDUCTIVITY FOR
            # THIS LAYER
            # ----------------------------------------------------------------------
            WDFCND(WDF2, WCND2, SH2OA[NSOIL - 1], SMCMAX, BEXP, DKSAT,
                   DWSAT, SICEMAX)
            # ----------------------------------------------------------------------
            # CALC A PARTIAL PRODUCT FOR LATER USE IN CALC'NG RHSTT
            # ----------------------------------------------------------------------
            DSMDZ2 = 0.0
            # ----------------------------------------------------------------------
            # SET MATRIX COEF CI TO ZERO
            # ----------------------------------------------------------------------
            CI[K] = 0.0
        # endIF
        # ----------------------------------------------------------------------
        # CALC RHSTT FOR THIS LAYER AFTER CALC'NG ITS NUMERATOR
        # ----------------------------------------------------------------------
        NUMER = (WDF2 * DSMDZ2) + SLOPX * WCND2 - (WDF * DSMDZ) \
                - WCND + ET[K]
        RHSTT[K] = NUMER / (-DENOM2)

        # ----------------------------------------------------------------------
        # CALC MATRIX COEFS, AI, AND BI FOR THIS LAYER
        # ----------------------------------------------------------------------
        AI[K] = -WDF * DDZ / DENOM2
        BI[K] = -(AI[K] + CI[K])

        # ----------------------------------------------------------------------
        # RESET VALUES OF WDF, WCND, DSMDZ, AND DDZ FOR LOOP TO NEXT LYR
        # RUNOFF2:  SUB-SURFACE OR BASEFLOW RUNOFF
        # ----------------------------------------------------------------------
        if K == NSOIL:
            RUNOFF2 = SLOPX * WCND2
        # endIF

        if K != NSOIL:
            WDF = WDF2
            WCND = WCND2
            DSMDZ = DSMDZ2
            DDZ = DDZ2

    return


# ----------------------------------------------------------------------
# END SUBROUTINE SRT
# ----------------------------------------------------------------------

# -- 21. PHYSICS SUBROUTINE ==>  SUBROUTINE SSTEP ----------------------

def SSTEP(SH2OOUT, SH2OIN, CMC, RHSTT, RHSCT, DT,
          NSOIL, SMCMAX, CMCMAX, RUNOFF3, ZSOIL, SMC, SICE,
          AI, BI, CI):
    import numpy as np
    RHSTTin = np.zeros(NSOIL, dtype=float)

    # ----------------------------------------------------------------------
    # SUBROUTINE SSTEP
    # ----------------------------------------------------------------------
    # CALCULATE/UPDATE SOIL MOISTURE CONTENT VALUES AND CANOPY MOISTURE
    # CONTENT VALUES.
    # ----------------------------------------------------------------------
    # ----------------------------------------------------------------------
    # CREATE 'AMOUNT' VALUES OF VARIABLES TO BE INPUT TO THE
    # TRI-DIAGONAL MATRIX ROUTINE.
    # ----------------------------------------------------------------------
    NSOLD = 20
    CIin = np.zeros(NSOLD, dtype=float)

    for K in range(0, NSOIL):
        RHSTT[K] = RHSTT[K] * DT
        AI[K] = AI[K] * DT
        BI[K] = 1. + BI[K] * DT
        CI[K] = CI[K] * DT

    # ----------------------------------------------------------------------
    # COPY VALUES FOR INPUT VARIABLES BEFORE CALL TO ROSR12
    # ----------------------------------------------------------------------
    for K in range(0, NSOIL):
        RHSTTin[K] = RHSTT[K]

    for K in range(0, NSOLD):
        CIin[K] = CI[K]

    # ----------------------------------------------------------------------
    # CALL ROSR12 TO SOLVE THE TRI-DIAGONAL MATRIX
    # ----------------------------------------------------------------------
    ROSR12(CI, AI, BI, CIin, RHSTTin, RHSTT, NSOIL)

    # ----------------------------------------------------------------------
    # SUM THE PREVIOUS SMC VALUE AND THE MATRIX SOLUTION TO GET A
    # NEW VALUE.  MIN ALLOWABLE VALUE OF SMC WILL BE 0.02.
    # RUNOFF3: RUNOFF WITHIN SOIL LAYERS
    # ----------------------------------------------------------------------
    WPLUS = 0.0
    RUNOFF3 = 0.
    DDZ = -ZSOIL[0]

    for K in range(0, NSOIL):
        if K != 0:
            DDZ = ZSOIL[K - 1] - ZSOIL[K]
        SH2OOUT[K] = SH2OIN[K] + CI[K] + WPLUS / DDZ
        SICE[K] = max(SICE[K], 0.)
        STOT = SH2OOUT[K] + SICE[K]
        if STOT > SMCMAX:
            if K == 0:
                DDZ = -ZSOIL[0]
            else:
                KK11 = K - 1
                DDZ = -ZSOIL[K] + ZSOIL[KK11]
            # endIF
            WPLUS = (STOT - SMCMAX) * DDZ
        else:
            WPLUS = 0.
        # endIF
        SMC[K] = max(min(STOT, SMCMAX), 0.02)
        SH2OOUT[K] = max((SMC[K] - SICE[K]), 0.0)

    RUNOFF3 = WPLUS

    # ----------------------------------------------------------------------
    # UPDATE CANOPY WATER CONTENT/INTERCEPTION (CMC).  CONVERT RHSCT TO
    # AN 'AMOUNT' VALUE AND ADD TO PREVIOUS CMC VALUE TO GET NEW CMC.
    # ----------------------------------------------------------------------
    CMC = CMC + DT * RHSCT
    if CMC < 1.E-20:
        CMC = 0.0
    CMC = min(CMC, CMCMAX)

    return


# ----------------------------------------------------------------------
# END SUBROUTINE SSTEP
# ----------------------------------------------------------------------

# -- 22. PHYSICS SUBROUTINE ==>  SUBROUTINE TBND  ----------------------

def TBND(TU, TB, ZSOIL, ZBOT, K, NSOIL):
    # ----------------------------------------------------------------------
    # SUBROUTINE TBND
    # ----------------------------------------------------------------------
    # CALCULATE TEMPERATURE ON THE BOUNDARY OF THE LAYER BY INTERPOLATION OF
    # THE MIDDLE LAYER TEMPERATURES
    # ----------------------------------------------------------------------
    T0 = 273.15

    # ----------------------------------------------------------------------
    # USE SURFACE TEMPERATURE ON THE TOP OF THE FIRST LAYER
    # ----------------------------------------------------------------------
    if K == 0:
        ZUP = 0.
    else:
        ZUP = ZSOIL[K - 1]
    # endIF

    # ----------------------------------------------------------------------
    # USE DEPTH OF THE CONSTANT BOTTOM TEMPERATURE WHEN INTERPOLATE
    # TEMPERATURE INTO THE LAST LAYER BOUNDARY
    # ----------------------------------------------------------------------
    if K == NSOIL - 1:
        ZB = 2. * ZBOT - ZSOIL[K]
    else:
        ZB = ZSOIL[K + 1]

    # ----------------------------------------------------------------------
    # LINEAR INTERPOLATION BETWEEN THE AVERAGE LAYER TEMPERATURES
    # ----------------------------------------------------------------------
    TBND1 = TU + (TB - TU) * (ZUP - ZSOIL[K]) / (ZUP - ZB)

    return TBND1


# ----------------------------------------------------------------------
# END SUBROUTINE TBND
# ----------------------------------------------------------------------

# -- 23. PHYSICS SUBROUTINE ==>  SUBROUTINE TDFCND ---------------------

def TDFCND(SMC, QZ, SMCMAX, SH2O):
    # ----------------------------------------------------------------------
    # SUBROUTINE TDFCND
    # ----------------------------------------------------------------------
    # CALCULATE THERMAL DIFFUSIVITY AND CONDUCTIVITY OF THE SOIL FOR A GIVEN
    # POINT AND TIME.
    # ----------------------------------------------------------------------
    # PETERS-LIDARD APPROACH (PETERS-LIDARD ET al., 1998)
    # June 2001 CHANGES: FROZEN SOIL CONDITION.
    # ----------------------------------------------------------------------
    # WE NOW GET QUARTZ AS AN INPUT ARGUMENT (SET IN ROUTINE REDPRM):
    #      DATA QUARTZ /0.82, 0.10, 0.25, 0.60, 0.52, &
    #                  0.35, 0.60, 0.40, 0.82/
    # ----------------------------------------------------------------------
    # IF THE SOIL HAS ANY MOISTURE CONTENT COMPUTE A PARTIAL SUM/PRODUCT
    # OTHERWISE USE A CONSTANT VALUE WHICH WORKS WELL WITH MOST SOILS
    # ----------------------------------------------------------------------
    #  THKW ......WATER THERMAL CONDUCTIVITY
    #  THKQTZ ....THERMAL CONDUCTIVITY FOR QUARTZ
    #  THKO ......THERMAL CONDUCTIVITY FOR OTHER SOIL COMPONENTS
    #  THKS ......THERMAL CONDUCTIVITY FOR THE SOLIDS COMBINED(QUARTZ+OTHER)
    #  THKICE ....ICE THERMAL CONDUCTIVITY
    #  SMCMAX ....POROSITY (= SMCMAX)
    #  QZ .........QUARTZ CONTENT (SOIL TYPE DEPENDENT)
    # ----------------------------------------------------------------------
    # USE AS IN PETERS-LIDARD, 1998 (MODIF. FROM JOHANSEN, 1975).
    #
    #                                  PABLO GRUNMANN, 08/17/98
    # REFS.:
    #      FAROUKI, O.T.,1986: THERMAL PROPERTIES OF SOILS. SERIES ON ROCK
    #              AND SOIL MECHANICS, VOL. 11, TRANS TECH, 136 PP.
    #      JOHANSEN, O., 1975: THERMAL CONDUCTIVITY OF SOILS. PH.D. THESIS,
    #              UNIVERSITY OF TRONDHEIM,
    #      PETERS-LIDARD, C. D., ET AL., 1998: THE EFFECT OF SOIL THERMAL
    #              CONDUCTIVITY PARAMETERIZATION ON SURFACE ENERGY FLUXES
    #              AND TEMPERATURES. JOURNAL OF THE ATMOSPHERIC SCIENCES,
    #              VOL. 55, PP. 1209-1224.
    # ----------------------------------------------------------------------
    # NEEDS PARAMETERS
    # POROSITY(SOIL TYPE):
    #      POROS = SMCMAX
    # SATURATION RATIO:
    import math
    SATRATIO = SMC / SMCMAX

    # PARAMETERS  W/(M.K)
    THKICE = 2.2
    THKW = 0.57
    THKO = 2.0
    #      IF (QZ .LE. 0.2) THKO = 3.0
    THKQTZ = 7.7
    # SOLIDS' CONDUCTIVITY
    THKS = math.pow(float(THKQTZ), float(QZ)) * math.pow(float(THKO), (1. - float(QZ)))

    # UNFROZEN FRACTION (FROM 1., i.e., 100%LIQUID, TO 0. (100% FROZEN))
    XUNFROZ = (SH2O + 1.E-9) / (SMC + 1.E-9)

    # UNFROZEN VOLUME FOR SATURATION (POROSITY*XUNFROZ)
    XU = XUNFROZ * SMCMAX
    # SATURATED THERMAL CONDUCTIVITY
    THKSAT = THKS ** (1. - SMCMAX) * THKICE ** (SMCMAX - XU) * THKW ** (XU)

    # DRY DENSITY IN KG/M3
    GAMMD = (1. - SMCMAX) * 2700.

    # DRY THERMAL CONDUCTIVITY IN W.M-1.K-1
    THKDRY = (0.135 * GAMMD + 64.7) / (2700. - 0.947 * GAMMD)

    if (SH2O + 0.0005) < SMC:
        # FROZEN
        AKE = SATRATIO
    else:
        # UNFROZEN
        # RANGE OF VALIDITY FOR THE KERSTEN NUMBER (AKE)
        if SATRATIO > 0.1:

            # KERSTEN NUMBER (USING "FINE" FORMULA, VALID FOR SOILS CONTAINING AT
            # LEAST 5% OF PARTICLES WITH DIAMETER LESS THAN 2.E-6 METERS.)
            # (FOR "COARSE" FORMULA, SEE PETERS-LIDARD ET AL., 1998).
            AKE = math.log(SATRATIO) + 1.0
        else:
            # USE K = KDRY
            AKE = 0.0
        #  THERMAL CONDUCTIVITY

    DF = AKE * (THKSAT - THKDRY) + THKDRY

    return DF


# ----------------------------------------------------------------------
# END SUBROUTINE TDFCND
# ----------------------------------------------------------------------


# -- 24. PHYSICS SUBROUTINE ==>  SUBROUTINE TMPAVG ---------------------

def TMPAVG(TAVG, TUP, TM, TDN, ZSOIL, NSOIL, K):
    # ----------------------------------------------------------------------
    # SUBROUTINE TMPAVG
    # ----------------------------------------------------------------------
    # CALCULATE SOIL LAYER AVERAGE TEMPERATURE (TAVG) IN FREEZING/THAWING
    # LAYER USING UP, DOWN, AND MIDDLE LAYER TEMPERATURES (TUP, TDN, TM),
    # WHERE TUP IS AT TOP BOUNDARY OF LAYER, TDN IS AT BOTTOM BOUNDARY OF
    # LAYER.  TM IS LAYER PROGNOSTIC STATE TEMPERATURE.
    # ----------------------------------------------------------------------
    T0 = 2.7315E2

    # ----------------------------------------------------------------------
    if K == 1:
        DZ = -ZSOIL[1]
    else:
        DZ = ZSOIL[K - 1] - ZSOIL[K]
    # endIF

    DZH = DZ * 0.5

    if TUP < T0:
        if TM < T0:
            if TDN < T0:
                # ----------------------------------------------------------------------
                # TUP, TM, TDN < T0
                # ----------------------------------------------------------------------
                TAVG = (TUP + 2.0 * TM + TDN) / 4.0
            else:
                # ----------------------------------------------------------------------
                # TUP & TM < T0,  TDN >= T0
                # ----------------------------------------------------------------------
                X0 = (T0 - TM) * DZH / (TDN - TM)
                TAVG = 0.5 * (TUP * DZH + TM * (DZH + X0) + T0 * (2. * DZH - X0)) / DZ
        # endIF
        else:
            if TDN < T0:
                # ----------------------------------------------------------------------
                # TUP < T0, TM >= T0, TDN < T0
                # ----------------------------------------------------------------------
                XUP = (T0 - TUP) * DZH / (TM - TUP)
                XDN = DZH - (T0 - TM) * DZH / (TDN - TM)
                TAVG = 0.5 * (TUP * XUP + T0 * (2. * DZ - XUP - XDN) + TDN * XDN) / DZ
            else:
                # ----------------------------------------------------------------------
                # TUP < T0, TM >= T0, TDN >= T0
                # ----------------------------------------------------------------------
                XUP = (T0 - TUP) * DZH / (TM - TUP)
                TAVG = 0.5 * (TUP * XUP + T0 * (2. * DZ - XUP)) / DZ
    # endIF
    # endIF
    else:
        if TM < T0:
            if TDN < T0:
                # ----------------------------------------------------------------------
                # TUP >= T0, TM < T0, TDN < T0
                # ----------------------------------------------------------------------
                XUP = DZH - (T0 - TUP) * DZH / (TM - TUP)
                TAVG = 0.5 * (T0 * (DZ - XUP) + TM * (DZH + XUP) + TDN * DZH) / DZ
            else:
                # ----------------------------------------------------------------------
                # TUP >= T0, TM < T0, TDN >= T0
                # ----------------------------------------------------------------------
                XUP = DZH - (T0 - TUP) * DZH / (TM - TUP)
                XDN = (T0 - TM) * DZH / (TDN - TM)
                TAVG = 0.5 * (T0 * (2. * DZ - XUP - XDN) + TM * (XUP + XDN)) / DZ
        # endIF
        else:
            if TDN < T0:
                # ----------------------------------------------------------------------
                # TUP >= T0, TM >= T0, TDN < T0
                # ----------------------------------------------------------------------
                XDN = DZH - (T0 - TM) * DZH / (TDN - TM)
                TAVG = (T0 * (DZ - XDN) + 0.5 * (T0 + TDN) * XDN) / DZ
            else:
                # ----------------------------------------------------------------------
                # TUP >= T0, TM >= T0, TDN >= T0
                # ----------------------------------------------------------------------
                TAVG = (TUP + 2.0 * TM + TDN) / 4.0
        # endIF
        # endIF
        # endIF

    return TAVG
    # ----------------------------------------------------------------------
    # END SUBROUTINE TMPAVG
    # ----------------------------------------------------------------------


# -- 25. PHYSICS SUBROUTINE ==>  SUBROUTINE TRANSP ---------------------

def TRANSP(ET1, NSOIL, ETP1, SMC, CMC, ZSOIL, SHDFAC, SMCWLT,
           CMCMAX, PC, CFACTR, SMCREF, SFCTMP, Q2, NROOT, RTDIS, GX=None):
    # ----------------------------------------------------------------------
    # SUBROUTINE TRANSP
    # ----------------------------------------------------------------------
    # CALCULATE TRANSPIRATION FOR THE VEG CLASS.
    # ----------------------------------------------------------------------
    # ----------------------------------------------------------------------
    # INITIALIZE PLANT TRANSP TO ZERO FOR ALL SOIL LAYERS.
    # ----------------------------------------------------------------------
    for K in range(0, NSOIL):
        ET1[K] = 0.

        # ----------------------------------------------------------------------
        # CALCULATE AN 'ADJUSTED' POTENTIAL TRANSPIRATION
        # IF STATEMENT BELOW TO AVOID TANGENT LINEAR PROBLEMS NEAR ZERO
        # NOTE: GX AND OTHER TERMS BELOW REDISTRIBUTE TRANSPIRATION BY LAYER,
        # ET(K), AS A FUNCTION OF SOIL MOISTURE AVAILABILITY, WHILE PRESERVING
        # TOTAL ETP1A.
        # ----------------------------------------------------------------------
        if CMC != 0.0:
            ETP1A = SHDFAC * PC * ETP1 * (1.0 - (CMC / CMCMAX) ** CFACTR)
        else:
            ETP1A = SHDFAC * PC * ETP1
    # endIF

    SGX = 0.0
    for I in range(0, NROOT):
        GX[I] = (SMC[I] - SMCWLT) / (SMCREF - SMCWLT)
        GX[I] = max(MIN(GX[I], 1.), 0.)
        SGX = SGX + GX[I]

    SGX = SGX / NROOT

    DENOM = 0.
    for I in range(0, NROOT):
        RTX = RTDIS[I] + GX[I] - SGX
        GX[I] = GX[I] * max(RTX, 0.)
        DENOM = DENOM + GX[I]

    if DENOM <= 0.0:
        DENOM = 1.

    for I in range(0, NROOT):
        ET1[I] = ETP1A * GX[I] / DENOM

    # ----------------------------------------------------------------------
    # ABOVE CODE ASSUMES A VERTIC
    # LY UNIFORM ROOT DISTRIBUTION
    # CODE BELOW TESTS A VARIABLE ROOT DISTRIBUTION
    # ----------------------------------------------------------------------
    #      ET(1) = ( ZSOIL(1) / ZSOIL(NROOT) ) * GX * ETP1A
    #      ET(1) = ( ZSOIL(1) / ZSOIL(NROOT) ) * ETP1A
    # ----------------------------------------------------------------------
    # USING ROOT DISTRIBUTION AS WEIGHTING FACTOR
    # ----------------------------------------------------------------------
    #      ET(1) = RTDIS(1) * ETP1A
    #      ET(1) = ETP1A * PART(1)
    # ----------------------------------------------------------------------
    # LOOP DOWN THRU THE SOIL LAYERS REPEATING THE OPERATION ABOVE,
    # BUT USING THE THICKNESS OF THE SOIL LAYER (RATHER THAN THE
    # ABSOLUTE DEPTH OF EACH LAYER) IN THE FINAL CALCULATION.
    # ----------------------------------------------------------------------
    #      DO K = 2,NROOT
    #        GX = ( SMC(K) - SMCWLT ) / ( SMCREF - SMCWLT )
    #        GX = MAX ( MIN ( GX, 1. ), 0. )
    # TEST CANOPY RESISTANCE
    #        GX = 1.0
    #        ET(K) = ((ZSOIL(K)-ZSOIL(K-1))/ZSOIL(NROOT))*GX*ETP1A
    #        ET(K) = ((ZSOIL(K)-ZSOIL(K-1))/ZSOIL(NROOT))*ETP1A
    # ----------------------------------------------------------------------
    # USING ROOT DISTRIBUTION AS WEIGHTING FACTOR
    # ----------------------------------------------------------------------
    #        ET(K) = RTDIS(K) * ETP1A
    #        ET(K) = ETP1A*PART(K)
    #      END DO
    return


# ----------------------------------------------------------------------
# END SUBROUTINE TRANSP
# ----------------------------------------------------------------------
# -- 26. PHYSICS SUBROUTINE ==>  SUBROUTINE WDFCND ---------------------

def WDFCND(WDF, WCND, SMC, SMCMAX, BEXP, DKSAT, DWSAT,
           SICEMAX):
    # ----------------------------------------------------------------------
    # SUBROUTINE WDFCND
    # ----------------------------------------------------------------------
    # CALCULATE SOIL WATER DIFFUSIVITY AND SOIL HYDRAULIC CONDUCTIVITY.
    # ----------------------------------------------------------------------
    #     CALC THE RATIO OF THE ACTUAL TO THE MAX PSBL SOIL H2O CONTENT
    # ----------------------------------------------------------------------
    SMC = SMC
    SMCMAX = SMCMAX
    FACTR1 = 0.2 / SMCMAX
    FACTR2 = SMC / SMCMAX

    # ----------------------------------------------------------------------
    # PREP AN EXPNTL COEF AND CALC THE SOIL WATER DIFFUSIVITY
    # ----------------------------------------------------------------------
    EXPON = BEXP + 2.0
    WDF = DWSAT * FACTR2 ** EXPON

    # ----------------------------------------------------------------------
    # FROZEN SOIL HYDRAULIC DIFFUSIVITY.  VERY SENSITIVE TO THE VERTICAL
    # GRADIENT OF UNFROZEN WATER. THE LATTER GRADIENT CAN BECOME VERY
    # EXTREME IN FREEZING/THAWING SITUATIONS, AND GIVEN THE RELATIVELY
    # FEW AND THICK SOIL LAYERS, THIS GRADIENT SUFFERES SERIOUS
    # TRUNCTION ERRORS YIELDING ERRONEOUSLY HIGH VERTICAL TRANSPORTS OF
    # UNFROZEN WATER IN BOTH DIRECTIONS FROM HUGE HYDRAULIC DIFFUSIVITY.
    # THEREFORE, WE FOUND WE HAD TO ARBITRARILY CONSTRAIN WDF
    # --
    # VERSION D_10CM: ........  FACTR1 = 0.2/SMCMAX
    # WEIGHTED APPROACH...................... PABLO GRUNMANN, 28_SEP_1999.
    # ----------------------------------------------------------------------
    if SICEMAX > 0.0:
        VKWGT = 1. / (1. + (500. * SICEMAX) ** 3.)
        WDF = VKWGT * WDF + (1. - VKWGT) * DWSAT * FACTR1 ** EXPON
    # endIF

    # ----------------------------------------------------------------------
    # RESET THE EXPNTL COEF AND CALC THE HYDRAULIC CONDUCTIVITY
    # ----------------------------------------------------------------------
    EXPON = (2.0 * BEXP) + 3.0
    WCND = DKSAT * FACTR2 ** EXPON

    return


# ----------------------------------------------------------------------
# END SUBROUTINE WDFCND
# ----------------------------------------------------------------------

# E(T)函数
def ETT123(T):
    import math
    LW = 0.0
    E = 0.0
    CPV = 1870.
    RV = 461.5
    CW = 4187.
    ESO = 611.2
    TO = 273.15
    LVH2O = 2.501000e6

    LW = LVH2O - (CW - CPV) * (T - TO)
    # print('8715', ESO, LW, TO, T, RV)
    ET = ESO * math.exp(LW * (1 / TO - 1 / T) / RV)
    return ET


def DQS(T):
    #  PURPOSE:  TO CALCULATE VALUES OF VAPOR PRESSURE (E)
    #            AND P * DQS/DT (P TIMES CHG IN SAT MXG RATIO WITH RESPECT
    #            TO THE CHG IN TEMP) IN SUBSTITUTION TO THE LOOK-UP TABLES.
    #
    #            FORMULAS AND CONSTANTS FROM ROGERS AND YAU, 1989.
    #                         ADDED BY PABLO J. GRUNMANN, 6/30/97.

    import math
    CPV = 1870.
    RV = 461.5
    CW = 4187.
    EPS = 0.622
    ESO = 611.2
    TO = 273.15
    LVH2O = 2.501000e+6

    #     ABOUT THE PARAMETERS:
    #     EPS ---------- WATER - DRY AIR MOLECULAR MASS RATIO, EPSILON
    # !   VALUES FOR SPECIFIC HEAT CAPACITY AND INDIVIDUAL GAS CONSTANTS

    LW = LVH2O - (CW - CPV) * (T - TO)
    ES = ESO * math.exp(LW * (1 / TO - 1 / T) / RV)
    DESDT = LW * ES / (RV * T * T)

    DQS1 = EPS * DESDT

    return DQS1


def DQSDT(SFCTMP, SFCPRS):
    #    PURPOSE:  TO RETRIEVE THE APPROPRIATE VALUE OF DQSDT (THE CHANGE
    #    =======   OF THE SATURATION MIXING RATIO WITH RESPECT TO THE
    #              CHANGE IN TEMPERATURE) FROM:
    #
    #              FORMULAS INTRODUCED IN FUNCTION DQS
    #                                  (MODIFIED BY PABLO GRUNMANN, 7/9/97).

    if (SFCTMP >= 173.0) and (SFCTMP <= 373.0):

        #  IF THE INPUT SFC AIR TEMP IS BTWN 173 K AND 373 K, USE
        #   FUNCTION DQS TO DETERMINE THE SLOPE OF SAT.MIX RATIO FUNCTION

        DQSDT = DQS(SFCTMP) / SFCPRS

    else:

        #  OTHERWISE, SET DQSDT EQUAL TO ZERO

        DQSDT = 0.0
    return DQSDT


# Noah_main
def noah_main():
    name11 = ''
    ierror = 0
    kk = 0
    NSOLD = 20
    K = 0
    ICE = 0
    NROOT = 0
    NSOIL = 0
    NSLAY = 0
    SOILTYP = 0
    SIBVEG = 0
    IJ = 0
    T = 0
    SLOPETYP = 0
    soilrz = 0.0

    DRIP = 0.0
    SLOPE = 0.0
    BETA = 0.0
    EC = 0.0
    EDIR = 0.0
    ET = [0.0 for i in range(NSOLD)]
    ETT = 0.0
    ESNOW = 0.0
    F = 0.0
    FXEXP = 0.0
    FLX1 = 0.0
    FLX2 = 0.0
    FLX3 = 0.0
    DEW = 0.0
    RIB = 0.0
    RUNOFF1 = 0.0
    RUNOFF2 = 0.0
    RUNOFF3 = 0.0
    SIGMA = 0.0
    Q1 = 0.0
    AET = 0.0
    ALB = 0.0
    ALBEDO = 0.0
    CH = 0.0
    CM = 0.0
    CMC = 0.0
    CCOND = 0.0
    DQDST = 0.0
    CZIL = 0.0
    REFKDT = 0.0
    REFDK = 0.0
    KDT = 0.0
    FRZK = 0.0
    FRZX = 0.0
    FRZFACT = 0.0
    DQSDT2 = 0.0
    DT = 0.0
    EMISS = 0.0
    ESAT = 0.0
    QSAT = 0.0
    E = 0.0
    ETA = 0.0
    EVP = 0.0
    ETP = 0.0
    FUP = 0.0
    SHTFLX = 0.0
    SOLDN = 0.0

    LWDN = 0.0
    PRCP = 0.0
    CPCP = 0.0
    ASNOW = 0.0
    PTU = 0.0

    Q2 = 0.0
    Q2SAT = 0.0
    RES = 0.0
    RNET = 0.0
    SFCSPD = 0.0
    UWIND = 0.0
    VWIND = 0.0
    SFCPRS = 0.0
    SFCTMP = 0.0
    SHDFAC = 0.0
    SHDMIN = 0.0

    SMC = [0.0 for i in range(NSOLD)]
    SNOMLT = 0.0
    SNOALB = 0.0
    MSTAVRZ = 0.0
    MSTAVTOT = 0.0
    SOILM = 0.0
    SOILRZ = 0.0
    SOIL1M = 0.0
    STC = [0.0 for i in range(NSOLD)]
    GFLX = 0.0

    T14 = 0.0
    T1V = 0.0
    T2V = 0.0
    TBOT = 0.0
    TH2 = 0.0
    TH2V = 0.0
    Z = 0.0
    ZO = 0.0
    SOLNET = 0.0

    PCPDAY = 0.0
    PCPSUM = 0.0
    ETADAY = 0.0
    ETSUM = 0.0
    RUNOFFSUM = 0.0
    RUNOFFDAY = 0.0
    SMCNW = 0.0

    SH20 = [0.0 for i in range(NSOLD)]
    SLDPTH = [0.0 for i in range(4)]
    SNOWH = 0.0
    SNEQV = 0.0
    SNCOVR = 0.0

    RSMIN = 0.0
    RGL = 0.0
    HS = 0.0
    SNUP = 0.0

    XLAI = 0.0
    RC = 0.0
    PC = 0.0

    RCS = 0.0
    RCT = 0.0
    RCQ = 0.0
    RCSOIL = 0.0

    SMCMAX = 0.0
    SMCREF = 0.0
    SMCREF1 = 0.0
    SMCWLT = 0.0
    SMCWLT1 = 0.0
    SMCDRY = 0.0
    SOILM = 0.0

    PSISAT = 0.0
    DKSAT = 0.0
    BEXP = 0.0
    DWSAT = 0.0
    F1 = 0.0
    QUARTZ = 0.0
    SMLOW = 0.0
    SMHIGH = 0.0
    MSTAVRZ, MSTAVTO = 0.0, 0.0

    R = 287.04
    CP = 1004.5
    T0 = 273.15
    LVH20 = 2.501000e6  # 科学计数法，表示2.5 *10的6次方。
    EPS = 0.622

    soilmtc = 0

    for t in range(0, int(lis.d.nch)):
        ierror = t
        DT = lis.t.ts
        TBOT = float(zd_varder['noah'][t].tempbot)  # fortran 中zd_varder['noah'][T].tempbot是不是这样改写？
        NROOT = int(zd_varder['noah'][t].VEGP[0])
        RSMIN = float(zd_varder['noah'][t].VEGP[1])
        RGL = float(zd_varder['noah'][t].VEGP[2])
        HS = float(zd_varder['noah'][t].VEGP[3])
        SNUP = float(zd_varder['noah'][t].VEGP[4])
        Z0 = float(zd_varder['noah'][t].VEGP[5])
        XLAI = float(zd_varder['noah'][t].VEGP[6])
        SHDFAC = float(zd_varder['noah'][t].VEGIP)
        SHDMIN = 0.0
        if zd_varder['noah'][t].vegt == 11:
            SHDFAC = 0.0

        NSOIL = 4
        SLDPTH[0] = 0.10
        SLDPTH[1] = 0.30
        SLDPTH[2] = 0.60
        SLDPTH[3] = 1.00

        SOILTYP = zd_varder['noah'][t].ZOBSOIL[0]
        SMCMAX = zd_varder['noah'][t].SOILP[0]
        PSISAT = zd_varder['noah'][t].SOILP[1]
        DKSAT = zd_varder['noah'][t].SOILP[2]
        BEXP = zd_varder['noah'][t].SOILP[3]
        QUARTZ = zd_varder['noah'][t].SOILP[4]
        SMCMAX = 0.3138
        PSISAT = 0.3433
        DKSAT = 198.4E-7
        BEXP = 11.6706

        SMLOW = 0.5

        SMHIGH = 6.0

        # if SMCMAX < 0:
        DWSAT = BEXP * DKSAT * (PSISAT / SMCMAX)
        F1 = math.log(PSISAT, 10) + BEXP * math.log(SMCMAX, 10) + 2.0  # 以10位底求对数。
        SMCREF1 = SMCMAX * (5.79E-9 / DKSAT) ** (1.0 / (2.0 * BEXP + 3.0))
        SMCREF = SMCREF1 + (SMCMAX - SMCREF1) / SMHIGH
        SMCWLT1 = SMCMAX * (200.0 / PSISAT) ** (-1.0 / BEXP)
        SMCWLT = SMCWLT1 - SMLOW * SMCWLT1
        SMCDRY = SMCWLT
        noahder.soildrymoist = SMCDRY

        REFDK = 41.5e-7
        REFKDT = 4.8614 * 1.3
        KDT = REFKDT * DKSAT / REFDK

        FRZK = 0.15

        FRZFACT = (SMCMAX / SMCREF) * (0.412 / 0.468)
        FRZX = FRZK * FRZFACT

        SLOPE = 1.0

        ALB = zd_varder['noah'][t].ALBSF

        SNOALB = zd_varder['noah'][t].MXSNALB

        SFCTMP = zd_varder['noah'][t].forcing[0]
        Q2 = zd_varder['noah'][t].forcing[1]
        SOLDN = zd_varder['noah'][t].forcing[2]
        LWDN = zd_varder['noah'][t].forcing[3]

        SFCSPD = zd_varder['noah'][t].forcing[4]
        SFCPRS = zd_varder['noah'][t].forcing[5]
        PRCP = zd_varder['noah'][t].forcing[6]

        Z = 2.0

        if SFCSPD <= 0.01:
            SFCSPD = 0.01
        if Q2 < 0.1e-5:
            Q2 = 0.1e-5
        # print('10124', t, SFCTMP, int(lis.d.nch))
        ESAT = ETT123(SFCTMP)

        Q2SAT = 0.622 * ESAT / (SFCPRS - (1.0 - 0.622) * ESAT)
        if Q2 >= Q2SAT:
            Q2 = Q2SAT * 0.99
        DQSDT2 = DQSDT(SFCTMP, SFCPRS)

        TH2 = SFCTMP + (0.0098 * Z)
        T2V = SFCTMP * (1.0 + 0.61 * Q2)

        T1V = zd_varder['noah'][t].T1 * (1.0 + 0.61 * Q2)
        TH2V = TH2 * (1.0 + 0.61 * Q2)

        PTU = 0.10

        SOILRZ = 0.0
        SOIL1M = 0.0

        # if zd_varder['noah'][t].vegt == 0:
        #     ICE = 1
        # else:
        ICE = 0.0
        if zd_varder['noah'][t].SNEQV == 0.0:
            ALBEDO = ALB
        else:
            ALBEDO = ALB + SNCOVR * (SNOALB - ALB)
            if ALBEDO > SNOALB:
                ALBEDO = SNOALB
        EC = 0.0
        # 该sflx函数在哪里？
        (zd_varder['noah'][t].CMC, zd_varder['noah'][t].T1, zd_varder['noah'][t].STC, zd_varder['noah'][t].SMC,
         zd_varder['noah'][t].SH2O,
         zd_varder['noah'][t].SNOWH, zd_varder['noah'][t].SNEQV, ALBEDO, zd_varder['noah'][t].CH,
         zd_varder['noah'][t].CM,
         EVP, ETA, SHTFLX,
         EC, EDIR, ET, ETT, ESNOW, DRIP, DEW,
         BETA, ETP, GFLX,
         FLX1, FLX2, FLX3,
         SNOMLT, SNCOVR,
         RUNOFF1, RUNOFF2, RUNOFF3,
         RC, PC, RCS, RCT, RCQ, RCSOIL, MSTAVRZ, MSTAVTOT, SOILM) = SFLX(T,
                                                                         ICE, DT, Z, NSOIL, SLDPTH,
                                                                         LWDN, SOLDN, SFCPRS, PRCP, SFCTMP, Q2, SFCSPD,
                                                                         TH2, Q2SAT, DQSDT2,
                                                                         SLOPE, SHDFAC, SHDMIN, PTU, ALB, SNOALB,
                                                                         RSMIN, RGL, HS, SNUP, Z0, XLAI, NROOT,
                                                                         PSISAT, BEXP, DKSAT, SMCMAX, QUARTZ, DWSAT,
                                                                         SMCWLT, SMCREF, SMCDRY, F1, KDT, FRZX, FRZFACT,
                                                                         TBOT,
                                                                         zd_varder['noah'][t].CMC,
                                                                         zd_varder['noah'][t].T1,
                                                                         zd_varder['noah'][t].STC,
                                                                         zd_varder['noah'][t].SMC,
                                                                         zd_varder['noah'][t].SH2O,
                                                                         zd_varder['noah'][t].SNOWH,
                                                                         zd_varder['noah'][t].SNEQV, ALBEDO,
                                                                         zd_varder['noah'][t].CH,
                                                                         zd_varder['noah'][t].CM,
                                                                         EVP, ETA, SHTFLX,
                                                                         EC, EDIR, ET, ETT, ESNOW, DRIP, DEW,
                                                                         BETA, ETP, GFLX,
                                                                         FLX1, FLX2, FLX3,
                                                                         SNOMLT, SNCOVR,
                                                                         RUNOFF1, RUNOFF2, RUNOFF3,
                                                                         RC, PC, RCS, RCT, RCQ, RCSOIL,
                                                                         MSTAVRZ, MSTAVTOT, SOILM)  # 83个参数

        T14 = zd_varder['noah'][t].T1 * zd_varder['noah'][t].T1 * zd_varder['noah'][t].T1 * zd_varder['noah'][t].T1
        FUP = 5.67E-8 * T14
        zd_varder['noah'][t].swnet = zd_varder['noah'][t].swnet + SOLDN * (1.0 - ALBEDO)

        zd_varder['noah'][t].lwnet = zd_varder['noah'][t].lwnet + (5.67E-8) * (zd_varder['noah'][t].T1 ** 4.0) - LWDN

        zd_varder['noah'][t].qle = zd_varder['noah'][t].qle + ETA

        zd_varder['noah'][t].qh = zd_varder['noah'][t].qh + SHTFLX
        zd_varder['noah'][t].qg = zd_varder['noah'][t].qg - GFLX

        if SFCTMP <= T0:
            zd_varder['noah'][t].snowf = zd_varder['noah'][t].snowf + PRCP
            zd_varder['noah'][t].rainf = zd_varder['noah'][t].rainf + 0.0
        else:
            zd_varder['noah'][t].snowf = zd_varder['noah'][t].snowf + 0.0
            zd_varder['noah'][t].rainf = zd_varder['noah'][t].rainf + PRCP
        zd_varder['noah'][t].evap = zd_varder['noah'][t].evap + EVP

        zd_varder['noah'][t].qs = zd_varder['noah'][t].qs + RUNOFF1 * 1000.0
        zd_varder['noah'][t].qsb = zd_varder['noah'][t].qsb + RUNOFF2 * 1000.0
        zd_varder['noah'][t].qsm = zd_varder['noah'][t].qsm + SNOMLT * 1000.0
        zd_varder['noah'][t].avgsurft = zd_varder['noah'][t].T1
        zd_varder['noah'][t].ALBEDO = ALBEDO
        zd_varder['noah'][t].swe = zd_varder['noah'][t].SNEQV * 1000.0

        zd_varder['noah'][t].soilmoist1 = zd_varder['noah'][t].SMC[0] * 1000.0 * SLDPTH[0]
        zd_varder['noah'][t].soilmoist2 = zd_varder['noah'][t].SMC[1] * 1000.0 * SLDPTH[1]
        zd_varder['noah'][t].soilmoist3 = zd_varder['noah'][t].SMC[2] * 1000.0 * SLDPTH[2]
        zd_varder['noah'][t].soilmoist4 = zd_varder['noah'][t].SMC[3] * 1000.0 * SLDPTH[3]

        zd_varder['noah'][t].soilwet = MSTAVTOT
        zd_varder['noah'][t].ecanop = zd_varder['noah'][t].ecanop + EC * 1000.0

        zd_varder['noah'][t].tveg = zd_varder['noah'][t].tveg + ETT * 1000.0
        zd_varder['noah'][t].esoil = zd_varder['noah'][t].esoil + EDIR * 1000.0

        zd_varder['noah'][t].canopint = zd_varder['noah'][t].CMC * 1000.0

        for k in range(0, int(NROOT)):
            soilrz = soilrz + (zd_varder['noah'][t].SMC[k] * SLDPTH[k] * 1000.0)
        zd_varder['noah'][t].rootmoist = soilrz
        zd_varder['noah'][t].COUNT = zd_varder['noah'][t].COUNT + 1

    return


# noahdrv
ierr = 0
nen1 = 0
win1 = 0
name = ''  # fortran 中common是什么类型？
ierror = 0
gridDesci = [0 for x in range(0, 50)]
LIS_config()
LIS_domain_init()  # createtiles_latlon in domains
LIS_lsm_init()  # noah_varder_ini
noah_varder_ini(lis.d.nch)
allocate_forcing_mem()
defnatncbinary(griddec)
LIS_dataassim_init()  #
LIS_setuplsm()
LIS_readrestart()
while not LIS_endofrun():
    if lis.a.daalg > 0:
        LIS_backuptime(0)
        evap_en4d_bkpe()
        getcyclencbinary()
        for nen in range(0, lis.a.nens):
            win1 = 0
            evap_enkf_setst(nen)
            LIS_backuptime(1)
            while win1 < lis.a.win:
                # LIS_ticktime()
                noah_dynsetup()
                getncbinary()
                time_interp_ncbinary()
                LIS_force2tile()
                LIS_lsm_main()
                evap_en4d_getst(win1, nen)
                win1 = win1 + 1
        read_evap()
        evap_enkf_update()
        LIS_backuptime(1)
    print('10276 lisdrv')
    if lis.a.daalg == 0:
        LIS_ticktime()
    noah_dynsetup()
    getncbinary()
    time_interp_ncbinary()
    LIS_force2tile()
    LIS_lsm_main()
    LIS_lsm_output()
    LIS_writerestart()
LIS_writerange()
LISdrv_clear()
baseforcing_clear()
freencbinary()
evap_en4d_clear()
noah_clear()
lis_log_blocked_msg('MSG: lisdrv -- Done')
