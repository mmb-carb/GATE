
from datetime import datetime
from numpy import arcsin, array, cos, isnan, pi, radians, sin, sqrt, tan
import numpy as np
from netCDF4 import Dataset
import os
from scipy.spatial import cKDTree
import sys
import time


# USER CONFIGURABLES
## RUN INFO
DATES = ['2012-07-18']
DATE_FORMAT = '%Y-%m-%d'
BASE_YEAR = 2012
REGIONS = range(1, 70)  # Can include OCS regions!
#REGIONS = [1, 16, 19, 25, 38, 63]  # California's extrema
## GRID INFO
GRID_DOT_FILE = 'input/grid/GRIDDOT2D.State_321x291'
MET_ZF_FILE = 'input/grid/METCRO3D_2012_01_extract_ZF_AVG'
NCOLS = 321
NROWS = 291
NLAYERS = 18
ABL_METERS = 1000
REGION_BOX_FILE = 'input/default/region_boxes.py'
## FLIGHT PATH INFO
TAKEOFF_ANGLES = [radians(10), radians(20), radians(30)]
LAND_ANGLES = [radians(2.5), radians(3), radians(3.5)]
RUNWAY_FILE = 'input/default/runway_info_cali.csv'
FLIGHT_FRACTS_FILE = 'input/default/flight_stage_fractions_20161004.csv'
## EMISSIONS INFO
EICS = [81080011400000, 81080211400000, 81080411400000, 81080611400000, 81080814000000,
        81080814300000, 81081014000000, 81081014500000, 81081214000000, 81081214500000]
AREA_FILES = ['input/emis/st_4k.ar.v0001.810.2012.2012.rf2095_snp20160627.SMOKEv4p0..ff10']
POINT_FILES = ['input/emis/st_4k.ps.v0001.810.2012.2012.rf2095_snp20160627.SMOKEv4p0.EIC14.ff10.csv']
GAI_CODES_FILE = 'input/default/gai_codes.py'
FACILITY_ID_FILE = 'input/default/facility_ids.py'
## TEMPORAL INFO
SMOKE_AREA_FILE = 'input/temporal/ATREF_pro2012_snp20160627_smk4.csv'  # TODO: Remove?
SMOKE_PNT_FILE = 'input/temporal/PTREF_pro2012_snp20160627_smk4.csv'
SMOKE_PROF_FILE = 'input/temporal/ARPTPRO_pro2012_snp20160627_smk3_smk4.csv'
## CMAQ OUTPUT INFO
VERSION = 'v0100'
GRID_CROS_FILE = 'input/grid/GRIDCRO2D.California_4km_291x321'
GSPRO_FILE = 'input/ncf/gspro.cmaq.saprc.31dec2015.all.csv'
GSREF_FILE = 'input/ncf/gsref_28july2016_2012s.txt'
WEIGHT_FILE = 'input/ncf/molecular.weights.txt'
OUT_DIR = 'output/'
SHOULD_ZIP = False


def main():
    # parse configurables
    config = {'DATES': DATES, 'DATE_FORMAT': DATE_FORMAT, 'BASE_YEAR': BASE_YEAR,
              'REGIONS': REGIONS, 'GRID_DOT_FILE': GRID_DOT_FILE, 'MET_ZF_FILE': MET_ZF_FILE,
              'NROWS': NROWS, 'NCOLS': NCOLS, 'NLAYERS': NLAYERS, 'ABL_METERS': ABL_METERS,
              'REGION_BOX_FILE': REGION_BOX_FILE,
              'TAKEOFF_ANGLES': TAKEOFF_ANGLES, 'LAND_ANGLES': LAND_ANGLES,
              'RUNWAY_FILE': RUNWAY_FILE, 'FLIGHT_FRACTS_FILE': FLIGHT_FRACTS_FILE,
              'EICS': EICS, 'AREA_FILES': AREA_FILES, 'POINT_FILES': POINT_FILES,
              'GAI_CODES_FILE': GAI_CODES_FILE, 'FACILITY_ID_FILE': FACILITY_ID_FILE,
              'SMOKE_AREA_FILE': SMOKE_AREA_FILE, 'SMOKE_PNT_FILE': SMOKE_PNT_FILE,
              'SMOKE_PROF_FILE': SMOKE_PROF_FILE, 'VERSION': VERSION,
              'GRID_CROS_FILE': GRID_CROS_FILE, 'GSPRO_FILE': GSPRO_FILE, 'GSREF_FILE': GSREF_FILE,
              'WEIGHT_FILE': WEIGHT_FILE, 'OUT_DIR': OUT_DIR, 'SHOULD_ZIP': SHOULD_ZIP}

    # parse command line
    a = 1
    while a < len(sys.argv):
        a += 1

    # run program
    gate = GATE(config)
    gate.run()


class GATE(object):

    def __init__(self, config):
        # build model
        self.emis_readr = EmissionsReader(config)
        self.temp_build = TemporalSurrogateBuilder(config)
        self.spat_build = SpatialSurrogateBuilder(config)
        self.emis_scale = EmissionsScaler(config)
        self.ncdf_write = DictToNcfWriter(config)

    def run(self):
        print('\nRunning GATE Model')
        emis = self.emis_readr.read()
        temp_surrs = self.temp_build.build(emis.keys())
        spat_surrs = self.spat_build.build(emis.keys())
        scaled_emis = self.emis_scale.scale(emis, spat_surrs, temp_surrs)
        self.ncdf_write.write(scaled_emis)


class EmissionsReader(object):

    def __init__(self, config):
        self.area_files = config['AREA_FILES']
        self.point_files = config['POINT_FILES']
        self.eics = config['EICS']
        self.regions = config['REGIONS']
        self.gai_codes = eval(open(config['GAI_CODES_FILE'], 'r').read())
        self.facility_ids = eval(open(config['FACILITY_ID_FILE'], 'r').read())
        self.airports = SpatialSurrogateBuilder.read_runways(config['RUNWAY_FILE'])
        self.airport_emis = {}

    def read(self):
        ''' generates emissions dict by: airport (w/ region), EIC, and pollutant
        '''
        print('\tReading Emissions Files')
        self._read_area_files()
        self._read_point_files()

        return self.airport_emis

    def _read_area_files(self):
        ''' read all the SMOKE-ready FF10 area files for aircraft emissions,
            and split those emissions by airport
        '''
        # read regional emissions from FF10s
        area_emis = {}
        for file_path in self.area_files:
            self._read_area_file(file_path, area_emis)

        # split emissions by airport
        self._split_area_to_airports(area_emis)

    def _read_area_file(self, file_path, region_emis):
        ''' Read a SMOKE-ready FF10 area file to get yearly aircraft Emissions (in tons)
            File Format:
            #DESC     FF10 Nonpoint format
            #DESC     AIR_BASIN,REGION_CODE,DISTRICT,CENSUS,ID,EIC,EMIS_TYPE,POLL,ANN_VALUE,,,...
            0NC,006008,NCU,,,81080011400000,,CO,4.80000019073,,,,,,,81080011400000,,,...
            0NC,006012,NCU,,,81080011400000,,CO,13.3999996185,,,,,,,81080011400000,,,...
        '''
        f = open(file_path, 'r')

        for line in f.xreadlines():
            if line.startswith('#'): continue
            ln = line.split(',')
            if len(ln) < 10: continue
            region = self.gai_codes[ln[0] + ln[1] + ln[2]]
            if region not in self.regions: continue
            eic = int(ln[5])
            if eic not in self.eics: continue
            pollutant = ln[7].upper()
            emis = float(ln[8]) / 365.0  # convert from annual to daily

            if region not in region_emis:
                region_emis[region] = {}
            if eic not in region_emis[region]:
                region_emis[region][eic] = {}
            if pollutant not in region_emis[region][eic]:
                region_emis[region][eic][pollutant] = 0.0
            region_emis[region][eic][pollutant] += emis

        f.close()

    def _split_area_to_airports(self, area_emis):
        ''' split area source aircraft emissions from regional to airport-specific
            using the number of flights at each airport
        '''
        for region, eic_emis in area_emis.iteritems():
            if region not in self.airports: continue
            if region not in self.airport_emis:
                self.airport_emis[region] = {}
            total_flights = float(sum([d['flights'] for d in self.airports[region].itervalues()]))

            for airport, airport_data in self.airports[region].iteritems():
                if airport not in self.airport_emis[region]:
                    self.airport_emis[region][airport] = {}
                fraction = self.airports[region][airport]['flights'] / total_flights

                for eic, poll_emis in eic_emis.iteritems():
                    for poll, emis in poll_emis.iteritems():
                        if eic not in self.airport_emis[region][airport]:
                            self.airport_emis[region][airport][eic] = {}
                        if poll not in self.airport_emis[region][airport][eic]:
                            self.airport_emis[region][airport][eic][poll] = 0.0
                        self.airport_emis[region][airport][eic][poll] += emis * fraction

    def _read_point_files(self):
        ''' read all the SMOKE-ready FF10 point files for aicraft emissions,
            and split those emissions by airport
        '''
        for file_path in self.point_files:
            self._read_point_file(file_path)

    def _read_point_file(self, file_path):
        ''' Read a SMOKE-ready FF10 point file to get yearly aircraft Emissions (in tons)
            File Format:
            #DESC     FF10 Point format
            #DESC     AIR_BASIN,REGION_CODE,DISTRICT,FACILITY_ID,POINT_ID,STACK_ID,SEGMENT_ID,AGY_FACILITY_ID,AGY_UNIT_ID,AGY_REL_POINT_ID,AGY_PROCESS_ID,EIC,POLL,ANN_TOTAL,ANN_PCT_RED,FACILITY_NAME,ERPTYPE,STKHGT,STKDIAM,STKTEMP,STKFLOW,STKVEL,NAICS,LONGITUDE,LATITUDE,LL_DATUM,HORIZ_COLL_MTHD,DESIGN_CAPACITY,DESIGN_CAPACITY_UNITS,SIC,FAC_SOURCE_TYPE,UNIT_TYPE_CODE,CONTROL_IDS,CONTROL_MEASURES,CURRENT_COST,CUMULATIVE_COST,PROJECTION_FACTOR,SUBMITTER_FAC_ID,CALC_METHOD,DATA_SET_ID,FACIL_CATEGORY,ORIS_FACILITY_CODE,ORIS_BOILER_ID,IPM_YN,CALC_YEAR,DATE_UPDATED,FUG_HEIGHT,FUG_WIDTH_YDIM,FUG_LENGTH_XDIM,FUG_ANGLE,ZIPCODE,ANNUAL_AVG_HR_YR,JAN_VALUE,FEB_VALUE,MAR_VALUE,APR_VALUE,MAY_VALUE,JUN_VALUE,JUL_VALUE,AUG_VALUE,SEP_VALUE,OCT_VALUE,NOV_VALUE,DEC_VALUE,JAN_PCTRED,FEB_PCTRED,MAR_PCTRED,APR_PCTRED,MAY_PCTRED,JUN_PCTRED,JUL_PCTRED,AUG_PCTRED,SEP_PCTRED,OCT_PCTRED,NOV_PCTRED,DEC_PCTRED,COMMENT
            0SC,006019,0SC,180002,3,0,1,,,,,81080411400000,CO,0.24699999392,,"Brackett Field",,121.4,11.1,699.5,,22.6,811420.0,-117.78167,34.091667,,,,,4581,,,,,,,,,,,,,,,,,,,,,91750.0,,,,,,,,,,,,,,,,,,,,,,,,,
            0SC,006019,0SC,180002,5,0,1,,,,,81080411400000,CO,194.68699646,,"Brackett Field",,121.4,11.1,699.5,,22.6,811420.0,-117.78167,34.091667,,,,,4581,,,,,,,,,,,,,,,,,,,,,91750.0,,,,,,,,,,,,,,,,,,,,,,,,,
        '''
        f = open(file_path, 'r')
        facs_not_found = set()

        for line in f.xreadlines():
            if line.startswith('#'): continue
            ln = line.split(',')
            if len(ln) < 14: continue
            fac_id = int(ln[3])
            if fac_id not in self.facility_ids:
                facs_not_found.add(str(fac_id))
                continue
            facility = self.facility_ids[fac_id]
            region = facility['gai']
            if region not in self.regions: continue
            eic = int(ln[11])
            if eic not in self.eics: continue
            poll = ln[12].upper()
            emis = float(ln[13]) / 365.0  # convert from annual to daily
            airport = facility['faa_lid']

            if region not in self.airport_emis:
                self.airport_emis[region] = {}
            if airport not in self.airport_emis[region]:
                self.airport_emis[region][airport] = {}
            if eic not in self.airport_emis[region][airport]:
                self.airport_emis[region][airport][eic] = {}
            if poll not in self.airport_emis[region][airport][eic]:
                self.airport_emis[region][airport][eic][poll] = 0.0
            self.airport_emis[region][airport][eic][poll] += emis

        if facs_not_found:
            print('\t\tThese facility IDs were not found. Their emissions will be dropped:\n\t\t' +
                  ' '.join(sorted(facs_not_found)))

        f.close()


class TemporalSurrogateBuilder(object):

    def __init__(self, config):
        self.base_year = int(config['BASE_YEAR'])
        self.date_format = config['DATE_FORMAT']
        self.dates = [datetime.strptime(d, self.date_format) for d in sorted(set(config['DATES']))]
        self.smoke_area_file = config['SMOKE_AREA_FILE']
        self.smoke_point_file = config['SMOKE_PNT_FILE']
        self.smoke_profile_file = config['SMOKE_PROF_FILE']
        self.eics = config['EICS']
        self.regions = config['REGIONS']
        self.gai_codes = eval(open(config['GAI_CODES_FILE'], 'r').read())
        self.temp_profs = self._default_profiles()

    def build(self, regions=False):
        ''' read temporal profile dict by: region (w/ default), month, DOW, and hour
        '''
        print('\tBuilding Temporal Profiles')
        # pull regions from emissions files, if you can
        if regions:
            self.regions = regions

        #profiles = self._read_profiles_file()
        #self._read_area_file(profiles)
        #self._read_point_file(profiles)

        return self._build_daily_profiles()

    def _build_daily_profiles(self):
        ''' generate the final daily, regional, and hourly profiles,
            after combining the monthly, weekly, and diurnal profiles.
        '''
        profs = {}

        # create a full set of scaling factors for each date
        for d in self.dates:
            d_str = datetime.strftime(d, self.date_format)
            profs[d_str] = {}
            dow = datetime(self.base_year, d.month, d.day).weekday()
            month = d.month - 1
            # create temporal factors for each region
            for region in self.regions:
                profs[d_str][region] = {}
                # and each airport, individually (-1 is default airport code)
                for airport in self.temp_profs[region]:
                    # monthly scaling factor
                    if 'monthly' in self.temp_profs[region][airport]:
                        factor_month = self.temp_profs[region][airport]['monthly'][month]
                    else:
                        factor_month = self.temp_profs[region][-1]['monthly'][month]

                    # dow scaling factor
                    if 'weekly' in self.temp_profs[region][airport]:
                        factor_dow = self.temp_profs[region][airport]['weekly'][dow]
                    else:
                        factor_dow = self.temp_profs[region][-1]['weekly'][dow]

                    # 24-hr diurnal scaling factors
                    if dow < 5:
                        if 'duirnal_weekday' in self.temp_profs[region][airport]:
                            factors_diurnal = self.temp_profs[region][airport]['duirnal_weekday']
                        else:
                            factors_diurnal = self.temp_profs[region][-1]['duirnal_weekday']
                    else:
                        if 'duirnal_weekend' in self.temp_profs[region][airport]:
                            factors_diurnal = self.temp_profs[region][airport]['duirnal_weekend']
                        else:
                            factors_diurnal = self.temp_profs[region][-1]['duirnal_weekend']

                    # combine scaling factors to a resultant 24-hr cycle
                    profs[d_str][region][airport] = [f * factor_month * factor_dow for f in factors_diurnal]

        return profs

    def _default_profiles(self):
        ''' Build the default temporal profiles for all regions.
            (BoTS = Beauro of Transportation Statistics)
            Monthly: found using BoTS 2015 California flight data
            Weekly:  found using BoTS 2015 California flight data
            Diurnal: flat between 6AM and 11PM (guessing: hard to find real data by hour)
        '''
        # TODO: Create profiles using 2002-2015, not just 2015!
        profiles = {}

        for r in self.regions:
            # lay out profile dict
            profiles[r] = {}
            profiles[r][-1] = {}

            # values for default airports
            profiles[r][-1]['monthly'] = [0.940019, 0.940760, 0.983350, 0.997757, 1.000840, 1.054673,
                                          1.064091, 1.058083, 0.997962, 1.001267, 0.980646, 0.980553]
            profiles[r][-1]['weekly'] = [1.043168, 1.013231, 1.027018, 1.033452, 1.039303, 0.854394, 0.989434]
            profiles[r][-1]['duirnal_weekday'] = [0, 0, 0, 0, 0, 1/18.,
                                                  1/18., 1/18., 1/18., 1/18., 1/18., 1/18.,
                                                  1/18., 1/18., 1/18., 1/18., 1/18., 1/18.,
                                                  1/18., 1/18., 1/18., 1/18., 1/18., 0]
            profiles[r][-1]['duirnal_weekend'] = [0, 0, 0, 0, 0, 1/18.,
                                                  1/18., 1/18., 1/18., 1/18., 1/18., 1/18.,
                                                  1/18., 1/18., 1/18., 1/18., 1/18., 1/18.,
                                                  1/18., 1/18., 1/18., 1/18., 1/18., 0]

        return profiles

    def _read_profiles_file(self):
        ''' read the SMOKE-ready ARPTPRO file to get profile codes
            for different temporal profiles:
            MONTHLY, WEEKLY, DIURNAL WEEKDAY, and DIURNAL WEEKEND
        '''
        section = 'MONTHLY'
        profiles = {}
        f = open(self.smoke_profile_file, 'r')
        for line in f.xreadlines():
            ln = line.strip()
            # check which section of the file we're in
            if ln.startswith('/'):
                if 'END' in ln.upper():
                    continue
                else:
                    section = ln.strip('/').replace(' ', '_')
                    profiles[section] = {}
                    continue

            # parse line
            prof_num = int(line[:5])
            data = [int(v) for v in line[5:-4].strip().split()]
            total = float(line[-4:])

            if section == 'MONTHLY':
                profiles[section][prof_num] = [12.0 * (d / total) for d in data]
            elif section == 'WEEKLY':
                profiles[section][prof_num] = [7.0 * (d / total) for d in data]
            else:
                profiles[section][prof_num] = [d / total for d in data]

        f.close()
        return profiles

    def _read_area_file(self, profiles):
        ''' Read a SMOKE-ready ATREF file to determine the temporal profiles for each region
            File Format:
            EIC,            monthly code, DOW code, diurnal code, pollutant, region, ...
            81080611400000, 4218,         22,       52,           0,         0SF0060490BA,...
            81081214000000,4219,21,37,0,0SF0060490BA,,0,,81081214000000
        '''
        f = open(self.smoke_area_file, 'r')
        for line in f.xreadlines():
            ln = line.split(',')
            eic = int(ln[0])
            if eic not in self.eics: continue
            region_code = ln[5]
            if region_code not in self.gai_codes: continue
            region = self.gai_codes[region_code]
            if region not in self.regions: continue
            month_code = int(ln[1])
            week_code = int(ln[2])
            diurn_code = int(ln[3])

            try:
                self.temp_profs[region]['monthly'][-1] = profiles['MONTHLY'][month_code]
                self.temp_profs[region]['weekly'][-1] = profiles['MONTHLY'][month_code]
                self.temp_profs[region]['monthly'][-1] = profiles['MONTHLY'][month_code]
                self.temp_profs[region]['monthly'][-1] = profiles['MONTHLY'][month_code]
            except KeyError as e:
                print('\t\tERROR: temporal code not found.' + e)

        f.close()

    def _read_point_file(self, profiles):
        '''  A unique point source is identified by a set of coabdis, facility ID, device ID, process ID, and stack ID
            File Format:
            eic, monthly temporal code, day of week code, diurnal code, pollutant, coabdis,      facid,  stk, dev, proid
            0,   3849,                  7,                81,           0,         0SC0060330SC, 145258, 4,   0,   1
            0,3849,7,81,0,0SC0060330SC,145258,3,0,1
        '''
        pass


class SpatialSurrogateBuilder(object):

    RAD_FACTOR = np.float32(pi / 180.0)  # need angles in radians

    def __init__(self, config):
        self.nrows = config['NROWS']
        self.ncols = config['NCOLS']
        self.nlayers = config['NLAYERS']
        self.abl_meters = config['ABL_METERS']
        self.zf_file = config['MET_ZF_FILE']
        self.corners_file = config['GRID_DOT_FILE']
        self.zf = None
        self._read_grid_heights()
        self.lat_dot = None
        self.lon_dot = None
        self._read_grid_corners_file()
        self.surrogates = dict((r, {}) for r in config['REGIONS'])
        self.region_boxes = eval(open(config['REGION_BOX_FILE'], 'r').read())
        self.kdtree = self.create_kdtrees()
        self.takoff_angles = config['TAKEOFF_ANGLES']
        self.landing_angles = config['LAND_ANGLES']
        self.airports = self.read_runways(config['RUNWAY_FILE'])
        self.flight_fracts_file = config['FLIGHT_FRACTS_FILE']
        self.flight_fracts = self.read_flight_fracts()

    def build(self, regions=None):
        ''' build spatial surrogate by: region, airport, pollutant, and grid cell
        '''
        print('\tBuilding Spatial Surrogates')
        # pull regions from emissions files, if you can
        if regions:
            self.regions = regions

        # build spatial surrogates
        for region in self.regions:
            if region not in self.airports:
                print('\t\tNo airports given for region #' + str(region) + '. Skipping.')
                continue
            for airport, airport_data in self.airports[region].iteritems():
                if airport not in self.surrogates[region]:
                    self.surrogates[region][airport] = {}

                for lat0, lon0, lat1, lon1 in airport_data['runways']:
                    # TODO: For efficiency, find the grid cell of both ends of the runway and memoize?
                    cell0 = tuple(self.find_grid_cell((0, lon0, lat0), region))
                    cell1 = tuple(self.find_grid_cell((0, lon1, lat1), region))

                    # landing
                    land = {}
                    for angle in self.landing_angles:
                            self.add_dict(land, self._gen_surrogate_1runway(region, lat0, lon0, lat1, lon1, angle))
                    self.scale_dict(land, len(self.landing_angles))

                    # take-off
                    toff = {}
                    for angle in self.takoff_angles:
                        self.add_dict(toff, self._gen_surrogate_1runway(region, lat1, lon1, lat0, lon0, angle))
                    self.scale_dict(toff, len(self.takoff_angles))

                    # taxi-ing
                    taxi = {cell0: 0.5, cell1: 0.5}

                    # fill surrogates by eic and pollutant
                    for eic, poll_fracts in self.flight_fracts.iteritems():
                        if eic not in self.surrogates[region][airport]:
                            self.surrogates[region][airport][eic] = {}
                        for poll, fractions in poll_fracts.iteritems():
                            # copy surrogates, so they can be reused
                            surr = land.copy()
                            taxi1 = taxi.copy()
                            toff1 = toff.copy()
                            # scale surrogate by flight phase fractions
                            self.scale_dict(surr, fractions['landing'])
                            self.scale_dict(taxi1, fractions['taxiing'])
                            self.scale_dict(toff1, fractions['takeoff'])
                            # sum three flight phase surrogates together
                            self.add_dict(surr, taxi1)
                            self.add_dict(surr, toff1)
                            # add to final spatial surrogate collection
                            self.surrogates[region][airport][eic][poll] = surr

        return self.surrogates

    @staticmethod
    def add_dict(orig, new):
        ''' sum the element in two flat dictionaries
        '''
        for key in new:
            if key not in orig:
                orig[key] = new[key]
            else:
                orig[key] += new[key]

    @staticmethod
    def scale_dict(d, factor):
        ''' scaled the float elements in a flat dictionary
        '''
        for key in d:
            d[key] *= factor

    def pprint_surrogate(self):
        ''' Test print method to inspect spatial surrogates
        '''
        print('\nSpatial Surrogates (Just for testing purposes)\n')
        # color scheme
        bounds = [0.0001, 0.001, 0.01, 0.025, 0.05, 0.06, 0.08, 0.1, 1.0]
        COLOR_SCHEME = ['#FF99FF', '#AA54FF', '#0000FF', '#0080FF',  # Rainbow: Red to Purple
                        '#58FAF4', '#00FF00', '#FFFF00', '#FF8000', '#FF0000']

        # print with colors
        for region, airports in self.surrogates.iteritems():
            for airport, eic_data in airports.iteritems():
                for eic, poll_data in eic_data.iteritems():
                    for poll, surr_data in poll_data.iteritems():
                        s = {}
                        for key,value in surr_data.iteritems():
                            for i in xrange(len(COLOR_SCHEME)):
                                if value <= bounds[i]:
                                    s[key] = COLOR_SCHEME[i]
                                    break

                        print(region, airport, eic, poll)
                        print([(coord, color) for coord,color in s.iteritems()])

    def _read_grid_heights(self):
        '''Read the heights of all the grid layers in the modeling domain.
        NOTE: Layer heights are presumed to be in units of Meters.
        NOTE: This function produces only time-independent grid layer heights.
        '''
        # read in grid cell heights
        data = Dataset(self.zf_file, 'r')

        # units must be Meters
        if data.variables[u'ZF'].units.strip() != 'M':
            raise ValueError('Grid file is not in units of meters: ' + file_path)

        self.zf = data.variables[u'ZF'][0]
        data.close()

        # validate dimensions
        if self.zf.shape != (self.nlayers, self.nrows, self.ncols):
            raise ValueError('Grid file has wrong number of vertical dimensions: ' + self.zf_file)

    def _read_grid_corners_file(self):
        '''Read the NetCDF-formatted, CMAQ-ready grid definition file "DOT file" to read
        the corners of each grid cell.
        The results should be of dimensions one more than the grid dimensions.
        '''
        # read in gridded lat/lon
        data = Dataset(self.corners_file, 'r')
        self.lat_dot = data.variables['LATD'][0][0]
        self.lon_dot = data.variables['LOND'][0][0]
        data.close()

        # validate dimensions
        if (self.lat_dot.shape[0] != self.nrows + 1) or (self.lon_dot.shape[0] != self.nrows + 1):
            raise ValueError('The grid file has the wrong number of columns: ' + self.corners_file)
        elif (self.lat_dot.shape[1] != self.ncols + 1) or (self.lon_dot.shape[1] != self.ncols + 1):
            raise ValueError('The grid file has the wrong number of rows: ' + self.corners_file)

    def _gen_surrogate_1runway(self, region, lat0, lon0, lat1, lon1, angle):
        ''' generate a sparse-matrix 3D spatial surrogate
            for a single runway, going one direction
        '''
        surrogate = {}

        # how long is the runway?
        runway_length = self.haversine(lon0, lat0, lon1, lat1)

        # build trajectory line (layers/vertical, cols/lon, rows/lat))
        h = runway_length * tan(angle)
        p1 = array([0.0, lon0, lat0], dtype=float)
        p2 = array([h, lon1, lat1], dtype=float)
        p_end = self._find_end_point(p1, p2, self.abl_meters)

        # subset grid
        start_bottom = self.find_grid_cell(p1, region)
        start_bottom[0] = 0
        end_top = self.find_grid_cell(p_end, region)
        end_top[0] = self.abl_meters

        # intersect the trajectory with the grid
        cells_by_meter = self.bresenham_line_3d(start_bottom, end_top)
        surrogate = self._convert_vertical_to_grid(cells_by_meter)
        self._normalize_surrogate(surrogate)

        return surrogate

    def _convert_vertical_to_grid(self, cells_by_meter):
        ''' The vertical grid cells are given in feet, but need to
            be coverted to grid cell number. This will create a lot of
            duplication, so the cell list is converted to a dict.
        '''
        cells = {}
        for cell in cells_by_meter:
            z_meters, x, y = cell
            z = self._find_vertical_grid_cell(z_meters, x, y)
            p = (z, x, y)
            if p not in cells:
                cells[p] = 0
            cells[p] += 1

        return cells

    @staticmethod
    def _normalize_surrogate(surrogate):
        ''' ensure that all the grid cells in the surrogate add to 1.0.
        '''
        total = float(sum(surrogate.values()))
        for cell in surrogate:
            surrogate[cell] /= total

    def _find_vertical_grid_cell(self, z_meters, x, y):
        '''Given the x and y grid cell, and the height (z)
        in meters, calculate what vertical grid cell the point
        lay in.
        '''
        z = 0
        layers = [self.zf[i][y][x] for i in xrange(len(self.zf))]
        for i, layer in enumerate(layers):
            if layer > z_meters:
                z = i
                break
            z = i

        return i

    def find_grid_cell(self, p, region):
        ''' Find the grid cell location of a single point in our 3D grid.
            (Point given as a tuple (height in meters, lon in degrees, lat in degrees)
        '''
        lat_min, lat_max = self.region_boxes[region]['lat']
        lon_min, lon_max = self.region_boxes[region]['lon']

        # define parameters
        lon0 = p[1] * self.RAD_FACTOR
        lat0 = p[2] * self.RAD_FACTOR

        # run KD Tree algorithm
        clat0,clon0 = cos(lat0),cos(lon0)
        slat0,slon0 = sin(lat0),sin(lon0)
        dist_sq_min, minindex_1d = self.kdtree[region].query([clat0*clon0, clat0*slon0, slat0])
        y, x = np.unravel_index(minindex_1d, (lat_max - lat_min, lon_max - lon_min))

        y = lat_min + y + 1
        x = lon_min + x + 1

        # truncate values that have gone past the grid boundaries
        if y < 0:
            y = 0
        elif y >= self.nrows:
            y = self.nrows - 1
        if x < 0:
            x = 0
        elif x >= self.ncols:
            x = self.ncols - 1

        # find vertical grid cell
        z = 0
        layers = [self.zf[i][y][x] for i in xrange(len(self.zf))]
        for i, layer in enumerate(layers):
            if layer > p[0]:
                z = i
                break
            z = i

        return array([z, x, y], dtype=int)

    def _is_point_in_2d_cell(self, p, x, y):
        ''' Test if the point "p" is inside the grid cell at x,y.
            Return a tuple of the shift you will need to make to find the correct grid cell.
            Returns (0, 0) when you are in the correct cell.
        '''
        x_shift = 0
        y_shift = 0

        # test the Y/lat coordinate
        if y > 0 and p[2] < self.lat_dot[x][y]:
            y_shift = -1
        elif y < (self.nrows[2] - 1) and p[2] > self.lat_dot[x][y + 1]:
            y_shift = 1

        # test the X/lon coordinate
        if x < (self.ncols[1] - 1) and p[1] < self.lon_dot[x + 1][y]:
            x_shift = 1
        elif x > 0 and p[1] > self.lon_dot[x][y]:
            x_shift = -1

        return (x_shift, y_shift)

    def create_kdtrees(self):
        """ Create a KD Tree for the entire state """
        lat_vals = self.lat_dot[:] * self.RAD_FACTOR
        lon_vals = self.lon_dot[:] * self.RAD_FACTOR

        kdtrees = {}
        for region in self.surrogates.iterkeys():
            # find the grid cell bounding box for the region in question
            lat_min, lat_max = self.region_boxes[region]['lat']
            lon_min, lon_max = self.region_boxes[region]['lon']

            # slice grid down to this region
            latvals = lat_vals[lat_min:lat_max, lon_min:lon_max]
            lonvals = lon_vals[lat_min:lat_max, lon_min:lon_max]

            # create tree
            clat,clon = cos(latvals),cos(lonvals)
            slat,slon = sin(latvals),sin(lonvals)
            triples = list(zip(np.ravel(clat*clon), np.ravel(clat*slon), np.ravel(slat)))
            kdtrees[region] = cKDTree(triples)

        return kdtrees

    @staticmethod
    def bresenham_line_3d(p1, p2):
        """ Bresenham's line algorithm, extended to 3D """
        points = []
        z0, x0, y0 = tuple(p1)
        z1, x1, y1 = tuple(p2)
        dx = abs(x1 - x0)
        dy = abs(y1 - y0)
        dz = abs(z1 - z0)
        z, x, y = z0, x0, y0
        sx = -1 if x0 > x1 else 1
        sy = -1 if y0 > y1 else 1
        sz = -1 if z0 > z1 else 1

        if dz > dx and dz > dy:
            err_x = dz / 2.0
            err_y = dz / 2.0
            while z != z1:
                points.append((z, x, y))
                err_x -= dx
                if err_x < 0:
                    x += sx
                    err_x += dz
                err_y -= dy
                if err_y < 0:
                    y += sy
                    err_y += dz
                z += sz
        elif dx > dy:
            err_z = dx / 2.0
            err_y = dx / 2.0
            while x != x1:
                points.append((z, x, y))
                err_y -= dy
                if err_y < 0:
                    y += sy
                    err_y += dx
                err_z -= dz
                if err_z < 0:
                    z += sz
                    err_z += dx
                x += sx
        else:
            err_x = dy / 2.0
            err_z = dy / 2.0
            while y != y1:
                points.append((z, x, y))
                err_x -= dx
                if err_x < 0:
                    x += sx
                    err_x += dy
                err_z -= dz
                if err_z < 0:
                    z += sz
                    err_z += dy
                y += sy

        points.append(p2)
        return points

    @staticmethod
    def haversine(lon0, lat0, lon1, lat1):
        """ Calculate the great circle distance between two points
            on the earth (specified in decimal degrees).
            Source: http://stackoverflow.com/questions/4913349/
            haversine-formula-in-python-bearing-and-distance-between-two-gps-points
        """
        # convert decimal degrees to radians
        lon0, lat0, lon1, lat1 = map(radians, [lon0, lat0, lon1, lat1])

        # haversine formula
        dlon = lon1 - lon0
        dlat = lat1 - lat0
        a = sin(dlat / 2.0) ** 2 + cos(lat0) * cos(lat1) * sin(dlon / 2.0) ** 2
        c = 2.0 * arcsin(sqrt(a))
        r = 6.371e6  # radius of Earth in meters (use 2.088768e7 for feet).

        return c * r

    @staticmethod
    def _find_end_point(p1, p2, z_end):
        ''' Given two points defining a 3D line, find the X and Y coordinates
            for a given Z coordinate. Using the eqn of a 3D line:
            P = p1 + t(p2 - p1)
            where:  P = (z_end, x_end, y_end)
            thus:   t = (z_end - p1[0]) / (p2[0] - p1[0])
                    x_end = p1[1] + t * (p2[1] - p1[1])
                    y_end = p1[2] + t * (p2[2] - p1[2])
        '''
        t = (z_end - p1[0]) / abs(p2[0] - p1[0])
        x_end = p1[1] + t * (p2[1] - p1[1])
        y_end = p1[2] + t * (p2[2] - p1[2])

        return SpatialSurrogateBuilder._nan_to_zero(array([z_end, x_end, y_end], dtype=float))

    @staticmethod
    def _nan_to_zero(a):
        ''' Change all the NaN/nan values in a numpy array to 0. '''
        a[isnan(a)] = 0
        return a

    def read_flight_fracts(self):
        ''' read the GATE-custom fractions file that divides emissions between
            the 3 flight stages by pollutant and EIC
            File format:
            EIC,Pollutant,Landing,Taxiing,Takeoff
            81080011400000,PM,0.213454075,0.420439845,0.36610608
            81080211400000,PM,0.213454075,0.420439845,0.36610608
        '''
        # open file for reading
        f = open(self.flight_fracts_file, 'r')

        # parse header
        header = f.readline().rstrip().lower().split(',')
        eics_col = header.index('eic') if 'eic' in header else 0
        poll_col = header.index('pollutant') if 'pollutant' else 1
        land_col = header.index('landing') if 'landing' else 2
        taxi_col = header.index('taxiing') if 'taxiing' else 3
        take_col = header.index('takeoff') if 'takeoff' else 4

        # read file line-by-line
        fracts = {}
        for line in f.xreadlines():
            # parse line
            ln = line.rstrip().split(',')
            if len(ln) < 5: continue
            eic = int(ln[eics_col])
            poll = ln[poll_col].upper()
            f_land = abs(float(ln[land_col]))
            f_taxi = abs(float(ln[taxi_col]))
            f_take = abs(float(ln[take_col]))
            # normalize fractions, just in case
            total = f_land + f_taxi + f_take
            f_land /= total
            f_taxi /= total
            f_take /= total

            # fill fraction dict
            if eic not in fracts:
                fracts[eic] = {}
            fracts[eic][poll] = {'landing': f_land, 'taxiing': f_taxi, 'takeoff': f_take}

        f.close()
        return fracts

    @staticmethod
    def read_runways(file_path):
        ''' Read custom runways file,
            to build a dictionary of all runways by region
            File Format:
            airport,region,runway,flights,land_lat,land_lon,takeoff_lat,takeoff_lon
            LAX,59,06L/24R,158967.0,33.9491124722,-118.431159861,33.9521039167,-118.401948917
            LAX,59,06R/24L,158967.0,33.9467474722,-118.435327222,33.9501944444,-118.401668667
        '''
        airports = {}
        f = open(file_path, 'r')

        # parse header for column numbers
        header = f.readline().rstrip().lower().split(',')
        airport_col = header.index('airport') if 'airport' in header else 0
        regions_col = header.index('region') if 'region' else 1
        flights_col = header.index('flights') if 'flights' else 3
        landlat_col = header.index('land_lat') if 'land_lat' else 4
        landlon_col = header.index('land_lon') if 'land_lon' else 5
        takelat_col = header.index('takeoff_lat') if 'takeoff_lat' else 6
        takelon_col = header.index('takeoff_lon') if 'takeoff_lon' else 7

        # read file, line by line
        for line in f.xreadlines():
            # parse line
            ln = line.rstrip().split(',')
            if len(ln) < 7: continue
            airport = ln[airport_col]
            region = int(ln[regions_col])
            flights = int(float(ln[flights_col]))
            land_lat = float(ln[landlat_col])
            land_lon = float(ln[landlon_col])
            take_lat = float(ln[takelat_col])
            take_lon = float(ln[takelon_col])

            # fill output dict
            if region not in airports:
                airports[region] = {}
            if airport not in airports[region]:
                airports[region][airport] = {'flights': 0, 'runways': []}
            airports[region][airport]['flights'] += flights
            airports[region][airport]['runways'].append((land_lat, land_lon, take_lat, take_lon))

        return airports


class EmissionsScaler(object):

    def __init__(self, config):
        self.emis = {}

    def scale(self, emis, spat_surrs, temp_surrs):
        ''' Create daily, gridded aircraft emissions
        Inputs:
            Emissions - multi-layer dictionary
                keys: region -> airport -> EIC -> pollutant -> tons/day
            Spatial Surrogates - multi-layer dictionary
                keys: region -> airport -> EIC -> pollutant -> grid cell -> fraction
            Temporal Surrogates - multi-layer dictionary (-1 is default airport code)
                keys: date_string -> region -> airport -> 24-hourly fractions
        Output:
            Gridded Emissions - multi-layer dictionary
                keys: date_string -> EIC -> hr -> poll -> grid cell -> tons/day
        '''
        print('\tScaling Emissions to Hourly 3D')

        self.emis = {}
        for d, temporal in temp_surrs.iteritems():
            self.emis[d] = {}
            for region, region_emis in emis.iteritems():

                for airport, airport_emis in region_emis.iteritems():
                    surrs = spat_surrs[region][airport]

                    diurnal = temporal[region][airport] if airport in temporal[region] else temporal[region][-1]
                    for eic, polls in airport_emis.iteritems():
                        if eic not in self.emis[d]:
                            self.emis[d][eic] = dict((hr, {}) for hr in range(24))

                        for hr in xrange(24):
                            fraction_hr = diurnal[hr]
                            if fraction_hr == 0.0:
                                continue

                            for poll, val in polls.iteritems():
                                if poll not in self.emis[d][eic][hr]:
                                    self.emis[d][eic][hr][poll] = {}
                                val0 = val * fraction_hr

                                for cell, fraction_cell in surrs[eic][poll].iteritems():
                                    if cell not in self.emis[d][eic][hr][poll]:
                                        self.emis[d][eic][hr][poll][cell] = 0.0
                                    self.emis[d][eic][hr][poll][cell] += val0 * fraction_cell

        return self.emis


class DictToNcfWriter(object):

    STONS_HR_2_G_SEC = 251.99583333333334
    POLLS = ['CO', 'NH3', 'NOX', 'SOX', 'PM', 'TOG']

    def __init__(self, config):
        self.directory = config['OUT_DIR']
        self.eics = config['EICS']
        self.nrows = config['NROWS']
        self.ncols = config['NCOLS']
        self.nlayers = config['NLAYERS']
        self.version = config['VERSION']
        self.grid_file = config['GRID_CROS_FILE']
        self.gspro_file = config['GSPRO_FILE']
        self.gsref_file = config['GSREF_FILE']
        self.weight_file = config['WEIGHT_FILE']
        truthy = [True, 1, 'True', 'true', 'TRUE', 'yes', 'zip']
        self.should_zip = True if config['SHOULD_ZIP'] in truthy else False
        self.gspro = {}
        self.gsref = {}
        self.groups = {}
        self.num_species = 0
        self.base_year = int(config['BASE_YEAR'])
        self.date_format = config['DATE_FORMAT']
        self.dates = [datetime.strptime(d, self.date_format) for d in sorted(set(config['DATES']))]
        # default NetCDF header for on-road emissions on California's 4km modeling domain
        self.header = {'IOAPI_VERSION': "$Id: @(#) ioapi library version 3.1 $" + " "*43,
                       'EXEC_ID': "????????????????" + " "*64,
                       'FTYPE': 1,               # file type ID
                       'STIME': 80000,           # start time    e.g. 80000 (for GMT)
                       'TSTEP': 10000,           # time step     e.g. 10000 (1 hour)
                       'NTHIK': 1,               # Domain: perimeter thickness (boundary files only)
                       'NCOLS': self.ncols,      # Domain: number of columns in modeling domain
                       'NROWS': self.nrows,      # Domain: number of rows in modeling domain
                       'NLAYS': 1,               # Domain: number of vertical layers
                       'GDTYP': 2,               # Domain: grid type ID (lat-lon, UTM, RADM, etc...)
                       'P_ALP': 30.0,            # Projection: alpha
                       'P_BET': 60.0,            # Projection: betha
                       'P_GAM': -120.5,          # Projection: gamma
                       'XCENT': -120.5,          # Projection: x centroid longitude
                       'YCENT': 37.0,            # Projection: y centroid latitude
                       'XORIG': -684000.0,       # Domain: -684000 for CA_4k, -84000 for SC_4k
                       'YORIG': -564000.0,       # Domain: -564000 for CA_4k, -552000 for SC_4k
                       'XCELL': 4000.0,          # Domain: x cell width in meters
                       'YCELL': 4000.0,          # Domain: y cell width in meters
                       'VGTYP': 2,               # Domain: grid type ID (lat-lon, UTM, RADM, etc...)
                       'VGTOP': 10000,           # Domain: Top Vertical layer at 10km
                       'VGLVLS': [1.0, 0.9958],  # Domain: Vertical layer locations
                       'GDNAM': "CMAQ Emissions  ",
                       'UPNAM': "combineEmis_wdwe",
                       'FILEDESC': "",
                       'HISTORY': ""}
        # Read speciation profiles & molecular weight files
        self._load_weight_file()
        self._load_gsref()
        self._load_gspro()

    def write(self, emis):
        ''' create CMAQ-ready NetCDF from results
        '''
        print('\tWriting Output NetCDF Files')
        for d in self.dates:
            d_str = datetime.strftime(d, self.date_format)
            print('\t\tWriting: ' + d_str)
            self._write_one_ncf(d, emis[d_str], d==self.dates[-1])

    def _write_one_ncf(self, date, emis, is_last_date):
        ''' Write a CMAQ-ready NetCDF file for a single day
        '''
        # get Julian date
        jdate = int(str(date.year) + datetime(self.base_year, date.month, date.day).strftime('%j'))

        # create empty netcdf file (including file path)
        out_path = self._build_state_file_path(date)
        rootgrp, gmt_shift = self._create_netcdf(out_path, date, jdate)

        # fill netcdf file with data
        self._fill_grid(emis, date, rootgrp, gmt_shift)

        # compress output file
        if self.should_zip:
            if is_last_date:
                os.system('gzip -1 ' + out_path)
            else:
                os.system('gzip -1 ' + out_path + ' &')

    def _fill_grid(self, scaled_emissions, date, rootgrp, gmt_shift):
        ''' Fill the entire modeling domain with a 3D grid for each pollutant.
            Fill the emissions values in each grid cell, for each polluant.
            Create a separate grid set for each date.

            Old Emis format: region -> date -> hr -> EIC -> poll_grid
            New Emis format: EIC -> hr -> poll -> grid cell -> tons/day
        '''
        species = {}
        for group in self.groups:
            for i in xrange(len(np.atleast_1d(self.groups[group]['species']))):
                species[self.groups[group]['species'][i]] = {'group': group, 'index': i}

        # some mass fractions are not EIC dependent
        nox_fraction = self.gspro['DEFNOX']['NOX']
        sox_fraction = self.gspro['SOX']['SOX']

        for hour in xrange(24):
            # adjust hr for DST
            if gmt_shift == '19':
                hr = (hour + 1) % 24
            else:
                hr = hour

            for poll in self.POLLS:
                for spec in self.groups[poll]['species']:
                    # get species information
                    ind = species[spec]['index']

                    # build default emissions grid, for the sum of all EICs
                    grid = np.zeros((self.nlayers, self.nrows, self.ncols))

                    for eic, eic_data in scaled_emissions.iteritems():
                        if poll not in eic_data[hour]: continue

                        # TOG and PM fractions are EIC-dependent
                        if int(eic) in self.gsref:
                            tog_fraction = self.gspro[self.gsref[int(eic)]['TOG']]['TOG']
                            pm_fraction = self.gspro[self.gsref[int(eic)]['PM']]['PM']
                        else:
                            tog_fraction = []
                            pm_fraction = []

                        # species fractions
                        fraction = (self.STONS_HR_2_G_SEC / self.groups[poll]['weights'][ind])
                        if poll == 'TOG' and len(tog_fraction):
                            fraction *= tog_fraction[ind]
                        elif poll == 'PM' and len(pm_fraction):
                            fraction *= pm_fraction[ind]
                        elif poll == 'NOX':
                            fraction *= nox_fraction[ind]
                        elif poll == 'SOX':
                            fraction *= sox_fraction[ind]

                        self._add_grid_cells(grid, eic_data[hour][poll], fraction)  # TODO: Efficient?

                    # write data block to file
                    rootgrp.variables[spec][hr,:,:,:] = grid
                    # last hour is the same as the first
                    if hr == 0:
                        rootgrp.variables[spec][24,:,:,:] = grid

        rootgrp.close()

    def _add_grid_cells(self, grid, grid_cells, fraction):
        ''' Given a dictionary of (layer, row, col) -> float,
            create a 3D grid to store the emissions.
        '''
        for (z, x, y), value in grid_cells.iteritems():
            grid[(z, y, x)] += value * fraction

    def _create_netcdf(self, out_path, d, jdate):
        ''' Creates a blank CMAQ-ready NetCDF file, including all the important
            boilerplate and header information. But does not fill in any emissions data.
        '''
        date = d.strftime(self.date_format)

        # define some header variables
        current_date = int(time.strftime("%Y%j"))
        current_time = int(time.strftime("%H%M%S"))

        # create and outline NetCDF file
        rootgrp = Dataset(out_path, 'w', format='NETCDF3_CLASSIC')
        TSTEP = rootgrp.createDimension('TSTEP', None)
        DATE_TIME = rootgrp.createDimension('DATE-TIME', 2)
        LAY = rootgrp.createDimension('LAY', self.nlayers)
        VAR = rootgrp.createDimension('VAR', self.num_species)  # number of variables/species
        ROW = rootgrp.createDimension('ROW', self.nrows)        # Domain: number of rows
        COL = rootgrp.createDimension('COL', self.ncols)        # Domain: number of columns

        # define TFLAG Variable
        TFLAG = rootgrp.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME',), zlib=False)
        TFLAG.units = '<YYYYDDD,HHMMSS>'
        TFLAG.long_name = 'TFLAG'
        TFLAG.var_desc = 'Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS'

        # define variables and attribute definitions
        varl = ''
        for group in self.groups:
            for species in self.groups[group]['species']:
                rootgrp.createVariable(species, 'f4', ('TSTEP', 'LAY', 'ROW', 'COL'), zlib=False)
                rootgrp.variables[species].long_name = species
                rootgrp.variables[species].units = self.groups[group]['units']
                rootgrp.variables[species].var_desc = 'emissions'
                varl += species.ljust(16)

        # global attributes
        rootgrp.IOAPI_VERSION = self.header['IOAPI_VERSION']
        rootgrp.EXEC_ID = self.header['EXEC_ID']
        rootgrp.FTYPE = self.header['FTYPE']    # file type ID
        rootgrp.CDATE = current_date            # current date  e.g. 2013137
        rootgrp.CTIME = current_time            # current time  e.g. 50126
        rootgrp.WDATE = current_date            # current date  e.g. 2013137
        rootgrp.WTIME = current_time            # current time  e.g. 50126
        rootgrp.SDATE = jdate                   # scenario date e.g. 2010091
        rootgrp.STIME = self.header['STIME']    # start time    e.g. 80000 (for GMT)
        rootgrp.TSTEP = self.header['TSTEP']    # time step     e.g. 10000 (1 hour)
        rootgrp.NTHIK = self.header['NTHIK']    # Domain: perimeter thickness (boundary files only)
        rootgrp.NCOLS = self.header['NCOLS']    # Domain: number of columns in modeling domain
        rootgrp.NROWS = self.header['NROWS']    # Domain: number of rows in modeling domain
        rootgrp.NLAYS = self.header['NLAYS']    # Domain: number of vertical layers
        rootgrp.NVARS = self.num_species        # number of variables/species
        rootgrp.GDTYP = self.header['GDTYP']    # Domain: grid type ID (lat-lon, UTM, RADM, etc...)
        rootgrp.P_ALP = self.header['P_ALP']    # Projection: alpha
        rootgrp.P_BET = self.header['P_BET']    # Projection: betha
        rootgrp.P_GAM = self.header['P_GAM']    # Projection: gamma
        rootgrp.XCENT = self.header['XCENT']    # Projection: x centroid longitude
        rootgrp.YCENT = self.header['YCENT']    # Projection: y centroid latitude
        rootgrp.XORIG = self.header['XORIG']    # Domain: -684000 for CA_4k, -84000 for SC_4k
        rootgrp.YORIG = self.header['YORIG']    # Domain: -564000 for CA_4k, -552000 for SC_4k
        rootgrp.XCELL = self.header['XCELL']    # Domain: x cell width in meters
        rootgrp.YCELL = self.header['YCELL']    # Domain: y cell width in meters
        rootgrp.VGTYP = self.header['VGTYP']    # Domain: grid type ID (lat-lon, UTM, RADM, etc...)
        rootgrp.VGTOP = self.header['VGTOP']    # Domain: Top Vertical layer at 10km
        rootgrp.VGLVLS = self.header['VGLVLS']  # Domain: Vertical layer locations
        rootgrp.GDNAM = self.header['GDNAM']
        rootgrp.UPNAM = self.header['UPNAM']
        rootgrp.FILEDESC = self.header['FILEDESC']
        rootgrp.HISTORY = self.header['HISTORY']
        rootgrp.setncattr('VAR-LIST', varl)     # use this command b/c of python not liking hyphen '-'

        # seconds since epoch
        secs = time.mktime(time.strptime("%s 12" % jdate, "%Y%j %H"))
        gmt_shift = time.strftime("%H", time.gmtime(secs))
        secs -= (int(gmt_shift) - 8) * 60 * 60

        # build TFLAG variable
        tflag = np.ones((25, self.num_species, 2), dtype=np.int32)
        for hr in xrange(25):
            gdh = time.strftime("%Y%j %H0000", time.gmtime(secs + hr * 60 * 60))
            a_date,ghr = map(int, gdh.split())
            tflag[hr,:,0] = tflag[hr,:,0] * a_date
            tflag[hr,:,1] = tflag[hr,:,1] * ghr
        rootgrp.variables['TFLAG'][:] = tflag

        return rootgrp, gmt_shift

    def _load_gsref(self):
        ''' load the gsref file
            File Format: eic,profile,group
            0,CO,CO
            0,NH3,NH3
            0,SOx,SOX
            0,DEFNOx,NOX
            0,900,PM
        '''
        self.gsref = {}

        f = open(self.gsref_file, 'r')
        for line in f.xreadlines():
            ln = line.rstrip().split(',')
            if len(ln) != 3:
                continue
            eic = int(ln[0])
            if eic not in self.eics:
                continue
            profile = ln[1].upper()
            group = ln[2].upper()
            if eic not in self.gsref:
                self.gsref[eic] = {}
            self.gsref[eic][group] = profile

        f.close()

    def _load_weight_file(self):
        """ load molecular weight file
            File Format:
            NO          30.006      NOX    moles/s
            NO2         46.006      NOX    moles/s
            HONO        47.013      NOX    moles/s
        """
        # read molecular weight text file
        fin = open(self.weight_file,'r')
        lines = fin.read()
        fin.close()

        # read in CSV or Fortran-formatted file
        lines = lines.replace(',', ' ')
        lines = lines.split('\n')

        self.groups = {}
        # loop through file lines and
        for line in lines:
            # parse line
            columns = line.rstrip().split()
            if not columns:
                continue
            species = columns[0].upper()
            weight = np.float32(columns[1])
            group = columns[2].upper()

            # file output dict
            if group not in self.groups:
                units = columns[3]
                self.groups[group] = {'species': [], 'weights': [], 'units': units}
            self.groups[group]['species'].append(species)
            self.groups[group]['weights'].append(weight)

        # convert weight list to numpy.array
        for grp in self.groups:
            self.groups[grp]['species'] = np.array(self.groups[grp]['species'], dtype=np.dtype('a8'))
            self.groups[grp]['weights'] = np.array(self.groups[grp]['weights'], dtype=np.float32)

        # calculate the number of species total
        self.num_species = 0
        for group in self.groups:
            self.num_species += len(self.groups[group]['species'])

    def _load_gspro(self):
        ''' load the gspro file
            File Format:  group, pollutant, species, mole fraction, molecular weight=1, mass fraction
            1,TOG,CH4,3.1168E-03,1,0.0500000
            1,TOG,ALK3,9.4629E-03,1,0.5500000
            1,TOG,ETOH,5.4268E-03,1,0.2500000
        '''
        self.gspro = {}

        f = open(self.gspro_file, 'r')
        for line in f.xreadlines():
            # parse line
            ln = line.rstrip().split(',')
            profile = ln[0].upper()
            group = ln[1].upper()
            if group not in self.groups:
                sys.exit('ERROR: Group ' + group + ' not found in molecular weights file.')
            pollutant = ln[2].upper()
            try:
                poll_index = list(self.groups[group]['species']).index(pollutant)
            except ValueError:
                # we don't care about that pollutant
                pass
            # start filling output dict
            if profile not in self.gspro:
                self.gspro[profile] = {}
            if group not in self.gspro[profile]:
                self.gspro[profile][group] = np.zeros(len(self.groups[group]['species']),
                                                      dtype=np.float32)
            self.gspro[profile][group][poll_index] = np.float32(ln[5])

        f.close()

    def _build_state_file_path(self, date):
        """ Build output file directory and path for a daily, multi-region NetCDF file.
            NOTE: This method uses an extremely detailed file naming convention.
                  For example:
            st_4k.mv.v0938..2012.203107d18..e14..ncf
            [statewide]_[4km grid].[mobile source].[version 938]..[base year 2012].
            [model year 2031][month 7]d[day 18]..[EIC 14 categories]..ncf
        """
        yr, month, day = date.strftime(self.date_format).split('-')

        out_dir = os.path.join(self.directory, 'ncf')
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        # define the grid size string
        grid_size = '4k'
        grid_name = os.path.basename(self.grid_file)
        if '12km' in grid_name:
            grid_size = '12k'
        elif '36km' in grid_name:
            grid_size = '36k'
        elif '1km' in grid_name:
            grid_size = '1k'
        elif '250m' in grid_name:
            grid_size = '250m'

        # TODO: "st" = state, "ac" = aircraft, and "e14" = EIC-14 All can change
        file_name = 'st_' + grid_size + '.ac.' + self.version + '..' + str(self.base_year) + '.' + \
                    yr + month + 'd' + day + '..e14..ncf'

        return os.path.join(out_dir, file_name)


if __name__ == '__main__':
    main()
