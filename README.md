# GATE Model

> Griding of Aircraft Trajectory Emissions Model

GATE is a command-line tool for processing raw aircraft emissions data into spatially and temporally-allocated emissions inventories, suitable for photochemical modeling or other analysis. GATE is an open-source, Python-based tool designed by the AQPSD branch of the [California EPA][CalEPA]'s [ARB][ARB].

## How to Run GATE

The GATE model is run from the Linux command line by typing:

    python GATE.py

Of course, this will only run the model with the default settings, as appear in the header of the file.  The GATE model is very flexible, and as such it has a lot of configurable variables (at the top of the script):

    ## RUN INFO
    DATES = ['2012-07-18', '...', '2012-07-20']
    DATE_FORMAT = '%Y-%m-%d'
    BASE_YEAR = 2012
    REGIONS = range(1, 70)

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

    ## TEMPORAL INFO (This section may be removed.)
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
    SHOULD_ZIP = True

You can change these manually in the script to customize your run. Or you can change them from the command line when you execute the program.

For instance, if you just want to change the date(s) of your run, you would modify the above configurable variable by doing:

    python GATE.py -DATES 2012-07-18

Notice that the command line flags are simply the name of the configurable variable prepended with a minus `-` sign. Thus, if you wanted to just change the version number, you might do:

    python GATE.py -VERSION v1234

Or you can mix and match multiple flags. If you wanted to just run the corner counties in California for one day, with a custom version number, and you wanted to leave the output file unzipped, you might do:

    python GATE.py -VERSION v0900 -REGIONS 1,16,19,25,38,63 -DATES 2012-07-18 -SHOULD_ZIP False


## Open-Source Licence

As GATE was developed by the California State Government for use in air quality modeling, it is part of the public domain. It is openly licensed under the GNU GPLv3 license, and free for all to use.

* [GNU GPLv3 License](LICENSE)


[ARB]: http://www.arb.ca.gov/homepage.htm
[CalEPA]: http://www.calepa.ca.gov/

