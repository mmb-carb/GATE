# GATE Model Usage Guide

> This article is a stub. It is a work in progress.

What follows is a description of all the configurable variables in the GATE model.  Though the GATE model is a single, monolythic Python script, it has quite a few configurable variables.  Each variable can be set in two ways:

1. In the script itself, in the "USER CONFIGURABLES" section.
2. From the command line.

Most Linux users typically prefer to set run parameters via command line.  However, setting parameters inside the script will allow the user to save the script itself as a logged record of the run.

## Run Parameters

It is useful to note when running from the command line that all command line parsing done in GATE is done using spaces as separates between flags and options.  So you cannot put a space in a file path or string without breaking the command-line parsing.

#### DATES

> dates to model aircraft emissions

This is a list of date strings.  The default format for each date will be `2012-12-31`, but the `DATE_FORMAT` parameter will allow the user to change that.  These dates are given with the model year, not the base year.

When configured inside the script, this variable is a list of date strings.  When configured from the command line, multiple dates are comma separated.

#### DATE_FORMAT

> Python datetime format string for the above

This is a single string in the Python `datetime` format.  The default value is `%Y-%m-%d`.  For more information on Python `datetime` formats see the [official documentation](https://docs.python.org/2/library/datetime.html#strftime-and-strptime-behavior).

#### THREE_DAY_MONTH

> True if each month can be represented by 3 days

If `True`, any month you want to run will be reduce to only three days: the second Wednesday, Saturday, and Sunday. Set it to `False` to run for the exact days in your range.

This may seem like a strange option, but it exists to help make GATE more uniform with certain other inventory sectors.

#### BASE_YEAR           

> base year

Frequently, the base year and the model year of an inventory are different.  If so, this can be set to a base year, to be differentiated from the year in the `DATES` variable above.

NOTE: The `DATES` and `BASE_YEAR` variable do not allow for GATE to be run for multiple years in one pass.

#### REGIONS             

> numerical region list

"Regions" can be counties, GAIs, states, or anything else.  This is kept vague, and regions are denoted by numbers so that GATE can more easily be used in different scenarios, like counties and states in different runs.

When configured inside the script, this variable is a list of integers.  When configured from the command line, multiple regions are comma separated.

#### NUM_PROCS           

> number of parallel processes to run (1 per day)

Each day in a GATE run is independent of the days around it.  So it is possible to run GATE in parallel; one parallel run for each day at maximum.  But the number of parallel runs that your system has resources for is the decision of the user.  The default value for this parameter is one.

#### GRID_DOT_FILE       

> path to CMAQ GRIDDOT2D file

This is the path to the CMAQ-formatted `GRIDDOT2D` file. This file defines the lat/lon coordinates of the corners of the grid cells in the modeling domain.  For the most current CMAQ documentation, look [here](https://www.epa.gov/cmaq/cmaq-documentation).

#### MET_ZF_FILE         

> path to CMAQ METCRO3D file

This is the path to the CMAQ-formatted `METCRO3d` file. Among other things, this file defines the height of all the z-layer grid cells at each I/J location around the modeling domain.  For the most current CMAQ documentation, look [here](https://www.epa.gov/cmaq/cmaq-documentation).

#### NCOLS               

> number of columns in modeling domain

An integer representing the number of West-to-East columns in the modeling domain.

#### NROWS               

> number of rows in modeling domain

An integer representing the number of North-to-South rows in the modeling domain.

#### NLAYERS             

> total number of vertical layers

An integer representing the number of vertical layers in the modeling domain.

#### NUM_NONZERO_LAYERS  

> number of vertical layers with emissions

An integer representing the number of non-zero vertical layers in the modeling domain.  This is a performance bump, because frequently the aircraft emissions will not stretch the entire vertical length of the modeling domain.  If unsure, the best default value here is the number of vertical layers.

#### ABL_METERS          

> height of the ABL, in meters

This is a float value of the approximate height of the Atmospheric Boundary Layer (ABL) in meters.

#### REGION_BOX_FILE     

> path to Python file with I/J box for each region

#### TAKEOFF_ANGLES      

> take-off angles to model

#### LAND_ANGLES         

> landing angles to model

#### RUNWAY_FILE         

> path to CSV with lat/lons for all runways

#### FLIGHT_FRACTS_FILE  

> path to CSV for species fractions by flight stage

#### CATEGORIES_FILE     

> path to Python file with aircraft EIC codes

#### AREA_FILES          

> path to FF10 file with area source emissions

#### POINT_FILES         

> path to CSV file with point source emissions

#### GAI_CODES_FILE      

> path to Python file with region code information

#### FACILITY_ID_FILE    

> path to Python file with airport FAA codes

#### TEMPORAL_FILE       

> path to CSV file with airport temporal profiles

#### VERSION             

> string used to identify the run

#### GSPRO_FILE          

> path to SMOKE-style GSPRO file

#### GSREF_FILE          

> path to SMOKE-style GSREF file

#### WEIGHT_FILE         

> path to file with molecular weights

#### OUT_DIR             

> path to output directory

#### SHOULD_ZIP          

> True if you want to gzip outputs, False otherwise

#### PRINT_TOTALS        

> True if you want to print totals to stdout

