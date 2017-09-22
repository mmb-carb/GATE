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

An integer representing the number of South-to-North rows in the modeling domain.


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

> path to CSV file with I/J box for each region

This is the path to a CSV file with the I/J bounding boxes for each region in the modeling domain.  For instance, the California County region bounding box starts like this:

    REGION,LAT_MIN,LAT_MAX,LON_MIN,LON_MAX
    1,176,195,177,194
    2,107,158,208,277

The example above starts with regions `1` (Alameda county), which stretches from grid cell 176 to 195 in the South-to-North direction and grid cells 177 to 194 in the West-to-East direction.

This file is provided to GATE to greatly improve the efficiency of doing one particular piece of math: determining which grid cell a given lat/lon belongs in.  This search is done extensively when trying to relate runway positions to the modeling domain. If the user does not know (and cannot calculate) the I/J bounding boxes of each region in their domain, these boxes can be replaced with the entire modeling domain.  For instance, in the California 4km statewide modeling domain this would look like:

    REGION,LAT_MIN,LAT_MAX,LON_MIN,LON_MAX
    1,0,291,0,321
    2,0,291,0,321

GATE will run with these domain-wide boxes, but it will probably run slower.


#### TAKEOFF_ANGLES

> take-off angles to model

Commercial aircraft are typical required by air traffic control to take off within a narrow range of angles.  Our research currently shows that range is between 10 and 30 degrees.

When modeling the take off and landing patterns of aircraft from a given runway, several different angles are selected within a range to represent aircraft flight patterns.  While any arbitrarily large number of angles can be chosen to improve accuracy, often just three angles in a small range will suffice.  Also, selecting more angles may slow down the model.

What the user provide is a list of one or more angles (in radians).  Inside the script, this might look like:

    TAKEOFF_ANGLES = [radians(10), radians(20), radians(30)]

Or, if the user wants to improve accuracy at the cost of model performance:

    TAKEOFF_ANGLES = [radians(10), radians(15), radians(20), radians(25), radians(30)]

And from the command line that might look like:

    -TAKEOFF_ANGLES 0.175,0.349,0.524


#### LAND_ANGLES

> landing angles to model

Commercial aircraft are typical required by air traffic control to land within a narrow range of angles.  Our research currently shows that range is between 2.5 and 3.5 degrees.

When modeling the take off and landing patterns of aircraft from a given runway, several different angles are selected within a range to represent aircraft flight patterns.  While any arbitrarily large number of angles can be chosen to improve accuracy, often just three angles in a small range will suffice.  Also, selecting more angles may slow down the model.

What the user provide is a list of one or more angles (in radians).  Inside the script, this might look like:

    LAND_ANGLES = [radians(2.5), radians(3), radians(3.5)]

And from the command line that might look like:

    -LAND_ANGLES 0.0436,0.0524,0.0611


#### RUNWAY_FILE

> path to CSV with lat/lons for all runways

This is a path to a CSV file with a complete description of the location of every runway (and helipad) in the modeling domain.  An example CSV might look like:

    airport,region,runway,flights,land_lat,land_lon,takeoff_lat,takeoff_lon
    LAX,59,06L/24R,158967.0,33.9491124722,-118.431159861,33.9521039167,-118.401948917
    LAX,59,06R/24L,158967.0,33.9467474722,-118.435327222,33.9501944444,-118.401668667
    LAX,59,07L/25R,158967.0,33.9358305833,-118.41934175,33.9398771944,-118.379776944
    LAX,59,07R/25L,158967.0,33.9336493889,-118.419018333,33.9373630278,-118.382713917
    SFO,44,01L/19R,107991.5,37.607897,-122.382928,37.62648025,-122.370609028
    SFO,44,01R/19L,107991.5,37.6063290833,-122.381040361,37.6273412222,-122.367110361
    SFO,44,10L/28R,107991.5,37.6287386111,-122.393391139,37.6135323611,-122.357141028
    SFO,44,10R/28L,107991.5,37.6262900278,-122.393105417,37.6117091389,-122.358342
    FAUX8009H,57,H,1.0,34.464612,-120.039316,H,H
    FAUX8013H,81,H,10.0,34.469136,-120.680817,H,H

It is worth explaining the columns in the above CSV file:

* **airport** - Official FAA airport code.  In the example above two official FAA airport codes are listed: LAX in Los Angeles and SFO in San Francisco. There are also two fake (faux) FAA codes listed, which were generated ad-hoc for a couple of helicopter landing pads: FAUX8009H and FAUX8013H.
* **region** - This is a single number to match the given airport/helipad with the region.  This region code must match the region code used elsewhere in the modeling.
* **runway** - This is a string representing the runway.  Many airports, and all major airports, have several runways.  In order to model aircraft emissions, we need to model each runway at an airport individually.  These strings can be any ad hoc string, though it is preferable to use official FAA code where possible. Notice that this field is meaningless for helicopter landing pads and so there the default value `H` is used.
* **flights** - Average number of flights per year.  This can be a decimal because it is an average.  Also, it should be noticed that in practice the average number of flights per year is not well known for smaller airports and helipads.
* **land_lat** - This is the latitude of the landing side of the runway.  A runway is straight line it is defined by two lat/lon points.  In the case of helipads, there is only one point, which we here call the landing side.
* **land_lon** - This is the longitude of the landing side of the runway.  A runway is straight line it is defined by two lat/lon points.  In the case of helipads, there is only one point, which we here call the landing side.
* **takeoff_lat** - This is the latitude of the take-off side of the runway.  A runway is straight line it is defined by two lat/lon points.  In the case of helipads, this takes the default value `H`.
* **takeoff_lon** - This is the longitude of the take-off side of the runway.  A runway is straight line it is defined by two lat/lon points.  In the case of helipads, this takes the default value `H`.

A good source of data for the lat/lon locations of the each end of every runway in America can be found on the commercial [Airport IQ 5010](http://www.gcr1.com/5010web/) website.


#### FLIGHT_FRACTS_FILE

> path to CSV for species fractions by flight stage

This is the path to the flight fractions CSV file.  This CSV exists to define what percentage of emissions goes into the landing, taxiing, and take-off flight stages.  The CSV contains columns for these three fractions for each EIC and pollutant.

This file exists so that emissions from aircraft can be placed into the correct spatial locations with reference to each runway.


#### CATEGORIES_FILE

> path to Python file with aircraft EIC codes

This is a single Python dictionary saved to a Python file, to explicitly list all of the EIC/SCC categories that will be modeled for aircraft.  The GATE model reads input non-spatial emissions files, which may have any number of emissions categories. This Python file exists so that only aircraft emissions are selected from the input emissions files (listed under the `AREA_FILES` and `POINT_FILES` variables).

This Python dictionary has two keys, that correspond to two aircraft categories:

* **eics** - This is a simple list of all the aircraft EICs.
* **scc2eic** - This is a mapping from SCC to EIC.

This file needs to incorporate both EIC and SCC codes because the CARB inventory uses EIC for area sources and SCC for point sources.


#### AREA_FILES

> path to FF10 file with area source emissions

This FF10 inventory file is the input for the region-wide, pre-gridded aircraft emissions.  That is, if the yearly inventory do not include aircraft emissions for a specific airport or helipad, but just have aircraft emissions for an entire county, those emissions need to go into this file.  It should be noted, this option will accept a list of such files, not just one.

The FF10 format looks much like a CSV, but has a multi-line comment header.  The FF10 format supported by GATE is not an original file format, but is taken direclty from the SMOKE model, version 3.7.  For a complete manual on this file format, see [this section](https://www.cmascenter.org/smoke/documentation/3.7/html/ch08s02.html) of the SMOKE manual.


#### POINT_FILES

> path to FF10 file with point source emissions

This FF10 inventory file is the input for the point-source, pre-gridded aircraft emissions.  That is, if the yearly inventory include aircraft emissions for specific airports / helipads, those emissions need to go into this file.  It should be noted, this option will accept a list of such files, not just one.

The FF10 format looks much like a CSV, but has a multi-line comment header.  The FF10 format supported by GATE is not an original file format, but is taken direclty from the SMOKE model, version 3.7.  For a complete manual on this file format, see [this section](https://www.cmascenter.org/smoke/documentation/3.7/html/ch08s02.html) of the SMOKE manual.


#### REGION_STRINGS_FILE

> path to CSV file with region code information

The path to a simple CSV file mapping an integer region code to the regional information in the input FF10 emissions files.  The SMOKE FF10 format provides three CSV columns of data to completely define a region.  Depending on your personal FF10 format, you may want to re-create this file with a different second column.

GATE comes provided with a couple different versions of the `region_strings.csv` file. These examples are used to support California modeling. However, this file could easily be replaced to match the user's modeling.

For instance, if the user is creating their own FF10 files from scratch, they may want to entirely forego the first and third columns and just make the second column their region code:

    #FORMAT   FF10_NONPOINT
    #...
    #DESC     AIR_BASIN,REGION_CODE,DISTRICT,...
    ,1,,...
    ,2,,...
    ...

Then their `region_strings.csv` file would be trivial:

    REGION,STRING
    1,1
    2,2
    3,3
    ...


#### FACILITY_ID_FILE

> path to CSC file with airport FAA codes

This CSV file contains flight code information about airports and helipads.  In specific, it maps a numerical ID for a given facility in the point source inventory (FF10) files, to an official `FFA_LID` code.  The user can include such a mapping for all the airports and helipads in their model, but the FF10 file format for area sources does not include specific airports/helipads, so this file will need be used for area source emissions.

The provide example file that comes with GATE looks like this:

    ID,FAA_LID,REGION,NAME
    180001,LAX,59,Los Angeles Int Airport
    180009,SNA,60,John Wayne Airport

The `ID` codes in the first column have to match the point source IDs in the FF10 inventory files.  The `FAA_LID` codes must match the codes used in the `RUNWAY_FILE`. The third column is unused, and exists purely to make this CSV easier to read.


#### TEMPORAL_FILE

> path to CSV file with airport temporal profiles

This CSV-like file includes all of the temporal profiles need to convert daily/annual inventories into hourly inventories. Each row has four comma-separated pieces of metadata:

* **region** - Region numbers, with the default value `default`.
* **airport** - Airport codes, with the default value `default`.
* **eic** - EIC vehicle codes, with the default value `default`.
* **type** - Must be one of the four temporal profile types (shown below).

After these four metadata identifiers, comes a comma-separated list of numbers that define the temporal profile, defined by (`|`).

The third column identifies which one of the four temporal profiles is being assigned:

* **diurnal_weekday** - 24 columns, giving hourly fractions for a weekday diurnal pattern, sums to 1.0
* **diurnal_weekend** - 24 columns, giving hourly fractions for a weekend diurnal pattern, sums to 1.0
* **monthly** - 12 columns, giving month-of-year fractions, sums to 12.0
* **weekly** - 7 columns, giving day-of-week fractions, sums to 7.0

The example data provided with the model is for California, and comes from a data pull done on the U.S. Beauro of Transportation Statistics [website](https://www.transtats.bts.gov/databases.asp?Mode_ID=1&Mode_Desc=Aviation&Subject_ID2=0) on December 20, 2016.


#### VERSION

> string used to identify the run

This is a short string with a version number or identifying note.  It is put into the output file names, as an identifier.


#### GSPRO_FILE

> path to SMOKE-style GSPRO file

The GSPRO and GSREF files together form a complete mapping from EIC-source category to chemical speciation.  The GSPRO and GSREF file formats were chosen from SMOKE v3.7 so that they would be standard formats that people in the community would already be familiar with.  You can find a detailed explaination of the GSPRO file format in the official documentation [here](https://www.cmascenter.org/smoke/documentation/3.7/html/ch08s05s02.html).


#### GSREF_FILE

> path to SMOKE-style GSREF file

The GSPRO and GSREF files together form a complete mapping from EIC-source category to chemical speciation.  The GSPRO and GSREF file formats were chosen from SMOKE v3.7 so that they would be standard formats that people in the community would already be familiar with.  You can find a detailed explaination of the GSREF file format in the official documentation [here](https://www.cmascenter.org/smoke/documentation/3.7/html/ch08s05s04.html).


#### WEIGHT_FILE

> path to file with molecular weights

The molecular weights text file simply provides the molecular weights for all the model species the user wants to output. It also provides the model a little metadata that is needed to speciation criteria pollutants.  The file can either be fixed format:

    NO          30.006      NOX     moles/s
    NO2         46.006      NOX     moles/s
    HONO        47.013      NOX     moles/s

Or in a headerless CSV format:

    NO,30.006,NOX,moles/s
    NO2,46.006,NOX,moles/s
    HONO,47.013,NOX,moles/s

Either format is acceptable.  While the CSV format is more in keeping with the rest of the GATE model, the fixed format will be more familiar to SMOKE users.  Either way, both are supported.

There are four columns in this file:

* **1** - species name
* **2** - molecular weight
* **3** - criteria pollutant this species belongs to
* **4** - units, for display in final NetCDF file


#### OUT_DIR

> path to output directory

This is the path to a directory for output files.  If this directory does not exist, GATE will try and create it.


#### SHOULD_ZIP

> True if you want to gzip outputs, False otherwise

The user may want output files compressed with gzip, or not.

This is a simple boolean flag, defined in the script with a Python boolean value:

    SHOULD_ZIP = True
    SHOULD_ZIP = False

But from the command line, it must be configured using a string:

    -SHOULD_ZIP True
    -SHOULD_ZIP False


#### PRINT_TOTALS

> True if you want to print totals to stdout

The user may want to print quick model-wide, daily emissions totals to standard out, or not.

This is a simple boolean flag, defined in the script with a Python boolean value:

    PRINT_TOTALS = True
    PRINT_TOTALS = False

But from the command line, it must be configured using a string:

    -PRINT_TOTALS True
    -PRINT_TOTALS False
