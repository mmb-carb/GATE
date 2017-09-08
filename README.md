# GATE Model

> Gridded Aircraft Trajectory Emissions Model

GATE is a command-line tool for processing raw aircraft emissions data into spatially and temporally-allocated emissions inventories, suitable for photochemical modeling. GATE is an open-source, Python-based tool designed by the AQPSD branch of the [California EPA][CalEPA]'s [ARB][ARB].

## How to Run GATE

The GATE model is run from the command line by typing:

    python GATE.py

Or it can be executed as a program:

    ./GATE.py

Of course, this will only run the model with the default settings, as appear in the header of the file.  The GATE model is very flexible, and has a lot of configurable variables. You can change these manually in the script to customize your run. Or you can change them from the command line when you execute the program.

For instance, if you just want to change the date(s) of your run, you would modify the above configurable variable by doing:

    python GATE.py -DATES 2012-07-18

Notice that the command line flags are simply the name of the configurable variable prepended with a minus `-` sign. Thus, if you wanted to just change the version number associated with your run, you might do:

    python GATE.py -VERSION v1234

Or you can mix and match multiple flags. If you wanted to just run the Los Angeles regions for one day, with a custom version number, and you wanted to leave the output file unzipped, you might do:

    python GATE.py -VERSION v0900 -REGIONS 59,68 -DATES 2012-07-18 -SHOULD_ZIP False

For a more detailed explaination of the config options in GATE see the [Usage Guide](USAGE.md).

## Open-Source Licence

As GATE was developed by the California State Government for use in air quality modeling, it is part of the public domain. It is openly licensed under the GNU GPLv3 license, and free for all to use.

* [GNU GPLv3 License](LICENSE)


[ARB]: http://www.arb.ca.gov/homepage.htm
[CalEPA]: http://www.calepa.ca.gov/

