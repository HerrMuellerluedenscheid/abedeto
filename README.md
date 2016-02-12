## Array Beam Depth Tool

### Prerequisites:

* [pyrocko](http://emolch.github.io/pyrocko/)

If you don't have appropriate Green's function databases you also need to install the
modelling codes as described in the [Fomosto Tutorial](http://emolch.github.io/pyrocko/v0.3/fomosto.html) in the
"Creating a new Green's function store" paragraph.

### Installation

    git clone https://github.com/HerrMuellerluedenscheid/ArrayBeamDepthTool.git
    cd ArrayBeamDepthTool
    sudo python setup.py install

### Processing
In general: If you need help add a ``--help`` to the command call in order to get additional information.

Initialize a project:

    abedeto init <catalog>

where <catalog> is a [pyrocko](http://emolch.github.io/pyrocko/) compatible <catalog> of one or several events. Have a look at the
[iquique example](https://github.com/HerrMuellerluedenscheid/ArrayBeamDepthTool/blob/master/examples) to see an example of such a file.
This will create project folders for each event within the catalog.
Change into one of the created project directories and run

    abedeto download

to start querying IRIS, Geofon and BGR data centers for available array data.

*Abedeto* can do beamforming. Run

    abedeto beam

*abedeto* can propose suitable Green's function stores based on [Crust2.0](http://igppweb.ucsd.edu/~gabi/crust2.html) profiles:

    abedeto stores

Set the depth range to test by appending

        --depths z_min:z_max:z_delta

to the previous command. Values are to be given in kilometers. Default is 0:15:1
km.
The proposed stores' config files contain a source and a receiver site model. These are 
combinations of the crust2 models at the top and beneath the AK135 model. 
You can modify those models as you please.
*abedeto* will set some parameters depending on the penetration depth of the
defined phase. E.g. it will remove everything beneath the turning point of the P ray
path (plus 10 %) from the earth model and set a narrow slowness taper (see
sub-folder: stores/'SOME\_STORE\_ID'/extra/qseis) depending on the P arrival. This
will decrease computational effort a lot.
After that you can process them as it is explained in the 
[Fomosto Tutorial](http://emolch.github.io/pyrocko/v0.3/fomosto.html).
Most likely, you want to run the commands

    fomosto ttt			# Interpolate travel time tables
    fomosto build		# Start generating data bases

Having finished this, run

    abedeto process [options]

to generate figures which might help to judge about the depth of the event. They are created within the project directory in PNG format.
Probably, synthetic and recorded traces are not well aligned. This can be corrected by appending a ``--correction [some_seconds]`` to the last command.
Currently, no restitution is performed. However, the quantity can be changed using ``--quantity [displacement|velocity]`` which either integrates the beam or differentiates the synthetic traces, respectively. It is possible to retrieve response functions for each trace using ``--get-repsponses``. These can be found afterwards in the underlying data directory and can be used to restitute, manually.

### General Information
Running the init command on events where the *name* is not specified will fail. In this case a name needs to be specified, manually:

    abedeto init catalog.pf --name

In order to avoid amiguity between projects this works only when one event is present in the given file.
By default, existing directories will not be overwritten unless appending ``--force`` to the command.

The hierarchy within the directory looks as follows::

    ProjectDir/				# Project directory
        |--array_data
           |--"SOME_ID1"		# Some Array ID
           |--"SOME_ID2"
               |--array_center.pf	# Array center location used for beam forming
               |--beam.mseed		# Beam
               |--stations.pf		# Station meta information
               |--traces.mseed		# Raw traces
           :
           :

        |--event.pf			# Event file
        |--store-mapping		# Maps store ids to array ids
        |--request.yaml			# Information concerning data selection
        |--stores
           |--StoreID1			# Green's function stores
           |--StoreID2			# The name combines the array ID and 
           :				# the ID of the Crust2x2 tile at the
           :				# source and receiver site


