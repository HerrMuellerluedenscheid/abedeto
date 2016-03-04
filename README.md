## Array Beam Depth Tool

### Prerequisites:

* [pyrocko](http://emolch.github.io/pyrocko/)
* [progressbar](https://pypi.python.org/pypi/progressbar) (optional)

If you don't have appropriate Green's function databases you also need to install the
modelling codes as described in the [Fomosto Tutorial](http://emolch.github.io/pyrocko/v0.3/fomosto.html) in the
"Creating a new Green's function store" paragraph.

### Download and Installation

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

This step is done based on synthetic back-azimuth and slownesses calculated using pyrocko's [cake](http://emolch.github.io/pyrocko/v0.3/cake_doc.html) module. In order to verify correct delay times used to stack traces use ``--plot``. After processing finished, there is a file *beam_shifts.png* inside each of the ``array_data`` subfolders showing the applied delay times.

Also, *abedeto* can propose suitable Green's function stores based on [Crust2.0](http://igppweb.ucsd.edu/~gabi/crust2.html) profiles:

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

### Further Information on Applications
You can specify the array-id you are going to process using e.g. ``--array-id=GERES``.

Change the y-position of the beam (blue traces) using e.g. ``--depth=1.6. This will plot the beam at the y-axis intersept coinciding with a depth of 1.6 km.

The number of synthetic traces to display can be modified using the ``--depths`` option. E.g. ``--depths=0.2:12:0.2`` will draw synthetic traces of sources between 0.2 and 12 km depths with an increment 0.2 km.

Given that you have a GF store stored in a different location other then in the *stores* sub-directory of your project, you can specify that location: ``--store-superdirs=[comma separated list of directories]``. Futhermore, if the store is not given the default name which is "ArrayID-xx-yy" (xx and yy are Crust2.0 tile IDs) the specific store to use can be defined: ``--store=STORE_ID``.

Synthetic and recorded traces might not be well aligned. This can be corrected by appending a ``--correction=[some_seconds]`` to the last command. Notice that the dashed line indicating the synthetic onset might also be shifted with respect to the synthetic traces. The reason for this is that the onset calculation is (at the moment) only performed using the source site model.

In order to change the visible time range, define a time window with respect to the synthetic onset. E.g: ``--zoom=-2.0:12`` will plot traces within a time window of 14 seconds length, starting 2 seconds before the synthetic P phase arrival.

All traces are normalized using their minmum/maximum amplitude within the visible time-range. Additionally, a gain factor can be applied using ``--gain=[X]`` (X is a float).

Filters are adjustable: e.g. ``--filter=1.1:9`` will apply a band-pass filter with corner frequencies between 1.1 and 9.0 Hz.

Currently, no restitution is performed. However, the quantity can be changed using ``--quantity [displacement|velocity]`` which either integrates the beam or differentiates the synthetic traces, respectively. It is possible to retrieve response functions for each trace using ``--get-repsponses``. These can be found afterwards in the underlying data directory and can be used to restitute, manually.

Some rudimentary information can be added as a caption to the figure by settings the ``--auto-caption`` flag.

Running the init command on events where the *name* is not specified will fail. In this case a name needs to be specified, manually:

    abedeto init catalog.pf --name

In order to avoid amiguity between projects this works only when one event is present in the given file.
By default, existing directories will not be overwritten unless appending ``--force`` to the command.

For quicker feedback ``--show`` will visualize all results right away.

Finally, the complete line to produce a meaningfull image might look like:

    abedeto process --array-id AliceSprings --depth 0.2 --depths=0.2:12:0.2 --quantity displacement --filter 0.9:9 --gain 3. --out-filename="%(array-id)s_%(quantity)s" --correction=0.4 --zoom=-2:12 --title="%(array-id)s" --show --auto-caption 

### General Information

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


### Notes on building stores
#### Handling ``number of layers (lmax) defined too small``

Especially at large distances, the ray's turning point can be deep resulting in an earthmodel (see store's config file) of great depth. In this case QSeis will fail with the exception mentioned above. In this case you can increase the parameter ``lmax`` in you QSeis' installation's ``qsglobal.h``, recompile and try again.
