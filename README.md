# Array Beam Depth Tool

## Installation

    sudo python setup.py install

### Prerequisites:

* [pyrocko](http://emolch.github.io/pyrocko/) 


## Processing
If you need help add a *--help* to the command call in order to get additional information.
Initialize a project using the following command

    abedeto init <catalog>

where using a [pyrocko](http://emolch.github.io/pyrocko/) compatible <catalog> of events.

Changing into one of the created project directories and running

    abedeto download

should start querying the IRIS data center for available array data.

*Abedeto* can do beamforming. Run

    abedeto beam

You can let *abedeto* propose suitable greens function stores based on Crust2.0 profiles by running

    abedeto stores --depths z_min:z_max:z_delta

The proposed stores' config files contain a source and a receiver site model. These are 
combinations of the crust2 models at the top and beneath the AK135 model. 
You can modify those models as you please.
*abedeto* will set some parameters depending on the penetration depth of the
defined phase. E.g. it will remove everything beneath the turning point of the ray
path (plus 10 %) from the earth model and set a narrow slowness taper (see
sub-folder: stores/'SOME\_STORE\_ID'/extra/qseis) depending on the P arrival. This
will decrease computational effort a lot.
After that you can process them as it is explained in the 
[Fomosto Tutorial](http://emolch.github.io/pyrocko/v0.3/fomosto.html)
Most likely, you want to run the commands

    fomosto ttt
    fomosto build

Having finished this, run

    abedeto process [options]

to generate some figures which might help to judge about the depth of the event.

