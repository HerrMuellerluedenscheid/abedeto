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

*Abedeto* can beamforming. Run

    abedeto beam

You can let *abedeto* propose suitable greens function stores based on Crust2.0 profiles by running

    abedeto stores

Have a look into the stores and modify them if needed. After that you can process them as it is explained in the 
[Fomosto Tutorial](http://emolch.github.io/pyrocko/v0.3/fomosto.html)

Having finished this, run

    abedeto process

to generate some figures which might help to judge about the depth of the event.

