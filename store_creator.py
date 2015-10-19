import os
import sys
import numpy as num
import math
from pyrocko import crust2x2, orthodrome as ortho
from pyrocko.fomosto import qseis
from pyrocko.gf.store import Store
from pyrocko.gf.meta import Timing, ConfigTypeA, TPDef
from pyrocko import cake
from util import create_directory
import copy
import tempfile
import logging
km = 1000.
pjoin = os.path.join
logging.basicConfig(level="INFO")
logger = logging.getLogger("propose-store")

def propose_store(station, event, superdir, source_depth_min=0., source_depth_max=15.,
                  source_depth_delta=1., sample_rate=10., force_overwrite=False):
    ''' Propose a fomosto store configuration for P-pP Array beam forming.
    :param event: Event instance
    :param station: Station instance
    :param superdir: where to create the store (default, current directory)
    :param source_depth_min: minimum source depth (default 0)
    :param source_depth_max: maximum source deoth (default 15)
    :param source_depth_delta: increment
    :param sample_rate: in Hz
    :param force_overwrite: overwrite potentially existent store'''
    modelling_code_id = 'qseis.2006a'
    distance = ortho.distance_accurate50m(station, event)
    depth = event.depth
    config = ConfigTypeA(id='[placeholder]',
                         source_depth_min=source_depth_min*km,
                         source_depth_max=source_depth_max*km,
                         source_depth_delta=source_depth_delta*km,
                         distance_min=0.,
                         distance_max=0.,
                         distance_delta=0.,
                         sample_rate=sample_rate,
                         ncomponents=10)

    config.distance_min = math.floor(distance-5*km)
    config.distance_max = math.ceil(distance+5*km)
    config.distance_delta = 10*km
    configid = ''
    for item in zip([station, event],['earthmodel_receiver_1d', 'earthmodel_1d']):
        location, model_id = item
        if model_id=='earthmodel_receiver_1d':
            mod = cake.LayeredModel.from_scanlines(
                cake.from_crust2x2_profile(
                    crust2x2.get_profile(location.lat, location.lon)))
            configid += '%s_%s_' % (station.station, crust2x2.get_profile(location.lat, location.lon)._ident)
        else:
            mod = cake.load_model()
            mod = mod.replaced_crust((location.lat, location.lon))
            configid += '%s' % (crust2x2.get_profile(location.lat, location.lon)._ident)
        setattr(config, model_id, mod)

    config.id = configid
    dest_dir = pjoin(superdir, config.id) 
    create_directory(dest_dir, force_overwrite)
    logger.info('Created store under: %s' % dest_dir)

    mean_z = num.mean([config.source_depth_min, config.source_depth_max])


    arrivals = config.earthmodel_1d.arrivals(phases=cake.PhaseDef('P'),
                                             distances=[distance*cake.m2d],
                                             zstart=mean_z)
    if len(arrivals)==0:
        slow = 0.1
    else:
        slow = arrivals[0].p/(cake.r2d*cake.d2m/km)

    config.modelling_code_id=modelling_code_id
    config.tabulated_phases=[
        TPDef(
            id='begin',
            definition='P,P\\,Pv_(cmb)p'),
        TPDef(
            id='end',
            definition='2.5'),
        TPDef(
            id='P',
            definition='!P')]

    qs = qseis.QSeisConfig()
    qs.qseis_version = config.modelling_code_id.split('.')[1]
    qs.time_region = (Timing('begin-50'), Timing('begin+100'))
    qs.cut = (Timing('begin-50'), Timing('begin+100'))
    qs.slowness_window = (0., 0., slow+slow*0.5, slow+slow*0.8)
    qs.filter_shallow_paths = 1
    qs.filter_shallow_paths_depth = round(distance/100000.)
    qs.validate()
    config.validate()
    Store.create_editables(dest_dir, config=config, extra={'qseis': qs})

    return dest_dir



if __name__=="__main__":
    import argparse
    from pyrocko.model import Station, Event, load_stations

    parser = argparse.ArgumentParser('suggest a store for P phases only')
    parser.add_argument('--stations',
                        help='stations file')
    parser.add_argument('--events',
                        help='event file')
    parser.add_argument('--force',
                        action='store_true',
                        help='force_overwrite')
    parser.add_argument('--superdir',
                        default='.',
                      help='directory where to put the store')
    parser.add_argument('--number_of_distances',
                      help='number of distances between outer grid nodes in GFDB',
                      default=2)

    args = parser.parse_args()

    stations = load_stations(args.stations)
    if len(stations)==1:
        s = stations[0]

    events = list(Event.load_catalog(args.events))
    if len(events)==1:
        e = events[0]

    propose_store(s, e, superdir=args.superdir, force_overwrite=args.force)
