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
from collections import defaultdict
import copy
import tempfile
import logging

km = 1000.
pjoin = os.path.join
logging.basicConfig(level="INFO")
logger = logging.getLogger("propose-store")

def propose_store(station, events, superdir, source_depth_min=0., source_depth_max=15.,
                  source_depth_delta=1., sample_rate=10., force_overwrite=False, numdists=2,
                  run_ttt=False, simplify=False):
    ''' Propose a fomosto store configuration for P-pP Array beam forming.
    :param event: Event instance
    :param station: Station instance
    :param superdir: where to create the store (default, current directory)
    :param source_depth_min: minimum source depth (default 0)
    :param source_depth_max: maximum source deoth (default 15)
    :param source_depth_delta: increment
    :param sample_rate: in Hz
    :param force_overwrite: overwrite potentially existent store
    :param run_ttt: generate travel time tables right away'''
    modelling_code_id = 'qseis.2006a'

    earthmodels_1d = {}
    dists = defaultdict(list)
    for e in events:
        prof = crust2x2.get_profile(e.lat, e.lon)
        dists[prof._ident].append(ortho.distance_accurate50m(e, station))
        earthmodels_1d[prof._ident] = prof

    _config = ConfigTypeA(id='[placeholder]',
                         source_depth_min=source_depth_min*km,
                         source_depth_max=source_depth_max*km,
                         source_depth_delta=source_depth_delta*km,
                         distance_min=0,
                         distance_max=0,
                         distance_delta=0,
                         sample_rate=sample_rate,
                         ncomponents=10)

    station_crust = crust2x2.get_profile(station.lat, station.lon)
    mod = cake.LayeredModel.from_scanlines(cake.from_crust2x2_profile(station_crust))
    setattr(_config, 'earthmodel_receiver_1d', mod)

    configs = []
    for ident, prof in earthmodels_1d.items():
        configid = '%s_%s_%s' % (station.station, station_crust._ident, ident)
        config = copy.copy(_config)

        config.distance_min = math.floor(min(dists[ident]))
        config.distance_max = math.ceil(max(dists[ident]))
        if len(events)==1:
            distance_delta = 1.
        else:
            distance_delta = (config.distance_max-config.distance_min)/(numdists-1)
        config.distance_delta = distance_delta
        mod = cake.load_model(crust2_profile=prof)
        if simplify:
            mod = mod.simplify(max_rel_error=0.002)
        setattr(config, 'earthmodel_1d', mod)
        configs.append(config)
        config.id = configid
        dest_dir = pjoin(superdir, config.id) 
        create_directory(dest_dir, force_overwrite)
        logger.info('Created store under: %s' % dest_dir)

        mean_z = num.mean([config.source_depth_min, config.source_depth_max])

        mean_dist = config.distance_min+(config.distance_max-config.distance_min)/2.
        arrivals = config.earthmodel_1d.arrivals(phases=cake.PhaseDef('P'),
                                                 distances=[mean_dist*cake.m2d],
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
        half_lapse_time = 70
        qs.time_region = (Timing('begin-%s' % (half_lapse_time*1.1)), Timing('begin+%s' % (half_lapse_time*1.1)))
        qs.cut = (Timing('begin-%s' % half_lapse_time), Timing('begin+%s' % half_lapse_time))
        qs.slowness_window = (0., 0., slow+slow*0.3, slow+slow*0.5)
        qs.wavelet_duration_samples = 0.001
        qs.sw_flat_earth_transform = 1
        qs.filter_shallow_paths = 1
        qs.filter_shallow_paths_depth = round(mean_dist/50000.)
        qs.validate()
        config.validate()
        Store.create_editables(dest_dir, config=config, extra={'qseis': qs})
        if run_ttt:
            st = Store(dest_dir)
            st.make_ttt()

    config_ids = [c.id for c in configs]
    return config_ids


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

    propose_store(s, events, superdir=args.superdir, force_overwrite=args.force)
