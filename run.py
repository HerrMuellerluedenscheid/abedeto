import os
from request import DataProvider
from beam_stack import BeamForming
from pyrocko import model
from util import create_directory
from request import DataProvider, CakeTiming
import logging

pjoin = os.path.join

logging.basicConfig(level='INFO')
logger = logging.getLogger('run')

def init(args):
    create_directory(args.name, args.force)
    create_directory(pjoin(args.name, 'meta'), args.force)

    length = 1000.
    e = list(model.Event.load_catalog(args.events))
    if len(e)>1:
        raise Exception('more than one event in catalog. Not implemented')
    else:
        e = e[0]
    model.Event.dump_catalog([e], pjoin(args.name, 'event.pf'))
    provider = DataProvider()
    tmin = CakeTiming(phase_selection='first(p|P|PP)-40', fallback_time=0.)
    tmax = CakeTiming(phase_selection='first(p|P|PP)+40', fallback_time=1000.)
    provider.download(e, timing=(tmin, tmax), prefix=args.name, dump_config=True)
    logger.info('.'*30)
    logger.info('Go to %s' % args.name)
    logger.info('modify get_parameters.pf')

def getagain(args):
    raise Exception('Not implemented')

def beam(args):
    provider = DataProvider.load(filename='available_arrays.yaml')
    for array_id in provider.available_arrays:
        traces = io.load_traces(pjoin('array_data', array_id, 'traces.mseed'))
        traces = [tr for trs in traces for tr in trs]
        stations = model.load_stations(pjoin('array_data', 'stations.pf'))
        BeamForming(statoins, event, traces, normalize=args.normalize)
    raise Exception('Not implemented')

def propose_stores(args):
    e = list(model.Event.load_catalog(pjoin(args.name, 'event.pf')))
    if len(e)>1:
        raise Exception('more than one event in catalog. Not implemented')
    else:
        e = e[0]
    for s in stations:
        propose_store(s, e, superdir=args.store_dir,
                      source_depth_min=args.sdmin,
                      source_depth_max=args.sdmax,
                      source_depth_delta=args.sddelta,
                      sample_rate=args.samplerate,
                      force_overwrite=args.force_overwrite)
    create_directory(pjoin(args.name, 'meta'), args.force)
    raise Exception('Not implemented')

def process(args):
    raise Exception('Not implemented')



if __name__=='__main__':
    import argparse

    parser = argparse.ArgumentParser('What was the depth, again?', add_help=True)

    group_init = parser.add_argument_group('Initialization')
    group_init.add_argument('--init',
                            action='store_true',
                            default=False,
                            help='create new job')
    group_init.add_argument('--events',
                        help='Event you don\'t know the depth of')
    group_init.add_argument('--name',
                        help='name')
    group_init.add_argument('--force',
                            action='store_true',
                            default=False,
                            help='force overwrite')

    group_getagain = parser.add_argument_group('Re-retrieve the data')
    group_getagain.add_argument('--getagain',
                                help='',
                                default=False,
                                action='store_true')

    group_getagain = parser.add_argument_group('Create stores')
    group_getagain.add_argument('--storify',
                                help='',
                                default=False,
                                action='store_true')
    group_getagain.add_argument('--super-dir',
                                dest='store_dir',
                                help='super directory where to search/create stores',
                                default='stores')
    group_getagain.add_argument('--source-depth-min',
                                dest='sdmin',
                                help='minimum source depth of store [km]. Default 0',
                                default=0.)
    group_getagain.add_argument('--source-depth-max',
                                dest='sdmax',
                                help='minimum source depth of store [km]. Default 15',
                                default=0.)
    group_getagain.add_argument('--source-depth-delta',
                                dest='sddelta',
                                help='delte source depth of store [km]. Default 1',
                                default=1.)
    group_getagain.add_argument('--sampling-rate',
                                dest='sample_rate',
                                help='samppling rate store [Hz]. Default 10',
                                default=10.)
    group_getagain.add_argument('--force-store',
                                dest='force_overwrite',
                                help='overwrite existent stores',
                                action='store_false')

    group_beam = parser.add_argument_group('Beam Forming')
    group_beam.add_argument('--beam',
                                help='run beamforming',
                                action='store_true')
    group_beam.add_argument('--normalize',
                                help='normlize by standard deviation of trace',
                                action='store_true',
                            default=True)

    args = parser.parse_args()

    if args.init:
        init(args)

    if args.beam:
        beam(args)
    # init
    # -download data, store download infos in subdir

    # getagain
    # first modify get options
    # run getagain

    # storify
    # create possible stores

    # process


