import os
import shutil
from collections import defaultdict
from pyrocko.fdsn import ws
from pyrocko import model
from pyrocko import cake
from pyrocko import orthodrome as ortho
from pyrocko.gf.meta import Timing as PyrockoTiming
from pyrocko.guts import Object, String, Float, List, Dict
import logging
from util import create_directory

pjoin = os.path.join
logging.basicConfig(level='INFO')
logger = logging.getLogger('data-request')


class CakeTiming(Object):
    '''Calculates and caches phase arrivals.
    :param fallback_time: returned, when no phase arrival was found for the
                        given depth-distance-phase-selection-combination

    E.g.:
    definition = 'first(p,P)-20'
    CakeTiming(definition)'''
    phase_selection = String.T()
    fallback_time = Float.T(optional=True)
    def __init__(self, phase_selection, fallback_time=None):
        self.arrivals = defaultdict(dict)
        self.fallback_time = fallback_time
        self.which = None
        self.phase_selection = phase_selection
        _phase_selection = phase_selection
        if '+' in _phase_selection:
            _phase_selection, self.offset = _phase_selection.split('+')
            self.offset = float(self.offset)
        elif '-' in _phase_selection:
            _phase_selection, self.offset = _phase_selection.split('-')
            self.offset = float(self.offset)
            self.offset = -self.offset

        if 'first' in _phase_selection:
            self.which = 'first'
        if 'last' in _phase_selection:
            self.which = 'last'
        if self.which:
            _phase_selection = self.strip(_phase_selection)

        self.phases = _phase_selection.split('|')

    def return_time(self, ray):
        if ray == None:
            return self.fallback_time
        else:
            return ray.t + self.offset

    def t(self, mod, z_dist, get_ray=False):
        ''':param phase_selection: phase names speparated by vertical bars
        :param z_dist: tuple with (depth, distance)
        '''
        z, dist = z_dist
        if (dist, z) in self.arrivals.keys():
            return self.return_time(self.arrivals[(dist, z)])

        phases = [cake.PhaseDef(pid) for pid in self.phases]
        arrivals = mod.arrivals(distances=[dist*cake.m2d], phases=phases, zstart=z)
        if arrivals==[]:
            logger.warn('none of defined phases at d=%s, z=%s. (return fallback time)'  % (dist, z))
            want = None
        else:
            want = self.phase_selector(arrivals)
        self.arrivals[(dist, z)] = want
        if get_ray:
            return want
        else:
            return self.return_time(want)

    def phase_selector(self, _list):
        if self.which=='first':
            return min(_list, key=lambda x: x.t)
        if self.which=='last':
            return max(_list, key=lambda x: x.t)

    def strip(self, ps):
        ps = ps.replace(self.which, '')
        ps = ps.rstrip(')')
        ps = ps.lstrip('(')
        return ps

    def pyrocko_timing(self):
        return Timing()

class Timings(Object):
    timings = List.T(CakeTiming.T())
    def __init__(self, timings):
        self.timings = timings


class DataProvider(Object):
    use = List.T(String.T())
    timings = Dict.T(String.T(), Timings.T())
    #receiver_crusts = Dict.T(String.T(), String.T())
    def __init__(self, channels='SHZ', use=None, timings=None):
        self.use = use or []
        self.timings = timings or {}
        iris_arrays = {'YKA': ('CN', 'YKA*', '', channels),
                       'ESK': [('IM', 'EKB?', '', channels),
                               ('IM', 'EKR*', '', channels)],
                       'ESK1': ('IM', 'EKA?', '', channels),
                       'ILAR': ('IM', 'IL*', '', channels),
                       'IMA': ('IM', 'IM0?', '', channels),
                       'NIA': ('IM', 'I56H?', '', channels),
                       'PFIA': [('IM', 'I57H?', '', channels),
                                ('IM', 'I57L?', '', channels)],
                       'BMA' : ('IM', 'BM0?', '', channels),
                       'BCA' : ('IM', 'BC0?', '', channels),
                       'HIA':  ('IM', 'I59H?', '', channels),
                       'NVAR': ('IM', 'NV*', '', channels),
                       'PDAR': ('IM', 'PD*', '', channels),
                       'TXAR': ('IM', 'TX*', '', channels),
                       'Pilbara': ('AU', 'PSA*', '', channels),
                       'AliceSprings': ('AU', 'AS*', '', channels),
                       'GERES': [('IM', 'GEA?', '', channels),
                                 ('IM', 'GEB?', '', channels),
                                 ('IM', 'GEC?', '', channels),
                                 ('IM', 'GED?', '', channels)],
                       # Diego Garcia Hydroacoustic array
                       'DGHAland': ('IM', 'I52H?', '', channels),
                       'DGHAS': ('IM', 'H08S?', '', channels),
                       'DGHAN': ('IM', 'H08N?', '', channels),
                       #Tristan da Cunha. channels: BDF.
                       'TDC': [('IM', 'H09N?', '', channels),
                                ('IM', 'I49H?', '', channels)],
                        'NarroginIA': ('IM', 'I04H?', '', channels),
                        'CocosIslands': ('IM', 'I06H?', '' , channels),
                        'Warramunga': ('IM', 'I07H?', '', channels),
                        'BermudaIA': ('IM', 'I51H?', '', channels),
                        'FairbanksIA': ('IM', 'I53H?', '', channels)}

        geofon_arrays = { 'ROHRBACH':('6A', 'V*', '', channels),
                              'AntaOffshore': ('GR', 'I27L?', '*', channels),
                              'AntaOnshore': ('AW', 'VNA*', '*', channels),
                              'NORES': [('NO', 'NA*', '*', channels),
                                       ('NO', 'NB*', '*', channels),
                                       ('NO', 'NC*', '*', channels)]}
        bgr_arrays = {'GERES':[('GR', 'GEA?', '*', channels),
                                     ('GR', 'GEB?', '*', channels),
                                     ('GR', 'GEC?', '*', channels),
                                     ('GR', 'GED?', '*', channels)],}

        self.providers = {'iris': iris_arrays,
                          'geofon': geofon_arrays,
                          'bgr': bgr_arrays}

    def download(self, event, directory='array_data', timing=None, length=None,
                 want='all', force=False, prefix=False, dump_config=False):
        """:param want: either 'all' or ID as string or list of IDs as strings
        """
        use = []
        ts = {}
        if all([timing, length]) == None:
            raise Exception('Define one of "timing" and "length"')
        prefix = prefix or ''
        directory = pjoin(prefix, directory)
        create_directory(directory, force)
        for site, array_data_provder in self.providers.items():
            logger.info('requesting data from site %s' % site)
            for array_id, codes in array_data_provder.items():
                if array_id not in want and want!='all':
                    continue
                sub_directory = pjoin(directory, array_id)
                logger.info("fetching %s" % array_id)
                codes = array_data_provder[array_id]
                if not isinstance(codes, list):
                    codes = [codes]
                selection = [c + tuple((event.time, event.time+1000.)) for c in codes]
                logger.debug('selection: ', selection)
                try:
                    st = ws.station(site=site, selection=selection)
                except ws.EmptyResult as e:
                    logging.error('%s on %s. skip' %(e, array_id))
                    continue

                stations = st.get_pyrocko_stations()
                min_dist = min(
                    [ortho.distance_accurate50m(s, event) for s in stations])
                max_dist = max(
                    [ortho.distance_accurate50m(s, event) for s in stations])

                mod = cake.load_model(crust2_profile=(event.lat, event.lon))
                if length:
                    tstart = 0.
                    tend = length
                elif timing:
                    tstart = timing[0].t(mod, (event.depth, min_dist))
                    tend = timing[1].t(mod, (event.depth, max_dist))
                selection = [c + tuple((event.time + tstart, event.time + tend)) for c in codes]
                try:
                    d = ws.dataselect(site=site, selection=selection)
                    create_directory(sub_directory, force)
                    fn = pjoin(sub_directory, 'traces.mseed')
                    with open(fn, 'w') as f:
                        f.write(d.read())
                        f.close()
                    model.dump_stations(stations, pjoin(sub_directory, 'stations.pf'))
                    #st.dump(filename=pjoin(sub_directory, 'fdsn_request.yaml'))
                    if dump_config and timing:
                        t = Timings(list(timing))
                        ts[array_id] = t
                        t.validate()
                        #t.dump(filename=pjoin(sub_directory, 'request.conf'))
                    if array_id not in use:
                        use.append(array_id)
                except ws.EmptyResult as e:
                    logging.error('%s on %s' %(e, array_id))
        if dump_config:
            self.timings = ts
            self.use = use
            self.dump(filename=pjoin(prefix, 'request.yaml'))



if __name__=="__main__":
    import argparse

    parser = argparse.ArgumentParser(description='download waveforms')
    parser.add_argument('--events', help='event file. uses only first event in list!')
    parser.add_argument('--site', help='geofon|iris', default='iris')
    parser.add_argument('--want', help='all|id1,id2,id3....', default='all')
    args = parser.parse_args()
    length = 1000.
    e = list(model.Event.load_catalog(args.events))[0]
    provider = DataProvider(channels='*Z')
    tmin = CakeTiming(phase_selection='first(p|P|PP)-100', fallback_time=0.001)
    tmax = CakeTiming(phase_selection='first(p|P|PP)+520', fallback_time=1000.)
    want = args.want
    if want!='all':
        want = want.split(',')
    provider.download(e, timing=(tmin, tmax), dump_config=True, force=True, site=args.site, want=want)
