import os
import shutil
from collections import defaultdict
from pyrocko.fdsn import ws
from pyrocko import model
from pyrocko import cake
from pyrocko import orthodrome as ortho
import logging
from util import create_directory

pjoin = os.path.join
logging.basicConfig(level='INFO')
logger = logging.getLogger('data-request')


class CakeTiming():
    '''Calculates and caches phase arrivals.
    :param fallback_time: returned, when no phase arrival was found for the
                        given depth-distance-phase-selection-combination

    E.g.:
    definition = 'first(p,P)-20'
    CakeTiming(definition)'''
    def __init__(self, phase_selection, fallback_time=None):
        self.arrivals = defaultdict(dict)
        self.fallback = fallback_time
        self.which = None
        if '+' in phase_selection:
            phase_selection, self.offset = phase_selection.split('+')
            self.offset = float(self.offset)
        elif '-' in phase_selection:
            phase_selection, self.offset = phase_selection.split('-')
            self.offset = float(self.offset)
            self.offset = -self.offset

        if 'first' in phase_selection:
            self.which = 'first'
        if 'last' in phase_selection:
            self.which = 'last'
        if self.which:
            phase_selection = self.strip(phase_selection)

        self.phases = phase_selection.split('|')


    def t(self, mod, z_dist):
        ''':param phase_selection: phase names speparated by vertical bars
        :param z_dist: tuple with (depth, distance)
        '''
        z, dist = z_dist
        if (dist, z) in self.arrivals.keys():
            return self.arrivals[(dist, z)]

        phases = [cake.PhaseDef(pid) for pid in self.phases]
        arrivals = mod.arrivals(distances=[dist*cake.m2d], phases=phases, zstart=z)
        if arrivals==[]:
            logger.warn('none of defined phases at d=%s, z=%s. (return fallback)'  % (dist, z))
            want = self.fallback
        else:
            want = self.phase_selector(arrivals)
            want = want.t + self.offset
        self.arrivals[(dist, z)] = want
        return want

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



class DataProvider():
    def __init__(self, channels='SHZ'):

        self.arrays = {'YKA': ('CN', 'YKA*', '', channels),
                       'ESK': [('IM', 'EKB?', '', channels),
                               ('IM', 'EKR*', '', channels)],
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


                       #'GERES': ()}

    def download(self, event, directory='array_data', timing=None, length=None,
                 want='all', force=False):
        """:param want: either 'all' or ID as string or list of IDs as strings
        """
        if all([timing, length]) == None:
            raise Exception('Define one of "timing" and "length"')

        if want=='all':
            wanted_ids = self.arrays.keys()
        elif isinstance(want, str):
            wanted_ids = want
        else:
            wanted_ids = want

        create_directory(directory, force)
        for array_id in wanted_ids:
            sub_directory = pjoin(directory, array_id)
            logger.info("fetching %s" % array_id)
            codes = self.arrays[array_id]
            if not isinstance(codes, list):
                codes = [codes]
            selection = [c + tuple((event.time, event.time+1000.)) for c in codes]
            logger.debug('selection: ', selection)
            try:
                st = ws.station(site='iris', selection=selection)
            except ws.EmptyResult as e:
                logging.error('%s on %s' %(e, array_id))

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
                d = ws.dataselect(site='iris', selection=selection)
                create_directory(sub_directory, force)
                fn = pjoin(sub_directory, 'traces.mseed')
                with open(fn, 'w') as f:
                    f.write(d.read())
                    f.close()
                model.dump_stations(stations, pjoin(sub_directory, 'stations.pf'))
            except ws.EmptyResult as e:
                logging.error('%s on %s' %(e, array_id))




if __name__=="__main__":
    import argparse

    parser = argparse.ArgumentParser(description='download waveforms')
    parser.add_argument('--events', help='event file. uses only first event in list!')
    args = parser.parse_args()
    length = 1000.
    e = list(model.Event.load_catalog(args.events))[0]
    provider = DataProvider()
    tmin = CakeTiming(phase_selection='first(p|P|PP)-40', fallback_time=0.)
    tmax = CakeTiming(phase_selection='first(p|P|PP)+40', fallback_time=1000.)
    provider.download(e, timing=(tmin, tmax))
