import numpy as num
from scipy.signal import detrend
import matplotlib
import math
matplotlib.use = 'QtAgg4'
import matplotlib.pyplot as plt
import logging
from pyrocko import io
from pyrocko import model, cake
from pyrocko.trace import IntegrationResponse, FrequencyResponse
from pyrocko.trace import ButterworthResponse, DifferentiationResponse
from pyrocko.gf import DCSource, Target, LocalEngine, seismosizer, meta
from pyrocko.util import str_to_time, match_nslc
from pyrocko.gui_util import load_markers
from pyrocko.guts import Object, Float, String, List, Bool
from pyrocko import orthodrome as ortho
km = 1000.
logging.basicConfig(loglevel="DEBUG")
logger = logging.getLogger('guesstimate')
arglist = ['station_filename', 'trace_filename', 'store_id', 'event_filename',
            'gain', 'correction', 'store_superdirs', 'depth', 'depths', 'zoom']
class PlotSettings(Object):
    trace_filename = String.T(help='filename of beam or trace to use for '
                              'plotting, incl. path.',
                              optional=True)
    station_filename = String.T(help='filename containing station meta '
                               'information related to *trace_filename*.',
                                optional=True)
    event_filename = String.T(help='filename containing event information '
                              'including the expected moment tensor.',
                              default='event.pf')
    store_id = String.T(help='Store ID to use for generating the synthetic '
                        'traces.', optional=True)
    store_superdirs = List.T(String.T(), optional=True)
    depth = Float.T(help='Depth [km] where to put the trace.', default=10.)
    depths = String.T(help='Synthetic source depths [km]. start:stop:delta. '
                      'default: 0:15:1', optional=True, default='0:15:1')
    filters = List.T(FrequencyResponse.T(help='List of filters used to filter '
                                         'the traces'))
    zoom = List.T(Float.T(), help='Window to visualize with reference to the P '
                    'phase onset [s].', default=[-7, 15])
    correction = Float.T(help='time shift, to move beam trace.', default=0.)
    normalize = Bool.T(help='normalize by maximum amplitude', default=True)
    save_as = String.T(default='depth_estimate.png', help='filename')
    force_nearest_neighbor = Bool.T(help='Handles OutOfBounds exceptions. '
                                  'applies only laterally!',
                                  default=False)
    auto_caption = Bool.T(help='add a caption giving basic information. Needs '
                               'implementation.',
                          default=False)
    title = String.T(default='%(array_id)s - %(event_name)s', help='Add default title.')
    quantity = String.T(default='velocity', help='velocity-> differentiate synthetic.'
                        'displacement-> integrate recorded')
    gain = Float.T(default=1., help='Gain factor')


    def update_from_args(self, args):
        kwargs = {}
        try:
            hp, lp = args.filter.split(':')
            filters = [
                ButterworthResponse(corner=float(lp), order=4, type='low'),
                ButterworthResponse(corner=float(hp), order=4, type='high')]
            kwargs.update({'filters': filters})
        except:
            pass

        for arg in arglist:
            try:
                val = getattr(args, arg)
                if val:
                    kwargs.update({arg: val})
            except AttributeError:
                logger.debug('%s not defined' % arg)
                continue
        for arg, v in kwargs.items():
            setattr(self, arg, v)

    @classmethod
    def from_argument_parser(cls, args):
        try:
            hp, lp = args.filter.split(':')
        except AttributeError:
            hp, lp = (0.7, 4.5)
        filters = [
            ButterworthResponse(corner=float(lp), order=4, type='low'),
            ButterworthResponse(corner=float(hp), order=4, type='high')]

        kwargs = {}
        for arg in arglist:
            try:
                val = getattr(args, arg)
                if val:
                    kwargs.update({arg: val})
            except AttributeError:
                logger.debug('%s not defined' % arg)
                continue

        return cls(filters=filters,
                   **kwargs)

    def do_filter(self, tr):
        for f in self.filters:
            tr = tr.transfer(transfer_function=f,
                             tfade=20,
                             cut_off_fading=False)
        return tr

class GaussNotch(FrequencyResponse):
    def __init__(self, center, fwhm):
        self.center = center
        self.fwhm = fwhm

    def evaluate(self, freqs):
        denom = self.fwhm / (2.*math.sqrt(math.log(2.)))
        return 1.-num.exp(1-((freqs-self.center)/denom)**2)

def integrate_differentiate(tr, do):
    if do == 'integrate':
        response = IntegrationResponse()
    elif do == 'differentiate':
        response = DifferentiationResponse()
    tr = tr.transfer(transfer_function=response,
                     tfade=20,
                     cut_off_fading=True,
                     freqlimits=(0.10, 0.2, 10., 20.))
    return tr

def notch_filter(tr, f_center, bandwidth):
    indi = num.arange(tr.data_len(), dtype=num.float)
    tr.set_ydata(detrend(tr.get_ydata()))
    freqs, FT_data = tr.spectrum(pad_to_pow2=True, tfade=None)
    fft_length = FT_data.size
    FFF = num.ones((fft_length), dtype=num.complex) * GaussNotch(f_center, bandwidth).evaluate(freqs)
    filtered_data = num.fft.irfft(FT_data*FFF)[:tr.data_len()]
    tr.set_ydata(filtered_data)
    return tr

def station_to_target(s, quantity, store_id):
    return Target(codes=s.nsl()+tuple('Z'),
                 lat=s.lat,
                 lon=s.lon,
                 elevation=s.elevation,
                 quantity=quantity,
                 store_id=store_id)


def plot(settings, show=False):

    #align_phase = 'P(cmb)P<(icb)(cmb)p'
    align_phase = 'P'
    zoom_window = settings.zoom

    print zoom_window
    ampl_scaler = '4*standard deviation'

    quantity = settings.quantity
    zstart, zstop, inkr = settings.depths.split(':')
    test_depths = num.arange(float(zstart)*km, float(zstop)*km, float(inkr)*km)

    traces = io.load(settings.trace_filename)

    event = model.load_events(settings.event_filename)
    assert len(event)==1
    event = event[0]
    event.depth = float(settings.depth) * 1000.
    base_source = DCSource.from_pyrocko_event(event)

    test_sources = []
    for d in test_depths:
        s = base_source.clone()
        s.depth = float(d)
        test_sources.append(s)

    if settings.store_superdirs:
        engine = LocalEngine(store_superdirs=settings.store_superdirs)
    else:
        engine = LocalEngine(use_config=True)
    try:
        store = engine.get_store(settings.store_id)
    except seismosizer.NoSuchStore as e:
        logger.info('%s ... skipping.' % e)
        return

    stations = model.load_stations(settings.station_filename)
    station = filter(lambda s: match_nslc('%s.%s.%s.*' % s.nsl(), traces[0].nslc_id), stations)
    assert len(station) == 1
    station = station[0] 
    targets = [station_to_target(station, quantity=quantity, store_id=settings.store_id)]
    try:
        request = engine.process(targets=targets, sources=test_sources)
    except seismosizer.NoSuchStore as e:
        logger.info('%s ... skipping.' % e)
        return
    except meta.OutOfBounds as error:
        if settings.force_nearest_neighbor:
            logger.warning('%s  Using nearest neighbor instead.' % error)
            mod_targets = []
            for t in targets:
                closest_source = min(test_sources, key=lambda s: s.distance_to(t))
                farthest_source = max(test_sources, key=lambda s: s.distance_to(t))
                min_dist_delta = store.config.distance_min - closest_source.distance_to(t)
                max_dist_delta = store.config.distance_max - farthest_source.distance_to(t)
                if min_dist_delta < 0:
                    azi, bazi = closest_source.azibazi_to(t)
                    newlat, newlon = ortho.azidist_to_latlon(t.lat, t.lon, azi, min_dist_delta*cake.m2d)
                elif max_dist_delta < 0:
                    azi, bazi = farthest_source.azibazi_to(t)
                    newlat, newlon = ortho.azidist_to_latlon(t.lat, t.lon, azi, max_dist_delta*cake.m2d)
                t.lat, t.lon = newlat, newlon
                mod_targets.append(t)
            request = engine.process(targets=mod_targets, sources=test_sources)
        else:
            raise error

    alldepths = list(test_depths)
    depth_count = dict(zip(sorted(alldepths), range(len(alldepths))))

    target_count = dict(zip([t.codes[:3] for t in targets], range(len(targets))))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    maxz = max(test_depths)
    minz = min(test_depths)
    relative_scale = (maxz-minz)*0.02
    for s, t, tr in request.iter_results():
        if quantity=='velocity':
            tr = integrate_differentiate(tr, 'differentiate')

        #ax = axs[target_count[t.codes[:3]]][depth_count[s.depth]]
        onset = engine.get_store(t.store_id).t(
            'begin', (s.depth, s.distance_to(t)))

        tr = settings.do_filter(tr)
        if settings.normalize:
            tr.set_ydata(tr.get_ydata()/num.max(abs(tr.get_ydata())))
            ax.tick_params(axis='y', which='both', left='off', right='off',
                           labelleft='off')

        y_pos = s.depth
        xdata = tr.get_xdata()-onset-s.time
        tr_ydata = tr.get_ydata()
        if ampl_scaler == 'trace min/max':
            ampl_scale = float(max(abs(tr_ydata)))
        elif ampl_scaler == '4*standard deviation':
            ampl_scale = 4*float(num.std(tr_ydata))
        else:
            ampl_scale = 1.
        ydata = (tr_ydata/ampl_scale * settings.gain)*relative_scale + y_pos
        ax.plot(xdata, ydata, c='blue', linewidth=1., alpha=0.5)
        #if self.fill_between:
        #    ax.fill_between(xdata, y_pos, ydata, where=ydata>y_pos, color='black', alpha=0.5)
        ax.text(zoom_window[0], y_pos, '%s km' % (s.depth/1000.), horizontalalignment='right', fontsize=12.)
        if True:
            mod = store.config.earthmodel_1d
            label = 'pP'
            arrivals = mod.arrivals(phases=[cake.PhaseDef(label)],
                                      distances=[s.distance_to(t)*cake.m2d],
                                      zstart=s.depth)

            try:
                t = arrivals[0].t
                ydata_absmax = num.max(num.abs(tr.get_ydata()))
                marker_length = 0.5
                x_marker = [t-onset]*2
                y = [y_pos-(maxz-minz)*0.025, y_pos+(maxz-minz)*0.025]
                ax.plot(x_marker, y, linewidth=1, c='blue')

                ax.text(x_marker[1]-x_marker[1]*0.005, y[1], label,
                        fontsize=12,
                        color='blue',
                        verticalalignment='top',
                        horizontalalignment='right')

            except IndexError:
                logger.warning('no pP phase at d=%s z=%s stat=%s' % (s.distance_to(t)*cake.m2d,
                                                                     s.depth, station.station))
                pass

    if len(traces)==0:
        raise Exception('No Trace found!')
    if len(traces)>1:
        raise Exception('More then one trace provided!')
    else:
        onset = 0
        tr = traces[0]
        correction = float(settings.correction)
        if quantity=='displacement':
            tr = integrate_differentiate(tr, 'integrate')
        tr = settings.do_filter(tr)
        onset = engine.get_store(targets[0].store_id).t(
            'begin', (event.depth, s.distance_to(targets[0]))) + event.time
        if settings.normalize:
            tr.set_ydata(tr.get_ydata()/max(abs(tr.get_ydata())))
            ax.tick_params(axis='y', which='both', left='off', right='off',
                           labelleft='off')

        y_pos = event.depth
        xdata = tr.get_xdata()-onset
        tr_ydata = tr.get_ydata()
        if ampl_scaler == 'trace min/max':
            ampl_scale = float(max(abs(tr_ydata)))
        elif ampl_scaler == '4*standard deviation':
            ampl_scale = 4*float(num.std(tr_ydata))
        else:
            ampl_scale = 1.
        ydata = (tr_ydata/ampl_scale * settings.gain)*relative_scale + y_pos
        ax.plot(xdata, ydata, c='red', linewidth=1.)
        ax.set_xlim(zoom_window)
    if settings.title:
        params = {'array_id': '.'.join(station.nsl()), 'event_name':event.name}
        fig.suptitle(settings.title % params)
    plt.subplots_adjust(hspace=-.4, bottom=0.1)
    if settings.save_as:
        logger.info('save as: %s ' % settings.save_as)
        fig.savefig(settings.save_as)
    if show:
        plt.show()

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Find depth.')
    parser.add_argument('--trace',
                        help='name of file containing trace',
                        dest='trace_filename',
                        required=True)
    parser.add_argument('--station',
                        dest='station_filename',
                        help='name of file containing station information',
                        required=True)
    parser.add_argument('--event',
                        help='name of file containing event catalog',
                        dest='event_filename',
                        required=True)
    parser.add_argument('--store',
                        dest='store_id',
                        help='name of store id',
                        required=True)
    parser.add_argument('--store-superdir',
                        help='store superdir',
                        nargs='*', dest='store_superdirs',
                        required=False)
    parser.add_argument('--pick',
                        help='name of file containing p marker',
                        required=False)
    parser.add_argument('--depth',
                        help='assumed source depth [km]',
                        default=10.,
                        required=False)
    parser.add_argument('--depths',
                        help='testing depths in km. zstart:zstop:delta',
                        required=True)
    parser.add_argument('--quantity',
                        help='velocity|displacement',
                        default='velocity',
                        required=False)
    parser.add_argument('--filter',
                        help='4th order butterw. default: "0.7:4.5"',
                        required=False)
    parser.add_argument('--correction',
                        help='correction in time [s]',
                        dest='correction',
                        default=0,
                        required=False)
    parser.add_argument('--normalize',
                        help='normalize traces to 1',
                        action='store_true',
                        required=False)
    parser.add_argument('--skip_true',
                        help='if true, do not plot recorded and the assigned synthetic trace on top of each other',
                        action='store_true',
                        required=False)
    parser.add_argument('--save-as',
                        dest='save_as',
                        help='file to store image',
                        required=False)
    parser.add_argument('--print_parameters',
                        help='creates a text field giving the used parameters',
                        required=False)
    parser.add_argument('--title',
                        help='title of figure',
                        required=False)
    parser.add_argument('--no-y-axis',
                        help='Do not plot depths',
                        dest='no_y_axis',
                        action='store_true',
                        required=False)
    parser.add_argument('--show',
                        help='Show results right away',
                        action='store_true',
                        required=False)
    parser.add_argument('--settings',
                        help='filename defining settings. Parameters defined '
                        'defined parameters will overwrite those.',
                        required=False)
    parser.add_argument('--gain', help='Gain factor', type=float, required=False)
    parser.add_argument('--zoom', help='Zoom window like t1:t2', nargs=2, type=float, required=False)

    args = parser.parse_args()
    settings = PlotSettings.from_argument_parser(args)
    plot(settings, show=args.show)
