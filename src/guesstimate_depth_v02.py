import numpy as num
from scipy.signal import detrend
import matplotlib

font = {'family' : 'normal',
        'size'   : 16}
matplotlib.rc('font', **font)
import math
matplotlib.use = 'QtAgg4'
import matplotlib.pyplot as plt
import logging
from pyrocko import io
from pyrocko import model, cake
from pyrocko.trace import IntegrationResponse, FrequencyResponse
from pyrocko.trace import ButterworthResponse, DifferentiationResponse
from pyrocko.gf import DCSource, Target, LocalEngine, seismosizer, meta
from pyrocko.util import str_to_time, time_to_str, match_nslc
from pyrocko.gui_util import load_markers
from pyrocko.guts import Object, Float, String, List, Bool
from pyrocko import orthodrome as ortho
km = 1000.
logging.basicConfig(loglevel="DEBUG")
logger = logging.getLogger('guesstimate')
arglist = ['station_filename', 'trace_filename', 'store_id', 'event_filename',
            'gain', 'gain_record', 'correction', 'store_superdirs', 'depth', 'depths', 'zoom',
           'title', 'save_as', 'color', 'auto_caption', 'quantity']
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
    save_as = String.T(default='depth_%(array_id)s.png', help='filename')
    force_nearest_neighbor = Bool.T(help='Handles OutOfBounds exceptions. '
                        'applies only laterally!', default=False)
    auto_caption = Bool.T(help='Add a caption giving basic information.',
                          default=False)
    title = String.T(default='%(array_id)s - %(event_name)s', help='Add default title.')
    quantity = String.T(default='velocity', help='velocity-> differentiate synthetic.'
                        'displacement-> integrate recorded')
    gain = Float.T(default=1., help='Gain factor')
    gain_record = Float.T(default=1., help='Gain factor')
    color = String.T(help='Trace color', default='blue')

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
                     freqlimits=(0.05, 0.1, 10., 20.))
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
    with_onset_line = True
    fill = True
    align_phase = 'P'
    zoom_window = settings.zoom
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

        onset = engine.get_store(t.store_id).t(
            'begin', (s.depth, s.distance_to(t)))

        tr = settings.do_filter(tr)
        if settings.normalize:
            tr.set_ydata(tr.get_ydata()/num.max(abs(tr.get_ydata())))
            ax.tick_params(axis='y', which='both', left='off', right='off',
                           labelleft='off')

        y_pos = s.depth
        xdata = tr.get_xdata()-onset-s.time
        tr_ydata = tr.get_ydata() * -1
        visible = tr.chop(tmin=event.time+onset+zoom_window[0],
                          tmax=event.time+onset+zoom_window[1])
        if ampl_scaler == 'trace min/max':
            ampl_scale = float(max(abs(visible.get_ydata())))
        elif ampl_scaler == '4*standard deviation':
            ampl_scale = 4*float(num.std(visible.get_ydata()))
        else:
            ampl_scale = 1.
        ampl_scale /= settings.gain
        ydata = (tr_ydata/ampl_scale)*relative_scale + y_pos
        ax.plot(xdata, ydata, c='black', linewidth=1., alpha=1.)
        if False:
            ax.fill_between(xdata, y_pos, ydata, where=ydata<y_pos, color='black', alpha=0.5)
        ax.text(zoom_window[0]*1.09, y_pos, '%i' % (s.depth/1000.), horizontalalignment='right') #, fontsize=12.)
        if False:
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
                        #fontsize=12,
                        color='black',
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
        xdata = tr.get_xdata()-onset+correction
        tr_ydata = tr.get_ydata() *-1
        visible = tr.chop(tmin=onset+zoom_window[0]+correction,
                          tmax=onset+zoom_window[1]+correction)
        if ampl_scaler == 'trace min/max':
            ampl_scale = float(max(abs(visible.get_ydata())))
        elif ampl_scaler == '4*standard deviation':
            ampl_scale = 4*float(num.std(visible.get_ydata()))
        else:
            ampl_scale = 1.
        ydata = (tr_ydata/ampl_scale * settings.gain*settings.gain_record)*relative_scale + y_pos
        ax.plot(xdata, ydata, c=settings.color, linewidth=1.)
        ax.set_xlim(zoom_window)
        zmax = max(test_depths)
        zmin = min(test_depths)
        zrange = zmax - zmin
        ax.set_ylim((zmin-zrange*0.2, zmax+zrange*0.2))
        ax.set_xlabel('Time [s]')
        ax.text(-0.08, 0.6, 'Source depth [km]',
                rotation=90,
                horizontalalignment='right',
                transform=ax.transAxes) #, fontsize=12.)

    if fill:
        ax.fill_between(xdata, y_pos, ydata, where=ydata<y_pos, color=settings.color, alpha=0.5)
    if with_onset_line:
        ax.text(0.08, zmax+zrange*0.1, align_phase, fontsize=14)
        vline = ax.axvline(0., c='black')
        vline.set_linestyle('--')
    if settings.title:
        params = {'array_id': '.'.join(station.nsl()),
                  'event_name': event.name,
                  'event_time': time_to_str(event.time)}
        ax.text(0.5, 1.05, settings.title % params,
                horizontalalignment='center', 
                transform=ax.transAxes)
    if settings.auto_caption:
        cax = fig.add_axes([0., 0., 1, 0.05], label='caption')
        cax.axis('off')
        cax.xaxis.set_visible(False)
        cax.yaxis.set_visible(False)
        captions = {'filters':''}
        for f in settings.filters:
            captions['filters'] += '%s pass, order %s, f$_c$=%s Hz, '%(f.type, f.order, f.corner)
        captions['store_sampling'] = 1./store.config.deltat
        cax.text(0, 0, 'Filters: %(filters)s GFDB sampled at %(store_sampling)s Hz.' % captions,
                 fontsize=12, transform=cax.transAxes)
        plt.subplots_adjust(hspace=.4, bottom=0.15)
    else:
        plt.subplots_adjust(bottom=0.1)

    ax.invert_yaxis()
    if settings.save_as:
        logger.info('save as: %s ' % settings.save_as)
        fig.savefig(settings.save_as%{'array_id': '.'.join(station.nsl())}, dpi=160, bbox_inches='tight')
    if show:
        plt.show()


class Inverter():
    def __init__(self, provider, args):
        self.provider = provider
        self.args = args

    def invert(self, args):
        align_phase = 'P'
        ampl_scaler = '4*standard deviation'

        for array_id in self.provider.use:
            try:
                if args.array_id and array_id != args.array_id:
                    continue
            except AttributeError:
                pass
            subdir = pjoin('array_data', array_id)
            settings_fn = pjoin(subdir, 'plot_settings.yaml')
            if os.path.isfile(settings_fn):
                settings = PlotSettings.load(filename=pjoin(settings_fn))
                settings.update_from_args(self.args)
            else:
                logger.warn('no settings found: %s' % array_id)
                continue
            if settings.store_superdirs:
                engine = LocalEngine(store_superdirs=settings.store_superdirs)
            else:
                engine = LocalEngine(use_config=True)
            try:
                store = engine.get_store(settings.store_id)
            except seismosizer.NoSuchStore as e:
                logger.info('%s ... skipping.' % e)
                return
            try:
                store = engine.get_store(settings.store_id)
            except seismosizer.NoSuchStore as e:
                logger.info('%s ... skipping.' % e)
                return

            if not settings.trace_filename:
                settings.trace_filename = pjoin(subdir, 'beam.mseed')
            if not settings.station_filename:
                settings.station_filename = pjoin(subdir, 'array_center.pf')
            zoom_window = settings.zoom
            mod = store.config.earthmodel_1d

            zstart, zstop, inkr = settings.depths.split(':')
            test_depths = num.arange(float(zstart)*km, float(zstop)*km, float(inkr)*km)
            traces = io.load(settings.trace_filename)
            event = model.load_events(settings.event_filename)
            assert len(event)==1
            event = event[0]
            event.depth = float(settings.depth) * 1000.
            base_source = MTSource.from_pyrocko_event(event)

            test_sources = []
            for d in test_depths:
                s = base_source.clone()
                s.depth = float(d)
                test_sources.append(s)

            stations = model.load_stations(settings.station_filename)
            station = filter(lambda s: match_nslc('%s.%s.%s.*' % s.nsl(), traces[0].nslc_id), stations)
            if len(station) != 1:
                logger.error('no matching stations found. %s %s' % []) 
            else:
                station = station[0]
            targets = [station_to_target(station, quantity=settings.quantity, store_id=settings.store_id)]
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

            candidates = []
            for s, t, tr in request.iter_results():
                tr.deltat = regularize_float(tr.deltat)
                if True:
                    tr = integrate_differentiate(tr, 'differentiate')
                tr = settings.do_filter(tr)
                candidates.append((s, tr))
            assert len(traces)==1
            ref = traces[0]
            ref = settings.do_filter(ref)
            dist = ortho.distance_accurate50m(event, station)
            tstart = self.provider.timings[array_id].timings[0].t(mod, (event.depth, dist)) + event.time
            tend = self.provider.timings[array_id].timings[1].t(mod, (event.depth, dist)) + event.time
            ref = ref.chop(tstart, tend)
            misfits = []

            center_freqs = num.arange(1., 9., 4.)
            num_f_widths = len(center_freqs)

            mesh_fc = num.zeros(len(center_freqs)*num_f_widths*len(candidates))
            mesh_fwidth = num.zeros(len(center_freqs)*num_f_widths*len(candidates))
            misfits_array = num.zeros((len(center_freqs), num_f_widths, len(candidates)))
            depths_array = num.zeros((len(center_freqs), num_f_widths, len(candidates)))
            debug = False
            pb = ProgressBar(maxval=max(center_freqs)).start()
            i = 0
            for i_fc, fc in enumerate(center_freqs):
                if debug:
                    fig = plt.figure()

                fl_min = fc-fc*2./5.
                fr_max = fc+fc*2./5.
                widths = num.linspace(fl_min, fr_max, num_f_widths)

                for i_width, width in enumerate(widths):
                    i_candidate = 0
                    mesh_fc[i] = fc
                    mesh_fwidth[i] = width
                    i += 1
                    for source, candidate in candidates:
                        candidate = candidate.copy()
                        tstart = self.provider.timings[array_id].timings[0].t(mod, (source.depth, dist)) + event.time
                        tend = self.provider.timings[array_id].timings[1].t(mod, (source.depth, dist)) + event.time
                        filters = [
                            ButterworthResponse(corner=float(fc+width*0.5), order=4, type='low'),
                            ButterworthResponse(corner=float(fc-width*0.5), order=4, type='high')]
                        settings.filters = filters
                        candidate = settings.do_filter(candidate)
                        candidate.chop(tmin=tstart, tmax=tend)
                        candidate.shift(float(settings.correction))
                        m, n, aproc, bproc = ref.misfit(candidate=candidate, setup=settings.misfit_setup, debug=True)
                        aproc.set_codes(station='aproc')
                        bproc.set_codes(station='bproc')
                        if debug:
                            ax = fig.add_subplot(len(test_depths)+1, 1, i+1)
                            ax.plot(aproc.get_xdata(), aproc.get_ydata())
                            ax.plot(bproc.get_xdata(), bproc.get_ydata())
                        mf = m/n
                        #misfits.append((source.depth, mf))
                        misfits_array[i_fc][i_width][i_candidate] = mf
                        i_candidate += 1
                pb.update(fc)

            pb.finish()
            fig = plt.figure()
            ax = fig.add_subplot(111)
            i_best_fits = num.argmin(misfits_array, 2)
            print 'best fits: \n', i_best_fits
            best_fits = num.min(misfits_array, 2)
            #cmap = matplotlib.cm.get_cmap()
            xmesh, ymesh = num.meshgrid(mesh_fc, mesh_fwidth)
            #c = (best_fits-num.min(best_fits))/(num.max(best_fits)-num.min(best_fits))
            ax.scatter(xmesh, ymesh, best_fits*100)
            #ax.scatter(mesh_fc, mesh_fwidth, c)
            #ax.scatter(mesh_fc, mesh_fwidth, s=best_fits)
            ax.set_xlabel('fc')
            ax.set_ylabel('f_width')
        plt.legend()
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
    parser.add_argument('--auto-caption', help='caption basic info',
                        action='store_true', required=False)

    parser.add_argument('--gain', help='Gain factor', type=float, required=False)
    parser.add_argument('--gain-record', dest='gain_record', help='Gain factor', type=float, required=False)
    parser.add_argument('--zoom', help='Zoom window like t1:t2', nargs=2, type=float, required=False)
    parser.add_argument('--color', help='color of trace', required=False)

    args = parser.parse_args()
    settings = PlotSettings.from_argument_parser(args)
    plot(settings, show=args.show)
