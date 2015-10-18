import numpy as num
from scipy.signal import detrend
import matplotlib
import math
matplotlib.use = 'QtAgg4'
import matplotlib.pyplot as plt
from pyrocko import io
from pyrocko import model, cake
from pyrocko.trace import IntegrationResponse, FrequencyResponse
from pyrocko.trace import ButterworthResponse, DifferentiationResponse
from pyrocko.gf import DCSource, Target, LocalEngine
from pyrocko.util import str_to_time
from pyrocko.gui_util import load_markers
from pyrocko.guts import Object, Float, String, List, Bool
km = 1000.


import argparse


class PlotSettings(Object):
    trace_filename = String.T(help='filename of beam or trace to use for '
                              'plotting, incl. path.')
    station_filename = String.T(help='filename containing station meta '
                               'information related to *trace_filename*.')
    event_filename = String.T(help='filename containing event information '
                              'including the expected moment tensor.')
    store_id = String.T(help='Store ID to use for generating the synthetic '
                        'traces.')
    depth = Float.T(help='Depth [km] where to put the trace.')
    depths = String.T(help='Synthetic source depths [km]. start:stop:delta. '
                      'default: 0:15:1',
                      optional=True, default='0:15:1')
    filters = List.T(FrequencyResponse.T(help='List of filters used to filter '
                                         'the traces'))
    zoom = List.T(Float.T(), help='Window to visualize with reference to the P '
                  'phase onset [s].', default=[-7, 15])
    onset_correction = Float.T(help='time shift, to move beam trace.', default=0.)
    normalization = Bool.T(help='normalize by maximum amplitude', default=True)
    # do I need this, here?
    save_as = String.T(default='depth_estimate.png', help='filename')
    auto_caption = Bool.T(help='add a caption giving basic information. Needs '
                               'implementation.')
    title = String.T(default='%{array_id}s - %{event_name}s', help='Add default title.')

    @classmethod
    def from_arguement_parser(cls, args):
        hp, lp = args.filter.split(':')
        filters = [
            ButterworthResponse(corner=float(lp), order=4, type='low'),
            ButterworthResponse(corner=float(hp), order=4, type='high')]
        cls(trace_filename=args.trace,
            station_filename=args.station,
            event_filename=args.event,
            store_id=args.store,
            depth=args.depth,
            depths=args.depths,
            filters=filters)
        raise Exception('needs implementation')


class GaussNotch(FrequencyResponse):
    def __init__(self, center, fwhm):
        self.center = center
        self.fwhm = fwhm

    def evaluate(self, freqs):
        denom = self.fwhm / (2.*math.sqrt(math.log(2.)))
        return 1.-num.exp(1-((freqs-self.center)/denom)**2)

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

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Find depth.')
    parser.add_argument('--trace',
                        help='name of file containing trace',
                        required=True)
    parser.add_argument('--station',
                        help='name of file containing station information',
                        required=True)
    parser.add_argument('--event',
                        help='name of file containing event catalog',
                        required=True)
    parser.add_argument('--store',
                        help='name of store id',
                        required=True)
    parser.add_argument('--pick',
                        help='name of file containing p marker',
                        required=True)
    parser.add_argument('--depth',
                        help='assumed source depth [km]',
                        default=10.,
                        required=False)
    parser.add_argument('--depths',
                        help='testing depths in km. zstart:zstop:delta',
                        default=0,
                        required=False)
    parser.add_argument('--quantity',
                        help='velocity|displacement',
                        default='velocity',
                        required=False)
    parser.add_argument('--filter',
                        help='4th order butterw. default: "0.7:4.5"',
                        default="0.7:4.5",
                        required=False)
    parser.add_argument('--correction',
                        help='correction in time [s]',
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
    parser.add_argument('--out_filename',
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
    parser.add_argument('--settings',
                        help='filename defining settings. Parameters defined '
                        'defined parameters will overwrite those.',
                        required=False)

    args = parser.parse_args()

    settings = PlotSettings().from_arguement_parser(args)


def plot(settings):

    # use test depth:
    # sind das die richtigen strike dip rake kobinationen?
    # zeitbereich, den man betrachten moechte relativ zur p phase
    zoom_window = settings.zoom
    # zoom_window = [-7, 15]

    #notch = 0.15 # use this for GERES
    #notch = False
    quantity = args.quantity
    lp, hp = args.filter.split(':')
    bandpass = {'order': 4, 'corner_hp': float(lp), 'corner_lp': float(hp) }

    if args.depths:
        zstart, zstop, inkr = args.depths.split(':')
        test_depths = num.arange(float(zstart)*km, float(zstop)*km, float(inkr)*km)
    else:
        test_depths = num.arange(0*km, 40*km, 1*km)

    traces = io.load(args.trace)

    event = model.load_events(args.event)[0]
    event.depth = float(args.sdepth) * 1000.
    base_source = DCSource.from_pyrocko_event(event)

    test_sources = []
    # setup sources:
    for d in test_depths:
        #if d==base_source.depth:
        #    continue
        s = base_source.clone()
        s.depth = float(d)
        test_sources.append(s)

    e = LocalEngine(store_superdirs=['/media/usb/stores'])
    station = model.load_stations(args.station)
    targets = [station_to_target(station[0], quantity=args.quantity, store_id=args.store)]
    print targets[0]
    print 'distance: ', targets[0].distance_to(base_source)
    request = e.process(targets=targets,
                         sources=test_sources)
    alldepths = list(test_depths)
    #alldepths.append(base_source.depth)
    depth_count = dict(zip(sorted(alldepths), range(len(alldepths))))
    print depth_count

    target_count = dict(zip([t.codes[:3] for t in targets], range(len(targets))))

    fig, axs = plt.subplots(len(test_sources),
                            len(targets),
                            sharex=True)
                            #frameon=True)
    if len(targets)==1:
        axs = [axs]

    for s, t, tr in request.iter_results():
        if s.depth == base_source.depth and args.skip_true:
            continue
        tr.bandpass(**bandpass)
        if notch:
            notch_filter(tr, 2*num.pi*notch, 1.5)
        if quantity=='velocity':
            diff_response = DifferentiationResponse()
            tr = tr.transfer(transfer_function=diff_response, tfade=20, freqlimits=(0.1, 0.2, 10., 20.))

        ax = axs[target_count[t.codes[:3]]][depth_count[s.depth]]
        onset = e.get_store(t.store_id).t('P', (s.depth, s.distance_to(t)))
        store = e.get_store(args.store)
        mod = store.config.earthmodel_1d
        #onset = mod.arrivals(phases=[cake.PhaseDef('P')], 
        #                              distances=[s.distance_to(t)*cake.m2d],
        #                              zstart=s.depth)[0].t
        ydata = tr.get_ydata()
        if args.normalize:
            ydata = ydata/num.max(abs(ydata))
            ax.tick_params(axis='y',
                           which='both',
                           left='off',
                           right='off',
                           labelleft='off')

        ax.plot(tr.get_xdata()-s.time-onset, ydata)
        ax.set_xlim(zoom_window)
        if not args.no_y_axis:
            ax.text(-0.01, 0.5,'%s km' % (s.depth/1000.),
                    transform=ax.transAxes,
                    horizontalalignment='right')
        ax.axes.patch.set_visible(False)
        if True:
            arrivals = mod.arrivals(phases=[cake.PhaseDef('pP')], 
                                      distances=[s.distance_to(t)*cake.m2d],
                                      zstart=s.depth)

            t = arrivals[0].t
            ydata_absmax = num.max(num.abs(ydata))
            marker_length = 0.5
            ax.plot([t-onset]*2,
                    [-ydata_absmax*marker_length, ydata_absmax*marker_length],
                    linewidth=1, 
                    c='red')
            print arrivals
        if s.depth==max(test_depths):
            ax.xaxis.set_ticks_position('bottom')
            for pos in ['left', 'top','right']:
                ax.spines[pos].set_visible(False)
            ax.set_xlabel('Time [s]')
        else:
            ax.axes.get_xaxis().set_visible(False)
            for item in ax.spines.values():
                item.set_visible(False)

    for tr in traces:
        try:
            correction = float(args.correction)
        except KeyError:
            correction = 0
        marker = load_markers(args.pick)[0]
        ponset = marker.tmin

        if quantity=='displacement':
            integration_response = IntegrationResponse()
            tr = tr.transfer(transfer_function=integration_response,
                             tfade=20,
                             freqlimits=(0.10, 0.2, 10., 20.))
        if notch:
            notch_filter(tr, 2*num.pi*notch, 1.5)
        tr.bandpass(**bandpass)
        ax = axs[target_count[tr.nslc_id[:3]]][depth_count[base_source.depth]]
        ydata = tr.get_ydata()
        if args.normalize:
            ydata = ydata/max(abs(ydata))
            ax.tick_params(axis='y',
                           which='both',
                           left='off',
                           right='off',
                           labelleft='off')

        ax.plot(tr.get_xdata()-ponset+correction, ydata, c='black', linewidth=2)
        if not args.no_y_axis:
            ax.text(-0.01, 0.5,'%s km' % (base_source.depth/1000.),
                    transform=ax.transAxes,
                    horizontalalignment='right')
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.patch.set_visible(False)
        for item in ax.spines.values():
            item.set_visible(False)
        ax.set_xlim(zoom_window)
    if args.title:
        fig.suptitle(args.title)
    bottom = 0.1
    plt.subplots_adjust(hspace=-.4,
                        bottom=bottom)
    if args.out_filename:
        fig.savefig(args.out_filename)
    plt.show()

