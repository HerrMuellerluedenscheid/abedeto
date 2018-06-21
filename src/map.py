import numpy as num 
import math
from pyrocko import automap
from pyrocko import model, gmtpy
from pyrocko.guts import Object, Float, List, String, Bool, Dict, Tuple, Int
from pyrocko.model import Event
from pyrocko.moment_tensor import magnitude_to_moment as mag2mom
from pyrocko.moment_tensor import to6


size_factor = 4.
cm = 10.


def uniso(m):
    xx = num.trace(m) / 3.
    mc = num.matrix([[xx, 0., 0.],[0., xx, 0.], [0., 0., xx]])
    return m - mc


class MapParameters(Object):
    width = Float.T(help='width of map', optional=True, default=16)
    height = Float.T(help='heigth of map', optional=True, default=16)
    lat = Float.T(help='center latitude')
    lon = Float.T(help='center longitude')
    radius = Float.T(help='radius of map')
    outfn = String.T(help='output filename', optional=True, default='map.png')
    stations = List.T(model.Station.T(), optional=True)
    station_label_mapping = List.T(String.T(), optional=True)
    station_colors = Dict.T(String.T(), String.T(), optional=True)
    events = List.T(model.Event.T(), optional=True)
    show_topo = Bool.T(optional=True, default=True)
    show_grid = Bool.T(optional=True, default=True)
    color_wet = Tuple.T(3, Int.T(), default=(200, 200, 200))
    color_dry = Tuple.T(3, Int.T(), default=(253, 253, 253))

    @classmethod
    def process_args(cls, *args, **kwargs):
        return cls(**kwargs)


class JMap(automap.Map):
    def __init__(self, *args, **kwargs):
        automap.Map.__init__(self, *args, **kwargs)

    def _setup_gmt(self):
            w, h = self.width, self.height
            scaler = self._scaler

            gmtconf = dict(
                TICK_PEN='1.25p',
                TICK_LENGTH='0.2c',
                ANNOT_FONT_PRIMARY='1',
                ANNOT_FONT_SIZE_PRIMARY='12p',
                LABEL_FONT='1',
                LABEL_FONT_SIZE='12p',
                CHAR_ENCODING='ISOLatin1+',
                BASEMAP_TYPE='fancy',
                PLOT_DEGREE_FORMAT='D',
                PAPER_MEDIA='Custom_%ix%i' % (
                    w*gmtpy.cm,
                    h*gmtpy.cm),
                GRID_PEN_PRIMARY='thinnest/0/50/0',
                DOTS_PR_INCH='1200',
                OBLIQUE_ANNOTATION='6')

            gmtconf.update(
                (k.upper(), v) for (k, v) in self.gmt_config.iteritems())

            gmt = gmtpy.GMT(config=gmtconf)

            layout = gmt.default_layout()
            layout.set_fixed_margins(*[x*cm for x in self._expand_margins()])

            widget = layout.get_widget()
            widget['P'] = widget['J']
            widget['J'] = ('-JE%g/%g' % (self.lon, self.lat)) + '/%(width)gp'
            scaler['R'] = '-R%g/%g/%g/%gr' % self._corners

            aspect = gmtpy.aspect_for_projection(*(widget.J() + scaler.R()))
            widget.set_aspect(aspect)

            self._gmt = gmt
            self._layout = layout
            self._widget = widget
            self._jxyr = self._widget.JXY() + self._scaler.R()
            self._pxyr = self._widget.PXY() + [
                '-R%g/%g/%g/%g' % (0, widget.width(), 0, widget.height())]
            self._have_drawn_axes = False
            self._have_drawn_labels = False


def make_map(lat=None, lon=None, radius=None, outfn=None, stations=None, events=None, stations_label_mapping=None, map_parameters=None, show_topo=True):
    if map_parameters:
        lat = map_parameters.lat
        lon = map_parameters.lon
        radius=map_parameters.radius
        outfn = map_parameters.outfn
        stations = map_parameters.stations
        station_label_mapping = map_parameters.station_label_mapping
        events = map_parameters.events
        station_colors = map_parameters.station_colors
        show_topo = map_parameters.show_topo
        color_wet = map_parameters.color_wet
        color_dry = map_parameters.color_dry
        show_grid = map_parameters.show_grid
        height = map_parameters.height
        width = map_parameters.width

    _map = JMap(width=width, height=height, lat=lat, lon=lon,
        radius=radius,
        topo_resolution_max=200,
        topo_resolution_min=40.,
        show_topo=show_topo,
        show_grid=show_grid,
        illuminate=True,
        illuminate_factor_land=0.5,
        illuminate_factor_ocean=0.25,
        color_wet=color_wet,
        color_dry=color_dry)

    #_map.draw_cities()
    if stations:
        for i_station, s in enumerate(stations):
            lats = [s.lat]
            lons = [s.lon]
            if station_colors:
                color = "%s" % station_colors[s.nsl()]
            else:
                color = "black"

            _map.gmt.psxy(
                in_columns=(lons, lats),
                S='t20p',
                G=color,
                *_map.jxyr)

            if station_label_mapping:
                label = station_label_mapping[i_station]
            else:
                label = "%s"%s.station
            _map.add_label(s.lat, s.lon, label)
    if events:
        shifts = [(-0.6, -0.6), (-0.6, 0.)]
        for i_e, e in enumerate(events):
            if e.moment_tensor:
                _map.gmt.psxy(
                    in_columns=([e.lon], [e.lat]),
                    S='a25p',
                    G='red',
                    *_map.jxyr)

                break

                size_cm = math.sqrt(math.sqrt(mag2mom(e.moment_tensor.magnitude) / 10e7)) * size_factor
                m = e.moment_tensor
                # mc = uniso(m)
                # mc = mc / e.moment_tensor.scalar_moment() * mag2mom(5.0)
                # m6 = to6(c)
                m6 = m.m6_up_south_east()
                # data = (e.lon), (e.lat) + (10,) + m6 + (1, 0, 0)
                # data = lonlat_to_en1_km(e.lon, e.lat) + (10,) + m6 + (1, 0, 0)
                # data = [e.lon, e.lat, 5, m6, 1,0,0]
                eshift, nshift = shifts[i_e]
                # data = (e.lon, e.lat, 10, m.strike1, m.dip1, m.rake1, 1, e.lon+eshift, e.lat+nshift, 'Test site')
                data = (e.lon, e.lat, 10, 1,1,1,0,0,0, 1, e.lon, e.lat, 'Test site')
                _map.gmt.psmeca(
                    S='%s%g' % ('m', size_cm*2.0),
                    # G = gmtpy.color(colors[e.cluster]),
                    # G = colors[i_e],
                    G='red',
                    C='3p,0/0/0',
                    # W = 'thinnest,%i/%i/%i' % (255, 255, 255),
                    # L = 'thinnest,%i/%i/%i' % (255, 255, 255),
                    in_rows=[data],
                    *_map.jxyr)
    _map.save(outpath=outfn)


if __name__=='__main__':
    e = list(Event.load_catalog(filename='event.pf'))[0]
    stations = model.load_stations('array_center.pf')
    color_wet = [200, 200, 200]
    color_dry = [253, 253, 253]
    params = MapParameters(lat=e.lat, lon=e.lon, radius=8000000, outfn='array-map-new.pdf',stations=stations, events=[e],
                           show_topo=False,
                           show_grid=False,
                           color_wet=color_wet,
                           color_dry=color_dry)
    make_map(map_parameters=params)