from pyrocko import automap
import math
from pyrocko import model, gmtpy
from pyrocko.moment_tensor import magnitude_to_moment as mag2mom
from pyrocko.moment_tensor import to6
from pyrocko.guts import Object, Float, String, List


class ArrayMap(automap.Map):
    stations = List.T(model.Station.T(), optional=True)
    station_label_mapping = List.T(String.T(), optional=True)
    event = model.Event.T(optional=True)

    def __init__(self, stations=None, station_label_mapping=None, event=None,
                 *args, **kwargs):
        automap.Map.__init__(self)
        outjfn = String.T(default='map.pdf')
        self.stations = stations
        self.station_label_mapping = station_label_mapping
        self.event = event

    def save(self, outpath, resolution=75., oversample=2., size=None,
             width=None, height=None):
        gmt = self.gmt
        self.draw_labels()
        self.draw_axes()
        self.draw_cities()
        size_factor = 9.
        if self.show_topo and self.show_topo_scale:
            self._draw_topo_scale()

        if self.stations:
            lats = [s.lat for s in stations]
            lons = [s.lon for s in stations]

            self.mt.psxy(
                in_columns=(lons, lats),
                S='t10p',
                G='black',
                *self.jxyr)

            for i_station, s in enumerate(stations):
                if self.station_label_mapping:
                    label = self.station_label_mapping[i_station]
                else:
                    label = "%s"%s.station
                self.add_label(s.lat, s.lon, label)
        if self.event:
            e = self.event
            if e.moment_tensor:
                size_cm = math.sqrt(math.sqrt(mag2mom(e.magnitude) / 100e17)) * size_factor
                m = e.moment_tensor
                m6 = m.m6_up_south_east()
                data = (e.lon, e.lat, 10, m.strike1, m.dip1, m.rake1, 1, 0, 0, 'M %1.1f' % e.magnitude)
                if True:
                    self.gmt.psmeca(
                        S='%s%g' % ('a', size_cm*2.0),
                        G='red',
                        #W='thinnest,%i/%i/%i' % darken(gmtpy.color_tup(colors[e.cluster])),
                        #L='thinnest,%i/%i/%i' % darken(gmtpy.color_tup(colors[e.cluster])),
                        in_rows=[data],
                        *self.jxyr)
        gmt = self._gmt
        if outpath.endswith('.eps'):
            tmppath = gmt.tempfilename() + '.eps'
        elif outpath.endswith('.ps'):
            tmppath = gmt.tempfilename() + '.ps'
        else:
            tmppath = gmt.tempfilename() + '.pdf'

        gmt.save(tmppath)

        if any(outpath.endswith(x) for x in ('.eps', '.ps', '.pdf')):
            shutil.copy(tmppath, outpath)
        else:
            convert_graph(tmppath, outpath,
                          resolution=resolution, oversample=oversample,
                          size=size, width=width, height=height)


