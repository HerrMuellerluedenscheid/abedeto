import numpy as num 
from pyrocko import automap
import math
from pyrocko import model, gmtpy
from pyrocko.guts import Object, Float, List, String
from pyrocko.moment_tensor import magnitude_to_moment as mag2mom
from pyrocko.moment_tensor import to6

def uniso(m):
    xx = num.trace(m) / 3.
    mc = num.matrix([[xx, 0., 0.],[0., xx, 0.], [0., 0., xx]])
    return m - mc

#width = 20.
#height = 18.
#lat = 45.
#lon = 8
#radius = 1000000
#outfn = 'map.pdf'
#stations_fn = 'GERES/array_center.pf'
#station_label_mapping = ['GERES']
#events_fn = 'castor_event_IGN.dat'
#size_factor = 9.


class MapParameters(Object):
    width = Float.T(help='width of map', optional=True, default=20)
    height = Float.T(help='heigth of map', optional=True, default=16)
    lat = Float.T(help='center latitude')
    lon = Float.T(help='center longitude')
    radius = Float.T(help='radius of map')
    outfn = Float.T(help='output filename', optional=True, default='map.png')
    stations = List.T(model.Station.T(), optional=True)
    stations_label_mapping = List.T(String.T(), optional=True, default=False)
    events = List.T(model.Event.T(), optional=True)


def make_map(lat=None, lon=None, radius=None, outfn=None, stations=None, events=None, stations_label_mapping=None, map_parameters=None):
    if map_parameters:
        width =map_parameters.width
        height=map_parameters.height
        lat=map_parameters.lat
        lon=map_parameters.lon
        radius=map_parameters.radius

    map = automap.Map( width=width, height=height, lat=lat, lon=lon,
        radius=radius,
        topo_resolution_max=200,
        topo_resolution_min=40.,
        show_topo=True,
        show_grid=True,
        illuminate=True,
        illuminate_factor_land=0.5,
        illuminate_factor_ocean=0.25)

    map.draw_cities()

    if stations:
        lats = [s.lat for s in stations]
        lons = [s.lon for s in stations]

        map.gmt.psxy(
            in_columns=(lons, lats),
            S='t10p',
            G='black',
            *map.jxyr)

        for i_station, s in enumerate(stations):
            if station_label_mapping:
                label = station_label_mapping[i_station]
            else:
                label = "%s"%s.station
            map.add_label(s.lat, s.lon, label)

    if events:
        for e in events:
            if e.moment_tensor:
                size_cm = math.sqrt(math.sqrt(mag2mom(e.magnitude) / 100e17)) * size_factor
                m = e.moment_tensor
                #mc = uniso(m)
                #mc = mc / e.moment_tensor.scalar_moment() * mag2mom(5.0)
                #m6 = to6(c)
                m6 = m.m6_up_south_east()
                #data = (e.lon), (e.lat) + (10,) + m6 + (1, 0, 0)
                #data = lonlat_to_en1_km(e.lon, e.lat) + (10,) + m6 + (1, 0, 0)
                #data = [e.lon, e.lat, 5, m6, 1,0,0]
                data = (e.lon, e.lat, 10, m.strike1, m.dip1, m.rake1, 1, 0, 0, 'M %1.1f' % e.magnitude)
                if True:
                    map.gmt.psmeca(
                        S='%s%g' % ('a', size_cm*2.0),
                        #G=gmtpy.color(colors[e.cluster]),
                        G='red',
                        #W='thinnest,%i/%i/%i' % darken(gmtpy.color_tup(colors[e.cluster])),
                        #L='thinnest,%i/%i/%i' % darken(gmtpy.color_tup(colors[e.cluster])),
                        in_rows=[data],
                        *map.jxyr)

    map.save(outfn)

if __name__=='__main__':
    params = MapParameters(lat=10, lon=10, radius=10, outfn='test-map.pdf')
    make_map(map_parameters=params)
