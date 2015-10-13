import numpy as num 
from pyrocko import automap
import math
from pyrocko import model, gmtpy
from pyrocko.moment_tensor import magnitude_to_moment as mag2mom
from pyrocko.moment_tensor import to6

def uniso(m):
    xx = num.trace(m) / 3.
    mc = num.matrix([[xx, 0., 0.],[0., xx, 0.], [0., 0., xx]])
    return m - mc

width = 20.
height = 18.
lat = 45.
lon = 8
radius = 1000000
outfn = 'map.pdf'
stations_fn = 'GERES/array_center.pf'
station_label_mapping = ['GERES']
events_fn = 'castor_event_IGN.dat'
size_factor = 9.




map = automap.Map(
    width=width,
    height=height,
    lat=lat,
    lon=lon,
    radius=radius,
    topo_resolution_max=200,
    topo_resolution_min=40.,
    show_topo=True,
    show_grid=True,
    illuminate=True,
    illuminate_factor_land=0.5,
    illuminate_factor_ocean=0.25)

map.draw_cities()

if stations_fn:
    stations = model.load_stations(stations_fn)
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

if events_fn:
    events = model.load_events(events_fn)
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
            print m6
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



    #lats = [e.lat for e in events]
    #lons = [e.lon for e in events]

    #depths = [e.depth for e in events]
    #strikes = [e.moment_tensor.strike1 for e in events]
    #dips = [e.moment_tensor.dip1 for e in events]
    #rakes = [e.moment_tensor.rake1 for e in events]
    #mags = [e.moment_tensor.magnitude for e in events]
    #newx = [0] * len(events)
    #newy = [0] * len(events)
    #titles = [''] * len(events)
    #print titles
    #print newx 
    #print newy 
    #print mags
    #print rakes 

    #map.gmt.psmeca(
    #    in_columns=(lons, lats, depths, strikes, dips, rakes, mags, newx, newy, titles),
    #    #S='a',
    #    #G='red',
    #    *map.jxyr)
    #map.gmt.psxy(
    #    in_columns=(lons, lats),
    #    S='a20p',
    #    G='red',
    #    *map.jxyr)

map.save(outfn)

