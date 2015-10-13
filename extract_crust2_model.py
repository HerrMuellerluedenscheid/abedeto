import sys
import numpy as num
from pyrocko import crust2x2, orthodrome as ortho
from pyrocko.fomosto import qseis
from pyrocko.gf.store import Store
from pyrocko.gf.meta import Timing, ConfigTypeA, TPDef
from pyrocko import cake
import copy
import tempfile
km = 1000.


def propose_store(station, event, source_depth_min=0., source_depth_max=15., source_depth_delta=1., sample_rate=10.):
    ''' Propose a fomosto store configuration for P-pP Array beam forming.'''
    tempdir = tempfile.mkdtemp(prefix='arraybeam')
    modelling_code_id = 'qseis.2006a'
    distance = ortho.distance_accurate50m(station, event)
    depth = event.depth
    config = ConfigTypeA(id=tempdir.split('/')[-1],
                         source_depth_min=source_depth_min*km,
                         source_depth_max=source_depth_max*km,
                         source_depth_delta=source_depth_delta*km,
                         distance_min=0.,
                         distance_max=0.,
                         distance_delta=0.,
                         sample_rate=sample_rate)

    config.distance_min = distance-5*km
    config.distance_max = distance+5*km
    config.distance_delta = 10*km
    for item in zip([station, event],['earthmodel_1d', 'earthmodel_receiver_1d']):
        location, model_id = item
        if model_id=='earthmodel_1d':
            mod = cake.load_model()
            mod = mod.replaced_crust((location.lat, location.lon))
        else:
            mod = cake.LayeredModel.from_scanlines(
                cake.from_crust2x2_profile(
                    crust2x2.get_profile(location.lat, location.lon)))
        setattr(config, model_id, mod)
    mean_z = num.mean([config.source_depth_min, config.source_depth_max])



    arrivals = config.earthmodel_1d.arrivals(phases=cake.PhaseDef('P'),
                                             distances=[distance*cake.m2d],
                                             zstart=mean_z)
    slow = arrivals[0].p/(cake.r2d*cake.d2m/km)

    config.modelling_code_id=modelling_code_id
    config.tabulated_phases=[
        TPDef(
            id='begin',
            definition='p,P,p\\,P\\,Pv_(cmb)p'),
        TPDef(
            id='end',
            definition='2.5'),
        TPDef(
            id='P',
            definition='!P'),
        TPDef(
            id='S',
            definition='!S'),
        TPDef(
            id='p',
            definition='!p'),
        TPDef(
            id='s',
            definition='!s')]

    qs = qseis.QSeisConfig()
    qs.time_region = (Timing('P-50'), Timing('P+60'))
    qs.cut = (Timing('P-50'), Timing('P+60'))
    qs.slowness_window = (0., 0., slow+slow*0.2, slow+slow*0.25)
    qs.validate()
    config.validate()
    Store.create_editables(tempdir, config=config, extra={'qseis': qs})

    return dempdir



if __name__=="__main__":
    from pyrocko.model import Station, Event
    # for testing:

    s = Station()
    e = Event(lat=90)
    propose_store(s, e)
