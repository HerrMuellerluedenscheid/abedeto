import os
import numpy as num
import math
import bisect
from pyrocko import crust2x2, orthodrome as ortho
from pyrocko.fomosto import qseis
from pyrocko.gf.store import Store, remake_dir
from pyrocko.gf.meta import Timing, ConfigTypeA, TPDef
from pyrocko import cake
import logging


km = 1000.0
pjoin = os.path.join
logger = logging.getLogger("propose-store")


class NoRay(Exception):
    def __init__(self, context):
        Exception.__init__(self, "NoRay at " + context)


def model_has_cmb(mod):
    disks = mod.discontinuities()
    return "cmb" in [d.name for d in disks]


def model_has_icb(mod):
    disks = mod.discontinuities()
    return "icb" in [d.name for d in disks]


def adjust_earthmodel_receiver_depth(config):
    rmod = config.earthmodel_receiver_1d
    smod = config.earthmodel_1d
    last_layer = list(rmod.elements())[-1]
    ll_zbot = last_layer.zbot
    s_z = smod.profile("z")
    if ll_zbot in s_z:
        return
    else:
        last_layer.zbot = s_z[bisect.bisect(s_z, ll_zbot)]


def setup_distances_crusts(stations, events, do_round=True):
    """Propose a fomosto store configuration for P-pP Array beam forming.
    :param event: Event instance
    :param station: Station instance."""
    distances = {}
    models = {}
    for s in stations:
        for e in events:
            event_profile = crust2x2.get_profile(e.lat, e.lon)
            station_profile = crust2x2.get_profile(s.lat, s.lon)
            k1 = station_profile._ident
            k2 = event_profile._ident
            models[k1] = station_profile
            models[k2] = event_profile
            key = (s.station, k1, k2)
            distance = ortho.distance_accurate50m(e, s)
            if not key in distances:
                if do_round:
                    distances[key] = (math.floor(distance - 1), math.ceil(distance + 1))
                else:
                    distances[key] = (distance - 1, distance + 1)
            else:
                d1, d2 = distances[key]
                if do_round:
                    distances[key] = (
                        min(d1, math.floor(distance - 1)),
                        max(d2, math.ceil(distance + 1)),
                    )
                else:
                    distances[key] = (min(d1, distance - 1), max(d2, distance + 1))

    return distances, models


# def propose_store(station, events, superdir, source_depth_min=0.,
#                  source_depth_max=15., source_depth_delta=1., sample_rate=10.,
#                  force=False, numdists=2, run_ttt=False, simplify=False,
#                  phases=['P'], classic=True):
def propose_stores(
    distances,
    models,
    superdir,
    source_depth_min=0.0,
    source_depth_max=15.0,
    source_depth_delta=1.0,
    sample_rate=10.0,
    force=False,
    numdists=2,
    run_ttt=False,
    simplify=False,
    phases=["P"],
    classic=True,
    distance_delta_max=None,
):
    """Propose a fomosto store configuration for P-pP Array beam forming.
    :param event: Event instance
    :param superdir: where to create the store (default, current directory)
    :param source_depth_min: minimum source depth (default 0)
    :param source_depth_max: maximum source deoth (default 15)
    :param source_depth_delta: increment
    :param sample_rate: in Hz
    :param force_overwrite: overwrite potentially existent store
    :param run_ttt: generate travel time tables right away"""

    modelling_code_id = "qseis.2006b"

    configs = []
    if classic:
        define_method = cake.PhaseDef
    else:
        define_method = cake.PhaseDef.classic

    wanted = [define_method(ph) for ph in phases]

    global_model = cake.load_model()
    remake_dir(superdir, force)
    for (station_id, key_station, key_event), (dist_min, dist_max) in distances.items():

        configid = "%s_%s_%s" % (station_id, key_station, key_event)
        distance_delta = dist_max - dist_min
        if distance_delta_max is not None:
            while distance_delta > distance_delta_max:
                distance_delta /= 2.0
        config = ConfigTypeA(
            id=configid,
            source_depth_min=source_depth_min * km,
            source_depth_max=source_depth_max * km,
            source_depth_delta=source_depth_delta * km,
            distance_min=dist_min,
            distance_max=dist_max,
            distance_delta=distance_delta,
            sample_rate=sample_rate,
            ncomponents=10,
        )

        station_crust = models[key_station]
        config.earthmodel_receiver_1d = cake.LayeredModel.from_scanlines(
            cake.from_crust2x2_profile(station_crust)
        )

        config.earthmodel_1d = global_model.replaced_crust(
            crust2_profile=models[key_event]
        )

        if simplify:
            config.earthmodel_1d = config.earthmodel_1d.simplify(max_rel_error=0.002)
        adjust_earthmodel_receiver_depth(config)
        configs.append(config)
        dest_dir = pjoin(superdir, config.id)
        remake_dir(dest_dir, force)
        logger.info("Created store: %s" % dest_dir)

        mean_z = num.mean((config.source_depth_min, config.source_depth_max))
        mean_dist = num.mean((config.distance_min, config.distance_max))
        arrivals = config.earthmodel_1d.arrivals(
            phases=wanted, distances=[mean_dist * cake.m2d], zstart=mean_z
        )
        if len(arrivals) == 0:
            logger.warning(
                NoRay(
                    "d: %s, z: %s, %s phases: %s"
                    % (
                        mean_dist * cake.m2d,
                        mean_z,
                        "classic" if classic else "",
                        "|".join(phases),
                    )
                )
            )
            slow = 0.1
            slowness_taper = (0.0, 0.0, 1.3 * slow, 1.5 * slow)
            z_turn = num.max(config.earthmodel_1d.profile("z"))
        else:
            slow = arrivals[0].p / (cake.r2d * cake.d2m / km)
            slowness_taper = (0.3 * slow, 0.5 * slow, 1.5 * slow, 1.7 * slow)
            z_turn = num.max(arrivals[0].zxt_path_subdivided()[0])

        zmax = max(z_turn * 1.1, config.earthmodel_receiver_1d.profile("z")[-1])

        config.earthmodel_1d = config.earthmodel_1d.extract(depth_max=zmax)
        begin_phase_defs = "P,P\\,PP"
        if model_has_icb(config.earthmodel_1d):
            begin_phase_defs += ",P(cmb)P(icb)P(icb)p(cmb)p,P(cmb)P<(icb)(cmb)p"
        elif model_has_cmb(config.earthmodel_1d):
            begin_phase_defs += ",Pv_(cmb)p"
        config.modelling_code_id = modelling_code_id
        config.tabulated_phases = [
            TPDef(id="begin", definition=begin_phase_defs),
            TPDef(id="end", definition="2.5"),
            TPDef(id="PP", definition="PP"),
            TPDef(id="P", definition="P"),
        ]

        qs = qseis.QSeisConfig()
        qs.qseis_version = config.modelling_code_id.split(".")[1]
        half_lapse_time = 55
        qs.time_region = (
            Timing("begin-%s" % (half_lapse_time * 1.1)),
            Timing("begin+%s" % (half_lapse_time * 1.1)),
        )
        qs.cut = (
            Timing("begin-%s" % half_lapse_time),
            Timing("begin+%s" % half_lapse_time),
        )
        qs.slowness_window = slowness_taper
        qs.wavelet_duration_samples = 0.001
        qs.sw_flat_earth_transform = 1
        qs.filter_shallow_paths = 1
        qs.filter_shallow_paths_depth = float(z_turn * 0.2)
        qs.sw_algorithm = 1
        Store.create_editables(dest_dir, config=config, extra={"qseis": qs})
        if run_ttt:
            st = Store(dest_dir)
            st.make_ttt()

    config_ids = [c.id for c in configs]
    return config_ids


if __name__ == "__main__":
    import argparse
    from pyrocko.model import Event, load_stations

    parser = argparse.ArgumentParser("suggest a store for P phases only")
    parser.add_argument("--stations", help="stations file")
    parser.add_argument("--events", help="event file")
    parser.add_argument("--force", action="store_true", help="force_overwrite")
    parser.add_argument(
        "--superdir", default=".", help="directory where to put the store"
    )
    parser.add_argument(
        "--number_of_distances",
        help="number of distances between outer grid nodes in GFDB",
        default=2,
    )

    args = parser.parse_args()

    stations = load_stations(args.stations)
    if len(stations) == 1:
        s = stations[0]

    events = list(Event.load_catalog(args.events))

    propose_store(s, events, superdir=args.superdir, force=args.force)
