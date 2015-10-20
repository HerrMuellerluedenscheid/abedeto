#!/usr/bin/env python
import argparse
import numpy as num
import matplotlib.pyplot as plt
from pyrocko import model, orthodrome, util, gui_util


def get_bounds(stations, events=None, usestations=False, printall=False, show_fig=False, km=False):
    scale = 1./1000. if km else 1
    if printall:
        print 'using stations: '
        for s in stations:
            print s
            print '.......................................\n\n\n'
    maxdistances = []
    mindistances = []
    alldists = {}
    for e in events:
        dists = [orthodrome.distance_accurate50m(e, s)*scale for s in stations]
        alldists[e] = dists
        maxdistances.append(max(dists))
        mindistances.append(min(dists))

    depths = [e.depth for e in events]
    if printall:
        for e, dists in alldists.items():
            i = 0
            for s in stations:
                print '%s \n %s DISTANCE: %s \n\n' % (e, s, dists[i])
                print '='*20
                i += 1
        if len(maxdistances)==1:
            print 'maximum distances: ', maxdistances[0]
            print 'minimum distances: ', mindistances[0]

    if show_fig:
        f, axs = plt.subplots(2)

        axs[0].hist(maxdistances, bins=20, color='r')
        axs[0].hist(mindistances, bins=20, color='g')
        axs[0].set_title('Minimum (green) and maximum (red) distances [km]')
        axs[1].hist(depths, bins=20)
        axs[1].set_title('Depths [km]')
        axs[0].get_xaxis().get_major_formatter().set_useOffset(False)
        axs[1].get_xaxis().get_major_formatter().set_useOffset(False)
        plt.show()

    else:
        return mindistances, maxdistances, depths


"""maximum and minimum distances of events with respect to stations. Requires python >= 2.7 
due to usage of argparse."""
if __name__=="__main__":

    parser = argparse.ArgumentParser(description='Find minimum and maximum distances/depth.')
    parser.add_argument('--markers',
                        help='file name containing event markers',
                        required=False,
                        default=False)
    parser.add_argument('--usestations',
                        help='regular expression to filter stations',
                        required=False,
                        default=False)
    parser.add_argument('--stations',
                        help='name of file containing station information',
                        required=True)
    parser.add_argument('--events',
                        help='name of file containing event catalog',
                        default=False,
                        required=False)
    parser.add_argument('--printall',
                        help='Print all results to terminal',
                        default=True,
                        required=False,
                        action='store_true')
    parser.add_argument('--show',
                        help='show figure at the end',
                        default=False,
                        required=False,
                        action='store_true')
    args = parser.parse_args()

    stations = model.load_stations(args.stations)

    if args.usestations:
        stations = [s for s in stations if util.match_nslc(args.usestations, s.nsl())]

    events = []
    if args.events:
        events.extend(model.load_events(args.events))
    if args.markers:
        markers = gui_util.load_markers(args.markers)
        events.extend([m.get_event() for m in markers])
    get_bounds(stations, events=events, usestations=args.usestations, printall=args.printall, show_fig=args.show)
