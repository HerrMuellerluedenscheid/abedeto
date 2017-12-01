from pyrocko import model, orthodrome as ortho
stations = model.load_stations('GERES/array_center.pf')
events = model.load_events('castor_event_IGN.dat')
for s in stations:
    print(s)
    print('distance in deg ', ortho.distance_accurate50m(s, events[0])/110.54/1000.)
