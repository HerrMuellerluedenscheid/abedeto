from pyrocko.fdsn import ws
from pyrocko import model

if __name__=="__main__":
    import argparse
    arrays = {'YKA': ('CN', 'YKA*', '*', 'SHZ')}

    parser = argparse.ArgumentParser(description='download waveforms')
    parser.add_argument('--events')
    args = parser.parse_args()
    e = list(model.Event.load_catalog(args.events))[0]
    tmin = ws.sdatetime(e.time)
    tmax = ws.sdatetime(e.time+300)
    for array_id, codes in arrays.items():
        selection = [codes+tuple((e.time, e.time+300))]
        d = ws.dataselect(site='iris', selection=selection)
        fn = '%s.mseed' % array_id
        with open(fn, 'w') as f:
            f.write(d.read())
            f.close()
