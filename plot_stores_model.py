import matplotlib.pyplot as plt
from pyrocko.gf import LocalEngine


def make_plot(which, store_ids, engine, parameters, out_filename,
              max_depth, xlim, store_id_mapping=None, combine=True):
    '''create pretty model plots'''
    if store_id_mapping and not len(store_id_mapping)==len(store_ids):
        raise Exception('If id mapping is used, it needs to have the same length '
                        'as store_ids')

    if combine:
        wantrows = 1
    else:
        wantrows = len(store_ids)

    fig, axs = plt.subplots(len(parameters),
                            wantrows,
                            figsize=(wantrows*5, len(parameters)*3))
    if combine:
        axs = [axs]
    if len(parameters)==1:
        axs = [axs]


    for i_store, store_id in enumerate(store_ids):
        store = engine.get_store(store_id)
        if which=='source':
            mod = store.config.earthmodel_1d
        elif which=='receiver':
            mod = store.config.earthmodel_receiver_1d
        else:
            raise Exception('which can only be receiver or source')
        if combine:
            i_row = 0
        else:
            i_row = i_store

        for i_param, param in enumerate(parameters):
            ax = axs[i_row][i_param]
            if store_id_mapping:
                label = store_id_mapping[i_store]
            else:
                label = store_id
            ax.plot(mod.profile(param)/1000., mod.profile('z')/1000., label=label)
            ax.set_ylim([0, max_depth/1000.])
            if xlim:
                ax.set_xlim(xlim)
            ax.set_ylabel('Depth [km]')
            ax.set_xlabel('$v_p$ [km/s]')
            if combine and len(store_ids)>1:
                ax.legend(fontsize=9)
            else:
                ax.set_title(label)
            ax.invert_yaxis()
    fig.savefig(out_filename,
                pad_inches=0.04,
                bbox_inches='tight')


engine = LocalEngine(use_config=True)
parameters = ['vp']
store_ids = ['castor_geres_4']
max_depth = 15000.
xlim = (2,9)


# Can be receiver or source:
store_id_mapping = ['Source Site']
out_filename = 'velocity_model_source.pdf'
which = 'source'

make_plot(which, store_ids, engine, parameters, out_filename, max_depth, xlim, store_id_mapping)

##############################################
out_filename = 'velocity_models.pdf'
store_id_mapping = ['GERES']
which = 'receiver'
make_plot(which, store_ids, engine, parameters, out_filename, max_depth, xlim, store_id_mapping)
