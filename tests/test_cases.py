import copy
from collections import OrderedDict
from subprocess import run

import pytest

import numpy as np
import matplotlib.pyplot as plt
import hissw


def run_idl(ebtel_idl_path, config):
    # Flags and keywords
    flags = []
    if 'dem' not in config or not config['dem']['use_new_method']:
        flags += ['dem_old']
    if not config['use_flux_limiting']:
        flags += ['classical']
    keys = [
        'flux_nt',
        'energy_nt',
        'a_tr',
        'a_c',
        'a_0',
        'l_fact',
    ]
    keywords = [f'{k}={config[k]}' for k in keys if k in config]
        
    # Heating
    time = np.arange(0, config['total_time']-config['tau'], config['tau'])
    heat = np.ones(time.shape) * config['heating']['background']
    for _e in config['heating']['events']:
        e = _e['event']
        # Rise
        i = np.where(np.logical_and(time >= e['rise_start'], time < e['rise_end']))
        heat[i] += e['magnitude'] * (time[i] - e['rise_start']) / (e['rise_end'] - e['rise_start'])
        # Plateau
        i = np.where(np.logical_and(time >= e['rise_end'], time < e['decay_start']))
        heat[i] += e['magnitude']
        # Decay
        i = np.where(np.logical_and(time >= e['decay_start'], time <= e['decay_end']))
        heat[i] += e['magnitude'] * (e['decay_end'] - time[i])/(e['decay_end'] - e['decay_start'])

    args = {
        'time': time.tolist(),
        'loop_length': config['loop_length'],
        'heat': heat.tolist(),
        'flags': flags,
        'keywords': keywords,
        'return_vars': config['return_vars'],
    }
    idl = hissw.Environment(extra_paths=[ebtel_idl_path])
    script = """time = {{ time }}
heat = {{ heat }}
loop_length = {{ loop_length }}
ebtel2,time,heat,loop_length,{{ return_vars | join(',') }}{% if flags %}, /{{ flags | join(', /') }}{% endif %}{% if keywords %}, {{ keywords | join(',') }}{% endif %}
    """
    
    return idl.run(script, args=args, save_vars=config['return_vars'])


def plot_results(result):
    fig = plt.figure(figsize=(10,10))
    for i,v in enumerate(['temperature', 'density', 'pressure']):
        ax = fig.add_subplot(2,2,i+1)
        ax.plot(result['time'], result[v])
        ax.set_xlim(result['time'][[0,-1]])
    ax = fig.add_subplot(224)
    ax.plot(result['temperature'], result['density'],)
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.show()


# Set up configuration dictionaries for all the cases we want to consider
# When adding a new test case, put it here and then add it to the list of 
# test cases below
base_config = {
    'total_time': 5e3,
    'tau': 1.0,
    'loop_length': 4e9,
    'use_flux_limiting': False,
    'heating': OrderedDict({
        'background': 3e-5,
        'events': []
    }),
    'return_vars': [
        'temperature',
        'density',
        'pressure',
        'velocity',
    ]
}
# Case 1: No events
no_events = copy.deepcopy(base_config)
# Case 2: 1 event
one_event = copy.deepcopy(base_config)
one_event['heating']['events'] = [
    {'event': {
        'rise_start':0,
        'rise_end':100,
        'decay_start':100,
        'decay_end':200,
        'magnitude':0.1,
    }}
]
# Case 3: Non-uniform cross-sections
varying_xs = copy.deepcopy(one_event)
varying_xs['a_c'] = 3.0
varying_xs['a_0'] = 2.0
varying_xs['a_tr'] = 1.0
varying_xs['l_fact'] = 0.85
# Case 4: 
mulitple_returns = copy.deepcopy(base_config)
mulitple_returns['return_vars'] = mulitple_returns['return_vars'] + [
    'ta', 
    'na', 
    'pa', 
    'c11', 
    'dem_tr',
    'dem_cor',
    'logtdem', 
    'f_ratio', 
    'rad_ratio', 
    'cond', 
    'rad_cor',
]


@pytest.mark.parametrize('config', [
    no_events,
    one_event,
    varying_xs,
    mulitple_returns,
])
def test_all_cases(ebtel_idl_path, show_plots, config):
    r = run_idl(ebtel_idl_path, config)
    for k in config['return_vars']:
        assert k in r
    if show_plots:
        plot_results(r)
