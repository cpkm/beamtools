# -*- coding: utf-8 -*-
"""
Created on Sun Apr 2 15:30 2017

@author: cpkmanchee

Creates a dictionary of file formats and exports the dictionary using pickle.
The .pkl file can be imported and used by other modules.

"""
'''
import csv
import pickle

import os

path = os.path.dirname((os.path.abspath(__file__)))
'''

file_formats = {'bt_regen_monitor': 
                {'alias': ['regen_monitor','regen','monitor','mon'],
                 'header_lines': 9,
                 'number_data_columns': 6,
                 'column_labels': ['time','current','power','crossover','t2','t2'],
                 'column_units': ['', 'A','W','ratio','degC','degC'],
                 'delimiter': '\t',
                 },
            'thorlabs_pm': 
                {'alias': ['thorlabs_pm','thor','thorlabs','pm100','pm'],
                 'header_lines': 3,
                 'number_data_columns': 3,
                 'column_labels': ['time','power','units'],
                 'column_units': ['', 'W', ''],
                 'delimiter': '\t',
                 },
            'oceanoptics_spectrometer':
                {'alias': ['oceanoptics_spectrometer','ocean_optics_spectrometer','oo_spectrometer','oospec','oo_spec','oo'],
                 'header_lines': 0,
                 'number_data_columns': 2,
                 'column_labels': ['wavelength','intensity'],
                 'column_units': ['nm', 'units'],
                 'delimiter': ',',
                 },
            'bt_autocorrelator': 
                {'alias':['bt_autocorrelator','autocorrelator','ac','auto_correlator','auto','bt_ac','btac'],
                 'header_lines': 0,
                 'number_data_columns': 2,
                 'column_labels': ['position','power'],
                 'column_units': ['mm', 'W'],
                 'delimiter': '\t',
                 },
            'bt_beamprofiler':
                {'alias':['bt_beamprofiler','bt_beam_profiler','bt_bp','beamprofiler','bt_beampointing','beampointing'],
                 'header_lines': 2,
                 'number_data_columns': 2,
                 'column_labels': ['time','x0','y0'],
                 'column_units': ['','um','um'],
                 'delimiter': '\t',
                 },


            }

'''
with open(os.path.join(path,'file_formats.pkl'), 'wb') as f:
    pickle.dump(file_formats,f,protocol=pickle.HIGHEST_PROTOCOL)
'''