# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 09:32:19 2016

@author: cpkmanchee
"""

import numpy as np
import csv
import pickle


file_formats = {'regen_monitor': 
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
                 'column_units': ['mm', 'W', ''],
                 'delimiter': '\t',
                 },
            'ocean_optics_spectrometer':
                {'alias': ['ocean_optics_spectrometer','oo_spectrometer','oospec','oo_spec','oo'],
                 'header_lines': 0,
                 'number_data_columns': 2,
                 'column_labels': ['wavelength','intensity'],
                 'column_units': ['nm', 'units'],
                 'delimiter': ',',
                 },
            'autocorrelator': 
                {'alias':['autocorrelator','ac','auto_correlator','auto'],
                 'header_lines': 0,
                 'number_data_columns': 2,
                 'column_labels': ['position','power'],
                 'column_units': ['mm', 'W'],
                 'delimiter': '\t',
                 }

            }



#w = csv.writer(open("file_formats.csv", "w"))
#for key, val in file_formats.items():
#    w.writerow([key, val]) 
with open('file_formats.pkl', 'wb') as f:
    pickle.dump(file_formats,f)
