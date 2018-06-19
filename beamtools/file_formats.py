# -*- coding: utf-8 -*-
"""
Created on Sun Apr 2 15:30 2017

@author: cpkmanchee

Dictionary of file formats
"""

file_formats = {
            'bt_regen_monitor': 
                {'alias': ['bt_regen_monitor','regen_monitor','regen'],
                 'header_lines': 9,
                 'number_data_columns': 6,
                 'column_labels': ['time','current','power','crossover','t1','t2'],
                 'column_units': ['', 'A','W','ratio','degC','degC'],
                 'delimiter': '\t',
                 'time_format': '%Y/%m/%d %H:%M:%S.%f',
                 },
            'bt_rod_monitor': 
                {'alias': ['bt_rod_monitor','rod_monitor','rod'],
                 'header_lines': 9,
                 'number_data_columns': 5,
                 'column_labels': ['time','current','power','t1','t2'],
                 'column_units': ['', 'A','W','degC','degC'],
                 'delimiter': '\t',
                 'time_format': '%Y/%m/%d %H:%M:%S.%f',
                 },
            'bt_autocorrelator': 
                {'alias':['bt_autocorrelator','autocorrelator','ac','auto_correlator','auto','bt_ac','btac'],
                 'header_lines': 0,
                 'number_data_columns': 3,
                 'column_labels': ['position', 'delay', 'power'],
                 'column_units': ['um', 'ps', 'W'],
                 'delimiter': '\t',
                 },
            'bt_autocorrelator_old': 
                {'alias':['bt_autocorrelator_old','autocorrelator_old','ac_old','auto_correlator_old','auto_old','bt_ac_old','btac_old'],
                 'header_lines': 0,
                 'number_data_columns': 2,
                 'column_labels': ['position','power'],
                 'column_units': ['um', 'W'],
                 'delimiter': '\t',
                 },
            'bt_beamprofiler':
                {'alias':['bt_beamprofiler','bt_beam_profiler','bt_bp','beamprofiler','bt_beampointing','beampointing'],
                 'header_lines': 2,
                 'number_data_columns': 3,
                 'column_labels': ['time','x0','y0'],
                 'column_units': ['','um','um'],
                 'delimiter': '\t',
                 },
            'bt_knifeedge':
                {'alias':['bt_knifeedge','knifeedge','bt_ke'],
                 'header_lines': 1,
                 'number_data_columns': 2,
                 'column_labels': ['position','power'],
                 'column_units': ['um','W'],
                 'delimiter': ',',
                 },
            'thorlabs_pm': 
                {'alias': ['thorlabs_pm','thor','thorlabs','pm100','pm'],
                 'header_lines': 2,
                 'number_data_columns': 3,
                 'column_labels': ['time','power','units'],
                 'column_units': ['', 'W', ''],
                 'delimiter': '\t',
                 'time_format': '%m/%d/%Y %I:%M:%S.%f %p',
                 },
            'oceanoptics_spectrometer':
                {'alias': ['oceanoptics_spectrometer','ocean_optics_spectrometer','oo_spectrometer','oospec','oo_spec','oo'],
                 'header_lines': 0,
                 'number_data_columns': 2,
                 'column_labels': ['wavelength','intensity'],
                 'column_units': ['nm', 'units'],
                 'delimiter': ',',
                 },
            'picoscope': 
                {'alias': ['picoscope','pico','pico_scope'],
                 'header_lines': 3,
                 'number_data_columns': 2,
                 'column_labels': ['time','volt'],
                 'column_units': ['s', 'mV'],
                 'delimiter': ',',
                 },
            'picoscope_2ch': 
                {'alias': ['picoscope_2ch','pico_2ch','pico_scope_2ch','pico2'],
                 'header_lines': 3,
                 'number_data_columns': 3,
                 'column_labels': ['time','voltA','voltB'],
                 'column_units': ['s', 'mV','mV'],
                 'delimiter': ',',
                 },
            'picoscope_spectrum': 
                {'alias': ['picoscope_spectrum','pico_spec','pico_scope_spectrum'],
                 'header_lines': 3,
                 'number_data_columns': 2,
                 'column_labels': ['freq','volt'],
                 'column_units': ['MHz', 'mV'],
                 'delimiter': ',',
                 },
            'bt_crosssection': 
                {'alias': ['bt_crosssection','bt_cs','crosssection'],
                 'header_lines': 7,
                 'number_data_columns': 3,
                 'column_labels': ['wavelength','absorption','emission'],
                 'column_units': ['m','m**2','m**2'],
                 'delimiter': ',',
                 }
            }
