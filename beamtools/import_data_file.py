"""
Created on Fri Apr 7 21:49:16 2017

@author: cpkmanchee
"""

import numpy as np
import os
import csv
import pickle
import warnings

from beamtools.file_formats import file_formats
from beamtools.common import DataObj


__all__ = ['import_data_file', 'list_atr','list_filetypes']


def list_filetypes():
    '''Display all filetypes in dictionary
    '''
    [print(k,v) for k,v in file_formats['filetype'].items()]
    return

def list_atr(given_filetype):
    '''List the attributes of resultant object from data import.
    '''
    filetype = filetype_lookup(file_formats,given_filetype.lower())
    column_labels = file_formats.get(filetype).get('column_labels')
    print(column_labels)
    return

                
def filetype_lookup(file_dict, given_type):
    '''Identify file type for given input. Only first found match is returned.
    '''
    for k,v in file_dict.items():
        if given_type in file_dict.get(k).get('alias'):
            return(k)
        
    raise RuntimeError('File type lookup failed. File type "%s" not found' %(given_type))
    return(None)

def import_data_file(file, given_filetype):
    '''Imports data of given filetype.
    Data returned as object with appropriate attributes
    '''

    filetype = filetype_lookup(file_formats,given_filetype.lower())

    header_lines = file_formats.get(filetype).get('header_lines')
    delimiter = file_formats.get(filetype).get('delimiter')
    column_labels = file_formats.get(filetype).get('column_labels')

    #initialize header and output dictionary
    header=[]
    output={}
    [output.update({c:[]}) for c in column_labels]
    with open(file, 'r') as f:
        #extract header information only
        data = csv.reader(f, delimiter = delimiter)
        for i in range(header_lines):
            header.append(data.__next__())
        #write rest of data to dictionary, keys are column_labels, values = data 
        [[(output[c].append(row[c_ind].strip())) for c_ind,c in enumerate(column_labels)] for row in data]

    #convert data to float
    v = []
    for c in output.keys():
        try:
            output[c] = np.asarray(output[c],dtype=np.float)
            v.append(~np.isnan(output[c]))

        except ValueError:
            warnings.warn('Unable to convert to float')
            try:
                output[c] = np.asarray(output[c])
            except ValueError:
                warning.warn('Unable to cast as array')
                
    #nan correction
    try:
        v = np.asarray(v)
        if not (v.size == 0):
            valid = np.prod(v,0).astype(bool)
            for c in output.keys():
                output[c] = output[c][valid]
    except:
        warnings.warn('Nan processing failure')

    output.update({'header': header})
    output.update({'filetype': filetype})
    output_obj = DataObj(output)

    return output_obj