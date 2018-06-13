# beamtools
The *beamtools* package is a collection of methods to assist with several aspects of optical system design and analysis.

## Importing data
One of the most useful methods in the package is `import_data_file` which can be used to seamlessly import a variety of data types.

### Importing data with `beamtools.import_data_file()`
Once of the most useful methods in the package is `import_data_file` which can be used to seamlessly import a variety of data types.
```
import beamtools as bt

#data from an ocean optics spectrometer
file = 'testdata.csv'
type = 'oo_spec'

data = bt.import_data_file(file,type)
```
Data imported by `import_data_file` is returned as a DataObject. In this example, `data.wavelength` and `data.intensity` returns numpy arrays of the wavelength and intensity data imported from the file.
Any header information is stored in `data.header`. The header syntax for each file is defined in `file_formats.py`.

### Data Objects
All raw data in the package is handle via a DataObj where the data is stored as object attributes. Data objects are created from dictionaries (set of key/value pairs) and the Data Object attributes are defined by the dictionary keys. The class definition can be found in `.import_data_file.py`.
All `DataObj` instances posses the class method `fields()`. Calling `DataObj.fields()` returns a list of the onjects attributes. This can be used to view the names available data arrays.

### File formats
A list of file formats can be found in `file_formats.py`. I have already included many common file types. More can be added to the package in the future by creating additional dictionary entries following the format and keywords given. Each file format includes a key named `alias`. These values comprise a set of alternate names which can be used anywhere in the package when specifying a file type.

## Common methods

## `beamtools.pulse` module
The `beamtools.pulse` module contains many methods for analysing optical pulse. This includes the analysis of spectra and autocorrelation traces.

## `beamtools.upp` module
The ultrafast pulse propagation module, imported as `beamtools.ufpp` possess numerous methods pertaining to the numerical propogation of ultrafast pulses. This includes nonlinear fiber-optic equations, fiber amplification, and several simulated fiber-optic components. 

## `beamtools.grating` module
