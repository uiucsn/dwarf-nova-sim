

# Extract Outburst, model outbursts, and select good candidatdes

## Related Scripts

- `cand_all.py`

- `light_curve_fitting.py`

>> - `Light_curve_bazin.py`

## Preparation

Install packages:
```sh
python3 -mpip install pandas scipy matplotlib numpy light_curve
```

Data: OGLE dwarf novae data in `./phot`

## Results

outburst extracted

outburst extracted plots

fitting plots, which is by default saved under `./pictures` directory.

## Run Code

### Extract outburst

`cand_all.py` extracts candidate outburst and return the time arrays and magnitudes arrays. It also plots the extracted data against the original OGLE data.
>> Doesn't produce all figures

After manually check each selected outburst, couple of good candidates have been selected and their information is stored in a csv file. 

### Fitting model

`light_curve_fitting.py` fits the extracted outbursts with piecewise function and uses `scipy.optimize` to fit the data. It also plots the model against data per light curve and is by default saved under `./pictures` directory.
>> fails with `FileNotFoundError: [Errno 2] No such file or directory: 'Dwarf nova.csv'`

>> `Light_curve_bazin.py` used Bazin function to fit the data and plot it. 

Further steps use fitting model from `light_curve` package

## Provided results

Selected objects: 

> selected_obj = ['OGLE BLG-DN-0001', 'OGLE BLG-DN-0002', 'OGLE BLG-DN-0036', 'OGLE BLG-DN-0087',  'OGLE BLG-DN-0168', 'OGLE BLG-DN-0174', 'OGLE BLG-DN-0233', 'OGLE BLG-DN-0275',  'OGLE BLG-DN-0286', 'OGLE BLG-DN-0305', 'OGLE BLG-DN-0373', 'OGLE BLG-DN-0376', 'OGLE BLG-DN-0421', 'OGLE BLG-DN-0444', 'OGLE BLG-DN-0458', 'OGLE BLG-DN-0531', 'OGLE BLG-DN-0588', 'OGLE BLG-DN-0595', 'OGLE BLG-DN-0690', 'OGLE BLG-DN-0783', 'OGLE BLG-DN-0826', 'OGLE BLG-DN-0899']

Model Fitting information: saved as a csv files `analysis/fitting_info.csv`, with OGLE ID, outburst index, starting time, ending time, and model parameters. 



# Retrieve Distance and extinction

## Related scripts

>> - `Get_distance.py`

- `get_extinction.py`

>> ## Preperation
>> 
>> Packages: cand_all, numpy, matplotlib, os, light_curve
>> 
>> Data: phot
>>
>> **Remove?**

>> Describe dustmaps data file fetching
```sh
python3 -c 'from dustmaps import sfd, bayestar; sfd.fetch(); beystar.fetch()'
```

## Run Code

Searching scope for Gaia is 1 arcsec. If no objects is found in Gaia, then there will be no distane data for this OGLE object.

If distance is not provided from Gaia, extinction is given by SFD

If distance is provided from Gaia, extinction is given by Bayestar

## Results:

Get_distance.py provides distance information for all extracted objects (not the manually selected objects) and is by default saved as "analysis/distance_1arcsec.csv". Objects without data from Gaia will be skipped.

`get_extinction.py` provides extinction information for all selected objects. Columns: OGLE ID (name), galactic coordinates (l,b), extinction from Bayestar, SFD, and Marshall.

## Provided Results:

Objects name, coordinates(l, b), distances, extinctions for all selected objects are given in analysis/extinction_1arces.csv

# Find Accretion rate and luminosities in other passbands

## Related Codes

Find_Mdot_DerivedLuminosity.py

## Preparation:

Packages: spicy, satrapy, functors, numpy, (matplotlib), os, light_curve, cand_all

Data:phot, analysis/fitting_info.csv, analysis/extinction_1arcese.csv

## Run code 

The code will extract outburst from original OGLE data, use Basin function to model the extracted outbursts, use standard disk model to derive accretion rate M_dot, and derive luminosity at LSST passbands. 

## Results

--Two sets of results are provided: 

In analysis_Mdot directory, there are individual csv files for each outburst with columns of time (t), accretion rate (Mdot), and luminosites in each LSST passband (L_u, L_g, etc.). The name for each file has a format of OGLE ID + outburst index

In analysis_luminosities directory, there are individual csv fiels fore ach outburst with columns of lumimosities and magnitudes. The name for each file has a format of OGLE ID + outburst index + luminosity. 

--When running the code, basic information will be printed: the method used to find Mdot (method 2: no distance provided and use 1kpc), whether the mode involves normalization, Mdot before normalization(Mdot_old), Mdot after normalization(Mdot_norm)

## Provided Results

The code has already been run and results have already been saved under analysis_Mdot and analysis_Luminosity, which can be used directly for generator.py

# Generating realsitic objects and write LCLIB file

## Related Code

Generator.py

## Preparation

packages: lumpy, astropy, matplotlib, pandas, dustmaps, ch_vars

Data: analysis_Mdot, dustmaps

## Run Code

Running generator.py will generate 10k realistic outbursts in all LSST passbands and output them into LCLIB file. 

It used numpy random generator with starting index of 42

Outburst instances are the selected outbursts. Coordinates instances follow milkyway density model.

It filters out instances with any luminosity number outside the range [5, 99]

## Results

The code will generate a txt file containing the header for the file and 10k outburst instances. 

## Provided results

The LCLIB txt file result is provided and is compressed as LCLIB_dwarf_nova_sim.txt.zip

