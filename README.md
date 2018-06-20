# High Proper Motion Calculator

Two Python scripts which can be used to calculate the proper motion of a potential high proper motion candidate. This is a script which was presumed lost and originally written for a dissertation project. The script was first written in 2015 in Python 2 (when I didn't know much about programming!) and I decided to revive/re-write it in Python 3.

## Usage

Example data can be found in the directory `Sample Data`. To use the script, invoke Python in the usual way

```bash
$ python HPM_Py3.py
```

You will be asked to fill in the path to data to be used. By default, this assumes that the input data is somewhere relative to the current working directory, but this is easily changeable in the script by setting the variable `work_dir=False`.

To use the sample data, input `UWISH2.csv` when prompted for UWISH2 data, `UKIDSS.csv` for the UKIDSS GPS data. The coordinates for the potential HPM object are found in the file `HPM_coords.txt`, as well as the date the UWISH2 image was taken. Note that in the current Python 3 version, the HPM candidate will not look significant compared to the data sample. This will be updated in future versions when script handles the noise in the images better.

The script will output all matched objects and their proper motions to file as well as plotting all the proper motions. The potential HPM object is indicated as a red diamond, and the rings on the plot correspond to standard deviation levels around the mean of the sample.

## Differences between the Python 2 and 3 versions

`HPM_Py2.py` should work as intended for the project. However, not all of the features in this version have been implemented into the Python 3 version. Currently, the colour limiting is not working as intended. I also elected to output the J, H and K magnitudes rather than the colour excesses in the Python 3 version. The Python 2 version also removes outliers when calculating the standard deviation and outputs the potential HPM object plus any sigma-3 detections, which the Python 3 version currently does not do.
