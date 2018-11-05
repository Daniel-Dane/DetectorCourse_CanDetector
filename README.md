Detector Course Report
============

This is the LaTeX document and Python code for the can gas detector and semiconductor detector from the Nerdic Detector Course.

## Run python scripts

Below is a hint on how to run the python scripts to produce desired output files:

### HV scan plots
The script produces the high voltage (HV) scan plots using the `.mca` files in `CanDetector/data/mca/`. The scripts are hardcoded
as they will not change. To produce the voltage scan plots and the overview of different gains and voltages one can simply run

```
./plot_HVscan.py
```

from `CanDetector/scripts`, which will produce the output files in `CanDetector/graphics`. With the additional flag

```
./plot_HVscan.py -s fe
```

only the iron plots are produced. Alternatively `am` can be used as a command line argument for the americium ones.


## Compile report
The subdirectories `CanDetector` and `SemiconductorDetector` contain the files needed to compile the corresponding reports.

### Can Detector
To compile the can detector a makefile is provided, therefore simply running:
```
make
```
in the directory is sufficient. For later iterations one needs to call
```
make cleanpdf
```
to remove the target PDF or `make cleanall` to remove all temporary LaTeX files as well.
