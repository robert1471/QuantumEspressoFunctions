import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import pandas as pd
import re
import os
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)


def pdos(atom, orbital, spin):
    """Function takes all PDOS files of a given atomic species and extracts the PDOS energy values
    returns standard df file that can be used for plotting"""
    df_master = pd.DataFrame()
    # df_master = pd.read_csv('./PDOS/MgAl2O4.pdos_tot', delimiter='\s+', header=None, usecols=[0],
    #                        skiprows=1, names=["E"])
    for i, file in enumerate(os.listdir('./')):
        match = re.search("\({}\).*\({}\)".format(atom, orbital), file)
        if match:
            if spin == "up":
                df = pd.read_csv(file, delimiter='\s+', header=None, usecols=[1], skiprows=1,
                                 names=["{}".format(file)])
            else:
                df = pd.read_csv(file, delimiter='\s+', header=None, usecols=[1], skiprows=1,
                                 names=["{}".format(file)])


            df_master = pd.concat([df_master, df], axis=1)
        else:
            continue
    return df_master.iloc[:, 0]


def xaxis(tot_pdos_file):
    axis = pd.read_csv(tot_pdos_file, delimiter='\s+', header=None, usecols=[0],
                       skiprows=1, names=["E"])
    return axis


def manual_data(datafile, column):
    axis = pd.read_csv(datafile, delimiter='\s+', header=None, usecols=[column],
                       skiprows=1, names=["E"])
    return axis
