import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import pandas as pd
import re
import os
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)

def PDOS(atom, orbital):
    """Function takes all PDOS files of a given atomic species and extracts the PDOS energy values
    returns standard df file that can be used for plotting"""
    df_master = pd.DataFrame()
    # df_master = pd.read_csv('./PDOS/MgAl2O4.pdos_tot', delimiter='\s+', header=None, usecols=[0],
    #                        skiprows=1, names=["E"])
    for i, file in enumerate(os.listdir('./PDOS')):
        match = re.search("\({}\).*\({}\)".format(atom, orbital), file)
        if match:
            df = pd.read_csv('./PDOS/' + file, delimiter='\s+', header=None, usecols=[1], skiprows=1,
                             names=["{}".format(file)])
            if atom == 'Ni1' or atom == 'Ni2':
                df = df * 2
            if atom == 'Ni2':
                df = df * -1
            df_master = pd.concat([df_master, df], axis=1)
        else:
            continue
    print(df_master)
    return df_master


def PDOS_plot():
    Os = PDOS("O", "s")
    Op = PDOS("O", "p")
    Od = PDOS("O", "d")

    Ni1s = PDOS("Ni1", "s")
    Ni1p = PDOS("Ni1", "p")
    Ni1d = PDOS("Ni1", "d")

    Ni2s = PDOS("Ni2", "s")
    Ni2p = PDOS("Ni2", "p")
    Ni2d = PDOS("Ni2", "d")

    x_axis = pd.read_csv('./PDOS/NiO.pdos_tot', delimiter='\s+', header=None, usecols=[0],
                         skiprows=1, names=["E"])

    y_axis1 = pd.read_csv('./PDOS/NiO.pdos_tot', delimiter='\s+', header=None, usecols=[1],
                         skiprows=1, names=["U"])

    y_axis2 = pd.read_csv('./PDOS/NiO.pdos_tot', delimiter='\s+', header=None, usecols=[1],
                         skiprows=1, names=["D"])

    plt.plot(x_axis, y_axis1)
    plt.plot(x_axis, y_axis2)
    '''plt.subplot(3, 1, 1)

    plt.plot(x_axis, Os.iloc[:, 0], label='O$_s$', lw=1)
    plt.plot(x_axis, Op.iloc[:, 0], label='O$_p$', lw=1)
    #plt.plot(x_axis, Od.iloc[:, 0], label='O$_d$', lw=1)
    plt.axvline(10.837, color="red", label="E$_{Fermi}$", lw=1)
    plt.ylim(0, 5)
    plt.title("Partial DOS plot for NiO", fontsize=20)
    plt.legend()

    plt.subplot(3, 1, 2)

    plt.plot(x_axis, Ni1s.iloc[:, 0], label='Ni1$_s$', lw=1)
    #plt.plot(x_axis, Ni1p.iloc[:, 0], label='Ni1$_p$', lw=1)
    plt.plot(x_axis, Ni1d.iloc[:, 0], label='Ni1$_d$', lw=1)
    plt.axvline(10.837, color="red", label="E$_{Fermi}$", lw=1)
    plt.legend()

    plt.subplot(3, 1, 3)

    plt.plot(x_axis, Ni2s.iloc[:, 0], label='Ni2$_s$', lw=1)
    #plt.plot(x_axis, Ni2p.iloc[:, 0], label='Ni2$_p$', lw=1)
    plt.plot(x_axis, Ni2d.iloc[:, 0], label='Ni2$_d$', lw=1)'''

    plt.axvline(10.837, color="red", label="E$_{Fermi}$", lw=1)
    plt.xlabel("Energy / eV", fontsize=16)
    plt.legend()
    plt.show()

PDOS_plot()