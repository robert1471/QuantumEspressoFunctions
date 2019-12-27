import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import pandas as pd
import re
import os
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)

"""This script contains functions to calculate, save, process and plot output from the QE epsilon.x postprocessing
 module. Properties considered are:
 
 Band Structure
 DOS
 PDOS
 Real and Imaginary parts of the dielectric function
 Refractive Index
 Reflectivity
 Extinction Coefficient
 
 
 TO DO:
    • Band Structure
    • DOS
 """

system = ""
e_1 = "real part of epsilon"
e_2 = "imaginary part of epsilon"


# i = "imaginary number"
# dielectric_function = e_1 + ie_2
# reflectivity = np.square(1 - np.sqrt(np.square(e_re[i]) + np.square(e_im[i]))) / np.square(1 + e_1) + np.square(e_2)

def energy(file):
    """Gets the energy values as a pd dataframe for use in plotting"""
    energies = pd.read_csv(file, delimiter="\s+", header=None, usecols=["E"], names=["E"], skiprows=2)
    return energies


def magnitude(file):
    """This function takes as an input the diagonal parts of a 3 x 3 x 3 matrix
    from the epsilon.x output files and returns the average

    Files consist of 4 columns with headings of:

    energy grid [eV]     x  y  z
    """
    df = pd.read_csv(file, delimiter="\s+", header=None, names=["E", "x", "y", "z", "Mag"], skiprows=2)
    for i, row in df.iterrows():
        val = (row["x"] + row["y"] + row["z"]) / 3
        df["Mag"].loc[i] = val
    df["Mag"].to_csv("vals_" + file)
    mags = df["Mag"].values.tolist()
    return mags


epsr_file = "epsr_MgAl2O4.dat"
epsi_file = "epsi_MgAl2O4.dat"


def refractive_index(epsr_file, epsi_file):
    """This function takes as an input the real and imaginary parts of
    the dielectric function and calculates the refractive index"""
    e_re = magnitude(epsr_file)
    e_im = magnitude(epsi_file)
    df = pd.DataFrame(columns=["V"])
    for i, val in enumerate(e_re):
        n_ij = (1 / np.sqrt(2)) * np.sqrt(np.sqrt(np.square(e_re[i]) + np.square(e_im[i])) + e_re[i])
        df.loc[i] = [n_ij]
    return df


def extinction(epsr_file, epsi_file):
    """This function takes as an input the real and imaginary parts of
        the dielectric function and calculates the extinction coefficient"""
    e_re = magnitude(epsr_file)
    e_im = magnitude(epsi_file)
    df = pd.DataFrame(columns=["V"])
    for i, val in enumerate(e_re):
        n_ij = (1 / np.sqrt(2)) * np.sqrt(np.sqrt(np.square(e_re[i]) + np.square(e_im[i])) - e_re[i])
        df.loc[i] = [n_ij]
    return df


def reflectivity(epsr_file, epsi_file):
    """This function takes as an input the refractive index and extinction coefficient
     and calculates the reflectivity"""
    refrac = refractive_index(epsr_file, epsi_file)
    ext = extinction(epsr_file, epsi_file)
    df = pd.DataFrame(columns=["V"])
    for i, entry in refrac.iterrows():
        # print(refrac.iloc[i])
        reflect = (((refrac.iloc[i] - 1) ** 2) + (ext.iloc[i])) / (((refrac.iloc[i] + 1) ** 2) + (ext.iloc[i]))
        df.loc[i] = [reflect]
    return df


def energy_loss(epsr_file, epsi_file):
    """This function takes as an input the real and imaginary parts of
         the dielectric function and calculates the extinction coefficient"""
    e_re = magnitude(epsr_file)
    e_im = magnitude(epsi_file)
    df = pd.DataFrame(columns=["V"])
    for i, val in enumerate(e_re):
        n_ij = (e_im[i] / (e_im[i] * e_im[i] + e_re[i] * e_re[i]))
        df.loc[i] = [n_ij]
    return df


plot_tile = "Absorption Spectrum for MgAl$_2$O$_4$"


def plotter(xlabel, ylabel, plot_title):
    """A general plotting function"""

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(plot_title)

    plt.plot(energy("epsr_MgAl2O4.dat"), magnitude("epsr_MgAl2O4.dat"), color="red")
    # plt.hlines(y=0, xmin=0, xmax=25, colors="red", linestyles="--")
    # plt.savefig("EnergyLoss.pdf")
    plt.show()


os.chdir("./MgAl2O4/EPSILON")
fig, ax = plt.subplots(3, 2, sharex=True)
ax[0, 0].plot(energy("epsr_MgAl2O4.dat"), magnitude("epsr_MgAl2O4.dat"), color='Red')  # Real
ax[0, 0].title.set_text('Real part of Dielectric Function')
ax[0, 0].set_ylabel("$\epsilon_{1}(\omega)$")

ax[0, 1].plot(energy("epsr_MgAl2O4.dat"), magnitude("epsi_MgAl2O4.dat"), color='Blue')  # Imaginary
ax[0, 1].title.set_text('Imaginary part of Dielectric Function')
ax[0, 1].set_ylabel("$\epsilon_{2}(\omega)$")

ax[1, 0].plot(energy("epsr_MgAl2O4.dat"), refractive_index("epsr_MgAl2O4.dat", "epsi_MgAl2O4.dat"), color='Green')  # Refractive Index
ax[1, 0].title.set_text('Refractive Index')
ax[1, 0].set_ylabel("$n(\omega)$")

ax[1, 1].plot(energy("epsr_MgAl2O4.dat"), extinction("epsr_MgAl2O4.dat", "epsi_MgAl2O4.dat"), color='Orange')  # Absorption Coefficient
ax[1, 1].title.set_text('Absorption Coefficient')
ax[1, 1].set_ylabel("a.u.")

ax[2, 0].plot(energy("epsr_MgAl2O4.dat"), reflectivity("epsr_MgAl2O4.dat", "epsi_MgAl2O4.dat"), color='Gray')  # Reflectivity
ax[2, 0].title.set_text('Reflectivity')
ax[2, 0].set_ylabel("$R(\omega)$")
ax[2, 0].set_xlabel("Energy / eV")

ax[2, 1].plot(energy("epsr_MgAl2O4.dat"), energy_loss("epsr_MgAl2O4.dat", "epsi_MgAl2O4.dat"), color='Black')  # Energy Loss Spectrum
ax[2, 1].title.set_text('Energy Loss Spectrum')
ax[2, 1].set_ylabel("$L(\omega)$")
ax[2, 1].set_xlabel("Energy / eV")
plt.show()


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
            df_master = pd.concat([df_master, df], axis=1)
        else:
            continue
    return df_master


# def PDOS_axis(subplotx, subploty, subploti, system, species, orbitals):


def PDOS_plot():
    Os = PDOS("O", "s")
    Op = PDOS("O", "p")
    Od = PDOS("O", "d")

    Als = PDOS("Al", "s")
    Alp = PDOS("Al", "p")
    Ald = PDOS("Al", "d")

    Mgs = PDOS("Mg", "s")
    Mgp = PDOS("Mg", "p")
    Mgd = PDOS("Mg", "d")

    x_axis = pd.read_csv('./PDOS/MgAl2O4.pdos_tot', delimiter='\s+', header=None, usecols=[0],
                         skiprows=1, names=["E"])

    plt.subplot(3, 1, 1)

    plt.plot(x_axis, Os.iloc[:, 0], label='O$_s$', lw=1)
    plt.plot(x_axis, Op.iloc[:, 0], label='O$_p$', lw=1)
    plt.plot(x_axis, Od.iloc[:, 0], label='O$_d$', lw=1)
    plt.axvline(5.5935, color="red", label="E$_{Fermi}$", lw=1)
    plt.title("Partial DOS plot for MgAl$_2$O$_4$", fontsize=20)
    plt.legend()

    plt.subplot(3, 1, 2)

    plt.plot(x_axis, Als.iloc[:, 0], label='Al$_s$', lw=1)
    plt.plot(x_axis, Alp.iloc[:, 0], label='Al$_p$', lw=1)
    plt.plot(x_axis, Ald.iloc[:, 0], label='Al$_d$', lw=1)
    plt.axvline(5.5935, color="red", label="E$_{Fermi}$", lw=1)
    plt.legend()

    plt.subplot(3, 1, 3)

    plt.plot(x_axis, Mgs.iloc[:, 0], label='Mg$_s$', lw=1)
    plt.plot(x_axis, Mgp.iloc[:, 0], label='Mg$_p$', lw=1)
    plt.plot(x_axis, Mgd.iloc[:, 0], label='Mg$_d$', lw=1)

    plt.axvline(5.5935, color="red", label="E$_{Fermi}$", lw=1)
    plt.xlabel("Energy / eV", fontsize=16)
    plt.legend()
    plt.show()


def dos_plot():
    df = pd.read_csv("MgAl2O4.dos", delimiter='\s+')
    print(df[['E', 'dos']])
    df['E'] = df['E'] - 5.593
    plt.plot(df['E'], df['dos'], lw=0.75, color='green')
    plt.axvline(x=0, ymin=0, ymax=45, linestyle='-', lw=0.55, color='red')
    plt.xlabel("Energy / eV")
    plt.title("DOS of MgAl$_2$O$_4$")
    plt.savefig("dos.pdf")
