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
