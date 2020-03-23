import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
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


def refractive_index(real, imaginary):
    """This function takes as an input the real and imaginary parts of
    the dielectric function and calculates the refractive index"""
    df = pd.DataFrame(columns=["V"])
    for i, val in enumerate(real):
        n_ij = (1 / np.sqrt(2)) * np.sqrt(np.sqrt(np.square(real[i]) + np.square(imaginary[i])) + real[i])
        df.loc[i] = n_ij
    return df


def extinction(real, imaginary):
    """This function takes as an input the real and imaginary parts of
        the dielectric function and calculates the extinction coefficient"""
    df = pd.DataFrame(columns=["V"])
    for i, val in enumerate(real):
        n_ij = (1 / np.sqrt(2)) * np.sqrt(np.sqrt(np.square(real[i]) + np.square(imaginary[i])) - real[i])
        df.loc[i] = n_ij
    return df


def reflectivity(refrac, extinct):
    """This function takes as an input the refractive index and extinction coefficient
     and calculates the reflectivity"""
    df = pd.DataFrame(columns=["V"])
    for i, entry in refrac.iterrows():
        # print(refrac.iloc[i])
        reflect = (((refrac.iloc[i] - 1) ** 2) + (extinct.iloc[i])) / (((refrac.iloc[i] + 1) ** 2) + (extinct.iloc[i]))
        reflect = reflect[0]
        df.loc[i] = reflect
    return df


def absorption_coefficient(extinct, x_axis=None):
    h = 6.63 * math.pow(10, -34)
    c = 3.00 * math.pow(10, 8)
    j = 1.60218 * math.pow(10, -19)
    df = pd.DataFrame(columns=["V"])

    x_axis_eV = pd.DataFrame(data=x_axis, columns=["E"])
    combo = x_axis_eV.join(extinct)
    mult = pd.DataFrame(combo["E"] * combo["V"])
    for i, entry in mult.iterrows():
        alpha = (4 * np.pi * j * entry) / (h * c)
        alpha = alpha[0]
        df.loc[i] = [alpha]
    return df


def energy_loss(real, imaginary):
    """This function takes as an input the real and imaginary parts of
         the dielectric function and calculates the extinction coefficient"""
    df = pd.DataFrame(columns=["V"])
    for i, val in enumerate(real):
        n_ij = (imaginary[i] / (imaginary[i] * imaginary[i] + real[i] * real[i]))
        df.loc[i] = [n_ij]
    return df


def eps_plot_all_diagram(system, data_loc="./", owd=os.getcwd(), rows=3, cols=2, ax_a=None, ax_b=None,
                         color=(0, 0, 0), label="None", title=None, program="espresso", vasp_x_axis=None,
                         vasp_real=None, vasp_imag=None):
    if type(data_loc) == list:
        os.chdir(data_loc[0])
    else:
        os.chdir(data_loc)
    data_file_r = "epsr_{}.dat".format(system)
    data_file_i = "epsi_{}.dat".format(system)

    global fig
    global axs

    if "fig" not in globals():
        fig, axs = plt.subplots(rows, cols, figsize=(6.4, 6.4 * 1.25), dpi=250)

    # initial set up of rows; if statement to allow for single figure, 1D and 2D subplots
    if rows == 1 and cols == 1:
        a = axs
    elif rows == 1 or cols == 1:
        a = axs[ax_a]
    else:
        a = axs

    # program choice
    if program == "espresso":
        x_axis = energy(data_file_r)
        e_re = magnitude("epsr_{}.dat".format(system))
        e_im = magnitude("epsi_{}.dat".format(system))
    elif program == "VASP":
        x_axis = vasp_x_axis
        e_re = vasp_real
        e_im = vasp_imag
    else:
        raise NameError("Program Name specified is not valid: only \"espresso\" and \"VASP\" implemented")

    refrac = refractive_index(e_re, e_im)
    extinct = extinction(e_re, e_im)
    reflect = reflectivity(refrac=refrac, extinct=extinct)
    epsilon = absorption_coefficient(extinct=extinct, x_axis=x_axis)
    loss = energy_loss(e_re, e_im)

    features = [e_re, e_im, epsilon, refrac, reflect, loss]
    labels = ["Real Dielectric", "Imaginary Dielectric", "Absorption",
              "Refractive Index", "Reflectivity", "Energy Loss"]

    # plotting
    for i, ax in enumerate(a.reshape(-1)):
        ax.plot(x_axis, features[i], color=color, label=label, lw=0.75)
        ax.set_title(labels[i], fontsize=10)
        ax.set_xlim(0, 5)
        # needed as e_re and e_im are lists not dfs
        if i < 2:
            ax.set_ylim(0, max(features[i]) * 1.1)
        else:
            ax.set_ylim(0, features[i]["V"].max() * 1.1)

        ax.grid(linestyle="--", color=(0.85, 0.85, 0.85, 0.1))
        ax.axvline(1.8, ls="-", color="black", lw=0.75, zorder=0)
        ax.axvline(3.1, ls="-", color="black", lw=0.75, zorder=0)

    # subtitles

    a[0, 0].legend(loc=3, frameon=False, fontsize=8)

    a[2, 0].set_xlabel("Energy / eV", fontsize=10)
    a[2, 1].set_xlabel("Energy / eV", fontsize=10)

    a[0, 0].set_ylabel("$\epsilon^{(1)}$($\omega$)", fontsize=10)
    a[0, 1].set_ylabel("$\epsilon^{(2)}$($\omega$)", fontsize=10)
    a[1, 0].set_ylabel("I($\omega$)", fontsize=10)
    a[1, 1].set_ylabel("n($\omega$)", fontsize=10)
    a[2, 0].set_ylabel("R($\omega$)", fontsize=10)
    a[2, 1].set_ylabel("L($\omega$)", fontsize=10)

    plt.tight_layout(rect=(0, 0, 1, .925))
    fig.suptitle(title)
    os.chdir(owd)
