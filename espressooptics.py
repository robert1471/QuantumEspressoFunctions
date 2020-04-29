import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
import os
import warnings
from matplotlib import rc

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
class EspressoOptics(object):
    def __init__(self, file_path, system, labels=None, owd=os.getcwd()):

        self.owd = owd
        os.chdir(file_path)
        self.file_path = file_path
        self.system = system

        self.real_file = "epsr_{}.dat".format(self.system)
        self.imag_file = "epsi_{}.dat".format(self.system)

        self.x_axis = self.energy(self.real_file)

        self.real = self.magnitude(self.real_file)
        self.imag = self.magnitude(self.imag_file)
        self.ext = self.extinction()
        self.abs = self.absorption_coefficient()
        self.refrac = self.refractive_index()
        self.reflect = self.reflectivity()
        self.eels = self.energy_loss()

        self.h = 6.63 * math.pow(10, -34)
        self.c = 3.00 * math.pow(10, 10)  # in cm
        self.j = 1.60218 * math.pow(10, -19)
        
        self.x_axis_nm = (self.h * self.c) * math.pow(10, 7) / (self.x_axis * self.j)

        if labels:
            self.labels = labels
        else:
            self.labels = {"real" : "Real Dielectric",
                           "imag" : "Imaginary Dielectric",
                           "extinct" : "Extinction Coefficient",
                           "abs" : "Absorption",
                           "reflect" : "Reflectivity",
                           "refrac" : "Refractive Index",
                           "eels" : "Energy Loss"}

        self.features = {"real" : self.real,
                         "imag" : self.imag,
                         "extinct" : self.ext,
                         "abs" : self.abs,
                         "refrac": self.refrac,
                         "reflect" : self.reflect,
                         "eels" : self.eels}

        self.units = {"real" : "$\epsilon^{(1)}$($\omega$)",
                         "imag" : "$\epsilon^{(2)}$($\omega$)",
                         "extinct" : "k($\omega$)",
                         "abs" : "$\epsilon$($\omega$)",
                         "refrac": "n($\omega$)",
                         "reflect" : "R($\omega$)",
                         "eels" : "L($\omega$)"}

        rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
        rc('text', usetex=True)

        os.chdir(self.owd)

    def energy(self, file):
        """Gets the energy values as a pd dataframe for use in plotting"""
        energies = pd.read_csv(file, delimiter="\s+", header=None, usecols=["E"], names=["E"], skiprows=2)
        return energies

    def magnitude(self, file):
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

    def refractive_index(self):
        """This function takes as an input the real and imaginary parts of
        the dielectric function and calculates the refractive index"""
        df = pd.DataFrame(columns=["V"])
        for i, val in enumerate(self.real):
            n_ij = (1 / np.sqrt(2)) * \
                   np.sqrt(np.sqrt(np.square(self.real[i]) + np.square(self.imag[i])) + self.real[i])
            df.loc[i] = n_ij
        return df

    def extinction(self):
        """This function takes as an input the self.real and imaginary parts of
            the dielectric function and calculates the extinction coefficient"""
        df = pd.DataFrame(columns=["V"])
        for i, val in enumerate(self.real):
            n_ij = (1 / np.sqrt(2)) * np.sqrt(np.sqrt(np.square(self.real[i]) + np.square(self.imag[i])) - self.real[i])
            df.loc[i] = n_ij
        return df

    def reflectivity(self):
        """This function takes as an input the refractive index and extinction coefficient
         and calculates the reflectivity"""
        df = pd.DataFrame(columns=["V"])
        for i, entry in self.refrac.iterrows():
            # print(refrac.iloc[i])
            reflect = (((self.refrac.iloc[i] - 1) ** 2) + (self.ext.iloc[i])) / (((self.refrac.iloc[i] + 1) ** 2)
                                                                                     + (self.ext.iloc[i]))
            reflect = reflect[0]
            df.loc[i] = reflect
        return df

    def absorption_coefficient(self):

        h = 6.63 * math.pow(10, -34)
        c = 3.00 * math.pow(10, 10)
        j = 1.60218 * math.pow(10, -19)

        df = pd.DataFrame(columns=["V"])

        x_axis_eV = pd.DataFrame(data=self.x_axis, columns=["E"])
        combo = x_axis_eV.join(self.ext)
        mult = pd.DataFrame(combo["E"] * combo["V"])
        for i, entry in mult.iterrows():
            alpha = (4 * np.pi * j * entry) / (h * c)
            alpha = alpha[0]
            df.loc[i] = [alpha]
        return df

    def energy_loss(self):
        """This function takes as an input the real and imaginary parts of
             the dielectric function and calculates the extinction coefficient"""
        df = pd.DataFrame(columns=["V"])
        for i, val in enumerate(self.real):
            n_ij = (self.imag[i] / (self.imag[i] * self.imag[i] + self.real[i] * self.real[i]))
            df.loc[i] = [n_ij]
        return df

    def eps_plot_all_diagram(self, rows=3, cols=2, ax_a=None, ax_b=None,
                             color=(0, 0, 0), label="None", title=None, program="espresso", vasp_x_axis=None,
                             vasp_real=None, vasp_imag=None):

        global fig
        global axs

        if "fig" not in globals():
            fig, axs = plt.subplots(rows, cols)

        # initial set up of rows; if statement to allow for single figure, 1D and 2D subplots
        if rows == 1 and cols == 1:
            a = axs
        elif rows == 1 or cols == 1:
            a = axs[ax_a]
        else:
            a = axs

        # plotting
        for feature, ax, label in zip(self.features, a.reshape(-1), self.labels):
            ax.plot(self.x_axis, self.features[feature], color=color, lw=0.75)
            ax.set_title(self.labels[label], fontsize=16)
            ax.set_xlim(0, 5)
            ax.set_ylabel(self.units[feature], fontsize=12)
            # needed as e_re and e_im are lists not dfs
            if feature == "real" or feature == "imag":
                ax.fill_between(self.x_axis["E"].tolist(), self.features[feature], color=color, alpha=0.4)
                ax.set_ylim(0, max(self.features[feature]) * 1.1)
            else:
                ax.fill_between(self.x_axis["E"].tolist(), self.features[feature]["V"].tolist(), color=color, alpha=0.4)
                ax.set_ylim(0, self.features[feature]["V"].max() * 1.1)

            """ax.grid(linestyle="--", color=(0.85, 0.85, 0.85, 0.1))
            ax.axvline(1.8, ls="-", color="black", lw=0.75, zorder=0)
            ax.axvline(3.1, ls="-", color="black", lw=0.75, zorder=0)"""

        # subtitles

        for column in range(cols):
            a[-1, column].set_xlabel("Energy / eV", fontsize=12)

        plt.tight_layout()
        fig.suptitle(title)
        fig.set_size_inches(12 * cols/rows, 12)
        fig.subplots_adjust(left=0.1, wspace=0.25, hspace=0.25)
        os.chdir(self.owd)

    def absorbance_plot(self, color, label=None):
        ax = plt.gca()
        ax.plot(self.x_axis_nm, self.abs, color=color, label=label, lw=1)
        ax.fill_between(self.x_axis_nm["E"], self.abs["V"], color=color, alpha=0.4)
        ax.set_xlim(200, 800)
        ax.set_ylabel("Absorption (a.u)", fontsize=16)
        ax.set_xlabel("Wavelength (nm)", fontsize=16)
        ax.set_title("Calculated Absorption Spectrum for {}".format(self.system), fontsize=20)
