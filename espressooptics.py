import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
import os
import re
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
 
 
 """
#
class EspressoOptics(object):
    """
    A class for calculating and plotting the optical properties from the QuantumEspresso epsilon.x module

    Attributes
    ----------
    Methods
    -------
    label_offset(self, ax, axis="y") : axis_title
        changes the scientific notation of an axis to "axis title (power)"
    energy : df
        loads in the energies over which the dielectric function is calculated
    magnitude: list
        returns a list of the average trace of a 3 x 3 matrix
    refractive_index: df
        returns the frequency dependent refractive index
    extinction: df
        returns the frequency dependent extinction coefficient
    reflectivity: df
        returns the frequency dependent reflectivity
    absorption_coefficient: df
        returns the frequency dependent absorption coefficient
    energy_loss: df
        returns the frequency dependent electron energy loss spectrum
    eps_plot_all_diagram: plot
        plots the dielectric (real and imaginary), extinction coefficient, reflectivity, refractive index and
        energy loss spectrum
    absorbance_plot: plot
        plots the absorbance spectrum against wavelength
    """
    def __init__(self, file_path, system, labels=None, owd=os.getcwd()):
        """
        Parameters
        ----------
        file_path: str
            file path to data
        system: str
            chemical system as formatted in the output file names
        labels: str
            labels of dielectric derived functions plotted
        owd: str
            home directory from which script is run
        """
        os.chdir(file_path)

        self.owd = owd
        self.file_path = file_path
        self.system = system
        pattern = re.compile(r'([0-9])')
        self.formatted_system = pattern.sub(r"$_{\1}$", self.system)

        self.real_file = "epsr_{}.dat".format(self.system)
        self.imag_file = "epsi_{}.dat".format(self.system)

        self.x_axis = self.energy(self.real_file)

        # dielectric function and derivatives
        self.real = self.magnitude(self.real_file)
        self.imag = self.magnitude(self.imag_file)
        self.ext = self.extinction()
        self.abs = self.absorption_coefficient()
        self.refrac = self.refractive_index()
        self.reflect = self.reflectivity()
        self.eels = self.energy_loss()

        #constants
        self.h = 6.63 * math.pow(10, -34)
        self.c = 3.00 * math.pow(10, 10)  # in cm s-1
        self.j = 1.60218 * math.pow(10, -19)

        self.x_axis_nm = (self.h * self.c) * math.pow(10, 7) / (self.x_axis * self.j)

        if labels:
            self.labels = labels
        else:
            self.labels = {"real" : "Real Dielectric",
                           "imag" : "Imaginary Dielectric",
                           "extinct" : "Extinction Coefficient",
                           "reflect" : "Reflectivity",
                           "refrac" : "Refractive Index",
                           "eels" : "Energy Loss"}

        self.features = {"real" : self.real,
                         "imag" : self.imag,
                         "extinct" : self.ext,
                         "refrac": self.refrac,
                         "reflect" : self.reflect,
                         "eels" : self.eels}

        self.units = {"real" : "$\epsilon^{(1)}$($\omega$)",
                      "imag" : "$\epsilon^{(2)}$($\omega$)",
                      "extinct" : "k($\omega$)",
                      "refrac": "n($\omega$)",
                      "reflect" : "R($\omega$)",
                      "eels" : "L($\omega$)"}

        #rc params
        rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
        rc('text', usetex=True)

        os.chdir(self.owd)

    def label_offset(self, ax, axis="y"):
        """
        This function changes label of an axis so that the order of magnitude is included in label

        Taken from https://peytondmurray.github.io/coding/fixing-matplotlibs-scientific-notation/#

        Parameters
        ----------
        ax: ax object
            The axis object to which the offset is applied
        axis: axis
            The specific axis (x or y) to which the offset is applied

        Returns
        -------

        """
        if axis == "y":
            fmt = ax.yaxis.get_major_formatter()
            ax.yaxis.offsetText.set_visible(False)
            set_label = ax.set_ylabel
            label = ax.get_ylabel()

        elif axis == "x":
            fmt = ax.xaxis.get_major_formatter()
            ax.xaxis.offsetText.set_visible(False)
            set_label = ax.set_xlabel
            label = ax.get_xlabel()

        def update_label():
            offset = fmt.get_offset()
            if offset == '':
                set_label("{}".format(label))
            else:
                set_label("{} ({})".format(label, offset))
            return

        ax.callbacks.connect("ylim_changed", update_label)
        ax.callbacks.connect("xlim_changed", update_label)
        ax.figure.canvas.draw()
        update_label(None)
        return

    def energy(self, file):
        """
        Gets the energy values as a pd dataframe for use in plotting
        Parameters
        ----------
        file: str
            file used for obtaining energy
        Returns
        -------
        energies

        """
        energies = pd.read_csv(file, delimiter="\s+", header=None, usecols=["E"], names=["E"], skiprows=2)
        return energies

    def magnitude(self, file):
        """
        This function takes as an input the diagonal parts of a 3 x 3 x 3 matrix
        from the epsilon.x output files and returns the average

        Files consist of 4 columns with headings of:
        energy grid [eV]     x  y  z

        Parameters
        ----------
        file:
            file used for obtaining magnitude

        Returns
        -------
        mags
        """
        df = pd.read_csv(file, delimiter="\s+", header=None, names=["E", "x", "y", "z", "Mag"], skiprows=2)
        for i, row in df.iterrows():
            val = (row["x"] + row["y"] + row["z"]) / 3
            df["Mag"].loc[i] = val
        df["Mag"].to_csv("vals_" + file)
        mags = df["Mag"].values.tolist()
        return mags

    def refractive_index(self):
        """
        This function calculates the refractive index from the complex dielectric function

        Returns
        -------
        df: pandas dataframe
        """
        df = pd.DataFrame(columns=["V"])
        for i, val in enumerate(self.real):
            n_ij = (1 / np.sqrt(2)) * \
                   np.sqrt(np.sqrt(np.square(self.real[i]) + np.square(self.imag[i])) + self.real[i])
            df.loc[i] = n_ij
        return df

    def extinction(self):
        """This function calculates the extinction coefficient from the complex dielectric function

        Returns
        -------
        df: pandas dataframe
        """
        df = pd.DataFrame(columns=["V"])
        for i, val in enumerate(self.real):
            n_ij = (1 / np.sqrt(2)) * np.sqrt(np.sqrt(np.square(self.real[i]) + np.square(self.imag[i])) - self.real[i])
            df.loc[i] = n_ij
        return df

    def reflectivity(self):
        """This calculates the reflectivity from the extinction coefficient and refractive index

        Returns
        -------
        df: pandas dataframe
        """
        df = pd.DataFrame(columns=["V"])
        for i, entry in self.refrac.iterrows():
            # print(refrac.iloc[i])
            reflect = (((self.refrac.iloc[i] - 1) ** 2) + (self.ext.iloc[i])) / (((self.refrac.iloc[i] + 1) ** 2)
                                                                                 + (self.ext.iloc[i]))
            reflect = reflect[0]
            df.loc[i] = reflect
        return df

    def absorption_coefficient(self):
        """This function calculates the absorption coefficient from the extinction coefficient

        Returns
        -------
        df: pandas dataframe
        """
        h = 6.63 * math.pow(10, -34)
        c = 3.00 * math.pow(10, 10)  # cm s-1
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
        """
        This calculates the electron energy loss spectrum from the complex dielectric function

        Returns
        -------
        df: pandas dataframe
        """
        df = pd.DataFrame(columns=["V"])
        for i, val in enumerate(self.real):
            n_ij = (self.imag[i] / (self.imag[i] * self.imag[i] + self.real[i] * self.real[i]))
            df.loc[i] = [n_ij]
        return df

    def eps_plot_all_diagram(self, rows=3, cols=2, ax_a=None, color=(0, 0, 0), title=None):
        """
        plots the optical properties

        Parameters
        ----------
        rows: int
            defines the number of rows of plots
        cols: int
            defined the number of columns of plots
        ax_a: int
            the pyplot axis object
        color: str
            pyplot colour of plot
        title: str
            defined the axis title

        Returns
        -------
        """

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
            self.label_offset(ax, "y")
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

    def absorbance_plot(self, color, title=None, x_lim=(200, 800), y_lim=None, label=None):
        """
        plots the frequency dependent absorption coefficient against wavelength

        Parameters
        ----------
        color: str
            pyplot colour of plot
        title:
            plot title
        x_lim: tuple
            range of x axis
        y_lim: tuple
            range of y axis
        label: str
            plot label for legend

        Returns
        -------

        """
        ax = plt.gca()
        ax.plot(self.x_axis_nm, self.abs, color=color, label=label, lw=1)

        if y_lim:
            ax.set_ylim(y_lim)
        if x_lim:
            ax.set_xlim(x_lim)

        ax.fill_between(self.x_axis_nm["E"], self.abs["V"], color=color, alpha=0.4)

        ax.set_ylabel("Absorption (a.u)", fontsize=16)
        ax.set_xlabel("Wavelength (nm)", fontsize=16)

        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        self.label_offset(ax, "y")

        if title:
            ax.set_title(title, fontsize=20)
        else:
            ax.set_title("Calculated Absorption Spectrum for {}".format(self.formatted_system), fontsize=20)
