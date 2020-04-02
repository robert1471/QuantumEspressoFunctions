import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import pandas as pd
import re
import os
import warnings
import collections

warnings.simplefilter(action='ignore', category=FutureWarning)

class Dos(object):
    def __init__(self, file_path, system, species, orbitals, fermi):

        self.file_path = file_path
        self.system = system
        self.species = species
        self.orbitals = orbitals
        self.fermi = fermi
        self.owd = os.getcwd()

        self.pdos_file = "{}.pdos_tot".format(system)

        self.d_orbitals_up = {
            "dz2" : 3,
            "dzx" : 5,
            "dzy" : 7,
            "dx2-y2" : 9,
            "dxy" : 11
        }
        self.d_orbitals_down = {
            "dz2" : 4,
            "dzx" : 6,
            "dzy" : 8,
            "dx2-y2" : 10,
            "dxy" : 12
        }
        self.d_orbitals_colors = {
            "dz2" : "Red",
            "dzx" : "Lime",
            "dzy" : "DeepPink",
            "dx2-y2" : "orange",
            "dxy" : "cyan"
        }
        self.p_orbitals_up = {
            "pz" : 3,
            "px" : 5,
            "py" : 7,
        }
        self.p_orbitals_down = {
            "pz" : 4,
            "px" : 6,
            "py" : 8,
        }
        self.p_orbitals_colors = {
            "pz" : "Red",
            "px" : "Dark Orange",
            "py" : "Green",
        }

    """Orbital Order
    Order of m-components for each l in the output:

        1, cos(phi), sin(phi), cos(2*phi), sin(2*phi), .., cos(l*phi), sin(l*phi)

    where phi is the polar angle:x=r cos(theta)cos(phi), y=r cos(theta)sin(phi)
    This is determined in file Modules/ylmr2.f90 that calculates spherical harmonics.

    for l=1:
      1 pz     (m=0)
      2 px     (real combination of m=+/-1 with cosine)
      3 py     (real combination of m=+/-1 with sine)

    for l=2:
      1 dz2    (m=0)
      2 dzx    (real combination of m=+/-1 with cosine)
      3 dzy    (real combination of m=+/-1 with sine)
      4 dx2-y2 (real combination of m=+/-2 with cosine)
      5 dxy    (real combination of m=+/-2 with sine)"""

    def pdos(self, species, orbital, spin=None):
        """Function takes all PDOS files of a given atomic species and extracts the PDOS energy values
        returns standard df file that can be used for plotting"""
        df_master = pd.DataFrame()

        for i, file in enumerate(os.listdir('./')):
            match = re.search("\({}\).*\({}\)".format(species, orbital), file)
            if match:
                if spin == 1 or 2:
                    df = pd.read_csv(file, delimiter='\s+', header=None, usecols=[spin], skiprows=1,
                                     names=["{}".format(file)])
                elif spin is None:
                    print("No Spin Specified")
                    df = pd.read_csv(file, delimiter='\s+', header=None, usecols=[1], skiprows=1,
                                     names=["{}".format(file)])
                df_master = pd.concat([df_master, df], axis=1)
            else:
                continue
        df_master = df_master.sum(axis=1)
        df_master = df_master.to_frame(name="E")
        df_master["E"] = df_master["E"].astype(float)

        return df_master["E"].to_numpy()


    def orbital_pdos(self, species, orbital, projected_orbitals=None, spin=None):

        # use all orbitals if none specified
        if projected_orbitals == None and orbital == "d":
            projected_orbitals = self.d_orbitals_up.keys()
        elif projected_orbitals == None and orbital == "p":
            projected_orbitals = self.p_orbitals_up.keys()

        # Define spin and orbitals to be looked at
        if spin == 1:
            if orbital == "d":
                orbital_dict = {k: self.d_orbitals_up[k] for k in projected_orbitals}
                print(orbital_dict)
            else:
                orbital_dict = {k: self.p_orbitals_up[k] for k in projected_orbitals}
        elif spin == 2:
            if orbital == "d":
                orbital_dict = {k: self.d_orbitals_down[k] for k in projected_orbitals}
            else:
                orbital_dict = {k: self.p_orbitals_down[k] for k in projected_orbitals}

        # sorting into df for and sum of each specified orbital
        x = 0
        df_master = pd.DataFrame(columns = [value for value in orbital_dict])
        for i, file in enumerate(os.listdir('./')):
            match = re.search("\({}\).*\({}\)".format(species, orbital), file)
            if match:
                df = pd.read_csv(file, delimiter='\s+', header=None,
                                 usecols=[orbital_dict[value] for value in orbital_dict], skiprows=1,
                                 names=[value for value in orbital_dict])

                # Either append to df for first time or add to existing values after
                if x == 0:
                    df_master = df_master.append(df)
                else:
                    df_master.add(df, fill_value=1)
                x += 1

        return df_master


    def axis(self, row=0):
        # gets the energy axis from any pdos file
        data = pd.read_csv(self.pdos_file, delimiter='\s+', header=None, usecols=[row], skiprows=1, names=["E"])
        data = data[data.E != "*******"].reset_index()
        data["E"] = data["E"].astype(float)

        return data["E"].to_numpy()


    def total_dos(self):
        # total dos data
        y_total_axis_up = self.axis(1)
        y_total_axis_down = self.axis(2)
        df = pd.concat([y_total_axis_up, y_total_axis_down], axis=1, names=["Spin up", "Spin down"])
        return df

    """a.plot(x_total_axis, y_total_axis_up, label="Total", color=tot_color, lw=0.5)
        a.plot(x_total_axis, -y_total_axis_down, label=None, color=tot_color, lw=0.5)
        a.fill_between(x_total_axis["E"], y_total_axis_up["E"], color="black", alpha=0.1)
        a.fill_between(x_total_axis["E"], -y_total_axis_down["E"], color="black", alpha=0.1)"""

    def pdos_plot_diagram(self, spin_polarised=True, ax_title="Untitled", data_loc="./",
                          labela="Spin Up", labelb="Spin Down", scalar=1, totaldos=False, orbital_projected=False,
                          projected_orbitals=None):

        fig, axs = plt.subplots()

        # change to data location
        os.chdir(data_loc)

        # x axis
        x_axis = self.axis(row=0) - self.fermi

        #axs.axvline(0, lw=0.5, color="Red", zorder=10, ls="-")
        for orbital, element in zip(self.orbitals, self.species):
            if spin_polarised == True:
                if orbital_projected == True:

                    # colours
                    if orbital == "d":
                        colors = self.d_orbitals_colors
                    elif orbital == "p":
                        colors = self.p_orbitals_colors

                    # use all orbitals if none specified
                    if projected_orbitals == None and orbital == "d":
                        projected_orbitals = self.d_orbitals_up.keys()
                    elif projected_orbitals == None and orbital == "p":
                        projected_orbitals = self.p_orbitals_up.keys()


                    # data
                    val_up = self.orbital_pdos(species=element, orbital=orbital,
                                               projected_orbitals=projected_orbitals, spin=1)
                    val_down = self.orbital_pdos(species=element, orbital=orbital,
                                                 projected_orbitals=projected_orbitals, spin=2)

                    # catch out any missing values
                    if len(x_axis) < len(val_up):
                        val = len(val_up) - len(x_axis)
                        val_up = val_up[:][val:]
                        val_down = val_down[:][val:]

                    # plot all orbitals
                    for entry in projected_orbitals:
                        axs.plot(x_axis, scalar * val_up[entry], label=entry, color=colors[entry], lw=0.5)
                        axs.plot(x_axis, scalar * -val_down[entry], color=colors[entry], lw=0.5)

                        axs.fill_between(x_axis, scalar * val_up[entry], color=colors[entry], alpha=0.4)
                        axs.fill_between(x_axis, scalar * -val_down[entry], color=colors[entry], alpha=0.4)

                else:
                    val_up = self.pdos(species=element, orbital=orbital, spin=1)
                    val_down = self.pdos(species=element, orbital=orbital, spin=2)

                    # catch out any missing values
                    if len(x_axis) < len(val_down):
                        val = len(val_up) - len(x_axis)
                        val_up = val_up[val:]
                        val_down = val_down[val:]

                    #data = pd.concat([x_axis, val_up, val_down], axis=1, names=["x_axis", "val_up", "val_down"])
                    axs.plot(x_axis, val_up)
                    axs.plot(x_axis, scalar * -val_down,  label=labelb, color="red", lw=0.5)

                    axs.fill_between(x_axis, scalar * val_up, color="red", alpha=0.4)
                    axs.fill_between(x_axis, scalar * -val_down, color="Red", alpha=0.4)
            else:
                val_up = self.pdos(element, orbital)
                axs.plot(x_axis, val_up, label=labela, color=color, lw=0.4)
                axs.fill_between(x_axis["E"], scalar * val_up, color=color, alpha=0.5)

        #axs.set_ylim(-self.y_total_axis_down.max(), self.y_total_axis_up.max())
        axs.set_title(ax_title)
        axs.set_xlim(-20, 20)
        box = axs.get_position()
        axs.set_position([box.x0, box.y0, box.width * 0.95, box.height])

        # Put a legend to the right of the current axis
        axs.legend(loc='center left', bbox_to_anchor=(1, 0, 1, 1))

        # change directory home
        os.chdir(self.owd)
        return axs
