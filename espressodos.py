import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import pandas as pd
import re
import os
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)


def pdos(atom, orbital, spin=None):
    """Function takes all PDOS files of a given atomic species and extracts the PDOS energy values
    returns standard df file that can be used for plotting"""
    df_master = pd.DataFrame()
    # df_master = pd.read_csv('./PDOS/MgAl2O4.pdos_tot', delimiter='\s+', header=None, usecols=[0],
    #                        skiprows=1, names=["E"])
    for i, file in enumerate(os.listdir('./')):
        match = re.search("\({}\).*\({}\)".format(atom, orbital), file)
        if match:
            if spin == 1:
                df = pd.read_csv(file, delimiter='\s+', header=None, usecols=[1], skiprows=1,
                                 names=["{}".format(file)])
            elif spin == 2:
                df = pd.read_csv(file, delimiter='\s+', header=None, usecols=[2], skiprows=1,
                                 names=["{}".format(file)])
            elif spin is None:
                print("No Spin Specified")
                df = pd.read_csv(file, delimiter='\s+', header=None, usecols=[1], skiprows=1,
                                 names=["{}".format(file)])
            else:
                print("No Spin Specified; Are you sure this is right?")
                df = pd.read_csv(file, delimiter='\s+', header=None, usecols=[1], skiprows=1,
                                 names=["{}".format(file)])
            df_master = pd.concat([df_master, df], axis=1)
        else:
            continue
    df_master = df_master.sum(axis=1)
    print(df_master)
    return df_master


def axis(pdos_file, row=0):
    # gets the energy axis from any pdos file
    data = pd.read_csv(pdos_file, delimiter='\s+', header=None, usecols=[row],
                       skiprows=1, names=["E"])
    return data


def pdos_plot_diagram(system, element=None, orbitals=[], spin_polarised=True, ax_title="Untitled", data_loc="./",
                      owd=os.getcwd(), fermi=0, rows=1, cols=1, ax_a=None, ax_b=None, color=(0, 0, 0),
                      tot_color=(0, 0, 0),
                      labela="Spin Up", labelb="Spin Down", total_dos=False):
    global fig
    global axs

    # change to data location
    os.chdir(data_loc)

    # total dos data
    x_total_axis = axis("{}.pdos_tot".format(system)) - fermi
    y_total_axis_up = axis("{}.pdos_tot".format(system), 1)
    y_total_axis_down = axis("{}.pdos_tot".format(system), 2)

    # figure set up or check if file already defined to exist and to plot total dos only once
    if "fig" not in globals():
        fig, axs = plt.subplots(rows, cols, figsize=(6.4 * 2, 6.4), dpi=250)

        # initial set up of rows; if statement to allow for single figure, 1D and 2D subplots
        if rows == 1 and cols == 1:
            a = axs
        elif rows == 1 or cols == 1:
            a = axs[ax_a]
        else:
            a = axs[ax_a, ax_b]

        a.axvline(0, lw=0.5, color="Red", zorder=10, ls="-")
        # plot total dos
        if total_dos:
            a.plot(x_total_axis, y_total_axis_up, label="Total", color=tot_color, lw=0.5)
            a.plot(x_total_axis, -y_total_axis_down, label=None, color=tot_color, lw=0.5)
            a.fill_between(x_total_axis["E"], y_total_axis_up["E"], color="black", alpha=0.1)
            a.fill_between(x_total_axis["E"], -y_total_axis_down["E"], color="black", alpha=0.1)

    else:
        # if statement to allow for single figure, 1D and 2D subplots
        if rows == 1 and cols == 1:
            a = axs
        elif rows == 1 or cols == 1:
            a = axs[ax_a]
        else:
            a = axs[ax_a, ax_b]
        a.axvline(0, lw=0.5, color="Red", zorder=10, ls="-")

    for orbital in orbitals:
        if spin_polarised is True:
            val_up = pdos(element, orbital, 1)
            val_down = pdos(element, orbital, 2)
            a.plot(x_total_axis, val_up, label=labela, color=color, lw=0.5)
            a.plot(x_total_axis, -val_down,  label=labelb, color=color, lw=0.5)
            a.fill_between(x_total_axis["E"], val_up, color=color, alpha=0.4)
            a.fill_between(x_total_axis["E"], -val_down, color=color, alpha=0.4)
        else:
            val_up = pdos(element, orbital)
            a.plot(x_total_axis, val_up, label=labela, color=color, lw=0.4)
            a.fill_between(x_total_axis["E"], val_up, color=color, alpha=0.5)

    a.set_ylim(-y_total_axis_down["E"].max(), y_total_axis_up["E"].max())
    a.set_title(ax_title)
    a.set_xlim(-20, 20)
    box = a.get_position()
    a.set_position([box.x0, box.y0, box.width * 0.95, box.height])

    # Put a legend to the right of the current axis
    a.legend(loc='center left', bbox_to_anchor=(1, 0, 1, 1))
    # change directory home
    os.chdir(owd)
    return a
