# This was written by Levi Lentz for the Kolpak Group at MIT
# This is distributed under the MIT license
import numpy as np
import os
from matplotlib import pyplot as plt


# This function extracts the high symmetry points from the output of bandx.out
from numpy.core.multiarray import ndarray


def Symmetries(fstring):
    f = open(fstring, 'r')
    x = np.zeros(0)
    for i in f:
        if "high-symmetry" in i:
            x = np.append(x, float(i.split()[-1]))
    f.close()
    return x


# This function takes in the datafile, the fermi energy, the symmetry file, a subplot, and the label
# It then extracts the band data, and plots the bands, the fermi energy in red, and the high symmetry points


def band_plot(datafile, fermi, symmetryfile, subplot, colour, highsympoints, label, labelloc):
    z = np.loadtxt(datafile)  # This loads the bandx.dat.gnu file
    x = np.unique(z[:, 0])  # This is all the unique x-points
    bands = []

    bndl = len(z[z[:, 0] == x[1]])  # This gives the number of bands in the calculation
    Fermi = float(fermi)
    axis = [min(x), max(x), Fermi - Fermi, Fermi + Fermi]

    for i in range(0, bndl):
        bands.append(np.zeros([len(x), 2]))  # This is where we store the bands

    for i in range(0, len(x)):
        sel = z[z[:, 0] == x[i]]  # Here is the energies for a given x

        for j in range(0, bndl):  # This separates it out into a single band
            bands[j][i][0] = x[i]
            bands[j][i][1] = np.multiply(sel[j][1], 1)

    for j, i in enumerate(bands):  # Here we plots the bands
        if j == 0:  # if statement needed so that legend only shows one entry not one for each band
            subplot.plot(i[:, 0], i[:, 1] - Fermi, lw=0.75, color=colour, label=label)
        else:
            subplot.plot(i[:, 0], i[:, 1] - Fermi, lw=0.75, color=colour)

    temp = Symmetries(symmetryfile)
    val = 0
    for j in temp:  # This is the high symmetry lines
        x1 = [j, j]
        x2 = [axis[2] - 50, axis[3] + 50]
        subplot.plot(x1, x2, '--', lw=0.55, color='blue', alpha=0.9)
        subplot.text(j, labelloc, highsympoints[val], va='center', ha='center', fontsize=10)
        val += 1

    subplot.plot([min(x), max(x)], [Fermi - Fermi, Fermi - Fermi], '--', lw=1, color="red", zorder=0)
    subplot.set_xticklabels([])
    subplot.set_ylim(-Fermi, Fermi)
    subplot.set_xlim([axis[0], axis[1]])
    subplot.set_ylabel(r'E - E\textsubscript{f} / eV', fontsize=16)


def bandgap(datafile, fermi):
    z = np.loadtxt(datafile)  # This loads the bandx.dat.gnu file
    x = np.unique(z[:, 0])  # This is all the unique x-points
    bands = []

    bndl = len(z[z[:, 0] == x[1]])  # This gives the number of bands in the calculation
    Fermi = float(fermi)

    for i in range(0, bndl):
        bands.append(np.zeros([len(x), 2]))  # This is where we store the bands

    for i in range(0, len(x)):
        sel = z[z[:, 0] == x[i]]  # Here is the energies for a given x

        for j in range(0, bndl):  # This separates it out into a single band
            bands[j][i][0] = x[i]
            bands[j][i][1] = np.multiply(sel[j][1], 1)

    maxs = []
    mins = []

    for j, i in enumerate(bands):  # Here we plots the bands
        if max(i[:, 1] - Fermi) < 0:  # This calculates the band gap (but does not determine direct or indirect!!)
            maxs.append(max(i[:, 1]) - Fermi)
        if min(i[:, 1] - Fermi) > 0:
            mins.append(min(i[:, 1]) - Fermi)
    bandgap = min(mins) - max(maxs)
    return bandgap


def band_plot_diagram(system, ax_title="Untitled", data_loc="./", high_sym_points=None, fermi=0,
                      color1="red", color2="darkorange", owd=os.getcwd(), rows=1, cols=1, ax_a=0, ax_b=None):
    global fig
    global axs

    # figure set up or check if file already defined to exist
    if "fig" not in globals():
        fig, axs = plt.subplots(rows, cols, figsize=(6.4 * 2, 6.4), dpi=250)

    # if statement to allow for single figure, 1D and 2D subplots
    if rows == 1 and cols == 1:
        a = axs
    elif rows == 1 or cols == 1:
        a = axs[ax_a]
    else:
        a = axs[ax_a, ax_b]

    # change to data location
    os.chdir(data_loc)

    # get band gaps
    band_gap_1 = bandgap("1.{}.bands.dat.gnu".format(system), fermi)
    band_gap_2 = bandgap("2.{}.bands.dat.gnu".format(system), fermi)

    # plot bands
    band_plot("1.{}.bands.dat.gnu".format(system), fermi, "1.bandx.out", a, color1, high_sym_points,
              r"Spin up (E\textsubscript{{g}} = {} eV)".format(round(band_gap_1, 2)), -(fermi + 1))
    band_plot("2.{}.bands.dat.gnu".format(system), fermi, "2.bandx.out", a, color2, high_sym_points,
              r"Spin down (E\textsubscript{{g}} = {} eV)".format(round(band_gap_2, 2)), -(fermi + 1))

    # style
    a.set_ylim((-10.5, 5))
    a.set_xlabel("Wave Vector", fontsize=16, y=-1)
    a.xaxis.set_ticks_position('none')
    a.xaxis.set_label_coords(0.5, -0.065)
    a.legend(framealpha=1)
    a.set_title(ax_title, fontsize=20)

    # change directory home
    os.chdir(owd)
