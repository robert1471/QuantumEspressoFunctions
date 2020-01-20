# This was written by Levi Lentz for the Kolpak Group at MIT
# This is distributed under the MIT license
import numpy as np


# This function extracts the high symmetry points from the output of bandx.out

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


def bandplot(datafile, fermi, symmetryfile, subplot, colour, highsympoints, label):
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
        subplot.text(j, - Fermi - 0.5, highsympoints[val], va='center', ha='center', fontsize=12)
        val += 1

    subplot.plot([min(x), max(x)], [Fermi - Fermi, Fermi - Fermi], '--', lw=1, color="red", zorder=0)
    subplot.set_xticklabels([])
    subplot.set_ylim(-Fermi, Fermi)
    subplot.set_xlim([axis[0], axis[1]])
    subplot.set_ylabel('E - E$_\mathregular{f}$ / eV', fontsize=16)


def bandgap(datafile, fermi, subplot):
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
