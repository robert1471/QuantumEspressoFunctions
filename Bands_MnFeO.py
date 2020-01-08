# This was written by Levi Lentz for the Kolpak Group at MIT
# This is distributed under the MIT license
import numpy as np
import matplotlib.pyplot as plt
import os


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

syms = ["$\Gamma$", "X", "W", "K", "$\Gamma$", "L", "U", "W", "L", "K", "U", "X", "$\Gamma$"]


def bndplot(datafile, fermi, symmetryfile, subplot, species, colour):
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

    for i in bands:  # Here we plots the bands
        print(i[:, 0])
        subplot.plot(i[:, 0], i[:, 1] - Fermi, lw=0.75, color=colour)

    temp = Symmetries(symmetryfile)
    val = 0
    for j in temp:  # This is the high symmetry lines
        x1 = [j, j]
        x2 = [axis[2] - 50, axis[3] + 50]
        subplot.plot(x1, x2, '--', lw=0.55, color='blue', alpha=0.9, label=" ")
        #subplot.text(j, axis[2] - 2 * Fermi, syms[val], va='center', ha='center', fontsize=12)
        val += 1

    subplot.plot([min(x), max(x)], [Fermi - Fermi, Fermi - Fermi], '--', lw=1, color="red", zorder=0)
    subplot.set_xticklabels([])
    subplot.set_ylim([axis[2] - (2 * Fermi) + 1, axis[3] + 2 * Fermi])
    subplot.set_xlim([axis[0], axis[1]])
    subplot.text((axis[1] - axis[0]) / 2.0, axis[3] + 2 * Fermi + 1, "Band Structure of {}".format(species),
                 va='center', ha='center', fontsize=20)
    subplot.set_ylabel('E - E$_{fermi}$ / eV')
    subplot.set_ylim(5, 12)


fig, axs = plt.subplots(1, 1)

os.chdir("./MnFe2O4")
bndplot("MnFeO22.bands.dat.gnu", "5.5934", "bandx.out", axs, "NiO", "black")

'''print([x[0] for x in os.walk("./")])

for i, x in enumerate(os.walk("./")):
    if x[0] == './':
        continue
    else:
        print(x[0])
        os.chdir(x[0])
        bndplot("NiO.bands.dat.gnu", "5.5934", "bandx.out", axs[i], "NiO", "black")
        os.chdir("../")'''
#bndplot("MgAl2O4.bands.dat.gnu", "5.5934", "bandx.out", axs, "MgAl$_2$O$_4$", "black", 1)
#bndplot("MgAl2O4.bands.dat.gnu.14", "5.5934", "bandx.out", axs, "MgAl$_2$O$_4$", "red", 1)


plt.show()