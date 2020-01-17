from matplotlib import pyplot as plt
#import espressooptics as eo
import espressobands as eb
import os



os.chdir("MnFeO")

fig, axs = plt.subplots(1, 2)

high_symmetry_points = ["$\Gamma$", "X", "W", "K", "$\Gamma$", "L", "U", "W", "L", "K", "U", "X"]

bandgap = eb.bandgap("NC.1.MnFeO22.bands.dat.gnu", "10.5424", axs[0])
eb.bandplot("NC.1.MnFeO22.bands.dat.gnu", "10.5424", "NC.1.bandx.out", axs[0], "Green", high_symmetry_points, "Spin up")
eb.bandplot("NC.2.MnFeO22.bands.dat.gnu", "10.5424", "NC.2.bandx.out", axs[0], "Red", high_symmetry_points, "Spin down")

eb.bandplot("USPP.1.MnFeO22.bands.dat.gnu", "10.2743", "USPP.1.bandx.out", axs[1], "Green", high_symmetry_points,
            "Spin Up (E$_\mathnormal{Cut}$ = {})".format(bandgap))
eb.bandplot("USPP.2.MnFeO22.bands.dat.gnu", "10.2743", "USPP.2.bandx.out", axs[1], "Red", high_symmetry_points, "Spin down")
#axs.legend(fontsize=16)
#axs.set_title("Spin polarised Band Structures of MnFeO", fontsize=20)
axs[0].set_title("Norm-Conserving", fontsize=20)
axs[1].set_title("Ultra-soft", fontsize=20)
plt.show()