from matplotlib import pyplot as plt
#import espressooptics as eo
import espressobands as eb
import os

os.chdir("MnFeO")

fig, axs = plt.subplots()

high_symmetry_points = ["$\Gamma$", "X", "W", "K", "$\Gamma$", "L", "U", "W", "L", "K", "U", "X"]


eb.bandplot("1.MnFeO22.bands.dat.gnu", "10.2743", "1.bandx.out", axs, "Green", high_symmetry_points, "Spin up")
eb.bandplot("2.MnFeO22.bands.dat.gnu", "10.2743", "2.bandx.out", axs, "Red", high_symmetry_points, "Spin down")
axs.legend(fontsize=16)
axs.set_title("Spin polarised Band Structures of MnFeO", fontsize=20)
plt.show()
