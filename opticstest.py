from matplotlib import pyplot as plt
import espressooptics as eo
import espressobands as eb
import os

os.chdir("results.hubbard/hub-4")

fig, axs = plt.subplots()

high_symmetry_points = []
for i in range(1, 18):
    high_symmetry_points.append(str(i))
print(high_symmetry_points)

eb.bandplot("NiO.bands.dat.gnu", "5", "bandx.out", axs, "Green", high_symmetry_points
            , "If this worked it would be a miracle", )
plt.show()
#x = eo.energy("epsi_MgAl2O4.dat")
#y = eo.extinction("epsr_MgAl2O4.dat", "epsi_MgAl2O4.dat")

#print(x, y)


#plt.plot(x, y)

#plt.show()
