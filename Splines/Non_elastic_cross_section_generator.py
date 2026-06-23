import math
import matplotlib.pyplot as plt
import numpy as np
import scipy
import scipy.interpolate

from endf_parserpy import EndfParser

parser = EndfParser()
Element = "argon"
CS_file = "Splines/" + Element + "_cs_ne.txt"
enang_file = "Splines/" + Element + "_enang_ne.txt"
endf_dict = parser.parsefile(enang_file) |# Open the energy-angle distribution file and parse it into a dictionary
densities_out = {}
for key, item in endf_dict[6][5]["subsection"][2]["E"].items():
    goinginen = item
    denstemp = 0
    denstemp2 = []
    energytemp = []
    r = []
    for key2 in endf_dict[6][5]["subsection"][2]["b"][key]:
        denstemp += endf_dict[6][5]["subsection"][2]["b"][key][key2][0]
        energytemp.append(endf_dict[6][5]["subsection"][2]["Ep"][key][key2])
        denstemp2.append(denstemp)
        r.append(endf_dict[6][5]["subsection"][2]["b"][key][key2][1])
    for j in range(len(denstemp2)):
        denstemp2[j] *= 1 / denstemp
    densities_out[item] = [energytemp, denstemp2, r]
endf_dict = parser.parsefile(CS_file) # Open the cross section file and parse it into a dictionary
avaenergies = []
avaenergies2 = []
for key in densities_out.keys():
    avaenergies.append(key) # ? never used
    avaenergies2.append(key / 1e6) # Convert energies to MeV for output

# Write out the non-elastic scattering rate to file
f = open("Splines/" + Element + "_ne_rate.txt", "w")
# First line contains the energy values in MeV
f.write(" ".join([str(z / 1e6) for z in endf_dict[3][5]["xstable"]["E"]]) + "\n")
# Second line contains the scattering rates corresponding to the energy values
tmp = " ".join([str(max(z, 0)) for z in endf_dict[3][5]["xstable"]["xs"]])
f.write(tmp + "\n")
f.close()

# Write out the non-elastic scattering rate and angle distribution to file
f = open("Splines/" + Element + "_ne_energyangle_cdf.txt", "w")
f.write(" ".join([str(z) for z in avaenergies2]) + "\n")
for keys in densities_out.keys():
    tmp = " ".join([str(z / 1e6) for z in densities_out[keys][0]])
    f.write(tmp + "\n")
    tmp = " ".join([str(max(z, 0)) for z in densities_out[keys][1]])
    f.write(tmp + "\n")
    tmp = " ".join([str(z) for z in densities_out[keys][2]])
    f.write(tmp + "\n")
f.close()
