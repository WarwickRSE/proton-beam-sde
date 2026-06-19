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
endf_dict = parser.parsefile(enang_file)
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
endf_dict = parser.parsefile(CS_file)
avaenergies = []
avaenergies2 = []
for key in densities_out.keys():
    avaenergies.append(key)
    avaenergies2.append(key / 1e6)
f = open("Splines/" + Element + "_ne_rate.txt", "w")
f.write(" ".join([str(z / 1e6) for z in endf_dict[3][5]["xstable"]["E"]]) + "\n")
tmp = " ".join([str(max(z, 0)) for z in endf_dict[3][5]["xstable"]["xs"]])
f.write(tmp + "\n")
f.close()
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

""" 
Claude Haiku 4.5:
This script processes nuclear physics data from ENDF (Evaluated Nuclear Data File) format files to extract and organize non-elastic cross-section information for argon. Here's what it does:

Main steps:

Parse ENDF data — Uses EndfParser to read two argon-related ENDF files:

Energy-angle correlation data (argon_enang_ne.txt)
Cross-section data (argon_cs_ne.txt)
Process energy-angle distributions — For each incident energy, it:

Accumulates secondary particle branch ratios from the b field
Collects corresponding outgoing energies and angular distributions
Normalizes the branch ratios to create a probability distribution (denstemp2)
Generate output files — Writes three formatted text files:

argon_ne_rate.txt — Cross-section values (in MeV) for each energy bin
argon_ne_energyangle_cdf.txt — Cumulative density functions (CDFs) for outgoing particle energies and angles at each incident energy, used for Monte Carlo sampling in simulations
The data is converted to MeV units (dividing by 1e6) and negative values are clamped to zero. This preprocessed data likely feeds into a particle transport simulation (like the Geant4 code in your workspace).

 """