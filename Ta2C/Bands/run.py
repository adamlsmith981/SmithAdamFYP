#!/usr/bin/env python
from Bands import *

datafile='Bands_results/Ta2C.Bandx.dat.gnu'
fermi = -0.2033 # get from nscf / scf step out
symmetryfile='Bands_results/Ta2C.bandx.out'
bool_shift_efermi= True
fig, ax = plt.subplots()

#bndplot(datafile,fermi,symmetryfile,ax)
bndplot(datafile,fermi,symmetryfile,ax,shift_fermi=0,\
color='black',linestyle='solid',name_k_points=['Γ','M','K','Γ'],legend='Ta2C, PBE')


fig.savefig("ta2c_bands.png")
plt.show()