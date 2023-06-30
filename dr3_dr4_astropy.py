import healpy as hp
import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.constants import c
from astropy.constants import h
from astropy.constants import k_B
import astropy.units as u
import cosmoglobe as cg

T_cmb = 2.7255*u.K

# Load data                                                                                                                                                                                                        
#--------------------------------------------                                                                                                                                                                      
print("Load data")
map353_dr3 = hp.read_map('/mn/stornext/d5/data/daniher/dang_runs/hi_fit/data/planck/dr3/HFI_map_353-psb_IQU_n1024_uKcmb.fits')*u.uK
map545_dr3 = hp.read_map('/mn/stornext/d5/data/daniher/dang_runs/hi_fit/data/planck/dr3/HFI_map_545_I_n1024_MJy.fits')*u.MJy/u.sr
map857_dr3 = hp.read_map('/mn/stornext/d5/data/daniher/dang_runs/hi_fit/data/planck/dr3/HFI_map_857_I_n1024_MJy.fits')*u.MJy/u.sr

map353_dr4 = hp.read_map('/mn/stornext/d5/data/daniher/dang_runs/hi_fit/data/planck/dr4/HFI_map_353-BPassCorrected_IQU_n1024_uK_CMB.fits')*u.uK
map545_dr4 = hp.read_map('/mn/stornext/d5/data/daniher/dang_runs/hi_fit/data/planck/dr4/HFI_map_545-BPassCorrected_IQU_n1024_uK_CMB.fits')*u.uK
map857_dr4 = hp.read_map('/mn/stornext/d5/data/daniher/dang_runs/hi_fit/data/planck/dr4/HFI_map_857-BPassCorrected_IQU_n1024_uK_CMB.fits')*u.uK

cmb_map = hp.read_map('/mn/stornext/u3/hke/BeyondPlanck/data/delivery/v10.00/v2/BP_cmb_resamp_I_n1024_v2.fits')*u.uK

bp_353_dr4 = np.loadtxt('/mn/stornext/d5/data/daniher/dang_runs/hi_fit/data/planck/dr4/bp_RIMO_R4.00_353.dat').T
bp_545_dr4 = np.loadtxt('/mn/stornext/d5/data/daniher/dang_runs/hi_fit/data/planck/dr4/bp_RIMO_R4.00_545.dat').T
bp_857_dr4 = np.loadtxt('/mn/stornext/d5/data/daniher/dang_runs/hi_fit/data/planck/dr4/bp_RIMO_R4.00_857.dat').T


# Find monopoles for each map                                                                                                                                                                                      
#--------------------------------------------                                                                                                                                                                      
res353_dr3,mono353_dr3,dip353_dr3 = hp.remove_dipole(map353_dr3,fitval=True)
res545_dr3,mono545_dr3,dip545_dr3 = hp.remove_dipole(map545_dr3,fitval=True)
res857_dr3,mono857_dr3,dip857_dr3 = hp.remove_dipole(map857_dr3,fitval=True)

res353_dr4,mono353_dr4,dip353_dr4 = hp.remove_dipole(map353_dr4,fitval=True)
res545_dr4,mono545_dr4,dip545_dr4 = hp.remove_dipole(map545_dr4,fitval=True)
res857_dr4,mono857_dr4,dip857_dr4 = hp.remove_dipole(map857_dr4,fitval=True)

mono353_dr3 = mono353_dr3*u.uK
mono545_dr3 = mono545_dr3*u.MJy/u.sr
mono857_dr3 = mono857_dr3*u.MJy/u.sr

mono353_dr4 = mono353_dr4*u.uK
mono545_dr4 = mono545_dr4*u.uK
mono857_dr4 = mono857_dr4*u.uK

# Raw difference                                                                                                                                                                                                   
diff_353_mono = (map353_dr3 - mono353_dr3) - (map353_dr4-cmb_map-mono353_dr4)
diff_545_mono = ((map545_dr3-mono545_dr3).to(u.uK, equivalencies=u.thermodynamic_temperature(545*u.GHz)) - (map545_dr4-cmb_map-mono545_dr4))
diff_857_mono = ((map857_dr3-mono857_dr3).to(u.uK, equivalencies=u.thermodynamic_temperature(857*u.GHz)) - (map857_dr4-cmb_map-mono857_dr4))

diff_353 = ((map353_dr3) - (map353_dr4-cmb_map))
diff_545 = ((map545_dr3).to(u.uK, equivalencies=u.thermodynamic_temperature(545*u.GHz)) - (map545_dr4-cmb_map))
diff_857 = ((map857_dr3).to(u.uK, equivalencies=u.thermodynamic_temperature(857*u.GHz)) - (map857_dr4-cmb_map))

# Ratios                                                                                                                                                                                                           
ratio_353_mono = ((map353_dr3 - mono353_dr3)/ (map353_dr4-cmb_map-mono353_dr4))
ratio_545_mono = ((map545_dr3-mono545_dr3).to(u.uK, equivalencies=u.thermodynamic_temperature(545*u.GHz)) / (map545_dr4-cmb_map-mono545_dr4))
ratio_857_mono = ((map857_dr3-mono857_dr3).to(u.uK, equivalencies=u.thermodynamic_temperature(857*u.GHz)) / (map857_dr4-cmb_map-mono857_dr4))

ratio_353 = (map353_dr3 / (map353_dr4-cmb_map))
ratio_545 = ((map545_dr3.to(u.uK, equivalencies=u.thermodynamic_temperature(545*u.GHz))) / (map545_dr4-cmb_map))
ratio_857 = ((map857_dr3.to(u.uK, equivalencies=u.thermodynamic_temperature(857*u.GHz))) / (map857_dr4-cmb_map))

# Write out the ratio and difference maps
hp.write_map('diff_353_mono.fits',diff_353_mono)
hp.write_map('diff_545_mono.fits',diff_545_mono)
hp.write_map('diff_857_mono.fits',diff_857_mono)

hp.write_map('diff_353.fits',diff_353)
hp.write_map('diff_545.fits',diff_545)
hp.write_map('diff_857.fits',diff_857)

hp.write_map('ratio_353_mono.fits',ratio_353_mono)
hp.write_map('ratio_545_mono.fits',ratio_545_mono)
hp.write_map('ratio_857_mono.fits',ratio_857_mono)

hp.write_map('ratio_353.fits',ratio_353)
hp.write_map('ratio_545.fits',ratio_545)
hp.write_map('ratio_857.fits',ratio_857)

# Plot the ratio maps
cg.plot(ratio_353_mono)
plt.savefig('dr3_dr4_353_ratio_monopole_corrected',dpi=100,bbox_inches='tight')
plt.close()
cg.plot(ratio_545_mono)
plt.savefig('dr3_dr4_545_ratio_monopole_corrected',dpi=100,bbox_inches='tight')
plt.close()
cg.plot(ratio_857_mono)
plt.savefig('dr3_dr4_857_ratio_monopole_corrected',dpi=100,bbox_inches='tight')
plt.close()
cg.plot(ratio_353)
plt.savefig('dr3_dr4_353_ratio_monopole',dpi=100,bbox_inches='tight')
plt.close()
cg.plot(ratio_545)
plt.savefig('dr3_dr4_545_ratio_monopole',dpi=100,bbox_inches='tight')
plt.close()
cg.plot(ratio_857)
plt.savefig('dr3_dr4_857_ratio_monopole',dpi=100,bbox_inches='tight')
plt.close()
