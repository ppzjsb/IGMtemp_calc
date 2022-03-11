

import numpy as np
import matplotlib.pyplot as plt

# File directory and name
base = '../output/'
file = base+'puchwein2019.dat'

# Open the binary file
readdata = open(file,"rb")

# Header data
omegaM  = np.fromfile(readdata,dtype=np.double,count=1) # Omega_m (matter density)
omegaL  = np.fromfile(readdata,dtype=np.double,count=1) # Omega_L (Lambda density)
omegab  = np.fromfile(readdata,dtype=np.double,count=1) # Omega_b (baryon density)
h100    = np.fromfile(readdata,dtype=np.double,count=1) # Hubble constant, H0 / 100 km/s/Mpc
Xh      = np.fromfile(readdata,dtype=np.double,count=1) # Hydrogen fraction by mass
logdens = np.fromfile(readdata,dtype=np.double,count=1) # log10 density of gas parcel
nbins   = np.fromfile(readdata,dtype=np.int32,count=1)  # Number of redshift points

# Line of sight scale
ztime  = np.fromfile(readdata,dtype=np.double,count=nbins[0]) # redshift
ttime  = np.fromfile(readdata,dtype=np.double,count=nbins[0]) # time [s]
temp   = np.fromfile(readdata,dtype=np.double,count=nbins[0]) # temperature [K]
nH0    = np.fromfile(readdata,dtype=np.double,count=nbins[0]) # nHI/nH
nHep   = np.fromfile(readdata,dtype=np.double,count=nbins[0]) # nHeII/nH
nHepp  = np.fromfile(readdata,dtype=np.double,count=nbins[0]) # nHeIII/nH
nelec  = np.fromfile(readdata,dtype=np.double,count=nbins[0]) # ne/nH
heat   = np.fromfile(readdata,dtype=np.double,count=nbins[0]) # specific heat rate [erg s^-1 g^-1]
gJH0   = np.fromfile(readdata,dtype=np.double,count=nbins[0]) # HI photoionisation rate [s^-1]
gJHe0  = np.fromfile(readdata,dtype=np.double,count=nbins[0]) # HeI photoionisation rate [s^-1]
gJHep  = np.fromfile(readdata,dtype=np.double,count=nbins[0]) # HeII photoionisation rate [s^-1]
epsH0  = np.fromfile(readdata,dtype=np.double,count=nbins[0]) # HI photoheating rate [erg s^-1]
epsHe0 = np.fromfile(readdata,dtype=np.double,count=nbins[0]) # HeI photoheating rate [erg s^-1]
epsHep = np.fromfile(readdata,dtype=np.double,count=nbins[0]) # HeII photoheating rate [ers s^-1]




# Close the binary file
readdata.close()


yhe = (1.0-Xh)/(4.0*Xh) # nHe/nH
ev  = 1.602e-12 # erg

# Gas parcel temperature 
plt.figure(1)
plt.plot(ztime,temp/1.0e3,color='black')
plt.xlim(0,8)
plt.ylim(3,17)
plt.xlabel(r'z',fontsize=12)
plt.ylabel(r'$T/10^{3}\rm\, K$',fontsize=12)

# Fractional abundances of H0, He+ and He2+
plt.figure(2)
plt.plot(ztime,np.log10(nH0),color='blue')
plt.plot(ztime,np.log10(nHep/yhe),color='red',linestyle='dashed')
plt.plot(ztime,np.log10(nHepp/yhe),color='orange',linestyle='dotted')
plt.xlim(0,8)
plt.ylim(-7,1)
plt.xlabel(r'z',fontsize=12)
plt.ylabel(r'$n_{\rm HI}/n_{\rm H},\,n_{\rm HeII}/n_{\rm He},\,n_{\rm HeIII}/n_{\rm He}$',fontsize=12)

# Energy per photoionisation for H0, He+ and He2+
plt.figure(3)
plt.plot(ztime,np.log10(epsH0/gJH0/ev),color='blue')
plt.plot(ztime,np.log10(epsHe0/gJHe0/ev),color='red',linestyle='dashed')
plt.plot(ztime,np.log10(epsHep/gJHep/ev),color='orange',linestyle='dotted')
plt.xlim(0,8)
plt.ylim(0.5,1.8)
plt.xlabel(r'z',fontsize=12)
plt.ylabel(r'$\epsilon_{\rm i}/\Gamma_{\rm i}\, [\rm eV]$',fontsize=12)






