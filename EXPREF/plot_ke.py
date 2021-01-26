import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

files = ['/work/bb0519/b350090/models/release-4.0.1/tests/DW_CHANNEL_SGS/EXP00/JURICKE', \
         '/work/bb0519/b350090/models/release-4.0.1/tests/DW_CHANNEL_SGS/EXP01/JURICKE' ] 

names = ['EXP00','EXP01']

ufiles = []
vfiles = []
tfiles = []

for f in files:
   ufile = '%s_grid_U.nc' % (f,)
   vfile = '%s_grid_V.nc' % (f,)
   tfile = '%s_grid_T.nc' % (f,) 
   
   ufiles.append( xr.open_dataset(ufile) )
   vfiles.append( xr.open_dataset(vfile) )
   tfiles.append( xr.open_dataset(tfile) )

for i in range(0,len(tfiles)):
   print(tfiles[i])

def integrate_spectrum(kx,ky,ket,dx,dy):
   """
   Integrate 2D spectrum to 1D 
   For example, KE(kx,ky) to just KE(k)
   
   Input: 
     kx  - 1D array of wavenumbers in x direction
     ky  - 1D array of wavenumbers in y direction
     ket - 2D array of data in spectral space
   
   Output: 
     vk  - 1D array of wavenumbers
     ek  - 1D array of energy (power)
     psd - 1D array of energy per wavenumber (power density)
   """
   
   # 2D arrays of kx,ky 
   wx,wy = np.meshgrid(kx,ky)
   # Construct  
   wk = wx**2 + wy**2
   kmax = np.max( np.sqrt(wk) )
   nx = kx.shape[0]
   ny = ky.shape[0]
   dk = min(1./(nx*dx),1./(ny*dy))
   k = np.arange(dk, kmax + dk, dk)
   print( ' Min max kx, ky ',kx.min(),kx.max(),ky.min(),ky.max())
   print( ' Min max k, dk ',k.min(),k.max(),dk)
   # 
   # Integrate KE around lines of constant wave number  
   # PSD = d/dk E [m3/s2]   
   # E = PSD * dk [m2/s2] 
   nk = k.shape[0]
   ek = np.zeros((nk))
   for jk in range(0,nk):
      indices = np.where(wk >= k[jk]**2)
      ek[jk]  = np.sum(ket[indices])
   
   # PSD is derivate with k 
   psd = -(ek[1:] - ek[0:-1]) / dk 
   # Remove last element in ek so it matches 
   ek  = ek[0:-1]
   # Average to make array match in size 
   vk  = 0.5 * (k[1:] + k[0:-1])
   
   return vk,ek,psd


def pwke(u,v,dx,dy):
   """
   Calculate spectral KE from u,v in grid-point space
   
   Input: 
     u  - 2D array of zonal velocity 
     v  - 2D array of meridional velocity 
     dx - grid spacing in x
     dy - grid spacing in y
   
   Output:
     vk  - 1D array of wavenumbers
     ek  - 1D array of power
     psd - 1D array of power density 
   """
   
   # Wavenumbers
   nx = u.shape[-1]
   ny = u.shape[-2]
   print(' Build wavenumbers: nx, ny, dx, dy ',nx,ny,dx,dy)
   kx = np.fft.fftfreq(nx, d=dx)
   ky = np.fft.fftfreq(ny, d=dy)   
   
   # Calculate two-dim FFTs
   ut = np.fft.fft2(u)
   vt = np.fft.fft2(v)
   
   # Calculate KE
   zz = 0.5 * (np.conj(ut) * ut + np.conj(vt) * vt)
   nn = (nx**2 * ny**2) 
   # Normalise by number of grid points, squared
   ket = zz/float(nn)
   
   vk,ek,psd = integrate_spectrum(kx,ky,ket,dx,dy)
   
   return vk,ek,psd


fig1, ax1 = plt.subplots(1,2)
fig2, ax2 = plt.subplots(1,2)
fig3, ax3 = plt.subplots(1,2)

# 
for i in range(0,len(ufiles)):
   ufile = ufiles[i]
   vfile = vfiles[i]
   tfile = tfiles[i]
   name  = names[i]
   
   print( ' File: ' )
   print( ufile )
   
   # Plot SGS KE
   fig_sgs, ax_sgs = plt.subplots(1,1)
   tfile['sgske'][-1,0,:,:].plot(ax=ax_sgs)
   fig_sgs.savefig('%s_sgs_ke_last.png' % (names[i],), format='png')
   
   uoce  = ufile['uoce']
   if i==0:
      uoce0 = uoce
   voce  = vfile['voce']
   keoce = tfile['ke_zint'] / 1600.
   
   nt = len(uoce['time_counter'])
   nk = len(uoce['depthu'])
   print(' time_counter size ',nt)
   print(' depth size ',nk)
   
   all_t = []
   for jn in range(0,nt):
      all_d = []
      for jk in range(0,1):
         u2d = uoce[jn,jk,:,:].values
         v2d = voce[jn,jk,:,:].values
         times = [uoce['time_counter'][jn]]
         depth = [uoce['depthu'][jk]]
         print(' working at depth, time ',depth[0],times[0])
         vk,ek,psd = pwke(u2d,v2d,30e3,30e3) 
         
         ds = xr.Dataset( {
                            "Ek"   : (["time", "depth", "k"], ek[np.newaxis,np.newaxis,:]),
                            "PSD"  : (["time", "depth", "k"], psd[np.newaxis,np.newaxis,:]),
                          },
                          coords={
                            "k"      : (["k"], vk),
                            "depth"  : (["depth"], depth),
                            "time"   : times,
                                 },
                          )         
         #print(ds)
         all_d.append(ds)
         
      ds_d = xr.concat(all_d,dim="depth")
      #print(ds_d)
      all_t.append(ds_d)
      
   ds_t = xr.concat(all_t,dim="time")
   print(ds_t)      
   
   #fig, ax = plt.subplots(1,1)
   #ax.loglog(vk,ek)
   ds = ds_t.mean("time")
   #xr.plot.line(ds.k, ds.Ek, ax=ax)
   #ax.set_xscale("log")
   #ax.set_yscale("log")
   ax1[0].loglog(ds.k.values, ds.Ek.values[0,:], label=names[i])
   
   # Grid-point space stats
   um = uoce.mean('x').mean('time_counter') 
   km = keoce.mean('time_counter')
   print(um)
   if i == 0:
      km.plot(ax=ax2[0])
      km0 = km
      ax2[0].set_title(names[i])
   if i == 1:
      (km-km0).plot(ax=ax2[1])
      ax2[1].set_title('%s - %s' % (names[1],names[0]))
   
   for var in ['utrd_hpg','utrd_spg','utrd_udx','utrd_pvo','utrd_zad','utrd_ldf','utrd_keb']:
      if i==0:
         ufiles[0][var].isel(depthu=0).mean('time_counter').mean('x').plot(ax=ax3[0],label=var)
      if i==1:
         (ufiles[1][var] - ufiles[0][var]).isel(depthu=0).mean('time_counter').mean('x').plot(ax=ax3[1],label=var)
         
ax1[0].legend()
ax3[0].legend()

fig1.savefig('ke_spectrum.pdf',format='pdf')
fig2.savefig('ke_forcing.pdf',format='pdf')
fig3.savefig('ke_zint.png',format='png')

plt.show()
