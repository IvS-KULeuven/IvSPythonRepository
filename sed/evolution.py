# -*- coding: utf-8 -*-
import numpy as np
import pylab as pl
from scipy.interpolate import bisplrep,bisplev,Rbf
import os
from divers import conversions
from ivs.units import constants
from ivs.catalogs import vizier

dir_schaller1992 = os.path.join(os.path.dirname(__file__),'schaller1992')
dir_driebe1998 = os.path.join(os.path.dirname(__file__),'driebe1998')
dir_pamyatnykh = os.path.join(os.path.dirname(__file__),'pamyatnykh')
schaller_masses = [1,1.25,1.5,1.7,2,2.5,3,4,5,7,9,12,15,20,25,40,60]

def schaller1992():
    masses = [1,1.25,1.5,1.7,2,2.5,3,4,5,7,9,12,15,20,25,40,60][:-1]
    tables = ['table20','table18','table17','table16','table15','table14',
              'table13','table12','table11','table10','table9','table8',
              'table7','table6','table5','table4','table3'][:-1]
    all_teffs = []
    all_loggs = []
    all_radii = []
    for mass,table in zip(masses,tables):
        data,comms,units = vizier.search('J/A+AS/96/269/%s'%(table))
        all_teffs.append(10**data['logTe'])
        all_radii.append(conversions.getRadius(logL=data['logL'],Teff=10**data['logTe'])[0])
        all_loggs.append(conversions.getLogg(mass*constants.Msol_cgs,all_radii[-1]*constants.Rsol_cgs)[0])
        pl.plot(all_teffs[-1],all_loggs[-1],'k-')
    
    all_teffs = np.hstack(all_teffs)
    all_radii = np.hstack(all_radii)
    all_loggs = np.hstack(all_loggs)
    
    keep = all_teffs>5000
    all_teffs = all_teffs[keep]
    all_radii = all_radii[keep]
    all_loggs = all_loggs[keep]
    
    #mygrid = bisplrep(np.log10(all_teffs),all_loggs,all_radii,s=10000)
    #print np.sort(np.log10(all_teffs))
    
    res = 1
    mygrid = Rbf(np.log10(all_teffs)[::res],all_loggs[::res],all_radii[::res],function='linear')
    
    rteffs = np.random.uniform(low=5000,high=20000,size=1000)
    rloggs = np.random.uniform(low=2.0,high=4.5,size=1000)
    #interpl = np.array([bisplev(x,y,mygrid) for x,y in zip(np.log10(rteffs),rloggs)])
    interpl = mygrid(np.log10(rteffs),rloggs)
    
    #help(bisplev)
    
    
    pl.scatter(all_teffs,all_loggs,c=all_radii,cmap=pl.cm.spectral,vmin=0.5,vmax=25,edgecolors='none')
    pl.scatter(rteffs,rloggs,c=interpl,s=100,cmap=pl.cm.spectral,vmin=0.5,vmax=25)
    pl.colorbar()
    pl.xlim(pl.xlim()[::-1])
    pl.ylim(pl.ylim()[::-1])
    pl.gca().set_xscale('log')
    
    pl.figure()
    pl.scatter(all_teffs,all_loggs,c=all_radii,cmap=pl.cm.spectral,vmin=0.5,vmax=5,edgecolors='none')
    pl.scatter(rteffs,rloggs,c=interpl,s=100,cmap=pl.cm.spectral,vmin=0.5,vmax=5)
    pl.colorbar()
    pl.xlim(pl.xlim()[::-1])
    pl.ylim(pl.ylim()[::-1])
    pl.gca().set_xscale('log')
    
    pl.figure()
    pl.scatter(all_teffs,all_loggs,c=all_radii,cmap=pl.cm.spectral,vmin=5,vmax=15,edgecolors='none')
    pl.scatter(rteffs,rloggs,c=interpl,s=100,cmap=pl.cm.spectral,vmin=5,vmax=15)
    pl.colorbar()
    pl.xlim(pl.xlim()[::-1])
    pl.ylim(pl.ylim()[::-1])
    pl.gca().set_xscale('log')
    pl.show()

if __name__=="__main__":
    schaller1992()