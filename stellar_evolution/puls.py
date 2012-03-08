# -*- coding: utf-8 -*-
"""
Interfaces to LOSC and ADIPLS pulsation codes.
"""
import os
import logging
import subprocess
import functools
import json
import itertools

import pyosclib
from ivs.io import ascii
from ivs.stellar_evolution import fileio
from ivs.stellar_evolution import adipls
from ivs.units import constants
from ivs.units import conversions
from ivs.sigproc import filtering
from ivs.aux import numpy_ext as ne
from ivs.sed import model

import numpy as np
from numpy import pi,sqrt,sin,cos
from scipy import interpolate
from scipy import integrate
from matplotlib import pyplot as pl

logger = logging.getLogger("IVS.STAR.PULS")

#-- necessary constants
np0max = 6001

class StellarModel:
    """
    Class facilitating reading in of stellar models from different codes,
    and calculating frequencies with different codes,
    
    Warning: try to first do the computations, then the plotting. It seems
    that you need to be careful with the order in which you call the OSCLIB
    routines.
    """
    #{ Built-in functions
    def __init__(self,star_model,unit='d-1',codes=('losc','adipls'),**kwargs):
        """
        Read in the stellar model.
        
        We set the model in physical units.
        
            - .mod files are OSC files
            - .dat files are CLES files 
            - .log files are MESA files
        
        Global parameters should be:
            - boundary condition
        """
        model = 'model'
        self.master = codes[0]
        self.codes = codes
        if isinstance(star_model,str):
            logger.info('Received file %s'%(star_model))
        else:
            logger.info('Received ext %s from FITS HDUList'%(kwargs['ext']))
        
        starg,starl,mtype = fileio.read_starmodel(star_model,do_standardize=True,**kwargs)
                
        self.constants = conversions.get_constants(units='cgs',values='standard')
        
        self.starg = starg
        self.starl = starl
        self.eigenfreqs = {}
        self.eigenfuncs = {}
        self.unit = unit
    
    def __len__(self):
        return len(self.starl)
    
    def __getitem__(self,key):
        return self.starl[key]   
    
    def _init_losc(self):
        """
        Convert the stellar profiles to input that LOSC understands.
        
        This effectively initiates LOSC, except for calling the meshing
        """
        #r,m,P,rho,G1 = self.starl['radius'],self.starl['mass'],self.starl['pressure'],\
        #               self.starl['Rho'],self.starl['gamma1']
        #for arr in [r,m,P,rho,G1]:
        #    if np.any(np.isnan(arr)): raise ValueError
        #    if np.any(np.isinf(arr)): raise ValueError
        #output,R,M = physical2osc(r,m,P,rho,G1,G=self.constants['GG'])
        output,R,M = star2losc(self.starg,self.starl,G=self.constants['GG'])
        #if np.isnan(output[1][0]):
        #    output[1][0] = 0
        shapenow = output.shape[1]
        if shapenow>=np0max:
            raise ValueError,'Mesh is too fine... now what?'
        else:
            logger.info('Initiated LOSC model with %d meshpoints'%(shapenow))
            #bla = self._get_current_losc_model()
            #print bla.shape
        shapeadd = np0max-shapenow
        output = np.hstack([output,np.zeros((6,shapeadd))])
        output = np.array(output,float)
        #print 'output',output.shape,shapenow
        #pyosclib.readmodpyg(np.zeros_like(output),0,6001,R,M,'model',self.constants['GG'],
                                                            #self.constants['Msol'])
        #pyosclib.mesh(0,6001,0)
        pyosclib.readmodpyg(output,0,shapenow+1,R,M,'model',self.constants['GG'],
                                                            self.constants['Msol'])
        pyosclib.mesh(0,shapenow+1,0)
        bla = self._get_current_losc_model()
        mass,radius,taudyn,npts,rc = pyosclib.modinfo()
        if rc!=0:
            raise IOError('Error in reading stellar model for estimating frequencies (%s)'%(rc))
        pyosclib.setboundary(2,rc)
        return rc
    
    def _get_current_losc_model(self):
        """
        Return model currently stored in LOSC.
        
        S1 = x
        S2 = q/x**3
        S3 = R*P_r/(G*M*rho_r)
        S4 = 4*pi*R**3*rho_r/M
        S5 = G1
        S6 = np.zeros(len(r))
        """
        star = np.zeros((6,6001))
        star = pyosclib.getmodel(star)
        keep = star[0]!=0
        if np.any(keep):
            keep[0] = True
            return star[:,keep]
        else:
            return star
    
    def _init_adipls(self):
        """
        Convert the stellar profiles to input that ADIPLS understands
        """
        names = ['x','q_x3','Vg','gamma1','brunt_A','U']
        starl = np.rec.fromarrays([self.starl[name] for name in names],
                                  names=names)
        names = ['star_mass','photosphere_r','P_c','rho_c','D5','D6','mu']
        starg = {}
        for name in names:
            starg[name] = self.starg[name]
        return starg,starl
    #}
    
    #{ Code-independent interface
    
    def compute_eigenfreqs(self,degrees,codes=None,
                           f0=(1.,'d-1'),fn=(12.,'d-1'),nscan=1000,
                           fspacing='p',boundary=0):
        """
        Compute eigenfrequencies for all available methods.
        
        Remember the values of the keyword arguments, they will be duplicated
        when computing eigenfunctions.
        """
        if codes is None:
            codes = self.codes
        self.kwargs_losc = {}
        self.kwargs_adip = {}
        self.kwargs_losc['ires'] = [None,'p','g'].index(fspacing)
        self.kwargs_losc['n'] = nscan
        self.kwargs_adip['fspacing'] = fspacing
        self.kwargs_adip['nf'] = nscan
        #self.kwargs_adip['mdintg'] = 1
        
        #-- translate boundary condition for LOSC to 1, 2 or raise ValueError
        # translate boundary condition for ADIPLS to istbc
        if boundary==0:
            pyosclib.setboundary(2,0) # vanishing pressure at surface
            self.kwargs_adip['istsbc'] = 0
        else:
            pyosclib.setboundary(1,0) # nonvanishing pressure at surface
            self.kwargs_adip['istsbc'] = 1
       
        if 'losc' in codes:
            self.compute_eigenfreqs_losc(degrees,f0=f0,fn=fn,**self.kwargs_losc)
        if 'adipls' in codes:
            self.compute_eigenfreqs_adipls(degrees,f0=f0,fn=fn,**self.kwargs_adip)
    
    def compute_eigenfuncs(self,degrees,modes=None,codes=None):
        """
        Compute eigenfunctions for all available methods.
        
        Remember the values of the keyword arguments, they will be duplicated
        when computing eigenfunctions.
        """
        if codes is None:
            codes = self.codes
        if not hasattr(self,'kwargs_losc') or not hasattr(self,'kwargs_adip'):
            raise ValueError,'please compute eigenfunctions with specialised funcs'
        for degree in degrees:
            if 'losc' in codes:
                self.compute_eigenfuncs_losc(degree,modes=modes,ires=self.kwargs_losc['ires'])
            if 'adipls' in codes:
                self.compute_eigenfuncs_adipls(degree,modes=modes,**self.kwargs_adip)
    def match(self,template='losc',tol=1.):
        """
        Match frequencies.
        
        Tolerance in percentage (unitless).
        """
        #-- now do the match for frequencies and eigenfunctions
        for code in self.eigenfreqs.keys():
            if code==template: continue
            for degree in self.eigenfreqs[code].keys():
                #-- we will match all frequencies against these ones:
                ftemplate = self.eigenfreqs[template][degree]['frequency']
                #-- match
                indices = ne.match_arrays(self.eigenfreqs[code][degree]['frequency'],ftemplate)
                mismatch = np.abs(self.eigenfreqs[code][degree]['frequency'][indices]-ftemplate)/ftemplate>=(tol/100.)
                #-- sort eigenfreqs
                self.eigenfreqs[code][degree] = self.eigenfreqs[code][degree][indices]
                self.eigenfreqs[code][degree]['frequency'][mismatch] = np.nan
                #-- sort eigenfuncs (optional)
                if code in self.eigenfuncs and degree in self.eigenfuncs[code]:
                    self.eigenfuncs[code][degree] = [self.eigenfuncs[code][degree][i] for i in indices]
    
    def get_photometry(self,photbands,distance=(10.,'pc'),**kwargs):
        """
        Retrieve colors for this stellar model.
        
        To retrieve absolute bolometric magnitude, set distance to 10 pc.
        
        @return: fluxes and bolometric magnitude at 10 pc
        """
        teff = self.starg['Teff']
        logg = conversions.convert('cm/s2','[cm/s2]',self.starg['g'])
        R = conversions.convert('cgs','Rsol',self.starg['photosphere_r'])
        distance = conversions.convert(distance[1],'Rsol',distance[0])
        scale = (R/distance)**2
        logger.info('Retrieving photometry for model teff=%.0fK logg=%.3f, R=%.3fRsol in bands %s'%(teff,logg,R,', '.join(photbands)))
        try:
            fluxes,Lbol = model.get_itable(teff=teff,logg=logg,photbands=photbands,**kwargs)
            Lbol = Lbol*scale
        except IOError:
            wave,flux = model.get_table(teff=teff,logg=logg,**kwargs)
            fluxes = model.synthetic_flux(wave,flux,photbands=photbands)
            Lbol = model.luminosity(wave,flux*scale)
        Lbol = -2.5*np.log10(Lbol)-38.406612779
        return fluxes,Lbol   
    #}
    
    #{ Resampling
    
    def resample(self,increase=2,method='simple',mode_type=None,n_new=5000,int_method='spline'):
        """
        Resample the current mesh.
        
        @parameter increase: increase the resolution of the mesh with this factor.
        """
        r = self.starl['radius']
        myR = self.starl['radius'].max()
        myC = self.starl['radius'].min()
        R = self.starg['photosphere_r']
        
        #-- check if our input model is OK
        for name in self.starl.dtype.names:
            weird = np.isnan(self.starl[name]) | np.isinf(self.starl[name])
            if np.any(weird):
                logger.warning('NaN or inf encountered in %s (before resampling), setting to zero'%(name))            
        
        n_init = len(r)
        #-- my own simple method
        if method.lower()=='simple':
            newr = r.copy()
            while increase>1:
                increase-=1
                newr = sorted(np.hstack([newr[:-1],newr[:-1]+np.diff(newr)/2,[newr[-1]]]))
            new_starl = np.zeros(len(newr),dtype=self.starl.dtype)
            for name in new_starl.dtype.names:
                if name=='radius':
                    new_starl['radius'] = newr
                    continue
                if int_method=='linear':
                    f = interpolate.interp1d(r,self.starl[name])
                    new_starl[name] = f(newr)
                elif int_method=='spline':
                    #new_starl[name] = interpolate.barycentric_interpolate(r,self.starl[name],newr)
                    ius = interpolate.InterpolatedUnivariateSpline(r,self.starl[name])
                    new_starl[name] = ius(newr)
                    new_starl[name][0] = self.starl[name][0]
        
            
        #-- ADIPLS method
        elif method.lower()=='adipls':
            #   initialize ADIPLS
            adig,adil = self._init_adipls()        
            #   write adipls variables to a GONG file, and convert it to the
            #   binary ADIPLS format
            model = 'tempmodel_int'
            fileio.write_gong(adig,adil,model+'.gong')
            subprocess.call('form-amdl.d 2 %s %s'%(model+'.gong',model+'.amdl'),stdout=open('/dev/null','w'),shell=True)
            #   make a control file for the the remeshing tool in ADIPLS
            control_file = adipls.make_redistrb_inputfile(model_name=model,nn=n_new,mode_type=mode_type)
            #   call the remeshing tool, read the results back in and remove
            #   temporary files
            subprocess.call('redistrb.c.d %s'%(control_file),stdout=open('/dev/null','r'),shell=True)
            adig_,adil_ = fileio.read_amdl(model+'.amdlr')
            #os.unlink(model+'.amdlr')
            #os.unlink(model+'.amdl')
            #os.unlink(model+'.gong')
            #os.unlink(control_file)
            #   now reinterpolate all the quantities onto this mesh. Remark
            #   that we loose the interpolation made by ADIPLS and only retain
            #   its mesh structure.
            newr = adil_['x']*self.starg['photosphere_r']
            new_starl = np.zeros(len(newr),dtype=self.starl.dtype)
            for name in new_starl.dtype.names:
                if name=='radius':
                    new_starl['radius'] = newr
                    continue
                if name=='x':
                    new_starl['x'] = newr/R
                    continue
                if int_method=='linear':
                    new_starl[name] = np.interp(newr,r,self.starl[name])
                elif int_method=='spline':
                    #new_starl[name] = interpolate.barycentric_interpolate(r,self.starl[name],newr)
                    ius = interpolate.InterpolatedUnivariateSpline(r,self.starl[name])
                    new_starl[name] = ius(newr)
                    new_starl[name][0] = self.starl[name][0]
                
                
            cols = new_starl.dtype.names
            for name in adil_.dtype.names:
                if not name in cols:
                    new_starl = pl.mlab.rec_append_fields(new_starl,[name],[adil_[name]])            
            
        #-- LOSC method
        elif method.lower()=='losc':
            self._init_losc()
            ires = 0
            if mode_type is not None:
                if 'p' in mode_type.lower(): ires += 1
                if 'g' in mode_type.lower(): ires += 2
            nres = n_new#max(len(self.starl),n_new)
            pyosclib.mesh(ires,nres,0)
            newr = self._get_current_losc_model()[0]*self.starg['photosphere_r']
            new_starl = np.zeros(len(newr),dtype=self.starl.dtype)
            #-- interpolate *all* quantities on the grid. This might introduce
            #   some artifacts!
            for name in new_starl.dtype.names:
                if name=='radius':
                    new_starl[name] = newr
                    continue
                if name=='x':
                    new_starl['x'] = newr/R
                    continue
                if int_method=='spline':
                    ius = interpolate.InterpolatedUnivariateSpline(r,self.starl[name])
                    new_starl[name] = ius(newr)
                    new_starl[name][0] = self.starl[name][0]
                else:
                    f = interpolate.interp1d(r,self.starl[name],kind=int_method)
                    new_starl[name] = f(newr)
                    
                
            if len(new_starl)>=np0max:
                raise ValueError,'LOSC: too many points (%d>=%d)'%(len(new_starl),np0max)
        else:
            raise ValueError,'cannot interpret interpolation method %s'%(method)
        
        
        logger.info('Interpolated (method=%s) from %d --> %d meshpoints'%(method,n_init,len(new_starl)))
        for name in new_starl.dtype.names:
            weird = np.isnan(new_starl[name]) | np.isinf(new_starl[name])
            if np.any(weird):
                logger.warning('NaN or inf encountered in %s, setting to zero'%(name))            
                new_starl[name][weird] = 0
        
        #-- check for double entries and check if the radius
        #   hasn't changed
        if np.any(np.diff(newr)<=0):
            index = np.arange(len(newr)-1)[np.diff(newr)<=0]
            print 'Difference:',np.diff(newr)[np.diff(newr)<=0]
            print 'Indices:',index
            print 'Region:',newr[index-1:index+2]
            pl.figure()
            pl.plot(newr,'o-')
            pl.show()
            raise ValueError,'interpolation non-increasing'
        if not myC==new_starl['radius'].min():
            raise ValueError,'new center mesh point'
        #if not myR==new_starl['radius'].max():
        #    raise ValueError,'new boundary %s-->%s'%(myR/constants.Rsol_cgs,
        #                              new_starl['radius'].max()/constants.Rsol_cgs)
        self.starl = new_starl
        
        #-- if the method was not LOSC, we have to initiate it anyway
        if method.lower()!='losc':
            #-- we need to initiate the mesh in LOSC one way or
            #   the other.
            self._init_losc()
            pyosclib.mesh(0,len(self.starl),0)
    
    
    def resample_adipls(self,n_new=5000.,mode_type=None,cleanup=True):
        #   initialize ADIPLS
        adig,adil = self._init_adipls()        
        #   write adipls variables to a GONG file, and convert it to the
        #   binary ADIPLS format
        model = 'tempmodel_int'
        fileio.write_gong(adig,adil,model+'.gong')
        subprocess.call('form-amdl.d 2 %s %s'%(model+'.gong',model+'.amdl'),stdout=open('/dev/null','w'),shell=True)
        #   make a control file for the the remeshing tool in ADIPLS
        control_file = adipls.make_redistrb_inputfile(model_name=model,nn=n_new,mode_type=mode_type)
        #   call the remeshing tool, read the results back in and remove
        #   temporary files
        subprocess.call('redistrb.c.d %s'%(control_file),stdout=open('/dev/null','r'),shell=True)
        adig_,adil_ = fileio.read_amdl(model+'.amdlr')
        #-- file cleanup
        if cleanup:
            os.unlink(model+'.amdlr')
            os.unlink(model+'.amdl')
            os.unlink(model+'.gong')
            os.unlink(control_file)
        return adig_,adil_
    
    #}
    #{ Plotting routines
    
    def plot_model(self):
        r,m,P,rho,G1,dPdr = self.starl['radius'],self.starl['mass'],self.starl['pressure'],\
                       self.starl['Rho'],self.starl['gamma1'],self.starl['dP_dr']
        pl.subplot(231);pl.title('Mesh density')
        pl.plot(r)
        ax1 = pl.subplot(232);pl.title('mass')
        pl.plot(r,m)
        pl.subplot(233,sharex=ax1);pl.title('pressure')
        pl.plot(r,P)
        pl.subplot(234,sharex=ax1);pl.title('Rho')
        pl.plot(r,rho)
        pl.subplot(235,sharex=ax1);pl.title('gamma1')
        pl.plot(r,G1)
        pl.subplot(236,sharex=ax1);pl.title('dP_dr')
        pl.plot(r,dPdr)
    
    def plot_pulsation_model(self,x='x',color_cycle=None):
        star = self._get_current_losc_model()
        starg,starl = self.starg,self.starl
        #['Vg','gamma1','A']
        
        if color_cycle is None:
            color_cycle = itertools.cycle(['k','r'])
        
        c1 = color_cycle.next()
        c2 = color_cycle.next()
        
        
        #pl.figure()
        #pl.plot(self.starl[x],(star[1]-star[1][0])/self.starl[x]**2,'k-')
        #pl.show()
        
        
        
            
        if len(self.codes)>1 and not np.all(star[0]==0):
            if not star.shape[1]==starl.shape[0]:
                pl.figure()
                pl.plot(star[0],star[1],'+-')
                pl.plot(starl['x'],starl['q_x3'],'x-')
                pl.show()
                raise ValueError,'sanity check failed (LOSC (%d) and ADIPLS (%d) model have unequal number of mesh points)'%(star.shape[1],starl.shape[0])
            if not np.allclose(star[0],starl['x']):
                #wrong = np.abs(star[0]-starl['x'])>=1e-3
                print zip(star[0],starl['x'])
                raise ValueError,'sanity check failed (difference in LOSC and ADIPLS model)'
                
        elif not 'losc' in self.codes:
            star = np.zeros((6,len(starl)))
        
        ax1 = pl.subplot(231)
        pl.plot(star[0],star[1],'-',lw=2,label='LOSC',color=c1)
        pl.plot(self.starl[x],starl['q_x3'],'--',lw=2,label='ADIPLS',color=c2)
        if not ax1.get_ylabel()=='q_x3':
            leg = pl.legend(fancybox=True)
            leg.get_frame().set_alpha(0.5)
        pl.ylabel(r'q/x$^3$')
        pl.xlabel('r/R')
        
        pl.subplot(232,sharex=ax1)
        pl.ylabel(r'RP/(GM$\rho$)')
        pl.plot(star[0],star[2],'-',lw=2,color=c1)
        pl.twinx(pl.gca())
        pl.ylabel('Vg')
        pl.xlabel('r/R')
        pl.plot(self.starl[x],starl['Vg'],'--',lw=2,color=c2)
        
        pl.subplot(233,sharex=ax1)
        pl.ylabel(r'4$\pi$r$^3\rho$/m')
        pl.xlabel('r/R')
        pl.plot(star[0],star[3]/star[1],'-',lw=2,color=c1)
        pl.plot(self.starl[x],starl['U'],'--',lw=2,color=c2)
        
        pl.subplot(234,sharex=ax1)
        pl.ylabel('gamma1')
        pl.xlabel('r/R')
        pl.plot(star[0],star[4],'-',lw=2,color=c1)
        pl.plot(self.starl[x],starl['gamma1'],'--',lw=2,color=c2)
        
        pl.subplot(235,sharex=ax1)
        pl.ylabel('Brunt A')
        pl.xlabel('r/R')
        pl.plot(star[0],-star[5]*star[0]**2,'-',lw=2,color=c1)
        pl.plot(self.starl[x],starl['brunt_A'],'--',lw=2,color=c2)
        
        pl.subplot(236,sharex=ax1)
        pl.ylabel('Mesh point')
        pl.xlabel('r/R')
        pl.plot(star[0],range(len(star[0])),'-',lw=2,color=c1)
        pl.plot(self.starl[x],range(len(starl['x'])),'--',lw=2,color=c2)    
    
    def plot_frequency_matching(self,master='losc'):
    
        rows,cols = 2,2
        
        color_cycle = itertools.cycle([pl.cm.spectral(j) for j in np.linspace(0,1,len(self.eigenfuncs[master]['l0']))])
        color_cycle_degrees = itertools.cycle([pl.cm.spectral(j) for j in np.linspace(0,1,len(self.eigenfreqs[master]))])
        
        #-- make some plots
        fig = pl.figure(figsize=(16,9))
        fig.canvas.set_window_title('Frequency matching')
        
        ax1 = pl.subplot(rows,cols,1)
        ax2 = pl.subplot(rows,cols,3,sharex=ax1)
        ax3 = pl.subplot(rows,cols,2)
        ax4 = pl.subplot(rows,cols,4,sharex=ax3)
        #-- frequency differences between matched frequencies and separate
        #   frequencies
        #-- run over all degrees in the master
        comp_keys = [key for key in sorted(self.eigenfreqs.keys()) if not key==master]
        for d,degree in enumerate(sorted(self.eigenfreqs[master].keys())):
            color = color_cycle_degrees.next()
            freq_master = self.eigenfreqs[master][degree]['frequency']
            
            
                    #x2,ar2,br2 = self.eigenfuncs['adipls']['l%d'%(l)][i]
                #except:
                    #continue
                #pl.plot(x1,np.abs(ar1),'-',color=color)
            
            #-- compare them with all the degrees of the comparison codes
            for i,comp in enumerate(comp_keys):
                if not degree in self.eigenfreqs[comp]: continue
                freq_comp = self.eigenfreqs[comp][degree]['frequency']
                
                pl.axes(ax1)
                pl.plot(freq_master,(freq_master-freq_comp)/freq_master*100.,'o-',color=color,label=degree)
                
                pl.axes(ax2)
                for fr in freq_comp:
                    pl.plot([fr,fr],[i+1,i+2],'-',color=color)
                
                pl.axes(ax3)
                pl.plot(1/freq_master,(1/freq_master-1/freq_comp)*freq_master*100.,'o-',color=color,label=degree)
                
                pl.axes(ax4)
                for fr in freq_comp:
                    pl.plot([1/fr,1/fr],[i+1,i+2],'-',color=color)
                
                
            
            pl.axes(ax2)
            for fr in freq_master:
                pl.plot([fr,fr],[0,1],'-',color=color)
            pl.axes(ax4)
            for fr in freq_master:
                pl.plot([1/fr,1/fr],[0,1],'-',color=color)
            
            ##-- plot eigenfunctions
            #ax3 = pl.subplot(rows,cols,3+d)
            #for i in range(len(freq_master)):
                #color_n = color_cycle.next()
                ##-- plot master eigenfunctions
                #if self.eigenfuncs[master][degree][i] is None:
                    #continue
                #x1,ar1,br1 = self.eigenfuncs[master][degree][i]['x'],\
                             #self.eigenfuncs[master][degree][i]['ar'],\
                             #self.eigenfuncs[master][degree][i]['br']
                #p, =pl.plot(x1,np.abs(ar1),'x-',color=color_n)
                #norm = np.abs(ar1).max()
                ##-- plot the other eigenfunctions, but only if they have been
                ##   matched
                #matched = False
                #if not comp_keys: #-- if no other comparison, plot them anyway
                    #matched = True
                #for comp in comp_keys:
                    #if not degree in self.eigenfuncs[comp]: continue
                    #out = self.eigenfuncs[comp][degree][i]
                    #if out is not None:
                        #matched = True
                        #x2,ar2,br2 = out['x'],out['ar'],out['br']
                        #pl.plot(x2,np.abs(ar2)/np.abs(ar2).max()*norm,'+--',color=color_n)
                ##-- don't clutter the plot with eigenfunction that are not
                ##   matched
                #if not matched:
                    #p.remove()
            
            #pl.axes(ax3)
            #pl.xlabel('Position in star [r/R]')
            #pl.ylabel('Eigenfunction amplitude [arbitrary units]')
            #pl.gca().set_yscale('log')
                    
        
        pl.axes(ax1)
        pl.xlabel('Frequency [muHz]')
        pl.ylabel('$f_\mathrm{losc}-f_\mathrm{adipls}/f_\mathrm{losc}$ [%]')
        pl.grid()
        pl.legend(loc='best',numpoints=1)
        
        pl.axes(ax2)
        pl.xlabel('Frequency [muHz]')
        pl.legend(loc='best')
        
        pl.axes(ax3)
        pl.xlabel('Period [1/muHz]')
        pl.ylabel('$P_\mathrm{losc}-P_\mathrm{adipls}/P_\mathrm{losc}$ [%]')
        pl.grid()
        pl.legend(loc='best',numpoints=1)
        
        pl.axes(ax4)
        pl.xlabel('Period [1/muHz]')
        pl.legend(loc='best')
    
    def plot_weight_functions(self,prefix=''):
        #-- compute weight functions
        for degree in self.eigenfuncs['losc']:
            for i,mode in enumerate(self.eigenfuncs['losc'][degree]):
            
                pl.figure(figsize=(16,9))
                order = self.eigenfreqs['losc'][degree]['modeLee'][i]
                sigma = self.eigenfreqs['losc'][degree]['sigma'][i]
                print degree,i,order
                x,a,b = mode['x'],mode['ar'],mode['br']
                Pp,T,C,N,G = mode['Pprime'],mode['T'],mode['C'],mode['N'],mode['G']
                total,charpinet = mode['weight'],mode['charpinet']
                Ccontr = np.trapz(C,x=x)
                Ncontr = np.trapz(N,x=x)
                Gcontr = np.trapz(G,x=x)
                tcontr = Ccontr+Ncontr+Gcontr
                varsig = np.sqrt(integrate.simps(total,x=self['radius'])/integrate.simps(T,x=self['radius']))
                dsig2 = (sigma**2-varsig**2)/(sigma**2)*100
                dsig = (sigma-varsig)/sigma*100
                
                pl.axes([0.1,0.3,0.4,0.65])
                pl.title(r'$\ell=%s, n=%d$'%(degree[1:],order))
                pl.plot(x,T/T.max(),'r-',label='T (Kinetic)')
                norm = total.max()
                pl.plot(x,C/norm,'c-',label='C (Acoustic %.3g%%)'%(Ccontr/tcontr*100))
                pl.plot(x,N/norm,'b-',label='N (Buoyancy %.3g%%)'%(Ncontr/tcontr*100))
                pl.plot(x,G/norm,'g-',label='G (Gravity %.3g%%)'%(Gcontr/tcontr*100))
                pl.plot(x,total/norm,'k-',label='Total weight')
                pl.plot(x,charpinet/norm,'m--',label='Charpinet')
                leg = pl.legend(loc='best',fancybox=True)
                leg.get_frame().set_alpha(0.5)
                pl.grid()
                pl.ylabel('Normalised weight function')
                
                pl.axes([0.1,0.1,0.4,0.17])
                pl.plot(x,abs(a)/abs(a).max(),'-',lw=2,color='0.5')
                pl.gca().set_yscale('log')
                pl.grid()
                pl.xlabel('Radial position [r/R]')
                
                pl.axes([0.55,0.1,0.4,0.85])
                perc = integrate.cumtrapz(total,x=x)/integrate.trapz(total,x=x)*100
                ninenine_perc = x[np.argmin(np.abs(perc-99.))]
                ninenine_perc_index = np.argmin(np.abs(ninenine_perc-self['x']))
                pl.plot(x[:-1],perc,'k-')
                pl.plot(ninenine_perc,99,'ro')
                pl.annotate('x(99%%)=%.3g'%(ninenine_perc),(0.5,0.95),xycoords='axes fraction')
                pl.annotate('q(99%%)=%.3g'%(self['q'][ninenine_perc_index]),(0.5,0.9),xycoords='axes fraction')
                pl.xlabel('Radial position [r/R]')
                pl.ylabel('Cumulative contribution of $\sigma^2$ (%)')
                pl.ylim(0,100)
                pl.xlim(0,1)
                pl.title('$\Delta\sigma^2$=%.3g%%, $\Delta\sigma$=%.3g%%'%(dsig2,dsig))
                pl.twinx(pl.gca())
                N = conversions.convert('rad/s',self.unit,np.sqrt(self['brunt_N2']))
                #pl.figure()
                #pl.plot(self['radius'],N,'+-')
                #pl.show()
                #raise SystemExit
                degree_ = int(degree[1:])
                Sl = conversions.convert('rad/s',self.unit,np.sqrt(self['csound']**2/self['radius']**2*(degree_*(degree_+1))))
                pl.plot(self['x'],N,'g-',label=r'$N$')
                pl.plot(self['x'],Sl,'r-',label=r'$S_{\ell=%s}$'%(degree[1:]))
                for fi,freq in enumerate(self.eigenfreqs['losc'][degree]['frequency']):
                    pl.plot([0,1],[freq,freq],'-',lw=2,color='0.5',alpha=0.5)
                    if i==fi: pl.plot([0,1],[freq,freq],'k--',lw=2)
                pl.gca().set_yscale('log')
                pl.xlim(0,1)
                pl.ylim(1,5000)
                N[np.isnan(N)] = 1e-10
                if int(degree_)>0:
                    maxNSl = np.where(Sl>=N,Sl,N)
                    maxNSl[np.isnan(maxNSl)] = pl.ylim()[1]
                    minNSl = np.where(Sl<=N,Sl,N)
                    minNSl[np.isnan(minNSl)] = 1e-10
                    pl.fill_between(self['x'],maxNSl,pl.ylim()[1],alpha=0.2,color='m')
                else:
                    minNSl = N
                pl.fill_between(self['x'],minNSl,np.ones_like(minNSl)*pl.ylim()[0],alpha=0.2,color='y',linewidth=0)
                pl.grid()
                pl.ylabel(r"Frequency (%s)"%(self.unit))
                leg = pl.legend(loc='best',fancybox=True)
                leg.get_frame().set_alpha(0.5)
                pl.savefig(prefix+'weight_function_%sn%+03d'%(degree,order))
                pl.close()
  
    def plot_density_temperature_profile(self):
        return None
    #}
    
    
    
    #{ Code-dependent interface
    def compute_eigenfreqs_losc(self,degrees,f0=(1.,'d-1'),fn=(12.,'d-1'),
                                     ires=0,n=500,ind=1,mmax=50,
                                     typecomp=1):
        """
        Compute eigenfrequencies with the LOSC code.
        
        @parameter f0: approximate start frequency
        @parameter fn: approximate end frequency
        @parameter ires: type of mesh for computation of the eigenfrequencies,
        1 to optimize for p modes, 2 for g modes and 3 for p and g modes
        @parameter n: number of subintervals in scanomaga
        @parameter ind: type of grid: 1 is equidistant in frequency space (for
        p modes), 2 is equidistant in period space (for g modes)
        @parameter mmax: maximum number of approximate eigenfrequencies to compute.
        @parameter typecomp: type of computations, 1 for normal computation, which is
        recommended, and type=2 if convergence is difficult to obtain.
        @parameter unit: units to put in the 'frequency' output column
        """
        old_settings = np.seterr(invalid='ignore')
        modekeys = ['l','nz','mode','modeLee','parity','sigma','beta','ev','xm',
                'delta','boundary','frequency','freqinit']
        dtypes = [int,int,int,int,int,float,float,float,float,float,float,float,float]
        #-- convert frequency to dimensionless angular frequencies for OSCL
        f0 = conversions.convert(f0[1],'rad/s',f0[0])
        fn = conversions.convert(fn[1],'rad/s',fn[0])
        R = self.starg['photosphere_r']
        M = self.starg['star_mass']
        t_dynamic = np.sqrt(R**3/(self.constants['GG']*M))
        f0 = f0*t_dynamic
        fn = fn*t_dynamic
        
        rc = 0
        if not 'losc' in self.eigenfreqs:
            self.eigenfreqs['losc'] = {}
        if not 'losc' in self.eigenfuncs:
            self.eigenfuncs['losc'] = {}
        #-- run over all the degrees
        for l in degrees:
            #-- estimate the frequencies
            pyosclib.setdegree(l)
            omegas = np.zeros(100)
            mm = len(omegas)
            omegas,mm = pyosclib.scanomega(f0,fn,n,ind,mmax,omegas,mm,rc)
            if rc!=0:
                raise ValueError('Error in scanning for frequencies')
            omegas = omegas[:mm]
            if 4*mm>n:
                raise ValueError('Frequency spectrum is dense, increase "n"')
            #logger.info("Estimated %d frequencies (l=%d) between %.3fd-1 and %.3fd-1"%(mm,l,omegas.min(),omegas.max()))
            #-- refine the frequencies and compute mode information
            freqinfo = np.zeros((len(omegas),len(modekeys)))
            for i,omega in enumerate(omegas):
                pyosclib.computemode(typecomp,omega,rc)
                freqinfo[i][:11] = pyosclib.oscinfo(omega,rc)
                freqinfo[i][12]  = omega
            #-- make the array into a convenient record array
            freqinfo = [np.array(i,dtype) for i,dtype in zip(freqinfo.T,dtypes)]
            freqinfo = np.rec.fromarrays(freqinfo,names=modekeys)
            #-- put the frequencies in cy/d in the array
            freqcd = conversions.convert('rad/s',self.unit,freqinfo['sigma'])
            freqinfo['frequency'] = freqcd
            self.eigenfreqs['losc']['l%d'%(l)] = freqinfo
            #-- make room for eigenfunctions
            self.eigenfuncs['losc']['l%d'%(l)] = [None for i in range(len(self.eigenfreqs['losc']['l%d'%(l)]))]
            
            logger.info('LOSC: found %d frequencies (l=%d) between %.3g%s and %.3g%s'%(len(freqcd),l,freqcd.min(),self.unit,freqcd.max(),self.unit))
        np.seterr(**old_settings)
    
    def compute_eigenfuncs_losc(self,degree,modes=None,ires=0,typecomp=1):
        """
        Compute eigenfunctions with the LOSC code.
        
        Also compute horizontal and vertical displacements. In the book of
        Aerts et al., these are known as $\tilde{\ksi}_r$ and $\tilde{\ksi}_h$.
            
        WARNING: the OSC manual notes that this is valid 'near the centre', but is
        actually a pretty fuzzy sentence...
        
        @parameter ires: type of mesh for computation of the eigenfrequencies,
        1 to optimize for p modes, 2 for g modes and 3 for p and g modes
        @parameter nres: the number of added points increases with n. This parameter
        should be increased as the complexity of the model increases
        @parameter typecomp: type of computations, 1 for normal computation, which is
        recommended, and type=2 if convergence is difficult to obtain.
        """
        old_settings = np.seterr(invalid='ignore',divide='ignore')
        if modes is None:
            modes = range(len(self.eigenfreqs['losc']['l%d'%(degree)]))
        rc = 0
        R = self.starg['photosphere_r']
        M = self.starg['star_mass']
        #-- below, we will need the dynamical time scale, the radial position
        #   in the star, the mass position and the density
        t_dynamic = np.sqrt(R**3/(self.constants['GG']*M))
        x_ = self.starl['radius']/R
        q_ = self.starl['mass']/M
        rho_ = self.starl['Rho']
        P_ = self.starl['pressure']
        G1_= self.starl['gamma1']
        dP_dr_ = self.starl['dP_dr']
        l = degree
        
        #-- run over all the degrees
        for mode in modes:
            #-- compute the eigenfrequency
            mymode = self.eigenfreqs['losc']['l%d'%(l)][mode]
            pyosclib.setdegree(l)
            omega = mymode['freqinit']
            if np.isnan(omega):
                logger.info('Mode %d is not matched (nan)'%(mode))
                continue
            pyosclib.computemode(typecomp,omega,rc)
            #-- maybe we should cal oscinfo here and check if results match
            #   the previous calculations, just as a sanity test
            info = np.array(pyosclib.oscinfo(omega,rc))
            if not np.allclose(list(mymode)[:11],info):
                print zip(mymode.dtype.names[:11],list(mymode)[:11],info)
                raise ValueError,'sanity check in eigenfunc computation failed' 
            #-- retrieve the eigenfunctions and the stellar (dimensionless)
            #   quantities from LOSC and remove zeros
            star,oscl = pyosclib.osc08py()
            clip = [oscl[0]!=0]
            star = [col[clip] for col in star]
            oscl = [col[clip] for col in oscl]
            Y,Z,U,V = oscl
            #-- get the physical quantities
            x = star[0]                              # normalised radius
            r = x*R                                  # radius
            rho = star[3]/(4*np.pi*R**3)*M           # density
            q = star[1]*x**3                         # normalised mass
            m = q*M                                  # mass
            P = star[2]/R*self.constants['GG']*M*rho # pressure
            g = self.constants['GG']*m/r**2          # gravity
            G1 = star[4]                             # Adiabatic exponent
            N2 = -g*star[5]*x/R                      # Brunt-Vaisalla frequency
            dP_dr = -rho*g                           # Pressure gradient using hydrostatic equilibrium
            U_ = r/m*4*pi*r**2*rho                   # r/m*dm/dr = dlnm/dlnr and dm/dr = 4*pi*r^2*rho using conservation of mass
            #-- compute the horizontal and vertical components
            #   first, we need the angular dimensionless frequency
            omega = mymode['sigma']*t_dynamic
            sigma = mymode['sigma']
            if l==0:#   in the case of a radial mode
                a_r = R*x*Y
                b_r = np.zeros(len(a_r))
            else:#   in the case of a nonradial mode
                a_r = R*x**(l-1)*Y
                b_r = R*x**(l-1)/(omega**2) * (U + R*P/(self.constants['GG']*M*rho)*Z + q/x**3*Y)
            #-- now compute the lagrangian pressure perturbation
            if l==0:#   in the case of a radial mode:
                dP_P = Z
            else:#   in the case of nonradial mode:
                dP_P = x**l*Z
            #  and convert it to Eulerian
            Pprime = dP_P*P - a_r*dP_dr
            #-- compute the Eulerian gravity potential perturbation and its
            #   derivative
            Phiprime = self.constants['GG']*M/R*x**l*U
            dPhiprime_dr = self.constants['GG']*M/R**2*x**(l-1)*(V-4*np.pi*R**3*rho/M*Y)
            #-- compute density perturbation ()
            #Rhoprime = see page 52 of JCD lecture notes
            #-- compute temperature perturbation ()
            #Tprime = see page 208 of JCD lecture notes
            #-- Dziembowski variables
            dr = a_r
            #dr = np.sqrt(a_r**2 + b_r**2)
            y1 = dr/r
            y2 = 1./(g*r) * (Pprime/rho + Phiprime)
            y3 = 1./(g*r) * Phiprime
            y4 = 1./g * dPhiprime_dr
            #-- weight functions
            Sl2 = G1*P/(r**2*rho)   # remove the l*(l+1) thing, it is cancelled in C anyway
            #-- compute the weight functions
            T = dr**2 + l*(l+1)/(r**2*sigma**2)*(Pprime/rho+Phiprime) # proportional to kinetic energy density
            C = g**2/Sl2*(y2-y3)**2                              # acoustic energy
            N = N2*dr**2                                        # buoyancy energy
            G = -1/(4*np.pi*rho*self.constants['GG'])*(dPhiprime_dr + l*(l+1)*Phiprime/r) # gravitational energy
            T[0] = 0
            C[0] = 0
            N[0] = 0
            G[0] = 0
            total = (C+N+G)*rho*r**2
            charpinet = (a_r**2*N2 + Pprime**2/(G1*P*rho) + (Pprime/(G1*P) + a_r*N2/g))*rho*r**2
            varsig = np.sqrt(integrate.simps(total,x=r)/integrate.simps((T*rho*r**2),x=r))
            logger.info("Difference between freq2 and varfreq2: %.3g%%"%(np.abs(varsig**2-sigma**2)/sigma**2*100))
            
            #-- we check here explicitly Masao's (Takata) dipolar condition
            if l==1:
                zero = Pprime + g/(4*np.pi*self.constants['GG'])*(dPhiprime_dr+2/r*Phiprime) - omega**2*r*(rho*dr + 1/(4*np.pi*self.constants['GG'])*(dPhiprime_dr-Phiprime/r))
                logger.info("Masao's check: %.3g"%(np.sqrt((zero**2).sum())))
            
            self.eigenfuncs['losc']['l%d'%(l)][mode] = np.rec.fromarrays([x,a_r,b_r,Pprime,Phiprime,dPhiprime_dr,
                 T*rho*r**2,C*rho*r**2,N*rho*r**2,G*rho*r**2,total,charpinet],
                 names=['x','ar','br','Pprime','Phiprime','dPhiprime_dr','T','C','N','G','weight','charpinet'])
        np.seterr(**old_settings)
        logger.info('LOSC: computed eigenfunctions for l=%d modes %s'%(degree,modes))
                
    
    def compute_eigenfreqs_adipls(self,degrees,f0=(1.,'d-1'),fn=(12.,'d-1'),
                                  nf=1000,verbose=False,**kwargs):
        """
        Compute eigenfrequencies with ADIPLS.
        
        Possible keyword arguments are described in ADIPLS documentation.
        """
        #-- in ADIPLS, we work with files, so we need to keep track of the names
        #   of these files
        model = 'tempmodel_frq'
        #-- initialize ADIPLS
        if 'n_new' in kwargs:
            adig,adil = resample_adipls(n_new=kwargs['n_new'],mode_type=kwargs.get('mode_type',None))
        else:
            adig,adil = self._init_adipls()        
        kwargs.setdefault('cgrav',self.constants['GG'])
        kwargs.setdefault('fspacing','p')
        #-- write adipls variables to file
        #ascii.write_array(adil,'tempmodel_frq.ascii',header=True,auto_width=True,use_float='%.15e',comments=['#'+json.dumps(adig)])
        #special.write_amdl(adig,adil,model+'.amdl')
        fileio.write_gong(adig,adil,model+'.gong')
        subprocess.call('form-amdl.d 2 %s %s'%(model+'.gong',model+'.amdl'),stdout=open('/dev/null','w'),shell=True)
        
        
        if not 'adipls' in self.eigenfreqs:
            self.eigenfreqs['adipls'] = {}
        if not 'adipls' in self.eigenfuncs:
            self.eigenfuncs['adipls'] = {}
        
        #-- run over all the degrees
        for degree in degrees:
            ckwargs = kwargs.copy()
            #-- use Takata ordering scheme for dipolar modes
            if degree==1:
                ckwargs.setdefault('irsord',20)
            #-- use Lee ordering scheme for other modes
            else:
                ckwargs.setdefault('irsord',11)
            #-- construct adipls control file
            control_file = adipls.make_adipls_inputfile(model_name=model,degree=degree,
                                               f0=f0,fn=fn,**ckwargs)
            #-- call adipls, perhaps with throwing away the output
            if not verbose:
                stdout = open('/dev/null', 'w')
            else:
                stdout = None
            subprocess.call('adipls.c.d %s'%(control_file),stdout=stdout,shell=True)
            #-- read small summary file and add contents to output
            starg,starl = fileio.read_assm(model+'.ssm')
            #os.unlink(model+'.ssm')
            #os.unlink(model+'.gsm')
            #-- convert frequencies into whatever unit
            starl['frequency'] = conversions.convert('mHz',self.unit,starl['frequency'])
            self.eigenfreqs['adipls']['l%d'%(degree)] = starl
            #-- make room for eigenfunctions
            self.eigenfuncs['adipls']['l%d'%(degree)] = [None for i in range(len(self.eigenfreqs['adipls']['l%d'%(degree)]))]
            logger.info('ADIPLS: found %d frequencies (l=%d) between %.3g%s and %.3g%s'%(len(starl),degree,starl['frequency'].min(),self.unit,starl['frequency'].max(),self.unit))
        #os.unlink(model+'.amdl')
        #os.unlink(model+'.gong')
        #mymod = self._get_current_losc_model()
        #ascii.write_array(mymod,model+'.tosc')
        #raise SystemExit
        
    def compute_eigenfuncs_adipls(self,degree,modes=None,verbose=False,**kwargs):
        """
        Compute eigenfunctions with adipls.
        """
        if modes is None:
            modes = range(len(self.eigenfreqs['adipls']['l%d'%(degree)]))
        #-- in ADIPLS, we work with files, so we need to keep track of the names
        #   of these files
        model = 'tempmodel_eig'
        #-- initialize ADIPLS
        adig,adil = self._init_adipls()        
        kwargs.setdefault('cgrav',self.constants['GG'])
        kwargs.setdefault('nfmode',2)
        #-- write adipls variables to file
        fileio.write_gong(adig,adil,model+'.gong')
        subprocess.call('form-amdl.d 2 %s %s'%(model+'.gong',model+'.amdl'),stdout=open('/dev/null','w'),shell=True)
        #special.write_amdl(adig,adil,model+'.amdl')
        #print adig
        #ascii.write_array(adil,'tempmodel_eig.ascii',header=True,auto_width=True,use_float='%.15e')
        l = degree        
        
        #-- run over all the degrees
        for mode in modes:
            mymode = self.eigenfreqs['adipls']['l%d'%(l)][mode]
            if np.isnan(mymode['frequency']):
                logger.info('Mode %d is not matched (nan)'%(mode))
                continue
            #-- use Takata ordering scheme for dipolar modes
            if degree==1:
                kwargs.setdefault('irsord',20)
            #-- use Lee ordering scheme for other modes
            else:
                kwargs.setdefault('irsord',11)
            #-- construct adipls control file
            df = mymode['frequency']/100.
            control_file = adipls.make_adipls_inputfile(model_name=model,degree=degree,
                                               f0=(mymode['frequency']-df,self.unit),
                                               fn=(mymode['frequency']+df,self.unit),**kwargs)
            #-- call adipls
            if not verbose:
                stdout = open('/dev/null', 'w')
            else:
                stdout = None
            subprocess.call('adipls.c.d %s'%(control_file),shell=True,stdout=stdout)
            #-- read eigenfunction file and add contents to output
            try:
                x,a_r,b_r = fileio.read_aeig(model+'.aeig')
            except:
                x,a_r,b_r = [np.linspace(0,1,2000),np.zeros(2000),np.zeros(2000)]
            #os.remove(model+'.aeig')
            os.remove(model+'.gsm')
            os.remove(model+'.ssm')
            #os.remove(control_file)
            Pprime = np.zeros_like(x)
            self.eigenfuncs['adipls']['l%d'%(degree)][mode] = np.rec.fromarrays([x,a_r,b_r,Pprime],names=['x','ar','br','Pprime'])
        #os.unlink(model+'.amdl')
        #os.unlink(model+'.gong')
        
        #raise SystemExit
        logger.info('ADIPLS: computed eigenfunctions for l=%d modes %s'%(degree,modes))
        
    #}

              
#{ General


def match_losc_adipls(losc,adip):
    """
    Match output from LOSC and ADIPLS.
    """
    #-- make unique codes per frequency for easy matching
    code_losc = losc['l']*100+losc['modeLee']
    code_adip = adip['l']*100+adip['nz']
    indices = ne.match_arrays(code_losc,code_adip)
    #-- we found matches for these frequencies
    losc = losc[indices]
    keep = (code_losc[indices]==code_adip)
    return losc[keep],adip[keep]


#}

#{ LOSC    
    
def physical2osc(r,m,P_r,rho_r,G1,G=None):
    """
    Convert the physical stellar model to dimonsionless variants
    which can be used in OSC.
    """
    if G is None:
        G = constants.GG_cgs
    R = r[-1] #-- despite what you think, this should be the total radius, not
              #   the photospheric one (priv comm Andrea Miglio)
    M = m[-1]
    x = r/R
    q = m/M
    
    old_settings = np.seterr(invalid='ignore',divide='ignore')
    S1 = x
    S2 = q/x**3
    S2[0] = np.polyval(np.polyfit(x[1:4],S2[1:4],2),0)
    S3 = R*P_r/(G*M*rho_r)
    S4 = 4*pi*R**3*rho_r/M
    S5 = G1
    S6 = np.zeros(len(r))
    np.seterr(**old_settings)
    
    logger.debug('Converted CGS model to dimensionless LOSC model')
    
    return np.array([S1,S2,S3,S4,S5,S6]),R,M

def star2losc(starg,starl,G=None):
    """
    Convert a star model to LOSC quantities.
    """
    old_settings = np.seterr(invalid='ignore',divide='ignore')
    R,M = starg['photosphere_r'],starg['star_mass']
    if G is None: G = constants.GG_cgs
    if 'x' in starl.dtype.names:
        x = starl['x']
    else:
        x = starl['radius']/R
    S1 = x
    if 'q' in starl.dtype.names:
        q = starl['q']
    else:
        q = starl['mass']/M
    if 'q_x3' in starl.dtype.names:
        S2 = starl['q_x3']
    else:
        S2 = q/x**3
        S2[0] = np.polyval(np.polyfit(x[1:4],S2[1:4],2),0)
    P = starl['pressure']
    rho = starl['Rho']
    S3 = R*P/(G*M*rho)
    S4 = 4*pi*R**3*rho/M
    S5 = starl['gamma1']
    S6 = np.zeros(len(S1))
    np.seterr(**old_settings)
    
    logger.debug('Converted CGS model to dimensionless LOSC model')
    return np.array([S1,S2,S3,S4,S5,S6]),R,M
        
    
    
def osc2physical(columns,R,M,G=None):
    """
    Convert the stellar model from dimensionless variants as used
    in OSC to real physical parameters.
    """
    if G is None:
        G = constants.GG_cgs
    S1,S2,S3,S4,S5,S6 = columns[:6]
    G1 = S5
    r = S1*R
    m = S2*S1**3*M
    rho_r = S4*M/(4*pi*R**3)
    P_r = S3*G*M*rho_r/R
    A = S6/S1**2
    logger.info('Converted dimensionless OSC model to CGS model')
    
    return r,m,P_r,rho_r,G1,A


def fgong2cles(starg,starl,G=None):
    """
    Convert FGONG file to CLES dat files naming convections.
    """
    if G is None: G = constants.GG_cgs
    #-- interpolate values of center points
    for name in starl.dtype.names:
        if name in ['radius']: continue
        starl[name][0] = np.polyval(np.polyfit(starl['radius'][1:4],starl[name][1:4],2),starl['radius'][0])
    mass = np.exp(starl['ln(m/M)'])*starg['star_mass']
    old_settings = np.seterr(invalid='ignore',divide='ignore')
    dP_dr = -starl['Rho']*G*mass/starl['radius']**2
    dP_dr[0] = 0
    np.seterr(**old_settings)
    starl = pl.mlab.rec_append_fields(starl,['mass','dP_dr'],[mass,dP_dr])
    #-- set the values for the D5 and D6 in adipls in the global parameters
    starg['D5'] = -starg['ddP_drr_c']/starl['gamma1'][0]
    starg['D6'] = -starg['ddrho_drr_c']
    return starg,starl
    

    
    

def get_pressure_perturbation(star_cles,star_osc,oscl,mode):
    """
    Calculate pressure perturbation for radial and nonradial modes.
    
    This is independent of theta and phi, and should thus be inserted in the
    right equation to get the 'real' perturbation at a specific point inside
    the star.
    """
    R = star_cles['radius'][-1]
    M = star_cles['mass'][-1]
    r,m,P,rho,G1 = osc2physical(star_osc,R,M)
    G = constants.G_cgs
    q = m/M
    tdyn = sqrt(R**3/(G*M))
    omega = mode['sigma']*tdyn # angular dimensionless frequency
    sigma = mode['sigma']
    l = mode['l']
    Y,Z,U,V = oscl
    x,S2,S3,S4,S5 = star_osc
    #-- in the case of a radial mode:
    if l == 0:
        dP_P = Z
    #-- in the case of nonradial mode:
    else:
        dP_P = x**l*Z
    return dP_P

def get_density_perturbation(star_cles,star_osc,oscl,mode):
    """
    Compute density perturbation in the quasi-adiabatic approx.
    
    See course notes of C.~Aerts (Stellar pulsations, p. 72) or C.~Aerts et al,
    Asteroseismology, 2009.
    """
    G1 = star_osc[4]
    dP_P = get_pressure_perturbation(star_cles,star_osc,oscl,mode)
    drho_rho = 1./G1   * dP_P
    return drho_rho

def get_temperature_perturbation(star_cles,star_osc,oscl,mode):
    """
    Compute temperature perturbation in the quasi-adiabatic approx.
    """
    G3_1 = star_cles['G3-1']
    R = star_cles['radius'][-1]
    x,S2,S3,S4,S5 = star_osc
    #-- interpolate CLES G_3-1 variable onto OSC grid
    G3_1 = interpol.linear_interpolation(star_cles['radius']/R,G3_1,x)
    drho_rho = get_density_perturbation(star_cles,star_osc,oscl,mode)
    dT_T = G3_1 * drho_rho
    return dT_T

def work_integral(star_cles,star_osc,oscl,mode):
    """
    Calculate integrated work integral in the quasi-adiabatic approx.
    
    We now interpolate onto the OSC grid, maybe we should do it the other way
    around, onto the CLES grid?
    
    WARNING: there might be a mix-up with dimensional and dimensionless units,
    but I guess this will only introduce constant factors, which will propagate
    outside of the integrals and derivatives. I guess the opacity derivatives
    are calculated in cgs, other stuff in dimensionless...
    """
    R = star_cles['radius'][-1]
    x,S2,S3,S4,S5 = star_osc
    lumi = star_cles['L']
    rho,T,X,Z = star_cles['Rho'],star_cles['T'],star_cles['X'],star_cles['Z']
    #-- calculate opacity
    kappa,kappa_rho,kappa_T = opacity(rho,T,X,Z)
    #-- calculate pressure perturbation
    drho_rho = get_density_perturbation(star_cles,star_osc,oscl,mode)
    #-- interpolate G3-1  and opacities/luminosity onto OSC grid
    G3_1  = star_cles['G3-1']
    G3_1  = interpol.linear_interpolation(star_cles['radius']/R,G3_1,x)
    kappa = interpol.linear_interpolation(star_cles['radius']/R,kappa,x)
    kappa_rho = interpol.linear_interpolation(star_cles['radius']/R,kappa_rho,x)
    kappa_T = interpol.linear_interpolation(star_cles['radius']/R,kappa_T,x)
    lumi = interpol.linear_interpolation(star_cles['radius']/R,lumi,x)
    #-- calculate psi_rad and the appropriate derivative
    psi_rad = (kappa_T - 4) * G3_1 + kappa_rho
    dpsi_rad = arman.deriv(x,psi_rad*drho_rho)
    #-- calculate work integral
    Wg = zeros_like(psi_rad)
    for i in range(len(Wg))[1:]:
        Wg[i] = lumi[i] * trapz(G3_1*drho_rho * dpsi_rad,x=x)
    return Wg
   

#}

def test():
    """
        >>> from pylab import show
        >>> p=show()
    """
    import doctest
    doctest.testmod()

if __name__ == "__main__":
    #test()
    from ivs.aux import loggers
    logger = loggers.get_basic_logger()
    
    model_cles = '/STER/100/pieterd/workspace/mesa/microphysics/Bstar/Bstar-0200.dat'
    model_mesa = ''
    
    fileio.read_starmodel(model)
    
    raise SystemExit
    sm = StellarModel(model)
    sm.resample()
    sm.compute_eigenfreqs([0],codes=['losc'])
