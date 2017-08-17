# -*- coding: utf-8 -*-
"""
Read and write files from stellar evolution codes.

We use the MESA naming conventions for the physical quantities as the defaults.
"""
import logging
import struct
import os
import glob
import pyfits
import numpy as np
import pylab as pl
#import pycles
from ivs.inout import ascii
from ivs.inout import fits
from ivs.inout import hdf5
from ivs.units import conversions
from ivs.units import constants
from ivs.aux import numpy_ext as ne
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

logger = logging.getLogger("IVS.STAR.FILEIO")

#{ Mesa

def read_mesa(filename='star.log',only_first=False):
    """
    Read star.log and .data files from MESA.
    
    This returns a record array with the global and local parameters (the latter
    can also be a summary of the evolutionary track instead of a profile if
    you've given a 'star.log' file.
    
    The stellar profiles are given from surface to center.
    
    @param filename: name of the log file
    @type filename: str
    @param only_first: read only the first model (or global parameters)
    @type only_first: bool
    @return: list of models in the data file (typically global parameters, local parameters)
    @rtype: list of rec arrays
    """
    models = []
    new_model = False
    header = None
    #-- open the file and read the data
    with open(filename,'r') as ff:
        #-- skip first 5 lines when difference file
        if os.path.splitext(filename)[1]=='.diff':
            for i in range(5):
                line = ff.readline()
            models.append([])
            new_model = True
        while 1:
            line = ff.readline()
            if not line: break # break at end-of-file
            line = line.strip().split()
            if not line: continue
            #-- begin a new model
            if line[0]=='1' and line[1]=='2':
                #-- wrap up previous model
                if len(models):
                    model = np.array(models[-1],float).T
                    models[-1] = np.rec.fromarrays(model,names=header)
                    if only_first: break
                models.append([])
                new_model = True
                continue
            #-- next line is the header of the data, remember it
            if new_model:
                header = line
                new_model = False
                continue
            models[-1].append(line)
    if len(models)>1:
        model = np.array(models[-1],float).T
        models[-1] = np.rec.fromarrays(model,names=header)
    
    if not only_first:
        logger.info('MESA log %s read'%(filename))
    
    return models

def list_mesa_data(filename='profiles.index'):
    """
    Return a chronological list of *.data files in a MESA LOG directory
    """
    number,priority,lognr = ascii.read2array(filename,skip_lines=1).T
    logfiles = [os.path.join(os.path.dirname(filename),'profile%d.data'%(nr)) for nr in lognr]
    return number,logfiles
            

def read_mesa_fits(filename,ext=2):
    """
    Read MESA FITS file containing tracks and summary.
    
    @param filename: name of the file
    @type filename: str
    @param ext: extension to read
    @type ext: int
    @return: global parameters, local parameters
    @rtype: dict, recarray
    """
    data,header = fits.read2recarray(filename,ext=ext,return_header=True)
    return header,data


def mesalog2fits(directory,filename=None,order='star_age'):
    """
    Convert the contents of a MESA log directory to one FITS file.
    
    @param directory: directory that contains 'star.log'
    @type directory: str
    @param filename: name of FITS file to append/write the data to
    @type filename: str
    @param order: order model via this keyword from the global parameters.
    @type order: str
    """
    #-- collect all the files
    files = np.sort(glob.glob('*.data'))
    sa = np.argsort([read_mesa(ff,only_first=True)[0][order][0] for ff in files])
    files = files[sa]
    if filename is None:
        filename = os.path.abspath(directory).split(os.sep)[-1]+'.fits'
    #-- first write summary, then all the models
    files = ['star.log']+list(files)
    for i,ff in enumerate(files):
        starg,starl = read_mesa(ff)
        header_dict = {}
        for name in starg.dtype.names:
            header_dict['HIERARCH '+name] = starg[name][0]
        if i==0:
            header_dict['EXTNAME'] = 'track'
        else:
            header_dict['EXTNAME'] = '%d'%(starg['model_number'])
        filename = fits.write_recarray(starl,filename,header_dict=header_dict,
                                         ext=header_dict['EXTNAME'],close=False)
    filename.close()


#}

#{ CLES

def read_cles_binary(filename):
    """
    Python version of FORTRAN routine.
    
    This one works on 64-bit system too...
    
    Returns a dictionary of global parameters and a record array of local
    parameters.
    
    @param filename: filename
    @type filename: str
    @return: global parameters, local parameters
    @rtype: dict, recarray
    """
    ff = open(filename,'rb')
    content = ff.read()
    ff.close()
    #-- read in globals
    starg = '<'+'i'+24*'d'
    head = struct.calcsize(starg)
    header = content[:head]
    content = content[head:]
    starg = struct.unpack(starg,header)[1:]
    #-- read in locals
    N = int(starg[8])
    star = np.zeros((44,N))
    for i in range(N):
        fmt = '<'+47*'d'
        size = struct.calcsize(fmt)
        line = content[:size]
        content = content[size:]
        star[:,i] = struct.unpack(fmt,line)[3:]
    return starg,star

def read_cles_dat(filename,do_standardize=False):
    """
    Read a CLES dat model (version 18).
    
    @param filename: name of the file
    @type filename: str
    @param do_standardize: standardize output
    @type do_standardize: bool
    @return: global parameters, local parameters
    @rtype: dict, recarray
    """
    #-- definition of the local columns
    cols_local = ['r','m','rho','T','L','P','Cv','G1','G3-1',\
               'Prho','PT','Cp','kappa','enuc','X','Z','dP/dr',\
               'dT/dr','(dlnT/dlnP)rad','(dlnT/dlnP)ad','dlnT/dlnP',\
               'GammaConv','alphaConv','dm1','dm2','xH1','xH2',\
               'xHe3','xHe4','xLi7','xC12','xC13','xN14','xN15',\
               'xO16','xO17','xO18','xNe20','Others','ZOthers',\
               'ZZOthers','MOthers','egrav','mix']
    #-- definition of the global columns
    cols_global = ['M','MMsun','R','Teff','L','LLsun','RRsun',
                   'age','np','logg','Lg/L','agey','atm',
                   'version','alphaConv','AlphaOvershooting',
                   'X0','Z0','TimeStep','DiffusionParameters',
                   'EOS type','Opacity table','Metal table','Step number']
    #-- read in model (you can exchange this with the FORTRAN program if you
    #   like)
    starg_,star_ = read_cles_binary(filename)#pycles.readmodel(filename)
    starg = {}
    star = []
    keep = star_[3,:]!=0
    for i,col in enumerate(cols_local):
        star.append(star_[i,keep])
    for i,col in enumerate(cols_global):
        starg[col] = starg_[i]
    
    star = np.rec.fromarrays(star,names=cols_local)
    logger.info("CLES model %s read"%(filename))
    
    if do_standardize:
        starg = standardize(starg,'cles','global')
        starl = standardize(starl,'cles','local')
    
    return starg,star


def read_cles_sum(sumfile):
    """
    Read the input from a CLES '-sum.txt' file.
    
    Keys: C{logTeff}, C{logg}, C{M}, C{R}, C{logL}, C{age}(, C{num}, C{Xc})
    
    @param sumfile: filename of the '-sum.txt' file.
    @type param:str
    @return: model track
    @rtype: recarray
    """
    c = {'summary':{'logTeff':7,'logg':11,'M':2,'R':10,'logL':8,'age':6},
         'sum':    {'logTeff':2,'logg': 8,'M':6,'R': 7,'logL':3,'age':1,'num':0,'Xc':4}}
    
    data = ascii.read2array(sumfile)
    
    if 'summary' in sumfile:
        t = 'summary'
    else:
        t = 'sum'
        #-- not all info readily available, build it
        data_ = np.zeros((len(data[:,0]),3))
        mass = float(os.path.basename(sumfile).split('_')[0][1:])
        data_[:,0] = mass
        radi = conversions.derive_radius((data[:,c[t]['logL']],'[Lsol]'),\
                                         (data[:,c[t]['logTeff']],'[K]'),units='Rsol')[0]
        data_[:,1] = radi
        logg_ = conversions.derive_logg((float(mass),'Msol'),(data_[:,1],'Rsol'))
        data_[:,2] = logg_
        data = np.hstack([data,data_])
    
    keys = np.array(c[t].keys())
    cols = np.array([c[t][key] for key in keys])
    sa = np.argsort(cols)
    keys = keys[sa]
    cols = cols[sa]
    mydata = np.rec.fromarrays([data[:,col] for col in cols],names=list(keys))
    
    logger.debug('CLES summary %s read'%(sumfile))
    
    return mydata
#}
#{ ADIPLS

def read_assm(ssm_file):
    """
    Read in a short summary file from ADIPLS.
    
    @param ssm_file: name of the short summary file
    @type ssm_file: str
    @return: global parameters, local parameters
    @rtype: dict, rec array
    """
    with open(ssm_file,'rb') as ff:
        contents = ff.read()
    
    #-- global parameters
    fmt = '<i' + 6*'d'
    size = struct.calcsize(fmt)
    out = struct.unpack(fmt,contents[:size])
    contents = contents[size:]
    M,R,P_c,rho_c = out[3:7]
    starg = dict(M=M,R=R,P_c=P_c,rho_c=rho_c)
    
    #-- local parameters
    fmt = '<'+4*'i'+5*'d'+'ii'
    size = struct.calcsize(fmt)
    starl = np.zeros((5,len(contents)/size))
    for i in range(starl.shape[1]):
        out = struct.unpack(fmt,contents[:size])
        contents = contents[size:]
        starl[:,i] = out[4:-2]
    
    starl = np.rec.fromarrays(starl,names=['l','nz','sigma2','E','frequency'])
    return starg,starl
    
def read_agsm(ssm_file):
    """
    Read in a grand summary file from ADIPLS.
    
    @param ssm_file: name of the grand summary file
    @type ssm_file: str
    @return: global parameters, local parameters
    @rtype: dict, rec array
    """
    with open(ssm_file,'rb') as ff:
        contents = ff.read()
    
    fmt = '<i'
    size = struct.calcsize(fmt)
    out = struct.unpack(fmt,contents[:size])
    contents = contents[size:]
    
    pars = []
    while 1:
        #-- global parameters
        fmt = '<' + 38*'d'
        size = struct.calcsize(fmt)
        out = struct.unpack(fmt,contents[:size])
        contents = contents[size:]
        
        mode = list(out)
        
        fmt = '<'+8*'i'
        size = struct.calcsize(fmt)
        out = struct.unpack(fmt,contents[:size])
        contents = contents[size:]
        mode += list(out)[:-1]
        pars.append(mode)
        
        if len(contents)<72:
            break
        fmt = '<'+18*'i'
        size = struct.calcsize(fmt)
        out = struct.unpack(fmt,contents[:size])
        contents = contents[size:]
    
    columns = np.array(['xmod','M','R','P_c','Rho_c','D5','D6','mu','D8',
               'A2_xs','A5_xs','x1','sigma_Omega','xf','fctsbc','fcttbc',
               'lambda','l','n','sigma2','sigmac2','y1max','xmax',
               'E','Pe','PV','varfreq','ddsig','ddsol','y1_xs','y2_xs','y3_xs','y4_xs',
               'z1max','xhatmax','beta_nl','freqR','m','in','nn','rich','ivarf',
               'icase','iorign','iekinr'])
    keep = np.array([17,18,19,23,26,35,36,37])
    data = np.rec.fromarrays(np.array(pars).T[keep],names=list(columns[keep]))
    
    fmt = '<'+17*'i'
    size = struct.calcsize(fmt)
    out = struct.unpack(fmt,contents[:size])
    contents = contents[size:]
    
    return data
    
def read_aeig(aeig_file,nfmode=2):
    """
    Read in a simple eigenfunction file written by ADIPLS.
    
    @param aeig_file: name of the short summary file
    @type aeig_file: str
    @return: x, ar, br
    @rtype: 3xarray
    """
    with open(aeig_file,'rb') as ff:
        contents = ff.read()
    
    if nfmode==2:
        #-- read in the size of the arrays in the file
        fmt = '<'+2*'i'
        size = struct.calcsize(fmt)
        N1,N2 = struct.unpack(fmt,contents[:size])
        contents = contents[size:]
        #-- read in the 'x' axis
        fmt = '<'+N2*'d'
        size = struct.calcsize(fmt)
        x = struct.unpack(fmt,contents[:size])
        contents = contents[size:]
        x = np.array(x)
        #-- read in the size of the following arrays in the file
        fmt = '<'+2*'i'
        size = struct.calcsize(fmt)
        out = struct.unpack(fmt,contents[:size])
        contents = contents[size:]
        #-- some weird intermezzo that I don't understand: global parameters?
        fmt = '<'+50*'d'
        size = struct.calcsize(fmt)
        y = struct.unpack(fmt,contents[:size])
        contents = contents[size:]
        #-- the actual eigenfunctions
        N = 2
        fmt = '<'+N*N2*'d'
        size = struct.calcsize(fmt)
        z = struct.unpack(fmt,contents[:size])
        contents = contents[size:]
        z = np.array(z)
        y1 = z[0::N]
        y2 = z[1::N]
        #-- and an end integer
        fmt = '<'+'i'
        size = struct.calcsize(fmt)
        y = struct.unpack(fmt,contents[:size])
        contents = contents[size:]
        output = x,y1,y2
    
    elif nfmode==1:
        fmt = '<'+1*'i'
        size = struct.calcsize(fmt)
        N1 = struct.unpack(fmt,contents[:size])
        contents = contents[size:]
        fmt = '<'+50*'d'
        size = struct.calcsize(fmt)
        y = struct.unpack(fmt,contents[:size])
        contents = contents[size:]
        #-- read in the size of the following arrays in the file
        fmt = '<'+1*'i'
        size = struct.calcsize(fmt)
        out = struct.unpack(fmt,contents[:size])
        contents = contents[size:]
        N2, = out
        #-- the actual eigenfunctions
        N = 7
        fmt = '<'+N*N2*'d'
        size = struct.calcsize(fmt)
        z = struct.unpack(fmt,contents[:size])
        contents = contents[size:]
        z = np.array(z)
        x = z[0::N]
        y1 = z[1::N]
        y2 = z[2::N]
        y3 = z[3::N]
        y4 = z[4::N]
        z1 = z[5::N]
        z2 = z[6::N]
        #-- and an end integer
        fmt = '<'+'i'
        size = struct.calcsize(fmt)
        y = struct.unpack(fmt,contents[:size])
        contents = contents[size:]
        output = x,y1,y2,y3,y4,z1,z2
    
    if len(contents):
        raise IOError,"Error reading ADIPLS eigenfrequency file: uninterpreted bytes remaining"
    return output

def read_amdl(filename):
    """
    Read unformatted binary AMDL files from ADIPLS.
    
    @param filename: name of the file
    @type filename: str
    @return: global parameters, local parameters
    @rtype: dict, rec array
    """
    with open(filename,'rb') as ff:
        #-- read in all the bytes at once
        content = ff.read()
    #-- have a look at the header of the file
    header_fmt = '<iii'
    line_fmt = '<'+8*'d'
    head = struct.calcsize(header_fmt)
    line = struct.calcsize(line_fmt)
    header = content[:head]
    content = content[head:]
    print struct.unpack(header_fmt,header)
    model_name,icase,nn = struct.unpack(header_fmt,header)
    print filename,'model_name,nn',model_name,nn
    
    #-- read in the first line of the data. This contains
    #   global information on the model, and also the number
    #   of columns in the datafile.
    record = content[:line]
    content = content[line:]
    M,R,P_c,rho_c,D5,D6,mu,D8 = struct.unpack(line_fmt,record)
    
    #-- print out some information
    M_ = conversions.convert('g','Msol',M)
    R_ = conversions.convert('cm','Rsol',R)
    print "Model information: (%d meshpoints)"%(nn)
    print "  M=%.3gMsol, R=%.3gRsol"%(M_,R_)
    print "  P_c=%.3g g/cm/s2, rho_c=%.3g g/cm3"%(P_c,rho_c)
    print '  Polytropic index: %.2g'%(mu)
    print '  %.3g %3g'%(D5,D6)
    
    #-- derive the number of variables in this model file
    if D8>=100: nvars = 9
    elif D8>=10:nvars = 7
    else:       nvars = 6
    
    #-- then read in the rest, line by line
    line_fmt = '<'+nvars*'d'
    line = struct.calcsize(line_fmt)
    starl = []
    for i in range(nn):
        record = content[:line]
        content = content[line:]
        data = struct.unpack(line_fmt,record)
        starl.append(data)
    
    #-- end of the file
    line_fmt = '<i'
    record = content[:line]
    content = content[line:]
    data = struct.unpack(line_fmt,record)
    print 'endoffile',data
    
    #-- wrap up the local and global information on the star
    starl = np.rec.fromarrays(np.array(starl).T,names=['x','q/x3','Vg','G1','A','U'])
    starg = {'M':M,'R':R,'P_c':P_c,'rho_c':rho_c,
             'D5':D5,'D6':D6,'mu':mu}
    if len(content):
        raise ValueError,'uninterpreted data in file',fileanem
    
    return starg,starl



def write_gong(adig,adil,filename):
    """
    Write star info to gong file
    
    @param adig: global parameters of the star
    @type adig: dict
    @param adil: local parameters of the star
    @type adil: rec array.
    @param filename: name of the file
    @type filename: str
    @return: global parameters, local parameters
    @rtype: dict, rec array
    """
    contents = np.ravel(np.array([adil[name] for name in adil.dtype.names]).T)
    with open(filename,'w') as ff:
        ff.write('%10d%10d%10d\n'%(1,len(adil),len(adil.dtype.names)-1))
        ff.write('%20.13e'%(adig['star_mass']))
        ff.write('%20.13e'%(adig['photosphere_r']))
        ff.write('%20.13e'%(adig['P_c']))
        ff.write('%20.13e\n'%(adig['rho_c']))
        ff.write('%20.13e'%(adig['D5']))
        ff.write('%20.13e'%(adig['D6']))
        ff.write('%20.13e'%(adig['mu']))
        ff.write('%20.13e'%(adig.get('D8',0)))
        for i in xrange(len(contents)):
            if i%4==0: ff.write('\n')
            ff.write('%20.13e'%(contents[i]))


def write_gyre(starg,starl,filename):
    """
    Write a file to the GYRE format.
    """
    gyre_model = {}
    
    gyre_model['L_star'] = float(conversions.convert('erg/s',"SI",starg['photosphere_L']))
    gyre_model['R_star'] = float(conversions.convert('cm',"SI",starg['photosphere_r']))
    gyre_model['X_star'] = -1.00
    gyre_model['Z_star'] = float(starg['initial_z'])
    gyre_model['M_star'] = float(conversions.convert('g','SI',starg['star_mass']))
    gyre_model['t_star'] = float(starg['star_age'])
    gyre_model['n_shells'] = len(starl)
    gyre_model['B3_VERSION'] = 1.0
    gyre_model['B3_SUBCLASS'] = ''
    gyre_model['B3_CLASS'] = 'STAR'
    gyre_model['B3_DATE'] = ''
    #for key in gyre_model:
        #print key,gyre_model[key]
    
    #-- local parameters
    #inner = starl['mass']<starg['star_mass']
    #w = -starl['mass'][inner]/(starl['mass'][inner]-starg['star_mass'])
    inner = starl['mass']<starl['mass'].max()
    w = -starl['mass'][inner]/(starl['mass'][inner]-starl['mass'].max())
    gyre_model['r'] = conversions.convert('cm','SI',starl['radius'])
    gyre_model['w'] =  np.hstack([w,1e16*np.ones(sum(-inner))])
    gyre_model['L_r'] = conversions.convert('erg/s','SI',starl['luminosity'])
    gyre_model['L_r'][0] = 0.
    gyre_model['T'] = starl['temperature']
    gyre_model['p'] = conversions.convert('ba','Pa',starl['pressure'])
    gyre_model['c_V'] = conversions.convert('erg/s/g/K','SI',starl['cv'])
    gyre_model['c_p'] = conversions.convert('erg/s/g/K','SI',starl['cp'])
    gyre_model['X'] = starl['h1']
    gyre_model['rho'] = conversions.convert('g/cm3','SI',starl['Rho'])
    gyre_model['chi_rho'] = starl['chiRho']
    gyre_model['chi_T'] = starl['chiT']
    gyre_model['N2'] = starl['brunt_N2']
    gyre_model['nabla'] = starl['grad_temperature']
    gyre_model['kappa'] = conversions.convert('cm2/g','SI',starl['opacity'])
    gyre_model['kappa_rho'] = starl['dkap_dlnrho']/starl['opacity']
    gyre_model['kappa_T'] = starl['dkap_dlnT']/starl['opacity']
    gyre_model['epsilon'] = np.zeros(len(starl))
    gyre_model['epsilon_rho'] = np.zeros(len(starl))
    gyre_model['epsilon_T'] = np.zeros(len(starl))
    
    #-- global parameters
    gyre_model['L_star'] = float(gyre_model['L_r'].max())
    gyre_model['R_star'] = float(gyre_model['r'].max())
    gyre_model['M_star'] = float(starl['mass'].max()/1000.)
    #gyre_model['L_star'] = float(gyre_model['L_r'].max())#float(conversions.convert('erg/s',"SI",starg['photosphere_L']))
    #gyre_model['R_star'] = float(gyre_model['r'].max())#float(conversions.convert('cm',"SI",starg['photosphere_r']))
    #gyre_model['X_star'] = -1.00
    #gyre_model['Z_star'] = float(starg['initial_z'])
    #gyre_model['M_star'] = float(starl['mass'].max())#float(conversions.convert('g','SI',starg['star_mass']))
    #gyre_model['t_star'] = float(starg['star_age'])
    #gyre_model['n_shells'] = len(starl)
    #gyre_model['B3_VERSION'] = 1.0
    #gyre_model['B3_SUBCLASS'] = ''
    #gyre_model['B3_CLASS'] = 'STAR'
    #gyre_model['B3_DATE'] = ''
    
    #keys = ['w','L_r','T','p','c_V','c_p','X','rho','chi_rho','chi_T','N2',\
            #'nabla','kappa','kappa_rho','kappa_T','epsilon','epsilon_rho','epsilon_T']
    ##for key in keys:
        #gyre_model[key] = np.array(gyre_model[key],np.float32)
        #print key,gyre_model[key].min(),gyre_model[key].max()
    
    hdf5.write_dict(gyre_model,filename,update=False,attr_types=[float,int,str,unicode])
    
#}
#{ General

def read_starmodel(filename,do_standardize=True,units='cgs',values='auto',**kwargs):
    """
    Automatically detect the format of the stellar evolution file.
    
    This function also adds additional information that can be derived from the
    given information, which is sometimes handy for plotting reasons or to
    compute frequencies.
    
    @param filename: name of the file
    @type filename: str
    @param do_standardize: standardize output
    @type do_standardize: bool
    @param units: work in this convention's units (cgs, SI...)
    @type units: str
    @param values: use this set of values for the fundamental constants
    @type values: str (mesa, standard, cles...)
    @return: global parameters, local parameters
    @rtype: dict, rec array
    """
    add_center_point = kwargs.get('add_center_point',False)
    old_settings = np.seterr(invalid='ignore',divide='ignore')
    if isinstance(filename,str):
        basename,ext = os.path.splitext(filename)
    elif isinstance(filename,pyfits.HDUList):
        ext = '.fits'
    else:
        starg,starl = filename
        ext = None
        _from = None
    #-- MESA
    if ext in ['.log','.data']:
        starg,starl = read_mesa(filename)
        _from = 'mesa'
    elif ext in ['.fits']:
        starg,starl = read_mesa_fits(filename,ext=kwargs.get('ext',2))
        _from = 'mesa'
    elif ext in ['.dat']:
        starg,starl = read_cles_dat(filename)
        _from = 'cles'
    
    if do_standardize:
        #-- if values is set to automatic, set the values for the fundamental
        #   constants to the set defined in the pulsation/stellar evolution code
        if values=='auto' and _from is not None:
            values=_from
        old = conversions.set_convention(units,values)
        starg = standardize(starg,_from,'global',units=units)
        starl = standardize(starl,_from,'local',units=units,add_center_point=add_center_point)
        if ext in ['.log']:
            starl = starl[::-1]
        #-- now just add some other parameters for which we need both the
        #   global and local parameters:
        available = starl.dtype.names
        if not 'q' in available and 'mass' in available and 'star_mass' in starg.keys():
            q = starl['mass']/starg['star_mass']
            starl = pl.mlab.rec_append_fields(starl,['q'],[q])
            logger.info('Added local q')
        if not 'x' in available and 'radius' in available and 'photosphere_r' in starg.keys():
            x = starl['radius']/starg['photosphere_r']
            starl = pl.mlab.rec_append_fields(starl,['x'],[x])
            logger.info('Added local x')
        if not 'q_x3' in available and 'mass' in available and 'star_mass' in starg.keys():
            q_x3 = starl['q']/starl['x']**3
            #-- some models contain the center values, others do not
            if ext not in ['.data']:
                q_x3[0] = starl['Rho'][0]/starg['star_mass']*starg['photosphere_r']**3*4./3.*np.pi
            else:
                q_x3[0] = np.polyval(np.polyfit(starl['x'][1:5]**2,q_x3[1:5],2),0)
                #starl['Rho'][0] = q_x3[0]*starg['star_mass']/starg['photosphere_r']**3/4.*3./np.pi
                U = starl['Rho']/q_x3
                starl['Rho'][0] = np.polyval(np.polyfit(starl['x'][1:5]**2,U[1:5],1),0)*q_x3[0]
                
            starl = pl.mlab.rec_append_fields(starl,['q_x3'],[q_x3])
            logger.info('Added local q/x3')
        #-- add some quantities that are useful for adipls
        if not 'P_c' in starg.keys() and 'pressure' in available:
            if ext in ['.data']:
                starl['pressure'][0] = starl['pressure'][1] + starl['Rho'][1]*starl['g'][1]*starl['radius'][1]
            starg['P_c'] = starl['pressure'][0]
            logger.info('Added global P_c')
        if not 'rho_c' in starg.keys() and 'Rho' in available:
            starg['rho_c'] = starl['Rho'][0]
            logger.info('Added global rho_c')
        if (not 'D5' in starg.keys() or starg['D5']==0) and 'gamma1' in available:
            
            #dp_dx = ne.deriv(x,starl['pressure'])
            #d2p_dx2 = ne.deriv(x,dp_dx)
            #drho_dx = ne.deriv(x,starl['Rho'])
            #d2rho_dx2 = ne.deriv(x,drho_dx)
            #starg['D5'] = (-1./(starl['gamma1']*starl['pressure'])*d2p_dx2)[0]
            #starg['D6'] = (-1./starl['Rho']*d2rho_dx2)[0]
            
            p1 = np.polyfit(starl['x'][:5]**2,starl['Vg'][:5],1)
            p2 = np.polyfit(starl['x'][:5]**2,starl['brunt_A'][:5],1)
            starg['D5'] = p1[0]
            starg['D6'] = p2[0]+p1[0]
            
            starg['mu'] = -1
            logger.info('Set global parameter D5=%.6g for ADIPLS'%(starg['D5']))
            logger.info('Set global parameter D6=%.6g for ADIPLS'%(starg['D6']))
            logger.info('Set global parameter mu=%.6g for ADIPLS'%(starg['mu']))            
        logger.info("Standardized model")
    
        #-- reset convention to what it was before
        if values=='auto' and _from is not None:
            conversions.set_convention(old[0],old[1])
    
    np.seterr(**old_settings)
    return starg,starl,_from

def read_cols(filename):
    """
    Read a .cols file.
    
    @param filename: name of the file
    @type filename: str
    @return: names, units and description of columns in a profile file
    @rtype: 3xlist
    """
    name = []
    unit = []
    desc = []
    with open(filename,'r') as ff:
        while 1:
            line = ff.readline().strip()
            if not line: break
            if line[0]=='!': continue
            contents = line.split('!')
            name.append(contents[0].strip())
            unit.append(contents[1].strip())
            desc.append(contents[2].strip())
    return name,unit,desc









def standardize(star,_from,info='local',units='cgs',add_center_point=False):
    """
    Standardize units and naming of a stellar model.
    
    @param star: information to standardize
    @type star: dict or rec array
    @param _from: name of the type of file (mesa, cles...)
    @type _from: str
    @param info: is this local or global information?
    @type info: str
    @param units: units to work with
    @type units: str
    @param add_center_point: if True, a center point will be added
    @type add_center_point: bool
    @return: standardized star
    @rtype: dict or recarray
    """
    basedir = os.path.dirname(__file__)
    if _from=='mesa' and info=='global' and not hasattr(star,'keys'):
        #-- convert global recarray into dictionary
        star_ = {}
        for name in star.dtype.names:
            star_[name] = star[name][0]
        star = star_
    
    if hasattr(star,'keys'):
        keys = star.keys()
    else:
        keys = star.dtype.names
    
    if _from is not None:
        mesa_onames,mesa_names,mesa_units = read_cols(os.path.join(basedir,'mesa_{0}.cols'.format(info)))
        mnames,onames,ounits = read_cols(os.path.join(basedir,_from+'_{0}.cols'.format(info)))        
        
        new_star = {}
        for key in keys:
            #-- maybe we cannot convert the column (yet):
            try:
                index = onames.index(key)
            except ValueError:
                logger.debug('Cannot convert %s from (%s,%s) to mesa name (skipping)'%(key,_from,info))
                continue
            #-- maybe information is not in 'standard' list. If it is, convert to
            #   standard units and standard name. Otherwise, just change the name
            try:
                if ounits[index]=='ignore':
                    continue
                mindex = mesa_names.index(mnames[index])
                if ounits[index]=='n' and mesa_units[mindex]=='n kg':
                    value = star[key]*getattr(constants,mnames[index].title())
                elif mesa_units[mindex] or ounits[index]:
                    value = conversions.convert(ounits[index],units,star[key])
                else:
                    value = star[key]
            except ValueError:
                #logger.warning('Cannot convert %s to mesa units (keeping original units)'%(key))
                value = star[key]
            
            new_star[mnames[index]] = value
        #-- if input was not a dictionary, convert it to a record array
        if not hasattr(star,'keys'):
            names = new_star.keys()
            data = [new_star[name] for name in names]
            new_star = np.rec.fromarrays(data,names=names)
            keys = names
        else:
            keys = star.keys()
    else:
        if not hasattr(star,'keys'):
            names = star.dtype.names
            data = [star[name] for name in names]
            new_star = np.rec.fromarrays(data,names=names)
            keys = list(names)
        else:
            new_star = star
            keys = star.keys()
    
    
    
    #-- some profiles start at the surface...
    if info=='local' and _from in ['mesa']:
        new_star = new_star[::-1]
        add_center_point = True
        logger.info('Reversed profile (first point is now center)')
    #-- some profiles have no inner point
    if add_center_point:
        if info=='global':
            raise ValueError,'cannot add center point in global variables'
        new_star = np.hstack([new_star[:1],new_star])
        #-- fill-in some initial values
        for key in ['mass','q','radius']:
            if key in keys:
                new_star[key][0] = 0.
        logger.info('Added center point')
    #-- sometimes it is easy to precompute additional information that can be
    #   derived from the other information available. We might need this
    #   information for pulsation analysis, profile studies etc...
    #-- (surface) gravity
    if not 'g' in keys:
        if info=='local' and 'mass' in keys and 'radius' in keys:
            g = constants.GG*new_star['mass']/new_star['radius']**2
            g[0] = 0.
            new_star = pl.mlab.rec_append_fields(new_star,['g'],[g])
        elif info=='global' and 'star_mass' in keys and 'photosphere_r' in keys:
            new_star['g'] = constants.GG*new_star['star_mass']/new_star['photosphere_r']**2
        logger.info('Derived g from model info (%s, %s)'%(_from,info))
    #-- dP/dr
    if info=='local' and not 'dP_dr' in keys and 'Rho' in keys and 'mass' in keys and 'radius' in keys:
        dP_dr = -new_star['Rho']*constants.GG*new_star['mass']/new_star['radius']**2 #dP/dr=-rho g
        dP_dr[0] = 0
        new_star = pl.mlab.rec_append_fields(new_star,['dP_dr'],[dP_dr])
        logger.info('Derived dP_dr from model info (%s, %s)'%(_from,info))
        keys.append('dP_dr')
    #-- Brunt A
    if info=='local' and not 'brunt_A' in keys and 'brunt_N2' in keys:
        A = new_star['brunt_N2']*new_star['radius']/new_star['g']
        new_star = pl.mlab.rec_append_fields(new_star,['brunt_A'],[A])
        logger.warning('Brunt A via brunt N2')
        keys.append('brunt_A')
    if info=='local' and not 'brunt_A' in keys and 'radius' in keys and 'Rho' in keys \
        and 'gamma1' in keys and 'pressure' in keys and 'dP_dr' in keys:
        r,rho,G1,P = new_star['radius'],new_star['Rho'],new_star['gamma1'],new_star['pressure']
        dP_dr = new_star['dP_dr']
        x = new_star['mass']/new_star['mass'][-1]
        A = -ne.deriv(r,rho)/rho*r + 1./G1 * dP_dr/P*r
        N2 = new_star['g']*A/r
        N2[0] = 0
        A[:2] = np.polyval(np.polyfit(x[2:5]**2,A[2:5],2),x[:2]**2)
        new_star = pl.mlab.rec_append_fields(new_star,['brunt_A','brunt_N2'],[A,N2])
        logger.warning('Brunt A/N via explicit derivatives')
        logger.info('Extrapolated inner two points of Brunt A')
    #-- override brunt_N2_dimensionless:
    if info=='local' and 'grad_density' in keys and 'radius' in keys and 'Rho' in keys \
        and 'gamma1' in keys and 'pressure' in keys and 'mass' in keys:
        g,rho,P,G1 = new_star['g'],new_star['Rho'],new_star['pressure'],new_star['gamma1']
        dlnrho_dlnP = new_star['grad_density']
        N2 = rho*g**2/P*(-1./G1 + dlnrho_dlnP)
        #factor = 3*constants.GG*new_star['mass'].max()/new_star['radius'].max()**3
        #new_star['brunt_N2_dimensionless'] = N2/factor
        x = new_star['mass']/new_star['mass'][-1]
        A = new_star['radius']*N2/g
        A[:2] = np.polyval(np.polyfit(x[2:5]**2,A[2:5],2),x[:2]**2)
        new_star['brunt_A'] = A
        if not 'brunt_N2' in new_star.dtype.names:
            new_star = pl.mlab.rec_append_fields(new_star,['brunt_N2'],[N2])
        else:
            new_star['brunt_N2'] = N2
        logger.info('Brunt A via precomputed partial derivatives (%s,%s)'%(_from,info))
        logger.info('Extrapolated inner two points of Brunt A')
    #-- ADIPLS quantity Vg
    if info=='local' and not 'Vg' in keys and 'gamma1' in keys and 'radius' in keys and 'Rho' in keys \
        and 'gamma1' in keys and 'pressure' in keys and 'mass' in keys:
        m,rho = new_star['mass'],new_star['Rho']
        G1,P,r = new_star['gamma1'],new_star['pressure'],new_star['radius']
        Vg = constants.GG*m*rho/(G1*P*r)
        Vg[0] = 0.
        new_star = pl.mlab.rec_append_fields(new_star,['Vg'],[Vg])
    elif info=='local' and not 'Vg' in keys and not 'gamma1' in keys:
        logger.warning('Cannot compute Vg because gamma1 is not available')
    #-- ADIPLS quantitiy U
    if info=='local' and not 'U' in keys and 'gamma1' in keys and 'radius' in keys and 'Rho' in keys \
        and 'gamma1' in keys and 'pressure' in keys and 'mass' in keys:
        rho,m,r = new_star['Rho'],new_star['mass'],new_star['radius']
        U = 4*np.pi*rho*r**3/m
        x = new_star['mass']/new_star['mass'][-1]
        U[0] = np.polyval(np.polyfit(x[1:4],U[1:4],2),x[0])
        new_star = pl.mlab.rec_append_fields(new_star,['U'],[U])
    #-- dynamical time scale
    if info=='global' and not 'dynamic_time' in keys and 'photosphere_r' in keys and 'star_mass' in keys:
        R = new_star['photosphere_r']
        M = new_star['star_mass']
        new_star['dynamic_time'] = np.sqrt(R**3/(constants.GG*M))
        logger.info('Derived dynamic_time from model info (%s, %s)'%(_from,info))
    if info=='local' and not 'csound' in keys and 'gamma1' in keys and 'radius' in keys and 'Rho' in keys \
        and 'gamma1' in keys and 'pressure' in keys and 'mass' in keys:
        csound = np.sqrt(new_star['gamma1']*new_star['pressure']/new_star['Rho'])
        new_star = pl.mlab.rec_append_fields(new_star,['csound'],[csound])
    #-- some profiles don't have a center value
    if _from in []:
        pass
    
    #-- that's it!
    return new_star
        
#}


#{ General plotting routines

def plot_logRho_logT(starl):
    """
    Plot Density -Temperature diagram for one given profile.
    """
    #-- list all burning regions
    pl.xlabel('log (Density [g cm$^{-1}$]) [dex]')
    pl.ylabel('log (Temperature [K]) [dex]')
    for species in ['hydrogen','helium','carbon','oxygen']:
        logRho,logT = ascii.read2array(os.path.join(os.path.dirname(__file__),'plot_info','%s_burn.data'%(species))).T
        pl.plot(logRho,logT,'k--',lw=2)
        pl.annotate(species,(logRho[-1],logT[-1]),ha='right')
    
    bounds = ['elect','gamma_4_thirds','kap_rad_cond_eq','opal_clip','psi4','scvh_clip']
    bounds = ['gamma_4_thirds','opal_clip','scvh_clip']
    names = ['$\Gamma$<4/3','OPAL','SCVH']
    for limits,name in zip(bounds,names):
        logRho,logT = ascii.read2array(os.path.join(os.path.dirname(__file__),'plot_info','%s.data'%(limits))).T
        pl.plot(logRho,logT,'--',lw=2,color='0.5')
        xann,yann = logRho.mean(),logT.mean()
        pl.annotate(name,(xann,yann),color='0.5')
    #-- plot stellar profile
    color_radiative = 'g'
    color_convective = (0.33,0.33,1.00)

    x = np.log10(starl['Rho'])
    y = np.log10(starl['temperature'])
    
    # first convective regions and radiative regions
    stab_type = starl['stability_type']
    # Create a colormap for convective and radiative regions
    cmap = ListedColormap([color_radiative,color_convective])
    norm = BoundaryNorm([0, 1, 5], cmap.N)
    # Create a set of line segments so that we can color them individually
    # This creates the points as a N x 1 x 2 array so that we can stack points
    # together easily to get the segments. The segments array for line collection
    # needs to be numlines x points per line x 2 (x and y)
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    # Create the line collection object, setting the colormapping parameters.
    # Have to set the actual values used for colormapping separately.
    lc = LineCollection(segments, cmap=cmap, norm=norm)
    lc.set_array(stab_type)
    lc.set_linewidth(10)
    pl.gca().add_collection(lc)
    
    # burning regions
    limits = [0,1,1e3,1e7,np.inf]
    burn = starl['eps_nuc'][:-1]
    colors = np.zeros((len(segments),4))
    colors[:,0] = 1
    colors[:,-1] = 1
    colors[(burn<=1e0),-1] = 0
    colors[(1e0<burn)&(burn<=1e3),1] = 1
    colors[(1e3<burn)&(burn<=1e7),1] = 0.5
    
    # Create the line collection object, setting the colormapping parameters.
    # Have to set the actual values used for colormapping separately.
    lc = LineCollection(segments,colors=colors)
    lc.set_linewidth(4)
    pl.gca().add_collection(lc)
    
    p1 = pl.Line2D([0,0],[1,1], color=(1,1.0,0),lw=3)
    p2 = pl.Line2D([0,0],[1,1], color=(1,0.5,0),lw=3)
    p3 = pl.Line2D([0,0],[1,1], color=(1,0,0),lw=3)
    p4 = pl.Line2D([0,0],[1,1], color='g',lw=3)
    p5 = pl.Line2D([0,0],[1,1], color=(0.33,0.33,1.0),lw=3)
    leg = pl.legend([p1,p2,p3,p4,p5], [">1 erg/s/g",">10$^3$ erg/s/g",">10$^7$ erg/s/g",'Radiative','Convective'],loc='best',fancybox=True)
    leg.get_frame().set_alpha(0.5)
    
    
    pl.xlim(-10,11)
    pl.ylim(4.0,9.4)

def plot_Kiel_diagram(starl):
    """
    Plot Kiel diagram.
    """
    x = starl['temperature']
    y = starl['g']
    age = starl['age']/1e6
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    
    cmap = pl.cm.spectral
    norm = pl.Normalize(age.min(), age.max())
    lc = LineCollection(segments, cmap=cmap,norm=norm)
    lc.set_array(age)
    lc.set_linewidth(3)
    pl.gca().add_collection(lc)
    pl.xlim(x.max(), x.min())
    pl.ylim(y.max(), y.min())
    pl.xlabel('Effective temperature [K]')
    pl.ylabel('log(surface gravity [cm s$^{-2}$]) [dex]')
    
    ax0 = pl.gca()
    ax1 = pl.mpl.colorbar.make_axes(ax0)[0]
    norm = pl.mpl.colors.Normalize(age.min(), age.max())
    cb1 = pl.mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
                                   norm=norm,orientation='vertical')
    cb1.set_label('Age [Myr]')
    pl.axes(ax0)
    
def plot_abundances(starl,x='q',species=('h1','he4','c12','o16','n14')):
    """
    Plot abundance profile
    """
    colors = dict(h1='b',he4='g',c12='k',o16='r',n14=(0.58,0.,0.83))
    for spec in species:
        pl.plot(starl[x],starl[spec],color=colors[spec],lw=2,label=spec)
    leg = pl.legend(loc='best',fancybox=True)
    leg.get_frame().set_alpha(0.5)
    pl.xlabel(x)
    pl.ylabel('Element mass fraction')



def plot_add_convection(starl,x='q'):
    stabtype = starl['stability_type']
    x = starl[x]
    #-- search for the limits of convection types
    limits = [[stabtype[0],x[0]]]
    for i in range(1,len(x)):
        if stabtype[i]!=limits[-1][0]:
            limits[-1].append(x[i])
            limits.append([stabtype[i],x[i]])
    limits[-1].append(x[-1])
    alphas = [0.1,0.25,0.1,0.5,0.5]
    colors = ['k','y','k','g','r'] # no conv / conv /  / only ledou / only sch
    alphas = [1.0,0.5,0.5,0.5,0.5]
    colors = ['k',(0.33,0.33,1.0),(0.33,0.33,1.0),(0.33,0.33,1.0),(0.33,0.33,1.0)]
    for limit in limits:
        if limit[0]==0: continue
        index = int(limit[0])
        logger.info('Convection region between %g-%g (%s,color=%s)'%(limit[1],limit[2],
                     ['Stable','Unstable','nan','Ledoux stable','Schwarzschild stable'][index],
                     colors[index]))
        out = pl.axvspan(limit[1],limit[2],color=colors[index],alpha=alphas[index])


def plot_add_burn(starl,x='q'):
    burn = starl['eps_nuc']
    burntype = np.zeros(len(burn))
    burntype[(burn<=1e0)] = 0
    burntype[(1e0<burn)&(burn<=1e3)] = 1
    burntype[(1e3<burn)&(burn<=1e7)] = 2
    burntype[(1e7<burn)] = 3
    x = starl[x]
    #-- search for the limits of convection types
    limits = [[burntype[0],x[0]]]
    for i in range(1,len(x)):
        if burntype[i]!=limits[-1][0]:
            limits[-1].append(x[i])
            limits.append([burntype[i],x[i]])
    limits[-1].append(x[-1])
    alphas = [1.0,0.5,0.5,0.5,0.5]
    colors = ['k',(1.00,1.00,0.0),(1.00,0.50,0.0),(1.00,0.00,0.0)]
    for limit in limits:
        if limit[0]==0: continue
        index = int(limit[0])
        logger.info('Burning between %g-%g (%s,color=%s)'%(limit[1],limit[2],
                     ['Low','Mid','High'][index],
                     colors[index]))
        out = pl.axvspan(limit[1],limit[2],color=colors[index],alpha=alphas[index])

#}
if __name__=="__main__":
    data = read_agsm('tempmodel_frq.gsm')
    ascii.write_array(data,'test.dat',header=True,auto_width=True)