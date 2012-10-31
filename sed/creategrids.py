"""
Create photometric and spectral grid.
"""
import pyfits
import logging
import numpy as np
from multiprocessing import cpu_count,Manager,Process
import os
import shutil
from ivs.sed import model
from ivs.sed import filters
from ivs.sed import reddening
from ivs.io import ascii

logger = logging.getLogger('IVS.SED.CREATE')


def calc_integrated_grid(threads=1,ebvs=None,law='fitzpatrick2004',Rv=3.1,
           units='Flambda',responses=None,update=False,add_spectrophotometry=False,**kwargs):
    """
    Integrate an entire SED grid over all passbands and save to a FITS file.
    
    The output file can be used to fit SEDs more efficiently, since integration
    over the passbands has already been carried out.
    
    WARNING: this function can take a day to complete on a 8-core processor!
    
    Extra keywords can be used to specify the grid.
    
    @param threads: number of threads
    @type threads; integer, 'max', 'half' or 'safe' 
    @param ebvs: reddening parameters to include
    @type ebvs: numpy array
    @param law: interstellar reddening law to use
    @type law: string (valid law name, see C{reddening.py})
    @param Rv: Rv value for reddening law
    @type Rv: float
    @param units: choose to work in 'Flambda' or 'Fnu'
    @type units: str, one of 'Flambda','Fnu'
    @param responses: respons curves to add (if None, add all)
    @type responses: list of strings
    @param update: if true append to existing FITS file, otherwise overwrite
    possible existing file.
    @type update: boolean
    """    
    if ebvs is None:
        ebvs = np.r_[0:4.01:0.01]
        
    #-- select number of threads
    if threads=='max':
        threads = cpu_count()
    elif threads=='half':
        threads = cpu_count()/2
    elif threads=='safe':
        threads = cpu_count()-1
    threads = int(threads)
    
    if threads > len(ebvs):
        threads = len(ebvs)
    logger.info('Threads: %s'%(threads))
    
    #-- set the parameters for the SED grid
    model.set_defaults(**kwargs)
    #-- get the dimensions of the grid: both the grid points, but also
    #   the wavelength range
    teffs,loggs = model.get_grid_dimensions()
    wave,flux = model.get_table(teff=teffs[0],logg=loggs[0])
    #-- get the response functions covering the wavelength range of the models
    #   also get the information on those filters
    if add_spectrophotometry:
        bands = filters.add_spectrophotometric_filters()
    if responses is None:
        responses = filters.list_response(wave_range=(wave[0],wave[-1]))
    else:
        responses_ = []
        if not any(['BOXCAR' in i for i in responses]) and add_spectrophotometry:
            responses.append('BOXCAR')
        for resp in responses:
            responses_ += filters.list_response(resp)
        responses = responses_
    filter_info = filters.get_info(responses)
    responses = filter_info['photband']
    responses = [resp for resp in responses if not (('ACS' in resp) or ('WFPC' in resp) or ('STIS' in resp) or ('ISOCAM' in resp) or ('NICMOS' in resp))]

    logger.info('Integrating {0} filters'.format(len(responses)))
    
    def do_ebv_process(ebvs,arr,responses):
        logger.debug('EBV: %s-->%s (%d)'%(ebvs[0],ebvs[-1],len(ebvs)))
        for ebv in ebvs:
            flux_ = reddening.redden(flux,wave=wave,ebv=ebv,rtype='flux',law=law,Rv=Rv)
            #-- calculate synthetic fluxes
            synflux = model.synthetic_flux(wave,flux_,responses,units=units)
            arr.append([np.concatenate(([ebv],synflux))])
        logger.debug("Finished EBV process (len(arr)=%d)"%(len(arr)))
    
    #-- do the calculations
    c0 = time.time()
    output = np.zeros((len(teffs)*len(ebvs),4+len(responses)))
    start = 0
    logger.info('Total number of tables: %i'%(len(teffs)))
    exceptions = 0
    exceptions_logs = []
    for i,(teff,logg) in enumerate(zip(teffs,loggs)):
        if i>0:
            logger.info('%s %s %s %s: ET %d seconds'%(teff,logg,i,len(teffs),(time.time()-c0)/i*(len(teffs)-i)))
        
        #-- get model SED and absolute luminosity
        wave,flux = model.get_table(teff,logg)
        Labs = model.luminosity(wave,flux)
        
        #-- threaded calculation over all E(B-V)s
        processes = []
        manager = Manager()
        arr = manager.list([])
        all_processes = []
        for j in range(threads):
            all_processes.append(Process(target=do_ebv_process,args=(ebvs[j::threads],arr,responses)))
            all_processes[-1].start()
        for p in all_processes:
            p.join()
        
        try:
            #-- collect the results and add them to 'output'
            arr = np.vstack([row for row in arr])
            sa = np.argsort(arr[:,0])
            arr = arr[sa]
            output[start:start+arr.shape[0],:3] = teff,logg,Labs
            output[start:start+arr.shape[0],3:] = arr
            start += arr.shape[0]
        except:
            logger.warning('Exception in calculating Teff=%f, logg=%f'%(teff,logg))
            logger.debug('Exception: %s'%(sys.exc_info()[1]))
            exceptions = exceptions + 1
            exceptions_logs.append(sys.exc_info()[1])
    
    #-- make FITS columns
    gridfile = get_file()
    if os.path.isfile(os.path.basename(gridfile)):
        outfile = os.path.basename(gridfile)
    else:
        outfile = os.path.join(os.path.dirname(gridfile),'i{0}'.format(os.path.basename(gridfile)))
    outfile = 'i{0}'.format(os.path.basename(gridfile))
    outfile = os.path.splitext(outfile)
    outfile = outfile[0]+'_law{0}_Rv{1:.2f}'.format(law,Rv)+outfile[1]
    logger.info('Precaution: making original grid backup at {0}.backup'.format(outfile))
    if os.path.isfile(outfile):
        shutil.copy(outfile,outfile+'.backup')
    output = output.T
    if not update or not os.path.isfile(outfile):
        cols = [pyfits.Column(name='teff',format='E',array=output[0]),
                pyfits.Column(name='logg',format='E',array=output[1]),
                pyfits.Column(name='ebv',format='E',array=output[3]),
                pyfits.Column(name='Labs',format='E',array=output[2])]
        for i,photband in enumerate(responses):
            cols.append(pyfits.Column(name=photband,format='E',array=output[4+i]))
    #-- make FITS columns but copy the existing ones
    else:
        hdulist = pyfits.open(outfile,mode='update')
        names = hdulist[1].columns.names
        cols = [pyfits.Column(name=name,format='E',array=hdulist[1].data.field(name)) for name in names]
        for i,photband in enumerate(responses):
            cols.append(pyfits.Column(name=photband,format='E',array=output[4+i]))
        
    #-- make FITS extension and write grid/reddening specifications to header
    table = pyfits.new_table(pyfits.ColDefs(cols))
    table.header.update('gridfile',os.path.basename(gridfile))
    for key in sorted(defaults.keys()):
        key_ = (len(key)>8) and 'HIERARCH '+key or key
        table.header.update(key_,defaults[key])
    for key in sorted(kwargs.keys()):
        key_ = (len(key)>8) and 'HIERARCH '+key or key
        table.header.update(key_,kwargs[key])
    table.header.update('FLUXTYPE',units)
    table.header.update('REDLAW',law,'interstellar reddening law')
    table.header.update('RV',Rv,'interstellar reddening parameter')
    
    #-- make/update complete FITS file
    if not update or not os.path.isfile(outfile):
        if os.path.isfile(outfile):
            os.remove(outfile)
            logger.warning('Removed existing file: %s'%(outfile))
        hdulist = pyfits.HDUList([])
        hdulist.append(pyfits.PrimaryHDU(np.array([[0,0]])))
        hdulist.append(table)
        hdulist.writeto(outfile)
        logger.info("Written output to %s"%(outfile))
    else:
        hdulist[1] = table
        hdulist.flush()
        hdulist.close()
        logger.info("Appended output to %s"%(outfile))
    
    logger.warning('Encountered %s exceptions!'%(exceptions))
    for i in exceptions_logs:
        print 'ERROR'
        print i

def update_grid(gridfile,responses,threads=10):
    shutil.copy(gridfile,gridfile+'.backup')
    hdulist = pyfits.open(gridfile,mode='update')
    existing_responses = set(list(hdulist[1].columns.names))
    responses = sorted(list(set(responses) - existing_responses))
    if not len(responses):
        hdulist.close()
        print "No new responses to do"
        return None
    law = hdulist[1].header['REDLAW']
    units = hdulist[1].header['FLUXTYPE']
    teffs = hdulist[1].data.field('teff')
    loggs = hdulist[1].data.field('logg')
    ebvs = hdulist[1].data.field('ebv')
    zs = hdulist[1].data.field('z')
    rvs = hdulist[1].data.field('rv')
    vrads = hdulist[1].data.field('vrad')
    names = hdulist[1].columns.names
    
    N = len(teffs)
    index = np.arange(N)
    
    output = np.zeros((len(responses),len(teffs)))
    print N
    
    #--- PARALLEL PROCESS
    def do_process(teffs,loggs,ebvs,zs,rvs,index,arr):
        output = np.zeros((len(responses)+1,len(teffs)))
        c0 = time.time()
        N = len(teffs)
        for i,(teff,logg,ebv,z,rv,ind) in enumerate(zip(teffs,loggs,ebvs,zs,rvs,index)):
            if i%100==0:
                dt = time.time()-c0
                print "ETA",index[0],(N-i)/100.*dt/3600.,'hr'
                c0 = time.time()
            #-- get model SED and absolute luminosity
            model.set_defaults(z=z)
            wave,flux = model.get_table(teff,logg)
            Labs = model.luminosity(wave,flux)
            flux_ = reddening.redden(flux,wave=wave,ebv=ebv,rtype='flux',law=law,Rv=rv)
            #-- calculate synthetic fluxes
            output[0,i] = ind
            output[1:,i] = model.synthetic_flux(wave,flux_,responses,units=units)
        arr.append(output)
    #--- PARALLEL PROCESS
    c0 = time.time()
    
    manager = Manager()
    arr = manager.list([])
    
    all_processes = []
    for j in range(threads):
        all_processes.append(Process(target=do_process,args=(teffs[j::threads],\
                                                                loggs[j::threads],\
                                                                ebvs[j::threads],\
                                                                zs[j::threads],\
                                                                rvs[j::threads],\
                                                                index[j::threads],arr)))
        all_processes[-1].start()
    for p in all_processes:
        p.join()
    
    output = np.hstack([res for res in arr])
    del arr
    sa = np.argsort(output[0])
    output = output[:,sa][1:]
    ascii.write_array(np.rec.fromarrays(output,names=responses),'test.temp',header=True)
    #-- copy old columns and append new ones
    cols = []
    for i,photband in enumerate(responses):
        cols.append(pyfits.Column(name=photband,format='E',array=output[i]))
    #-- create new table
    table = pyfits.new_table(pyfits.ColDefs(cols))
    table = pyfits.new_table(hdulist[1].columns + table.columns,header=hdulist[1].header)
    hdulist[1] = table
    hdulist.close()


def fix_grid(grid):
    hdulist = pyfits.open(grid,mode='update')
    names = hdulist[1].columns.names
    for i,name in enumerate(names):
        if name.lower() in ['teff','logg','ebv','labs','vrad','rv','z']:
            names[i] = name.lower()
    cols = [pyfits.Column(name=name,format='E',array=hdulist[1].data.field(name)) for name in names]
    N = len(hdulist[1].data)

    keys = [key.lower() for key in hdulist[1].header.keys()]

    if not 'z' in names:
        z = hdulist[1].header['z']
        logger.info('Adding metallicity from header {}'.format(z))
        cols.append(pyfits.Column(name='z',format='E',array=np.ones(N)*z))
    else:
        logger.info("Metallicity already in there")
    if not 'vrad' in names:
        vrad = 0.
        logger.info('Adding radial velocity {}'.format(vrad))
        cols.append(pyfits.Column(name='vrad',format='E',array=np.ones(N)*vrad))
    else:
        logger.info("Radial velocity already in there")
        
    fix_rv = False
    if not 'rv' in names:
        if 'rv' in keys:
            rv = hdulist[1].header['Rv']
            logger.info("Adding interstellar Rv from header {}".format(rv))
        else:
            rv = 3.1
            logger.info("Adding default interstellar Rv {}".format(rv))
        cols.append(pyfits.Column(name='rv',format='E',array=np.ones(N)*rv))
    elif not hdulist[1].header['Rv']==hdulist[1].data.field('rv')[0]:
        rv = hdulist[1].header['Rv']
        fix_rv = rv
        logger.info('Correcting interstellar Rv with {}'.format(rv))
    else:
       logger.info("Interstellar Rv already in there")
        
    table = pyfits.new_table(pyfits.ColDefs(cols))
    if fix_rv:
        table.data.field('rv')[:] = rv
    fake_keys = [key.lower() for key in table.header.keys()]
    for key in hdulist[1].header.keys():
        if not key.lower() in fake_keys:
            if len(key)>8:
                key = 'HIERARCH '+key
            table.header.update(key,hdulist[1].header[key])
    hdulist[1] = table
    print "Axis:"
    for name in hdulist[1].columns.names:
        if name.islower() and not name=='labs':
            ax = np.unique(hdulist[1].data.field(name))
            print name,len(ax),min(ax),max(ax)

    teffs = hdulist[1].data.field('teff')
    loggs = hdulist[1].data.field('logg')
    ebvs  = hdulist[1].data.field('ebv')

    #for logg in np.unique(loggs):
        #keep = loggs==logg
        #plt.figure()
        #plt.title(logg)
        #plt.plot(teffs[keep],ebvs[keep],'ko',ms=2)
        
    keep = hdulist[1].data.field('teff')>0
    logger.info('Removing {}/{} false entries'.format(sum(-keep),len(keep)))
    #print len(ff[1].data),sum(keep)
    hdulist[1].data = hdulist[1].data[keep]
    hdulist.close()