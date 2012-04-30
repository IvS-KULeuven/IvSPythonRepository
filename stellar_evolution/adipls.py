import logging
import numpy as np
from ivs.units import constants
from ivs.units import conversions
from ivs.aux import numpy_ext as ne

logger = logging.getLogger('STAR.ADIPLS')

    
def make_adipls_inputfile(model_name='model',degree=0,f0=None,fn=None,nf=1000,
                          fspacing='square',**kwargs):
    """
    Prepare an input file for ADIPLS.
    
    fspacing is uniform in 'square'd frequency, 'frequency' or 'period'.
    """
    global_params = dict(input_file=model_name+'.amdl',model_name=model_name)
    #-- parameter names
    par_names = {}
    #-- equilibirum model controls (mod)
    par_names['mod'] = [['ifind','xmod','imlds','in','nprmod'],
                        ['xtrnct','ntrnsf','imdmod']]
    #-- controls for mode selection
    par_names['osc'] = [['el','nsel','els1','dels','dfsig1','dfsig2','nsig1','nsig2'],
                        ['itrsig','sig1','istsig','inomde','itrds'],
                        ['dfsig','nsig','iscan','sig2'],
                        ['eltrw1','eltrw2','sgtrw1','sgtrw2']]
    #-- controls for rotational solution
    #par_names['rot'] = ['irotsl','em','nsem','ems1','dems']
    #-- controls for fundamental constants
    par_names['cst'] = [['cgrav']]
    #-- controls for equations and integration
    par_names['int'] = [['iplneq','iturpr','icow','alb'],
                        ['istsbc','fctsbc','ibotbc','fcttbc'],
                        ['mdintg','iriche','xfit','fcnorm','eps','epssol','itmax','dsigre'],
                        ['fsig','dsigmx','irsevn','xmnevn','nftmax','itsord']]
    #-- controls for output
    par_names['out'] = [['istdpr','nout','nprcen','irsord','iekinr'],
                        ['iper','ivarf','kvarf','npvarf','nfmode'],
                        ['irotkr','nprtkr','igm1kr','npgmkr','ispcpr'],
                        ['icaswn','sigwn1','sigwn2','frqwn1','frqwn2','iorwn1','iorwn2','frlwn1','frlwn2']]
    par_names['dgn'] = [['itssol','idgtss','moddet','iprdet','npout'],
                        ['imstsl','imissl','imjssl','idgnrk']]
    
    #-- set degree of the mode
    if not hasattr(degree,'__iter__'):
        kwargs['el'] = degree
        kwargs['nsel'] = 0
    #-- set frequency range
    #   if f0 and fn are not given as a tuple, they are dimensionless
    if f0 is not None or fn is not None:
        kwargs['istsig'] = 2 # sig1 and sig2 are limits in cyclic frequency nu
        if f0 is not None: kwargs['sig1'] = conversions.convert(f0[1],'mHz',f0[0])
        if fn is not None: kwargs['sig2'] = conversions.convert(fn[1],'mHz',fn[0])
    if f0 is not None and fn is not None:
        #-- frequency range
         kwargs['itrsig'] = 1 # check for change of sign in the matching
                         # determinant and iterate for eigenfrequency
                         # at each change of sign
         #-- set density of trial frequency range
         kwargs['iscan'] = nf
    elif f0 is not None and fn is None:
        #-- explicit setting of frequency, e.g. when computing eigenfrequencies
        kwargs['itrsig'] = 0 # trial frequency taken from sig1
        raise ValueError,'istsig=2 not allowed with itrsig=0'
    
    #-- set frequency step
    if 's' in fspacing.lower():
        kwargs['nsig'] = 1 # square in frequency for low order p modes
    elif 'p' in fspacing.lower():
        kwargs['nsig'] = 2 # linear in frequency for p modes
    elif 'g' in fspacing.lower():
        kwargs['nsig'] = 3 # linear in period for g modes
    else:
        raise ValueError,"Cannot interpret fspacing value"
    
    #-- check for nans
    for sec in par_names.keys():
        for line in par_names[sec]:
            for par in line:
                if np.isnan(kwargs.get(par,0)):
                    raise ValueError,'%s:%s is nan'%(sec,par)
    #-- write the parameter control input file
    control_file = model_name+'.in'
    #os.system('cp tempmodel.amdlr %s'%(model_name+'.amdl')) # historical debugging
    #a,b = special.read_amdl('tempmodel.amdlr') # historical debugging
    #special.write_amdl(a,b,model_name+'.amdl') # historical debugging
    
    with open(control_file,'w') as ff:
        #-- first write global parameters
        ff.write("2 '{0}'    @\n".format(model_name+'.amdl'))
        #ff.write("2 '{0}'    @\n".format('tempmodel.amdlr')) # historical debugging
        ff.write("9 '0'    @\n")
        ff.write("11 '{0}'    @\n".format(model_name+'.gsm'))
        ff.write("15 '{0}'    @\n".format(model_name+'.ssm'))
        if 'nfmode' in kwargs and kwargs['nfmode']>0:
            ff.write("4 '{0}'    @\n".format(model_name+'.aeig'))
            
        ff.write("-1 ''    @\n")
        #-- then specific parameters
        #   the sections have to be given in a specific order
        sections = [sec for sec in ['mod','osc','cst','int','out','dgn'] if sec in par_names.keys()]
        
        ff.write('\n')
        ff.write('cntrd,\n%s     @\n'%('.'.join(sections)))
        for section in sections:
            ff.write('%s:\n'%(section))
            #-- write line per line:
            for line in par_names[section]:
                #   first the names of the parameters for clarity
                ff.write(','.join(line)+'\n')
                #   then the values
                values = [str(kwargs.get(par_name,'')) for par_name in line]
                ff.write(','.join(values)+',     @\n')
    #-- that's it!  
    return control_file






def make_redistrb_inputfile(model_name='model',nn=5000,mode_type='gp',**kwargs):
    """
    Prepare an input file for ADIPLS's redistrb tool.
    """
    global_params = dict(input_file=model_name+'.amdl',model_name=model_name)
    #-- parameter names
    par_names = [['nn','icnmsh'],
                 ['icase','icvzbn','nsmth','ndisc','dlxdsc','dlgrmx','cacvzb'],
                 ['cg','cx','ca','cdgr','cddgr','cdg','alphsf','adda4','accrmn'],
                 ['nout','cn','irsu','unew','iresa4'],
                 ['nmodel','kmodel','itsaml','ioldex','idsin','idsout']]
    icase = 0
    if 'g' in mode_type.lower():
        icase += 12
    if 'p' in mode_type.lower():
        icase += 11
    kwargs['icase'] = icase
    kwargs['nn'] = int(nn)
    
    #-- write the parameter control input file
    control_file = model_name+'.rin'
    with open(control_file,'w') as ff:
        #-- first write global parameters
        ff.write("2 '{0}'    @\n".format(model_name+'.amdl'))
        ff.write("3 '{0}'    @\n".format(model_name+'.amdlr'))
        ff.write("-1 ''    @\n")
        #-- then specific parameters
        #   the sections have to be given in a specific order
        
        ff.write('\n')
        #-- write line per line:
        for line in par_names:
            #   first the names of the parameters for clarity
            ff.write(','.join(line)+'\n')
            #   then the values
            values = [str(kwargs.get(par_name,'')) for par_name in line]
            ff.write(','.join(values)+',     @\n')
    #-- that's it!  
    return control_file    
    