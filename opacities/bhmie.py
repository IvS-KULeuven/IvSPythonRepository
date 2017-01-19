import numpy as np

"""
Calculate extinction, absorption and scattering efficiencies for dust grains.

Currently implemented:
    - Spherical grains
    
"""

def _times(a,b):
    """
    Multiply an array by a matrix
    
    a is a vector of length n
    b is a matrix of dimensions (n,m)
    """
    if not a.shape[0] == b.shape[0] and len(a) == 1 and len(b.shape) == 2.:
        raise ValueError("the dimensions are not wat they should be, check the manual")
       
    out = np.zeros_like(b)
    for i in np.arange(b.shape[1]):
        out[:,i] = a*b[:,i]
    return out    
            

def bhmie_single_point(wavelength,refrel,radius,nang):
    """
    This file is converted from mie.m, see http://atol.ucsd.edu/scatlib/index.htm
    Bohren and Huffman originally published the code in their book on light scattering
    
    Calculation based on Mie scattering theory
    
    @param wavelength: wavelength in micron
    @type newx: float
    @param refrel: refraction index (n in complex form for example:  1.5+0.02*i)
    @type refrel: complex
    @param radius: radius in micron
    @type radius: float
    @param nang: number of angles for S1 and S2 function in range from 0 to pi/2
    @type nang: integer
    
    @return: extinction efficiency, scattering efficiency, backscatter efficiency, asymmetry parameter
    @rtype: float, float, float, ndarray
    """
    x    = radius*2.*np.pi/wavelength
    nmxx = 150000
    
    s1_1 = np.zeros(nang,dtype=np.complex128)
    s1_2 = np.zeros(nang,dtype=np.complex128)
    s2_1 = np.zeros(nang,dtype=np.complex128)
    s2_2 = np.zeros(nang,dtype=np.complex128)
    pi   = np.zeros(nang,dtype=np.complex128)
    tau  = np.zeros(nang,dtype=np.complex128)
    
    if (nang > 1000):
        print ('error: nang > mxnang=1000 in bhmie')
        return
    
    # Require NANG>1 in order to calculate scattering intensities
    if (nang < 2):
        nang = 2
    
    pii = 4.*np.arctan(1.)
    dx = x
    drefrl = refrel
    y      = x*drefrl
    ymod   = abs(y)
    
    #    Series expansion terminated after NSTOP terms
    #    Logarithmic derivatives calculated from NMX on down
    xstop = x + 4.*x**0.3333 + 2.0
    nmx   = max(xstop,ymod) + 15.0
    nmx   = np.fix(nmx)
    nstop = int(xstop)
    
    if (nmx > nmxx):
        raise ValueError("nmx > nmxx=%f for |m|x=%f" %(nmxx, ymod))

    dang = .5*pii/ (nang-1)
    amu = np.arange(0.0,nang,1)
    amu = np.cos(amu*dang)

    pi0 = np.zeros(nang,dtype = np.complex128)
    pi1 = np.ones(nang,dtype  = np.complex128)

    # Logarithmic derivative D(J) calculated by downward recurrence
    # beginning with initial value (0.,0.) at J=NMX
    nn = int(nmx)-1
    d  = np.zeros(nn+1,dtype = np.complex128)
    for n in range(0,nn):
        en = nmx - n
        d[nn-n-1] = (en/y) - (1./ (d[nn-n]+en/y))


    # Riccati-Bessel functions with real argument X calculated by upward recurrence
    psi0 = np.cos(dx)
    psi1 = np.sin(dx)
    chi0 = -np.sin(dx)
    chi1 = np.cos(dx)
    xi1 = psi1-chi1*1j
    qsca = 0.
    gsca = 0.
    p = -1

    for n in range(0,nstop):
        en = n+1.0
        fn = (2.*en+1.)/(en* (en+1.))
        # Calculate psi_n and chi_n
        psi = (2.*en-1.)*psi1/dx - psi0
        chi = (2.*en-1.)*chi1/dx - chi0
        xi = psi-chi*1j
        # Store previous values of AN and BN for use in computation of g=<cos(theta)>
        if (n > 0):
            an1 = an
            bn1 = bn
        # Compute AN and BN:
        an = (d[n]/drefrl+en/dx)*psi - psi1
        an = an/ ((d[n]/drefrl+en/dx)*xi-xi1)
        bn = (drefrl*d[n]+en/dx)*psi - psi1
        bn = bn/ ((drefrl*d[n]+en/dx)*xi-xi1)
        # Augment sums for Qsca and g=<cos(theta)>
        qsca += (2.*en+1.)* (abs(an)**2+abs(bn)**2)
        gsca += ((2.*en+1.)/(en* (en+1.)))*(np.real(an)*np.real(bn)+np.imag(an)*np.imag(bn))
        if (n > 0):
            gsca += ((en-1.)*(en+1.)/en)*(np.real(an1)*np.real(an)+np.imag(an1)*np.imag(an)+np.real(bn1)*np.real(bn)+np.imag(bn1)*np.imag(bn))
        # Now calculate scattering intensity pattern, but first do angles from 0 to 90
        pi=0+pi1    # 0+pi1 because we want a hard copy of the values
        tau=en*amu*pi-(en+1.)*pi0
        s1_1 += fn*(an*pi+bn*tau)
        s2_1 += fn*(an*tau+bn*pi)
        # Now do angles greater than 90 using PI and TAU from angles less than 90.
        # P=1 for N=1,3,...% P=-1 for N=2,4,... remember that we have to reverse
        # the order of the elements of the second part of s1 and s2 after the calculation
        p = -p
        s1_2+= fn*p* (an*pi-bn*tau)
        s2_2+= fn*p* (bn*pi-an*tau)
        psi0 = psi1
        psi1 = psi
        chi0 = chi1
        chi1 = chi
        xi1 = psi1-chi1*1j
        # Compute pi_n for next value of n
        # For each angle J, compute pi_n+1 from PI = pi_n , PI0 = pi_n-1
        pi1 = ((2.*en+1.)*amu*pi- (en+1.)*pi0)/ en
        pi0 = 0+pi   # 0+pi because we want a hard copy of the values

    # Have summed sufficient terms. Now compute QSCA,QEXT,QBACK,and GSCA we have to
    # reverse the order of the elements of the second part of s1 and s2
    s1=np.concatenate((s1_1,s1_2[-2::-1]))
    s2=np.concatenate((s2_1,s2_2[-2::-1]))
    gsca = 2.*gsca/qsca
    qsca = (2./ (dx*dx))*qsca
    qext = (4./ (dx*dx))*np.real(s1[0])

    # more common definition of the backscattering efficiency,
    # so that the backscattering cross section really
    # has dimension of length squared
    qback = 4*(abs(s1[2*nang-2])/dx)**2
    return qext,qsca,qback,gsca

def bhmie_safe(wavelength, refrel, nang=2, radius=0.1):
    """
    Calculate absorption and scattering efficiencies based on Mie scattering theory
    
    This method to calculate bhmie for large arrays is slower, but the fast program
    will crash if the requirement "wavelength > grain size" is not fullfilled.
    
    @param wavelength: wavelength in micron
    @type newx: float
    @param refrel: refraction index (n in complex form for example:  1.5+0.02*i)
    @type refrel: complex
    @param radius: radius in micron
    @type radius: float
    @param nang: number of angles for S1 and S2 function in range from 0 to pi/2
    @type nang: integer
    
    @return: extinction efficiency, scattering efficiency, backscatter efficiency, asymmetry parameter
    @rtype: float, float, float, ndarray
    """
    qext  = []
    qsca  = []
    qback = []
    gsca  = []
    for i in xrange(len(wavelength)):
        qe,qs,qb,gs = bhmie_single_point(wavelength[i],refrel[i],radius,nang)
        qext.append(qe)
        qsca.append(qs)
        qback.append(qb)
        gsca.append(gs)
    return np.array(qext, dtype=float), np.array(qsca, dtype=float), np.array(qback, dtype=float), np.array(gsca, dtype=float)

def bhmie(wavelength, refrel, nang=2, radius=0.1):
    """
    Calculate absorption and scattering efficiencies based on Mie scattering theory
    
    This file is converted from mie.m, see http://atol.ucsd.edu/scatlib/index.htm
    Bohren and Huffman originally published the code in their book on light scattering
    
    Caution: if the requirement "wavelength > grain size", this method will likely
    return NaN values. In this case, we urge the user to use the bhmie_safe method
    
    @param wavelength: wavelength in micron
    @type newx: float
    @param refrel: refraction index (n in complex form for example:  1.5+0.02*i)
    @type refrel: complex
    @param radius: radius in micron
    @type radius: float
    @param nang: number of angles for S1 and S2 function in range from 0 to pi/2
    @type nang: integer
    
    @return: extinction efficiency, scattering efficiency, backscatter efficiency, asymmetry parameter
    @rtype: float, float, float, ndarray
    """
    
    x    = radius*2.*np.pi/wavelength
    nx   = len(x)
    nmxx = 150000

    s1_1 = np.zeros((nx, nang), dtype=np.complex128)
    s1_2 = np.zeros((nx, nang), dtype=np.complex128)
    s2_1 = np.zeros((nx, nang), dtype=np.complex128)
    s2_2 = np.zeros((nx, nang), dtype=np.complex128)
    pi   = np.zeros((nx, nang), dtype=np.complex128)
    tau  = np.zeros((nx, nang), dtype=np.complex128)

    if (nang > 1000):
      raise ValueError('Nang > mxnang=1000 in bhmie')

    # Require NANG>1 in order to calculate scattering intensities
    if (nang < 2):
        raise ValueError("Nang < 2 in bhmie")

    y      = x*refrel
    ymod   = abs(y)
    
    #    Series expansion terminated after NSTOP terms
    #    Logarithmic derivatives calculated from NMX on down
    xstop = x + 4.*x**0.3333 + 2.0
    nmx   = max(xstop.max(),ymod.max()) + 15.0
    nmx   = np.ceil(nmx)

    # BTD experiment 91/1/15: add one more term to series and compare resu<s
    #      NMX=AMAX1(XSTOP,YMOD)+16
    # test: compute 7001 wavelen>hs between .0001 and 1000 micron
    # for a=1.0micron SiC grain.  When NMX increased by 1, only a single
    # computed number changed (out of 4*7001) and it only changed by 1/8387
    # conclusion: we are indeed retaining enough terms in series!

    nstop = int(xstop.max())

    if (nmx > nmxx):
        raise ValueError("Nmx > nmxx=%f for |m|x=%f" % (nmxx, ymod))

    amu = np.cos(np.linspace(0., np.pi/2., nang))

    pi0 = np.zeros((nx, nang), dtype=np.complex128)
    pi1 = np.ones((nx, nang), dtype=np.complex128)

    # Logarithmic derivative D(J) calculated by downward recurrence
    # beginning with initial value (0.,0.) at J=NMX

    nn = int(nmx-1)
    d  = np.zeros((nx,nn+1),dtype=np.complex128)
    for n in range(0,nn):
        en = nmx - n
        d[:,nn-n-1] = (en/y) - (1./ (d[:,nn-n]+en/y))


    #*** Riccati-Bessel functions with real argument X
    #    calculated by upward recurrence

    psi0 = np.cos(x)
    psi1 = np.sin(x)
    chi0 = -np.sin(x)
    chi1 = np.cos(x)
    xi1  = psi1-chi1*complex(0,1)
    qsca = 0.
    gsca = 0.
    p    = -1

    for n in range(0,nstop):
        en = n+1.0
        fn = (2.*en+1.)/(en* (en+1.))
        # for given N, PSI  = psi_n        CHI  = chi_n
        #              PSI1 = psi_{n-1}    CHI1 = chi_{n-1}
        #              PSI0 = psi_{n-2}    CHI0 = chi_{n-2}
        # Calculate psi_n and chi_n
        psi = (2.*en-1.)*psi1/x - psi0
        chi = (2.*en-1.)*chi1/x - chi0
        xi = psi-chi*complex(0,1)
        
        # Store previous values of AN and BN for use in computation of g=<cos(theta)>
        if (n > 0):
            an1 = an
            bn1 = bn
            
        # Compute AN and BN:
        an = (d[:,n]/refrel+en/x)*psi - psi1
        an = an/ ((d[:,n]/refrel+en/x)*xi-xi1)
        bn = (refrel*d[:,n]+en/x)*psi - psi1
        bn = bn/((refrel*d[:,n]+en/x)*xi-xi1)
        
        # Augment sums for Qsca and g=<cos(theta)>
        qsca += (2.*en+1.)* (abs(an)**2+abs(bn)**2)
        gsca += ((2.*en+1.)/ (en*(en+1.)))*(np.real(an)*np.real(bn)+np.imag(an)*np.imag(bn))
        if (n > 0):
            gsca += ((en-1.)* (en+1.)/en)*(np.real(an1)*np.real(an)+np.imag(an1)*np.imag(an)+np.real(bn1)*np.real(bn)+np.imag(bn1)*np.imag(bn))
        

        
        # Now calculate scattering intensity pattern
        #    First do angles from 0 to 90
        pipi  = 0+pi1    # 0+pi1 because we want a hard copy of the values
        tau   = en*amu*pipi-(en+1.)*pi0
        s1_1 += fn*(_times(an,pipi)+_times(bn,tau))
        s2_1 += fn*(_times(an,tau)+_times(bn,pipi))
        
        # Now do angles greater than 90 using PI and TAU from angles less than 90.
        # P=1 for N=1,3,...% P=-1 for N=2,4,...
        # remember that we have to reverse the order of the elements
        # of the second part of s1 and s2 after the calculation
        p = -p
        s1_2+= fn*p*(_times(an,pipi)-_times(bn,tau))
        s2_2+= fn*p*(_times(bn,pipi)-_times(an,tau))
        
        psi0 = psi1
        psi1 = psi
        chi0 = chi1
        chi1 = chi
        xi1 = psi1-chi1*complex(0,1)
        
        # Compute pi_n for next value of n
        # For each angle J, compute pi_n+1 from PI = pi_n , PI0 = pi_n-1
        pi1 = ((2.*en+1.)*amu*pi- (en+1.)*pi0)/ en
        pi0 = 0+pipi   # 0+pi because we want a hard copy of the values

    #*** Have summed sufficient terms.
    #    Now compute QSCA,QEXT,QBACK,and GSCA

    #   we have to reverse the order of the elements of the second part of s1 and s2
    
    s1=np.concatenate((s1_1,s1_2[:,-2::-1]), axis=1)
    s2=np.concatenate((s2_1,s2_2[:,-2::-1]), axis=1)
    gsca = 2.*gsca/qsca
    qsca = (2./(x*x))*qsca
    qext = (4./(x*x))*np.real(s1[:,0])

    # more common definition of the backscattering efficiency,
    # so that the backscattering cross section really
    # has dimension of length squared
    qback = 4*(abs(s1[:,2*nang-2])/x)**2

    return np.array(qext, dtype=float), np.array(qsca, dtype=float), np.array(qback, dtype=float), np.array(gsca, dtype=float)

if __name__ == '__main__':
    from ivs.inout.ascii import read2array
    lnk        = read2array("/home/kristofs/python/IVSdata/optical_constants/alumino-silicates/A/ca2al2sio7_02.91_0000_001.lnk")
    wavelength = lnk[:,0]
    refrel     = lnk[:,1] + lnk[:,2]*complex(0,1)
    radius     = .01
    nang       = 5.
    qext1, qsca1, qback1, gsca1 = bhmie(wavelength, refrel, nang, radius)
    qext2, qsca2, qback2, gsca2 = bhmie_slow(wavelength, refrel, nang, radius)
