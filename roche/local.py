"""
Derive local quantities such as effective temperature and surface gravity for Roche potentials.
"""
import numpy as np
import pylab as pl
from numpy import pi,cos,sin,sqrt
from scipy.optimize import newton
try:
    from scipy.spatial import Delaunay
except ImportError:
    print('No Delaunay')

from ivs.sed import limbdark
from ivs.coordinates import vectors


#{ Gridding functions

def get_grid(*args,**kwargs):
    """
    Construct a coordinate grid

    If you give two resolutions, the first is for theta, the second for phi

    @param args: one or two integers indicating number of grid points in theta
    and phi direction
    @type args: integer
    @keyword gtype: grid type ('spher' or 'delaunay')
    @type gtype: str
    @return: theta,phi(,grid)
    """
    gtype = kwargs.get('gtype','spherical')
    full = kwargs.get('full',False)

    if 'spher' in gtype.lower():
        if len(args)==1: res1 = res2 = args[0] # same resolution for both coordinates
        else:            res1,res2 = args      # different resolution

        dtheta = pi/res1 # step size coordinate 1
        dphi = pi/res2     # step size coordinate 2

        #-- full grid or only one quadrant
        if full:
            theta0,thetan = dtheta, pi-dtheta
            phi0,phin = 0,2*pi-2*dphi
        else:
            theta0,thetan = dtheta/4,pi/2-dtheta/4
            phi0,phin = 0,pi-dphi

        theta,phi = np.mgrid[theta0:thetan:res1*1j,phi0:phin:res2*1j]
        return theta,phi

    elif 'cil' in gtype.lower():
        if len(args)==1: res1 = res2 = args[0] # same resolution for both coordinates
        else:            res1,res2 = args      # different resolution

        dsintheta = 2./res1 # step size sin(theta)
        dphi = pi/res2     # step size coordinate 2

        #-- full grid or only one quadrant
        if full:
            sintheta0,sinthetan = -1+dsintheta, 1.-dsintheta
            phi0,phin = 0,2*pi-dphi
        else:
            sintheta0,sinthetan = dsintheta,1-dsintheta
            phi0,phin = 0,pi-dphi

        sintheta,phi = np.mgrid[sintheta0:sinthetan:res1*1j,phi0:phin:res2*1j]
        return np.arccos(sintheta),phi

    elif 'delaunay' in gtype.lower():
        #-- resolution doesn't matter that much anymore
        if len(args)==1: res1 = res2 = args[0] # same resolution for both coordinates
        else:            res1,res2 = args      # different resolution

        u,v = np.random.uniform(size=res1*res2),np.random.uniform(size=res1*res2)
        phi = 2*pi*u
        theta = np.arccos(2*v-1)

        x = sin(theta)*sin(phi)
        y = sin(theta)*cos(phi)
        z = cos(theta)

        points = np.array([x,y,z]).T
        grid = Delaunay(points)
        centers = np.zeros((len(grid.convex_hull),3))
        for i,indices in enumerate(grid.convex_hull):
            centers[i] = [x[indices].sum()/3,y[indices].sum()/3,z[indices].sum()/3]
        theta,phi = np.arccos(centers[:,2]),np.arctan2(centers[:,1],centers[:,0])+pi

        return theta,phi,grid

    elif 'tri' in gtype.lower():
        #-- resolution doesn't matter that much anymore
        if len(args)==1: res1 = res2 = args[0] # same resolution for both coordinates
        else:            res1,res2 = args      # different resolution

        u,v = np.random.uniform(size=res1*res2),np.random.uniform(size=res1*res2)
        phi = 2*pi*u
        theta = np.arccos(2*v-1)
        return theta,phi







def stitch_grid(theta,phi,*quant,**kwargs):
    """
    Stitch a grid together that was originally defined on 1 quadrant.

    We add the three other quandrants.
    """
    seamless = kwargs.get('seamless',False)
    vtype = kwargs.get('vtype',['scalar' for i in quant])
    gtype = kwargs.get('gtype','spher')
    ravel = kwargs.get('ravel',False)

    if gtype == 'spher':
        #-- basic coordinates
        alltheta = np.vstack([np.hstack([theta,theta]),np.hstack([theta+pi/2,theta+pi/2])])
        allphi   = np.vstack([np.hstack([phi,phi+pi]),np.hstack([phi,phi+pi])])

        #-- all other quantities
        allquan = []
        for i,iquant in enumerate(quant):
            #-- if they are scalar values, they do not change direction
            top1,top2 = iquant,iquant[:,::-1]
            bot1,bot2 = iquant[::-1],iquant[::-1][:,::-1]
            if vtype[i]=='scalar':
                allquan.append(np.vstack([np.hstack([top1,top2]),np.hstack([bot1,bot2])]))
            #-- vector components do change direction
            elif vtype[i]=='x':
                allquan.append(np.vstack([np.hstack([top1,top2]),np.hstack([bot1,bot2])]))
            elif vtype[i]=='y':
                allquan.append(np.vstack([np.hstack([top1,-top2]),np.hstack([bot1,-bot2])]))
            elif vtype[i]=='z':
                allquan.append(np.vstack([np.hstack([top1,top2]),np.hstack([-bot1,-bot2])]))
            elif vtype[i]=='vx':
                allquan.append(np.vstack([np.hstack([top1,-top2]),np.hstack([bot1,-bot2])]))
            elif vtype[i]=='vy':
                allquan.append(np.vstack([np.hstack([top1,top2]),np.hstack([bot1,bot2])]))
            elif vtype[i]=='vz':
                allquan.append(np.vstack([np.hstack([top1,top2]),np.hstack([-bot1,-bot2])]))

        out = [alltheta,allphi]+allquan

        #-- for plotting reasons, remove the vertical seam and the the bottom hole.
        if seamless:
            #-- vertical seam
            out = [np.column_stack([i,i[:,0]]) for i in out]
            out[1][:,-1] += 2*pi
            #-- fill bottom hole
            out = [np.vstack([i,i[0]]) for i in out]
            out[0][-1] += pi

        #-- ravel to 1d arrays if asked for
        if ravel:
            out = [i.ravel() for i in out]
    else:
        out = [theta,phi] + list(quant)

    return out

#}
def surface_normals(r,phi,theta,grid,gtype='spher'):
    """
    Numerically compute surface normals of a grid (in absence of analytical alternative).

    Also computes the surface elements, making L{surface_elements} obsolete.
    """
    if gtype=='spher':
        raise NotImplementedError
    elif gtype=='delaunay':
        raise NotImplementedError
    elif gtype=='triangular':
        #-- compute the angle between the surface normal and the radius vector
        x,y,z = vectors.spher2cart_coord(r,phi,theta)

        centers = np.zeros((len(grid.convex_hull),3))
        normals = np.zeros((len(grid.convex_hull),3))
        sizes = np.zeros(len(grid.convex_hull))

        #vertx,verty,vertz = points.T

        #-- compute centers,normals and sizes
        for i,indices in enumerate(grid.convex_hull):
            #-- center is triangle's barycenter
            centers[i] = [x[indices].sum()/3,y[indices].sum()/3,z[indices].sum()/3]
            #-- size is size of triangle
            a = sqrt((x[indices[0]]-x[indices[1]])**2 + (y[indices[0]]-y[indices[1]])**2 + (z[indices[0]]-z[indices[1]])**2)
            b = sqrt((x[indices[0]]-x[indices[2]])**2 + (y[indices[0]]-y[indices[2]])**2 + (z[indices[0]]-z[indices[2]])**2)
            c = sqrt((x[indices[1]]-x[indices[2]])**2 + (y[indices[1]]-y[indices[2]])**2 + (z[indices[1]]-z[indices[2]])**2)
            s = 0.5*(a+b+c)
            sizes[i] = sqrt( s*(s-a)*(s-b)*(s-c))
            #-- normal is cross product of two sides
            side1 = [x[indices[1]]-x[indices[0]],y[indices[1]]-y[indices[0]],z[indices[1]]-z[indices[0]]]
            side2 = [x[indices[2]]-x[indices[0]],y[indices[2]]-y[indices[0]],z[indices[2]]-z[indices[0]]]
            normals[i] = np.cross(side1,side2)

        #-- make sure the normal is pointed outwards
        normal_r,normal_phi,normal_theta = vectors.cart2spher(centers.T,normals.T)
        normal_r = np.abs(normal_r)
        centers_sph = vectors.cart2spher_coord(*centers.T)
        normals = np.array(vectors.spher2cart(centers_sph,(normal_r,normal_phi,normal_theta)))

        #-- normalise and compute angles
        normals_T = normals.T
        normals = normals_T / vectors.norm(normals_T)
        #cos_gamma = vectors.cos_angle(a,normals)
        print(centers.shape,sizes.shape,normals.shape)
        return centers, sizes, normals#, cos_gamma

#{ Derivation of local quantities

def surface_elements(radius_and_mygrid, surface_normals_xyz, gtype='spher'):
    """
    Compute surface area of elements in a grid.

    theta,phi must be generated like mgrid(theta_range,phi_range)

    usually, the surfnormals are acquired via differentiation of a gravity potential,
    and is then equal to the *negative* of the local surface gravity.
    """
    (r,mygrid) = radius_and_mygrid
    (surfnormal_x,surfnormal_y,surfnormal_z) = surface_normals_xyz
    theta,phi = mygrid[:2]
    if gtype=='spher':
        #-- compute the grid size at each location
        dtheta = theta[1:]-theta[:-1]
        dtheta = np.vstack([dtheta,dtheta[-1]])

        dphi = phi[:,1:]-phi[:,:-1]
        dphi = np.column_stack([dphi,dphi[:,-1]])
        #-- compute the angle between the surface normal and the radius vector
        x,y,z = vectors.spher2cart_coord(r,phi,theta)

        a = np.array([x,y,z])
        b = np.array([surfnormal_x,surfnormal_y,surfnormal_z])

        cos_gamma = vectors.cos_angle(a,b)

        return r**2 * sin(theta) * dtheta * dphi / cos_gamma, cos_gamma

    elif gtype=='delaunay':
        #-- compute the angle between the surface normal and the radius vector
        x,y,z = vectors.spher2cart_coord(r,phi,theta)
        a = np.array([x,y,z])
        b = np.array([surfnormal_x,surfnormal_y,surfnormal_z])
        cos_gamma = vectors.cos_angle(a,b)

        delaunay_grid = mygrid[2]

        sizes = np.zeros(len(delaunay_grid.convex_hull))
        points = delaunay_grid.points
        vertx,verty,vertz = points.T

        #from enthought.mayavi import mlab
        #mlab.figure()
        #mlab.triangular_mesh(vertx,verty,vertz,delaunay_grid.convex_hull,scalars=np.ones_like(vertx),colormap='gray',representation='wireframe')
        #mlab.points3d(x/r,y/r,z/r,scale_factor=0.02)

        centers = np.zeros((len(delaunay_grid.convex_hull),3))
        for i,indices in enumerate(delaunay_grid.convex_hull):
            #centers[i] = [vertx[indices].sum()/3,verty[indices].sum()/3,vertz[indices].sum()/3]
            a = sqrt((vertx[indices[0]]-vertx[indices[1]])**2 + (verty[indices[0]]-verty[indices[1]])**2 + (vertz[indices[0]]-vertz[indices[1]])**2)
            b = sqrt((vertx[indices[0]]-vertx[indices[2]])**2 + (verty[indices[0]]-verty[indices[2]])**2 + (vertz[indices[0]]-vertz[indices[2]])**2)
            c = sqrt((vertx[indices[1]]-vertx[indices[2]])**2 + (verty[indices[1]]-verty[indices[2]])**2 + (vertz[indices[1]]-vertz[indices[2]])**2)
            s = 0.5*(a+b+c)
            sizes[i] = sqrt( s*(s-a)*(s-b)*(s-c))

        #theta,phi = np.arccos(centers[:,2]),np.arctan2(centers[:,1],centers[:,0])+pi
        #mlab.points3d(centers[:,0],centers[:,1],centers[:,2],sizes,scale_factor=0.05,scale_mode='none',colormap='RdBu')
        #mlab.show()

        #pl.show()

        return sizes*r**2, cos_gamma

def temperature(surface_gravity,g_pole,T_pole,beta=1.):
    """
    Calculate local temperature.

    beta is gravity darkening parameter.
    """
    Grav = abs(surface_gravity/g_pole)**beta
    Teff = Grav**0.25 * T_pole
    return Teff

def intensity(teff,grav,mu=None,photband='OPEN.BOL'):
    """
    Calculate local intensity.

    beta is gravity darkening parameter.
    """
    if mu is None:
        mu = np.ones_like(teff)
    if (teff<3500).any() or np.isnan(teff).any():
        print('WARNING: point outside of grid, minimum temperature is 3500K')
        teff = np.where((teff<3500) | np.isnan(teff),3500,teff)
    if (grav<0.01).any() or np.isnan(grav).any():
        print('WARNING: point outside of grid, minimum gravity is 0 dex')
        grav = np.where((np.log10(grav*100)<0.) | np.isnan(grav),0.01,grav)
    intens = np.array([limbdark.get_itable(teff=iteff,logg=np.log10(igrav*100),absolute=True,mu=imu,photbands=[photband])[0] for iteff,igrav,imu in zip(teff.ravel(),grav.ravel(),mu.ravel())])
    return intens.reshape(teff.shape)


def projected_intensity(teff,gravity,areas,line_of_sight,photband='OPEN.BOL'):
    """
    Compute projected intensity in the line of sight.

    gravity is vector directed inwards in the star
    line of sight is vector.
    """
    ones = np.ones_like(gravity[0])
    losx = line_of_sight[0]*ones
    losy = line_of_sight[1]*ones
    losz = line_of_sight[2]*ones
    angles = vectors.angle(-gravity,np.array([losx,losy,losz]))
    mus = cos(angles)
    grav_ = vectors.norm(gravity)
    #-- intensity is less if we look at the limb
    intens = intensity(teff,grav_,mu=mus,photband=photband)
    #-- intensity is less if the surface element area is small (it does not need
    #   to be projected anymore!)
    return intens*areas*mus,mus


def project(star,view_long=(0,0,0),view_lat=(pi/2,0,0),photband='OPEN.BOL',
            only_visible=False,plot_sort=False,scale_factor=1.):
    """
    Project and transform coordinates and vectors to align with the line-of-sight.

    Parameter C{star} should be a record array containing fields 'teff','gravx',
    'gravy','gravz','areas','vx','vy','vz'

    and either you suply ('r','theta','phi') or ('x','y','z')

    The XY direction is then the line-of-sight, and the YZ plane is the plane
    of the sky.

    An extra column 'projflux' and 'eyeflux' will be added. Projected flux
    takes care of limb darkening, and projected surface area. Eye flux only
    takes care of limbdarkening, and should only be used for plotting reasons.

    view_long[0] of 0 means looking in the XY line, pi means looking in the YX line.
    view_lat[0] of pi/2 means edge on, 0 or pi is pole-on.

    This function updates all Cartesian coordinates present in the star, but not
    the polar coordinates! The projected fluxes are added as a field 'projflux'
    to the returned record array.

    If you set 'only_visible' to True, only the information on the visible parts
    of the star will be contained.

    If you set 'plot_sort' to True, the arrays will be returned in a sorted order,
    where the areas at the back come first. This is especially handy for plotting.

    @parameters star: record array containing all necessary information on the
    star
    @type star: numpy record array
    @parameter view_long: longitude viewing angle (radians) and coordinate zeropoint
    @type view_long: tuple floats  (radians,x,y)
    @parameter view_lat: inclination viewing angle (radians) and coordinate zeropoint
    @type view_lat: tuple floats (radians,x,z)
    @parameter photband: photometric passband
    @type photband: string
    @parameter only_visible: flag to return only information on visible surface elements
    @type only_visible: boolean
    @parameter plot_sort: flag to sort the surface elements from back to front
    @type plot_sort: boolean
    """
    myshape = star['gravx'].shape
    gravx,gravy,gravz = np.array([star['gravx'].ravel(),star['gravy'].ravel(),star['gravz'].ravel()])
    areas = star['areas'].ravel()
    teff = star['teff'].ravel()
    vx,vy,vz = star['vx'].ravel(),star['vy'].ravel(),star['vz'].ravel()
    #-- if 'x' is not in the star's record array, we assume the polar coordinates
    #   are in there and convert them to Cartesian coordinates
    if not 'x' in star.dtype.names:
        x,y,z = vectors.spher2cart_coord(star['r'].ravel(),star['phi'].ravel(),star['theta'].ravel())
    else:
        x,y,z = star['x'].ravel(),star['y'].ravel(),star['z'].ravel(),

    #-- first we rotate in the XY plane (only for surface coordinates is the
    #   coordinate zeropoint important, the rest are vectors!):
    x,y = vectors.rotate(x,y,view_long[0],x0=view_long[1],y0=view_long[2])
    gravx,gravy = vectors.rotate(gravx,gravy,view_long[0])
    vx,vy = vectors.rotate(vx,vy,view_long[0])
    #-- then we rotate in the YZ plane:
    if view_lat[0]!=pi/2:
        rot_i = -(pi/2 - view_lat[0])
        x,z = vectors.rotate(x,z,rot_i)
        gravx,gravz = vectors.rotate(gravx,gravz,rot_i)
        vx,vz = vectors.rotate(vx,vz,rot_i)
    #-- ... and project the fluxes in the line of sight, which is now in the XY
    #   direction:
    view_vector = np.array([1.,0,0])#np.array([-sin(pi/2),0,-cos(pi/2)])
    grav_local = np.array([gravx,gravy,gravz])
    proj_flux,mus = projected_intensity(teff,grav_local,areas,view_vector,photband=photband)

    #-- we now construct a copy of the star record array with the changed
    #   coordinates
    new_star = star.copy()
    new_star['gravx'],new_star['gravy'],new_star['gravz'] = gravx,gravy,gravz
    new_star['vx'],new_star['vy'],new_star['vz'] = vx,vy,vz
    if 'x' in star.dtype.names:
        new_star['x'],new_star['y'],new_star['z'] = x*scale_factor,y*scale_factor,z*scale_factor
    else:
        new_star = pl.mlab.rec_append_fields(new_star,'x',x*scale_factor)
        new_star = pl.mlab.rec_append_fields(new_star,'y',y*scale_factor)
        new_star = pl.mlab.rec_append_fields(new_star,'z',z*scale_factor)
    new_star = pl.mlab.rec_append_fields(new_star,'projflux',proj_flux)
    new_star = pl.mlab.rec_append_fields(new_star,'eyeflux',proj_flux/areas)
    new_star = pl.mlab.rec_append_fields(new_star,'mu',mus)

    #-- clip visible areas and sort in plotting order if necessary
    if only_visible:
        new_star = new_star[-np.isnan(new_star['projflux'])]
    if plot_sort:
        new_star = new_star[np.argsort(new_star['x'])]
    return new_star
