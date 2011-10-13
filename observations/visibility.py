"""
Calculate the visibility of an object in the sky.

Examples:

Plot the visibility of an object for 10 days, every 10 minutes.

>>> eph = Ephemeris()
>>> eph.set_objects(objects=['HD50230'])
>>> eph.set_date(startdate='2011/09/27 12:00:00.0',dt=10.,days=10)
>>> eph.set_site(sitename='lapalma')
>>> eph.visibility()
>>> eph.plot(fmt='bo')

You can get the same behaviour in less commands:

>>> eph = Ephemeris(objects=['HD50230'],startdate='2011/09/27 12:00:00.0',dt=10.,days=10,sitename='lapalma')
>>> eph.visibility()

Or, equivalently,

>>> eph = Ephemeris()
>>> eph.visibility(objects=['HD50230'],startdate='2011/09/27 12:00:00.0',dt=10.,days=10,sitename='lapalma')

Plot the visibility of an object only for today:

>>> p = pl.figure()
>>> eph.visibility(days=1)
>>> eph.plot(fmt='y-',lw=2)

"""
import pylab as pl
from matplotlib.dates import date2num,DayLocator,HourLocator,DateFormatter
import numpy as np
from numpy import sqrt,pi,cos,sin
import sys
import time
import pytz
import logging
import ephem

from ivs.units import conversions
from ivs.catalogs import sesame
from ivs.observations import airmass as obs_airmass
from ivs.aux import decorators

logger = logging.getLogger("OBS.VIS")


class Ephemeris(object):
    def __init__(self,**kwargs):
        """
        Initialize an Ephemeris object.
        
        You can set the site, date and objects via additional keyword arguments during initialization.
        You can also give the same information when calling 'visibility'.
        """
        self.set_site(**kwargs)
        self.set_date(**kwargs)
        self.set_objects(**kwargs)
    
    def __get_object(self,objectname):
        """
        Retrieve coordinates for an object in 'ephem'-style.
        """
        jpos = sesame.search(objectname,db='S')['jpos']
        logger.info('Found object %s at %s'%(objectname,jpos))
        myobject = ephem.readdb("%s,f|M|A0,%s,8.0,2000"%(objectname,','.join(jpos.split())))
        return myobject
    
    #{ Set object, site and dates
    @decorators.filter_kwargs
    def set_objects(self,objects=None,**kwargs):
        if objects is not None:
            objects = [self.__get_object(name) for name in objects]
            self.objects = objects
    
    @decorators.filter_kwargs
    def set_site(self,sitename='lapalma',sitelat=None,sitelong=None,siteelev=None):
        """
        Set the observing site.
        
        Supported are: 'lapalma', 'lasilla' and a bunch of capitals from pyephem.
        
        @keyword sitename: the name of the site or city
        @keyword sitelong: longitude (EAST) of site
        @type sitelong: string
        @keyword sitelat: latitude (NORTH) of site
        @type sitelat: string
        @keyword siteelev: site elevation above sea level (m)
        @type siteelev: integer or float
        """
        if sitename == 'lapalma':
            #-- La Palma
            mysite = ephem.Observer()
            #mysite.long = '-17:%.3f'%(52+42/60.) # 17:52:42 WEST
            mysite.long = '342.1184' # 342.1184 EAST
            mysite.lat  = '28:%.3f'%(45+44/60.) # 28:45:44 NORTH
            mysite.elevation = 2326
            mysite.name = sitename
        elif sitename == 'lasilla':
            mysite = ephem.Observer()
            mysite.long = '-70:43.8' # WEST
            mysite.lat  = '-29:15.4' # SOUTH
            mysite.elevation = 2347
            mysite.name = sitename
        elif sitename == 'palomar':
            mysite = ephem.Observer()
            mysite.long = '-116:51:46.80' # WEST
            mysite.lat  = '33:21:21.6' # NORTH
            mysite.elevation = 1706
            mysite.name = sitename
        else:
            try:
                mysite = ephem.city(sitename)
            except:
                logger.critical('Site %s not known, using user definitions'%(sitename))
                mysite.name = sitename
                mysite.elevation = siteelev
                mysite.lat = sitelat
                mysite.long = sitelong
        self.thesite = mysite
        logger.info('Set site to %s (long=%s EAST, lat=%s NORTH, elev=%s m)'%(mysite.name,mysite.long,mysite.lat,mysite.elevation))
    
    @decorators.filter_kwargs
    def set_date(self,startdate=None,days=20,dt=10.):
        """
        Set the start date for ephemeris calculations.
        
        Format 1: Calendar date: '2010/08/16 12:00:00.0'
        Format 2: Julian Date:  2455832.25
        Format 3: None: current time
        
        Careful: all times are UTC, they do not take care of time zones.
        default values are also set in UTC (gmtime)
        """
        if startdate is None:
            startdate = time.strftime('%Y/%m/%d %H:%M:%S',time.gmtime())
        elif not isinstance(startdate,str):
            YYYY,MM,DD = conversions.convert('JD','CD',startdate)
            DD,frac = np.floor(DD),DD-np.floor(DD)
            hh = frac*24.
            hh,frac = np.floor(hh),hh-np.floor(hh)
            mm = frac*60.
            mm,frac = np.floor(mm),mm - np.floor(mm)
            ss = frac*60.
            ss,frac = np.floor(ss),ss - np.floor(ss)
            startdate = '%4d/%02d/%02d %02d:%02d:%.1f'%(YYYY,MM,DD,hh,mm,ss)
        self.thesite.date = ephem.Date(startdate)
        self.dt = dt
        self.days = days
        logger.info('Set start date to %s, covering %d days per %d minutes'%(startdate,days,dt))
    
    #}
    
    #{ Compute visibilities
    
    def visibility(self,midnight=None,**kwargs):
        """
        Calculate ephemeris.
        
        For the list of objects, this function returns:
            - altitudes
            - hours elapsed since startdate
            - sunrises
            - sunsets
            - night time / day time (boolean 1-0)
        """
        #-- set optional additional keywords
        self.set_site(**kwargs)
        self.set_date(**kwargs)
        self.set_objects(**kwargs)
        
        sun = ephem.Sun()
        #moon = ephem.Moon()
        #moon_theta = 34.1 # minutes
        
        #-- set stuff for calculations
        timestep_minutes = self.dt
        total_days = self.days
        if midnight is None:
            total_minutes = int(total_days*24*60./timestep_minutes)
        else:
            total_minutes = 365
        
        #-- set initial arrays
        alts = np.zeros((len(self.objects),total_minutes))
        hours = np.zeros(total_minutes)
        during_night = np.zeros(total_minutes)
        #moon_separation = zeros_like(alts)
        
        #-- run over all timesteps
        if midnight is None:
            for i in range(total_minutes):
                self.thesite.date += ephem.minute*timestep_minutes
                prev_set = float(self.thesite.previous_setting(sun))
                prev_rise = float(self.thesite.previous_rising(sun))
                hours[i] = float(self.thesite.date)
                if (prev_rise<=prev_set):
                    during_night[i] = 1
                for j,star in enumerate(self.objects):               
                    star.compute(self.thesite)
                    alts[j,i] = float(star.alt)
                    
        else:
            i = 0
            while i<365:
                next_rise = float(self.thesite.next_rising(sun))
                prev_set = float(self.thesite.previous_setting(sun))
                prev_rise = float(self.thesite.previous_rising(sun))
                if prev_rise<=prev_set: # now we're in a night
                    self.thesite.date = ephem.Date((next_rise+prev_set)/2.)
                else:
                    #-- set 4 hours forwards
                    self.thesite.date += ephem.minute*60*4
                    continue
                hours[i] = float(self.thesite.date)
                during_night[i] = 1
                for j,star in enumerate(self.stars):               
                    star.compute(self.thesite)
                    alts[j,i] = float(star.alt)
                i += 1
                #-- set 1 day forwards
                self.thesite.date += ephem.minute*60*24
                    
        #-- calculate airmass
        airmass = obs_airmass.airmass(90-alts/pi*180)
        #-- calculate dates for plotting reasons
        dates = np.array([date2num(ephem.date(h).datetime()) for h in hours])
        self.vis = dict(hours=hours,dates=dates,alts=alts,airmass=airmass,during_night=during_night)
        for i,obj in enumerate(self.objects):
            keep = (during_night==1) & (0<=airmass[i,:]) & (airmass[i,:]<=2.5)
            logger.info('Object %s: %s visible during night time (%.1f<airmass<%.1f)'%(obj.name,-np.any(keep) and 'not' or '',sum(keep) and airmass[i,keep].min() or np.nan,sum(keep) and airmass[i,keep].max() or np.nan))
        
    def plot(self,**kwargs):
        #-- plot
        #figure()
        hours = self.vis['hours']
        dates = self.vis['dates']
        alts = self.vis['alts']
        airmass = self.vis['airmass']
        during_night = self.vis['during_night']
        yaxis = kwargs.pop('yaxis','airmass')
        
        for i in range(len(airmass)):
            #-- only keep times during the night and with an airmass between 0 and 2.5
            keep = (during_night==1) & (0<=airmass[i,:]) & (airmass[i,:]<=2.5)

            if yaxis=='airmass':
                pl.plot_date(dates[keep],airmass[i,keep],xdate=True,tz=pytz.utc,**kwargs)
            elif yaxis=='degree':
                pl.plot_date(dates[keep],alts[i,keep]/pi*180,xdate=True,tz=pytz.utc,**kwargs)
        
        #-- take care of the x-axis labelling format
        pl.gca().xaxis.set_major_formatter(DateFormatter("%d %b '%y %H:%M"))
        #gca().fmt_xdata = DateFormatter('%H:%M')
        pl.gcf().autofmt_xdate()
        pl.xlabel('time (UTC)')
        
        #-- take care of labels and direction of the y-axis
        if yaxis=='airmass':
            pl.ylabel('Airmass')  
            if pl.ylim()[0]<pl.ylim()[1]:
                pl.ylim(pl.ylim()[::-1])
        elif yaxis=='degree':
            pl.ylabel('Altitude (degrees)')
        
        pl.gca().set_autoscale_on(False)
        start_nights = dates[:-1][np.diff(during_night)==1]
        end_nights = dates[:-1][np.diff(during_night)==-1]
        for start,end in zip(start_nights,end_nights):
            pl.fill([start,start,end,end],[pl.ylim()[0],pl.ylim()[1],pl.ylim()[1],pl.ylim()[0]],'k',alpha=0.5)
            
        #-- if only one night was asked, show the whole night
        if len(start_nights)==1:
            width = pl.xlim()[1] - pl.xlim()[0]
            pl.xlim(start-0.05*width,end+0.05*width)
        #-- add a grid
        pl.grid()
        
        
    #}
    
    
if __name__=="__main__":
    if len(sys.argv)<2:
        import doctest
        doctest.testmod()
        pl.show()
    else:
        from ivs.aux import loggers
        logger = loggers.get_basic_logger()
    
        eph = Ephemeris(objects=sys.argv[1].split(','))
        eph.visibility()
        eph.plot()
        pl.show()