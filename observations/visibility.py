"""
This module allows to calculate the altitude and airmass of an object in the sky, to determine whether it is 'observable' during day or night time,
and to get the separation to the moon.

Examples:

Plot the visibility of an object for 10 days, every 10 minutes.

>>> eph = Ephemeris()
>>> eph.set_objects(objects=['HD50230'])
>>> eph.set_date(startdate='2011/09/27 12:00:00.0',dt=10.,days=10)
>>> eph.set_site(sitename='lapalma')
>>> eph.visibility()
>>> eph.plot(fmt='bo')

]]include figure]]ivs_observations_visibility_period.png]

You can get the same behaviour in less commands:

>>> eph = Ephemeris(objects=['HD50230'],startdate='2011/09/27 12:00:00.0',dt=10.,days=10,sitename='lapalma')
>>> eph.visibility()

Or, equivalently,

>>> eph = Ephemeris()
>>> eph.visibility(objects=['HD50230'],startdate='2011/09/27 12:00:00.0',dt=10.,days=10)

Plot the visibility of 2 objects and the moon only for today:

>>> eph = Ephemeris(objects=['HD163506','HD143454'])
>>> eph.visibility(startdate=2456084.01922,days=1,sitename='palomar')
>>> eph.plot(moon=True,lw=2)

]]include figure]]ivs_observations_visibility_today.png]

One can also get the visibility of an object in a random capital at the middle of the night throughout the year:

>>> eph = Ephemeris(objects=['V* V441 Her'])
>>> eph.visibility(midnight=True,sitename='Paris')
>>> eph.plot(moon=False)

Calculate visibility for a list of objects, sites and times:

>>> eph = Ephemeris()
>>> days = np.array([2456084.0,2456085.02,2456085.29,2456085.55,2456095.6321,2456132.5478])
>>> objectnames = np.array(['HD163506','HD163506','HD163506','HD143454','HD143454','HD163506'])
>>> sitenames = np.array(['lapalma','lapalma','palomar','palomar','lapalma','palomar'])
>>> eph.visibility(multiple=True,startdate=days,objects=objectnames,sitename=sitenames)
>>> eph.plot(moon=True)

or set e.g. the site fixed:

>>> eph = Ephemeris()
>>> days = np.array([2456084.5,2456085.32,2456085.69,2456085.15,2456095.8321,2456132.1478])
>>> objectnames = np.array(['HD163506','HD163506','HD163506','HD143454','HD143454','HD163506'])
>>> eph.visibility(multiple=True,startdate=days,objects=objectnames,sitename='palomar')
>>> eph.plot(moon=True)


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
        
        You can set the site, date and objects via the keyword arguments in 'set_objects', 'set_site' and 'set_date'. Defaults are 'no objects', 'la palma' and
        'the current time'.
        You can also give the same information when calling 'visibility'.
        """
        self.set_site(**kwargs)
        self.set_date(**kwargs)
        self.set_objects(**kwargs)
    
    def __set_object(self,objectname):
        """
        Retrieve coordinates for an object in 'ephem'-style.
        """
        try:
	    jpos = sesame.search(objectname,db='S')['jpos']
	    logger.info('Found object %s at %s'%(objectname,jpos))
        except KeyError:
	    logger.warning('Object %s not found in SIMBAD, trying NED.'%(objectname))
	    try:
		jpos = sesame.search(objectname,db='N')['jpos']
	    except KeyError:
		logger.warning('Object %s not found in NED either.'%(objectname))
		raise IOError, 'No coordinates retrieved for object %s.'%(objectname)
        myobject = ephem.readdb("%s,f|M|A0,%s,8.0,2000"%(objectname,','.join(jpos.split())))
        return myobject
    
    #{ Set object, site and dates
    @decorators.filter_kwargs
    def set_objects(self,objects=None,**kwargs):
        """
        Initializes a list of objects by searching their coordinates on SIMBAD or NED.
        
        @keyword objects: a list of object names, resolved by SIMBAD or NED.
        @type objects: list
        """
        if objects is not None:
            objects = [self.__set_object(name) for name in objects]
            self.objects = objects
    
    def get_objectnames(self):
	"""
	Return the names of the objects.
	"""
	return [member.name for member in self.objects]
    
    @decorators.filter_kwargs
    def set_site(self,sitename='lapalma',sitelat=None,sitelong=None,siteelev=None):
        """
        Set the observing site.
        
        Supported are: 'lapalma', 'lasilla', 'palomar' and a bunch of capitals defined in pyephem.
        The other keywords allow to manually define a site.
        
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
        
        @keyword startdate: the exact time of the first visibility calculation 
        @type startdate: string or float
        @keyword days: total number of days on which to calculate visibilities (if visibility keyword 'multiple' = False)
        @type days: float
        @keyword dt: sampling time in minutes
        @type dt: float
        
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
        if self.dt is None and self.days is None:
	    logger.info('Set start date to %s'%(startdate))
	else:
	    logger.info('Set start date to %s, covering %d days per %d minutes'%(startdate,days,dt))
    
    #}

    #{ Get output objects and output object names
    def get_out_objects(self,indexes):
	"""
	Return a list of ephem object(s) corresponding to the rows in the output arrays defined by indexes.
	"""
	return self.objects[self._objectIndices[index]]
    
    def get_out_objectnames(self,indexes):
	"""
	Return an array of object names corresponding to the rows in the output arrays defined by indexes.
	"""
	return array([member.name for member in get_out_objects(indexes)])
    
    def get_out_sites(self,indexes):
	"""
	Return the ephem site(s) corresponding to the rows in the output arrays defined by indexes.
	"""
	return self._uniqueSites[self._siteIndices[indexes]]
    
    def get_out_sitenames(self,indexes):
	"""
	Return the names of the sites corresponding to the rows in the output arrays defined by indexes.
	"""
	return array([site.name for site in get_out_sites(indexes)])
    
    #}
    
    #{ Compute and plot visibilities
    
    def visibility(self,multiple=False,midnight=None,airmassmodel='Pickering2002',**kwargs):
        """
        Calculate ephemeri.
        
        *If 'multiple' = False, a single site definition and startdate (and/or values for the days and dt keywords) are expected. 
        If several astrophysical objects are given, ephemeri will be given for each object at the defined site and times. 
        If 'midnight' is not None, ephemeri are calculated in the middle of the night, for the whole year.
        
        *If 'multiple' = True, a list/array of site names and startdates can be given. Note: the length of keywords 'sitename', 'startdate' and 'objects' 
        can be either N or 1 and should be in the same format as previously. In this case, ephemeri are only calculated at the times given in 'startdate', 
        and thus the keywords 'days', 'dt' and 'midnight' are not used. Sites are currently only accessible throught the keyword 'sitename'.
           
        The result of this function is an attribute 'vis', a dictionary containing:
            - MJDs: Modified Julian Dates
            - dates calendar dates (Format 1)
            - alts: altitudes
            - moon_alts: altitudes of the moon
            - airmass
            - moon_airmass: airmasses of the moon
            - moon_separation: separation to moon
            - during_night: night time / day time (boolean 1-0)
            
        NOTE: use the 'get_out_objects', 'get_out_objectnames', 'get_out_sites' and 'get_out_sitenames' functions to retrieve the object, name of the object, 
        site and name of the site corresponding to a particular row in each of the arrays 'vis' contains.
        
        """
        #-- set optional additional keywords
        sun = ephem.Sun()
        moon = ephem.Moon()
        #moon_theta = 34.1 # minutes
        
        #-- keep track of which options were used for plotting purposes
        self._multiple = multiple
        self._midnight = midnight
        
        if not multiple:
	    #-- set the site, date and objects
	    self.set_site(**kwargs)
	    self.set_date(**kwargs)
	    self.set_objects(**kwargs)
	    
	    #-- set stuff for calculations
	    timestep_minutes = self.dt
	    total_days = self.days
	    if midnight is None:
		total_minutes = int(total_days*24*60./timestep_minutes)
	    else:
		total_minutes = 365
	    
	    #-- set initial arrays
	    alts = np.zeros((len(self.objects),total_minutes))
	    MJDs = np.zeros(total_minutes)
	    during_night = np.zeros(total_minutes)
	    moon_separation = np.zeros((len(self.objects),total_minutes))
	    moon_alts = np.zeros(total_minutes)
	    
	    #-- run over all timesteps
	    if midnight is None:
		for i in range(total_minutes):
		    self.thesite.date += ephem.minute*timestep_minutes
		    prev_set = float(self.thesite.previous_setting(sun))
		    prev_rise = float(self.thesite.previous_rising(sun))
		    MJDs[i] = float(self.thesite.date) + 15019.5              # the zeropoint of ephem is 1899/12/31 12:00:00.0
		    #-- compute the moon position
		    moon.compute(self.thesite)
		    moon_alts[i] = float(moon.alt)
		    if (prev_rise<=prev_set):
			during_night[i] = 1
		    for j,star in enumerate(self.objects):               
			star.compute(self.thesite)
			alts[j,i] = float(star.alt)
			moon_separation[j,i] = ephem.separation(moon,star)
			
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
		    MJDs[i] = float(self.thesite.date) + 15019.5
		    during_night[i] = 1
		    for j,star in enumerate(self.objects):               
			star.compute(self.thesite)
			alts[j,i] = float(star.alt)
		    i += 1
		    #-- set 1 day forwards
		    self.thesite.date += ephem.minute*60*24
	    
	    alts = alts.ravel()                                                                       
	    MJDs = np.outer(np.ones(len(self.objects)),MJDs).ravel()
	    during_night = np.outer(np.ones(len(self.objects),int),during_night).ravel()              
	    moon_separation = moon_separation.ravel()      #                                          
	    moon_alts = np.outer(np.ones(len(self.objects)),moon_alts).ravel()                        
	    self._objectIndices = np.outer(np.arange(len(self.objects)),np.ones(total_minutes)).ravel()
	    self._siteIndices = np.zeros_like(alts)
	    self._uniqueSites = [self.thesite]    

        else:
	    startdate = kwargs.get('startdate',None)
	    sitename = kwargs.get('sitename','lapalma')
	    objects = kwargs.get('objects',None)
	    if objects is None:
		objects = self.get_objectnames()
	    
	    # the other option is to iterate over given 'dates' and/or sites and objects
	    if not type(startdate) == np.ndarray: startdate = np.array([startdate])
	    if not type(sitename) == np.ndarray: sitename = np.array([sitename])
	    if not type(objects) == np.ndarray: objects = np.array([objects]).ravel()
	    
	    uniqueObjectNames = np.unique(objects)
	    uniqueSiteNames = np.unique(sitename)
	    self.set_objects(uniqueObjectNames)
	    self._objectIndices = np.zeros_like(objects,int)
	    self._siteIndices = np.zeros_like(sitename,int)
	    self._uniqueSites = []
	    for i in range(len(uniqueObjectNames)):
		self.objectIndices = np.where(objects == uniqueObjectNames[i],i,self._objectIndices)
	    for i in range(len(uniqueSiteNames)):
		self._siteIndices = np.where(sitename == uniqueSiteNames[i],i,self._siteIndices)
		self.set_site(uniqueSiteNames[i])
		self._uniqueSites.append(self.thesite)
	    
	    #-- define the iterator
	    it = np.broadcast(startdate,self._siteIndices,self._objectIndices)
	    
	    #-- set initial arrays
	    alts = np.zeros(it.shape[0])
	    MJDs = np.zeros_like(alts) 
	    during_night = np.zeros_like(alts)
	    moon_separation = np.zeros_like(alts)
	    moon_alts = np.zeros_like(alts)
	    
	    #-- run over all elements
	    for i,element in enumerate(it):
		#-- set the site and date
		self.thesite = self._uniqueSites[element[1]]
		self.set_date(startdate=element[0],days=None,dt=None)
		#-- determine whether the time corresponds to night or day
		prev_set = float(self.thesite.previous_setting(sun))
		prev_rise = float(self.thesite.previous_rising(sun))
		if (prev_rise<=prev_set):
		    during_night[i] = 1
		MJDs[i] = float(self.thesite.date)
		#-- compute the moon position
		moon.compute(self.thesite)
		moon_alts[i] = float(moon.alt)
		#-- compute the star's position
		star = self.objects[element[2]]
		star.compute(self.thesite)
		alts[i] = float(star.alt)
		moon_separation[i] = ephem.separation(moon,star)
	
	#-- calculate airmass
	airmass = obs_airmass.airmass(90-alts/pi*180,model=airmassmodel)
	moon_airmass = obs_airmass.airmass(90-moon_alts/pi*180,model=airmassmodel)
	
	#-- calculate dates for plotting and output
	self._plotdates = np.array([date2num(ephem.date(h).datetime()) for h in MJDs-15019.5])
	dates = np.array([ephem.date(h).datetime() for h in MJDs-15019.5])
	
	#-- the output dictionary
	self.vis = dict(MJDs=MJDs,dates=dates,alts=alts,airmass=airmass,during_night=during_night,moon_alts=moon_alts,moon_airmass=moon_airmass,moon_separation=moon_separation)
	for i,obj in enumerate(self.objects):
	    keep = (self._objectIndices == i) & (during_night==1) & (0<=airmass) & (airmass<=2.5)
	    logger.info('Object %s: %s visible during night time (%.1f<airmass<%.1f)'%(obj.name,-np.any(keep) and 'not' or '',sum(keep) and airmass[keep].min() or np.nan,sum(keep) and airmass[keep].max() or np.nan))


    def plot(self,**kwargs):
        #-- input
        moon = kwargs.pop('moon',True)
        alts = self.vis['alts']
        airmass = self.vis['airmass']
        moon_airmass = self.vis['moon_airmass']
        during_night = self.vis['during_night']
        yaxis = kwargs.pop('yaxis','airmass')          # or 'alts'
        
        factor = (yaxis=='alts') and 180./pi or 1.
        
        #-- plot
        fmt = kwargs.pop('fmt',False)
        plotcolors = ['b','g','r','c','m','k']
        plotsymbols = ['o','.','v','^','<','>','s','p','*','h','H','+','d','D']
        
        pl.figure(figsize=(16,8))
        
        #-- run over all objects and plot them
        for i in range(len(self.objects)):
            #-- only keep times during the night and with an airmass between 0 and 2.5
            keep = (self._objectIndices == i) & (during_night==1) & (0<=airmass) & (airmass<=3.5)
            ifmt = fmt and fmt or plotsymbols[i/6]+plotcolors[np.fmod(i,6)]
            pl.plot_date(self._plotdates[keep],self.vis[yaxis][keep]*factor,xdate=True,tz=pytz.utc,fmt=ifmt,label=self.objects[i].name,**kwargs)
            
            if moon:
		keep = (self._objectIndices == i) & (during_night==1) & (0<=moon_airmass) & (moon_airmass<=3.5)
		
		#-- check the moon separation!
		if keep.sum():
		    index = np.argmin(self.vis['moon_separation'][keep]*factor)
		    closest = (self.vis['moon_separation'][keep]*factor)[index]
		    if closest<40.:
			t = self._plotdates[keep][index]
			time_closest = pytz.datetime.datetime.fromordinal(int(t)) + pytz.datetime.timedelta(days=t-int(t))
			logger.warning('Object-moon distance is %.1f deg at %s'%(closest,time_closest))
            
        if moon:
	    keep = (during_night==1) & (0<=moon_airmass) & (moon_airmass<=3.5)
	    pl.plot_date(self._plotdates[keep],self.vis['moon_'+yaxis][keep]*factor,xdate=True,tz=pytz.utc,fmt='yo',label='Moon')

        #-- take care of labels and direction of the y-axis, and the legend
        if yaxis=='airmass':
            pl.ylabel('Airmass')  
            if pl.ylim()[0]<pl.ylim()[1]:
                pl.ylim(pl.ylim()[::-1])
        elif yaxis=='alts':
            pl.ylabel('Altitude (degrees)')
        pl.legend(loc='upper right',numpoints=1,borderpad=0.2,labelspacing=0.2,handlelength=0.2)
        xwidth = pl.xlim()[1] - pl.xlim()[0]
        pl.xlim(pl.xlim()[0],pl.xlim()[1]+xwidth*0.15)

        #-- take care of the x-axis labelling format
        pl.gca().xaxis.set_major_formatter(DateFormatter("%d %b '%y %H:%M"))
        #gca().fmt_xdata = DateFormatter('%H:%M')
        pl.gcf().autofmt_xdate()
        pl.xlabel('time (UTC)')
        
        #-- the following only works when visibility is called with multiple = False and midnight = None
        if not self._multiple and self._midnight == None:
	    pl.gca().set_autoscale_on(False)
	    start_nights = self._plotdates[:-1][np.diff(during_night)==1]
	    end_nights = self._plotdates[:-1][np.diff(during_night)==-1]
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