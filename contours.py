from __future__ import print_function
import numpy as np
import healpy as hp
from scipy.interpolate import LinearNDInterpolator
from skimage import measure
import matplotlib.pyplot as plt
import operator
import pandas as pd 
import math
import pickle
import numpy.ma as ma
from astropy.coordinates import SkyCoord
import astropy.units as u

def radius(ra, dec,i):
    """
    Calculate centre coordinates and diameter of circle to enclose contour.
    
    Parameters:
    -----------
    ra, dec :  list
        a list of right ascension and declination coordinates for a contour.
    i : int
        index corresponding to contour being analysed.
    Return:
    -------
    centrera : float 
        Centre RA coordinate for contour
    centredec : float
        Centre Dec coordinate for contour
    maxdist : float
        Diameter of contour
    
    
    """

    width=max(ra)-min(ra)
    height=max(dec)-min(dec)
    centrera=min(ra)+width/2
    centredec=min(dec)+height/2
    deltadec=max(dec)-min(dec)
    deltadec=deltadec/2
    minra_index, minra_value = min(enumerate(ra), key=operator.itemgetter(1))
    maxra_index, maxra_value = max(enumerate(ra), key=operator.itemgetter(1))
    mindec_index, mindec_value = min(enumerate(dec), key=operator.itemgetter(1))
    maxdec_index, maxdec_value = max(enumerate(dec), key=operator.itemgetter(1))

    dist1=np.sqrt((ra[minra_index]-centrera)**2+(dec[minra_index]-centredec)**2)
    dist2=np.sqrt((ra[mindec_index]-centrera)**2+(dec[mindec_index]-centredec)**2)
    dist3=np.sqrt((ra[maxra_index]-centrera)**2+(dec[maxra_index]-centredec)**2)
    dist4=np.sqrt((ra[maxdec_index]-centrera)**2+(dec[maxdec_index]-centredec)**2)
    distarray=np.array([dist1,dist2,dist3,dist4])
    maxdist=max(distarray)
    return centrera, centredec, maxdist


def hpix_contours(m,levels=[0.9],nest=True):
	"""
	Compute iso-lines for a healpix map
	
	Parameters:
	-----------
	m: 1D numpy array
		The input healpix map.
	levels: list of floats
		The values for which to compute the iso-lines. Default: [0.5,0.9]
	nest: boolean
		If True, nested ordering is assumed for the healpix map.
	
	Return:
	-------
	contours: a list of masked numpy arrays
		Each element in the list is a 2D numpy array containing the contour lines corresponding to a given level.
		Each contour c in the list has shape (2,N): c[0] represents the RA and c[1] the Dec coordinates of the 
		points constituting the contour.
	"""
	
	
	nside = hp.npix2nside(len(m))
	
	# define a grid over which to evaluate the density
	lat = np.linspace(-np.pi/2.,np.pi/2.,300)
	lon = np.linspace(0.,2.*np.pi,300)
	
	# evaluate the map and keep track of both coordinates and indices where
	# the evaluation is done
	values = np.zeros([300*300,2])
	points = np.zeros([300*300,2])
	MAP = np.zeros([300,300])

	
	k=0
	for i in range(300):
		for j in range(300):
			p = hp.ang2pix(nside,np.pi/2.-lat[j],lon[i],nest=nest)
			MAP[i,j] = m[p]
			values[k] = np.array([lon[i],lat[j]])
			points[k,0] = i
			points[k,1] = j
			k = k+1

	#construct a linear interpolator to get interpolated
	#coordinates from fractional indices
	lint = LinearNDInterpolator(points,values)

	#use skimage find_contours method to find fractional indices
	#defining the contours, then compute corresponding coordinates
	#by using the linear interpolator
	contours = []
	for l in levels:
		contour_components = measure.find_contours(MAP,l)
		#the find_contours function above returns a list of the
		#connected components of the contour. We unpack it,
		#compute the interpolated coordinates for each, and
		#store them in one array, separating them with nans
		
		whole_contour = np.array([[np.nan,np.nan]])

		for contour in contour_components:
			cont_coords = lint(contour)
			whole_contour = np.concatenate((whole_contour,cont_coords),axis=0)
			whole_contour = np.concatenate((whole_contour,[[np.nan,np.nan]]),axis=0)
		
		#then we mask the nans
		C = whole_contour.transpose()
		contours.append(np.ma.masked_invalid(C))
	
	return contours

def integrated_probability(skymap):
	"""
	Take a (healpix) skymap and return the corresponding integrated probability skymap.
	The result can be then used to compute the confidence regions.
	"""

	sort_idx = np.argsort(skymap)[::-1]

	csm = np.empty(len(skymap))

	csm[sort_idx] = np.cumsum(skymap[sort_idx])

	return csm

def split_contours(contours, percent,index):
    """
    Split masked array into list of individual arrays.
    
    Parameters:
    -----------
    contours: a list of masked numpy arrays
        Each element in the list is a 2D numpy array containing the contour lines corresponding to a given level.
        Each contour c in the list has shape (2,N): c[0] represents the RA and c[1] the Dec coordinates of the 
        points constituting the contour.
    index: int
        Index of contours array to split. 0 corresponds to 90%, 1 to 50%, 2 to 99$
    
    Return:
    -------
    split_ra, split_dec : list
        a list of arrays of right ascension and declination coordinates, each array in the list corresponding to a separate contour.
        
    """
    ra,dec=contours[index][0]*180/np.pi, contours[index][1]*180/np.pi
    split_ra=[]
    split_dec=[]

    m=ma.getmaskarray(ra)
    md=ma.getmaskarray(dec)
    d = np.diff(m)
    dd=np.diff(md)
    cuts = np.flatnonzero(d) +1
    cutsd=np.flatnonzero(dd)+1
    

    asplit = np.split(ra, cuts)
    dsplit=np.split(dec,cutsd)
    msplit = np.split(m, cuts)

    for i in range(0,len(asplit)):
        if len(asplit[i])>1:
            split_ra.append(asplit[i])
    for i in range(0,len(dsplit)):
        if len(dsplit[i])>1:
            split_dec.append(dsplit[i])
    return split_dec,split_ra

def join_0_360(split_ra, split_dec):
    """
    Join together contours which cross the 0/360 RA boundary.
    
    Parameters:
    -----------
    split_ra, split_dec :  list
        a list of arrays of right ascension and declination coordinates, each array in the list corresponding to a separate contour.
        
    Return:
    -------
    split_ra2, split_dec2 : list
        a list of arrays of right ascension and declination coordinates, each array in the list corresponding to a separate whole contour.
    """
    deletes=[]
    adds=[]
    newras=[]
    newdecs=[]
    for i in range(0,len(split_ra)):
            for j in range(0,len(split_ra)):
                if i!=j:
                    if np.abs(split_dec[i][0]-split_dec[j][0])<0.5 or np.abs(split_dec[i][-1]-split_dec[j][0])<0.5:
                        if j not in deletes and j not in adds:
                            deletes.append(j)
                            adds.append(i)

                            newra=np.concatenate((split_ra[i], split_ra[j]))
                            newdec=np.concatenate((split_dec[i], split_dec[j]))
                            newdecs.append(newdec)
                            newras.append(newra)

    split_ra2=[]
    split_dec2=[]         

    for i in range(0,len(adds)):
        split_ra[adds[i]]=newras[i]
        split_dec[adds[i]]=newdecs[i]

    for i in range(0,len(split_ra)):
        if i not in deletes:
            split_ra2.append(split_ra[i])
            split_dec2.append(split_dec[i])
    return split_ra2, split_dec2        

def contour_plots(split_ra2, split_dec2,graceid, prelim, percentage):
    """
    Plot contours.
    
    Parameters:
    -----------
    split_ra2, split_dec2 :  list
        a list of arrays of right ascension and declination coordinates, each array in the list corresponding to a separate whole contour.
    graceid :  string
        identifier for GW event eg. S190707q
    prelim :  string
        notice type eg. Preliminary
    percentage : int
        percentage confidence region eg. 90%.

    Return:
    -------
    Plot.
    """
    plt.clf()
    fig5 = plt.figure(figsize=(15, 7), dpi=100)
    ax = plt.axes(
    [0.05, 0.05, 0.9, 0.9],
    projection='astro degrees aitoff')
    ax.grid()

    for i in range(0,len(split_ra2)):
        f=i+1
        skycoordcontours=SkyCoord(ra=split_ra2[i]*u.degree, dec=split_dec2[i]*u.degree, frame='icrs')
        centrera,centredec, maxdist=radius(split_ra2[i], split_dec2[i],i)
        centrecoordcontours=SkyCoord(ra=centrera*u.degree, dec=centredec*u.degree, frame='icrs')
        ax.plot(skycoordcontours.ra.deg,skycoordcontours.dec.deg,'-',transform=ax.get_transform('world'))
        ax.text(centrecoordcontours.ra.deg,  centrecoordcontours.dec.deg, '%i' %f, transform=ax.get_transform('world'))


    plt.savefig(graceid+prelim+'contourplot' +str(percentage)+'.png')
    plt.show()
    return

def split_ra_dec(split_ra2, split_dec2, split_index_ra, i):
    """
    Split contours where they cross the 0/360 degree RA boundary.
    
    Parameters:
    -----------
    split_ra2, split_dec2 :  list
        a list of arrays of right ascension and declination coordinates, each array in the list corresponding to a separate whole contour.
    split_index_ra : list
        index at which the contour coordinates pass the 0/360 degree RA boundary.
    i : int
        contour index

    Return:
    -------
    split_ra_360, split_dec_360 : a list of arrays of right ascension and declination coordinates, each array in the list corresponding to a separate contour where contours that cross the RA 360/0 degree boundary are split into separate contour arrays in the list. 
    
    """
    split_ra_360 = np.split(split_ra2[i], split_index_ra)
    split_dec_360 = np.split(split_dec2[i], split_index_ra)
    return split_ra_360, split_dec_360

def get_centres(split_ra_360, split_dec_360,i):
    """
    Create list of centre coordinates and diameter of a circle that would enclose each contour.
    
    Parameters:
    -----------
    split_ra_360, split_dec_360  :  list
        a list of arrays of right ascension and declination coordinates, each array in the list corresponding to a separate contour where contours that cross the RA 360/0 degree boundary are split into separate contour arrays in the list. 
    i : int
        index corresponding to percentage confidence region to analyse.

    Return:
    -------
    centreras : list
        list of centre RA coordinates for each contour
    centredecs : list
        list of centre Dec coordinates for each contour
    maxdists : list
        lost of diameter of circles that can enclose each contour.
    """
    centreras=[]
    centredecs=[]
    maxdists=[]
    for j in range(0,len(split_ra_360)):
        centrera,centredec, maxdist=radius(split_ra_360[j], split_dec_360[j],i)
        centreras.append(centrera)
        centredecs.append(centredec)
        maxdists.append(maxdist)
    return centreras, centredecs, maxdists 
    
import healpy as hp
import numpy as np
import sys
sys.path.append('/anaconda3/lib/python3.6/site-packages/MOCPy-0.4.9-py3.6.egg')
from mocpy import MOC

from math import log

class MOC_confidence_region(object):
    """
    Multi-Order coverage map (MOC) of sky areas enclosed within a contour plot at a given confidence level.
    """

    def read_skymap(self, infile):
        """Reading healpix skymap.
        
        Input parameters
        ----------------
        infile : string
              LVC probability sky localization in healpix format
              
        Return
        -------
        hpx : list
            1D array of values (probability stored in each pixel)
        nside : int
           skymap resolution
        """      
        
        self.hpx = hp.read_map(infile, verbose = False)
        npix = len(self.hpx)
        self.nside = hp.npix2nside(npix)
        
        return self.hpx, self.nside
 
    def ipixs_in_percentage(self, percentage):
        """Finding ipix indices confined in a given percentage.
        
        Input parameters
        ----------------
        percentage : float
                 fractional percentage from 0 to 1  
        
        Return
        ------- 
        ipixs : list
              indices of pixels
        """
        
        sort = sorted(self.hpx, reverse = True)
        cumsum = np.cumsum(sort)
        index, value = min(enumerate(cumsum), key = lambda x: abs( x[1] - percentage ))

        index_hpx = range(0, len( self.hpx ))
        hpx_index = np.c_[self.hpx, index_hpx]

        sort_2array = sorted(hpx_index, key = lambda x: x[0], reverse = True)
        value_contour = sort_2array[0:index]

        j = 1 
        table_ipix_contour = [ ]

        for i in range (0, len(value_contour)):
            ipix_contour = int(value_contour[i][j])
            table_ipix_contour.append(ipix_contour)
    
        self.ipixs = table_ipix_contour
          
        return self.ipixs
     
    def sky_coords(self):
        """Creating an astropy.table with RA[deg] and DEC[deg] ipix positions
        
        Return
        ------- 
        contour_ipix : 
                    sky coords in degrees
        """
       
        # from index to polar coordinates
        theta, phi = hp.pix2ang(self.nside, self.ipixs)

        # converting these to right ascension and declination in degrees
        ra = np.rad2deg(phi)
        dec = np.rad2deg(0.5 * np.pi - theta)
        
        # creating an astropy.table with RA[deg] and DEC[deg]
        from astropy.table import Table

        self.contour_ipix = Table([ra, dec], names = ('RA[deg]', 'DEC[deg]'), 
                             meta = {'ipix': 'ipix table'})

        return self.contour_ipix
    
    def moc_order(self):
        """Setting MOC order.
        
        Return
        ------- 
        moc_order : int
              
        """       
        
        self.moc_order = int(log( self.nside, 2))
        print(self.moc_order)
     
        return self.moc_order

    def create_moc(self):
        """Creating a MOC map from the contour_ipix table."""
        import sys 
        sys.path.append('/anaconda3/lib/python3.6/site-packages/MOCPy-0.4.9-py3.6.egg')
        from mocpy import MOC
        import astropy.units as u
        self.moc = MOC.from_lonlat(self.contour_ipix['RA[deg]'].T*u.deg,self.contour_ipix['DEC[deg]'].T*u.deg,max_norder=self.moc_order)

        return self.moc

    def write_moc(self, percentage, short_name):
        """Writing MOC file in fits format.
        
        Input parameters
        ----------------
        percentage : float
                 fractional percentage from 0 to 1 converted into a string
        short_name : str
                 file output
        """
        
        return self.moc.write(short_name + '_MOC_' + str(percentage), format = 'fits')

    
    
    def contour_plot(self,infile, percentage, short_name=''):
        """Creating/Writing a MOC contour region at a fixed level of probability.
        
        Input parameters
        ---------------
        infile : string
              LVC probability sky localization in healpix format
        percentage : float
                 fractional percentage from 0 to 1
        """
        ipixes=[]
        self.read_skymap(infile)
        print('i read skymap')
        self.ipixs_in_percentage(percentage)
        print('i got percentages')
        ipixes.append(self.ipixs_in_percentage(percentage))
        self.sky_coords()
        self.moc_order()
        self.create_moc()

        return self.write_moc(percentage, short_name), ipixes
    
    
class LocalizeSource(MOC_confidence_region):
    """The class is designed to localize an astrophysical source inside a probability skymap."""
        
    def __init__(self):
        self.aladin = AladinScriptCommands()
        
    def in_skymap(self, infile, ra, dec, percentage, label=' ', show = True):
        """Checking if an object falls in a given probability level.
        
        Input parameters
        ---------------
        infile : string
                LVC probability sky localization in healpix format
        percentage : float
                fractional percentage from 0 to 1 
        ra, dec : float
                sky coordinates in degrees
        label : string
                name of the object (optional, by default = '')
        show = True
                show the MOC confidence region and the object in the Aladin planes;
                otherwise no (optional; by default = True )
        """
        
        self.read_skymap(infile)
        ipixs=self.ipixs_in_percentage(percentage)
        
        theta = 0.5 * np.pi - np.deg2rad(dec)
        phi = np.deg2rad(ra)
        ipix = hp.ang2pix(self.nside, theta, phi)
        
        is_there = ipix in ipixs

        if is_there == True:
            print ("The sky coord", "ra="+str(ra)+"°,","dec="+str(dec)+"°", "(label:" + label+")", \
                    "lies within the", str(percentage*100)+'%', "c.l.")
        else:
            print ("The sky coord", "ra="+str(ra)+"°,","dec="+str(dec)+"°", "(label:" + label+")", \
                    "is outside the", str(percentage*100)+'%', "c.l.")

        if show == True:
            self.aladin.draw_newtool("sources")
            self.aladin.draw_string(ra, dec, label)
            self.aladin.draw_circle(ra, dec, size = '10arcmin')
    
    def pinpoint(self, infile, ra, dec, from_percentage=0.1, to_percentage=1,
                 resolution_percentage=0.1, label=' ', show=True):
            
        """Find in which confidence level the source falls.
        
        Input parameters
        ---------------
        infile : string
            LVC probability sky localization in healpix format
        from_percentage : float
            fractional percentage from 0 to 1
        to_percentage : float
            fractional percentage from 0 to 1
        resolution_percentage : float
            fractional percentage from 0 to 1
        ra, dec : float
            sky coordinates in degrees
        label : string
            name of the object (optional, by default = '')
        show = True
            show the MOC confidence region and the object in the Aladin planes;
        otherwise no (optional; by default = True )
        """
            
        self.read_skymap(infile)

        # finding ipix
        theta = 0.5 * np.pi - np.deg2rad(dec)
        phi = np.deg2rad(ra)
        ipix = hp.ang2pix(self.nside, theta, phi)

        find="n"
        while from_percentage < to_percentage or find=="y":
            ipixs = self.ipixs_in_percentage(from_percentage)
            is_there = ipix in ipixs

            if is_there != True:
                print ("It is not localized within the", str(from_percentage*100)+'%', "c.l.")
                from_percentage = from_percentage + resolution_percentage
            else:
                find="y"      
                print ("The sky coord", "ra="+str(ra)+"°,","dec="+str(dec)+"°", "(label:" + label+")" \
                       "lies within the", str(from_percentage*100)+'%', "c.l.")
                    
                self.aladin.draw_newtool("sources")
                self.aladin.draw_string(ra, dec, label)
                self.aladin.draw_circle(ra, dec, size = '10arcmin')
                    
                return find

def in_hull(p, hull):
    """
        Test if points in `p` are in `hull`
        
        `p` should be a `NxK` coordinates of `N` points in `K` dimensions
        `hull` is either a scipy.spatial.Delaunay object or the `MxK` array of the
        coordinates of `M` points in `K`dimensions for which Delaunay triangulation
        will be computed
        """
    from scipy.spatial import Delaunay
    if not isinstance(hull,Delaunay):
        hull = Delaunay(hull)
    
    return hull.find_simplex(p)>=0

def pixels_in_region(p,pipeline_event):
    """
    Determine the pixels that lie within the percentage confidence region
    
    Parameters:
    -----------
    p : float
        percentage confidence region (eg. 0.9 corresponds to 90%)
    pipeline_event : param
        skymap

    Return:
    -------
    moc : multi-order map
        MOC map
    ipixes : array
        pixels within percentage confidence region

    """

    moc = MOC_confidence_region()
    _, ipixes=moc.contour_plot(infile =pipeline_event, 
             percentage = p, 
             short_name = pipeline_event)
    return moc, ipixes

def checkifinellipsemoc(moc,nside, ipixes, ra,dec,dist,Bmag,name,f):
    """
    Join together contours which cross the 0/360 boundary.
    
    Parameters:
    -----------
    moc :  
        a list of arrays, each one corresponding to a separate contour.
    nside :  int
        resolution of healpix map
    ipixes :  
        pixels within percentage confidence region
    ra, dec : list
        list of ra, dec coordinates of galaxies
    dist : list 
        distances to galaxies
    Bmag : list
        B magnitudes of galaxies
    name : list
        list of galaxy names
    f : int
        index corresponding to contour

    Return:
    -------
    keepgalaxyra, keepgalaxydec : list
        list of ra, dec coordinates within contour region
    keepgalaxydist : list
        list of distances corresponding to galaxies within contour region
    keepgalaxyBmag : list
        list of B magnitude corresponding to galaxies within contour region
    keepgalaxyname : list
        list of names corresponding to galaxies within contour region
    contours : list
        list of contours corresponding to galaxies within contour region
    
    """

    galipix=[]
    theta= 0.5 * np.pi - np.deg2rad(dec)
    phi = np.deg2rad(ra)
    galipix = hp.ang2pix(nside, theta, phi)
        
    keepgalaxyra=[]
    keepgalaxydec=[]
    keepgalaxydist=[]
    keepgalaxyBmag=[]
    keepgalaxyname=[]
    contours=[]
    import random
    import math

    s = set(ipixes[0])
    for i,x in enumerate(galipix):
        if x in s:
            keepgalaxyra.append(ra[i])
            keepgalaxydec.append(dec[i])
            keepgalaxydist.append(dist[i])
            keepgalaxyBmag.append(Bmag[i])
            keepgalaxyname.append(name[i])
            contours.append(f)

    return keepgalaxyra,keepgalaxydec, keepgalaxydist, keepgalaxyBmag, keepgalaxyname,contours
