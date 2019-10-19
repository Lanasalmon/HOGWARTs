from astroquery.simbad import Simbad
from astroquery.vizier import VizierClass
from astropy.coordinates import SkyCoord
import pandas as pd
import astropy.units as u
import numpy as np

def GLADEV2coordinates(distmax,distmin):

    vizier = VizierClass(
    row_limit=-1, columns=['HyperLEDA', '_RAJ2000', '_DEJ2000', 'Dist','Bmag'],column_filters={'Bmag':'!=null'})
    cat, = vizier.get_catalogs('VII/281/glade2')
    data=pd.DataFrame({'RA': cat['_RAJ2000'], 'Dec':cat['_DEJ2000'],'dist':cat['Dist'],'Bmag':cat['Bmag'],'HyperLEDA':cat['HyperLEDA']})
    msk1=data[['dist']]<=distmax
    msk2=data[['dist']]>=distmin
    msk3=data[['dist']]>0
    msk4=data[['dist']]!='NaN'
    msk5=data[['Bmag']]!='NaN' 
    msk6=data[['Bmag']]!='null' 
    msk7=data[['Bmag']]>0 
    msk=pd.concat((msk1,msk2,msk3,msk4,msk5,msk6,msk7),axis=1)
    slct=msk.all(axis=1)
    data=data.ix[slct]

    coordinates=SkyCoord(data['RA'].values*u.deg, data['Dec'].values*u.deg,data['dist'].values*u.Mpc)
    return coordinates,data



def Viziergalaxies(split_ra2, split_dec2,i,distest,diststd):
    diff_values_ra = np.diff(split_ra2[i])

    #if the contour crosses 0/360 degrees, splut it up again so we can perform query in Vizier
    if max(diff_values_ra)>350:
        

        split_index_ra = np.where(np.abs(diff_values_ra) > 350)[0]
        split_index_ra += 1

        #split at point where 0 goes to 360
        split_ra_360, split_dec_360=split_ra_dec(split_ra2, split_dec2, split_index_ra, i)
        
        #get centre coordinates and diameter of circle which encloses contour region
        centreras, centredecs, maxdists =get_centres(split_ra_360, split_dec_360,i)
        
        #query vizier using the circle and save in database for that contour
        galaxy_ra, galaxy_dec, galaxy_dist, galaxy_Bmag, galaxyname=make_database_list(centreras, centredecs, maxdists, distest, diststd)
        
        
    else: 
        #if the contour doesn't cross 0/360 degrees, get centre coordinates and diameter of circle
        centrera,centredec, maxdist=radius(split_ra2[i], split_dec2[i],i)
        
        #query vizier using the circle and save in database for that contour
        galaxy_ra, galaxy_dec, galaxy_dist, galaxy_Bmag,galaxyname=make_database( centrera, centredec, maxdist, distest, diststd)
    return galaxy_ra, galaxy_dec, galaxy_dist, galaxy_Bmag,galaxyname
