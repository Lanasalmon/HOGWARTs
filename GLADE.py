from astroquery.simbad import Simbad
from astroquery.vizier import VizierClass
from astropy.coordinates import SkyCoord
import pandas as pd
import astropy.units as u

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