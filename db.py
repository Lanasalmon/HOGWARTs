import sqlite3
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
import astropy.units as u
import numpy as np
import numpy as np
import operator
import psycopg2
import time
from psycopg2.extensions import AsIs
import pandas as pd
import numpy as np
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

def make_database_list(database, table, centrera, centredec, maxdist, distest, diststd):
    """
    Create and fill database with galaxies within Vizier queried circle 
    for contours which are split at 0/360 degree RA boundary.
    
    Parameters:
    -----------
    database  :  string
        name of database to save to
    table : string
        name of table to create
    centrera : float 
        Centre RA coordinate for contour
    centredec : float
        Centre Dec coordinate for contour
    maxdist : float
        Diameter of contour
    distest : float
        LIGO/Virgo estimated distance
    diststd : float
        Error on LIGO/Virgo estimated distance

    Return:
    -------
    ralist, declist : list
        galaxy coordinates within contour
    distlist : list
        galaxy distances within contour
    maglist : list
        galaxy B magnitudes within contour
    name : list
        galaxy names within contour
    """
    conn = psycopg2.connect("dbname='%s' user='lana' host='localhost' password='myPass'" %database)
    conn.autocommit = True
    cursor = conn.cursor()
    cursor.execute('SET SESSION synchronous_commit TO OFF;')
    distmax=distest+diststd
    distmin=distest-diststd
    cursor.execute('''CREATE TABLE %s (ra FLOAT,dec FLOAT, distance FLOAT, mag FLOAT,  name VARCHAR(100))''' % table)
    apassconfig = {'catname':'VII/281/glade2',
    'nicecatname':'GLADE',
    'ra':'RAJ2000',
    'dec':'DEJ2000',
    'distance':'Dist',
    "B":'Bmag',
    "maj":'maj',
    "nameGWGC":'GWGC',
    "nameHyperLEDA": 'HyperLEDA',
    }
    ralists=[]
    declists=[]
    distlists=[]
    maglists=[]
    namelists=[]
    for i in range(0,len(centrera)):
        coord = SkyCoord(ra=centrera[i]*u.degree, dec=centredec[i]*u.degree, frame='icrs')
        configdict = apassconfig
        
        v = Vizier(columns=['RAJ2000', 'DEJ2000','Bmag', 'maj', 'Dist','GWGC', 'HyperLEDA'],column_filters={"Bmag":"!= null"})
        v.ROW_LIMIT = -1
       
        result = v.query_region(coord, radius=maxdist[i]*u.deg, catalog=configdict['catname'])
        mask = (result[0]['Dist'] >= distmin) & (result[0]['Dist'] <=distmax)
        ralist=np.array(result[0][mask]['RAJ2000'])
        declist=np.array(result[0][mask]['DEJ2000'])
        distlist=np.array(result[0][mask]['Dist'])
        maglist=np.array(result[0][mask]['Bmag'])
        name=np.array(result[0][mask]['HyperLEDA'])
        ralists.append(ralist)
        declists.append(declist)
        distlists.append(distlist)
        maglists.append(maglist)
        namelists.append(name)

    ralist=np.concatenate((ralists[0], ralists[1])) 
    declist=np.concatenate((declists[0], declists[1]))   
    maglist=np.concatenate((maglists[0], maglists[1]))  
    distlist=np.concatenate((distlists[0], distlists[1]))  
    name=np.concatenate((namelists[0], namelists[1])) 
    
    df = pd.DataFrame({"ra" : ralist, "dec" : declist, "distance" : distlist, "mag" : maglist, "name": name})
    df=df[(df['distance'] >= distmin) & (df['distance'] <=distmax )]
    df = df[['ra', 'dec','distance','mag', 'name']]
    df.to_csv("test1.csv", header=False, index=False)
    f = open('test1.csv')
    cursor.copy_from(f, str(table), sep=",")

    return ralist, declist, distlist, maglist, name
def make_database(database, table, centrera, centredec, maxdist, distest, diststd):
    """
    Create and fill database with galaxies within Vizier queried circle.
    
    Parameters:
    -----------
    database  :  string
        name of database to save to
    table : string
        name of table to create
    centrera : float 
        Centre RA coordinate for contour
    centredec : float
        Centre Dec coordinate for contour
    maxdist : float
        Diameter of contour
    distest : float
        LIGO/Virgo estimated distance
    diststd : float
        Error on LIGO/Virgo estimated distance

    Return:
    -------
    ralist, declist : list
        galaxy coordinates within contour
    distlist : list
        galaxy distances within contour
    maglist : list
        galaxy B magnitudes within contour
    name : list
        galaxy names within contour
    """
    conn = psycopg2.connect("dbname='%s' user='lana' host='localhost' password='myPass'" %database)
    conn.autocommit = True
    cursor = conn.cursor()
    cursor.execute('SET SESSION synchronous_commit TO OFF;')
    distmax=distest+diststd
    distmin=distest-diststd
    cursor.execute('''CREATE TABLE %s (ra FLOAT,dec FLOAT, distance FLOAT, mag FLOAT, name VARCHAR(100))''' % table)
    apassconfig = {'catname':'VII/281/glade2',
    'nicecatname':'GLADE',
    'ra':'RAJ2000',
    'dec':'DEJ2000',
    'distance':'Dist',
    "B":'Bmag',
    "maj":'maj',
    "nameGWGC":'GWGC',
    "nameHyperLEDA": 'HyperLEDA',
    }
    coord = SkyCoord(ra=centrera*u.degree, dec=centredec*u.degree, frame='icrs')
    configdict = apassconfig
    
    v = Vizier(columns=['RAJ2000', 'DEJ2000','Bmag', 'maj', 'Dist','GWGC', 'HyperLEDA'],column_filters={"Bmag":"!= null"})
    v.ROW_LIMIT = -1
    print(coord)
    ralist=[]
    declist=[]
    distlist=[]
    maglist=[]
    majlist=[]
    cc=[]
    name=[]
    result = v.query_region(coord, radius=maxdist*u.deg, catalog=configdict['catname'])
    print('done query')
    if len(result)<1:
        ralist=[]
        declist=[]
        distlist=[]
        maglist=[]
        majlist=[]
        cc=[]
        name=[]
    else:

        mask = (result[0]['Dist'] >= distmin) & (result[0]['Dist'] <=distmax)
        ralist=np.array(result[0][mask]['RAJ2000'])
        declist=np.array(result[0][mask]['DEJ2000'])
        distlist=np.array(result[0][mask]['Dist'])
        maglist=np.array(result[0][mask]['Bmag'])
        name=np.array(result[0][mask]['HyperLEDA'])
    print('appended')
        

        
    df = pd.DataFrame({"ra" : ralist, "dec" : declist, "distance" : distlist, "mag" : maglist, "name": name})
    df=df[(df['distance'] >= distmin) & (df['distance'] <=distmax )]
    df = df[['ra', 'dec','distance','mag','name']]
    print('df made')
    df.to_csv("test1.csv", header=False, index=False)
    f = open('test1.csv')
    print('csv made')
    cursor.copy_from(f, str(table), sep=",")
    print('table made')
         

    return ralist, declist, distlist, maglist,name
def database_to_csv(database, table, ra,dec,dist, bmag, prob,distmu,distsigma,distnorm, galname,finalprob,contour):
    """
    Create and fill csv with galaxies from database.
    
    Parameters:
    -----------
    database  :  string
        name of database to save to
    table : string
        name of table to create
    ra, dec : list
        list of galaxy coordinates within contours
    dist : list
        list of galaxy distances within contours
    bmag : list
        list of galaxy B magnitudes within contours
    prob : list
        list of skymap probabilities at galaxy coordinates within contours
    distmu: list 
        list of skymap mean distances at galaxy coordinates within contours
    distssigma: list 
        list of skymap distance standard deviations at galaxy coordinates within contours
    distnorm: list 
        list of skymap normalisation factors at galaxy coordinates within contours
    galname: list 
        list of galaxy names within contours
    finalprob: list 
        list of computed galaxy probabilities within contours
    contour: list 
        list of contours corresponding to galaxies
    
    
    
    Return:
    -------
    
    """
    conn = psycopg2.connect("dbname='%s' user='lana' host='localhost' password='myPass'" %database)
    conn.autocommit = True
    cursor = conn.cursor()

    df = pd.DataFrame({"ra" : ra, "dec" : dec, "distance" : dist, "mag" : bmag,  "prob": prob,"distmu": distmu, "distsigma": distsigma, "distnorm": distnorm, "galname": galname, "contour": contour, "finalprob": finalprob[0]})
    df = df[['ra', 'dec','distance','mag','prob','distmu', 'distsigma','distnorm','galname','contour','finalprob']]
    df.to_csv("test2.csv", header=False, index=False)
    f = open('test2.csv')
    cursor.copy_from(f, str(table), sep=",")
    


    return
def initialise_database(database, table):
    """
    Create database.
    
    Parameters:
    -----------
    database  :  string
        name of database to save to
    table : string
        name of table to create
    

    Return:
    -------

    """
    conn = psycopg2.connect("dbname='%s' user='lana' host='localhost' password='myPass'" %database)
    conn.autocommit = True
    cursor = conn.cursor()
    cursor.execute('''CREATE TABLE %s ( ra FLOAT, dec FLOAT, distance FLOAT, mag FLOAT, prob FLOAT,distmu FLOAT, distsigma FLOAT, distnorm FLOAT, galname VARCHAR(40), contour FLOAT, finalprob FLOAT)''' % table)

    return
