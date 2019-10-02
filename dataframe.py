from astropy.coordinates import SkyCoord
from astroquery.simbad import Simbad
import astropy.units as u
import pandas as pd
import json

def create_dataframe(finaldictsorted, ra, dec,name, dist, bmag,contour,cumsumprobs):
    """
    Create sorted dataframe of galaxies.
    
    Parameters:
    -----------
    
    ra, dec : list
        list of galaxy coordinates within contours
    dist : list
        list of galaxy distances within contours
    bmag : list
        list of galaxy B magnitudes within contours
    
    name: list 
        list of galaxy names within contours
    finaldictsorted: dict
        dictionary of sorted probability values
    contour: list 
        list of contours corresponding to galaxies
    
    
    
    Return:
    -------
    
    """
    finalgalname=[]
    finalprob=[]
    finalra=[]
    finaldec=[]
    timetoset=[]
    finaldist=[]
    finalbmag=[]
    finalcontour=[]
    cumsum=[]
    for i in range(0,len(finaldictsorted)):
        r = 10* u.arcminute
        finalgalname.append(name[finaldictsorted[i][0]])
        finalprob.append(finaldictsorted[i][1])
        finalra.append(ra[finaldictsorted[i][0]])
        finaldec.append(dec[finaldictsorted[i][0]])
        finaldist.append(dist[finaldictsorted[i][0]])
        finalbmag.append(bmag[finaldictsorted[i][0]])
        finalcontour.append(contour[finaldictsorted[i][0]])
        cumsum.append(cumsumprobs[i])
        
    
    dataf = pd.DataFrame(
                     {'Galaxy name' : finalgalname,
                     'Galaxy probability': finalprob,
                     'RA (degrees)': finalra,
                     'Dec (degrees)': finaldec,
                     'Galaxy name': finalgalname,
                     'Distance (Mpc)': finaldist,
                     'B magnitude': finalbmag,
                     'Contour': finalcontour,
                     'Cumulative Probability':cumsum
                     
                     })
    
    dataf = dataf[['Galaxy name', 'Galaxy probability', 'RA (degrees)', 'Dec (degrees)','Distance (Mpc)', 'B magnitude', 'Contour', 'Cumulative Probability']]
    dataf.to_json('output.json')
    return dataf
