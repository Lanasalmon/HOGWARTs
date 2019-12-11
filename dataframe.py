from astropy.coordinates import SkyCoord
from astroquery.simbad import Simbad
import astropy.units as u
import pandas as pd
import json

def create_dataframe(finaldictsorted, ra, dec,probs,name, dist, pdist,bmag, Slum,contour,cumsumprobs):
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

    finalgalname=[name[finaldictsorted[i][0]] for i,x in enumerate(finaldictsorted)]
    finalprob=[finaldictsorted[i][1] for i,x in enumerate(finaldictsorted)]
    finalra=[ra[finaldictsorted[i][0]] for i,x in enumerate(finaldictsorted)]
    finaldec=[dec[finaldictsorted[i][0]] for i,x in enumerate(finaldictsorted)]
    finaldist=[dist[finaldictsorted[i][0]] for i,x in enumerate(finaldictsorted)]
    finalbmag=[bmag[finaldictsorted[i][0]] for i,x in enumerate(finaldictsorted)]
    finalcontour=[contour[finaldictsorted[i][0]] for i,x in enumerate(finaldictsorted)]
    finalSlum=[Slum[finaldictsorted[i][0]] for i,x in enumerate(finaldictsorted)]
    finalpdist=[pdist[finaldictsorted[i][0]] for i,x in enumerate(finaldictsorted)]
    cumsum=[cumsumprobs[i] for i,x in enumerate(finaldictsorted)]
    finalorigprob=[probs[finaldictsorted[i][0]] for i,x in enumerate(finaldictsorted)]
    
        
    
    dataf = pd.DataFrame(
                     {'Galaxy name' : finalgalname,
                     'Galaxy probability score': finalprob,
                     'RA (degrees)': finalra,
                     'Dec (degrees)': finaldec,
                     'Location probability score' : finalorigprob,
                     'Galaxy name': finalgalname,
                     'Distance (Mpc)': finaldist,
                     'Distance probability score': finalpdist,
                     'B magnitude': finalbmag,
                     'B luminosity probability score': finalSlum,
                     'Contour': finalcontour,
                     'Cumulative Score':cumsum
                     
                     })
    
    dataf = dataf[['Galaxy name', 'Galaxy probability score', 'RA (degrees)', 'Dec (degrees)','Location probability score','Distance (Mpc)', 'Distance probability score', 'B magnitude', 'B luminosity probability score','Contour', 'Cumulative Score']]
    dataf.to_json('output.json')
    return dataf







