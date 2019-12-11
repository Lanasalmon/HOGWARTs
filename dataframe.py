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
    finalgalname=[]
    finalprob=[]
    finalra=[]
    finaldec=[]
    timetoset=[]
    finaldist=[]
    finalbmag=[]
    finalcontour=[]
    finalpdist=[]
    finalSlum=[]
    finalorigprob=[]
    cumsum=[]
    # finalgalname=[name[finaldictsorted[i][0]] for i,x in finaldictsorted]
    # finalprob=[finaldictsorted[i][1] for i,x in finaldictsorted]
    # finalra=[ra[finaldictsorted[i][0]] for i,x in finaldictsorted]
    # finaldec=[dec[finaldictsorted[i][0]] for i,x in finaldictsorted]
    # finaldist=[dist[finaldictsorted[i][0]] for i,x in finaldictsorted]
    # finalbmag=[bmag[finaldictsorted[i][0]] for i,x in finaldictsorted]
    # finalcontour=[contour[finaldictsorted[i][0]] for i,x in finaldictsorted]
    # finalSlum=[Slum[finaldictsorted[i][0]] for i,x in finaldictsorted]
    # finalpdist=[pdist[finaldictsorted[i][0]] for i,x in finaldictsorted]
    # cumsum=[cumsumprobs[i] for i,x in finaldictsorted]
    # finalorigprob=[probs[finaldictsorted[i][0]] for i,x in finaldictsorted]
    for i in range(0,len(finaldictsorted)):
        r = 10* u.arcminute
        finalgalname.append(name[finaldictsorted[i][0]])
        finalprob.append(finaldictsorted[i][1])
        finalra.append(ra[finaldictsorted[i][0]])
        finaldec.append(dec[finaldictsorted[i][0]])
        finaldist.append(dist[finaldictsorted[i][0]])
        finalbmag.append(bmag[finaldictsorted[i][0]])
        finalcontour.append(contour[finaldictsorted[i][0]])
        finalSlum.append(Slum[finaldictsorted[i][0]])
        finalpdist.append(pdist[finaldictsorted[i][0]])
        cumsum.append(cumsumprobs[i])
        finalorigprob.append(probs[finaldictsorted[i][0]])
        
    
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





