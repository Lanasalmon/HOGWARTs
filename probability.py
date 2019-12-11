import functools
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
import numpy as np
import healpy as hp
from numpy.linalg import eig, inv
import operator


def unique_galaxies(contourlist, contourss,ra_incontourlist1, ra_incontour, dec_incontourlist1, dec_incontour, dist_incontourlist1, dist_incontour, probs_incontourlist1, name_incontour, finalgalname, probs,Bmag_incontourlist1, Bmag_incontour,mudists_incontourlist1, mudists,distssigma_incontourlist1, distssigma,distsnorm_incontourlist1, distsnorm):
    """
    Create lists of unique galaxies in contour.
    
    Parameters:
    -----------
    
    ra_incontourlist1: list
        empty list
    ra_incontour : list
        list of galaxy RAs within this contour
    dec_incontourlist1: list
        empty list
    dec_incontour : list
        list of galaxy declinations within this contour
    dist_incontourlist1: list
        empty list
    dist_incontour : list
        list of galaxy distances within this contour
    bmag_incontourlist1: list
        empty list
    bmag_incontour : list
        list of galaxy B magnitudes within this contour
    probs_incontourlist1: list
        empty list
    probs_incontour : list
        list of galaxy probabilities within this contour
    name_incontourlist1: list
        empty list
    finalgalname : list
        list of galaxy names within this contour
    mudists_incontourlist1: list
        empty list
    mudists : list
        list of galaxy mean distances within this contour
    distssigma_incontourlist1: list
        empty list
    distssigma : list
        list of galaxy standard deviation on distances within this contour
    distnorm_incontourlist1: list
        empty list
    distsnorm : list
        list of galaxy normalisation factors within this contour
    
    contourlist: list
        empty list
    contours : list
        list of contours corresponding to galaxies within this contour
    
    
    Return:
    -------
    indices : list
        list of indices of unique galaxies.

    ra_incontourlist1: list
        list of unique galaxy RAs within this contours,
   
    dec_incontourlist1: list
        list of unique galaxy declinations within all contours,

    dist_incontourlist1: list
        list of unique galaxy distances within all contours,
    
    bmag_incontourlist1: list
        list of unique galaxy b magnitudes within all contours,
    
    probs_incontourlist1: list
        list of unique galaxy probabilities within all contours,
    name_incontourlist1: list
        list of unique galaxy names within all contours,
    
    mudists_incontourlist1: list
        list of unique galaxy mean distances within all contours,
   
    distssigma_incontourlist1: list
        list of unique galaxy standard deviation on distances within all contours,
    
    distnorm_incontourlist1: list
        list of unique galaxy normalisation factors within all contours,
    
    contourlist: list
        list of  ontours corresponding to unique galaxies within all contours
    
    """
    ra_incontourlist1=np.append(ra_incontourlist1,ra_incontour)
    ra_incontourlist2, indices = np.unique(ra_incontourlist1, return_index=True)
    ra_incontourlist1= [ra_incontourlist1[i] for i in indices]
    dec_incontourlist1=np.append(dec_incontourlist1,dec_incontour)
    dec_incontourlist1= [dec_incontourlist1[i] for i in indices]
    dist_incontourlist1=np.append(dist_incontourlist1,dist_incontour)
    dist_incontourlist1= [dist_incontourlist1[i] for i in indices]
    Bmag_incontourlist1=np.append(Bmag_incontourlist1,Bmag_incontour)
    Bmag_incontourlist1= [Bmag_incontourlist1[i] for i in indices]
    probs_incontourlist1=np.append(probs_incontourlist1,probs)
    probs_incontourlist1= [probs_incontourlist1[i] for i in indices]
    mudists_incontourlist1=np.append(mudists_incontourlist1,mudists)
    mudists_incontourlist1= [mudists_incontourlist1[i] for i in indices]
    distssigma_incontourlist1=np.append(distssigma_incontourlist1,distssigma)
    distssigma_incontourlist1= [distssigma_incontourlist1[i] for i in indices]
    distsnorm_incontourlist1=np.append(distsnorm_incontourlist1,distsnorm)
    distsnorm_incontourlist1= [distsnorm_incontourlist1[i] for i in indices]
    finalgalname=np.append(finalgalname,name_incontour)
           
    finalgalname= [finalgalname[i] for i in indices]
    contourlist=np.append(contourlist,contourss)
           
    contourlist= [contourlist[i] for i in indices]
    return indices, finalgalname,ra_incontourlist1,dec_incontourlist1,dist_incontourlist1,Bmag_incontourlist1,probs_incontourlist1,mudists_incontourlist1,distssigma_incontourlist1,distsnorm_incontourlist1, contourlist

def makelists(finalprobss,ra_incontourlist,ra_incontourlist1,dec_incontourlist,dec_incontourlist1,probs, probs_incontourlist, finalprobslist,finalgalnamelist,finalgalname,dist_incontourlist,dist_incontourlist1,Bmag_incontourlist,Bmag_incontourlist1,contourlist,contourss, pdist, pdistlist, Slum, Slumlist):
    """
    Combine all parameters from all contours together into single lists

    Parameters:
    -----------

    ra_incontourlist: list
        list of galaxy RAs within all contours,
   
    dec_incontourlist1: list
        list of galaxy declinations within all contours,

    dist_incontourlist1: list
        list of galaxy distances within all contours,
    
    bmag_incontourlist1: list
        list of galaxy b magnitudes within all contours,
    
    probs_incontourlist1: list
        list of galaxy probabilities within all contours,

    coord_incontourlist1: list
        list of galaxy SkyCoord coordinates within all contours,

    name_incontourlist1: list
        list of galaxy names within all contours,
    
    mudists_incontourlist1: list
        list of galaxy mean distances within all contours,
   
    distssigma_incontourlist1: list
        list of galaxy standard deviation on distances within all contours,
    
    distnorm_incontourlist1: list
        list of galaxy normalisation factors within all contours,
    
    contourlist: list
        list of contours corresponding to galaxies within all contours
    
    Return:
    -------
    indices : list
        list of indices of unique galaxies.

    ra_incontourlist1: list
        list of galaxy RAs within all contours,
   
    dec_incontourlist1: list
        list of galaxy declinations within all contours,

    dist_incontourlist1: list
        list of galaxy distances within all contours,
    
    bmag_incontourlist1: list
        list of galaxy b magnitudes within all contours,
    
    probs_incontourlist1: list
        list of galaxy probabilities within all contours,

    coord_incontourlist1: list
        list of galaxy SkyCoord coordinates within all contours,

    name_incontourlist1: list
        list of galaxy names within all contours,
    
    mudists_incontourlist1: list
        list of galaxy mean distances within all contours,
   
    distssigma_incontourlist1: list
        list of galaxy standard deviation on distances within all contours,
    
    distnorm_incontourlist1: list
        list of galaxy normalisation factors within all contours,
    
    contourlist: list
        list of contours corresponding to galaxies within all contours
    
    """
    finalprobss=list(finalprobss)
    ra_incontourlist=np.append(ra_incontourlist,ra_incontourlist1)
    dec_incontourlist=np.append(dec_incontourlist,dec_incontourlist1)
    finalprobslist=np.append(finalprobslist,finalprobss)
    finalgalnamelist=np.append(finalgalnamelist,finalgalname)
    dist_incontourlist=np.append(dist_incontourlist,dist_incontourlist1)
    Bmag_incontourlist=np.append(Bmag_incontourlist,Bmag_incontourlist1)
    contourlist=np.append(contourlist,contourss)
    pdistlist=np.append(pdistlist,pdist)
    Slumlist=np.append(Slumlist,Slum)
    probs_incontourlist=np.append(probs_incontourlist,probs)
    return finalprobss,ra_incontourlist,dec_incontourlist,finalprobslist,probs_incontourlist,finalgalnamelist, dist_incontourlist,Bmag_incontourlist,contourlist, pdistlist, Slumlist


def extract_LIGO_probability(ra, dec, nside, distsigma, prob, distnorm, distmu):
    """
    Identify and save LIGO/Virgo probabilities at galaxy positions.
    
    Parameters:
    -----------
    ra, dec : list
        list of galaxy coordinates within contour
    distmu: list 
        list of skymap mean distances
    distssigma: list 
        list of skymap distance standard deviations
    distnorm: list 
        list of skymap normalisation factors 
    
    prob: list 
        list of LIGO/Virgo galaxy probabilities 
    nside: int
        resolution of skymap
    
    
    
    Return:
    ----
    pixel_prob : list 
        list of LIGO/Virgo galaxy probabilities at galaxy positions
    pixel_mudist : list 
        list of LIGO/Virgo galaxy mean distances at galaxy positions
    pixel_distsigma : list 
        list of LIGO/Virgo galaxy standard deviations on distance at galaxy positions
    pixel_distnorm : list 
        list of LIGO/Virgo galaxy normalisation factors at galaxy positions
    
    """
    pixel_prob = []
    pixel_mudist = []
    pixel_distsigma = []
    pixel_distnorm = []
    theta = 0.5 * np.pi - np.deg2rad(dec)
    phi = np.deg2rad(ra)
    ipixes = hp.ang2pix(nside, theta, phi)
    pixel_prob=prob[ipixes]
    pixel_mudist=distmu[ipixes]
    pixel_distsigma=distsigma[ipixes]
    pixel_distnorm=distnorm[ipixes]

    return pixel_prob, pixel_mudist, pixel_distsigma, pixel_distnorm

def sortbyprob(finalprobslist):

    Snumber=np.arange(0,len(finalprobslist),1)
    fullprobs = dict(zip(Snumber, finalprobslist))
    finaldictsorted = sorted(fullprobs.items(), key=operator.itemgetter(1), reverse=True)
    finalprobsorted=sorted(finalprobslist, reverse=True)
    cumsumprobs=np.cumsum(finalprobsorted)
    return finaldictsorted, cumsumprobs


def calculate_absolute_probability(distance, Bmag, mudist, distsigma, distnorm, probs):
    """
    Calculate probability score for each galaxy
    
    Parameters:
    -----------
    probs : list 
        list of LIGO/Virgo galaxy probabilities at galaxy positions
    mudist : list 
        list of LIGO/Virgo galaxy mean distances at galaxy positions
    distsigma : list 
        list of LIGO/Virgo galaxy standard deviations on distance at galaxy positions
    distnorm : list 
        list of LIGO/Virgo galaxy normalisation factors at galaxy positions
    
    
    distance: list 
        list of galaxy distances 
    Bmag: list
        list of galaxy B magnitudes
    
    
    
    Return:
    ----
    pixel_prob : list 
        list of LIGO/Virgo galaxy probabilities at galaxy positions
    pixel_mudist : list 
        list of LIGO/Virgo galaxy mean distances at galaxy positions
    pixel_distsigma : list 
        list of LIGO/Virgo galaxy standard deviations on distance at galaxy positions
    pixel_distnorm : list 
        list of LIGO/Virgo galaxy normalisation factors at galaxy positions
    
    """

    Lsun = 3.86e26
    Msun = 4.83
    Mknmin = -19
    Mknmax = -12
    Lblist = []
    for i in range(0, len(distance)):
        Mb = Bmag[i] - 5 * np.log10((distance[i] * 10 ** 6)) + 5
        L = Lsun * 2.512 ** (Msun - Mb)
        Lb = Lsun * 2.512 ** (Msun - Mb)
        Lblist.append(Lb)
        # if distance[i]==1432.19 and Bmag[i]==14.185:
        #     print('Lb',Lb)
        # if distance[i]==1393.16 and Bmag[i]==17.67:
        #     print('Lb2',Lb)
        # if distance[i]==1445.44 and Bmag[i]==18.442:
        #     print('Lb3',Lb)

    targets={}

    Sdet=[]
    Lsun = 3.86e26
    Msun = 4.83
    Mknmin = -17
    Mknmax = -12
    Slist = []
    pdist=np.array(distnorm) * np.exp(-((np.array(distance) - np.array(mudist)) ** 2) / (2 * np.array(distsigma) ** 2))
    Sloc=(np.array(probs)*np.array(pdist))
    Slum=np.array(Lblist)/np.sum(Lblist)
    


    SS=Sloc * Slum

    S = Sloc * Slum 
    
    # for s in range(0,len(pdist)):
    #     if distance[s]==1432.19 and Bmag[s]==14.185:
    #         print('pdist',pdist[s], Sloc[s], Slum[s],S[s])
    #     if distance[s]==1393.16 and Bmag[s]==17.67:
    #         print('pdist2',pdist[s],Sloc[s],Slum[s],S[s])
    #     if distance[s]==1445.44 and Bmag[s]==18.442:
    #         print('pdist3',pdist[s],Sloc[s],Slum[s],S[s])    
    Slist.append(S)
    

    print(sum(Slist/np.sum(Slist)))
   
    return Slist/np.sum(Slist), pdist, Slum
