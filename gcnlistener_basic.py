import sys
import lxml.etree
import gcn
import healpy as hp
import numpy as np
import astropy_healpix
import astropy.utils.data
from astropy.coordinates import SkyCoord
import astropy.units as u
import pandas as pd
from ligo.skymap.postprocess import crossmatch
from ligo.skymap.io import read_sky_map
import ligo.skymap.plot
from probability import calculate_absolute_probability,unique_galaxies,extract_LIGO_probability,makelists,sortbyprob
from dataframe import create_dataframe
from contours import contour_plots,join_0_360,split_contours,integrated_probability,hpix_contours
from createfiles import createtxt,createjsonfile,createasciifile
from GLADE import GLADEV2coordinates
from skymapio import readskymap


sys.setrecursionlimit(6000000)
# Function to call every time a GCN is received.
# Run only for notices of type
# LVC_PRELIMINARY, LVC_INITIAL, or LVC_UPDATE.
@gcn.handlers.include_notice_types(
    gcn.notice_types.LVC_PRELIMINARY,
    gcn.notice_types.LVC_INITIAL,
    gcn.notice_types.LVC_UPDATE)

def process_gcn(payload, root):
    
    # respond to only real 'observation' events.
    if root.attrib['role'] != 'observation':
        return

    # Read all of the VOEvent parameters from the "What" section.
    params = {elem.attrib['name']:
              elem.attrib['value']
              for elem in root.iterfind('.//Param')}
    
    
    # Respond only to 'CBC' events. Change 'CBC' to "Burst'
    # to respond to only unmodeled burst events.
    if params['Group'] != 'CBC':
        return

    graceid=params['GraceID']
    prelim=params['AlertType']+params['Pkt_Ser_Num']
    
    #read sky map parameters
    skymap, prob, distmu, distsigma, distnorm, npix, nside, distest, diststd = readskymap(params['skymap_fits'])
    
    #create integrated probability skymap
    csm=integrated_probability(prob)
    #calculate contours from csm
    contours=hpix_contours(csm,levels=[0.99],nest=False)
    
    #define 99% region
    levels=[0.99]
    levelsper=[99]

    #define distance limits
    distmax=distest+5*diststd
    distmin=distest-5*diststd
    
    #get coordinates of all GLADEV2 galaxies
    coordinates,data=GLADEV2coordinates(distmax,distmin)
    
    #crossmatch GLADE with multiorder skymap
    if 'v1' in params['skymap_fits']:
        version='.v1'
    if 'v0' in params['skymap_fits']:
        version='.v0'
    if 'v2' in params['skymap_fits']:
        version='.v2'
    else:
        version=''
    if ',0' in params['skymap_fits']:
        after=',0'
    if ',1' in params['skymap_fits']:
        after=',1'
    if ',2' in params['skymap_fits']:
        after=',2'
    if ',3' in params['skymap_fits']:
        after=',3'
    else:
        after=''

    
    if 'bayestar' in params['skymap_fits']:

        url='https://gracedb.ligo.org/api/superevents/'+graceid+'/files/bayestar.multiorder.fits'+after
    else:
        url='https://gracedb.ligo.org/api/superevents/'+graceid+'/files/LALInference'+version+'.multiorder.fits'+after
    skymap=read_sky_map(url, moc=True)
    result=crossmatch(skymap,coordinates)

    #for each contour region (Eg. 99%)
    for d in range(0,len(contours)):

        jsonlist=[]
        jsonlist2=[]
        tablenames=[]
        ra_incontourlist=[]
        contourlens=[]
        dec_incontourlist=[]
        finalprobslist=[]
        finalgalnamelist=[]
        dist_incontourlist=[]
        Bmag_incontourlist=[]
        ra_incontourlist1=[]
        dec_incontourlist1=[]
        probs_incontourlist1=[]
        probs_incontourlist=[]
        finalgalnamelist1=[]
        dist_incontourlist1=[]
        Bmag_incontourlist1=[]
        finalgalname=[]
        mudists_incontourlist1=[]
        distssigma_incontourlist1=[]
        distsnorm_incontourlist1=[]
        pdist_incontourlist1=[]
        Slum_incontourlist1=[]
        contourlist1=[]
        contourlist=[]
        contourss=[]
        ccc=[]
        
        #separate masked array into separate contours
        split_dec, split_ra = split_contours(contours, levels[d],d)

        #retrieve galaxies in 99 percent regions
        results=data[result.searched_prob<0.99]
        ra_incontour=results['RA'].values
        dec_incontour=results['Dec'].values
        dist_incontour=results['dist'].values
        Bmag_incontour=results['Bmag'].values
        name_incontour=results['HyperLEDA'].values
        
        # if the contour is split at 0/360 degrees, rejoin back together for plot
        split_ra2, split_dec2=join_0_360(split_ra, split_dec)
        
        #create a plot of contours and number them
        contour_plots(split_ra2, split_dec2,graceid, prelim, levelsper[d])
        contourss=np.ones(len(ra_incontour))
        
        # extract probability parameters at galaxy positions
       
        probs, mudists, distssigma, distsnorm = extract_LIGO_probability(ra_incontour, dec_incontour, nside, distsigma, prob, distnorm, distmu)
        
        # remove duplicates
        indices, finalgalname, ra_incontourlist1,dec_incontourlist1,dist_incontourlist1,Bmag_incontourlist1,probs_incontourlist1,mudists_incontourlist1,distssigma_incontourlist1,distsnorm_incontourlist1, contourlist=unique_galaxies(contourlist, contourss,ra_incontourlist1, ra_incontour, dec_incontourlist1, dec_incontour, dist_incontourlist1, dist_incontour, probs_incontourlist1, name_incontour, finalgalname, probs,Bmag_incontourlist1, Bmag_incontour,mudists_incontourlist1, mudists,distssigma_incontourlist1, distssigma,distsnorm_incontourlist1, distsnorm)
        
        # Calculate probability score
       
        finalprobs,pdist,Slum = calculate_absolute_probability(dist_incontourlist1, Bmag_incontourlist1, mudists_incontourlist1, distssigma_incontourlist1,distsnorm_incontourlist1,probs_incontourlist1)
        
        finalprobss=[]
        for j in range(0,len(finalprobs[0])):
            finalprobss.append(finalprobs[0,j])
        # make lists for dataframes
        
        finalprobss,ra_incontourlist,dec_incontourlist,finalprobslist,probs_incontourlist, finalgalnamelist, dist_incontourlist,Bmag_incontourlist,contourlist1,pdist_incontourlist1, Slum_incontourlist1 = makelists(finalprobss,ra_incontourlist,ra_incontourlist1,dec_incontourlist,dec_incontourlist1,finalprobslist, probs_incontourlist1, probs_incontourlist, finalgalnamelist,finalgalname,dist_incontourlist,dist_incontourlist1,Bmag_incontourlist,Bmag_incontourlist1,contourlist1,contourlist,pdist, pdist_incontourlist1, Slum, Slum_incontourlist1)
        
        #sort by descending probability
                                                                                         
        finaldictsorted, cumsumprobs=sortbyprob(finalprobslist)
        
        #create dataframe for jsons
        
        dataf=create_dataframe(finaldictsorted, ra_incontourlist, dec_incontourlist, probs_incontourlist, finalgalnamelist, dist_incontourlist, pdist_incontourlist1, Bmag_incontourlist, Slum_incontourlist1, contourlist,cumsumprobs)
        
        
        jsonlist.append(dataf[['Galaxy name', 'Galaxy probability score', 'RA (degrees)', 'Dec (degrees)','Location probability score','Distance (Mpc)', 'Distance probability score', 'B magnitude', 'B luminosity probability score','Cumulative Score']].to_json())
        jsonlist2.append(dataf[['Galaxy name', 'Galaxy probability score', 'RA (degrees)', 'Dec (degrees)','Location probability score','Distance (Mpc)', 'Distance probability score','B magnitude','B luminosity probability score','Cumulative Score']].to_csv())
        
        #createtxt(dataf,finalgalnamelist, finaldictsorted,graceid,prelim,levelsper,d,ccc)                                                                
        createjsonfile(jsonlist,graceid,prelim,levelsper,d)
        createasciifile(jsonlist2,graceid,prelim,levelsper,d)
  
          
#payload = open('S190910d-1-Preliminary.xml', 'rb').read()
#root = lxml.etree.fromstring(payload)
#process_gcn(payload, root)

gcn.listen(handler=process_gcn,  iamalive_timeout=300)
