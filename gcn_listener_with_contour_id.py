import sys
import ligo.skymap.plot
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
from probability import calculate_absolute_probability,unique_galaxies,extract_LIGO_probability,makelists,sortbyprob
from dataframe import create_dataframe
from contours import contour_plots,join_0_360,split_contours,integrated_probability,hpix_contours,pixels_in_region,checkifinellipsemoc
from createfiles import createtxt,createjsonfile,createasciifile
from GLADE import Viziergalaxies
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

    # Read the HEALPix sky map and the FITS header.
    skymap, prob, distmu, distsigma, distnorm, npix, nside, distest, diststd=readskymap(params['skymap_fits'])
    
    #create integrated probability map and contours
    csm=integrated_probability(prob)
    contours=hpix_contours(csm,levels=[0.99],nest=False)

    #define percentage localisation regions
    levels=[0.99]
    levelsper=[99]
    
    #for each percentage localisation region
    for d in range(0,len(contours)):
        jsonlist=[]
        jsonlist2=[]
        ccc=[]
        
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
        finalgalnamelist1=[]
        dist_incontourlist1=[]
        Bmag_incontourlist1=[]
        finalgalname=[]
        mudists_incontourlist1=[]
        distssigma_incontourlist1=[]
        distsnorm_incontourlist1=[]
        contourlist1=[]
        contourlist=[]
        contourss=[]
        
        #separate masked array into separate contours
        split_dec, split_ra = split_contours(contours, levels[d],d)
        
        #retrieve pixels in 50, 90, 99 percent regions
        moc,ipixes=pixels_in_region(levels[d],params['skymap_fits'])
        
        # if the contour is split at 0/360 degrees, rejoin back together
        split_ra2, split_dec2=join_0_360(split_ra, split_dec)

        #create a plot of contours and number them
        contour_plots(split_ra2, split_dec2,graceid, prelim, levelsper[d])
       
        #for each separate contour
        for i in range(0,len(split_ra2)):
            #identify galaxies in circular regions surrounding contour
            galaxy_ra, galaxy_dec, galaxy_dist, galaxy_Bmag,galaxyname=Viziergalaxies(split_ra2,split_dec2,i,distest,diststd)  
            
            #identify galaxies within localisation region
            ra_incontour, dec_incontour, dist_incontour, Bmag_incontour, name_incontour,contourss=checkifinellipsemoc(moc,nside, ipixes, galaxy_ra, galaxy_dec,galaxy_dist, galaxy_Bmag, galaxyname,i+1)
            contourlens.append(len(ra_incontour))
            
            # extract probability parameters at galaxy positions
            probs, mudists, distssigma, distsnorm = extract_LIGO_probability(ra_incontour, dec_incontour, nside, distsigma, prob, distnorm, distmu)
            
            # remove duplicates
            indices, finalgalname, ra_incontourlist1,dec_incontourlist1,dist_incontourlist1,Bmag_incontourlist1,probs_incontourlist1,mudists_incontourlist1,distssigma_incontourlist1,distsnorm_incontourlist1, contourlist=unique_galaxies(contourlist, contourss,ra_incontourlist1, ra_incontour, dec_incontourlist1, dec_incontour, dist_incontourlist1, dist_incontour, probs_incontourlist1, name_incontour, finalgalname, probs,Bmag_incontourlist1, Bmag_incontour,mudists_incontourlist1, mudists,distssigma_incontourlist1, distssigma,distsnorm_incontourlist1, distsnorm)
           
            r = 10* u.arcminute
            if len(ra_incontour)==0:
                    print('nothing happened')
                    continue
          
        # Calculate probability score
        finalprobs = calculate_absolute_probability(dist_incontourlist1, Bmag_incontourlist1,mudists_incontourlist1, distssigma_incontourlist1,distsnorm_incontourlist1,probs_incontourlist1)
        
        finalprobss=[]
        for j in range(0,len(finalprobs[0])):
            finalprobss.append(finalprobs[0,j])
        # make lists for dataframes
        finalprobss,ra_incontourlist,dec_incontourlist,finalprobslist,finalgalnamelist, dist_incontourlist,Bmag_incontourlist,contourlist1 = makelists(finalprobss,ra_incontourlist,ra_incontourlist1,dec_incontourlist,dec_incontourlist1,finalprobslist,finalgalnamelist,finalgalname,dist_incontourlist,dist_incontourlist1,Bmag_incontourlist,Bmag_incontourlist1,contourlist1,contourlist)
        
        #sort by descending probability
        finaldictsorted, cumsumprobs=sortbyprob(finalprobslist)    
        
        #create dataframe for jsons
        
        dataf=create_dataframe(finaldictsorted, ra_incontourlist, dec_incontourlist, finalgalnamelist, dist_incontourlist, Bmag_incontourlist, contourlist,cumsumprobs)
        jsonlist.append(dataf[['Galaxy name', 'Galaxy probability', 'RA (degrees)', 'Dec (degrees)','Distance (Mpc)', 'B magnitude', 'Contour']].to_json())
        jsonlist2.append(dataf[['Galaxy name', 'Galaxy probability', 'RA (degrees)', 'Dec (degrees)','Distance (Mpc)', 'B magnitude','Contour']].to_csv())
        
        createtxt(dataf,finalgalnamelist, finaldictsorted,graceid,prelim,levelsper,d,ccc)                                                                           
        createjsonfile(jsonlist,graceid,prelim,levelsper,d)
        createasciifile(jsonlist2,graceid,prelim,levelsper,d)
    
#payload = open('S190910d-1-Preliminary.xml', 'rb').read()
#root = lxml.etree.fromstring(payload)
#process_gcn(payload, root)


gcn.listen(handler=process_gcn)
