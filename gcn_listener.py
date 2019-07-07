import astropy.utils.data
import sys
import lxml.etree
import gcn
import gcn.handlers
import gcn.notice_types
import healpy as hp
import numpy as np
import lxml.etree
import sys
import ligo.skymap.plot
import astropy_healpix
from astropy.io import fits
import numpy.ma as ma
from db import * 
from probability import *
from dataframe import * 
from contours import * 
import numpy as np
import time
import matplotlib
import pandas as pd
from astropy.io import fits
import psycopg2
import healpy as hp
import os
import json
import operator
from numpy.linalg import eig, inv
import functools
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.simbad import Simbad
import datetime

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
    prelim=params['AlertType']

    # Read the HEALPix sky map and the FITS header.
    skymap, header = hp.read_map(params['skymap_fits'],
                                 h=True, verbose=False)
    
    header = dict(header)
    hdulist = fits.open(params['skymap_fits'])
    distest = hdulist[1].header['DISTMEAN']
    diststd= hdulist[1].header['DISTSTD']

    prob, distmu, distsigma, distnorm = hp.read_map(params['skymap_fits'], field=[0, 1, 2, 3])
    npix = len(prob)
    nside=hp.npix2nside(npix)

    csm=integrated_probability(prob)
    contours=hpix_contours(csm,levels=[0.9,0.5,0.99],nest=False)

    levels=[0.9,0.5,0.99]
    levelsper=[90,50,99]
    
    for d in range(0,len(contours)):

        datafslong={}
        datafsshort={}
        datafslong99={}
        datafsshort99={}
        datafslong50={}
        datafsshort50={}
        
        numgals=[]
        numgals99=[]
        numgals50=[]

        tablenames=[]
        tablenames99=[]
        tablenames50=[]

        jsonlist=[]
        jsonlist2=[]
        jsonlist99=[]
        jsonlist299=[]
        jsonlist50=[]
        jsonlist250=[]

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
       

        database='new'
        table='ER4'+graceid+prelim+str(levelsper[d])+'probs230'
        tablenames.append(table)

        #make a database to put parameters for each galaxy into
        initialise_database(database, table)
        for i in range(0,len(split_ra2)):
            database='new'
            table='ER4'+graceid+prelim+str(levelsper[d])+'23%i' %i
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
                galaxy_ra, galaxy_dec, galaxy_dist, galaxy_Bmag, galaxyname=make_database_list(database, table, centreras, centredecs, maxdists, distest, diststd)
                
                
            else: 
                #if the contour doesn't cross 0/360 degrees, get centre coordinates and diameter of circle
                centrera,centredec, maxdist=radius(split_ra2[i], split_dec2[i],i)
                database='new'
                table='ER4'+graceid+prelim+str(levelsper[d])+'23%i' %i
                
                #query vizier using the circle and save in database for that contour
                galaxy_ra, galaxy_dec, galaxy_dist, galaxy_Bmag,galaxyname=make_database(database, table, centrera, centredec, maxdist, distest, diststd)
                
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
        finalprobs = calculate_absolute_probability(dist_incontourlist1, Bmag_incontourlist1,distsnorm_incontourlist1, mudists_incontourlist1, distssigma_incontourlist1,probs_incontourlist1)
        database='new'
        table='ER4'+graceid+prelim+str(levelsper[d])+'probs230'

        finalprobss=[]
        for j in range(0,len(finalprobs[0])):
            finalprobss.append(finalprobs[0,j])
        # make lists for dataframes
        finalprobss,ra_incontourlist,dec_incontourlist,finalprobslist,finalgalnamelist, dist_incontourlist,Bmag_incontourlist,contourlist1 = makelists(finalprobss,ra_incontourlist,ra_incontourlist1,dec_incontourlist,dec_incontourlist1,finalprobslist,finalgalnamelist,finalgalname,dist_incontourlist,dist_incontourlist1,Bmag_incontourlist,Bmag_incontourlist1,contourlist1,contourlist)
        
        #sort by descending probability
        Snumber=np.arange(0,len(finalprobslist),1)
        numgals.append(len(finalprobslist))
        fullprobs = dict(zip(Snumber, finalprobslist))
        finaldictsorted = sorted(fullprobs.items(), key=operator.itemgetter(1), reverse=True)
        

        database='new'
        table='ER4'+graceid+prelim+str(levelsper[d])+'probs230'
        
        #make/add to final database
        database_to_csv(database, table, ra_incontourlist1,dec_incontourlist1,dist_incontourlist1, Bmag_incontourlist1, probs_incontourlist1,mudists_incontourlist1,distssigma_incontourlist1,distsnorm_incontourlist1, finalgalname,finalprobs, contourlist1)
        
        #create dataframe for jsons
        dataf=create_dataframe(finaldictsorted, ra_incontourlist, dec_incontourlist, finalgalnamelist, dist_incontourlist, Bmag_incontourlist, contourlist)
        jsonlist.append(dataf[['Galaxy name', 'Galaxy probability', 'RA (degrees)', 'Dec (degrees)','Distance (Mpc)', 'B magntiude', 'Contour']].to_json())
        jsonlist2.append(dataf[['Galaxy name', 'Galaxy probability', 'RA (degrees)', 'Dec (degrees)','Distance (Mpc)', 'B magntiude','Contour']].to_csv())
       

        csv=str('[')
        

        for i in range(0,len(jsonlist)):
            if i<len(jsonlist)-1:
                csv= csv+ str(jsonlist[i])+','
            else:
                csv = csv + str(jsonlist[i])
        csv= csv+ ']'


        f = open("/users/lana/gxgwtest/static/"+graceid+prelim+str(levelsper[d])+".json", "w")
        f.write(csv)      # str() converts to string
        f.close()


        csv2=str('')
        for i in range(0,len(jsonlist2)):

            csv2= csv2+ str(jsonlist2[i])
        f = open("/users/lana/gxgwtest/static/"+graceid+prelim+str(levelsper[d])+".dat", "w")
        f.write( csv2 )      # str() converts to string
        f.close()

    
#payload = open('S190512at-2-Initial.xml', 'rb').read()
#root = lxml.etree.fromstring(payload)
#process_gcn(payload, root)


gcn.listen(handler=process_gcn)