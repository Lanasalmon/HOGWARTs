import astropy.utils.data
import sys
sys.path.append('/home/swalsh/anaconda3/lib/python3.6/site-packages')
import lxml.etree
import gcn
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
import boto3, botocore
from email.mime.base import MIMEBase
from email import encoders
import smtplib
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
sys.path.append('/home/swalsh/anaconda3/lib/python3.6/site-packages')
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
from astroquery.vizier import VizierClass
from astropy.coordinates import SkyCoord
import astropy.units as u
from ligo.skymap.io import read_sky_map
from ligo.skymap.postprocess import crossmatch
import h5py

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
    print('I am listening')

    #data.to_hdf('test.h5','data',mode='w',format='table',data_columns=True)
    
    print(root.attrib['role'])
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
    
    hdulist = fits.open(params['skymap_fits'])
    distest = hdulist[1].header['DISTMEAN']
    diststd= hdulist[1].header['DISTSTD']
    

    skymap=read_sky_map(params['skymap_fits'], moc=False,distances=True)
    
    prob=skymap[0][0]
    distmu=skymap[0][1]
    distsigma=skymap[0][2]
    distnorm=skymap[0][3]
    npix = len(prob)
    nside=hp.npix2nside(npix)
   

    
    csm=integrated_probability(prob)
    contours=hpix_contours(csm,levels=[0.99],nest=False)
    
    levels=[0.99]
    levelsper=[99]

    distmax=distest+diststd
    distmin=distest-diststd

    data=pd.read_hdf('test2.h5','df')
   
    data.columns=['PGC','GWGC','HyperLEDA','2MASS','SDSS','flag1','RA','Dec','dist','dist_err','z','Bmag','a','b','c','d','e','f','g','h','i','j']


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
    url='https://gracedb.ligo.org/api/superevents/'+graceid+'/files/bayestar.multiorder.fits'
    skymap=read_sky_map(url, moc=True)

    result=crossmatch(skymap,coordinates)

    
    
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
        ccc=[]

        
        #separate masked array into separate contours
        split_dec, split_ra = split_contours(contours, levels[d],d)

        #retrieve pixels in 50, 90, 99 percent regions
        results=data[result.searched_prob_vol<0.99]
        ra_incontour=results['RA'].values
        dec_incontour=results['Dec'].values
        dist_incontour=results['dist'].values
        Bmag_incontour=results['Bmag'].values
        name_incontour=results['HyperLEDA'].values
        
       
        # if the contour is split at 0/360 degrees, rejoin back together
        split_ra2, split_dec2=join_0_360(split_ra, split_dec)
        
        #create a plot of contours and number them
        
        contour_plots(split_ra2, split_dec2,graceid, prelim, levelsper[d])


        ra_incontour,dec_incontour, dist_incontour, Bmag_incontour, name_incontour,contourss=checkifinellipsemoc(moc,nside, ipixes, data['RA'].values,data['Dec'].values,data['dist'].values, data['Bmag'].values, data['HyperLEDA'].values,1)
        contourss=np.ones(len(ra_incontour))

        #if the contour crosses 0/360 degrees, splut it up 
        #identify galaxies within localisation region
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
       
        finalprobs = calculate_absolute_probability(dist_incontourlist1, Bmag_incontourlist1, mudists_incontourlist1, distssigma_incontourlist1,distsnorm_incontourlist1,probs_incontourlist1)
        
        
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
        finalprobsorted=sorted(finalprobslist, reverse=True)
        
        
        #make/add to final database

        cumsumprobs=np.cumsum(finalprobsorted)
        
        #create dataframe for jsons
        
        dataf=create_dataframe(finaldictsorted, ra_incontourlist, dec_incontourlist, finalgalnamelist, dist_incontourlist, Bmag_incontourlist, contourlist,cumsumprobs)
        
        
        jsonlist.append(dataf[['Galaxy name', 'Galaxy probability', 'RA (degrees)', 'Dec (degrees)','Distance (Mpc)', 'B magnitude', 'Contour', 'Cumulative Probability']].to_json())
        jsonlist2.append(dataf[['Galaxy name', 'Galaxy probability', 'RA (degrees)', 'Dec (degrees)','Distance (Mpc)', 'B magnitude','Contour','Cumulative Probability']].to_csv())
        
        c=SkyCoord(ra=dataf['RA (degrees)'].values*u.degree, dec=dataf['Dec (degrees)'].values*u.degree)

        newrah=[c.ra.hms.h]
        newram=[c.ra.hms.m]
        newras=[c.ra.hms.s]
        newdecs=[c.dec]

        newc=[str(int(newrah[0][c]))+ ' ' +str(int(newram[0][c]))+' '+str(newras[0][c].round(2))+ ' '+str(newdecs[0][c]) for c,x in enumerate(newrah)]
        
        newc = [w.replace('d', ' ') for w in newc]
        newc = [w.replace('h', ' ') for w in newc]
        newc = [w.replace('m', ' ') for w in newc]
        newc = [w.replace('s', ' ') for w in newc]
        ccc.append(newc)


        csv=str('[')

        

        for i in range(0,len(jsonlist)):
            if i<len(jsonlist)-1:
                csv= csv+ str(jsonlist[i])+','
            else:
                csv = csv + str(jsonlist[i])
        csv= csv+ ']'

        
        
        f = open(graceid+prelim+str(levelsper[d])+".json", "w")
        f.write(csv)      # str() converts to string
        f.close()
        

        csv2=str('')
        for i in range(0,len(jsonlist2)):

            csv2= csv2+ str(jsonlist2[i])
        f = open(+graceid+prelim+str(levelsper[d])+".dat", "w")
        f.write( csv2 )      # str() converts to string
        f.close()
        

        
        for i in range(0,len(ccc)):
            coords=[]
            name=[]
            jtwo=[]
            probtxt=[]
            galname=[]
            exclam=[]
            for k in range(0,len(ccc[i])):
                
                coords.append(ccc[i][k])
                jtwo.append('J2000')
                name.append('COORD0'+str(k))
                exclam.append('!')
                galname.append(str(finalgalnamelist[finaldictsorted[k][0]]))
                probtxt.append(str(finaldictsorted[k][1]))
            datafc=pd.DataFrame(
                         {'Galaxy Name' : name,
                         'Coordinates': coords,
                         'J2000': jtwo,
                         '!': exclam,
                         'Name' :galname,
                         'Prob': probtxt
                         })
        
            datafc = datafc[['Galaxy Name', 'Coordinates', 'J2000','!','Name','Prob']]
            tfile = open(graceid+prelim+str(levelsper[d])+".txt", 'a')
            tfile.write(datafc.to_string(header=False, index=False))
            tfile.close()
            
#payload = open('S190910d-1-Preliminary.xml', 'rb').read()
#root = lxml.etree.fromstring(payload)
#process_gcn(payload, root)

gcn.listen(handler=process_gcn,  iamalive_timeout=300)
