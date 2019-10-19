
import sys
import os
import json
import datetime
import lxml.etree
import gcn
import lxml.etree
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
import boto3, botocore
from email.mime.base import MIMEBase
from email import encoders
import smtplib
import healpy as hp
import numpy as np
import ligo.skymap.plot
import astropy_healpix
from astropy.io import fits
import astropy.utils.data
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.simbad import Simbad
from astroquery.vizier import VizierClass
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy.ma as ma
import time
import matplotlib
import pandas as pd
import psycopg2
import operator
from numpy.linalg import eig, inv
import functools
from ligo.skymap.io import read_sky_map
from ligo.skymap.postprocess import crossmatch
import h5py
from db import * 
from probability import *
from dataframe import * 
from contours import *

sys.setrecursionlimit(6000000)
# Function to call every time a GCN is received.
# Run only for notices of type
# LVC_PRELIMINARY, LVC_INITIAL, or LVC_UPDATE.
@gcn.handlers.include_notice_types(
    gcn.notice_types.LVC_PRELIMINARY,
    gcn.notice_types.LVC_INITIAL,
    gcn.notice_types.LVC_UPDATE)

def createtxtfile(dataf):
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
    return   
def GLADEV2coordinates():

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
    return coordinates
def readskymap(skymap):
        
    # Read the HEALPix sky map and the FITS header.

    hdulist = fits.open(skymap)
    distest = hdulist[1].header['DISTMEAN']
    diststd= hdulist[1].header['DISTSTD']


    skymap=read_sky_map(skymap, moc=False,distances=True)

    prob=skymap[0][0]
    distmu=skymap[0][1]
    distsigma=skymap[0][2]
    distnorm=skymap[0][3]
    npix = len(prob)
    nside=hp.npix2nside(npix)
    return skymap, prob, distmu, distsigma, distnorm, npix, nside, distest, diststd

def createjsonfile(jsonlist):
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
    return


def createasciifile(jsonlist2):

    csv2=str('')
    for i in range(0,len(jsonlist2)):

        csv2= csv2+ str(jsonlist2[i])
    f = open(graceid+prelim+str(levelsper[d])+".dat", "w")
    f.write( csv2 )      # str() converts to string
    f.close()
    return

def sortbyprob(finalprobslist):

    Snumber=np.arange(0,len(finalprobslist),1)
    fullprobs = dict(zip(Snumber, finalprobslist))
    finaldictsorted = sorted(fullprobs.items(), key=operator.itemgetter(1), reverse=True)
    finalprobsorted=sorted(finalprobslist, reverse=True)
    cumsumprobs=np.cumsum(finalprobsorted)
    return finaldictsorted, cumsumprobs

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
    distmax=distest+3*diststd
    distmin=distest-3*diststd
    
    #get coordinates of all GLADEV2 galaxies
    coordinates=GLADEV2coordinates()
    
    #crossmatch GLADE with multiorder skymap
    url='https://gracedb.ligo.org/api/superevents/'+graceid+'/files/bayestar.multiorder.fits'
    skymap=read_sky_map(url, moc=True)
    result=crossmatch(skymap,coordinates)

    #for each contour region (Eg. 99%)
    for d in range(0,len(contours)):

        jsonlist=[]
        jsonlist2=[]
        
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

        #retrieve galaxies in 99 percent regions
        results=data[result.searched_prob_vol<0.99]
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
       
        finalprobs = calculate_absolute_probability(dist_incontourlist1, Bmag_incontourlist1, mudists_incontourlist1, distssigma_incontourlist1,distsnorm_incontourlist1,probs_incontourlist1)
        
        finalprobss=[]
        for j in range(0,len(finalprobs[0])):
            finalprobss.append(finalprobs[0,j])
        # make lists for dataframes
        
        finalprobss,ra_incontourlist,dec_incontourlist,finalprobslist,finalgalnamelist, dist_incontourlist,Bmag_incontourlist,contourlist1 = makelists(finalprobss,ra_incontourlist,ra_incontourlist1,dec_incontourlist,dec_incontourlist1,finalprobslist,finalgalnamelist,finalgalname,dist_incontourlist,dist_incontourlist1,Bmag_incontourlist,Bmag_incontourlist1,contourlist1,contourlist)
        
        #sort by descending probability
                                                                                         
        finaldictsorted, cumsumprobs=sortbyprob(finalprobslist)    
        
        #create dataframe for jsons
        
        dataf=create_dataframe(finaldictsorted, ra_incontourlist, dec_incontourlist, finalgalnamelist, dist_incontourlist, Bmag_incontourlist, contourlist,cumsumprobs)
        
        #create files
        jsonlist.append(dataf[['Galaxy name', 'Galaxy probability', 'RA (degrees)', 'Dec (degrees)','Distance (Mpc)', 'B magnitude', 'Cumulative Probability']].to_json())
        jsonlist2.append(dataf[['Galaxy name', 'Galaxy probability', 'RA (degrees)', 'Dec (degrees)','Distance (Mpc)', 'B magnitude','Cumulative Probability']].to_csv())
        createtxtfile(dataf)                                                                           
        createjsonfile(jsonlist)
        createasciifile(jsonlist2)
        return
  
          
#payload = open('S190910d-1-Preliminary.xml', 'rb').read()
#root = lxml.etree.fromstring(payload)
#process_gcn(payload, root)

gcn.listen(handler=process_gcn,  iamalive_timeout=300)
