import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
import time
def createtxt(dataf, finalgalnamelist, finaldictsorted,graceid,prelim,levelsper,d,ccc):
    c=SkyCoord(ra=dataf['RA (degrees)'].values*u.degree, dec=dataf['Dec (degrees)'].values*u.degree)

    newrah=[c.ra.hms.h]
    newram=[c.ra.hms.m]
    newras=[c.ra.hms.s]
    newdecs=[c.dec]
    start=time.time()
    newc=[str(int(newrah[0][c]))+ ' ' +str(int(newram[0][c]))+' '+str(newras[0][c].round(2))+ ' '+str(newdecs[0][c]) for c,x in enumerate(newrah[0])]
    end=time.time()
    print('time',end-start)
    newc = [w.replace('d', ' ') for w in newc]
    newc = [w.replace('h', ' ') for w in newc]
    newc = [w.replace('m', ' ') for w in newc]
    newc = [w.replace('s', ' ') for w in newc]
    ccc.append(newc)
    coords=[x for x in ccc[0]]
    jtwo=['J2000']*len(ccc[0])
    name=['COORD'+str(k) for k in range(0,len(ccc[0]))]
    exclam=['!']*len(ccc[0])
    galname=[finalgalnamelist[finaldictsorted[d][0]] for d,x in enumerate(ccc[0])]
    probtxt=[finaldictsorted[d][1] for d,x in enumerate(ccc[0])]


        
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
def createjsonfile(jsonlist,graceid,prelim,levelsper,d):
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


def createasciifile(jsonlist2,graceid,prelim,levelsper,d):

    csv2=str('')
    for i in range(0,len(jsonlist2)):

        csv2= csv2+ str(jsonlist2[i])
    f = open(graceid+prelim+str(levelsper[d])+".dat", "w")
    f.write( csv2 )      # str() converts to string
    f.close()
    return
