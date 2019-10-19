import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
def createtxt(dataf, finalgalnamelist, finaldictsorted,graceid,prelim,levelsper,d,ccc):
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