from ligo.skymap.io import read_sky_map
from astropy.io import fits
import healpy as hp
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


