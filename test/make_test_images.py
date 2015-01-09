import numpy as np
import pyfits

imsize = (2048, 4096)
np.random.seed(6563)

# Image 1, passed from step to step
data = 2*np.random.random(imsize)-1
hdu = pyfits.PrimaryHDU(data)
hdu.header.set('FORALL', True)
hdu.header.set('TREE', 'oak')
hdulist = pyfits.HDUList([hdu])
hdulist.writeto('im1.fits')

# Image 2, used by foo only
data = 2*np.random.random(imsize)-1
hdu = pyfits.PrimaryHDU(data)
hdu.header['FORFOO'] = True
hdu.header['BOGUS'] = 'nope'
hdulist = pyfits.HDUList([hdu])
hdulist.writeto('im2.fits')

# Image 3, used by bar only
data = 2*np.random.random(imsize)-1
hdu = pyfits.PrimaryHDU(data)
hdu.header['FORBAR'] = True
hdu.header['ISLAND'] = 'Hawaii'
hdulist = pyfits.HDUList([hdu])
hdulist.writeto('im3.fits')
    


    
    

