import numpy as np
import fitsio

imsize = (2048, 4096)
np.random.seed(6563)

# Image 1, passed from step to step
data = 2*np.random.random(imsize)-1
hdr_dict = {'FORALL': True,
            'TREE': 'oak'}
fitsio.write('im1.fits', data, header=hdr_dict, clobber=True)

# Image 2, used by foo only
data = 2*np.random.random(imsize)-1
hdr_dict = {'FORFOO': True,
            'BOGUS': 'nope'}
fitsio.write('im2.fits', data, header=hdr_dict, clobber=True)

# Image 3, used by bar only
data = 2*np.random.random(imsize)-1
hdr_dict = {'FORFOO': True,
            'ISLAND': 'Hawaii'}
fitsio.write('im3.fits', data, header=hdr_dict, clobber=True)
    


    
    

