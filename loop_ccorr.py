import corrspline2
from astropy.io import fits
import ipdb
# need to figure how many orders there are in the image
# need to determine whether we are doing correlating many orders together
#     or if we are correlating all orders from target with a full spectrum

if __name__ == '__main__':
    ipdb.set_trace()
    # get number of orders
    target_fits = './blg0966red_multi.fits'
    hdu = fits.open(target_fits)
    num_orders = len(hdu[0].data[6,:,0])
    for order in range(num_orders):
        
