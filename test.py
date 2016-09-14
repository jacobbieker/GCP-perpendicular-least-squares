__author__ = 'jacob'
from astropy.table import Table
fits_table = Table.read("rxj1226allfit.fits")
for key, value in fits_table.meta.items():
    print('{0} = {1}'.format(key, value))
#for value in fits_table:
 #   print(value)