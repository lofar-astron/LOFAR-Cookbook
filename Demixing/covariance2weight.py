#!/usr/bin/python

import pyrap.tables
import numpy
import sys

msname = sys.argv[1]
if len(sys.argv) == 3 :
   correctedcol = sys.argc[2]
else:   
   correctedcol = "CORRECTED_DATA"

t = pyrap.tables.table(msname, readonly = False )
covarcolname = t.getcolkeyword( correctedcol, 'LOFAR_COVARIANCE_COLUMN')
for t1 in t.iter('TIME'):
   covar = t1.getcol(covarcolname)
   weight = 1.0/covar.diagonal(axis1=2, axis2=3).sum(axis=1)
   weight[numpy.nonzero(numpy.isnan(weight))] = 0.0
   t1.putcol('WEIGHT', weight)

t.close()
