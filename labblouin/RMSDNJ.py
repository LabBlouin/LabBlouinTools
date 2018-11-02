from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
import pandas as pd
from skbio import tree, DistanceMatrix
import numpy as np
import sys

m = pd.read_csv(sys.argv[1])
m[m.isnull()] = 0
arr = m.as_matrix()
M = arr + arr.T
dm = DistanceMatrix(M)
tree = tree.nj(dm)
