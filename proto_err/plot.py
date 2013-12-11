## Everything to do with plotting
import math, datetime
import rpy2.robjects.lib.ggplot2 as ggplot2
import rpy2.robjects as ro
from rpy2.robjects.packages import importr

## Import some R modules
base = importr('base')
stats = importr('stats')
datasets = importr('datasets')