import logging.config
import os

import yaml

with open(os.path.join(os.path.dirname(__file__), 'logging_conf.yaml'), 'rt') as f:
    config = yaml.safe_load(f.read())
logging.config.dictConfig(config)

class Semiring:
    REAL = 'real'
    LOG = 'log'
    TROPICAL = 'standard'

DEF_GAP_SYMBOL = '-'
MAX_INT32 = 2**31 - 1

from fstlib.cext.pywrapfst import *
from fstlib.cext.ops import DELTA, SHORTEST_DELTA
from fstlib.core import *
from fstlib.ext import *
from fstlib.paths import *
import fstlib.algos
import fstlib.factory
import fstlib.tools
