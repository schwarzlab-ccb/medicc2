import logging.config
import os

import yaml

import medicc.bootstrap
import medicc.core
import medicc.io
import medicc.nj
import medicc.plot
import medicc.sim
import medicc.stats
import medicc.tools
import medicc.tree_hash
import medicc.nni
import medicc.fst_utils
from medicc.ancestors import reconstruct_ancestors
from medicc.core import *
from medicc.factory import *
from medicc._version import __version__

with open(os.path.join(os.path.dirname(__file__), 'logging_conf.yaml'), 'rt') as f:
    config = yaml.safe_load(f.read())
logging.config.dictConfig(config)
