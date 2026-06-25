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
import medicc.spr
import medicc.nni
import medicc.tree_hash
from medicc.ancestors import reconstruct_ancestors, reconstruct_ancestors_incremental, reconstruct_ancestors_constrained_sankoff
from medicc.core import *
from medicc.factory import *
from medicc._version import __version__


with open(os.path.join(os.path.dirname(__file__), 'logging_conf.yaml'), 'rt') as f:
    config = yaml.safe_load(f.read())
logging.config.dictConfig(config)
