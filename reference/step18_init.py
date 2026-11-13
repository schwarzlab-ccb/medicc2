# REFERENCE for migration Step 18 (not part of REF; untracked). Goes to medicc/__init__.py.
# One added line: import medicc.fstengine (REF also adds import medicc.nni: that is Phase 7).

import logging.config
import os

import yaml

import medicc.bootstrap
import medicc.core
import medicc.fstengine
import medicc.io
import medicc.nj
import medicc.plot
import medicc.sim
import medicc.stats
import medicc.tools
from medicc.ancestors import reconstruct_ancestors
from medicc.core import *
from medicc.factory import *
from medicc._version import __version__

with open(os.path.join(os.path.dirname(__file__), 'logging_conf.yaml'), 'rt') as f:
    config = yaml.safe_load(f.read())
logging.config.dictConfig(config)
