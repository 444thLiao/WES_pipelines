import sys
from os.path import dirname
sys.path.insert(0,dirname(dirname(__file__)))

from luigi_pipelines.share_luigi_tasks import PrintReads, Annovar1, Annovar2
