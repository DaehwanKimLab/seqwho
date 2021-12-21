import sys, os
import pkg_resources

modulepath = os.path.dirname(__file__)
if modulepath not in sys.path:
    sys.path.append(os.path.dirname(__file__))

__version__ = pkg_resources.get_distribution('seqwho').version
