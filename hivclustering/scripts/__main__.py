from .hivnetworkannotate import main as hivnetworkannotate_main, parse_args as hivnetworkannotate_parse_args
from .hivnetworkcsv import make_hiv_network
from .TNS import main as tns_main


def hivnetworkcsv():
    # not sure why this would be necessary...
    __spec__ = None
    make_hiv_network()


def hivnetworkannotate():
    args = hivnetworkannotate_parse_args()
    hivnetworkannotate_main(args)


def TNS():
    tns_main()
