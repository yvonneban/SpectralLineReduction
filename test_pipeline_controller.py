import sys
from pipeline_controller import *

def main(argv):
    print(argv)
    H = HandlePipelineOptions()
    H.parse_options(argv, 'test', print_options=True)


main(sys.argv[1:])
