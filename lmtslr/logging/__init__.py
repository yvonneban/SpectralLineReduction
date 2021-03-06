"""
Logging Framework for lmtslr
TODO: Simple test cases
"""

import logging as pylogging
import logging.handlers as handlers
import os

DEBUG0 = 5
DEBUG1 = 6
DEBUG2 = 7
DEBUG3 = 8
DEBUG4 = 9
DEBUG  = 10
INFO   = 20
WARNING = 30
ERROR  = 40
CRITICAL = 50

def _get_lmtslr_logname():
    """
    Returns the path to the lmtslr log file
    """
    if 'HOME' in os.environ:
        home = os.environ['HOME']
        if os.path.exists(os.path.join(home, '.lmtslr')):
            if os.path.isdir(os.path.join(home, '.lmtslr')):
                fname = os.path.join(home, '.lmtslr', 'lmtslr.log')
            else:
                fname = os.path.join(home, 'lmtslr.log')
        else:
            os.makedirs(os.path.join(home, '.lmtslr'))
            fname = os.path.join(home, '.lmtslr', 'lmtslr.log')
        return fname
    else:
        return 'lmtslr.log'


class NullHandler(pylogging.Handler):
    """A simple NullHandler to be used by modules used
    inside lmtslr 
    Every module can use at least this one null handler.
    It does not do anything. The application that uses the
    lmtslr library can set up logging, and setup appropriate
    handlers, and then the lmtslr library will appropriately
    handle correct logging.
    A lmtslr module should do the following:
    from lmtslr.logging import logger
    logger.name = __name__  #to appropriately catch that module's logs
    And then you can use the logger object within the modules
    """
    def emit(self, record):
        pass

def will_debug():
    return logger.isEnabledFor(pylogging.DEBUG)

def add_file_handler(log, fname):
    handler = handlers.RotatingFileHandler(fname, maxBytes=1e6,
                                           backupCount=5)
    # create formatter
    formatter = pylogging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")

    handler.setLevel(pylogging.DEBUG)
    handler.setFormatter(formatter)
    #print "Setting up Rotating File Handler %s" % fname
    log.addHandler(handler)

def add_stream_handler(log):
    handler = pylogging.StreamHandler()
    #level = lmtslrParams['lmtslr'].get('loglevel', 10)
    handler.setLevel(pylogging.DEBUG)
    #handler.setLevel(level)
    log.addHandler(handler)

def debug_factory(logger, debug_level):
    def custom_debug(msg, *args, **kwargs):
        if logger.level >= debug_level:
           return
        logger._log(debug_level, msg, args, kwargs)
    return custom_debug    


logger = pylogging.getLogger('lmtslr')
logger.logging = pylogging
logger.addHandler(NullHandler())
logger.will_debug = will_debug
add_file_handler(logger, _get_lmtslr_logname())
add_stream_handler(logger)
logger.setLevel(DEBUG0)

for i in range(5, 10):
    pylogging.addLevelName(i, 'DEBUG%i' % (i-5))
    setattr(logger, 'debug%i' % (i-5), debug_factory(logger, i))
