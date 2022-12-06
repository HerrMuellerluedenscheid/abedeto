import shutil
import os
import sys
import logging

logger = logging.getLogger("util")


def fail(message):
    logger.critical(message)
    sys.exit(1)


def create_directory(directory, force):
    if not os.path.exists(directory):
        os.makedirs(directory)
        logger.info("creating new directory: %s" % directory)
    else:
        if force:
            logger.info("overwriting and creating new directory: %s" % directory)
            shutil.rmtree(directory)
            create_directory(directory, False)
        else:
            logger.info("directory exists: %s" % directory)
