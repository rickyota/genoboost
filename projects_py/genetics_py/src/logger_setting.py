import logging


def setting():
    # logging.basicConfig(level=logging.INFO,
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s:%(name)s:%(funcName)s:%(lineno)d[%(levelname)s] %(message)s')
    libs = ['matplotlib']
    for lib in libs:
        logging.getLogger(lib).setLevel(logging.INFO)
