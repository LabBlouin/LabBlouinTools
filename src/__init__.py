# below is boilerplate code in order to recognize all files in this directory

import os

thisdir = os.path.split(os.path.realpath(__file__))[0]
itlist = os.listdir(thisdir)
__all__ = [os.path.split(x)[-1].strip('.py') for x in itlist if x.endswith('.py') and not x.endswith('__init__.py')]
