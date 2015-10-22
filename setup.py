

import os
#from setuptools import setup
from distutils.core import setup

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "abedeto",
    version = "0.0.1",
    author = "Marius Kriegerowski",
    author_email = "marius@gfz-potsdam.de",
    description = ("Array Beam Depth Tool."),
    license = "BSD",
    keywords = "Array Beam Depth",
    packages=['abedeto'],
              #'abedeto.arclink'],
    package_dir={'abedeto': 'src'},
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 1"
        "Topic :: Utilities",
        "License :: OSI Approved :: BSD License",
    ],
    scripts=['src/abedeto'],
)


