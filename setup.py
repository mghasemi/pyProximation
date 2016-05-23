from distutils.core import setup
from distutils.extension import Extension

Description = """This package provides basic operations with measures,\
             orthogonal system of functions and solving integro-differential equations. """

setup(
    name = 'IntgDiff',
    version = '1.0.1',
    author = 'Mehdi Ghasemi',
    author_email = 'mehdi.ghasemi@gmail.com',
    packages = ['IntgDiff'],
    #url = 'https://github.com/mghasemi/CvxAlgGeo.git',
    license = 'MIT License',
    description = Description,
    #long_description = open('README.txt').read(),
)