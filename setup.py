from distutils.core import setup
from distutils.extension import Extension

Description = """This package provides basic operations with measures, orthogonal system of functions, interpolation and solving integro-differential equations. """

setup(
    name = 'ApproxPy',
    version = '1.0.0',
    author = 'Mehdi Ghasemi',
    author_email = 'mehdi.ghasemi@gmail.com',
    packages = ['ApproxPy'],
    #url = 'https://github.com/mghasemi/ApproxPy.git',
    license = 'MIT License',
    description = Description,
    long_description = open('README.rst').read(),
    keywords = "approximation, interpolation, measures, collocation",
    install_requires = ['sympy', 'numpy', 'scipy', 'matplotlib', 'itertools']
)