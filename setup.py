from distutils.core import setup
from distutils.extension import Extension

Description = """This package provides basic operations with measures, orthogonal system of functions, interpolation and solving integro-differential equations. """

setup(
    name='pyProximation',
    version='1.3.2',
    author='Mehdi Ghasemi',
    author_email='mehdi.ghasemi@gmail.com',
    packages=['pyProximation'],
    #url = 'https://github.com/mghasemi/pyProximation.git',
    license='MIT License',
    description=Description,
    long_description=open('README.rst').read(),
    keywords="approximation, interpolation, measures, collocation",
    install_requires=['sympy', 'numpy', 'scipy', 'matplotlib', 'itertools']
)
