from setuptools import setup
from Cython.Build import cythonize

setup(
    name='DiffusionC function',
    ext_modules=cythonize("diffusionC.pyx"),
    zip_safe=False,
)