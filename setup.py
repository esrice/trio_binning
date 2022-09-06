from setuptools import Extension, setup

setup(ext_modules=[Extension(name="trio_binning.kmers_c", sources=["c/kmers.c"])])
