from setuptools import setup, Extension
import numpy as np

try:
    from Cython.Build import cythonize
    ext_modules = cythonize(
        [
            Extension("pyatac.fragments", ["pyatac/fragments.pyx"], include_dirs=[np.get_include()]),
            Extension("nucleoatac.multinomial_cov", ["nucleoatac/multinomial_cov.pyx"], include_dirs=[np.get_include()]),
        ],
        compiler_directives={"language_level": "3"},
    )
except ImportError:
    ext_modules = []

setup(
    name='NucleoATAC',
    version='0.5.0',
    description='Python package for calling nucleosomes with ATAC-Seq',
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    keywords='ATAC-Seq sequencing bioinformatics',
    url='https://github.com/sjessa/NucleoATAC',
    author='Alicia Schep',
    author_email='aschep@stanford.edu',
    license='MIT',
    python_requires='>=3.9',
    packages=['pyatac', 'pyatac.pwm', 'nucleoatac', 'nucleoatac.vplot'],
    install_requires=[
        'cython >= 3.0',
        'numpy >= 1.22',
        'scipy >= 1.7',
        'pysam >= 0.20',
        'matplotlib >= 3.5',
    ],
    ext_modules=ext_modules,
    scripts=['bin/pyatac', 'bin/nucleoatac'],
    include_package_data=True,
    zip_safe=False,
)
