from setuptools import setup, find_packages
import versioneer

'''
This package provides tools for combining, genotyping and refining structural variant calls from LUMPY in a highly scalable way. The current package efficiently scales to process thousands of individuals.
'''


setup(
    name='svtools',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),

    description='Tools for processing and analyzing structural variants',
    long_description=__doc__,

    url='https://github.com/hall-lab/svtools',
    author='Ira Hall lab',
    author_email='dlarson@genome.wustl.edu',
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
    ],

    keywords='genomics structural variants sv bioinformatics',

    packages=find_packages(exclude=['tests']),
    include_package_data=True,

    install_requires=[
        'svtyper==0.7.1',
        'numpy==1.16.6',
        'scipy==1.2.3',
        'statsmodels==0.10.0',
        'pandas==0.21.1',
        'setuptools==44.1.1',
        'google-auth==1.35.0',
        'google-cloud-storage==1.44.0',
        'google-compute-engine==2.8.13',
        'crcmod==1.7',
        'logzero==1.7.0'
    ],
    scripts=['scripts/create_coordinates'],

    entry_points={
        'console_scripts': [
            'svtools=svtools.cli:main',
            ]
    },
)
