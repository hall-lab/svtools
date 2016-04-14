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

    install_requires=['pysam>=0.8.1', 'numpy', 'scipy', 'statsmodels', 'pandas', 'setuptools'],
    scripts=['scripts/create_coordinates'],
    
    entry_points={
        'console_scripts': [
            'svtools=svtools.cli:main',
            ]
    },
)
