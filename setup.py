from setuptools import setup, find_packages

__version__ = '0.0.1-dev'

setup(
    name = 'subsamplebam',
    version = __version__,
    py_modules = ['subsamplebam', 'subsample_mindepth'],
    setup_requires = [
        'nose',
        'python-coveralls'
    ],
    entry_points = {
        'console_scripts': [
            'subsamplebam = subsamplebam:main' 
        ]
    },
    author = 'Tyghe Vallard',
    author_email = 'vallardt@gmail.com',
    description = 'Subsample samtools view randomly',
    license = 'GPL v2',
    keywords = 'bam, samtools, subsample',
    url = 'https://github.com/necrolyte2/subsamplebam',
)
