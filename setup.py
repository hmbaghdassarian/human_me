# python setup.py develop
# python setup.py install
from setuptools import setup


CLASSIFIERS = '''\
License :: OSI Approved
Programming Language :: Python :: 3.6 :: 3.9
Topic :: Genome-Scale Modeling
Operating System :: Microsoft :: Windows
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
'''

DISTNAME = 'human_me'
AUTHOR = 'Hratch Baghdassarian'
AUTHOR_EMAIL = 'hmbaghdassarian@gmail.com'
DESCRIPTION = 'Python package to generate and analyze human ME Models.'
LICENSE = 'MIT'
README = 'Python package to generate and analyze human ME Models.'

VERSION = '0.1.0'
ISRELEASED = False

PYTHON_MIN_VERSION = '3.6'
PYTHON_MAX_VERSION = '3.9'
PYTHON_REQUIRES = f'>={PYTHON_MIN_VERSION}, <={PYTHON_MAX_VERSION}'

INSTALL_REQUIRES = [
    'pathos==0.2.7'
    'cmake==3.18.2'
    'cython==0.29.20'
    'more-itertools==8.4.0'
    'numpy==1.19.5'
    'pandas==1.1.5'
    'scipy==1.5.4'
    'six==1.16.0'
    'statsmodels==0.10.2'
    'sympy==1.9'
    'tqdm==4.62.3'
    'cobra==0.18.1'
    'xlrd==1.2.0'
    'requests==2.26.0'
    'tables==3.6.1'
    'matplotlib==3.3.4'
    'seaborn==0.11.2'
    'biopython==1.79'
    'Faker==8.5.1'
    'googledrivedownloader==0.4'
    'make==4.1 '
]

PACKAGES = [
    'human_me'
]

metadata = dict(
    name=DISTNAME,
    version=VERSION,
    long_description=README,
    packages=PACKAGES,
    python_requires=PYTHON_REQUIRES,
    install_requires=INSTALL_REQUIRES,
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    description=DESCRIPTION,
    classifiers=[CLASSIFIERS],
    license=LICENSE
)


def setup_package() -> None:
    setup(**metadata)


if __name__ == '__main__':
    setup_package()
