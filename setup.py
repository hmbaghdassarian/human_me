# python setup.py develop
# python setup.py install
from setuptools import setup
from setuptools import find_packages


CLASSIFIERS = '''\
License :: OSI Approved :: MIT license
Programming Language :: Python :: 3.6 :: 3.9
Topic :: Genome-Scale Modeling
Operating System :: Microsoft :: Windows
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
'''

DISTNAME = 'human_me'
AUTHOR = 'Hratch Baghdassarian'
AUTHOR_EMAIL = 'hmbaghdassarian@eng.ucsd.edu'
DESCRIPTION = 'Python package to generate and analyze human ME Models.'
LICENSE = 'MIT'

VERSION = '0.1.0'
ISRELEASED = False

PYTHON_MIN_VERSION = '3.8'
PYTHON_MAX_VERSION = '3.9'
PYTHON_REQUIRES = f'>={PYTHON_MIN_VERSION}, <={PYTHON_MAX_VERSION}'

INSTALL_REQUIRES = [
    'pathos==0.2.9',
    'cppy==1.2.0',
    'kiwisolver==1.3.1',
    'gdown==4.4.0',
    'cmake==3.18.2',
    'cython==0.29.20',
    'multiprocess==0.70.13',
    'numpy==1.19.5',
    'pandas==1.1.5',
    'scipy==1.5.4',
    'statsmodels==0.10.2',
    'sympy==1.12',
    'tqdm==4.62.3',
    'cobra==0.18.1',
    'matplotlib==3.3.4',
    'seaborn==0.11.2',
    'biopython==1.79',
    'Faker==8.5.1',
    'openpyxl==3.0.10',
    'tables==3.6.1',
    'numexpr==2.7.3'
    # 'swiglpk==5.0.5'
    # 'ruamel_yaml==0.17.4'
]

EXTRAS_REQUIRES = {'interactive': ['jupyter', 'ipykernel']
                  }

PACKAGES = [
    'human_me'
]

with open('README.md') as f:
    long_description = f.read()

metadata = dict(
    name=DISTNAME,
    version=VERSION,
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=long_description,
    url='https://github.com/hmbaghdassarian/human_me',  # homepage
    packages=find_packages(include=('human_me*'), exclude=('*test*',)),  # PACKAGES
#     scripts=['install_solver.py'],
    include_package_data=True,
    project_urls={'Documentation': 'https://hmbaghdassarian.github.io/human_me/'},
    # py_modules=['io'],
    python_requires=PYTHON_REQUIRES,
    install_requires=INSTALL_REQUIRES,
    extras_require=EXTRAS_REQUIRES,
    classifiers=[CLASSIFIERS],
    license=LICENSE
)


def setup_package() -> None:
    setup(**metadata)


if __name__ == '__main__':
    setup_package()
