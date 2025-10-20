from setuptools import setup, find_packages

CLASSIFIERS = [
    "License :: OSI Approved :: MIT License",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Operating System :: Microsoft :: Windows",
    "Operating System :: POSIX",
    "Operating System :: Unix",
    "Operating System :: MacOS",
]

# Dependency versions locked to the historically compatible window
INSTALL_REQUIRES = [
    # Core scientific + parallel
    "pathos==0.2.9",
    "cppy==1.2.0",
    "kiwisolver==1.3.1",
    "gdown==4.4.0",
    "cmake==3.18.2",
    "cython==0.29.20",
    "multiprocess==0.70.13",

    # Numpy/Scipy pinned for COBRA & QMINOSPy compatibility
    "numpy>=1.19.0,<1.22",   # <=1.21.x keeps qminospy and np.object happy
    "scipy>=1.5.0,<1.6",     # 1.5.x avoids vstack IndexError and ABI mismatch
    "pandas==1.1.5",
    "statsmodels>=0.13.0,<0.14",
    "sympy==1.12",

    # COBRA + plotting
    "cobra==0.18.1",
    "matplotlib==3.3.4",
    "seaborn==0.11.2",

    # Bio / I/O helpers
    "biopython==1.79",
    "Faker==8.5.1",
    "openpyxl==3.0.10",
    "tables==3.7.0",
    "numexpr==2.7.3",
    "tqdm==4.62.3",
]

EXTRAS_REQUIRE = {
    "interactive": ["jupyter", "ipykernel"],
}

with open("README.md", encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="human_me",
    version="0.1.0",
    author="Hratch Baghdassarian",
    author_email="hmbaghdassarian@gmail.com",
    description="Python package to generate and analyze human ME Models.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/hmbaghdassarian/human_me",
    packages=find_packages(include=["human_me*"], exclude=["*test*"]),
    include_package_data=True,
    python_requires=">=3.8,<3.10",
    install_requires=INSTALL_REQUIRES,
    extras_require=EXTRAS_REQUIRE,
    classifiers=CLASSIFIERS,
    license="MIT",
)
