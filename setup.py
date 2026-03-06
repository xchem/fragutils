from setuptools import setup, find_packages

# To use a consistent encoding
from codecs import open
from os import environ, path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, "README.rst"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="xchem-frag",
    version=environ.get("FRAG_VERSION", "0.0.0"),
    description="Library for fragment based analysis",
    long_description=long_description,
    # The project's main homepage.
    url="https://github.com/xchem/frag.git",
    # Author details
    author="Max Winokan",
    author_email="max.winokan@diamond.ac.uk",
    # Choose your license
    license="Apache 2.0",
    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Topic :: Software Development :: Build Tools",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    # What does your project relate to?
    keywords="",
    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(exclude=["contrib", "docs", "tests"]),
    # Alternatively, if you want to distribute just a my_module.py, uncomment
    # this:
    #   py_modules=["my_module"],
    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=[
        "neo4j-driver==4.4.11",
        "ipython>5.4.1",
        "tqdm>=4.65.0",
        "numpy>=1.25.2",
        "requests>=2.31.0",
    ],
    # List additional groups of dependencies here (e.g. development
    # dependencies). You can install these using the following syntax,
    # for example:
    # $ pip install -e .[dev,test]
    extras_require={"dev": ["check-manifest"], "test": ["coverage"]},
    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    package_data={"frag": ["utils/data/RDKitPh4.fdef"]},
    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    #    data_files=[('my_data', ['data/data_file'])],
    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
        "console_scripts": [
            "build_db=frag.network.scripts.build_db:main",
            "header_shred=frag.network.scripts.header_shred:main",
            "colate=frag.network.scripts.colate_all:main",
        ]
    },
)
