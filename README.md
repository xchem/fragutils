# Frag

[![test](https://github.com/xchem/frag/actions/workflows/test.yaml/badge.svg)](https://github.com/xchem/frag/actions/workflows/test.yaml)
[![release](https://github.com/xchem/frag/actions/workflows/release.yaml/badge.svg)](https://github.com/xchem/frag/actions/workflows/release.yaml)

[![License](http://img.shields.io/badge/license-Apache%202.0-blue.svg?style=flat)](https://github.com/xchem/frag/blob/master/LICENSE.txt)

![PyPI](https://img.shields.io/pypi/v/xchem-frag)

[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

Basic RDKit based Python tools for analysis of protein-ligand interactions.

>   This was originally the fragalysis GitHub repository, which is now being used for
    another purpose. All new releases of this package will come from here,
    with the new name `xchem-frag`.

Currently contains: -

1. Clustering - based on WONKA method - but separated from that code-base.
    Cluster waters, residues, ligands and pharmacophores. (Under development)
2. Astex Fragment Network - implementation on the basis of their recent paper
3. Conformer generation code - based on known X-ray structures
4. Support for the neo4j 4.4.2 graph database

## Pre-commit
The project uses [pre-commit] to enforce linting of files prior to committing
them to the upstream repository.

To get started review the pre-commit utility and then set-up your local clone
by following the **Installation** and **Quick Start** sections of the
pre-commit documentation.

Ideally from a Python environment...

    python -m venv venv
    source venv/bin/activate

    pip install --upgrade pip
    pip install -r build-requirements.txt
    pre-commit install -t commit-msg -t pre-commit

Now the project's rules will run on every commit and you can check the
state of the repository as it stands with...

    pre-commit run --all-files

## Publishing (to PyPI)
We rely on out **release** GitHub workflow to publish to PyPI, something that
is done automatically when the repository main branch is tagged.

---

[pre-commit]: https://pre-commit.com
