.. pyranges documentation master file, created by
   sphinx-quickstart on Fri Jun 23 11:27:13 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

The pyranges0 documentation
===========================
PyRanges is a Python library specifically designed for efficient and intuitive manipulation of genomics data,
particularly genomic intervals (like genes, genomic features, or reads).
The library is optimized for fast querying and manipulation of genomic annotations.

PyRanges is open source, and hosted at GitHub: https://github.com/pyranges/pyranges0


Deprecation note
~~~~~~~~~~~~~~~~

This repo currently hosts the "version 0" of pyranges. A redesigned, faster version 1 is available at https://github.com/pyranges/pyranges1
Version 1 will soon replace version 0. We encourage existing v0 users to migrate to v1 (see guide at https://pyranges1.readthedocs.io/en/latest/migration_guide.html), and new users to directly adopt v1. Read about v1 in our preprint at https://doi.org/10.64898/2025.12.11.693639


Citation
~~~~~~~~

Stovner EB, Sætrom P (2020) PyRanges: efficient comparison of genomic intervals in Python. *Bioinformatics 36(3):918-919*  http://dx.doi.org/10.1093/bioinformatics/btz615


Documentation outline
~~~~~~~~~~~~~~~~~~~~~

#. 🚀 :doc:`Installation instructions <./installation>`
#. 🚀 :doc:`The tutorial <./tutorial>`,  recommended for all new users
#. 🚀 :doc:`The how-to pages <./how_to_pages>`, further below, where functionalities are grouped by topic
#. 🚀 The `API reference <./autoapi/index.html>`_, where all methods are explained in detail.
#. 🚀 :doc:`The developer guide <./developer_guide>`, to follow in order to contribute to PyRanges.

Asking for help
~~~~~~~~~~~~~~~

If you encounter bugs, or the documentation is not enough and you cannot accomplish a specific task of interest, or if you'd like new features implemented, open an issue on GitHub: https://github.com/pyranges/pyranges0/issues


.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Contents:

   installation

   tutorial

   how_to_pages

   developer_guide

       

.. * :ref:`search`
.. * :ref:`modindex`
