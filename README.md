# catmaid-interface
Python tools for interacting with CATMAID projects

This is a slowly increasing toolbox for extracting CATMAID (http://www.catmaid.org) annotations, typically EM reconstructions of neurons and eventually visualizing and doing various kinds of analysis of neuronal morphology, connectivity, and connections between the two. See Quantitative neuroanatomy for connectomics in Drosophila (http://biorxiv.org/content/early/2016/02/17/026617) for more details of what sort of things will eventually be possible.

Currently it needs the following packages:
* json - for parsing and handling the JSON files that CATMAID trades in.
* requests - for happily interacting with a CATMAID server.
* scipy - For quickly handling graphs using sparse matrices.

A similar tool has been developed by Greg Jefferis for R, which can be found at https://github.com/jefferis/rcatmaid
