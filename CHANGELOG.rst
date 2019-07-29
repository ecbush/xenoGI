==========
Change Log
==========

-------------------
2.2.0_ - 2019-07-29
-------------------

- runBlast now checks to see if a particular comparison has already been run, and if so doesn't repeat. This is mostly needed for a development version working with large numbers of species.
- Replaced the geneNames object with a genes object more suitable to large numbers of strains (thousands).
- The string form of gene names now also contains the number.
- Now use string based names for nodes only, rather than numerical and string. As a part of these changes, we've created a strainInfo file, and a corresponding parameter strainInfoFN in the example params.py.
- parseGenbank can now optionlly create dna files (depending on the new parameter dnaBasedGeneTrees).
- added makeSpeciesTree option. This makes gene trees for the hard core set of gene families (using FastTree with GTR+CAT) and then combines them into a species tree using ASTRAL.
- Updated example assemblies so they mostly are from type strains.
- Changed the parameter numThreads to be called numProcesses (since that's a more acurate description of what's being done under the hood).

-------------------
2.1.0_ - 2019-05-29
-------------------

- Fixed a bug in island formation that arose when there were no islands on a branch.
- Changed naming for blast files to avoid conflicts with strain names. (In particular with the '-' character.)
- parseGenbank now throws an error when it encounters a gbff file that lacks protein annotations.
- Created a sharedScores class that shares memory between separate processes. This reduces memory usage significantly in calcScores.

-------------------
2.0.0_ - 2019-03-20
-------------------

This is a major release, meaning things like parameter files and output files have changed.

- Reduced RAM requirements. We've eliminated the norm score. Synteny scores are now based on raw scores, and thresholds for family formation calculated directly for each strain pair.
- Improved speed. This results from changes to island formation, and makes a particularly big difference for large data sets (50 or more strains).
- Simplified the user parameters file. Parameters users are unlikely to ever change have been moved elsewhere.
- Added plotScoreHists flag. Makes rawSc.pdf, synSc.pdf, coreSynSc.pdf.
- Introduced several additional classes (LocusFamily, LocusIsland). This is largely for the sake of future development, and should make it easier to capture things like gene duplication.
- Modified Score class to better organize scores from particular strain pairs. This means the current version of xenoGI will not read scores objects from older versions.

-------------------
1.1.2_ - 2018-10-06
-------------------

Fixed a bug in printAnalysis.

-------------------
1.1.1_ - 2018-06-11
-------------------

Interactive analysis is now part of the main xenoGI script.

-------------------
1.1.0_ - 2018-06-03
-------------------

This is the first pypi release. In setting it up to be distributed as a package, we've also changed the user interface (see README.rst).

Hereafter, tagged releases correspond to pypi releases. The master branch will have the latest stable code. 

-------------------
1.0.0_ - 2017-11-14
-------------------

Initial release, corresponding to our article: "xenoGI: reconstructing the history of genomic island insertions in clades of closely related bacteria".

.. _2.2.0:  https://github.com/ecbush/xenoGI/compare/v2.1.0...v2.2.0
.. _2.1.0:  https://github.com/ecbush/xenoGI/compare/v2.0.0...v2.1.0
.. _2.0.0:  https://github.com/ecbush/xenoGI/compare/v1.1.2...v2.0.0
.. _1.1.2:  https://github.com/ecbush/xenoGI/compare/v1.1.1...v1.1.2
.. _1.1.1:  https://github.com/ecbush/xenoGI/compare/v1.1.0...v1.1.1
.. _1.1.0:  https://github.com/ecbush/xenoGI/compare/v1.0.0...v1.1.0
.. _1.0.0:  https://github.com/ecbush/xenoGI/releases/tag/v1.0.0
