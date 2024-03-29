|GitHub release| |made-with-python|

.. |GitHub release| image:: https://img.shields.io/github/v/release/fabothch/TMprep
   :target: https://github.com/fabothch/TMprep/releases/latest

.. |made-with-python| image:: https://img.shields.io/badge/Made%20with-Python-1f425f.svg
   :target: https://www.python.org/

==============
TMprep Project
==============

Introduction
============

``tmprep`` is designed to create a usable control file for TURBOMOLE > version 7.5.
Useful settings are predefined. The files .CHRG and .UHF containing integer numbers
for charge and unpaired number of electrons are read. Symmetry information in Schoenflies 
notation can be read from the .SYM file. Settings as defined in the ~/.cefinerc are read 
and are compatible to ``tmprep``. The number of electrons which 
are printed out are only calculated within ``tmprep`` and do not stem from an EHT guess.

**WARNING**

    This work is in development and should not be used in production runs!
    (The warning mainly focuses on the creation of the basis and auxbasis files!)

Setup:
======

Download tmprep.py and make it executable:

.. code:: bash

    $ chmod u+x tmprep.py


Usage example:
==============

.. code:: bash

    $ tmprep.py -func r2scan-3c -scfconv 7 -radsize 10 -cosmo 80.0 -sym c1

The basis and auxbasis files can be generated by ``tmprep`` and this can be prone 
to errors so use with caution.

.. code:: bash

    $ tmpprep.py -gen_bas -gen_auxbas <other settings>


Otherwise ridft will generate the basis and auxbasis files.

Current options:
================

.. code:: bash

   tmprep.py - Command line input for TURBOMOLE V 0.1.4 FB 2021
   usage: tmprep.py [-h] [-chrg] [-uhf] (-func  | -wfunc  | -hf) [-basis] [-sym]
                    [-radsize] [-scfconv] [-grid] [-d3zero] [-d3] [-d3atm] [-d4]
                    [-donl] [-novdw] [-nori] [-cosmo] [-dcosmors] [-noopt]
                    [-gen_bas] [-gen_auxbas]

   optional arguments:
     -h, --help            show this help message and exit
     -func , --func        Density functional aproximation.
     -wfunc , --wfunc      Wavefunction based method. E.g. HF.
     -hf, --hf             Selecting Hartree Fock (HF).

   Options:
     -chrg , --chrg        Charge of the molecule.
     -uhf , --uhf          Integer number of unpaired electrons of the molecule.
     -basis , --basis      Basis set
     -sym , --sym          Symmetry in Schönflies lower case
     -radsize , --radsize 
                           Radsize
     -scfconv , --scfconv 
                           scfconv
     -grid , --grid        DFA grid
     -d3zero, --d3zero     D3(0)
     -d3, --d3             D3(BJ)
     -d3atm, --d3atm       D3(BJ)ATM
     -d4, --d4             D4
     -donl, --donl         NL dispersion correction, needs sym c1
     -novdw, --novdw       No dispersion correction
     -nori, --nori         No resolution of the identity
     -cosmo , --cosmo      Dielectric constant for COSMO
     -dcosmors , --dcosmors 
                           Add DCOSMO-RS to control. Usage: -dcosmors [solvent].
                           Options are acetone, chcl3, acetonitrile, ch2cl2,
                           dmso, h2o, methanol, thf, toluene, octanol, woctanol,
                           hexadecane. It is not necessary to additionally use
                           -cosmo
     -noopt, --noopt       Use cartesian coordinates.
     -gen_bas, --gen_bas   Create basis file.
     -gen_auxbas, --gen_auxbas
                           Create auxbasis file.

   tmprep is designed to create a usable control file for TURBOMOLE > version 7.5.
   Useful settings are predefined. .CHRG and .UHF files containing integer numbers
   for charge and unpaired number of electrons are read. Settings as defined in the 
   ~/.cefinerc are read and are compatible to tmprep. The number of electrons which 
   are printed out are only calculated within tmprep and do not stem from an EHT guess.

   Usage exmple:

   tmprep.py -func r2scan-3c -scfconv 7 -radsize 10 -cosmo 80.0 -sym c1







License
=======

``tmprep`` is free software: you can redistribute it and/or modify it under the terms
of the GNU Lesser General Public License as published by the Free Software 
Foundation, either version 3 of the License, or (at your option) any later version.

``tmprep`` is distributed in the hope that it will be useful, but without any 
warranty; without even the implied warranty of merchantability or fitness for 
a particular purpose. See the GNU Lesser General Public License for more details.

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in ``tmprep`` by you, as defined in the GNU Lesser General Public license, 
shall be licensed as above, without any additional terms or conditions.
