#!/usr/bin/env python3

# Copyright (C) 2021 Fabian Bohle
#
# tmprep is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# tmprep is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with tmprep. If not, see <https://www.gnu.org/licenses/>.

"""tmprep is designed to create a usable control file for TURBOMOLE > version 7.5 .
Useful settings are predefined. .CHRG and .UHF files containing integer numbers
for charge and unpaired number of electrons are read. Settings as defined in the
~/.cefinerc are read and are compatible to tmprep. The number of electrons which
are printed out are only calculated within tmprep and do not stem from an EHT guess.

Usage exmple:

tmprep.py -func r2scan-3c -scfconv 7 -radsize 10 -cosmo 80.0 -sym c1

"""

import os
import sys
import argparse
from collections import Counter

version = "0.1.5"

class Settings():
    """Settings for the preparation of a TM control file."""
    
    solvent_dcosmors = {
        "acetone": [20.7, "$dcosmo_rs file=propanone_25.pot"],
        "chcl3": [4.8, "$dcosmo_rs file=chcl3_25.pot"],
        "acetonitrile": [36.6, "$dcosmo_rs file=acetonitrile_25.pot"],
        "ch2cl2": [9.1, "$dcosmo_rs file=chcl3_25.pot"],
        "dmso": [47.2, "$dcosmo_rs file=dimethylsulfoxide_25.pot"],
        "h2o": [80.1, "$dcosmo_rs file=h2o_25.pot"],
        "methanol": [32.7, "$dcosmo_rs file=methanol_25.pot"],
        "thf": [7.6, "$dcosmo_rs file=thf_25.pot"],
        "toluene": [2.4, "$dcosmo_rs file=toluene_25.pot"],
        "octanol": [9.86, "$dcosmo_rs file=octanol_25.pot"],
        "woctanol": [8.1, "$dcosmo_rs file=wet-octanol_25.pot"],
        "hexadecane": [2.08, "$dcosmo_rs file=hexadecane_25.pot"],
    }

    def read_remote_configuration(self):
        """read ~/.censorc and update to user set defaults"""
        rc_defaults = {}
        path_cefinerc =os.path.expanduser('~/.cefinerc')
        if os.path.isfile(path_cefinerc):
            # read defaults concerning memory to internal_defaults dictionary
            with open(path_cefinerc, 'r') as inp:
                data = inp.readlines()
                for line in data:
                    if 'ricore' in line:
                        try:
                            rc_defaults['ricore'] = int(line.split()[1])
                        except (ValueError, TypeError):
                            pass
                    if 'maxcor' in line:
                        try:
                            rc_defaults['maxcor'] = int(line.split()[1])
                        except (ValueError, TypeError):
                            pass
                    if 'rpacor' in line:
                        try:
                            rc_defaults['rpacor'] = int(line.split()[1])
                        except (ValueError, TypeError):
                            pass
                    if 'func' in line:
                        try:
                            rc_defaults['functional'] = line.split()[1].strip()
                        except (ValueError, TypeError):
                            pass
                    if 'bas' in line:
                        try:
                            rc_defaults['basis'] = line.split()[1].strip()
                        except (ValueError, TypeError):
                            pass
                    if 'grid' in line:
                        try:
                            rc_defaults['grid'] = line.split()[1].strip()
                        except (ValueError, TypeError):
                            pass
                    if 'scfconv' in line:
                        try:
                            rc_defaults['scfconv'] = int(line.split()[1])
                        except (ValueError, TypeError):
                            pass
        self.rc_defaults = rc_defaults

    def __init__(self, **kwargs):
        """Constructor"""
        ### read data from .cefinerc
        self.read_remote_configuration()
        # adjust other settings
        self.charge = kwargs.get('charge', 0)
        self.unpaired = kwargs.get('unpaired', 0)
        self.functional = kwargs.get('functional', self.rc_defaults.get('functional', 'pbeh-3c'))
        self.wave_func = kwargs.get('wave_func', None)
        self.basis = kwargs.get('basis', self.rc_defaults.get('basis', 'def2-mSVP'))
        self.modbasis = kwargs.get('modbasis', False)
        self.gcp = kwargs.get('gcp', "")
        self.symmetry = kwargs.get('symmetry', None)
        self.radsize = kwargs.get('radsize', None)
        self.modradsize = kwargs.get('modradsize', False)
        self.scfconv = kwargs.get('scfconv', self.rc_defaults.get('scfconv', 7))
        self.scfiterlimit = kwargs.get('scfiterlimit', 125)
        self.grid = kwargs.get('grid', self.rc_defaults.get('grid', 'm4'))
        self.modgrid = kwargs.get('modgrid', False)
        self.disp = kwargs.get('disp', 'disp3 -bj') ###
        self.moddisp = kwargs.get('moddisp', False)
        self.novdw = kwargs.get('novdw', False)
        self.ri  = kwargs.get('ri', True)
        self.ricore = kwargs.get('ricore', self.rc_defaults.get('ricore' ,4000))
        self.maxcor = kwargs.get('maxcor', self.rc_defaults.get('maxcor' ,4000))
        self.rpacor = kwargs.get('rpacor', self.rc_defaults.get('rpacor' ,1000))
        self.twoint = kwargs.get('twoint', 1000)
        self.thime = kwargs.get('thime', 4)
        self.thize = kwargs.get('thize', 0.0000001)
        self.extol = kwargs.get('extol', None)
        self.cosmo = kwargs.get('cosmo', None)
        self.dcosmors = kwargs.get('dcosmors', None)
        self.gen_bas = kwargs.get('gen_bas', False)
        self.gen_auxbas = kwargs.get('gen_auxbas', False)
        self.noauxg = kwargs.get('noauxg', False)
        self.noopt = kwargs.get('noopt', False)
        self.xcfun = kwargs.get('xcfun', False)
    

def read_chrg(default=0):
    "read molecular charge from file: .CHRG"
    chrg_path = os.path.join(os.getcwd(), '.CHRG')
    if os.path.isfile(chrg_path):
        with open(chrg_path, 'r') as inp:
            try:
                charge = int(inp.readline().strip().split()[0])
            except (FileNotFoundError, ValueError):
                print("Can't read .CHRG file! Going to exit!")
                sys.exit(1)
    else:
        charge=default
    return charge

def read_uhf(default=0):
    "read number of unpaired electrons from file: .UHF"
    uhf_path = os.path.join(os.getcwd(), '.UHF')
    if os.path.isfile(uhf_path):
        with open(uhf_path, 'r') as inp:
            try:
                unpaired = int(inp.readline().strip().split()[0])
            except (FileNotFoundError, ValueError):
                print("Can't read .UHF file! Going to exit!")
                sys.exit(1)
    else:
        unpaired = default
    return unpaired

def read_sym(default=None):
    "read molecular symmetry in schoenfliess notation from file: .SYM"
    sym_path = os.path.join(os.getcwd(), '.SYM')
    if os.path.isfile(sym_path):
        with open(sym_path, 'r') as inp:
            try:
                symmetry = str(inp.readline().strip().split()[0])
            except (FileNotFoundError, ValueError):
                print("Can't read .SYM file! Going to exit!")
                sys.exit(1)
    else:
        symmetry = default
    return symmetry

def cml(settings, argv=None):
    """
    Process commandline arguments
    """

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""tmprep is designed to create a usable control file for TURBOMOLE > version 7.5.
Useful settings are predefined. .CHRG and .UHF files containing integer numbers
for charge and unpaired number of electrons are read. Settings as defined in the 
~/.cefinerc are read and are compatible to tmprep. The number of electrons which 
are printed out are only calculated within tmprep and do not stem from an EHT guess.

Usage exmple:

tmprep.py -func r2scan-3c -scfconv 7 -radsize 10 -cosmo 80.0 -sym c1

"""
    )
    group0 = parser.add_mutually_exclusive_group(required=True)
    group1 = parser.add_argument_group("Options")
    group1.add_argument(
        "-chrg",
        "--chrg",
        dest="charge",
        action="store",
        required=False,
        type=int,
        metavar="",
        help="Charge of the molecule.",
    )
    group1.add_argument(
        "-uhf",
        "--uhf",
        dest="unpaired",
        action="store",
        required=False,
        type=int,
        metavar="",
        help="Integer number of unpaired electrons of the molecule.",
    )
    group0.add_argument(
        "-func",
        "--func",
        dest="functional",
        action="store",
        required=False,
        metavar="",
        help="Density functional aproximation.",
    )
    group0.add_argument(
        "-wfunc",
        "--wfunc",
        dest="wave_func",
        action="store",
        required=False,
        metavar="",
        help="Wavefunction based method. E.g. HF.",
    )
    group0.add_argument(
        "-hf",
        "--hf",
        dest="hf",
        action="store_true",
        required=False,
        help="Selecting Hartree Fock (HF).",
    )
    group1.add_argument(
        "-basis",
        "--basis",
        dest="basis",
        action="store",
        required=False,
        metavar="",
        help="Basis set",
    )
    group1.add_argument(
        "-sym",
        "--sym",
        dest="symmetry",
        action="store",
        required=False,
        metavar="",
        help="Symmetry in Schönflies lower case",
    )
    group1.add_argument(
        "-radsize",
        "--radsize",
        dest="radsize",
        action="store",
        required=False,
        metavar="",
        type=int,
        help=("Radsize"),
    )
    group1.add_argument(
        "-scfconv",
        "--scfconv",
        dest="scfconv",
        action="store",
        required=False,
        metavar="",
        type=int,
        help=("scfconv"),
    )
    group1.add_argument(
        "-grid",
        "--grid",
        dest="grid",
        action="store",
        required=False,
        metavar="",
        help=("DFA grid"),
    )
    group1.add_argument(
        "-d3zero",
        "--d3zero",
        dest="d3zero",
        action="store_true",
        required=False,
        help=('D3(0)'),
    )
    group1.add_argument(
        "-d3",
        "--d3",
        dest="d3",
        action="store_true",
        required=False,
        help=("D3(BJ)"),
    )
    group1.add_argument(
        "-d3atm",
        "--d3atm",
        dest="d3atm",
        action="store_true",
        required=False,
        help=("D3(BJ)ATM"),
    )
    group1.add_argument(
        "-d4",
        "--d4",
        dest="d4",
        action="store_true",
        required=False,
        help=("D4"),
    )
    group1.add_argument(
        "-donl",
        "--donl",
        dest="donl",
        action="store_true",
        required=False,
        help=("NL dispersion correction, needs sym c1"),
    )
    group1.add_argument(
        "-novdw",
        "--novdw",
        dest="novdw",
        action="store_true",
        required=False,
        help=("No dispersion correction"),
    )
    group1.add_argument(
        "-nori",
        "--nori",
        dest="nori",
        action="store_true",
        required=False,
        help=("No resolution of the identity"),
    )
    group1.add_argument(
        "-cosmo",
        "--cosmo",
        dest="cosmo",
        action="store",
        required=False,
        type=float,
        metavar="",
        help=("Dielectric constant for COSMO"),
    )
    group1.add_argument(
        "-dcosmors",
        "--dcosmors",
        dest="dcosmors",
        action="store",
        required=False,
        choices=list(settings.solvent_dcosmors.keys()),
        metavar="",
        help=("Add DCOSMO-RS to control. Usage: -dcosmors [solvent]. "
            "Options are {}. It is not necessary to additionally use -cosmo".format(', '.join(list(settings.solvent_dcosmors.keys())))),
    )
    group1.add_argument(
        "-noopt",
        "--noopt",
        dest="noopt",
        action="store_true",
        default=False,
        required=False,
        help=("Optimization in Cartesian coordinates, instead of internal redundant coordinates."),
    )
    group1.add_argument(
        "-gen_bas",
        "--gen_bas",
        dest="gen_bas",
        action="store_true",
        default=False,
        required=False,
        help=("Create basis file."),
    )
    group1.add_argument(
        "-gen_auxbas",
        "--gen_auxbas",
        dest="gen_auxbas",
        action="store_true",
        default=False,
        required=False,
        help=("Create auxbasis file."),
    )
    group1.add_argument(
        "-noauxg",
        "--noauxg",
        dest="noauxg",
        action="store_true",
        default=False,
        required=False,
        help=("Remove g-functions from auxiliary basis set."),
    )
    args = parser.parse_args(argv)
    cml_input = {}
    if args.hf:
        args.wave_func ="hf"
    if args.wave_func:
        cml_input['wave_func'] = getattr(args, "wave_func").lower()
    if args.functional:
        cml_input['functional'] = getattr(args, "functional").lower()
    if args.basis:
        cml_input['basis'] = args.basis
        cml_input['modbasis'] = True
    if args.grid:
        cml_input['grid'] = args.grid
        cml_input['modgrid'] = True
    if not args.charge:
        cml_input['charge'] = read_chrg()
    else:
        cml_input['charge'] = args.charge
    if not args.unpaired:
        cml_input['unpaired'] = read_uhf()
    else:
        cml_input['unpaired'] = args.unpaired
    if not args.symmetry:
        cml_input['symmetry'] = read_sym()
    else:
        cml_input['symmetry'] = args.symmetry
    # dispersion corrections:
    if (args.d4 or
        args.d3zero or
        args.d3 or
        args.d3atm or
        args.donl or
        args.novdw
        ):
        cml_input['moddisp'] = True
    if args.d3zero:
        cml_input['disp'] = 'disp3'
    if args.d3:
        cml_input['disp'] = 'disp3 -bj'
    if args.d3atm:
        cml_input['disp'] = 'disp3 -bj -abc'
    if args.d4:
        cml_input['disp'] = 'disp4'
    if args.donl:
        cml_input['disp'] = 'donl '
    if args.novdw:
       cml_input['disp'] = ''
    if args.nori:
        cml_input['ri'] = False
    if args.radsize is not None:
        cml_input['radsize'] = args.radsize
        cml_input['modradsize'] = True
    if args.scfconv:
        cml_input['scfconv'] = args.scfconv
    if args.cosmo:
        cml_input['cosmo'] = args.cosmo
    if args.dcosmors and not args.cosmo:
        cml_input['cosmo'] =  settings.solvent_dcosmors.get(args.dcosmors)[0]
    if args.dcosmors:
        cml_input['dcosmors'] = settings.solvent_dcosmors.get(args.dcosmors)[1]
    if args.gen_auxbas:
        cml_input['gen_auxbas'] = args.gen_auxbas
    if args.gen_bas:
        cml_input['gen_bas'] = args.gen_bas
    if args.noauxg:
        cml_input['noauxg'] = args.noauxg
    if args.noopt:
        cml_input['noopt'] = args.noopt
    return cml_input

jbasis = 'universal'
ecp = 'def2-ecp'
# not considered :
#intmem = 1000
#twoint = 1000
# marij
# intmem
# fp

print("tmprep.py - Command line input for TURBOMOLE V {} FB 2021".format(version))

# read coord and get number of atoms, elements
coord_path = os.path.join(os.getcwd(), 'coord')
if not os.path.isfile(coord_path):
    print("coord file not found. Going to exit.")
    sys.exit(1)
else:
    elements = []
    with open(coord_path, 'r') as inp:
        coord_data = inp.readlines()
    if '$coord' not in coord_data[0]:
        print("Corrupt coord file! Going to exit.")
        sys.exit(1)
    for line in coord_data[1:]:
        if '$' in line:
            break
        elements.append(line.strip().split()[3].lower())

# remove old data
for file in ('mos', 'alpha', 'beta', 'basis', 'auxbasis', 'control'):
    file_path = os.path.join(os.getcwd(), file)
    if os.path.isfile(file_path):
        os.remove(file_path)

cml_input = cml(Settings())
do = Settings(**cml_input)

nat = len(elements)
element_occurence = Counter(elements)
element_electrons={'h':1, 'he':2, 'li':3, 'be':4, 'b':5, 'c':6,
                   'n':7, 'o':8, 'f':9, 'ne':10, 'na':11, 'mg':12,
                   'al':13, 'si':14, 'p':15, 's':16, 'cl':17, 'ar':18,
                   'k':19, 'ca':20, 'sc':21, 'ti':22, 'v':23, 'cr':24,
                   'mn':25, 'fe':26, 'co':27, 'ni':28, 'cu':29, 'zn':30,
                   'ga':31, 'ge':32, 'as':33, 'se':34, 'br':35, 'kr':36,
                   'rb':37, 'sr':38, 'y':39, 'zr':40, 'nb':41, 'mo':42,
                   'tc':43, 'ru':44, 'rh':45, 'pd':46, 'ag':47, 'cd':48,
                   'in':49, 'sn':50, 'sb':51, 'te':52, 'i':53, 'xe':54,
                   'cs':55, 'ba':56, 'la':57, 'ce':58, 'pr':59, 'nd':60,
                   'pm':61, 'sm':62, 'eu':63, 'gd':64, 'tb':65, 'dy':66,
                   'ho':67, 'er':68, 'tm':69, 'yb':70, 'lu':71, 'hf':72,
                   'ta':73, 'w':74,  're':75, 'os':76, 'ir':77, 'pt':78,
                   'au':79, 'hg':80, 'tl': 81, 'pb': 82, 'bi':83, 'po':84,
                   'at':85, 'rn': 86, 'fr':87, 'ra':88, 'ac':89, 'th':90,
                   'pa':91, 'u':92, 'np':93, 'pu':94, 'am':95, 'cm':96,
                   'bk': 97, 'cf':98, 'es':99, 'fm':100, 'md':101, 'no':102,
                   'lr': 103}
if not all(element in element_electrons.keys() for element in element_occurence.keys()):
    print("There are some elements that are unknown:")
    for element in element_occurence.keys():
        if element not in element_electrons.keys():
            print("NOT KNOWN: {}".format(element))
            error_logical = True
    if error_logical:
        sys.exit(1)

ecp28=['rb', 'sr', 'y', 'zr', 'nb', 'mo', 'tc', 'ru', 'rh', 'pd', 'ag', 'cd',
       'in', 'sn', 'sb', 'te', 'i', 'xe', 'ce', 'pr', 'nd', 'pm', 'sm', 'eu',
       'gd', 'tb', 'dy', 'ho', 'er', 'tm', 'yb', 'lu'
      ]
ecp46=['cs', 'ba', 'la']
ecp60=['hf', 'ta', 'w', 're', 'os', 'ir', 'pt', 'au', 'hg', 'tl', 'pb', 'bi',
       'po', 'at', 'rn'
       ]

all_ecp = ecp28 + ecp46 + ecp60
subscript = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
superscript = str.maketrans("0123456789+-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻")
sumformula = []
nall_electrons = 0
necp_electrons = 0
for key, value in element_occurence.items():
    sumformula.append(key[0].upper()+key[1:])
    sumformula.append(str(value))
    nall_electrons += value * element_electrons.get(key, 0)
    if key in ecp28:
        necp_electrons += value * (element_electrons.get(key, 0) - 28)
    elif key in ecp46:
        necp_electrons += value * (element_electrons.get(key, 0) - 46)
    elif key in ecp60:
        necp_electrons += value * (element_electrons.get(key, 0) - 60)
    else:
        necp_electrons += value * element_electrons.get(key, 0)

nall_electrons += -do.charge
necp_electrons += -do.charge

print('chemical formula: ', "".join(sumformula).translate(subscript))
print('number of atoms: {}'.format(nat))
if do.charge != 0:
    print('charge: {}'.format(do.charge))
if do.unpaired != 0:
    print('unpaired number of electrons: {}'.format(do.unpaired))
if do.symmetry is not None:
    print("symmetry: {}".format(do.symmetry))
print("number of all electrons: {}".format(nall_electrons))
print("number of considered electrons: {}".format(necp_electrons))
if nat == 1:
    print("TMPREP CHOICE: found one atom only, setting symmetry: c1")
    do.symmetry = 'c1'


# definition of composite methods:
#r2SCAN-3c
if do.functional in ('r2scan-3c', 'R2SCAN-3c', 'r2SCAN-3c'):
    # uses gcp D4 def2-mTZVP
    do.functional = 'r2scan-3c'
    if not do.modbasis:
        do.basis = 'def2-mTZVPP'
    if not do.modgrid:
        do.grid = 'm4'
    if not do.modradsize:
        do.radsize = 8
    if not do.moddisp:
        do.disp = 'disp4'
#PBEh-3c
elif do.functional in ('pbeh-3c', 'pbeh3c'):
    do.functional = 'pbeh-3c'
    if not do.modbasis:
        do.basis = 'def2-mSVP'
    if not do.modgrid:
        do.grid = 'm4'
    if not do.moddisp:
        disp = 'disp3 -bj'
    if not do.extol:
        do.extol = 2.5
#B97-3c
elif do.functional in ('b97-3c', 'b973c'):
    if do.functional == 'b973c':
        print("Using the slower b973c because of XCFun "
              "e.g. for hessian calculation")
    if not do.modbasis:
        do.basis = 'def2-mTVP'
    if not do.modgrid:
        do.grid = 'm4'
    if not do.moddisp:
        do.disp = 'disp3 -bj -abc'
#PBE-3c
elif do.functional in ('pbe3c', 'pbe-3c'):
    # uses gcp D3(BJ)ATM
    print('Setting up PBE-3c calculation (NOT PBEh-3c)')
    do.functional = 'pbe'
    do.gcp = 'gcp dft/sv(p)'
    if not do.modgrid:
        do.grid = 'm3'
    if not do.modbasis:
        do.basis = 'def2-mSVP'
    if not do.moddisp:
        do.disp = 'disp3 -bj -abc'
    if not do.extol:
        do.extol = 2.5
#B3-LYP-3c
elif do.functional in ('b3-lyp-3c', 'b3lyp-3c'):
    print('Setting up B3-LYP-3c calculation')
    do.functional = 'b3-lyp'
    do.gcp = 'gcp dft/sv(p)'
    if not do.modgrid:
        do.grid = 'm4'
    if not do.modbasis:
        do.basis = 'def2-mSVP'
    if not do.moddisp:
        do.disp = 'disp3 -bj -abc'
    if not do.extol:
        do.extol = 2.5
#PBE0-3c
elif do.functional in ('pbe0-3c', 'pbe03c'):
    print('Setting up PBEO-3c calculation')
    do.functional = 'pbe0'
    do.gcp = 'gcp dft/sv(p)'
    if not do.modgrid:
        do.grid = 'm4'
    if not do.modbasis:
        do.basis = 'def2-mSVP'
    if not do.moddisp:
        do.disp = 'disp3 -bj -abc'
    if not do.extol:
        do.extol = 2.5
# KT1 KT2
if do.functional in ('kt2', 'kt1'):
    setattr(do, 'xcfun', True)
  
if do.novdw:
    do.disp = ""

if do.functional and not do.wave_func:
    print("Settings: {}/{} scfconv {} grid {}".format(
        do.functional, do.basis, do.scfconv, do.grid)
        )
else:
    print("Settings: {}/{} scfconv {}".format(
        do.wave_func, do.basis, do.scfconv)
        )

#------- write basis file-------------------------------------------------------
if do.gen_bas:
    basis_element_data = {}
    ecp_element_data = {}
    for element in element_occurence.keys():
        do_ecp = False
        element = element.lower()
        if element in all_ecp:
            do_ecp = True
        basis_set = do.basis
        comment = 'comment'

        if os.path.isfile(os.path.join(os.path.expandvars("$TURBODIR"), "basen", element)):
            basis_file = os.path.join(os.path.expandvars("$TURBODIR"), "basen", element)
            print("Using basis ({}) from {}".format(element, os.path.join(os.path.expandvars('$TURBODIR'), 'basen', element)))
        else:
            print("ERROR: the basis set file can not be found!")

        with open(basis_file, 'r', encoding="ISO-8859-1") as inp:
            data = inp.readlines()

        # check if the basis set name can be found in the basis set file
        found = False
        for line in data:
            if do.basis in line:
                found = True
                break
        if not found:
            print("The basis set {} can not be found for element {}!".format(do.basis, element))
            print("ERROR: can not create the basis and auxbasis file!!!!")
            do.gen_auxbas = False
            do.gen_bas = False
            found = False
            break

        found = False
        start = False
        find_basis = []
        basis_set = do.basis
        # read first iteration
        for line in data: 
            if basis_set in line and not '#' in line and not '->' in line:
                found = True
            if found:
                if element in line.split()[0] and not '->' in line and basis_set in line:
                    initial = line.strip()
                if '*' in line and not start and element not in line:
                    start=True
                    continue
                elif '*' in line and start:
                    break
                if start and '->' in line:
                    find_basis.append(line.strip().strip('-> '))
        if not find_basis:
            find_basis.append(initial)
        found = False
        start = False

        # read second iteration (and insert at appropriate position)
        count = 0
        for basis_set in find_basis:
            for line in data: 
                if basis_set in line and not '#' in line and not '->' in line:
                    found = True
                if found:
                    if element in line.split()[0] and not '->' in line and basis_set in line:
                        initial = line.strip()
                    if '*' in line and not start and element not in line:
                        start=True
                        continue
                    elif '*' in line and start:
                        break
                    if start and '->' in line:
                        count+=1
                        find_basis.insert(count,line.strip().strip('-> '))
            found = False
            start = False
        basis_set = do.basis

        # read basis data for found basis sets
        found = False
        start = False
        basis_data = []
        for basis in find_basis:
            for line in data:
                if basis in line and not '#' in line and not '->' in line:
                    if basis == line.strip():
                        found = True
                        continue
                if found:
                    if '*' in line and not start and element not in line and not '#' in line:
                        start=True
                        continue
                    elif '*' in line and start:
                        found = False
                        start = False
                        break
                    if start and '#' not in line and '->' not in line:
                        basis_data.append(line.strip())
        basis_element_data[element] = basis_data
        
        found = False
        start = False
        ecp_data = []
        if do_ecp:
            for line in data:
                if ecp in line and not '#' in line and not '->' in line:
                    found = True
                if found:
                    if '*' in line and not start and element not in line:
                        start=True
                        continue
                    elif '*' in line and start:
                        found = False
                        start = False
                        break
                    if start and '#' not in line:
                        ecp_data.append(line.strip())
            ecp_element_data[element] = ecp_data
    
    if do.gen_bas:
        # write basis file:
        with open('basis', 'w') as out:
            out.write("$basis\n")
            for element in element_occurence.keys():
                out.write("*\n")
                out.write("{} {}\n".format(element, basis_set))
                out.write("# {}\n".format(comment))
                out.write("*\n")
                for line in basis_element_data[element]:
                    tmp = line.split()
                    if len(tmp) == 2:
                        try:
                            a = "{:.12}".format(float(tmp[0]))
                            b = "{:.12}".format(float(tmp[1]))
                            out.write("  {:<12}     {:<12}\n".format(a, b))
                        except (ValueError,Exception) as e:
                            out.write("   {}  {}\n".format(tmp[0],tmp[1]))
                    else:
                        out.write(line+'\n')
            out.write("*\n")
            if ecp_element_data:
                out.write("$ecp\n")
            for element in ecp_element_data.keys():
                out.write("*\n")
                out.write("{} {}\n".format(element, ecp))
                out.write("*\n")
                for line in ecp_element_data[element]:
                    if 'ncore' in line:
                        out.write(line+'\n')
                        out.write("#        coefficient   r^n          exponent\n")
                        continue
                    if len(line.split()) >= 3:
                        out.write("    {}\n".format(line))
                    else:
                        out.write(line+'\n')
            if ecp_element_data:
                out.write("*\n")
            out.write("$end\n")
#----END write basis file-------------------------------------------------------



# write auxbasis file ----------------------------------------------------------
if do.gen_auxbas:
    jbasis_element_data = {}
    for element in element_occurence.keys():
        jbasis_set = 'universal'
        comment = 'comment'

        if os.path.isfile(os.path.join(os.path.expandvars("$TURBODIR"), "jbasen", element)):
            jbasis_file = os.path.join(os.path.expandvars("$TURBODIR"), "jbasen", element)
            print("Using jbasis ({}) from {}".format(element, os.path.join(os.path.expandvars('$TURBODIR'), 'jbasen', element)))
        else:
            print("ERROR: the jbasis set file can not be found!")

        with open(jbasis_file, 'r', encoding="ISO-8859-1") as inp:
            data = inp.readlines()

        found = False
        start = False
        skip_x = 0
        jbasis_data = []
        for line in data:
            if skip_x > 0:
                skip_x += -1
                continue
            if jbasis_set in line and not '#' in line and not '->' in line:
                if line.strip().split()[1] == jbasis:
                    found = True
                    jbasis_data.append(line.strip())
                elif jbasis_set in line and 'ecp' in line:
                    found = True
                    jbasis_data.append(line.strip())
            if found:
                if '*' in line and not start and element not in line:
                    start=True
                    continue
                elif '*' in line and start:
                    found = False
                    start = False
                    break
                if start and '#' not in line:
                    if do.noauxg and 'g' in line and len(line.strip().split()) == 2:
                        skip_x = int(line.strip().split()[0])
                        continue
                    else:
                        jbasis_data.append(line.strip())
        jbasis_element_data[element] = jbasis_data

    # write auxbasis file:
    with open('auxbasis', 'w') as out:
        out.write("$jbas\n")
        for element in element_occurence.keys():
            out.write("*\n")
            if not jbasis_element_data[element]:
                print("ERROR: jbas not found for element: {}".format(element))
                continue
            out.write(jbasis_element_data[element][0]+'\n')
            out.write("# {}\n".format(comment))
            out.write("*\n")
            for line in jbasis_element_data[element][1:]:
                tmp = line.split()
                if len(tmp) == 2:
                    try:
                        a = "{:.12}".format(float(tmp[0]))
                        b = "{:.12}".format(float(tmp[1]))
                        out.write("  {:<12}    {:<12}\n".format(a, b))
                    except (ValueError,Exception) as e:
                        out.write("   {}  {}\n".format(tmp[0], tmp[1]))
                else:
                    out.write(line+'\n')
        out.write("*\n")
        out.write("$end\n")
#----END write auxbasis file-------------------------------------------------------


outputfile = os.path.join(os.getcwd(), 'control')
with open(outputfile, 'w', newline='\n') as out:
    out.write("$title \n")
    out.write("$coord file=coord\n")
    if do.unpaired !=0:
        out.write("$eht charge={} unpaired={}\n".format(do.charge, do.unpaired))
    else:
        out.write("$eht charge={}\n".format(do.charge))
    if do.symmetry is not None:
        out.write("$symmetry {}\n".format(str(do.symmetry).lower()))

    out.write("$atoms\n")

    for element in element_occurence.keys():
        # element position in coord:
        positions = []
        i = 1
        for atom in elements:
            if atom == element:
                positions.append(i)
            i+=1
        positions.sort()
        # shorten sequences of atoms
        short_positions = []
        first = last = positions[0]
        for n in positions[1:]:
            if n - 1 == last:
                last = n
            else:
                if first == last:
                    short_positions.append(str(first))
                else:
                    short_positions.append("{}-{}".format(first, last))
                first = last = n
        if first == last:
            short_positions.append(str(first))
        else:
            short_positions.append("{}-{}".format(first, last))
        # break long lines
        newpos = []
        tmp = []
        lenpos = 0
        for item in short_positions:
            lenitem = len(item)+1
            if (lenpos + lenitem) > 70:
                newpos.append(','.join(tmp))
                lenpos = lenitem
                tmp = [item,]
            else:
                tmp.append(item)
                lenpos+=lenitem
        if tmp:
            newpos.append(','.join(tmp))
        for count, item in enumerate(newpos,0):
            if count == 0:
                out.write('{} {} \ \n'.format(element, item))
            else:
                out.write('{} {} \ \n'.format(" ", item))
        #out.write('{} {} \ \n'.format(element, ','.join(short_positions)))
        out.write('  basis ={} {} \ \n'.format(element, do.basis))
        if element in all_ecp:
            if element in ecp28:
                out.write('  jbas  ={} {}\n'.format(element, 'universal-ecp-28'))
            elif element in ecp46:
                out.write('  jbas  ={} {}\n'.format(element, 'universal-ecp-46'))
            elif element in ecp60:
                out.write('  jbas  ={} {}\n'.format(element, 'universal-ecp-60'))
            out.write('  ecp   ={} {}\n'.format(element, ecp))
        else:
            out.write('  jbas  ={} {}\n'.format(element, 'universal'))
    if do.gen_bas:
        out.write('$basis    file=basis \n')
        if any(element in all_ecp for element in element_occurence.keys()):
            out.write('$ecp      file=basis \n')
    if do.gen_auxbas:
        out.write('$jbas    file=auxbasis \n')
    out.write("$scforbitalshift automatic=0.1 \n")
    if do.functional and not do.wave_func:
        out.write("$dft\n")
        if do.xcfun: # currently only kt1 kt2
            # can be improved for sure, currently just as a demonstration
            out.write("  functional xcfun {}\n".format('set-gga'))
            out.write("  functional xcfun {} {} \n".format(do.functional, '1.0'))
        else:
            out.write("  functional {}\n".format(do.functional))
        out.write("  gridsize {}\n".format(do.grid))
        if do.radsize is not None:
            out.write("  radsize {}\n".format(do.radsize))
    out.write("$scfconv {}\n".format(do.scfconv))
    if do.ri:
        out.write("$rij\n")
        out.write("$ricore {}\n".format(do.ricore))
    out.write("$maxcor {}\n".format(do.maxcor))
    if do.rpacor:
        out.write("$rpacor {}\n".format(do.rpacor))
    if do.disp:
        out.write("${}\n".format(do.disp))
    if do.gcp:
        out.write("${}\n".format(do.gcp))
    if do.extol:
        out.write("$extol  {}\n".format(do.extol))
    out.write("$scfiterlimit {}\n".format(do.scfiterlimit))
    out.write("$thize {} \n".format(do.thize))
    out.write("$thime {} \n".format(do.thime))
    if do.noopt:
        out.write("$optimize \n")
        out.write("  internal   off \n")
        out.write("  redundant  off \n")
        out.write("  cartesian  on  \n")
    out.write('$statpt \n')
    out.write('  itrvec 0 \n')
    out.write('  tradius 0.3 \n')
    out.write('  threchange  5.0d-7 \n')
    out.write('  thrrmsgrad  5.0d-5 \n')
    out.write('  thrmaxdispl 1.0d-3 \n')
    out.write('  thrrmsdispl 1.0d-3 \n')
    # COSMO:
    if do.cosmo:
        out.write("$cosmo \n")
        out.write("  epsilon= {:.2f} \n".format(float(do.cosmo)))
    # DCOSMO-RS:
    if do.dcosmors:
        out.write("{}\n".format(do.dcosmors))
    # redirect output
    out.write("$energy    file=energy \n")
    out.write("$grad    file=gradient \n")
    ### terminate control file
    out.write("$end\n")

print("\nTMPREP: All done!")
