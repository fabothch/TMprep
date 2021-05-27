#!/usr/bin/env python3

# Copyright (C) 2020 Fabian Bohle
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

# To do:
# clear error codes!


import os
import sys
import argparse
from collections import Counter

version = "0.1.4"

def read_chrg(default=0):
    # READ .CHRG
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
    # READ .UHF
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
    # READ .SYM
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

def cml(internal_defaults, solvent_dcosmors, argv=None):
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
        choices=list(solvent_dcosmors.keys()),
        metavar="",
        help=("Add DCOSMO-RS to control. Usage: -dcosmors [solvent]. "
            "Options are {}. It is not necessary to additionally use -cosmo".format(', '.join(list(solvent_dcosmors.keys())))),
    )
    group1.add_argument(
        "-noopt",
        "--noopt",
        dest="noopt",
        action="store_true",
        required=False,
        help=("Use cartesian coordinates."),
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
    args = parser.parse_args(argv)
    if args.hf:
        args.wave_func ="hf"
        try:
            delattr(args, "hf")
        except Exception as e:
            print(e)
    if args.wave_func:
        args.wave_func = getattr(args, "wave_func").lower()
    if args.functional:
        args.functional = getattr(args, "functional").lower()
    if args.basis:
        setattr(args, 'modbasis', True)
    if args.grid:
        setattr(args, 'modgrid', True)
    if not args.charge:
        setattr(args,'charge', read_chrg())
    if not args.unpaired:
        setattr(args,'unpaired', read_uhf())
    if not args.symmetry:
        setattr(args, 'symmetry', read_sym())
    if (args.d4 or
        args.d3zero or
        args.d3 or
        args.d3atm or
        args.donl or
        args.novdw
        ):
        setattr(args, 'moddisp', True)
    # dispersion corrections:
    if args.d3zero:
        setattr(args, 'disp', 'disp3')
    if args.d3:
        setattr(args, 'disp', 'disp3 -bj')
    if args.d3atm:
        setattr(args, 'disp', 'disp3 -bj -abc')
    if args.d4:
        setattr(args, 'disp', 'disp4')
    if args.donl:
        setattr(args, 'disp', 'donl ')
    if args.nori:
        setattr(args, 'ri', False)
    for key, value in internal_defaults.items():
        if getattr(args, key, None) is None:
            setattr(args, key, value)
    if args.novdw:
        args.disp = ""
    if args.radsize is not None:
        args.modradsize = True
    if args.dcosmors and not args.cosmo:
        setattr(args, 'cosmo', solvent_dcosmors.get(args.dcosmors)[0])
    return args

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


internal_defaults = {
    "symmetry": None, # not c1
    "basis": "def2-mSVP", 
    'functional': 'pbeh-3c',
    'grid': 'm4',
    'ri': True,
    'ricore': 4000,
    'maxcor': 4000,
    'radsize': None,
    'scfconv': 7,
    'charge': 0,
    'unpaired': 0,
    'gcp': "",
    'disp': 'disp3 -bj',
    'modgrid' : False,
    'modradsize': False,
    'modbasis': False,
    'moddisp': False,
    'novdw' : False,
    'rpacor': 1000,
    'twoint': 1000,
    'scfiterlimit': 125,
    'cosmo': None,
    'xcfun': False,
    'thime': 4,
    'thize': 0.0000001,
    'extol': None,
    'noauxg' : True,
}

path_cefinerc =os.path.expanduser('~/.cefinerc')
if os.path.isfile(path_cefinerc):
    # read defaults concerning memory to internal_defaults dictionary
    with open(path_cefinerc, 'r') as inp:
        data = inp.readlines()
        for line in data:
            if 'ricore' in line:
                try:
                    internal_defaults['ricore'] = int(line.split()[1])
                except (ValueError, TypeError):
                    pass
            if 'maxcor' in line:
                try:
                    internal_defaults['maxcor'] = int(line.split()[1])
                except (ValueError, TypeError):
                    pass
            if 'rpacor' in line:
                try:
                    internal_defaults['rpacor'] = int(line.split()[1])
                except (ValueError, TypeError):
                    pass
            if 'func' in line:
                try:
                    internal_defaults['functional'] = line.split()[1].strip()
                except (ValueError, TypeError):
                    pass
            if 'bas' in line:
                try:
                    internal_defaults['basis'] = line.split()[1].strip()
                except (ValueError, TypeError):
                    pass
            if 'grid' in line:
                try:
                    internal_defaults['grid'] = line.split()[1].strip()
                except (ValueError, TypeError):
                    pass
            if 'scfconv' in line:
                try:
                    internal_defaults['scfconv'] = int(line.split()[1])
                except (ValueError, TypeError):
                    pass

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

args = cml(internal_defaults, solvent_dcosmors)

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

nall_electrons += -args.charge
necp_electrons += -args.charge

print('chemical formula: ', "".join(sumformula).translate(subscript))
print('number of atoms: {}'.format(nat))
if args.charge != 0:
    print('charge: {}'.format(args.charge))
if args.unpaired != 0:
    print('unpaired number of electrons: {}'.format(args.unpaired))
if args.symmetry is not None:
    print("symmetry: {}".format(args.symmetry))
print("number of all electrons: {}".format(nall_electrons))
print("number of considered electrons: {}".format(necp_electrons))
if nat == 1:
    print("TMPREP CHOICE: found one atom only, setting symmetry: c1")
    args.symmetry = 'c1'


# definition of composite methods:
#r2SCAN-3c
if args.functional in ('r2scan-3c', 'R2SCAN-3c', 'r2SCAN-3c'):
    # uses gcp D4 def2-mTZVP
    args.functional = 'r2scan-3c'
    if not args.modbasis:
        args.basis = 'def2-mTZVPP'
    if not args.modgrid:
        args.grid = 'm4'
    if not args.modradsize:
        args.radsize = 8
    if not args.moddisp:
        args.disp = 'disp4'
#PBEh-3c
elif args.functional in ('pbeh-3c', 'pbeh3c'):
    args.functional = 'pbeh-3c'
    if not args.modbasis:
        args.basis = 'def2-mSVP'
    if not args.modgrid:
        args.grid = 'm4'
    if not args.moddisp:
        disp = 'disp3 -bj'
    if not args.extol:
        args.extol = 2.5
#B97-3c
elif args.functional in ('b97-3c', 'b973c'):
    if args.functional == 'b973c':
        print("Using the slower b973c because of XCFun "
              "e.g. for hessian calculation")
    if not args.modbasis:
        args.basis = 'def2-mTVP'
    if not args.modgrid:
        args.grid = 'm4'
    if not args.moddisp:
        args.disp = 'disp3 -bj -abc'
#PBE-3c
elif args.functional in ('pbe3c', 'pbe-3c'):
    # uses gcp D3(BJ)ATM
    print('Setting up PBE-3c calculation (NOT PBEh-3c)')
    args.functional = 'pbe'
    args.gcp = 'gcp dft/sv(p)'
    if not args.modgrid:
        args.grid = 'm3'
    if not args.modbasis:
        args.basis = 'def2-mSVP'
    if not args.moddisp:
        args.disp = 'disp3 -bj -abc'
    if not args.extol:
        args.extol = 2.5
#B3-LYP-3c
elif args.functional in ('b3-lyp-3c', 'b3lyp-3c'):
    print('Setting up B3-LYP-3c calculation')
    args.functional = 'b3-lyp'
    args.gcp = 'gcp dft/sv(p)'
    if not args.modgrid:
        args.grid = 'm4'
    if not args.modbasis:
        args.basis = 'def2-mSVP'
    if not args.moddisp:
        args.disp = 'disp3 -bj -abc'
    if not args.extol:
        args.extol = 2.5
#PBE0-3c
elif args.functional in ('pbe0-3c', 'pbe03c'):
    print('Setting up PBEO-3c calculation')
    args.functional = 'pbe0'
    args.gcp = 'gcp dft/sv(p)'
    if not args.modgrid:
        args.grid = 'm4'
    if not args.modbasis:
        args.basis = 'def2-mSVP'
    if not args.moddisp:
        args.disp = 'disp3 -bj -abc'
    if not args.extol:
        args.extol = 2.5
# KT1 KT2
if args.functional in ('kt2', 'kt1'):
    setattr(args, 'xcfun', True)
    
if args.novdw:
    args.disp = ""

if args.functional and not args.wave_func:
    tmp_method = args.functional
else:
    tmp_method = args.wave_func
print("Settings: {}/{} scfconv {} grid {}".format(
      tmp_method, args.basis, args.scfconv, args.grid)
      )

#------- write basis file-------------------------------------------------------
if args.gen_bas:
    basis_element_data = {}
    ecp_element_data = {}
    for element in element_occurence.keys():
        do_ecp = False
        element = element.lower()
        if element in all_ecp:
            do_ecp = True
        basis_set = args.basis
        comment = 'comment'

        if os.path.isfile(os.path.join(os.path.expandvars("$TURBODIR"), "basen", element)):
            basis_file = os.path.join(os.path.expandvars("$TURBODIR"), "basen", element)
            print(f"Using basis ({element}) from {os.path.join(os.path.expandvars('$TURBODIR'), 'basen', element)}")
        else:
            print(f"ERROR: the basis set file can not be found!")

        with open(basis_file, 'r', encoding="ISO-8859-1") as inp: 
            data = inp.readlines()

        found = False
        start = False
        find_basis = []
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
                    if start and '#' not in line:
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

    # write basis file:
    with open('basis', 'w') as out:
        out.write("$basis\n")
        for element in element_occurence.keys():
            out.write("*\n")
            out.write(f"{element} {basis_set}\n")
            out.write(f"# {comment}\n")
            out.write("*\n")
            for line in basis_element_data[element]:
                tmp = line.split()
                if len(tmp) == 2:
                    try:
                        float(tmp[1])
                        out.write(f"  {tmp[0]}    {tmp[1]}\n")
                    except (ValueError,Exception) as e:
                        out.write(f"   {tmp[0]}  {tmp[1]}\n")
                else:
                    out.write(line+'\n')
        out.write("*\n")
        if ecp_element_data:
            out.write("$ecp\n")
        for element in ecp_element_data.keys():
            out.write("*\n")
            out.write(f"{element} {ecp}\n")
            out.write("*\n")
            for line in ecp_element_data[element]:
                if 'ncore' in line:
                    out.write(line+'\n')
                    out.write(f"#        coefficient   r^n          exponent\n")
                    continue
                if len(line.split()) >= 3:
                    out.write(f"    {line}\n")
                else:
                    out.write(line+'\n')
        if ecp_element_data:
            out.write("*\n")
        out.write("$end\n")
#----END write basis file-------------------------------------------------------



# write auxbasis file ----------------------------------------------------------
if args.gen_auxbas:
    jbasis_element_data = {}
    for element in element_occurence.keys():
        jbasis_set = 'universal'
        comment = 'comment'

        if os.path.isfile(os.path.join(os.path.expandvars("$TURBODIR"), "jbasen", element)):
            jbasis_file = os.path.join(os.path.expandvars("$TURBODIR"), "jbasen", element)
            print(f"Using jbasis ({element}) from {os.path.join(os.path.expandvars('$TURBODIR'), 'jbasen', element)}")
        else:
            print(f"ERROR: the jbasis set file can not be found!")

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
                    if args.noauxg and 'g' in line and len(line.strip().split()) == 2:
                        skip_x = int(line.strip().split()[0])
                        continue
                    else:
                        jbasis_data.append(line.strip())
        jbasis_element_data[element] = jbasis_data

    # write basis file:
    with open('auxbasis', 'w') as out:
        out.write("$jbas\n")
        for element in element_occurence.keys():
            out.write("*\n")
            #out.write(f"{element} {basis_set}\n")
            if not jbasis_element_data[element]:
                print("ERROR: jbas not found for element: {}".format(element))
                continue
            out.write(jbasis_element_data[element][0]+'\n')
            out.write(f"# {comment}\n")
            out.write("*\n")
            for line in jbasis_element_data[element][1:]:
                tmp = line.split()
                if len(tmp) == 2:
                    try:
                        float(tmp[1])
                        out.write(f"  {tmp[0]}    {tmp[1]}\n")
                    except (ValueError,Exception) as e:
                        out.write(f"   {tmp[0]}  {tmp[1]}\n")
                else:
                    out.write(line+'\n')
        out.write("*\n")
        out.write("$end\n")
#----END write auxbasis file-------------------------------------------------------


outputfile = os.path.join(os.getcwd(), 'control')
with open(outputfile, 'w', newline='\n') as out:
    out.write("$title \n")
    out.write("$coord file=coord\n")
    if args.unpaired !=0:
        out.write("$eht charge={} unpaired={}\n".format(args.charge, args.unpaired))
    else:
        out.write("$eht charge={}\n".format(args.charge))
    if args.symmetry is not None:
        out.write("$symmetry {}\n".format(str(args.symmetry).lower()))

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
        out.write('  basis ={} {} \ \n'.format(element, args.basis))
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
    if args.gen_bas:
        out.write('$basis    file=basis \n')
        if any(element in all_ecp for element in element_occurence.keys()):
            out.write('$ecp      file=basis \n')
    if args.gen_auxbas:
        out.write('$jbas    file=auxbasis \n')
    out.write("$scforbitalshift automatic=0.1 \n")
    if args.functional and not args.wave_func:
        out.write("$dft\n")
        if args.xcfun: # currently only kt1 kt2
            # can be improved for sure, currently just as a demonstration
            out.write("  functional xcfun {}\n".format('set-gga'))
            out.write("  functional xcfun {} {} \n".format(args.functional, '1.0'))
        else:
            out.write("  functional {}\n".format(args.functional))
        out.write("  gridsize {}\n".format(args.grid))
        if args.radsize is not None:
            out.write("  radsize {}\n".format(args.radsize))
    out.write("$scfconv {}\n".format(args.scfconv))
    if args.ri:
        out.write("$rij\n")
        out.write("$ricore {}\n".format(args.ricore))
    out.write("$maxcor {}\n".format(args.maxcor))
    if args.rpacor:
        out.write("$rpacor {}\n".format(args.rpacor))
    if args.disp:
        out.write("${}\n".format(args.disp))
    if args.gcp:
        out.write("${}\n".format(args.gcp))
    if args.extol:
        out.write("$extol  {}\n".format(args.extol))
    out.write("$scfiterlimit {}\n".format(args.scfiterlimit))
    out.write("$thize {} \n".format(args.thize))
    out.write("$thime {} \n".format(args.thime))
    if args.noopt:
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
    if args.cosmo:
        out.write("$cosmo \n")
        out.write("  epsilon= {:.2f} \n".format(float(args.cosmo)))
    # DCOSMO-RS:
    if args.dcosmors:
        out.write("{}\n".format(solvent_dcosmors.get(args.dcosmors)[1]))
    ### terminate control file
    out.write("$end\n")

print("\nTMPREP: All done!")
