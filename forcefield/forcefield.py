# -*- coding: utf-8 -*-
"""Main class for handling forcefields"""

import logging
import molssi_util
import os.path
import pprint

logger = logging.getLogger(__name__)


metadata = {
    'quadratic_bond': {
        'equation': ['K2*(R-R0)^2'],
        'constants': [
            ('R0', 'angstrom'),
            ('K2', 'kcal/mol/angstrom^2'),
        ],
        'topology': {
            'n_atoms': 2,
            'symmetry': 'like_bond',
            'fill': 0,
            'flip': 0
        }
    },
    'quartic_bond': {
        'equation': ['K2*(R-R0)^2 + K3*(R-R0)^3 + K4*(R-R0)^4'],
        'constants': [
            ('R0', 'angstrom'),
            ('K2', 'kcal/mol/angstrom^2'),
            ('K3', 'kcal/mol/angstrom^3'),
            ('K4', 'kcal/mol/angstrom^4'),
        ],
        'topology': {
            'n_atoms': 2,
            'symmetry': 'like_bond',
            'fill': 0,
            'flip': 0
        }
    },
    'quadratic_angle': {
        'equation': ['K2*(Theta-Theta0)^2'],
        'constants': [
            ('Theta0', 'degree'),
            ('K2', 'kcal/mol/radian^2'),
        ],
        'topology': {
            'n_atoms': 3,
            'symmetry': 'like_angle',
            'fill': 0,
            'flip': 0
        }
    },
    'quartic_angle': {
        'equation': ['K2*(Theta-Theta0)^2 + K3*(Theta-Theta0)^3'
                     '+ K4*(Theta-Theta0)^4'],
        'constants': [
            ('Theta0', 'degree'),
            ('K2', 'kcal/mol/radian^2'),
            ('K3', 'kcal/mol/radian^3'),
            ('K4', 'kcal/mol/radian^4'),
        ],
        'topology': {
            'n_atoms': 3,
            'symmetry': 'like_angle',
            'fill': 0,
            'flip': 0
        }
    },
    'torsion_1': {
        'equation': ['KPhi * [1 + cos(n*Phi - Phi0)]'],
        'constants': [
            ('KPhi', 'kcal/mol'),
            ('n', ''),
            ('Phi0', 'degree'),
        ],
        'topology': {
            'n_atoms': 4,
            'symmetry': 'like_torsion',
            'fill': 0,
            'flip': 0
        }
    },
    'torsion_3': {
        'equation': ['V1 * [1 + cos(Phi - Phi0_1)]'
                     ' + V2 * [1 + cos(2*Phi - Phi0_2)]'
                     ' + V3 * [1 + cos(3*Phi - Phi0_3)]'],
        'constants': [
            ('V1', 'kcal/mol'),
            ('Phi0_1', 'degree'),
            ('V2', 'kcal/mol'),
            ('Phi0_2', 'degree'),
            ('V3', 'kcal/mol'),
            ('Phi0_3', 'degree'),
        ],
        'topology': {
            'n_atoms': 4,
            'symmetry': 'like_torsion',
            'fill': 0,
            'flip': 0
        }
    },
    'wilson_out_of_plane': {
        'equation': ['K*(Chi - Chi0)^2'],
        'constants': [
            ('K', 'kcal/mol/radian^2'),
            ('Chi0', 'degree'),
        ],
        'topology': {
            'n_atoms': 4,
            'symmetry': 'like_oop',
            'fill': 0,
            'flip': 0
        }
    },
    'nonbond(9-6)': {
        'equation': [
            'eps(ij) [2(r(ij)*/r(ij))**9 - 3(r(ij)*/r(ij))**6]',
            'r(ij) = [(r(i)**6 + r(j)**6))/2]**(1/6)',
            'eps(ij) = 2 * sqrt(eps(i) * eps(j)) * '
            'r(i)^3 * r(j)^3/[r(i)^6 + r(j)^6]'
        ],
        'constants': [
            ('r', 'angstrom'),
            ('eps', 'kcal/mol'),
        ],
        'topology': {
            'n_atoms': 1,
            'symmetry': 'none',
            'fill': 0,
            'flip': 0
        }
    },
    'bond-bond': {
        'equation': ["K*(R-R0)*(R'-R0')"],
        'constants': [
            ('K', 'kcal/mol/angstrom^2'),
        ],
        'topology': {
            'n_atoms': 3,
            'symmetry': 'like_angle',
            'fill': 0,
            'flip': 0
        }
    },
    'bond-bond_1_3': {
        'equation': ["K*(R-R0)*(R'-R0')"],
        'constants': [
            ('K', 'kcal/mol/angstrom^2'),
        ],
        'topology': {
            'n_atoms': 4,
            'symmetry': 'like_torsion',
            'fill': 0,
            'flip': 0
        }
    },
    'bond-angle': {
        'equation': ["K*(R-R0)*(Theta-Theta0)"],
        'constants': [
            ('K', 'kcal/mol/angstrom/radian'),
        ],
        'topology': {
            'n_atoms': 3,
            'symmetry': 'like_angle',
            'fill': 1,
            'flip': 1
        }
    },
    'angle-angle': {
        'equation': ["K*(Theta-Theta0)*(Theta'-Theta0')"],
        'constants': [
            ('K', 'kcal/mol/angstrom/radian'),
        ],
        'topology': {
            'n_atoms': 4,
            'symmetry': 'like_angle-angle',
            'fill': 0,
            'flip': 0
        }
    },
    'end_bond-torsion_3': {
        'equation': [
            '(R_L - R0_L) * (V1_L * [1 + cos(Phi - Phi0_1)]'
            ' + V2_L * [1 + cos(2*Phi - Phi0_2)]'
            ' + V3_L * [1 + cos(3*Phi - Phi0_3)])',
            '(R_R - R0_R) * (V1_R * [1 + cos(Phi - Phi0_1)]'
            ' + V2_R * [1 + cos(2*Phi - Phi0_2)]'
            ' + V3_R * [1 + cos(3*Phi - Phi0_3)])',
        ],
        'constants': [
            ('V1_L', 'kcal/mol'),
            ('V2_L', 'kcal/mol'),
            ('V3_L', 'kcal/mol'),
            ('V1_R', 'kcal/mol'),
            ('V2_R', 'kcal/mol'),
            ('V3_R', 'kcal/mol'),
        ],
        'topology': {
            'n_atoms': 4,
            'symmetry': 'like_torsion',
            'fill': 3,
            'flip': 3
        }
    },
    'middle_bond-torsion_3': {
        'equation': [
            '(R_M - R0_M) * (V1 * [1 + cos(Phi - Phi0_1)]'
            ' + V2 * [1 + cos(2*Phi - Phi0_2)]'
            ' + V3 * [1 + cos(3*Phi - Phi0_3)])',
        ],
        'constants': [
            ('V1', 'kcal/mol'),
            ('V2', 'kcal/mol'),
            ('V3', 'kcal/mol'),
        ],
        'topology': {
            'n_atoms': 4,
            'symmetry': 'like_torsion',
            'fill': 0,
            'flip': 0
        }
    },
    'angle-torsion_3': {
        'equation': [
            '(Theta_L - Theta0_L) * (V1_L * [1 + cos(Phi - Phi0_1)]'
            ' + V2_L * [1 + cos(2*Phi - Phi0_2)]'
            ' + V3_L * [1 + cos(3*Phi - Phi0_3)])',
            '(Theta_R - Theta0_R) * (V1_R * [1 + cos(Phi - Phi0_1)]'
            ' + V2_R * [1 + cos(2*Phi - Phi0_2)]'
            ' + V3_R * [1 + cos(3*Phi - Phi0_3)])',
        ],
        'constants': [
            ('V1_L', 'kcal/mol'),
            ('V2_L', 'kcal/mol'),
            ('V3_L', 'kcal/mol'),
            ('V1_R', 'kcal/mol'),
            ('V2_R', 'kcal/mol'),
            ('V3_R', 'kcal/mol'),
        ],
        'topology': {
            'n_atoms': 4,
            'symmetry': 'like_torsion',
            'fill': 3,
            'flip': 3
        }
    },
    'angle-angle-torsion_1': {
        'equation': [
            'K * (Theta_L - Theta0_L) * (Theta_R - Theta0_R) * '
            '(Phi - Phi0_1'
        ],
        'constants': [
            ('K', 'kcal/mol/degree^2/degree'),
        ],
        'topology': {
            'n_atoms': 4,
            'symmetry': 'like_torsion',
            'fill': 0,
            'flip': 0
        }
    },
    'torsion-torsion_1': {
        'equation': [
            'K * cos(Phi_L) * cos(Phi_R)'
        ],
        'constants': [
            ('K', 'kcal/mol'),
        ],
        'topology': {
            'n_atoms': 5,
            'symmetry': 'like_torsion-torsion',
            'fill': 0,
            'flip': 0
        }
    },

}


class Forcefield(object):
    def __init__(self, filename=None, fftype=None):
        """
        Read, write, and use a forcefield

        The Forcefield object is the main interface for working with
        forcefields. It provides methods to read and write
        forcefields, to assign the forcefield to a molecule, as well
        as to get parameters for bonds, angles, etc.

        Args:
            filename ('str', optional): An optional filename for the forcefield
            fftype ('str', optional): An optional type for the
                forcefield. If not given and a forcefield is read, the
                code will try to divine the type of forcefield.
        """
        # the extensions and types that can be handled
        self._ff_extensions = {
            '.frc': 'Biosym',
        }
        self._ff_readers = {
            'Biosym': self._read_biosym_ff,
        }

        self._fftype = None
        self._filename = None
        self.keep_lines = False
        self.data = {}

        self.fftype = fftype
        self.filename = filename

    @property
    def filename(self):
        """'str' name of file for this forcefield.

        When the filename is set, if the file exists it is read. If it
        does not exist, it is created and initialized as a forcefield
        file. The type of the file may be given by self.fftype; if
        not, the code tries the divine the type of the forcefield. The
        default type for new forcefields is the Biosym .frc format.

        If the filename is changed the object is reset.
        """
        return self._filename

    @filename.setter
    def filename(self, value):
        if not value:
            self.clear()
            self._filename = None
        else:
            if value == self._filename:
                return

            if os.path.isfile(value):
                self.clear()
                self._filename = value
                self._read()
            else:
                self._filename = value
                self._create_file()

    @property
    def fftype(self):
        """'str' the type of forcefield to handle

        When set, the type is checked to make sure it can be handled. If not
        a RunTimeError is raised.
        """
        return self._fftype

    @fftype.setter
    def fftype(self, value):
        if not value:
            self._fftype = None
        else:
            if value not in self._ff_readers:
                raise RuntimeError(
                    "Forcefield type '{}' not supported".format(value))
            self._fftype = value

    def clear(self):
        """
        Reset the object to its initial, empty, state
        """
        # self._fftype = None  # leave the type ????
        self._filename = None
        self.data = {}

    def _read(self):
        """Read the forcefield from the file self.filename

        self.fftype gives the type of forcefield file. If it is not set
        the code attempts to divine the type from the extension and the
        first lines.
        """
        if self.fftype:
            if self.fftype in self._ff_readers:
                reader = self._ff_readers[self.fftype]
            else:
                raise RuntimeError("Forcefield type '{}' not supported".format(
                    self.fftype))
        else:
            ext = molssi_util.splitext(self.filename)
            if ext in self._ff_extensions:
                reader = self._ff_readers[self._ff_extensions[ext]]
            else:
                raise RuntimeError(
                    "Don't recognize forcefield by extension '{}'".format(ext))

        with molssi_util.Open(self.filename, 'r') as fd:
            reader(fd)

        pprint.pprint(self.data['reference'])

    def _read_biosym_ff(self, fd):
        """
        Read and parse a forcefield in Biosym's format

        Args:
            fd (file object): the file handle
        """
        self.data = {
            'forcefield': {},
        }

        try:
            # Read and process the first line, which should say
            # what the file is e.g. '!BIOSYM forcefield 1'
            line = next(fd)
            if line[0] == '!' and len(line.split()) == 3:
                file_variant, file_type, version = line[1:].split()
                logger.info(
                    "reading '{}', a {} file from {}, version {}".format(
                        self.filename, file_type, file_variant, version))
            else:
                logger.warning(
                    "reading '{}', expected a header line but got\n\t'{}'".
                    format(self.filename, line))

            # Read the rest of the file, processing the '#'
            # delimited sections
            for line in fd:
                line = line.strip()

                # Empty and comment lines
                if line == '' or line[0] == '!':
                    continue

                if line[0] == '#':
                    # fd.push()
                    words = line[1:].split()
                    section = words[0]

                    # Just ignore #end sections, as they simply close a section
                    if section == 'end':
                        continue
                    elif section == 'version':
                        self._parse_biosym_version(words)
                        continue

                    if len(words) < 2:
                        print(line)
                        logger.warning(
                            section + ' section does not have a label!\n\t' +
                            '\n\t'.join(fd.stack()))
                        label = 'missing'
                        priority = 0
                    else:
                        label = words[1]
                        if len(words) > 2:
                            priority = float(words[2])
                        else:
                            priority = 0

                    result = self._read_biosym_section(fd)

                    result['section'] = section
                    result['label'] = label
                    result['priority'] = priority

                    # If we have metadata, we can automatically parse the
                    # section
                    if section in metadata:
                        self._parse_biosym_section(result)
                    else:
                        # There should be a specific parser!
                        method_name = '_parse_biosym_' + section
                        if method_name not in Forcefield.__dict__:
                            logger.warning('Cannot find parser for ' + section)
                        else:
                            method = Forcefield.__dict__[method_name]
                            method(self, result)

        except IOError as e:
            logger.exception("Encountered I/O error opening '{}'".format(
                self.filename))
            raise

    def _read_biosym_section(self, fd):
        """
        Read the body of a section of the forcefield

        Keeps tracks of comments ('!'), annotations ('>'), and modifiers ('@'),
        returning a dictionary with them plus tte raw lines of data
        """
        result = {
            'comments': [],
            'lines': [],
            'annotations': [],
            'modifiers': []
        }

        for line in fd:
            line = line.strip()

            # Empty and comment lines
            if line == '':
                continue

            if line[0] == '!':
                result['comments'].append(line[1:])
                continue

            if line[0] == '#':
                # At the end of the section, push the line back so the
                # main reader handles it and return the dict with the
                # data
                fd.push()
                return result

            if line[0] == '>':
                # An annotation
                result['annotations'].append(line[1:])
                continue

            if line[0] == '@':
                # A modifier such as units or form
                result['modifiers'].append(line[1:])
                continue

            # Must be a line of data! :-)
            result['lines'].append(line)

    def _parse_biosym_version(self, words):
        """
        Process the 'version' section, which looks like

        #version	pcff.frc	1.0	1-July-91
        """
        pass

    def _parse_biosym_define(self, data):
        """
        Process a forcefield definition section

        #define cff91

        !Ver Ref		Function	     Label
        !--- ---    ------------------------------   ------
         1.0  1     atom_types                       cff91
         1.0  1     equivalence                      cff91
        ...
        """
        section = 'forcefield'
        ff_name = data['label']

        if section not in self.data:
            self.data[section] = {}
        self.data[section][ff_name] = data
        sections = self.data[section][ff_name]['sections'] = {}

        for line in data['lines']:
            words = line.split()
            if len(words) < 4:
                logger.error(
                    "In a define section for {}, the line is too short:".
                    format(ff_name))
                logger.error("    " + line)
            else:
                version, reference, functional_form = words[0:3]
                labels = words[3:]
                if functional_form not in sections:
                    sections[functional_form] = {}
                sections[functional_form][version] = {
                    'version': version,
                    'reference': reference,
                    'sections': labels
                }

        if not self.keep_lines:
            del data['lines']

    def _parse_biosym_atom_types(self, data):
        """
        Process the atom types

        #atom_types           cff91

        > Atom type definitions for any variant of cff91
        > Masses from CRC 1973/74 pages B-250.

        !Ver Ref  Type     Mass      Element   connection   Comment
        !--- ---  -----  ----------  -------   ----------   ---------------------------
        2.1 11   Ag     107.86800     Ag          0        Silver metal
        2.1 11   Al      26.98200     Al          0        Aluminium metal
        ...
        """  # nopep8
        section = data['section']
        label = data['label']

        if section not in self.data:
            self.data[section] = {}
        if label in self.data[section]:
            msg = "'{}' already defined in section '{}'".format(label, section)
            logger.error(msg)
            raise RuntimeError(msg)
        self.data[section][label] = data
        atom_types = self.data[section][label]['atom_types'] = {}

        for line in data['lines']:
            words = line.split()
            version, reference, atom_type, mass, element, connections = words[
                0:6]
            comment = ' '.join(words[6:])
            if atom_type not in atom_types:
                atom_types[atom_type] = {}
            if version in atom_types[atom_type]:
                msg = "atom type '{}' defined more than ".format(atom_type) + \
                      "once in section '{}'!".format(section)
                logger.error(msg)
                raise RuntimeError(msg)
            atom_types[atom_type][version] = {
                'reference': reference,
                'mass': mass,
                'element': element,
                'connections': connections,
                'comment': comment
            }

    def _parse_biosym_equivalence(self, data):
        """
        Process the atom type equivalences

        #equivalence          cff91

        !                      Equivalences
        !       ------------------------------------------
        !Ver Ref  Type   NonB   Bond   Angle  Torsion  OOP
        !--- ---  -----  -----  -----  -----  -------  -----
        2.1 11   Ag     Ag     Ag     Ag     Ag       Ag
        2.1 11   Al     Al     Al     Al     Al       Al
        ...
        """  # nopep8
        section = data['section']
        label = data['label']

        if section not in self.data:
            self.data[section] = {}
        if label in self.data[section]:
            msg = "'{}' already defined in section '{}'".format(label, section)
            logger.error(msg)
            raise RuntimeError(msg)
        self.data[section][label] = data
        equivalences = self.data[section][label]['equivalences'] = {}

        for line in data['lines']:
            words = line.split()
            version, reference, atom_type, nonbond, bond, angle, \
                torsion, oop = words
            if atom_type not in equivalences:
                equivalences[atom_type] = {}
            if version in equivalences[atom_type]:
                msg = "atom type '{}' defined more than ".format(atom_type) + \
                      "once in section '{}'!".format(section)
                logger.error(msg)
                raise RuntimeError(msg)
            equivalences[atom_type][version] = {
                'reference': reference,
                'nonbond': nonbond,
                'bond': bond,
                'angle': angle,
                'torsion': torsion,
                'oop': oop
            }

    def _parse_biosym_auto_equivalence(self, data):
        """
        Process the atom type equivalences for automatic types

        #auto_equivalence     cff91_auto

        !                      Equivalences
        !       ------------------------------------------
        !Ver  Ref   Type  NonB Bond   Bond     Angle    Angle     Torsion   Torsion      OOP      OOP
        !                       Inct           End atom Apex atom End Atoms Center Atoms End Atom Center Atom
        !---- ---   ----  ---- ------ ----  ---------- --------- --------- -----------  -------- -----------
        2.0  1     Br    Br   Br     Br_   Br_        Br_       Br_       Br_          Br_      Br_
        2.0  1     Cl    Cl   Cl     Cl_   Cl_        Cl_       Cl_       Cl_          Cl_      Cl_
        ...
        """  # nopep8
        section = data['section']
        label = data['label']

        if section not in self.data:
            self.data[section] = {}
        if label in self.data[section]:
            msg = "'{}' already defined in section '{}'".format(label, section)
            logger.error(msg)
            raise RuntimeError(msg)
        self.data[section][label] = data
        equivalences = self.data[section][label]['equivalences'] = {}

        for line in data['lines']:
            words = line.split()
            version, reference, atom_type, nonbond, bond_increment, bond, \
                angle_end_atom, angle_center_atom, torsion_end_atom, \
                torsion_center_atom, oop_end_atom, oop_center_atom = words
            if atom_type not in equivalences:
                equivalences[atom_type] = {}
            if version in equivalences[atom_type]:
                msg = "atom type '{}' defined more than ".format(atom_type) + \
                      "once in section '{}'!".format(section)
                logger.error(msg)
                raise RuntimeError(msg)
            equivalences[atom_type][version] = {
                'reference': reference,
                'nonbond': nonbond,
                'bond_increment': bond_increment,
                'bond': bond,
                'angle_end_atom': angle_end_atom,
                'angle_center_atom': angle_center_atom,
                'torsion_end_atom': torsion_end_atom,
                'torsion_center_atom': torsion_center_atom,
                'oop_end_atom': oop_end_atom,
                'oop_center_atom': oop_center_atom
            }

    def _parse_biosym_bond_increments(self, data):
        """
        Process the bond increments

        #bond_increments      cff91_auto

        !Ver Ref    I     J     DeltaIJ   DeltaJI
        !--- ---  ----- -----   -------   -------
        2.1 11   Ag    Ag       0.0000   0.0000
        2.1 11   Al    Al       0.0000   0.0000
        ...
        """  # nopep8
        section = data['section']
        label = data['label']

        if section not in self.data:
            self.data[section] = {}
        if label in self.data[section]:
            msg = "'{}' already defined in section '{}'".format(label, section)
            logger.error(msg)
            raise RuntimeError(msg)
        self.data[section][label] = data
        parameters = data['parameters'] = {}
        data['constants'] = [
            ('deltaij', 'e'),
            ('deltaji', 'e')
        ]

        for line in data['lines']:
            words = line.split()
            version, reference, i, j, deltaij, deltaji = words
            # order canonically, i>j
            if i < j:
                i, j = j, i
                deltaij, deltaji = deltaji, deltaij
            key = (i, j)
            if key not in parameters:
                parameters[key] = {}
            if version in parameters[key]:
                msg = "bond increment '{}' '{}' defined more ".format(i, j) + \
                      "than once in section '{}'!".format(section)
                logger.error(msg)
                raise RuntimeError(msg)
            parameters[key][version] = {
                'reference': reference,
                'deltaij': deltaij,
                'deltaji': deltaji
            }

    def _parse_biosym_reference(self, data):
        """
        Process a 'reference' section, which looks like

        #reference 1
        @Author Biosym Technologies inc
        @Date 25-December-91
        cff91 forcefield created
        December 1991

        """
        section = data['section']
        label = data['label']

        if section not in self.data:
            self.data[section] = {}
        if label in self.data[section]:
            msg = "'{}' already defined in section '{}'".format(label, section)
            logger.error(msg)
            raise RuntimeError(msg)
        self.data[section][label] = data
        data['reference'] = data['lines']

        if not self.keep_lines:
            del data['lines']

    def make_canonical(self, symmetry, atom_types):
        """
        Using the symmetry, order the atom_types canonically
        """

        n = len(atom_types)
        flipped = False
        if n == 1:
            i = atom_types[0]
            return ((i), flipped)
        elif n == 2:
            i, j = atom_types
            if symmetry == 'like_bond':
                # order canonically, i<j
                if i > j:
                    i, j = j, i
                    flipped = True
                return ((i, j), flipped)
        elif n == 3:
            i, j, k = atom_types
            if symmetry == 'like_angle':
                # order canonically, i<k
                if i > k:
                    i, k = k, i
                    flipped = True
                return ((i, j, k), flipped)
        elif n == 4:
            i, j, k, l = atom_types
            if symmetry == 'like_torsion':
                # order canonically, j<k; i<l if j==k
                if j == k and i > l:
                    i, l = l, i
                    flipped = True
                elif j > k:
                    i, j, k, l = l, k, j, i
                    flipped = True
                return ((i, j, k, l), flipped)
            elif symmetry == 'like_oop':
                # j is central atom
                # order canonically, i<k<l; i=i<l or i<j=l
                i, k, l = sorted((i, k, l))
                flipped = [i, j, k, l] != atom_types
                return ((i, j, k, l), flipped)
            elif symmetry == 'like_angle-angle':
                # order canonically, i<l;
                if i > l:
                    i, l = l, i
                    flipped = True
                return ((i, j, k, l), flipped)

    def _parse_biosym_section(self, data):
        """
        Process the 1-term torsion parameters

        #torsion_1            cff91_auto

        > E = Kphi * [ 1 + cos(n*Phi - Phi0) ]

        !Ver Ref    I     J     K     L       KPhi     n     Phi0
        !--- ---  ----- ----- ----- -----   --------  ---  ---------
        2.0  2   *     c'_   c'_   *         0.4500    2   180.0000
        2.0  2   *     c'_   c=_   *         0.4500    2   180.0000
        ...
        """  # nopep8
        section = data['section']
        label = data['label']

        if section not in self.data:
            self.data[section] = {}
        if label in self.data[section]:
            msg = "'{}' already defined in section '{}'".format(label, section)
            logger.error(msg)
            raise RuntimeError(msg)

        self.data[section][label] = data

        # Copy in the metadata about this functional form
        data.update(metadata[section])

        parameters = data['parameters'] = {}

        for line in data['lines']:
            words = line.split()
            version, reference = words[0:2]
            symmetry = data['topology']['symmetry']
            n_atoms = data['topology']['n_atoms']
            key, flipped = self.make_canonical(symmetry, words[2:2+n_atoms])

            if key not in parameters:
                parameters[key] = {}
            if version in parameters[key]:
                msg = "value for '" + "' '".join(key) + " defined more " + \
                      "than once in section '{}'!".format(section)
                logger.error(msg)
                raise RuntimeError(msg)
            params = parameters[key][version] = {'reference': reference}
            values = words[2+n_atoms:]
            if 'fill' in data['topology']:
                n = data['topology']['fill']
                if n > 0:
                    if len(values) < 2*n:
                        values.extend(values[0:n])
            if flipped and 'flip' in data['topology']:
                n = data['topology']['flip']
                if n > 0:
                    first = values[0:n]
                    values = values[n:2*n]
                    values.extend(first)
            for constant, value in zip(data['constants'], values):
                params[constant[0]] = value

        if not self.keep_lines:
            del data['lines']
