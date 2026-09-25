=======
History
=======
2026.9.25 -- Internal: declare rdkit directly
    * RDKit is imported directly but was only obtained through molsystem. It is now
      a declared requirement, so the package's dependencies are complete on their own.

2026.9.20.1 -- Warnings about the charges now reach the step's output
    * When the charges from the forcefield do not add up to the charge of the
      system, they are adjusted to compensate and a warning is issued. The warning
      only went to the log, which in a running job is whatever its output happens
      to be -- a queueing system's output file, or nothing at all -- so in practice
      it was rarely seen. Assigning a forcefield now returns its warnings as well,
      so that a step can print them in its own output alongside everything else it
      reports.
    * The warning said how much each atom's charge was adjusted by, rounded to
      three decimals, which showed a real adjustment of a hundred-thousandth of an
      electron as '0.000'. It now reports enough figures to be meaningful, as does
      the total charge it quotes.

2026.9.20 -- Bugfix: an atom type with no equivalence hid the real error
    * A forcefield's equivalence table does not have to cover every atom type: a
      combined forcefield such as 'oplsaa+' takes its equivalences from one of the
      forcefields it is built from, while another contributes types of its own that
      have none -- and need none, because it gives an explicit parameter for every
      term they appear in. Looking for a parameter for such a type reported, e.g.,
      ``KeyError: 'SBLRHMKNNHXPHG_803'``, which named the atom type and read as
      though the type itself were unknown. The type is now recognized as simply
      having no equivalences to fall back on, so the search carries on and, if
      nothing matches, says which term has no parameters.
    * Internal: reformatted a file for the current version of black.

2025.8.1 -- Bugfix: tabulated angles cause a crash
   * The changes for tracking the original parameters introduced an error into the code
     handling tabulated angles. This is now fixed.
   * Added the element symbol corresponding to the atom type in the energy
     expression. This allows LAMMPS to add the correct element symbol to trajectories,
     if desired.

2025.5.26 -- Track more of the original parameters
   * Added tracking of the original parameters for cross terms like bond-bond terms.
   * Improved chacking for duplicate parameters to remove most duplicates, except for
     class 2 forcefields, where the cross terms must have the same terms as the
     diagonal time, i.e. angle and bond-bond terms must share a list of angles.

2025.5.23 -- Added charges in templates and tracking original parameters
   * Added charges in templates, which override any assignment from the bond increments
     or charges for the atom types.
   * Track the original parameters before any unit conversions or other transformations,
     so that they can be printed along with the actual parameters used by e.g. LAMMPS.
   * Improved checking for duplicate parameters in the energy expression.

2025.4.7 -- Added ability to handle ReaxFF forcefields
   * Added code to handle ReaxFF forcefield.
   * Added seamm-reaxff utility to import forcefields in the standard Reax/LAMMPS format
     into SEAMM's format.
   * Added handling of a metadata section in the forcefields to be able to support the
     various types of forcefields.

2025.3.16 -- Added handling of Dreiding forcefield.

2025.1.21 -- Bugfix: torsions in 3-membered rings
   * The code allowed the torsion around a 3-membered ring which had the same atom at
     each end of the torsion. This is not a valid torsion, and the code now checks for
     it and removes it.

2024.6.27 -- Support for local forcefield files
   * Added support for local forcefields files which can either be used directly
     or included by existing files.
   * Added URI handler to support local files
   * Added support for BibTex references in forcefield files, and automatically adding
     citations to the Reference Handler.
   * Add 'fragments' section to forcefields for atom-typing via a fragment or entire
     molecule. This supports using LigParGen for OPLS-AA forcefields.

2023.8.27 -- Added support for tabulated angle potentials

2023.4.6 -- Added support for Buckingham potentials
   * Also improved unit handling across all terms in forcefields.

2023.3.5 -- Added molecule numbers for LAMMPS input
   * Added the molecule number for each atom for when using LAMMPS

2023.2.6 -- Added handling of OPLS-AA forcefield
   * Added handling of the OPLS-AA forcefield
   * Moved documentation to new MolSSI theme and diátaxis layout
   * Cleaned up internal dependencies and workflows for GitHub

2022.5.29 -- Fixed bug typing larger systems
   * Fixed bug with atom typing due to limit in matches. by @paulsaxe in #59

2022.2.3 -- Fixed bug due to changing ordering of atoms.
   * Fixed bug with atom type assignment due to changed order of atoms. In the process,
     switch to using RDKit directly, which is both more direct and avoids the ordering
     problem.

0.1.0 -- (2017-12-05)
   * First release on PyPI.
