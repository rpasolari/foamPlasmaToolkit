/*---------------------------------------------------------------------------*\
  File: solveDielectric.H
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Solves the governing equations for the dielectric regions in plasmaFoam.
      Executed for each dielectric region within the main time loop.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

Info << "\nSolving for dielectric region " 
     << dielectricRegions[i].name() << endl;

#include "LaplaceEquation.H"
