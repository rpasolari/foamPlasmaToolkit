/*---------------------------------------------------------------------------*\
  File: solveFluid.H
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Derived from solveFluid.H in OpenFOAM v2412.
      Solves the governing equations for the fluid regions.
      Executed for each fluid region within the main time loop.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

Info << "\nSolving for fluid region "
     << fluidRegions[i].name() << endl;

