/*---------------------------------------------------------------------------*\
  File: setRegionFluidFields.H
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Derived from setRegionFluidFields.H in OpenFOAM v2412.
      Defines convenient field references for each fluid region.
      Associates region-specific variables with their corresponding
      data structures created in createFluidFields.H for use within plasmaFoam.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

volScalarField& ePotential = ePotentialFluid[i];
volScalarField& eCharge = eChargeFluid[i];

dimensionedScalar& epsilon = epsilonFluid[i];

fv::options& fvOptions = fvOptionsFluid[i];


