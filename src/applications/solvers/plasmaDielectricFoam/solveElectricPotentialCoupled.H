/*---------------------------------------------------------------------------*\
  File: solveElectricPotentialCoupled.H
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Solves the Poisson/Laplace equations for the electric potential in all
      regions (fluid + dielectric), in a coupled way (monolithically). 

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

forAll(fluidRegions, i)
{
    fvMesh& mesh = fluidRegions[i];
    #include "setRegionFluidFields.H"
    #include "PoissonEquation.H"
}

forAll(dielectricRegions, i)
{
    fvMesh& mesh = dielectricRegions[i];
    #include "setRegionDielectricFields.H"
    #include "LaplaceEquation.H"
}

Info<< "\nSolving for ePotential in coupled regions (monolithically)\n" << endl;

fvMatrixAssemblyPtr->solve();

forAll(fluidRegions, i)
{
    volScalarField& ePotential = ePotentialFluid[i];
    fv::options& fvOptions = fvOptionsFluid[i];

    fvOptions.correct(ePotential);
    ePotential.correctBoundaryConditions();
}

forAll(dielectricRegions, i)
{
    volScalarField& ePotential = ePotentialDielectric[i];

    ePotential.correctBoundaryConditions();
}

Info<< "\n>>> Solved for ePotential in coupled regions" << endl;
fvMatrixAssemblyPtr->clear();





