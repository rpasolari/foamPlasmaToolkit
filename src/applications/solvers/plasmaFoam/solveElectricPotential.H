/*---------------------------------------------------------------------------*\
  File: solveElectricPotential.H
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Solves the Poisson/Laplace equations for the electric potential in all
      regions (fluid + dielectric). 

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

if(!coupled)
{
    for (int outer = 1; outer <= ePotentialNonCoupledLoopMaxIter; ++outer)
    {
        Info<< "\nOuter electric solve iteration = " << outer << nl;

        // Reset residual lists each iteration
        forAll(ePotentialNonCoupledResidualFluid, i)
            ePotentialNonCoupledResidualFluid[i] = Foam::GREAT;

        forAll(ePotentialNonCoupledResidualDielectric, i)
            ePotentialNonCoupledResidualDielectric[i] = Foam::GREAT;

        // --- FLUID REGIONS ---
        forAll(fluidRegions, i)
        {
            fvMesh& mesh = fluidRegions[i];
            #include "setRegionFluidFields.H"
            #include "PoissonEquation.H"
        }

        // --- DIELECTRIC REGIONS ---
        forAll(dielectricRegions, i)
        {
            fvMesh& mesh = dielectricRegions[i];
            #include "setRegionDielectricFields.H"
            #include "LaplaceEquation.H"
        }

        // --- CHECK CONVERGENCE ---
        bool allOK = true;

        forAll(ePotentialNonCoupledResidualFluid, i)
        {
            if
            (
                ePotentialNonCoupledResidualFluid[i]
              >= ePotentialNonCoupledToleranceFluid[i]
            )
            {
                allOK = false;
            }
        }

        forAll(ePotentialNonCoupledResidualDielectric, i)
        {
            if
            (
                ePotentialNonCoupledResidualDielectric[i]
              >= ePotentialNonCoupledToleranceDielectric[i]
            )
            {
                allOK = false;
            }
        }

        if (allOK)
        {
            Info<< ">> ePotential converged in "
                << outer << " outer iterations." << nl;
            break;
        }
    }
}
else
{
    forAll(fluidRegions, i)
    {
        fvMesh& mesh = fluidRegions[i];
        #include "setRegionFluidFields.H"
        #include "PoissonEquation.H"
    }

    // --- DIELECTRIC REGIONS ---
    forAll(dielectricRegions, i)
    {
        fvMesh& mesh = dielectricRegions[i];
        #include "setRegionDielectricFields.H"
        #include "LaplaceEquation.H"
    }

    Info<< "\nSolving for electric potential in coupled regions" 
        << endl;

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

    Info<< "\nSolved for electric potential in coupled regions" << endl;
    fvMatrixAssemblyPtr->clear();
}




