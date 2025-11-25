/*---------------------------------------------------------------------------*\
  File: PoissonEquation.H
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Solves the Poisson equation for the electric potential in the plasma
      (gas) region:

          ∇·(ε∇φ) = -ρ_e

      where ε is the permittivity and ρ_e is the charge density. This
      equation couples the plasma charge distribution to the resulting
      electrostatic potential field.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

if (!coupled)
{
    Info<< "\nSolving for ePotential in region " 
        << gasRegionName << endl;

    for (int nonOrth = 0; nonOrth <= nNonOrthoCorr; ++nonOrth)
    {
        fvScalarMatrix PoissonEquation
        (
            fvm::laplacian(epsilon,ePotential)
            ==
            -chargeDensity
        );

        PoissonEquation.relax();

        fvOptions.constrain(PoissonEquation);

        auto sp = PoissonEquation.solve();

        if (nonOrth == 0)
        {
            ePotentialNonCoupledResidual = sp.initialResidual();
        }
    }

    fvOptions.correct(ePotential);
}
else
{
        fvScalarMatrix PoissonEquation
        (
            fvm::laplacian(epsilon,ePotential)
            ==
            -chargeDensity
        );

        PoissonEquation.relax();

        fvOptions.constrain(PoissonEquation);

        fvMatrixAssemblyPtr->addFvMatrix(PoissonEquation);
}
