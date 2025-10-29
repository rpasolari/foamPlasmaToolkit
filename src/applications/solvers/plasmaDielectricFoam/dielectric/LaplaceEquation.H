/*---------------------------------------------------------------------------*\
  File: LaplaceEquation.H
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Solves the Laplace equation for the electric potential in the dielectric
      regions:

          ∇·(ε∇φ) = 0

      where ε is the permittivity. This represents electrostatic equilibrium
      in a charge-free region (ρ_e = 0). 

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

{
    fvScalarMatrix LaplaceEquation
    (
        fvm::laplacian(epsilon,ePotential)
    );

    LaplaceEquation.relax();

    if(coupled)
    {
        fvMatrixAssemblyPtr->addFvMatrix(LaplaceEquation);
    }
    else
    {
        Info<< "\nSolving for ePotential in region " 
            << dielectricRegions[i].name() << endl;

        for (int nonOrth = 0; nonOrth <= nNonOrthoCorrDielectric[i]; ++nonOrth)
        {
            auto sp = LaplaceEquation.solve();

            if (nonOrth == 0)
            {
                ePotentialNonCoupledResidualDielectric[i] =
                    sp.initialResidual();
            }

            if (nonOrth == nNonOrthoCorrDielectric[i])
            {
                ePotential.correctBoundaryConditions();
            }
        }
    }    
}
