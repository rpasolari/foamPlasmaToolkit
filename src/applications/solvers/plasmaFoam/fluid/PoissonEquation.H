/*---------------------------------------------------------------------------*\
  File: PoissonEquation.H
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Solves the Poisson equation for the electric potential in the plasma
      (fluid) regions:

          ∇·(ε∇φ) = -ρ_e

      where ε is the permittivity and ρ_e is the charge density. This
      equation couples the plasma charge distribution to the resulting
      electrostatic potential field.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

{
    fvScalarMatrix ePotentialEquation
    (
          fvm::laplacian(epsilon,ePotential)
        ==
          -eCharge
    );

    if(coupled)
    {
        Info << "TAKA" << endl;
        fvMatrixAssemblyPtr->addFvMatrix(ePotentialEquation);
        Info <<"ELAFKI" << endl;
    }
    else
    {
        ePotentialEquation.solve();
    }
}
