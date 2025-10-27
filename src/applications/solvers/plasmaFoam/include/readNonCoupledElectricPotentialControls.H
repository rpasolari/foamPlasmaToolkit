/*---------------------------------------------------------------------------*\
  File: readNonCoupledElectricPotentialControls.H
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

Description:
    Reads the convergence tolerances for the uncoupled (segregated) electric
    potential solution. These values are used to perform an outer iteration
    loop over all regions and to stop the iterations as soon as *all* regions
    reach the specified residual tolerance.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

List<scalar> ePotentialNonCoupledToleranceFluid
(
    fluidRegions.size(),
    Foam::GREAT
);

List<scalar> ePotentialNonCoupledToleranceDielectric
(
    dielectricRegions.size(),
    Foam::GREAT
);

scalar ePotentialNonCoupledLoopMaxIter = 0;

bool foundAnyRegionMaxIter = false;

// --- FLUID REGIONS ---
forAll(fluidRegions, i)
{
    const fvMesh& mesh = fluidRegions[i];
    const dictionary rc =
        mesh.solutionDict().subOrEmptyDict
        (
            "electricPotentialNonCoupledResidualConvergence"
        );

    ePotentialNonCoupledToleranceFluid[i] =
        rc.getOrDefault<scalar>("tolerance", 1.0e-16);

    if (rc.found("maxIter"))
    {
        foundAnyRegionMaxIter = true;
        scalar regionMaxIter = rc.get<scalar>("maxIter");
        ePotentialNonCoupledLoopMaxIter =
            max(ePotentialNonCoupledLoopMaxIter, regionMaxIter);
    }
}

// --- DIELECTRIC REGIONS ---
forAll(dielectricRegions, i)
{
    const fvMesh& mesh = dielectricRegions[i];
    const dictionary rc =
        mesh.solutionDict().subOrEmptyDict
        (
            "electricPotentialNonCoupledResidualConvergence"
        );

    ePotentialNonCoupledToleranceDielectric[i] =
        rc.getOrDefault<scalar>("tolerance", 1.0e-16);

    if (rc.found("maxIter"))
    {
        foundAnyRegionMaxIter = true;
        scalar regionMaxIter = rc.get<scalar>("maxIter");
        ePotentialNonCoupledLoopMaxIter =
            max(ePotentialNonCoupledLoopMaxIter, regionMaxIter);
    }
}

if (!foundAnyRegionMaxIter)
{
    ePotentialNonCoupledLoopMaxIter = 100;
}
