/*---------------------------------------------------------------------------*\
License
    This file is part of the foamPlasmaToolkit.

    The foamPlasmaToolkit is not part of OpenFOAM but is developed using the
    OpenFOAM framework and linked against OpenFOAM libraries.

    Copyright (C) 2025 Rention Pasolari

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

Application
    plasmaDielectricFoam

Description
    Transient solver for coupled gas (plasma) and dielectric domains developed
    for plasma simulation purposes.

Usage
    \b plasmaDielectricFoam [OPTIONS]

    Example:
        plasmaDielectricFoam -case plasmaDielectricTest

Author
    Rention Pasolari
    Contact: r.pasolari@gmail.com
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "regionProperties.H"
#include "fvOptions.H"
#include "coordinateSystem.H"
#include "loopControl.H"
#include "fvSolution.H"
#include "solutionControl.H"
#include "mappedPatchBase.H"
#include "solidSurfaceFluxFvPatchScalarField.H"

#include "foamPlasmaToolkitConstants.H"
#include "plasmaSpecies.H"
#include "plasmaTransport.H"
#include "plasmaTimeControl.H"

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for coupled gas (plasma) and dielectric domains"
        " developed for plasma simulation purposes."
    );

    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMeshes.H"
    #include "createFields.H"

    plasmaTimeControl timeControl(runTime, gasMesh());
    timeControl.setInitialDeltaT(transport);

    #include "createCoupledRegions.H"

    #include "readElectricPotentialControls.H"

    #include "reportSimulationSummary.H"
    #include "updateChargeDensity.H"

    Info<< "\nStarting iteration loop\n" << endl;

    while (runTime.run())
    {
        ++runTime;

        Info << "Time = " << runTime.timeName() << nl << endl;

        gasMesh().update();

        timeControl.adjustDeltaT(transport);

        // Solve the Poisson/Laplace Equation (electric potential)
        if(coupled)
        {
            #include "solveElectricPotentialCoupled.H"
        }
        else
        {
            #include "solveElectricPotentialNonCoupled.H"
        }

        // Update the Electric field
        #include "calculateElectricField.H"

        transport.correct();
        // Info << "ELAFKI3" << endl;
        #include "updateChargeDensity.H"
        #include "updateSurfCharge.H"

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}
