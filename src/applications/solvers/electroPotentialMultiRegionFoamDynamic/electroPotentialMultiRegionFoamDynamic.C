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
    electroPotentialMultiRegionFoam

Description
    Steady-state multi-region electrostatic potential solver for coupled
    fluid and dielectric domains.

Usage
    \b electroPotentialMultiRegionFoam [OPTIONS]

    Example:
        electroPotentialMultiRegionFoam -case [CASE]

Author
    Rention Pasolari
    Contact: r.pasolari@gmail.com
\*---------------------------------------------------------------------------*/

#include "foamPlasmaToolkitConstants.H"
#include "dynamicFvMesh.H"
#include "fvCFD.H"
#include "regionProperties.H"
#include "coordinateSystem.H"
#include "loopControl.H"

#include "adaptiveFvMesh.H"
#include "mappedPolyPatch.H"
#include "fvSolution.H"
#include "solutionControl.H"

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Steady-state multi-region electrostatic potential solver for coupled"
        " fluid and dielectric domains."
    );

    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMeshes.H"
    #include "createFields.H"

    #include "createCoupledRegions.H"

    #include "readElectricPotentialControls.H"

    #include "AMRConfiguration.H"

    Info<< "\nStarting iteration loop\n" << endl;

// In Region A (Fluid)
volScalarField refineFieldA
(
    IOobject("refineField", runTime.timeName(), fluidRegions[0], IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
    fluidRegions[0],
    dimensionedScalar(dimless, 0.0)
);

// In Region B (Dielectric)
volScalarField refineFieldB
(
    IOobject("refineField", runTime.timeName(), dielectricRegions[0], IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
    dielectricRegions[0],
    dimensionedScalar(dimless, 0.0)
);


    while (runTime.run())
    {
        ++runTime;

// 1. Cast the regions to adaptiveFvMesh
adaptiveFvMesh& amrMeshA = refCast<adaptiveFvMesh>(fluidRegions[0]);
adaptiveFvMesh& amrMeshB = refCast<adaptiveFvMesh>(dielectricRegions[0]);

// 1. Calculate the Error in the Master (Region A)
amrMeshA.amr().error().update();
refineFieldA.primitiveFieldRef() = amrMeshA.amr().error().error();

// 2. Prepare Region B for a Fresh Signal
// Resetting to 0 prevents the "double refinement" at x=0.8

// refineFieldB.primitiveFieldRef() = 0.0;

// 3. Synchronize Boundaries (Mapping)
// This pushes 1.0 from A into the internal cells of B
refineFieldA.correctBoundaryConditions();
refineFieldB.correctBoundaryConditions();
// 2. Force the mapped patches to reconstruct their maps

// 4. Mesh Topology Update
// Region A updates its mesh
amrMeshA.update();

// Region B updates its mesh using the fresh 1.0 signals
amrMeshB.update();

// 5. Cleanup
// Clear fields so they don't carry over to the next physical time step
// refineFieldA.primitiveFieldRef() = 0.0;
// refineFieldB.primitiveFieldRef() = 0.0;

// --- 4. Re-sync after topology change (Optional but safer) ---
// refineFieldA.instance() = runTime.timeName();
// refineFieldB.instance() = runTime.timeName();

forAll(amrMeshA.boundaryMesh(), patchi)
{
    const polyPatch& pp = amrMeshA.boundaryMesh()[patchi];
    const mappedPatchBase* mpPtr = dynamic_cast<const mappedPatchBase*>(&pp);

    if (mpPtr)
    {
        // Use const_cast to allow calling the non-const clearOut()
        const_cast<mappedPatchBase*>(mpPtr)->clearOut(); 
    }
}
forAll(amrMeshB.boundaryMesh(), patchi)
{
    const polyPatch& pp = amrMeshB.boundaryMesh()[patchi];
    const mappedPatchBase* mpPtr = dynamic_cast<const mappedPatchBase*>(&pp);

    if (mpPtr)
    {
        const_cast<mappedPatchBase*>(mpPtr)->clearOut();
    }
}


if (coupled)
{
// if (topologyChanged && fvMatrixAssemblyPtr)
{
    Info << "Forcing mesh to forget old matrix assembly..." << endl;

    // This reaches into the mesh database and removes the cached addressing
    // 'lduAssembly4' is the default name OpenFOAM uses internally
    ePotentialFluid[0].mesh().thisDb().checkOut("lduAssembly4");

    fvMatrixAssemblyPtr.reset(nullptr);
    fvMatrixAssemblyPtr.reset
    (
        new fvMatrix<scalar>
        (
            ePotentialFluid[0],
            dimensionSet(0,0,1,0,0,1,0)
        )
    );
}
}



 


        Info << "Time = " << runTime.timeName() << nl << endl;
        
        // Solve the Poisson/Laplace Equation (electric potential)
        if(coupled)
        {
            #include "solveElectricPotentialCoupled.H"
        }
        else
        {
            #include "solveElectricPotentialNonCoupled.H"
        }

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}







//    {
//     // Array of region pairs to iterate through both directions
//     // Pair: {SourceMesh, TargetMesh, SourcePatchName}
//     struct DebugPair { const fvMesh& src; const fvMesh& tgt; word patchName; word targetPatchName; };
    
//     DebugPair pairs[2] = {
//         {amrMeshA, amrMeshB, "regionA_to_regionB", "regionB_to_regionA"},
//         {amrMeshB, amrMeshA, "regionB_to_regionA", "regionA_to_regionB"}
//     };

//     for (int pI=0; pI<2; ++pI)
//     {
//         const fvMesh& srcMesh = pairs[pI].src;
//         const fvMesh& tgtMesh = pairs[pI].tgt;
        
//         label patchID = srcMesh.boundaryMesh().findPatchID(pairs[pI].patchName);

//         if (patchID != -1)
//         {
//             const polyPatch& ppSrc = srcMesh.boundaryMesh()[patchID];
//             Info << "\n--- Debugging: [" << srcMesh.name() << "] -> [" << tgtMesh.name() << "] ---" << endl;
//             Info << "Patch: " << ppSrc.name() << " | Faces: " << ppSrc.size() << endl;

//             if (isA<mappedPatchBase>(ppSrc))
//             {
//                 const mappedPatchBase& mpp = refCast<const mappedPatchBase>(ppSrc);
//                 const mapDistribute& map = mpp.map();
                
//                 if (map.subMap().size() > 0)
//                 {
//                     const labelList& faceMap = map.subMap()[0];
//                     label tgtPatchID = tgtMesh.boundaryMesh().findPatchID(mpp.samplePatch());

//                     if (tgtPatchID != -1)
//                     {
//                         const polyPatch& ppTgt = tgtMesh.boundaryMesh()[tgtPatchID];
//                         const labelList& srcOwn = srcMesh.faceOwner();
//                         const labelList& tgtOwn = tgtMesh.faceOwner();

//                         forAll(ppSrc, facei)
//                         {
//                             label localFaceTgt = faceMap[facei];

//                             if (localFaceTgt >= 0)
//                             {
//                                 label globalFaceSrc = ppSrc.start() + facei;
//                                 label globalFaceTgt = ppTgt.start() + localFaceTgt;
                                
//                                 label cellSrc = srcOwn[globalFaceSrc];
//                                 label cellTgt = tgtOwn[globalFaceTgt];

//                                 Info<< "  SrcLocFace: " << facei 
//                                     << " (GlobFace: " << globalFaceSrc << ", Cell: " << cellSrc << ")"
//                                     << "\n    -> TgtLocFace: " << localFaceTgt 
//                                     << " (GlobFace: " << globalFaceTgt << ", Cell: " << cellTgt << ")"
//                                     << "\n    CoordSrc: " << srcMesh.faceCentres()[globalFaceSrc]
//                                     << " | CoordTgt: " << tgtMesh.faceCentres()[globalFaceTgt]
//                                     << " | Dist: " << mag(srcMesh.faceCentres()[globalFaceSrc] - tgtMesh.faceCentres()[globalFaceTgt])
//                                     << nl << endl;
//                             }
//                             else
//                             {
//                                 Info<< "  SrcLocFace: " << facei << " | *** UNMAPPED ***" << endl;
//                             }
//                         }
//                     }
//                 }
//                 else { Info << "  Map not yet initialized by solver." << endl; }
//             }
//         }
//     }
//     Info << "--- End Full Debugging ---\n" << endl;
// }
