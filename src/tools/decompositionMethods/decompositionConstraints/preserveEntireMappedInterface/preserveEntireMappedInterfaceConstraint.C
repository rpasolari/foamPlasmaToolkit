/*---------------------------------------------------------------------------*\
  File: preserveEntireMappedInterfaceConstraint.C
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::decompositionConstraints::preserveEntireMappedInterface

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "preserveEntireMappedInterfaceConstraint.H"
#include "addToRunTimeSelectionTable.H"
#include "syncTools.H"
#include "mappedPolyPatch.H"
#include "mappedWallPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace decompositionConstraints
{
    defineTypeName(preserveEntireMappedInterface);

    addToRunTimeSelectionTable
    (
        decompositionConstraint,
        preserveEntireMappedInterface,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decompositionConstraints::preserveEntireMappedInterface::preserveEntireMappedInterface
(
    const dictionary& dict
)
:
    decompositionConstraint(dict, typeName),
    patches_(coeffDict_.get<wordRes>("patches"))
{
    if (decompositionConstraint::debug)
    {
        Info<< type() << " : locking processor assignment for patches " 
            << patches_ << endl;
    }
}


Foam::decompositionConstraints::preserveEntireMappedInterface::preserveEntireMappedInterface
(
    const UList<wordRe>& patches
)
:
    decompositionConstraint(dictionary(), typeName),
    patches_(patches)
{
    if (decompositionConstraint::debug)
    {
        Info<< type() << " : locking processor assignment for patches " 
            << patches_ << endl;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::decompositionConstraints::preserveEntireMappedInterface::add
(
    const polyMesh& mesh,
    boolList& blockedFace,
    PtrList<labelList>& specifiedProcessorFaces,
    labelList& specifiedProcessor,
    List<labelPair>& explicitConnections
) const
{
    // BlockedFace size ensures all faces are unblocked for consideration
    // if they belong to the target patches.
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();
    const labelList patchIDs(patches_.matching(pbm.names()));

    label nUnblocked = 0;

    for (const label patchi : patchIDs)
    {
        const polyPatch& pp = pbm[patchi];

        forAll(pp, i)
        {
            label meshFacei = pp.start() + i;
            if (blockedFace[meshFacei])
            {
                blockedFace[meshFacei] = false;
                ++nUnblocked;
            }
        }
    }

    if (decompositionConstraint::debug & 2)
    {
        Info<< type() << " : unblocked "
            << returnReduce(nUnblocked, sumOp<label>()) << " patch faces" << endl;
    }

    syncTools::syncFaceList(mesh, blockedFace, andEqOp<bool>());
}


void Foam::decompositionConstraints::preserveEntireMappedInterface::apply
(
    const polyMesh& mesh,
    const boolList& blockedFace,
    const PtrList<labelList>& specifiedProcessorFaces,
    const labelList& specifiedProcessor,
    const List<labelPair>& explicitConnections,
    labelList& decomposition
) const
{
    Info << "Run the apply.. " << endl;
    // const polyBoundaryMesh& pbm = mesh.boundaryMesh();
    // const labelList patchIDs(patches_.matching(pbm.names()));

    // Info<< " [DEBUG] preserveEntireMappedInterface checking mesh: " << mesh.name() << endl;
    
    // if (patchIDs.empty())
    // {
    //     Info<< " [DEBUG] No matching patches found in mesh " << mesh.name() 
    //         << ". Check your decomposeParDict names!" << endl;
    // }

    // label nChanged = 0;
    // do
    // {
    //     nChanged = 0;
    //     for (const label patchi : patchIDs)
    //     {
    //         // Check for both generic mapped patches and mappedWall patches
    //         if (isType<mappedPolyPatch>(pbm[patchi]) || isType<mappedWallPolyPatch>(pbm[patchi]))
    //         {
    //             // Use polyPatch to get faceCells since it's common to both
    //             const polyPatch& pp = pbm[patchi];
    //             const labelList& faceCells = pp.faceCells();
                
    //             Info<< " [DEBUG] Found mapped patch: " << pp.name() 
    //                 << " with " << faceCells.size() << " faceCells." << endl;

    //             for (const label celli : faceCells)
    //             {
    //                 if (decomposition[celli] != Pstream::myProcNo())
    //                 {
    //                     decomposition[celli] = Pstream::myProcNo();
    //                     ++nChanged;
    //                 }
    //             }
    //         }
    //         else
    //         {
    //             Info<< " [DEBUG] Patch " << pbm[patchi].name() 
    //                 << " is NOT recognized as mapped. Type is: " << pbm[patchi].type() << endl;
    //         }
    //     }

    // } while (returnReduceOr(nChanged));
    
    // Info<< " [DEBUG] Total cells anchored on " << mesh.name() << ": " 
    //     << returnReduce(nChanged, sumOp<label>()) << endl;
}


// ************************************************************************* //
