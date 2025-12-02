/*---------------------------------------------------------------------------*\
  File: noneTransportModel.C
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::noneTransportModel.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "noneTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(noneTransportModel, 0);
addToRunTimeSelectionTable(plasmaTransportModel,noneTransportModel,dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

noneTransportModel::noneTransportModel
(
    const word& modelName,
    const dictionary& dict,
    const fvMesh& mesh,
    const plasmaSpecies& species,
    const label specieIndex,
    const volVectorField& E
)
:
    plasmaTransportModel(modelName, dict, mesh, species, specieIndex, E),
    driftVelocity_(
        IOobject(
            "driftVelocity",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "driftVelocity",
            dimensionSet(0, 1, -1, 0, 0, 0, 0),
            vector::zero
        )
    )
{}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void noneTransportModel::correct()
{
    Info << "Correct is called in noneTransportModel!!" << endl;
    //Do nothing here
}

const volVectorField& noneTransportModel::driftVelocity() const
{
    return driftVelocity_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
