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
    const plasmaSpecies& species
)
:
    plasmaTransportModel(modelName, dict, mesh, species)
{}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void noneTransportModel::correct()
{
    //Do nothing here
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
