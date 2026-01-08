/*---------------------------------------------------------------------------*\
  File: analyticMobility.C
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::analyticMobility.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "analyticMobility.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(analyticMobility, 0);
addToRunTimeSelectionTable(plasmaMobilityModel, analyticMobility, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

analyticMobility::analyticMobility
(
    const word& modelName,
    const dictionary& dict,
    const fvMesh& mesh,
    const plasmaSpecies& species,
    const label specieIndex
)
:
    plasmaMobilityModel(modelName, dict, mesh, species, specieIndex),
    muCoeff_(0.0),
    fieldExp_(0.0)
{
    if (!dict_.found("muCoeff"))
    {
        FatalIOErrorInFunction(dict_)
            << "Mobility model '" << modelName
            << "' requires entry 'muCoeff' in mobilityCoeffs dictionary."
            << exit(FatalIOError);
    }

    dict_.lookup("muCoeff") >> muCoeff_;

    if (!dict_.found("fieldExp"))
    {
        FatalIOErrorInFunction(dict_)
            << "Mobility model '" << modelName
            << "' requires entry 'fieldExp' in mobilityCoeffs dictionary."
            << exit(FatalIOError);
    }

    dict_.lookup("fieldExp") >> fieldExp_;
}


// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void analyticMobility::correct(volScalarField& mu) const
{
    const volVectorField& E = mesh_.lookupObject<volVectorField>("E");

    volScalarField dimensionlessE =
        mag(E) / dimensionedScalar("1", E.dimensions(), 1.0);

    mu = dimensionedScalar("mu", mu.dimensions(), muCoeff_)
       * pow(max(dimensionlessE, 1e-6), fieldExp_);

    mu.correctBoundaryConditions();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
