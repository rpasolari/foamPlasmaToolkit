/*---------------------------------------------------------------------------*\
  File: analyticDiffusivity.C
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::analyticDiffusivity.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "analyticDiffusivity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(analyticDiffusivity, 0);
addToRunTimeSelectionTable
(
    plasmaDiffusivityModel,
    analyticDiffusivity,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

analyticDiffusivity::analyticDiffusivity
(
    const word& modelName,
    const dictionary& dict,
    const fvMesh& mesh,
    const plasmaSpecies& species,
    const label specieIndex
)
:
    plasmaDiffusivityModel(modelName, dict, mesh, species, specieIndex),
    DCoeff_(0.0),
    fieldExp_(0.0)
{
    if (!dict_.found("DCoeff"))
    {
        FatalIOErrorInFunction(dict_)
            << "Diffusivity model '" << modelName
            << "' requires entry 'DCoeff' in diffusivityCoeffs dictionary."
            << exit(FatalIOError);
    }

    dict_.lookup("DCoeff") >> DCoeff_;

    if (!dict_.found("fieldExp"))
    {
        FatalIOErrorInFunction(dict_)
            << "Diffusivity model '" << modelName
            << "' requires entry 'fieldExp' in diffusivityCoeffs dictionary."
            << exit(FatalIOError);
    }

    dict_.lookup("fieldExp") >> fieldExp_;
}


// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void analyticDiffusivity::correct(volScalarField& D) const
{
    const volVectorField& E = mesh_.lookupObject<volVectorField>("E");

    volScalarField dimensionlessE =
        mag(E) / dimensionedScalar("1", E.dimensions(), 1.0);

    D = dimensionedScalar("D", D.dimensions(), DCoeff_)
       * pow(max(dimensionlessE, 1e-6), fieldExp_);

    D.correctBoundaryConditions();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
