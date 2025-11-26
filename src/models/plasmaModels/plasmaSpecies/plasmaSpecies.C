/*---------------------------------------------------------------------------*\
  File: plasmaSpecies.C
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::plasmaSpecies.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

//TODO: ADD ne for return of electrons density, nn for neutrals etc

#include "plasmaSpecies.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(plasmaSpecies, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

plasmaSpecies::plasmaSpecies(const fvMesh& mesh)
:
    IOdictionary
    (
        IOobject
        (
            "plasmaSpecies",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    speciesNames_(),
    nSpecies_(0),
    speciesCharge_(),
    numberDensity_(),
    flux_()
{
    if (!found("species"))
    {
        FatalIOErrorInFunction(*this)
            << "Required entry 'species' is missing in dictionary "
            << objectPath() << nl
            << exit(FatalIOError);
    }

    // Read species list
    lookup("species") >> speciesNames_;
    nSpecies_ = speciesNames_.size();

    speciesCharge_.setSize(nSpecies_);
    numberDensity_.setSize(nSpecies_);
    flux_.setSize(nSpecies_);

    // Read species properties dictionary
    if (!found("speciesProperties"))
    {
        FatalIOErrorInFunction(*this)
            << "Missing required dictionary 'speciesProperties' in "
            << objectPath() << nl << exit(FatalIOError);
    }

    const dictionary& propsDict = subDict("speciesProperties");

    // Read defaultProperties if present
    if (propsDict.found("defaultProperties"))
    {
        defaultDict_ = propsDict.subDict("defaultProperties");
    }
    else
    {
        defaultDict_.clear();
    }

    // Prepare storage for merged per-species dicts
    speciesDicts_.clear();
    speciesDicts_.resize(nSpecies_);

    // Loop over all species
    for (label i = 0; i < nSpecies_; ++i)
    {
        const word& sName = speciesNames_[i];

        // Each species must have a dictionary under speciesProperties
        if (!propsDict.found(sName))
        {
            FatalIOErrorInFunction(*this)
                << "Species '" << sName << "' is listed in 'species' but "
                << "has no sub-dictionary in " << objectPath() << nl
                << exit(FatalIOError);
        }

        // Build merged properties (defaults + overrides)
        dictionary mergedDict(defaultDict_);
        mergedDict.merge(propsDict.subDict(sName));
        speciesDicts_.insert(sName, mergedDict);

        // Read charge number (required)
        if (!mergedDict.found("charge"))
        {
            FatalIOErrorInFunction(*this)
                << "Species '" << sName << "' is missing required entry "
                << "'charge' in " << objectPath() << nl
                << exit(FatalIOError);
        }

        speciesCharge_[i] = readScalar(mergedDict.lookup("charge"));

        numberDensity_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "n_" + speciesNames_[i],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );

        flux_.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "Gamma_" + speciesNames_[i], 
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedVector
                (
                    "zero",
                    dimensionSet(0, -2, -1, 0, 0, 0, 0),
                    vector::zero
                )
            )
        );
    }

}  

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
