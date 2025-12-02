/*---------------------------------------------------------------------------*\
  File: driftDiffusion.C
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::driftDiffusion.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "fvc.H"

#include "driftDiffusion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(driftDiffusion, 0);
addToRunTimeSelectionTable(plasmaTransportModel, driftDiffusion, dictionary);

// * * * * * * * * * * * * * * Private Member Functions * * * * * * * * * *  //

void driftDiffusion::constructModels()
{
    // Construct mobility model

    if (!dict_.found("mobilityModel"))
    {
        FatalIOErrorInFunction(dict_)
            << "driftDiffusion requires a 'mobilityModel' entry\n"
            << "Please specify mobilityModel <modelName>;\n"
            << exit(FatalIOError);
    }

    word mobilityModelName;
    dict_.lookup("mobilityModel") >> mobilityModelName;

    if (!dict_.found("mobilityCoeffs"))
    {
        FatalIOErrorInFunction(dict_)
            << "Missing sub-dictionary 'mobilityCoeffs' for mobilityModel '"
            << mobilityModelName << "'\n"
            << "driftDiffusion requires mobilityCoeffs { ... }\n"
            << exit(FatalIOError);
    }

    const dictionary& mobilityCoeffs = dict_.subDict("mobilityCoeffs");

    mobilityModel_.reset
    (
        plasmaMobilityModel::New
        (
            mobilityModelName,
            mobilityCoeffs,
            mesh_,
            species_,
            specieIndex_
        )
    );

    // Construct diffusivity model

    if (!dict_.found("diffusivityModel"))
    {
        FatalIOErrorInFunction(dict_)
            << "driftDiffusion requires a 'diffusivityModel' entry\n"
            << "Please specify diffusivityModel <modelName>;\n"
            << exit(FatalIOError);
    }

    word diffusivityModelName;
    dict_.lookup("diffusivityModel") >> diffusivityModelName;

    if (!dict_.found("diffusivityCoeffs"))
    {
        FatalIOErrorInFunction(dict_)
            << "Missing sub-dictionary 'diffusivityCoeffs' for "
            << "diffusivityModel '" << diffusivityModelName << "'\n"
            << "driftDiffusion requires diffusivityCoeffs { ... }\n"
            << exit(FatalIOError);
    }

    const dictionary& diffusivityCoeffs = dict_.subDict("diffusivityCoeffs");

    diffusivityModel_.reset
    (
        plasmaDiffusivityModel::New
        (
            diffusivityModelName,
            diffusivityCoeffs,
            mesh_,
            species_,
            specieIndex_
        )
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

driftDiffusion::driftDiffusion
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
            "driftVelocity" + species.speciesNames()[specieIndex],
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
    ),

    mobility_
    (
        IOobject
        (
            "mu_" + species.speciesNames()[specieIndex],
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(-1, 0, 2, 0, 0, 1, 0), 0.0)
    ),
    
    diffusivity_
    (
        IOobject
        (
            "D_" + species.speciesNames()[specieIndex],
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0, 0, 0), 0.0)
    )
{
    constructModels();
}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void driftDiffusion::correct()
{
    mobilityModel_->correct(mobility_);
    diffusivityModel_->correct(diffusivity_);

    scalar qSign = (species_.speciesChargeNumber(specieIndex_) > 0)? 1.0 : -1.0;

    driftVelocity_ = qSign * mobility_ * E_;
}

const volVectorField& driftDiffusion::driftVelocity() const
{
    return driftVelocity_;
}

const volScalarField* driftDiffusion::diffusivity() const
{
    return &diffusivity_;
}

tmp<volScalarField> driftDiffusion::elecConductionCoeff() const
{
    return 
        mag(species_.speciesCharge(specieIndex_))
      * mobility_
      * species_.numberDensity(specieIndex_);

}

tmp<volScalarField> driftDiffusion::elecDiffusionCharge() const
{
    return species_.speciesCharge(specieIndex_)
           * fvc::laplacian(diffusivity_, species_.numberDensity(specieIndex_));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
