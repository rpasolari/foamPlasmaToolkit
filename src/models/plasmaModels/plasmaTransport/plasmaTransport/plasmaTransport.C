/*---------------------------------------------------------------------------*\
  File: plasmaTransport.C
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::plasmaTransport.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "ScharfetterGummel.H"
#include "plasmaTransport.H"
#include "plasmaTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(plasmaTransport, 0);

// * * * * * * * * * * * * * * Private Member Functions * * * * * * * * * *  //

void plasmaTransport::constructModels()
{
    // Loop over species and create a transport model for each one
    for (label i = 0; i < species_.nSpecies(); ++i)
    {
        const word& sName = species_.speciesNames()[i];
        const dictionary& sDict = species_.speciesDict(sName);

        if (!sDict.found("transportModel"))
        {
            FatalIOErrorInFunction(sDict)
                << "Species '" << sName
                << "' is missing required entry 'transportModel' in "
                << species_.dictName() << nl
                << exit(FatalIOError);
        }

        word modelName;
        sDict.lookup("transportModel") >> modelName;

        // Construct the model using the runtime selection system
        transportModels_.set
        (
            i,
            plasmaTransportModel::New(modelName, sDict, mesh_, species_, i, E_)
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

plasmaTransport::plasmaTransport
(
    plasmaSpecies& species,
    const fvMesh& mesh,
    const volVectorField& E
)
:
    mesh_(mesh),
    species_(species),
    E_(E),
    transportModels_(species.nSpecies()),
    flux_(),
    surfaceFlux_()
{
    constructModels();

    flux_.setSize(species.nSpecies());
    surfaceFlux_.setSize(species.nSpecies());
    for (label i = 0; i < species.nSpecies(); ++i)
    {
        flux_.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "Gamma_" + species.speciesNames()[i], 
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

        surfaceFlux_.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "phi_" + species.speciesNames()[i],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar
                (
                    "zero",
                    dimensionSet(0, 3, -1, 0, 0, 0, 0),
                    0.0)
            )
        );
    }
}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void plasmaTransport::correct()
{
    const label nSpecies = species_.nSpecies();

    for (label i = 0; i < nSpecies; ++i)
    {
        const word& sName = species_.speciesNames()[i];
        const dictionary& sDict = species_.speciesDict(i);

        transportModels_[i].correct();

        const volVectorField& U = transportModels_[i].driftVelocity();

        surfaceFlux_[i] = fvc::flux(U);

        volScalarField& n = species_.numberDensity(i);
        fvScalarMatrix nEqn(fvm::ddt(n));

        word modelName;
        sDict.lookup("transportModel") >> modelName;
        word advMode = sDict.lookupOrDefault<word>("advectionMode", "implicit");
        bool isExplicit = false;

        if (advMode == "explicit")
        {
            isExplicit = true;
        }
        else if (advMode != "implicit")
        {
            FatalIOErrorInFunction(sDict)
                << "Species '" << sName << "': Invalid advectionMode '" 
                << advMode << "'. Valid options: 'implicit', 'explicit'."
                << exit(FatalIOError);
        }

        word scheme = sDict.lookupOrDefault<word>
        (
            "driftDiffusionFluxScheme", 
            "standard"
        );

        if (modelName == "driftDiffusion")
        {
            const volScalarField* D = transportModels_[i].diffusivity();
            if (!D)
            {
                FatalErrorInFunction
                    << "DriftDiffusion model has null Diffusivity."
                    << exit(FatalError);
            }

            if (scheme == "ScharfetterGummel")
            {
                // SG Scheme: Forces Implicit
                if (isExplicit)
                {
                    WarningInFunction
                        << "Species '" << sName 
                        << "': advectionMode 'explicit' is ignored because "
                        << "Scharfetter-Gummel requires implicit coupling."
                        << endl;
                }
                nEqn += fvm::ScharfetterGummel(n, surfaceFlux_[i], *D);
            }
            else if (scheme == "standard")
            {
                // Standard Scheme
                if (isExplicit)
                {
                    nEqn -= fvc::div(surfaceFlux_[i], n);
                }
                else
                {
                    nEqn += fvm::div(surfaceFlux_[i], n);
                }
                // Diffusion is always implicit for DD
                nEqn -= fvm::laplacian(*D, n);
            }
            else
            {
                FatalIOErrorInFunction(sDict)
                    << "Species '" << sName 
                    << "': Invalid driftDiffusionFluxScheme '" << scheme 
                    << "'. Valid: 'standard', 'ScharfetterGummel'."
                    << exit(FatalIOError);
            }
        }
        else
        {
            // Check for invalid configuration
            if (scheme != "standard")
            {
                FatalIOErrorInFunction(sDict)
                    << "Species '" << sName 
                    << "': 'driftDiffusionFluxScheme' cannot be set to '" 
                    << scheme << "' when using transportModel '" 
                    << modelName << "'. Only 'standard' is allowed." 
                    << exit(FatalIOError);
            }

            if (isExplicit)
            {
                nEqn -= fvc::div(surfaceFlux_[i], n);
            }
            else
            {
                nEqn += fvm::div(surfaceFlux_[i], n);
            }
        }

        nEqn.solve();
    }
}

tmp<volScalarField> plasmaTransport::elecConductionCoeff() const
{
    tmp<volScalarField> elecConductionCoeff =
        tmp<volScalarField>::New
        (
            volScalarField
            (
                IOobject
                (
                    "elecConductionCoeff",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar
                (
                    "zeroConduction",
                    dimensionSet(-1,-3,3,0,0,2,0), 
                    0.0
                )
            )
        );

    const label nSpecies = species_.nSpecies();
    for (label i = 0; i < nSpecies; ++i)
    {
        elecConductionCoeff.ref() += transportModels_[i].elecConductionCoeff();
    }

    return elecConductionCoeff;
}

tmp<volScalarField> plasmaTransport::elecDiffusionCharge() const
{
    tmp<volScalarField> elecDiffusionCharge =
        tmp<volScalarField>::New
        (
            volScalarField
            (
                IOobject
                (
                    "elecDiffusionCharge",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar
                (
                    "zeroDiffusionCharge",
                    dimensionSet(0,-3,0,0,0,1,0), 
                    0.0
                )
            )
        );

    const label nSpecies = species_.nSpecies();
    for (label i = 0; i < nSpecies; ++i)
    {
        elecDiffusionCharge.ref() += transportModels_[i].elecDiffusionCharge();
    }

    return elecDiffusionCharge;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //












