/*---------------------------------------------------------------------------*\
  File: createGasFields.H
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Derived from createFluidFields.H in OpenFOAM v2412.
      Defines and initializes the gas-region fields.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

Info << "*** Reading fields for region "
     << gasRegionName << nl << endl;

Info << "   Adding to ePotential\n" << endl;
volScalarField ePotential
(
    IOobject
    (
        "ePotential",
        runTime.timeName(),
        gasMesh(),
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    gasMesh()
);

Info << "   Adding to E\n" << endl;
volVectorField E
(
    IOobject
    (
        "E",
        runTime.timeName(),
        gasMesh(),
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    gasMesh(),
    dimensionedVector("zero", dimensionSet(1, 1, -3, 0, 0, -1, 0), vector::zero)
);

autoPtr<surfaceScalarField> surfCharge;

{
    IOobject surfChargeHeader
    (
        "surfCharge",
        runTime.timeName(),
        gasMesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    );

    if (surfChargeHeader.typeHeaderOk<surfaceScalarField>())
    {
        surfCharge.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "surfCharge",
                    runTime.timeName(),
                    gasMesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                gasMesh()
            )
        );

        regIOobject::store(surfCharge.ptr());
        Info << "   Adding to surfCharge\n" << endl;
    }
}

IOdictionary electricPropertiesDict
(
    IOobject
    (
        "electricProperties",
        runTime.constant(),
        gasMesh(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

scalar epsilonR =
    electricPropertiesDict.lookupOrDefault<scalar>
    (
        "dielectricConstant",
        1.0
    );

Info << "   Adding material's permittivity\n" << endl;

dimensionedScalar epsilon
(
    "epsilon",
    epsilonR * constant::plasma::epsilon0
);

Info << "   Adding fvOptions\n" << endl;
fv::options fvOptions(gasMesh());

scalar ePotentialNonCoupledResidual = Foam::GREAT;

Info << "   Adding to chargeDensity\n" << endl;
volScalarField chargeDensity
(
    IOobject
    (
        "chargeDensity",
        runTime.timeName(),
        gasMesh(),
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    gasMesh()
);

Info<< "Constructing plasma species\n" << endl;

plasmaSpecies species(gasMesh());
Info << "Species loaded: " << species.speciesNames() << nl;
Info << "Number of species: " << species.nSpecies() << nl << endl;

for (label i = 0; i < species.nSpecies(); ++i)
{
    Info << "  Species " << i << ": "
            << species.speciesNames()[i]
            << ", charge = " << species.speciesCharge()[i] << nl;
}

Info << "\n--- Constructing plasmaTransport ---" << nl << endl;

plasmaTransport transport(species, gasMesh());

Info << "Transport models constructed." << nl << endl;

for (label i = 0; i < species.nSpecies(); ++i)
{
    Info << "  Model for " << species.speciesNames()[i]
            << " : " << transport.model(i).modelName() << nl;
}

Info << "\nSetup complete. Exiting..." << endl;
