/*---------------------------------------------------------------------------*\
  File: createFluidFields.H
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Derived from createFluidFields.H in OpenFOAM v2412.
      Defines and initializes the fluid-region fields for each region.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

PtrList<volScalarField> ePotentialFluid(fluidRegions.size());
PtrList<volScalarField> eChargeFluid(fluidRegions.size());
PtrList<surfaceScalarField> surfChargeFluid(fluidRegions.size());
PtrList<IOdictionary> propertiesDictFluid(fluidRegions.size());
PtrList<dimensionedScalar> epsilonFluid(fluidRegions.size());

List<scalar> epsilonRFluid(fluidRegions.size());

List<scalar> ePotentialNonCoupledResidualFluid
(
    fluidRegions.size(), Foam::GREAT
);

forAll(fluidRegions, i)
{
    Info << "*** Reading fields for region "
         << fluidRegions[i].name() << nl << endl;

    Info << "   Adding to ePotentialFluid\n" << endl;
    ePotentialFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "ePotential",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info << "   Adding to eChargeFluid\n" << endl;
    eChargeFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "eCharge",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    {
        IOobject surfChargeHeader
        (
            "surfCharge",
            runTime.timeName(),
            fluidRegions[i],
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        );

        if (surfChargeHeader.typeHeaderOk<surfaceScalarField>())
        {
            surfChargeFluid.set
            (
                i,
                new surfaceScalarField
                (
                    IOobject
                    (
                        "surfCharge",
                        runTime.timeName(),
                        fluidRegions[i],
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    fluidRegions[i]
                )
            );

            regIOobject::store(&surfChargeFluid[i]);
            Info << "   Adding to surfChargeFluid\n" << endl;
        }
    }

    propertiesDictFluid.set
    (
        i,
        new IOdictionary
        (
            IOobject
            (
                "properties",
                runTime.constant(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );

    epsilonRFluid[i] = propertiesDictFluid[i].lookupOrDefault<scalar>
    (
        "dielectricConstant", 1.0
    );

    Info << "   Adding material's permittivity\n" << endl;

    epsilonFluid.set
    (
        i,
        new dimensionedScalar
        (
            "epsilonFluid",
            dimensionSet(-1,-3,4,0,0,2,0),
            epsilonRFluid[i] * 8.8541878128e-12
        )
    );
}


