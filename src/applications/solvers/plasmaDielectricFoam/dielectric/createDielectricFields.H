/*---------------------------------------------------------------------------*\
  File: createDielectricFields.H
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Derived from createSolidFields.H in OpenFOAM v2412.
      Defines and initializes the dielectric-region fields for each region.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

PtrList<volScalarField> ePotentialDielectricList(dielectricRegions.size());
PtrList<volVectorField> EDielectricList(dielectricRegions.size());
PtrList<IOdictionary>
    electricPropertiesDictDielectricList(dielectricRegions.size());
PtrList<surfaceScalarField> surfChargeDielectricList(dielectricRegions.size());
PtrList<dimensionedScalar> epsilonDielectricList(dielectricRegions.size());

List<scalar> epsilonRDielectricList(dielectricRegions.size());

List<scalar> ePotentialNonCoupledResidualDielectricList
(
    dielectricRegions.size(), Foam::GREAT
);

forAll(dielectricRegions, i)
{
    Info << "*** Reading fields for region "
         << dielectricRegions[i].name() << nl << endl;

    Info << "   Adding to ePotentialDielectric\n" << endl;
    ePotentialDielectricList.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "ePotential",
                runTime.timeName(),
                dielectricRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            dielectricRegions[i]
        )
    );

    Info << "   Adding to EDielectric\n" << endl;
    EDielectricList.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "E",
                runTime.timeName(),
                dielectricRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            dielectricRegions[i],
            dimensionedVector
            (
                "zero",
                dimensionSet(1, 1, -3, 0, 0, -1, 0),
                vector::zero
            )
        )
    );

    {
        IOobject surfChargeHeader
        (
            "surfCharge",
            runTime.timeName(),
            dielectricRegions[i],
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        );

        if (surfChargeHeader.typeHeaderOk<surfaceScalarField>())
        {
            surfChargeDielectricList.set
            (
                i,
                new surfaceScalarField
                (
                    IOobject
                    (
                        "surfCharge",
                        runTime.timeName(),
                        dielectricRegions[i],
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    dielectricRegions[i]
                )
            );

            regIOobject::store(&surfChargeDielectricList[i]);
            Info << "   Adding to surfChargeDielectric\n" << endl;
        }
    }

    electricPropertiesDictDielectricList.set
    (
        i,
        new IOdictionary
        (
            IOobject
            (
                "electricProperties",
                runTime.constant(),
                dielectricRegions[i],
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );

    epsilonRDielectricList[i] =
        electricPropertiesDictDielectricList[i].lookupOrDefault<scalar>
        (
            "dielectricConstant",
            1.0
        );

    Info << "   Adding material's permittivity\n" << endl;
    epsilonDielectricList.set
    (
        i,
        new dimensionedScalar
        (
            "epsilonDielectric",
            epsilonRDielectricList[i] * constant::plasma::epsilon0
        )
    );
}
