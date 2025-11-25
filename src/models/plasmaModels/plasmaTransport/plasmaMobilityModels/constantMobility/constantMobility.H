/*---------------------------------------------------------------------------*\
License
    This file is part of the foamPlasmaToolkit.

    The foamPlasmaToolkit is not part of OpenFOAM but is developed using the
    OpenFOAM framework and linked against OpenFOAM libraries.

    Copyright (C) 2025
        Rention Pasolari

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program. If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::constantMobility

Description
    Constant mobility model for plasma species.

    This model returns a uniform mobility value μ throughout the domain,
    specified in the mobilityCoeffs sub-dictionary:

        mobilityCoeffs
        {
            mu   0.1;
        }

    It is selected at runtime via:

        mobilityModel  constant;

SourceFiles
    constantMobility.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_constantMobility_H
#define Foam_constantMobility_H

#include "mobilityModel.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class constantMobility Declaration
\*---------------------------------------------------------------------------*/

class constantMobility
:
    public mobilityModel
{
private:

    //- Constant mobility value
    scalar mu_;

public:

    //- Runtime type information
    TypeName("constant");

    //- Construct from dictionary and mesh
    constantMobility
    (
        const dictionary& coeffs,
        const fvMesh& mesh
    );

    //- Destructor
    virtual ~constantMobility() {}

    //- Return mobility μ
    virtual tmp<volScalarField> mu() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
