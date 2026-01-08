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

    Required dictionary entries:
    \verbatim
        mobilityCoeffs
        {
            mu   <scalar>;
        }
    \verbatim

    It is selected at runtime via:

        mobilityModel  constant;

SourceFiles
    constantMobility.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_constantMobility_H
#define Foam_constantMobility_H

#include "plasmaMobilityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class constantMobility Declaration
\*---------------------------------------------------------------------------*/

class constantMobility
:
    public plasmaMobilityModel
{
private:

    //- Constant mobility value [m2/(V·s)]
    scalar mu0_;

public:

    //- Runtime type information
    TypeName("constant");

    // Constructors

        //- Construct from components
        constantMobility
        (
            const word& modelName,
            const dictionary& dict,
            const fvMesh& mesh,
            const plasmaSpecies& species,
            const label specieIndex
        );

    // Member Functions
        
        //- Correct model
        virtual void correct(volScalarField& mu) const override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
