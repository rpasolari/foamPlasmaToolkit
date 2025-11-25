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
    Foam::noneTransportModel

Description
    A null transport model for plasma species.

    This model provides a do-nothing implementation of plasmaTransportModel.
    It is used when a species has no active transport behaviour and serves
    as a simple placeholder within the run-time selection framework.

    The model is activated by specifying:

        transportModel none;

    in the corresponding species dictionary.

SourceFiles
    noneTransportModel.C
    
\*---------------------------------------------------------------------------*/

#ifndef Foam_noneTransportModel_H
#define Foam_noneTransportModel_H

#include "plasmaTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class noneTransportModel Declaration
\*---------------------------------------------------------------------------*/

class noneTransportModel
:
    public plasmaTransportModel
{
public:

    //- Runtime type information
    TypeName("none");

    // Constructors

        //- Construct from components
        noneTransportModel
        (
            const word& modelName,
            const dictionary& dict,
            const fvMesh& mesh,
            const plasmaSpecies& species
        );

    // Member Functions
        
        //- Correct model
        virtual void correct() override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
