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
    Foam::constantDiffusivity

Description
    Constant diffusivity model for plasma species.

    This model returns a uniform diffusivity value D throughout the domain,
    specified in the diffusivityCoeffs sub-dictionary:

    Required dictionary entries:
    \verbatim
        diffusivityCoeffs
        {
            D   <scalar>;
        }
    \verbatim

    It is selected at runtime via:

        diffusivityModel  constant;

SourceFiles
    constantDiffusivity.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_constantDiffusivity_H
#define Foam_constantDiffusivity_H

#include "plasmaDiffusivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class constantDiffusivity Declaration
\*---------------------------------------------------------------------------*/

class constantDiffusivity
:
    public plasmaDiffusivityModel
{
private:

    //- Constant diffusivity value [m2/s]
    scalar D0_;

public:

    //- Runtime type information
    TypeName("constant");

    // Constructors

        //- Construct from components
        constantDiffusivity
        (
            const word& modelName,
            const dictionary& dict,
            const fvMesh& mesh,
            const plasmaSpecies& species,
            const label specieIndex
        );

    // Member Functions
        
        //- Correct model
        virtual void correct(volScalarField& D) const override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
