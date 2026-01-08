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
    Foam::analyticDiffusivity

Description
    Analytic diffusivity model for plasma species based on empirical power-law 
    functions of the electric field magnitude.

    This model calculates the species diffusivity field (\D) throughout the 
    domain using the local electric field strength |E|. 

    The diffusivity is computed as:
        D = DCoeff * |E|^fieldExp

        diffusivityCoeffs
        {
            D   0.1;
        }

    Required dictionary entries:
    \verbatim
        diffusivityCoeffs
        {
            DCoeff    <scalar>;
            fieldExp   <scalar>;
        }
    \endverbatim

    It is selected at runtime via:

        diffusivityModel  analytic;

SourceFiles
    analyticDiffusivity.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_analyticDiffusivity_H
#define Foam_analyticDiffusivity_H

#include "plasmaDiffusivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class analyticDiffusivity Declaration
\*---------------------------------------------------------------------------*/

class analyticDiffusivity
:
    public plasmaDiffusivityModel
{
private:

    //- Pre-exponential factor (Coefficient)
    scalar DCoeff_;

    //- Field dependency exponent
    scalar fieldExp_;

public:

    //- Runtime type information
    TypeName("analytic");

    // Constructors

        //- Construct from components
        analyticDiffusivity
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
