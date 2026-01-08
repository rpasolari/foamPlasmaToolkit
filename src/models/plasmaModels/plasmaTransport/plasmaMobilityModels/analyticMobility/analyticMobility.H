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
    Foam::analyticMobility

Description
    Analytic mobility model for plasma species based on empirical power-law 
    functions of the electric field magnitude.

    This model calculates the species mobility field (\mu) throughout the 
    domain using the local electric field strength |E|. 

    The mobility is computed as:
        mu = muCoeff * |E|^fieldExp

    Required dictionary entries:
    \verbatim
        mobilityCoeffs
        {
            muCoeff    <scalar>;
            fieldExp   <scalar>;
        }
    \endverbatim

    It is selected at runtime via:

        mobilityModel  analytic;

SourceFiles
    analyticMobility.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_analyticMobility_H
#define Foam_analyticMobility_H

#include "plasmaMobilityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class analyticMobility Declaration
\*---------------------------------------------------------------------------*/

class analyticMobility
:
    public plasmaMobilityModel
{
private:

    //- Pre-exponential factor (Coefficient)
    scalar muCoeff_;

    //- Field dependency exponent
    scalar fieldExp_;

public:

    //- Runtime type information
    TypeName("analytic");

    // Constructors

        //- Construct from components
        analyticMobility
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
