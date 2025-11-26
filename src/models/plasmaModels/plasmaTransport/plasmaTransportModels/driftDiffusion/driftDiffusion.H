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
    Foam::driftDiffusion
    
Description
    Drift–diffusion transport model for plasma species.

    This model implements the standard drift–diffusion formulation, combining
    mobility-driven drift in the electric field with concentration-driven
    diffusion. It updates the species fluxes and densities according to the
    provided mobility and diffusion coefficients.

    The model is activated by specifying:

        transportModel driftDiffusion;

    in the corresponding species dictionary.

SourceFiles
    driftDiffusion.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_driftDiffusion_H
#define Foam_driftDiffusion_H

#include "plasmaTransportModel.H"
#include "plasmaMobilityModel.H"
#include "plasmaDiffusivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class driftDiffusion Declaration
\*---------------------------------------------------------------------------*/

class driftDiffusion
:
    public plasmaTransportModel
{
private:

    //- Mobility field μ_i(x)
    volScalarField mobility_;

    //- Diffusivity field D_i(x)
    volScalarField diffusivity_;

    //- Mobility/Diffusivity submodels
    autoPtr<plasmaMobilityModel> mobilityModel_;
    autoPtr<plasmaDiffusivityModel> diffusivityModel_;

    //- Construct mobility and diffusivity models
    void constructModels();

public:

    //- Runtime type information
    TypeName("driftDiffusion");

    // Constructors

        //- Construct from components
        driftDiffusion
        (
            const word& modelName,
            const dictionary& dict,
            const fvMesh& mesh,
            const plasmaSpecies& species,
            const label specieIndex
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
