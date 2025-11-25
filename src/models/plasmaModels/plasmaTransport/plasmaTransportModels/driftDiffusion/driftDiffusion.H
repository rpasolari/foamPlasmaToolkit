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
    Drift–diffusion plasma transport model.

    This class implements a drift–diffusion–based transport formulation
    for plasma species. It is selected at runtime using the standard
    runtime-selection mechanism and operates directly on the species data
    provided by plasmaSpecies.

    The model participates in advancing the species densities by computing
    transport fluxes based on drift and diffusion processes, using the
    general interface defined in plasmaTransportModel.

SourceFiles
    driftDiffusion.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_driftDiffusion_H
#define Foam_driftDiffusion_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class driftDiffusion Declaration
\*---------------------------------------------------------------------------*/

class driftDiffusion
{
protected:



public:


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
