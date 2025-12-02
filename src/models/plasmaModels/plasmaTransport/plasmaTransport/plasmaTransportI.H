/*---------------------------------------------------------------------------*\
  File: plasmaTransportI.H
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Inline member function definitions for the plasmaTransport class.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#ifndef Foam_plasmaTransportI_H
#define Foam_plasmaTransportI_H

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Mesh and species metadata access
inline const Foam::fvMesh& 
Foam::plasmaTransport::mesh() const
{
    return mesh_;
}

inline const Foam::plasmaSpecies& 
Foam::plasmaTransport::species() const
{
    return species_;
}


//- Model access
inline Foam::plasmaTransportModel& 
Foam::plasmaTransport::model(const label i)
{
    return transportModels_[i];
}

inline const Foam::plasmaTransportModel&
Foam::plasmaTransport::model(const label i) const
{
    return transportModels_[i];
}

//- Flux access
inline const Foam::volVectorField&
Foam::plasmaTransport::flux(Foam::label i) const
{
    return flux_[i];
}

inline Foam::volVectorField&
Foam::plasmaTransport::flux(Foam::label i)
{
    return flux_[i];
}

inline const Foam::PtrList<Foam::volVectorField>&
Foam::plasmaTransport::flux() const
{
    return flux_;
}

inline Foam::PtrList<Foam::volVectorField>&
Foam::plasmaTransport::flux()
{
    return flux_;
}

//- Surface Flux access
inline const Foam::surfaceScalarField&
Foam::plasmaTransport::surfaceFlux(Foam::label i) const
{
    return surfaceFlux_[i];
}

inline Foam::surfaceScalarField&
Foam::plasmaTransport::surfaceFlux(Foam::label i)
{
    return surfaceFlux_[i];
}

inline const Foam::PtrList<Foam::surfaceScalarField>&
Foam::plasmaTransport::surfaceFlux() const
{
    return surfaceFlux_;
}

inline Foam::PtrList<Foam::surfaceScalarField>&
Foam::plasmaTransport::surfaceFlux()
{
    return surfaceFlux_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
