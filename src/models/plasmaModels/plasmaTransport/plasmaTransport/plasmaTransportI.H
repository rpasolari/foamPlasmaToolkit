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

//- Return mesh
inline const Foam::fvMesh& 
Foam::plasmaTransport::mesh() const
{
    return mesh_;
}

//- Return species database
inline const Foam::plasmaSpecies& 
Foam::plasmaTransport::species() const
{
    return species_;
}

//- Access non-const model
inline Foam::plasmaTransportModel& 
Foam::plasmaTransport::model(const label i)
{
    return transportModels_[i];
}

//- Access const model
inline const Foam::plasmaTransportModel&
Foam::plasmaTransport::model(const label i) const
{
    return transportModels_[i];
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
