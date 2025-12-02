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
    Foam::plasmaTimeControl

Description
    Scharfetter-Gummel finite volume operator for drift-diffusion.

    This function builds a coupled fvMatrix for the transport equation:
        div(J) = Source
    Where the total flux J is defined as:
        J = phiDrift * psi - D * grad(psi)

    This scheme solves the 1D convection-diffusion equation analytically 
    between cell centers (assuming constant E and D locally) to determine 
    the exponential profile of the scalar field. This ensures positivity 
    of the species density even in regions of high electric fields (sheaths) 
    where standard central or upwind schemes may fail or produce oscillations

SourceFiles
    
\*---------------------------------------------------------------------------*/

#ifndef ScharfetterGummel_H
#define ScharfetterGummel_H

#include "fvMatrix.H"
#include "fvm.H"
#include "fvc.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "cmath"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Function ScharfetterGummel
\*---------------------------------------------------------------------------*/

namespace fvm
{

    template<class Type>
    tmp<fvMatrix<Type>> ScharfetterGummel
    (
        const GeometricField<Type, fvPatchField, volMesh>& rho,
        const surfaceScalarField& phiDrift,
        const volScalarField& D
    )
    {
        const fvMesh& mesh = rho.mesh();

        // Initialize the matrix
        tmp<fvMatrix<Type>> tmat
            (
                new fvMatrix<Type>(rho, dimVol/dimTime * rho.dimensions())
            );
        fvMatrix<Type>& mat = tmat.ref();

        const surfaceScalarField& magSf = mesh.magSf();
        const surfaceScalarField& deltaCoeffs = mesh.deltaCoeffs();

        const labelList& owner = mesh.owner();
        const labelList& neighbour = mesh.neighbour();

        surfaceScalarField Df = fvc::interpolate(D);

        // Bernoulli Function Definition: B(x) = x / exp(x) - 1
        auto Bern = [](scalar x) -> scalar
        {
            if (mag(x) < 1e-4) return 1.0 - x/2.0;
            if (x > 80)        return 0.0;
            if (x < 80)        return -x;
            return x / (exp(x) - 1.0); 
        };

        // Internal Face loop
        for (label facei = 0; facei < owner.size(); ++facei)
        {
            label P = owner[facei];
            label N = neighbour[facei];
            scalar faceArea = magSf[facei];

            scalar diffCond = Df[facei] * faceArea * deltaCoeffs[facei];
            scalar driftFlux = phiDrift[facei];
            scalar Pe = driftFlux / (diffCond + VSMALL);

            scalar coeffP = diffCond * Bern(Pe);
            scalar coeffN = diffCond * Bern(-Pe);

            mat.diag()[P] += coeffP;
            mat.diag()[N] += coeffN;

            mat.upper()[facei] -= coeffN;
            mat.lower()[facei] -= coeffP;
        }

        // Boundary condition handling
        forAll(rho.boundaryField(), patchi)
        {
            const fvPatchField<Type>& ptf = rho.boundaryField()[patchi];

            if (!ptf.coupled())
            {
                mat.internalCoeffs()[patchi] = ptf.valueInternalCoeffs(rho);
                mat.boundaryCoeffs()[patchi] = ptf.valueBoundaryCoeffs(rho);
            }
        }

        return tmat;
    }

    // Overload for diffusivity as a dimensionedScalar
    template<class Type>
    tmp<fvMatrix<Type>> ScharfetterGummel
    (
        const GeometricField<Type, fvPatchField, volMesh>& rho,
        const surfaceScalarField& phiDrift,
        const dimensionedScalar& D    
    )
    {
        volScalarField DField
        (
            IOobject
            (
                "D_tmp", 
                rho.mesh().time().timeName(), 
                rho.mesh(), 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ),
            rho.mesh(),
            D
        );

        // Call the core implementation
        return ScharfetterGummel(rho, phiDrift, DField);
    }

    // Overload for drift velocity instead of flux
    template<class Type>
    tmp<fvMatrix<Type>> ScharfetterGummel
    (
        const GeometricField<Type, fvPatchField, volMesh>& rho,
        const volVectorField& vDrift,
        const volScalarField& D
    )
    {
        // Convert vDrift to phiDrift (Flux)
        surfaceScalarField phiDrift = fvc::flux(vDrift);
        
        // Call the core implementation
        return ScharfetterGummel(rho, phiDrift, D);
    }

    // Overload for Electric Field (E) and Mobility (mu) instead of flux
    template<class Type>
    tmp<fvMatrix<Type>> ScharfetterGummel
    (
        const GeometricField<Type, fvPatchField, volMesh>& rho,
        const volVectorField& E,
        const volScalarField& mu,
        const volScalarField& D
    )
    {
        volVectorField vDrift = mu * E;
        
        surfaceScalarField phiDrift = fvc::flux(vDrift);
        
        // Call the core implementation
        return ScharfetterGummel(rho, phiDrift, D);
    }

} // End namespace fvm

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
