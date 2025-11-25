/*---------------------------------------------------------------------------*\
  File: mobilityModel.C
  Part of: foamPlasmaToolkit
\*---------------------------------------------------------------------------*/

#include "mobilityModel.H"

namespace Foam
{

defineTypeNameAndDebug(mobilityModel, 0);
defineRunTimeSelectionTable(mobilityModel, dictionary);


// Selector
autoPtr<mobilityModel> mobilityModel::New
(
    const word& modelName,
    const dictionary& coeffs,
    const fvMesh& mesh
)
{
    auto ctorIter = dictionaryConstructorTablePtr_->cfind(modelName);

    if (!ctorIter.good())
    {
        FatalIOErrorInLookup
        (
            coeffs,
            mobilityModel::typeName,
            modelName,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<mobilityModel>(ctorIter()(coeffs, mesh));
}

} // End namespace Foam

// ************************************************************************* //
