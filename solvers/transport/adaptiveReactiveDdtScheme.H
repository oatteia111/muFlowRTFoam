#ifndef adaptiveReactiveDdtScheme_H
#define adaptiveReactiveDdtScheme_H

#include "fvMesh.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "EulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"

namespace Foam
{

template<class Type>
class adaptiveReactiveDdtScheme
{
    const fvMesh& mesh_;
    const volScalarField& reactIndicator_; // 0 = phase réactive, 1 = transport

    // Schémas internes
    EulerDdtScheme<Type> euler_;
    CrankNicolsonDdtScheme<Type> cn_;

public:

    // Constructeur
    adaptiveReactiveDdtScheme
    (
        const fvMesh& mesh,
        const word& fieldName,
        const volScalarField& reactIndicator
    )
    :
        mesh_(mesh),
        reactIndicator_(reactIndicator),
        euler_(mesh, fieldName),
        cn_(mesh, fieldName, 0.5) // facteur CN
    {}

    // Retourne le terme fvm::ddt en fonction du champ indicateur
    tmp<GeometricField<Type, fvPatchField, volMesh>>
    fvmDdt(GeometricField<Type, fvPatchField, volMesh>& vf) const
    {
        if (reactIndicator_.internalField().max() < 0.5)
        {
            // phase réactive → Euler
            return euler_.fvmDdt(vf);
        }
        else
        {
            // transport post-réaction → Crank-Nicolson
            return cn_.fvmDdt(vf);
        }
    }
};

}

#endif
