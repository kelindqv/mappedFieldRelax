/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::mappedFieldRelaxFvPatchField

Group
    grpGenericBoundaryConditions grpCoupledBoundaryConditions

Description
    A mappedField-derived BC which updates the values every initStep iterations
    with and optional (0<x<=1) decrement on each update.

Usage
    \table
        Property     | Description             | Required    | Default value
        fieldName    | name of field to be mapped | no       | this field name
        setAverage   | flag to activate setting of average value | yes |
        average      | average value to apply if \c setAverage = yes | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            mappedFieldRelax;
        sampleMode      nearestPatchFace;
        sampleRegion    fluid;             // For multi-region cases
        samplePatch     outlet;
        offsetMode      uniform;
        offset          (0.079 0 0);       // Distance to samplePatch
        field           U;
        setAverage      1;
        average         (1.0 0 0);
        interpolationScheme cell;
        period          1500;              // # of iterations btw updates
        value           uniform (1.0 0 0); 
    }
    \endverbatim

Note
    Since this condition can be applied on a per-field and per-patch basis,
    it is possible to duplicate the mapping information.  If possible, employ
    the \c mapped condition in preference to avoid this situation, and only
    employ this condition if it is not possible to change the underlying
    geometric (poly) patch type to \c mapped.

See also
    Foam::mappedPatchBase
    Foam::mappedPolyPatch
    Foam::mappedFvPatch
    Foam::fixedValueFvPatchField

SourceFiles
    mappedFieldRelaxFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef mappedFieldRelaxFvPatchField_H
#define mappedFieldRelaxFvPatchField_H

#include "mappedPatchBase.H"
//#include "mymappedPatchFieldBase.H"
#include "mappedPatchFieldBase.H"
#include "fixedValueFvPatchFields.H"
#include "interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class mappedFieldRelaxFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class mappedFieldRelaxFvPatchField
:
    public fixedValueFvPatchField<Type>,
    public mappedPatchBase,
//    public mymappedPatchFieldBase<Type>
    public mappedPatchFieldBase<Type>
{

protected:
        
//        //- Current step size
//        int step_;
//        
//        //- How much smaller should the sted size get every update?
//        const scalar decrement_;
//        
//        //- Which iteration is up for update?
//        int nextUpd_;
//        
//        //- Time lag between updates of U and updates of p
//        const int lag_;

        // Period between updates
        const int period_;

public:

    //- Runtime type information
    TypeName("mappedFieldRelax");


    // Constructors

        //- Construct from patch and internal field
        mappedFieldRelaxFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mappedFieldRelaxFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct from patch, internal field and distance for normal type
        //  sampling
        mappedFieldRelaxFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,

            // mappedPatchBase
            const word& sampleRegion,
            const sampleMode sampleMode,
            const word& samplePatch,
            const scalar distance,

            // My settings
            const word& fieldName,
            const bool setAverage,
            const Type average,
            const word& interpolationScheme,
            const int period_
//            const int step,
//            const scalar decrement,
//            const int nextUpd
        );

        //- Construct by mapping given
        //  mappedFieldRelaxFvPatchField
        //  onto a new patch
        mappedFieldRelaxFvPatchField
        (
            const mappedFieldRelaxFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        mappedFieldRelaxFvPatchField
        (
            const mappedFieldRelaxFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new mappedFieldRelaxFvPatchField<Type>
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        mappedFieldRelaxFvPatchField
        (
            const mappedFieldRelaxFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new mappedFieldRelaxFvPatchField<Type>
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mappedFieldRelaxFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
