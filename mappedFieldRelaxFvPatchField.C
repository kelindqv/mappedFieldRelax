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

\*---------------------------------------------------------------------------*/

#include "mappedFieldRelaxFvPatchField.H"
#include "volFields.H"
#include "interpolationCell.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::mappedFieldRelaxFvPatchField<Type>::mappedFieldRelaxFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    mappedPatchBase(p.patch()),
    mappedPatchFieldBase<Type>(*this, *this),
    step_(1000),
    decrement_(1),
    nextUpd_(1000)
{}


template<class Type>
Foam::mappedFieldRelaxFvPatchField<Type>::mappedFieldRelaxFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict),
    mappedPatchBase(p.patch(), dict),
    mappedPatchFieldBase<Type>(*this, *this, dict),
    step_(readInt(dict.lookup("initStep"))),
    decrement_(readScalar(dict.lookup("decrement")))
{
    step_ = max(100,step_);
    nextUpd_ = step_;
}


template<class Type>
Foam::mappedFieldRelaxFvPatchField<Type>::mappedFieldRelaxFvPatchField
(
    const mappedFieldRelaxFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    mappedPatchBase(p.patch(), ptf),
    mappedPatchFieldBase<Type>(*this, *this, ptf),
    step_(ptf.step_),
    decrement_(ptf.decrement_),
    nextUpd_(ptf.nextUpd_)
{}


template<class Type>
Foam::mappedFieldRelaxFvPatchField<Type>::mappedFieldRelaxFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,

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
    const int step,
    const scalar decrement,
    const int nextUpd
)
:
    fixedValueFvPatchField<Type>(p, iF),
    mappedPatchBase
    (
        p.patch(),
        sampleRegion,
        sampleMode,
        samplePatch,
        distance
    ),
    mappedPatchFieldBase<Type>
    (
        *this,
        *this,
        fieldName,
        setAverage,
        average,
        interpolationScheme
    ),
    step_(step),
    decrement_(decrement),
    nextUpd_(nextUpd)
{}


template<class Type>
Foam::mappedFieldRelaxFvPatchField<Type>::mappedFieldRelaxFvPatchField
(
    const mappedFieldRelaxFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    mappedPatchBase(ptf.patch().patch(), ptf),
    mappedPatchFieldBase<Type>(ptf),
    step_(ptf.step_),
    decrement_(ptf.decrement_),
    nextUpd_(ptf.nextUpd_)
{}


template<class Type>
Foam::mappedFieldRelaxFvPatchField<Type>::mappedFieldRelaxFvPatchField
(
    const mappedFieldRelaxFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    mappedPatchBase(ptf.patch().patch(), ptf),
    mappedPatchFieldBase<Type>(*this, *this, ptf),
    step_(ptf.step_),
    decrement_(ptf.decrement_),
    nextUpd_(ptf.nextUpd_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::mappedFieldRelaxFvPatchField<Type>::updateCoeffs()
{
    // Only update if needed and when iteration is evenly div. by 1000
    if
    (
        this->updated() || 
        static_cast<int>(
            this->patchField_.patch().boundaryMesh().mesh().
            time().value()
        ) != nextUpd_
//        // Update at regular intervals:
//        (
//            (
//                static_cast<int>
//                (
//                    this->patchField_.patch().boundaryMesh().mesh().
//                    time().value()
//                ) % period_
//            ) > 0
//        ) ||
//        // Don't update when Time==0
//        static_cast<int>
//        (
//            this->patchField_.patch().boundaryMesh().mesh().
//            time().value()
//        ) == 0
//        // Update only on first iteration (and do restarts with write to disk in between)
////        (
////            !firstTime_ // not implemented
////        )
    )
    {
        return;
    }

    this->operator==(this->mappedField());
//    this->setPeriod(this->period_ - this->decrement_);
//    period_ = period_ - decrement_;
    step_ = max(100,static_cast<int>(step_*decrement_));
    nextUpd_ = nextUpd_ + step_;

    if (debug)
    {
        Info<< "operating on field:" << this->internalField().name()
            << " patch:" << this->patch().name()
            << "  avg:" << gAverage(*this)
            << "  min:" << gMin(*this)
            << "  max:" << gMax(*this)
            << endl;
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::mappedFieldRelaxFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    mappedPatchBase::write(os);
    mappedPatchFieldBase<Type>::write(os);
    os.writeKeyword("initStep") << step_ << token::END_STATEMENT << nl;
    os.writeKeyword("decrement") << decrement_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// ************************************************************************* //
