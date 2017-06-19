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

#include "mappedFieldRelaxStagFvPatchField.H"
#include "volFields.H"
#include "interpolationCell.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::mappedFieldRelaxStagFvPatchField<Type>::mappedFieldRelaxStagFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    mappedPatchBase(p.patch()),
    mappedPatchFieldBase<Type>(*this, *this),
    period_(1000),
//    step_(1000),
//    decrement_(1),
//    nextUpd_(1000),
    lag_(500)
{}


template<class Type>
Foam::mappedFieldRelaxStagFvPatchField<Type>::mappedFieldRelaxStagFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict),
    mappedPatchBase(p.patch(), dict),
    mappedPatchFieldBase<Type>(*this, *this, dict),
    period_(readInt(dict.lookup("period"))),
    lag_(readInt(dict.lookup("lag")))
//    step_(readInt(dict.lookup("initStep"))),
//    decrement_(readScalar(dict.lookup("decrement")))
{
//    step_ = max(100,step_);
//    lag_ = step_/2;
//    nextUpd_ = step_;
}


template<class Type>
Foam::mappedFieldRelaxStagFvPatchField<Type>::mappedFieldRelaxStagFvPatchField
(
    const mappedFieldRelaxStagFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    mappedPatchBase(p.patch(), ptf),
    mappedPatchFieldBase<Type>(*this, *this, ptf),
    period_(ptf.period_),
//    step_(ptf.step_),
//    decrement_(ptf.decrement_),
//    nextUpd_(ptf.nextUpd_),
    lag_(ptf.lag_)
{}


template<class Type>
Foam::mappedFieldRelaxStagFvPatchField<Type>::mappedFieldRelaxStagFvPatchField
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
    const int period,
//    const int step,
//    const scalar decrement,
//    const int nextUpd,
    const int lag
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
    period_(period),
//    step_(step),
//    decrement_(decrement),
//    nextUpd_(nextUpd),
    lag_(lag)
{}


template<class Type>
Foam::mappedFieldRelaxStagFvPatchField<Type>::mappedFieldRelaxStagFvPatchField
(
    const mappedFieldRelaxStagFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    mappedPatchBase(ptf.patch().patch(), ptf),
    mappedPatchFieldBase<Type>(ptf),
    period_(ptf.period_),
//    step_(ptf.step_),
//    decrement_(ptf.decrement_),
//    nextUpd_(ptf.nextUpd_),
    lag_(ptf.lag_)
{}


template<class Type>
Foam::mappedFieldRelaxStagFvPatchField<Type>::mappedFieldRelaxStagFvPatchField
(
    const mappedFieldRelaxStagFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    mappedPatchBase(ptf.patch().patch(), ptf),
    mappedPatchFieldBase<Type>(*this, *this, ptf),
    period_(ptf.period_),
//    step_(ptf.step_),
//    decrement_(ptf.decrement_),
//    nextUpd_(ptf.nextUpd_),
    lag_(ptf.lag_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//// Overload the mapped field - if pressure is too big, scale it back
//template<class Type>
//tmp<Field<Type>> Foam::mappedFieldRelaxStagFvPatchField<Type>::mappedField() const
//{
//    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

//    // Since we're inside initEvaluate/evaluate there might be processor
//    // comms underway. Change the tag we use.
//    int oldTag = UPstream::msgType();
//    UPstream::msgType() = oldTag + 1;

//    const fvMesh& thisMesh = patchField_.patch().boundaryMesh().mesh();
//    const fvMesh& nbrMesh = refCast<const fvMesh>(mapper_.sampleMesh());

//    // Result of obtaining remote values
//    tmp<Field<Type>> tnewValues(new Field<Type>(0));
//    Field<Type>& newValues = tnewValues.ref();

//    switch (mapper_.mode())
//    {
//        case mappedPatchBase::NEARESTCELL:
//        {
//            const mapDistribute& distMap = mapper_.map();

//            if (interpolationScheme_ != interpolationCell<Type>::typeName)
//            {
//                // Send back sample points to the processor that holds the cell
//                vectorField samples(mapper_.samplePoints());
//                distMap.reverseDistribute
//                (
//                    (
//                        mapper_.sameRegion()
//                      ? thisMesh.nCells()
//                      : nbrMesh.nCells()
//                    ),
//                    point::max,
//                    samples
//                );

//                autoPtr<interpolation<Type>> interpolator
//                (
//                    interpolation<Type>::New
//                    (
//                        interpolationScheme_,
//                        sampleField()
//                    )
//                );
//                const interpolation<Type>& interp = interpolator();

//                newValues.setSize(samples.size(), pTraits<Type>::max);
//                forAll(samples, celli)
//                {
//                    if (samples[celli] != point::max)
//                    {
//                        newValues[celli] = interp.interpolate
//                        (
//                            samples[celli],
//                            celli
//                        );
//                    }
//                }
//            }
//            else
//            {
//                newValues = sampleField();
//            }

//            distMap.distribute(newValues);

//            break;
//        }
//        case mappedPatchBase::NEARESTPATCHFACE:
//        case mappedPatchBase::NEARESTPATCHFACEAMI:
//        {
//            const label nbrPatchID =
//                nbrMesh.boundaryMesh().findPatchID(mapper_.samplePatch());

//            if (nbrPatchID < 0)
//            {
//                FatalErrorInFunction
//                 << "Unable to find sample patch " << mapper_.samplePatch()
//                 << " in region " << mapper_.sampleRegion()
//                 << " for patch " << patchField_.patch().name() << nl
//                 << abort(FatalError);
//            }

//            const fieldType& nbrField = sampleField();

//            newValues = nbrField.boundaryField()[nbrPatchID];
//            mapper_.distribute(newValues);

//            break;
//        }
//        case mappedPatchBase::NEARESTFACE:
//        {
//            Field<Type> allValues(nbrMesh.nFaces(), Zero);

//            const fieldType& nbrField = sampleField();

//            forAll(nbrField.boundaryField(), patchi)
//            {
//                const fvPatchField<Type>& pf =
//                    nbrField.boundaryField()[patchi];
//                label faceStart = pf.patch().start();

//                forAll(pf, facei)
//                {
//                    allValues[faceStart++] = pf[facei];
//                }
//            }

//            mapper_.distribute(allValues);
//            newValues.transfer(allValues);

//            break;
//        }
//        default:
//        {
//            FatalErrorInFunction
//             << "Unknown sampling mode: " << mapper_.mode()
//             << nl << abort(FatalError);
//        }
//    }

//    if (setAverage_)
//    {
////        Type averagePsi =
////            gSum(patchField_.patch().magSf()*newValues)
////           /gSum(patchField_.patch().magSf());
//        Type tehMax = gMax(newValues)

//        if (tehMax > mag(average_))
//        {
//            newValues -= average_;
//        }
//    }

//    // Restore tag
//    UPstream::msgType() = oldTag;

//    return tnewValues;
//}

template<class Type>
void Foam::mappedFieldRelaxStagFvPatchField<Type>::updateCoeffs()
{
    // Only update if needed and when iteration is evenly div. by 1000
    if
    (
        this->updated() || 
//        static_cast<int>(
//            this->patchField_.patch().boundaryMesh().mesh().
//            time().value()
//        ) != (nextUpd_ + lag_)
        (
            (
                (
                    static_cast<int>
                    (
                        this->patchField_.patch().boundaryMesh().mesh().
                        time().value()
                    )
                    - lag_
                ) % period_
            ) > 0
        ) ||
        // Don't update when Time==0
        static_cast<int>
        (
            this->patchField_.patch().boundaryMesh().mesh().
            time().value()
        ) == 0
    )
    {
        return;
    }

    this->operator==(this->mappedField());
//    step_ = max(100,static_cast<int>(step_*decrement_));
//    nextUpd_ = nextUpd_ + step_;
//    lag_ = step_/2;

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
void Foam::mappedFieldRelaxStagFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    mappedPatchBase::write(os);
    mappedPatchFieldBase<Type>::write(os);
    os.writeKeyword("period") << period_ << token::END_STATEMENT << nl;
    os.writeKeyword("lag") << lag_ << token::END_STATEMENT << nl;
//    os.writeKeyword("initStep") << step_ << token::END_STATEMENT << nl;
//    os.writeKeyword("decrement") << decrement_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// ************************************************************************* //
