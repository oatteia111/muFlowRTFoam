/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

#include "myFixedValueConstraint.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "DimensionedField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fv::myFixedValueConstraint<Type>::myFixedValueConstraint
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::fv::myFixedValueConstraint<Type>::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        const dictionary& fieldValuesDict = coeffs_.subDict("fieldValues");
        fieldNames_.setSize(fieldValuesDict.size());
        fieldValues_.setSize(fieldNames_.size());

        label i = 0;
        forAllConstIter(dictionary, fieldValuesDict, iter)
        {
            fieldNames_[i] = iter().keyword();
            fieldValuesDict.lookup(iter().keyword()) >> fieldValues_[i];
            i++;
        }

        applied_.setSize(fieldNames_.size(), false);
		
        return true;
    }
    else
    {
        return false;
    }
}


template<class Type>
void Foam::fv::myFixedValueConstraint<Type>::constrain
(
	//const volScalarField& v0,
	fvMatrix<Type>& eqn,
    const label fieldi
)
{
    DebugInfo
        << "myFixedValueConstraint<"
        << pTraits<Type>::typeName
        << ">::constrain for source " << name_ << endl;

	const GeometricField<Type, fvPatchField, volMesh>& psi = eqn.psi();

	//Info << "cellSetName " << cellSetName_ << " data0 " << cellsData_[0]<< " field " << fieldi <<" value "<< fieldValues_[fieldi] << endl;
	// seach for the right time
	float time = mesh().time().value(); //const Time&
	//int nlay = int(cellsData_[0]);int ncell_lay = int(cellsData_[1]);
	int icount = 0; //Info <<"icount "<< icount << endl;//counter of lines within poption file
 	float tnow = cellsData_[2]*86400;
	while (time>=tnow)
		{
		icount += cells_.size(); 
		tnow = cellsData_[2+icount*4]*86400; //Info <<"icount "<< icount << " topt "<<topt<<endl;
		}
	if (icount>0) {icount -= cells_.size();};
	
 	//for (int i=0; i<cells_.size(); i++) {Info << "val of Cw "<< mag(psi()[cells_[i]])<<endl;}
	//find the correct values and prepare the b list
	//Info << "fix icount "<< icount << " cells_ size " << cells_.size() << endl;
	std::vector<double> a(cells_.size(),0.); //Info <<"popt size "<< poptions.size()<< endl;
	std::vector<double> aAdd(cells_.size(),0.); //Info <<"popt size "<< poptions.size()<< endl;
	List<Type> b;
	for (int i=0; i<cells_.size(); i++) { a[cIndex_[i]] = cellsData_[2+(icount+i)*4+2] ;}
	//if (cellSetName_.substr(0,1) == "c") {for (int i=0; i<cells_.size(); i++) { aAdd[cIndex_[i]] = cellsDataAdd_[2+(icount+i)*4+2] ;}}

	//set the values in the b list
	//for cFix, cWel and phreeqc 
	//Info << "variable name "<< psi().name()<< endl;
	if (cellSetName_.substr(0,1) == "c" && psi().name().length()>2)
		{
		int ispec = std::stoi(psi().name().substr(2,psi().name().length())); //the species number
		for (int i=0; i<cells_.size(); i++) { 
			double conc = cellsSolutions_[a[i]*cellsNspecies_+ispec]; 
			//if (i<100) {Info << "in opt ispec "<<ispec<<" data "<< cellsData_[2+i*4+2] << " conc " << conc << endl;}
			b.append(  conc * fieldValues_[fieldi] ); } //here fieldValues is just one, left for consistency
		}
	else if (cellSetName_.substr(0,1) == "g" && psi().name().length()>2) // gases
		{
		int ispec = std::stoi(psi().name().substr(2,psi().name().length())); //the species number
		for (int i=0; i<cells_.size(); i++) { 
			double conc = cellsGases_[a[i]*cellsGspecies_+ispec]; 
			b.append(  conc * fieldValues_[fieldi] ); //Info<<"spec g "<<ispec<<" conc "<<conc<<" b "<<b<<endl;
			} //here fieldValues is just one, left for consistency
		}
	else 
		{
		for (int i=0; i<cells_.size(); i++) { b.append(a[i]*fieldValues_[fieldi]); }
		}
	//for (int i=0; i<cells_.size(); i++) {Info << cells_[i] << " in fix const "<< b[i] << endl;}
	eqn.setValues(cells_, b);
    //eqn.setValues(cells_, List<Type>(cells_.size(), a[i]*fieldValues_[fieldi]));
}

// ************************************************************************* //
