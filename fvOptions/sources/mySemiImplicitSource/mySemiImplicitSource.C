/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
#include "mySemiImplicitSource.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "fvmSup.H"

#include <stdlib.h>  //OA
#include <set>
#include <vector>  //OA
#include <string>

/*
//"""""""""""""""" to be able to read options file """"""""""""""""""""""
#include <iostream>
#include <fstream>
#include <iterator>
#include <sstream>

#include <string>
#include <unistd.h>
#define GetCurrentDir getcwd
std::string get_current_dir() {
   char buff[FILENAME_MAX]; //create string buffer to hold path
   GetCurrentDir( buff, FILENAME_MAX );
   std::string current_working_dir(buff);
   return current_working_dir;
}
std::string cur_dir = get_current_dir();
*/

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

template<class Type>
const Foam::wordList Foam::fv::mySemiImplicitSource<Type>::volumeModeTypeNames_
(
    IStringStream("(absolute specific)")()
);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
typename Foam::fv::mySemiImplicitSource<Type>::volumeModeType
Foam::fv::mySemiImplicitSource<Type>::wordToVolumeModeType
(
    const word& vmtName
) const
{
    forAll(volumeModeTypeNames_, i)
    {
        if (vmtName == volumeModeTypeNames_[i])
        {
            return volumeModeType(i);
        }
    }

    FatalErrorInFunction
        << "Unknown volumeMode type " << vmtName
        << ". Valid volumeMode types are:" << nl << volumeModeTypeNames_
        << exit(FatalError);

    return volumeModeType(0);
}


template<class Type>
Foam::word Foam::fv::mySemiImplicitSource<Type>::volumeModeTypeToWord
(
    const volumeModeType& vmtType
) const
{
    if (vmtType > volumeModeTypeNames_.size())
    {
        return "UNKNOWN";
    }
    else
    {
        return volumeModeTypeNames_[vmtType];
    }
}


template<class Type>
void Foam::fv::mySemiImplicitSource<Type>::setFieldData(const dictionary& dict)
{
    fieldNames_.setSize(dict.toc().size());
    injectionRate_.setSize(fieldNames_.size());

    applied_.setSize(fieldNames_.size(), false);

    label i = 0;
 	forAllConstIter(dictionary, dict, iter)
    {
        fieldNames_[i] = iter().keyword();
        dict.lookup(iter().keyword()) >> injectionRate_[i];
        i++;
    }
    // Set volume normalisation
    if (volumeMode_ == vmAbsolute)
    {
        VDash_ = V_;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fv::mySemiImplicitSource<Type>::mySemiImplicitSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh),
    volumeMode_(vmAbsolute),
    VDash_(1.0),
	injectionRate_()

{
    read(dict);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Type>
void Foam::fv::mySemiImplicitSource<Type>::addSup
(
	const volScalarField& sw,
    fvMatrix<Type>& eqn,
    const label fieldi
)
{
const GeometricField<Type, fvPatchField, volMesh>& psi = eqn.psi();

typename GeometricField<Type, fvPatchField, volMesh>::Internal Su
    (
        IOobject
        (name_ + fieldNames_[fieldi] + "Su",mesh_.time().timeName(),mesh_,IOobject::NO_READ,IOobject::NO_WRITE),
        mesh_,
        dimensioned<Type> ("zeros",eqn.dimensions()/dimVolume,Zero), false
    );

volScalarField::Internal Sp
    (
        IOobject
        (name_ + fieldNames_[fieldi] + "Sp",mesh_.time().timeName(),mesh_,IOobject::NO_READ,IOobject::NO_WRITE),
        mesh_,
        dimensioned<scalar>("zero",eqn.dimensions()/dimVolume/psi.dimensions(),0.0), false
    );

//Info << "cellSetName " << cellSetName_ << " psi name " << psi().name() << endl;

	// seach for the right time
	float time = mesh().time().value(); //const Time&
	int icount = 0; //Info <<"icount "<< icount << endl;//counter of lines within poption file
	int nlay = int(cellsData_[0]);int ncell_lay = int(cellsData_[1]);
    float tnow = cellsData_[2]*86400; //time in files is in days
	while (time>=tnow)
	{
		icount += cells_.size(); 
		tnow = cellsData_[2+icount*4]*86400; 
	}

	if (icount>0) {icount -= cells_.size();};//Info <<"implicit icount "<< icount << endl;// we have to go to the previous vlaue, the last time<topt

	int mtype =0; //0 for flow, 1 transport, 3 chrmistry like for gwaterFoam
	if (cellSetName_.substr(0,1) == "c" && psi().name().length()<=2) {mtype=1;}
	if ((cellSetName_.substr(0,1) == "c" || cellSetName_.substr(0,1) == "g") && psi().name().length()>2) {mtype=3;} 
	//Info <<"icount "<< icount << " tnow "<< tnow <<" mtype "<<mtype<<endl;
	
	// set Su et Sp
	//to reorder the things (cIndex is calculated once in cellSetOptions)
	std::vector<double> a(cells_.size(),0.);
	std::vector<double> aAdd(cells_.size(),0.); 
	std::vector<double> b(cells_.size(),0.);
	std::vector<double> bAdd(cells_.size(),0.); 
	float dv = 1.;
	if (mtype==0) {dv=86400.;}
	for (int i=0; i<cells_.size(); i++) 
		{
		if (cellSetName_ == "hrch")
			{
			a[i] = cellsData_[2+(icount+i)*4+2]/dv; //fluxes  and conductance in files are in m3/d (17/8 removed cIndex_[i])
			b[i] = cellsData_[2+(icount+i)*4+3]/dv; //fluxes in files are in m3/d  cIndex_[i]
			}
		else
			{
			a[cIndex_[i]] = cellsData_[2+(icount+i)*4+2]/dv; 
			b[cIndex_[i]] = cellsData_[2+(icount+i)*4+3]/dv; 
			}
		}
		
	if (mtype>0)
		{
		for (int i=0; i<cells_.size(); i++)  {
			aAdd[cIndex_[i]] = cellsDataAdd_[2+(icount+i)*4+2]/86400; 
			bAdd[cIndex_[i]] = cellsDataAdd_[2+(icount+i)*4+3]/86400; 
			}
		}
	//Info<<"a et bAdd size "<<aAdd.size()<<" "<<bAdd.size()<<endl;
	//for recharge it needs to be applied to the 1st non-dry layer and divided by cell thikckness (thk imported in cellset)
	int ilay,ilay1,nc;
	double q = 0;double conc = 0;float sumRch=0;float sumDrn=0;float sumGhb=0;float sumWel=0;

	//////////////////////// FLOW   ////////////////////////
	if (mtype==0)
	{
	if (cellSetName_ == "hrch")  // ################## RECHARGE ##############"
		{ //Info << " sw "<< sw[400] <<endl;
		//Info <<"in rch "<<cells_.size()<<" cells"<<endl;
		for (int i=0; i<cells_.size(); i++) 
			{
			for (ilay=0;ilay<nlay;ilay++) { if (sw[cells_[i]-ilay*ncell_lay]>0.999) {break;} ;} //search for the first confined layer from the top
			ilay1 = max(0, ilay-1);
			nc = cells_[i]-ilay1*ncell_lay;
			Su[nc] = injectionRate_[fieldi].first()*a[i]/mesh_.V()[nc];//cells_[cIndex_[i]]-ilay1*ncell_lay];  // injeciton of recharge in the lowest unconfined layer ) /mag(psi()[cells_[i]])
			//if (i%5000==0) {Info<< "hrch cell "<< cells_[i] <<" lay " << ilay1 <<" new c "<<nc <<" h "<<mag(psi()[cells_[i]])<< " sw "<<sw[nc]<<"a "<<a[i]<<" Su "<<Su[nc]<<" vol "<<mesh_.V()[nc]<<" sum "<<sumRch<<endl;}
			sumRch += a[i];
			} // /VDash_ removed because it is in the source data file
		//Info<<"end opt rch, sum "<<sumRch*86400<<endl;
		}
	// for well we need to inject/pump also mass with a rate given in hwel (cellsDataAdd_) ##### Make two options fro transport and chemistry

	else if (cellSetName_ == "hdrn") //################  DRAINS  ########### a drain only discahrge so set to 0 if positive
		{
		//Info <<"in drn "<<cells_.size()<<" cells"<<endl;
		for (int i=0; i<cells_.size(); i++)
			{
			Su[cells_[i]] = injectionRate_[fieldi].first()*a[i]/mesh_.V()[cells_[i]];
			Sp[cells_[i]] = injectionRate_[fieldi].second()*b[i]/mesh_.V()[cells_[i]];
			scalar hprev = mag(psi()[cells_[i]]);
			scalar dff = a[i]+b[i]*hprev;
			//if (cells_[i]%ncell_lay==696) {Info << cells_[i]<< " h prev "<< hprev << " dff "<< dff <<" Su " << Su[cells_[i]]<< " Sp "<< Sp[cells_[i]]<<endl;}
			if (dff<0) {Su[cells_[i]]/=1000;Sp[cells_[i]]/=1000;}  //psi()[cells_[i];10/abs(dff)
			sumDrn += (mag(Su[cells_[i]])-hprev*Sp[cells_[i]])*mesh_.V()[cells_[i]];
			}
		//Info<<"end opt drn, sum "<<sumDrn*86400<<endl;
		}
	else if (cellSetName_ == "hghb")//##################  flow GHB ############"""
		{ 
		//Info <<"in ghb "<<cells_.size()<<" cells"<<endl;
		for (int i=0; i<cells_.size(); i++) 
			{
			Su[cells_[i]] = injectionRate_[fieldi].first()*a[i]/mesh_.V()[cells_[i]];// Info << cells_[i] <<" "<< a[i]<< endl;}
			Sp[cells_[i]] = injectionRate_[fieldi].second()*b[i]/mesh_.V()[cells_[i]]; // /VDash_ removed
			sumGhb += (mag(Su[cells_[i]])-mag(psi()[cells_[i]])*Sp[cells_[i]])*mesh_.V()[cells_[i]];
			} 
		//Info<<"end opt ghb, sum "<<sumGhb*86400<<endl;
		}
	else if (cellSetName_.substr(1,4) == "wel")//################# wells  (for h or p) ###################
		{ 
		//Info <<"in well "<<cells_.size()<<" cells"<<endl;
		for (int i=0; i<cells_.size(); i++) 
			{
			Su[cells_[i]] = injectionRate_[fieldi].first()*a[i]/mesh_.V()[cells_[i]]; Info << cells_[i] <<" "<< a[i]<< endl;
			sumWel += mag(Su[cells_[i]])*mesh_.V()[cells_[i]];
			} 
		//Info<<"end opt wel, sum "<<sumWel*86400<<" time "<<time<<endl;
		}
	}
	
	//////////////////////// TRANSPORT   ////////////////////////
	if (mtype == 1)
	{
	if (cellSetName_ == "cwel") //transport for wells
		{
		for (int i=0; i<cells_.size(); i++) 
			{ 
			q = aAdd[i]/mesh_.V()[cells_[i]]; // dataAdd is the rate (from hwel)
			scalar cloc = mag(psi()[cells_[i]]);
			if (q<0) {Sp[cells_[i]] = -injectionRate_[fieldi].second()*q;} //in Sp term as discharge prop to local concentration
			else {Su[cells_[i]] = injectionRate_[fieldi].first()*q*a[i];} //the injected mass is the water flow x the conc
			Info<<"cell "<<cells_[i]<<" q "<<q<<" cloc "<<cloc<<" Su "<<Su[cells_[i]]<<" Sp "<<Sp[cells_[i]]<<endl;
			}
		}
		
	if (cellSetName_ == "cghb") //transport for ghb, a is concentration
		{
		for (int i=0; i<cells_.size(); i++) 
			{ 
			scalar hprev = mag(psi()[cells_[i]]);
			q = (aAdd[i] + bAdd[i]*hprev)/mesh_.V()[cells_[i]];
			if (q<0) {Sp[cells_[i]] = -injectionRate_[fieldi].second()*q;} //in Sp term as discharge prop to local concentration
			else {Su[cells_[i]] = injectionRate_[fieldi].first()*q*a[i];} 
			}
		}

	if (cellSetName_ == "crch") //transport for recharge, a is concentration
		{ //int flg = 0;
		for (int i=0; i<cells_.size(); i++) 
			{ 
			for (ilay=0;ilay<nlay;ilay++) { if (sw[cells_[i]-ilay*ncell_lay]> 0.999) {break;} ;} //search for the first confined layer from the top
			ilay1 = max(0, ilay-1);
			nc = cells_[i]-ilay1*ncell_lay;
			q = aAdd[i]/mesh_.V()[nc]; // dataAdd is the rate (from hrch)
			Su[nc] = injectionRate_[fieldi].first()*q*a[i];
			//if (a[i]>0) {Info<< i <<" cell "<< nc <<" conc "<< a[i] <<" aAdd "<<aAdd[i]<<" q "<<q << " Su "<< Su[nc] << endl;}
			}
		}
	}
	
	////////////////////// CHEMISTRY  /////////////////////
	if (mtype==3) 
	{
	int ispec = std::stoi(psi().name().substr(2,psi().name().length())); //the species number
	if (cellSetName_ == "crch") //recharge chemistry
		{
		for (int i=0; i<cells_.size(); i++) 
			{ 
			for (ilay=0;ilay<nlay;ilay++) { if (sw[cells_[i]-ilay*ncell_lay] >0.999) {break;} ;} //search for the first confined layer from the top
			ilay1 = max(0, ilay-1);  //-1
			nc = cells_[i]-ilay1*ncell_lay;
			q = aAdd[i]/mesh_.V()[nc]; // dataAdd is the rate (from hrch)			
			conc = cellsSolutions_[(int)a[i]*cellsNspecies_+ispec]; //conc used for injection calc from solution number
			Su[nc] = injectionRate_[fieldi].first()*q*conc;
			//if (i%400==0) {Info<< i <<" spec "<< ispec <<" cell "<< nc <<" conc "<< conc <<" q "<<q << " Su "<< Su[nc] <<" vol "<<mesh_.V()[i]<< endl;}
			}
		}
	else if (cellSetName_ == "cwel") //chemistry for wells
		{
		for (int i=0; i<cells_.size(); i++) 
			{ 
			q = aAdd[i]/mesh_.V()[cells_[i]]; // dataAdd is the rate (from hwel)
			conc = cellsSolutions_[(int)a[i]*cellsNspecies_+ispec]; //conc used for injection calc from solution number
			scalar cloc = mag(psi()[cells_[i]]);
			if (q<0) {Sp[cells_[i]] = -injectionRate_[fieldi].second()*q;} //in Sp term as discharge prop to local concentration
			else {Su[cells_[i]] = injectionRate_[fieldi].first()*q*conc;} //the injected mass is the water flow x the conc
			//Info<<"cwel cell "<<cells_[i]<<" conc inj "<<conc<<" conc local "<<cloc<<" q "<<q<<" Su "<<Su[cells_[i]]<<" Sp "<<Sp[cells_[i]]<<endl;
			}
		}
	else if (cellSetName_ == "gwel") //chemistry for gas injection
		{
		//Info <<"in gwel "<<endl;
		for (int i=0; i<cells_.size(); i++) 
			{ 
			q = aAdd[i]/mesh_.V()[cells_[i]]; // dataAdd is the rate (from pwel)
			conc = cellsGases_[(int)a[i]*cellsGspecies_+ispec]; //conc used for injection calc from solution number
			scalar cloc = mag(psi()[cells_[i]]);
			if (q<0) {Sp[cells_[i]] = -injectionRate_[fieldi].second()*q;} //in Sp term as discharge prop to local concentration
			else {Su[cells_[i]] = injectionRate_[fieldi].first()*q*conc;} //the injected mass is the water flow x the conc
			//Info<<"gwel cell "<<cells_[i]<<" conc inj "<<conc<<" conc local "<<cloc<<" q "<<q<<" Su "<<Su[cells_[i]]<<" Sp "<<Sp[cells_[i]]<<endl;
			}
		}
	else
		{
		//Info<<"size cells "<<cells_.size()<<" size aAdd "<<aAdd.size()<<endl;
		for (int i=0; i<cells_.size(); i++) 
			{ 
			q = aAdd[i]/mesh_.V()[nc]; //Info<<i<<" cell "<<cells_[i]<<" q "<<q<<endl;// dataAdd is the pumping/injection rate
			if (q<0) {
				Sp[cells_[i]] = -injectionRate_[fieldi].second()*q; //this means discharge proportional to local conc (Sp)
				//Info<<i<<" cell "<<cells_[i]<<" ispec "<<ispec<<" q "<<q<<" local conc "<<psi()[cells_[i]]<<" Sp "<<Sp[cells_[i]]<<endl;
				} // in Sp -> mass discharge prop to local conc
			else { // injection given concentration
				//Info<<" a "<<a[i]<<" "<<(int)a[i]<<" nspec "<<cellsNspecies_<<endl;
				conc = cellsSolutions_[(int)a[i]*cellsNspecies_+ispec]; //conc used for injection calc from solution number
				Su[cells_[i]] = injectionRate_[fieldi].first()*q*conc;
				//Info<<i<<" cell "<<cells_[i]<<" ispec "<<ispec<<" q "<<q<<" conc "<<conc<<" Su "<<Su[cells_[i]]<<endl;
				} 
			}
		}
	}
		
	//eqn += Su - fvm::SuSp(-Sp, psi);
	//Info << "in opt Su "<<Su<<endl;
	int sgn1 = 1;
	if(mag(injectionRate_[fieldi].second())==1 && time>0) {sgn1 = -1;}  //&& cellSetName_.substr(0,1) == "h" 
	eqn += sgn1*Su+ sgn1*fvm::SuSp(Sp, psi); // with su<0 and sp>0 in opfwriter
	
	//const objectRegistry &db ;
	//const Time &runTime= db().time();
	//Info<< mesh().time().deltaT()<<endl;
	
	/* // straonge : if written here, they give wrong numbers for rch, ghb, drn
	if (mesh().time().write() && cellSetName_ == "hrch") {Info<<"time "<<time<<" tot rch "<<sumRch*86400<<endl;}
	if (mesh().time().write() && cellSetName_ == "hghb") {Info<<"time "<<time<<" tot ghb "<<sumGhb*86400<<endl;}
	if (mesh().time().write() && cellSetName_ == "hdrn") {Info<<"time "<<time<<" tot drn "<<sumDrn*86400<<endl;}
	if (mesh().time().write() && cellSetName_ == "hwel") {Info<<"time "<<time<<" tot wel "<<sumWel*86400<<endl;}
	*/
	
	/*
	if (mtype == 0 && time<=0) {eqn += Su + fvm::SuSp(-Sp, psi); }//sgn1*Su + sgn1*fvm... Info << "fvopt eqn set " <<endl;
	else if (mtype == 0 && time>0) {eqn += Su - fvm::SuSp(Sp, psi); }//sgn1*Su + sgn1*fvm... Info << "fvopt eqn set " <<endl;
	else {eqn += Su + fvm::SuSp(Sp, psi); }  // with su>0 and sp>0 (all negative doe snot work)
	*/
	/*
	Info<<"mtype "<<mtype<<endl;
	if (mtype == 0) {eqn += Su - fvm::SuSp(Sp, psi); }//sgn1*Su + sgn1*fvm... Info << "fvopt eqn set " <<endl;
	else {eqn += Su + fvm::SuSp(Sp, psi); }  // with su>0 and sp>0 (all negative doe snot work)
	*/
	}


// ************************************************************************* //
