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

#include "cellSetOption.H"
#include "volFields.H"

#include <iostream> // this nad below added for reading
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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
	
	namespace fv
    {
        defineTypeNameAndDebug(cellSetOption, 0);
    }

    template<> const char* NamedEnum
    <
        fv::cellSetOption::selectionModeType,
        4
        >::names[] =
    {
        "points",
        "cellSet",
        "cellZone",
        "all"
    };

    const NamedEnum<fv::cellSetOption::selectionModeType, 4>
        fv::cellSetOption::selectionModeTypeNames_;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::cellSetOption::setSelection(const dictionary& dict)
{
    switch (selectionMode_)
    {
        case smPoints:
        {
            dict.lookup("points") >> points_;
            break;
        }
        case smCellSet:
        {
            dict.lookup("cellSet") >> cellSetName_;
            break;
        }
        case smCellZone:
        {
            dict.lookup("cellZone") >> cellSetName_;
            break;
        }
        case smAll:
        {
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown selectionMode "
                << selectionModeTypeNames_[selectionMode_]
                << ". Valid selectionMode types are" << selectionModeTypeNames_
                << exit(FatalError);
        }
    }
}


void Foam::fv::cellSetOption::setCellSet()
{
    switch (selectionMode_)
    {
        case smPoints:
        {
            Info<< indent << "- selecting cells using points" << endl;

            labelHashSet selectedCells;

            forAll(points_, i)
            {
                label celli = mesh_.findCell(points_[i]);
                if (celli >= 0)
                {
                    selectedCells.insert(celli);
                }

                label globalCelli = returnReduce(celli, maxOp<label>());
                if (globalCelli < 0)
                {
                    WarningInFunction
                        << "Unable to find owner cell for point " << points_[i]
                        << endl;
                }

            }

            cells_ = selectedCells.toc();

            break;
        }
        case smCellSet:
        {
            Info<< indent
                << "- selecting cells using cellSet " << cellSetName_ << endl;

            cellSet selectedCells(mesh_, cellSetName_);
            cells_ = selectedCells.toc();

            break;
        }
        case smCellZone:
        {
            Info<< indent
                << "- selecting cells using cellZone " << cellSetName_ << endl;

            label zoneID = mesh_.cellZones().findZoneID(cellSetName_);
            if (zoneID == -1)
            {
                FatalErrorInFunction
                    << "Cannot find cellZone " << cellSetName_ << endl
                    << "Valid cellZones are " << mesh_.cellZones().names()
                    << exit(FatalError);
            }
            cells_ = mesh_.cellZones()[zoneID];

            break;
        }
        case smAll:
        {
            Info<< indent << "- selecting all cells" << endl;
            cells_ = identity(mesh_.nCells());

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown selectionMode "
                << selectionModeTypeNames_[selectionMode_]
                << ". Valid selectionMode types are" << selectionModeTypeNames_
                << exit(FatalError);
        }
    }

    // Set volume information (not used)
    V_ = 0.0;
    forAll(cells_, i)
    {
        V_ += mesh_.V()[cells_[i]];
    }
    reduce(V_, sumOp<scalar>());
	
	// get data, thickness, solutions
	std::string cur_dir = get_current_dir(); //Info << cur_dir <<endl;
	std::ifstream inputFcdata{cur_dir+"/constant/options/" + cellSetName_ };
	cellsData_ = {std::istream_iterator<float>{inputFcdata}, {}};
	//for (int i=0;i<cellsData_.size();i++) {Info<<"fv data "<<cellsData_[i]<<endl;}
	
	// if concentration (chem or transport) we need additionnal data (the h or p ones)
	
	if (cellSetName_.substr(0,1) == "c" || cellSetName_.substr(0,1) == "g") 
		{
		std::ifstream inputFcdata1{cur_dir+"/constant/options/h" + cellSetName_.substr(1,3) };
		cellsDataAdd_ = {std::istream_iterator<float>{inputFcdata1}, {}};
		}
		
	std::vector<char> str = {'0','1','2','3','4','5','6','7','8','9'}; // (2)
	char lastc = name_.back();
	bool chem = (std::find(str.begin(), str.end(), lastc) != str.end());
	//Info<<" in cellOpt name "<<name_<<" "<<lastc<<" chem "<<chem<<endl;
	if (chem) //chemistry
		{
		std::ifstream inputSolutions{cur_dir+"/constant/options/solutions"};
		cellsSolutions_ = {std::istream_iterator<float>{inputSolutions}, {}};
		std::ifstream inputGases{cur_dir+"/constant/options/gases"};
		cellsGases_ = {std::istream_iterator<float>{inputGases}, {}};

		std::ifstream inputPhq{cur_dir+"/phqfoam.txt"};
		a = {std::istream_iterator<int>{inputPhq}, {}};
		cellsNspecies_ = a[1];//Info<<" in cellOpt nspec "<<cellsNspecies_<<endl;
		cellsGspecies_ = a[2];
		}
	
	//std::ifstream inputThk{cur_dir+"/constant/options/thk"};
	//cellsThk_ = {std::istream_iterator<float>{inputThk}, {}};
	
	// create the index of cellSetData in the cells_ (they might have been reordered by openFoam)
 	float tnow = cellsData_[2];
	int icd; 
	cIndex_.resize(cells_.size(),0);
	for (int i=0; i<cells_.size();i++)  // reads the first ncells lines
		{
		icd = cellsData_[2+i*4+1];  
		auto iter = std::find(cells_.begin(), cells_.end(), icd);
		int a = {std::distance(cells_.begin(), iter)};  //Info << i << " icd "<< icd << " cll " << cells_[i] << " indx " << a << endl; // cell number in cellsData //index of cellsData in cells_
		cIndex_[i] = a;
		}
	//Info << " cdata in celSet " << cellsData_[6] << endl;
    //for (int i=0; i<cells_.size();i++) {Info << "cindex "<<cIndex_[i]<< endl;}
	Info<< indent
        << "- selected " << returnReduce(cells_.size(), sumOp<label>())
        << " cell(s) with volume " << V_ << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::cellSetOption::cellSetOption
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh),
    timeStart_(-1.0),
    duration_(0.0),
    selectionMode_
    (
        selectionModeTypeNames_.read(coeffs_.lookup("selectionMode"))
    ),
    cellSetName_("none"),
    V_(0.0)
{
    Info<< incrIndent;
    read(dict);
    setSelection(coeffs_);
    setCellSet();
    Info<< decrIndent;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::cellSetOption::~cellSetOption()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::cellSetOption::isActive()
{
    if (option::isActive() && inTimeLimits(mesh_.time().value()))
    {
        // Update the cell set if the mesh is changing
        if (mesh_.changing())
        {
            setCellSet();
        }

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
