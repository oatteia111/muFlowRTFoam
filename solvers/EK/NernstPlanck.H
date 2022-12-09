/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

dimensionedScalar cum("0", dimensionSet(0, -3, 0, 0, 1, 0, 0), 0);   
volScalarField rhoE( psiE * dimensionedScalar("norm", epsilonK.dimensions()/dimArea, 0.) );
volScalarField souE( psiE * dimensionedScalar("norm1",dimless/dimArea,0.));          
dimensionedScalar souI("SouI", dimless/dimArea, 0.);

/* for (i=1; i<ph_nspc;i++)
	{        
        cum += species_[specI].DebyeLengthP(relPerm, TK); 
    }     
    dimensionedScalar DebL ( sqrt (epsilonK* relPerm* kbK*TK / (cum * FK * eK) ) );   
    Info << "Debye length: " << DebL.value() << " m." << nl << endl; 
*/   

/*for (i=1; i<ph_nspc;i++)
    {
      rhoE  += ( z[i]*Cw[i]()*FK * exp(1.- eK*z[i]*psiE/(kbK*TK) ) );
    }

    // Fe? = -rhoE * ( fvc::grad(phiE+psi) - extraE) //seems in fact to be used only in more complex models
*/
//- Equation for the external potential (loop for the case
//  of non-orthogonal grids) 
int nIter=6;
float minRes=1e-12;
Info<<"solving external "<<endl;
 for (int j=0; j<nIter; j++) 
{
   scalar res=GREAT; 
   scalar iter=0;  
   while (res > minRes && iter < nIter)
		{ 
		fvScalarMatrix psiEEqn
		(
			fvm::laplacian(psiE)
		);
		psiEEqn.relax();
		res=psiEEqn.solve().initialResidual();
		iter++;
		}
    
  //- Equation for the intrinsic potential  
Info<<"solving intrinsic "<<endl;
  res=GREAT;
  iter=0; souE *= 0;
  forAll (Cw, i)
  {
   souE +=  -z[i]*Cw[i]()*1000*FK/(relPerm*epsilonK); //*1000 as in gw3 Cwi is in mol/L
   Info<<" i "<<i<<" z[i] "<<z[i]<<" Cwi "<<Cw[i]()[20]<<" souE "<<souE[20]<<endl;
  }

  while (res > minRes && iter < nIter)
  { 
    fvScalarMatrix psiIEqn
    (	   
      fvm::laplacian(psiI) == souE		           		    
    );
	     
    psiIEqn.relax();
    res=psiIEqn.solve().initialResidual();

    iter++;
  } 
         
  //- Nernst-Planck equation for each ionic specie  
    
  //  surfaceScalarField eMigFluxp( fvc::snGrad(phiE_+psi_) * phi().mesh().magSf() ); // Compute once, outside the loop 
Info<<"solving Cwi "<<endl;
        
	forAll (Cw, i)
	{         
		res=GREAT;
		iter=0;
		D25.value() = diff25[i];							
		dimensionedScalar cf(D25 * eK * z[i] / (kbK*TK));  
		volScalarField& Cwi = Cw[i]();
		word nameCi ("Cw" + name(i));
		
			   
		while (res > minRes && iter < nIter)
		{      
			fvScalarMatrix ciEqn 
			(
			 1000*fvm::ddt(Cwi)   == //+ fvm::div(phiw(), Cwi, "div(phiw,Cwi)")  ==  //no water movement at first
			 fvm::laplacian(1000*D25, Cwi, "laplacian(D,Cwi)") 	 
			 //    + fvm::div(eMigFluxp*cf, ci, "div(eMigFlux,ci)") 
			+ fvc::laplacian(1000*Cwi*cf, psiE+psiI, "laplacian(elecM)")	 
			);      

			if (psiE().mesh().relaxEquation("Cwi")) ciEqn.relax(psiE().mesh().equationRelaxationFactor("Cwi"));

			res=ciEqn.solve(psiE().mesh().solverDict("Cwi")).initialResidual();
			//Cwi = max(Cwi,0.); //Foam::max( dimensionedScalar("lowerLimit",Cwi.dimensions(), 0.), Cwi );
			Cwi.correctBoundaryConditions();
			iter++;
		} 
	}
              
}

// ************************************************************************* //
