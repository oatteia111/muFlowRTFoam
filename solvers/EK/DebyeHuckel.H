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
dimensionedScalar souE("SouE", psiE.dimensions()/dimArea, 0.);       
dimensionedScalar souI("SouI", dimless/dimArea, 0.);
dimensionedScalar TK("temp", dimensionSet(0,0,0,1,0,0,0), 0);   

for (i=1; i<ph_nspc;i++) {cw0[i] = mean(Cw[i]();} // OA added, trial to calc reference conc

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
   scalar res=GREAT; 
   scalar iter=0;  
   while (res > 1e-8 && iter < 4)
		{ 
		fvScalarMatrix phiEEqn
		(
			fvm::laplacian(phiE)
		);
		phiEEqn.relax();
		res=phiEEqn.solve().initialResidual();
		iter++;
		}
      
   //- Equation for the intrinsic potential
   res=GREAT;
   iter=0;  
   forAll (Cwi, i)
	 {
	   souE += z[i]*cw0[i]*FK/(relPerm*epsilonK);
	   souI += eK*z[i]*z[i]*cw0[i]*FK/(relPerm*epsilonK*kbK*TK);
	 }

	while (res > 1e-8 && iter < 4)
	{ 
		fvScalarMatrix psiEqn // eq.3.42 in rheoTool (ok)
		(	  
			 fvm::laplacian(psiE) == fvm::Sp( souI, psiE) - souE
		); 
		psiEqn.relax();
		res=psiEqn.solve().initialResidual();
		iter++;
	} 

}

// ************************************************************************* //
