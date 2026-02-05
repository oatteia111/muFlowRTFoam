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
int nIter=6;
float minRes=1e-9;
scalar res=GREAT;
dimensionedScalar JSMALL("",dimensionSet(1,-2,-1,0,0,0,0),1e-36);

//we iterate with picard to get the residual small
//we sum up all positive species to get JE+ and same for JE-, the objective is to find factors fJE+ and FJE- that makes fJE+*JE+=fJE-*abs(JE-)
//we cannot do it at once because there is diffusion mixed with EK in J

iterPicard = 0; resPicard = 1000.;
JE1 = JE;
scalar deltaTFact = 1.;
bcCwi.resize(ncell,0);

// get fixed solutions at h fix places
if (cellsShfix.size()>0)
	{
	if (time>=tShfix) {
		output = getCbuffer("shfix",itShfix,cellsShfix.size()); 
		tShfix = output.t; cellsShfixData=output.d; itShfix +=1;
		tnext = min(tnext,tShfix);
		}
	for (i=0; i<cellsShfix.size(); i++) 
		{condShfix[i] = static_cast<int>(cellsShfixData[i*4+2]);i1=cellsShfix[i];bcCwi[i1]=1;}
	Info<<" tShfix "<<tShfix<<" cond 0 "<<condShfix[0]<<endl;
	}
// get fixed solutions 
if (cellsSfix.size()>0)
	{
	if (time>=tSfix) {
		output = getCbuffer("sfix",itSfix,cellsSfix.size()); 
		tSfix = output.t; cellsSfixData=output.d; itSfix +=1;
		tnext = min(tnext,tSfix);
		}
	for (i=0; i<cellsSfix.size(); i++) 
		{condSfix[i] = static_cast<int>(cellsSfixData[i*4+2]);i1=cellsSfix[i];bcCwi[i1]=1;}
	}
//get the solution well data
if (cellsSwel.size()>0)
	{
	if (time>=tSwel) {
		output = getCbuffer("swel",itSwel,cellsSwel.size()); 
		tSwel = output.t; cellsSwelData=output.d; itSwel +=1;
		tnext = min(tnext,tSwel);
		}
	for (i=0; i<cellsSwel.size(); i++) 
		{condSwel[i] = static_cast<int>(cellsSwelData[i*4+2]);i1=cellsSwel[i];bcCwi[i1]=1;}
	Info<<" tSwell "<<tSwel<<" cond 0 "<<condSwel[0]<<endl;
	}

while (resPicard > tolPicard)
	{
	iterPicard++;
	JE = JE_tmp; 
	Info<<"solving E "<<endl;
	// external potential equation 	
	fvScalarMatrix psiEEqn ( fvm::laplacian(psiE) );
	psiEEqn.relax();
	//psiEEqn.solve();
	res = psiEEqn.solve().initialResidual();
	Info<<"residu psiE "<<res<<endl;
	
    /*
	souE *= 0;
	forAll (Cw, i)
    {
    souE +=  -invMolweight*freak.z[i]*Cw[i]()*FK/(relPerm*epsilonK); //invMol because C are in kg/m3 not moles/m3
    }
	*/
    //printing some results
	for (i=0;i<freak.nspc;i++) {if (freak.z[i]!=0) {Info<<freak.z[i]<<" ";} } Info<<endl;
    for (j=0;j<nxyz;j++) {
		for (i=0;i<freak.nspc;i++) {if (freak.z[i]!=0) {Info<<Cw[i]()[j]<<" ";} } Info<<souE[j]<<" "<<psiE[j]<<endl;
	}

	/*
	fvScalarMatrix psiIEqn
	(	   
	  fvm::laplacian(psiE) == souE		           		    
	); 
	psiIEqn.relax();
	psiIEqn.solve();
	res += psiIEqn.solve().initialResidual();
    */
    //- Nernst-Planck equation for each ionic specie  
    //  surfaceScalarField eMigFluxp( fvc::snGrad(phiE_+psi_) * phi().mesh().magSf() ); // Compute once, outside the loop 
    Info<<"solving Cwi "<<endl;//+ adding all fluxes in JEpos and JEneg
        
	for (int ic=1; ic<freak.nspc;ic++)		
	{   
		Info<<"species nb "<<i<<endl;
		// --------------------- internal conditions ------------------------
		if (cellsHwel.size()>0)  // solution in well (even if no Swl data we ned to pump the conc)
		{
			for (i=0; i<cellsHwel.size(); i++) 
				{ 
				i1=cellsHwel[i];v=mesh.V()[i1];unity[i1]=0;
				if (qHwel[i]<0) {CiSp[i1] = -qHwel[i]*(1+1e-6)/v;} //in Sp term as discharge prop to local concentration
				else { CiSu[i1] = -qHwel[i]/v*solu_species[condSwel[i]*freak.nspc+ic]; } 
				}
			Info<< ic << " Csu wel "<<solu_species[condSwel[0]*freak.nspc+ic]<<" "<<CiSu[i1]<<endl;
		} //Info<<"Cweli set"<<endl;
		if (cellsHfix.size()>0)
			{
			for (i=0; i<cellsHfix.size(); i++) 
				{ 
				i1=cellsHfix[i];v=mesh.V()[i1];//
				if (qHfix[i]<0) {CiSp[i1] = -qHfix[i]/v;} // outflow discharge prop to local concentration
				else {CiSu[i1] = -qHfix[i]/v*solu_species[condShfix[i]*freak.nspc+ic];}
				}
			}
		//------------------  end internal   -----------------------
		if (freak.z[i]>0) {fJE1=1.+fJE;} else {fJE1=1.-fJE;}//a unitless volscalarfield
		D25.value() = freak.diff25[i];	// Info<<"diff "<<D25<<endl;					
		cf = fJE1 * D25 * freak.z[i] *FK / (Rgas*TK) * eps; Info<<"cf "<<cf[0]<<endl;	
		volScalarField& Cwi = Cw[i]();
		word nameCi ("Cw" + name(i));
		
		fvScalarMatrix EciEqn 
		(
		 fvm::ddt(Cwi) //+ fvm::div(phiw, ECwi, "div(phiw,C)")  // water movement removed up to now
		 - fvm::laplacian(fJE1 * D25, Cwi, "laplacian(D,ECwi)") 	//!!!! we need to add tortuosity 
		 - fvc::laplacian(cf * Cwi , psiE, "laplacian(elecM)")	 
		); 
		
		for (i=0; i<cellsSfix.size(); i++) 
			{
			EciEqn.setReference(cellsSfix[i],solu_species[condSfix[i]*freak.nspc+ic]);
			}

		JE = mag(fJE1 * D25 * fvc::grad(Cwi) + cf*Cwi * fvc::grad(psiE) ); // flux of concentrations
		Info<<"JE 0-2 "<<JE[0]<<" "<<JE[1]<<" "<<JE[2]<<endl;
		if (freak.z[i]>0) {JEpos += mag(JE);} else {JEneg += mag(JE);}
		//if (psiE().mesh().relaxEquation("Cwi")) ciEqn.relax(psiE().mesh().equationRelaxationFactor("Cwi"));

		res += EciEqn.solve().initialResidual();
		//Cwi = max(Cwi,0.); //Foam::max( dimensionedScalar("lowerLimit",Cwi.dimensions(), 0.), Cwi );
		Cwi.correctBoundaryConditions();
	}
	Info<<"residual "<<res<<endl;
	Info<<"pos neg "; for (i=0;i<4;i++) {Info<<JEpos[i]<<" "<<JEneg[i]<<" ";} Info<<endl;
	fJE = JEpos*2./(JEpos+JEneg+JSMALL) -1.;
	
	resPicard = gMax((mag(JE-JE_tmp))->internalField());
	JE_tmp = JE;
	if (iterPicard == maxIterPicard)
		{ runTime.setTime(oldTime,tstep); break; }
	}  // end picard iter
	Info << "Picard nb iterations : "<<iterPicard<<endl;
	if (iterPicard == maxIterPicard) {deltaTFact = 0.2;JE_tmp=JE1;JE=JE1;rewind=1;}
	else if (iterPicard <= nIterStability) {deltaTFact = 1.2;}
	newDeltaT = min(deltaTFact*runTime.deltaTValue(),maxDeltaT);

// ************************************************************************* //
