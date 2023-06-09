thetag = max(eps*(1 - sw),SMALL);
if (diffusionEqn==0) {Deffg = Dg0 * pow(thetag,3.333) / pow(eps,2.) * tensor::I;}
if (diffusionEqn==1) {Deffg = Dg0 * pow(thetag,2.333) / pow(eps,2.) * tensor::I;}
Info<<"sw "<<sw[2]<<" thetag "<<thetag[2]<<" Deffg 2 "<<Deffg[2]<<endl;
int ig;
// -----------------  get internal conditions  ----------------------
//get the gas mixture well data
if (cellsGwl.size()>0)
	{
	float tnow = cellsHwData[2]*day; //time in files is in days
	while (time>=tnow) {icount += cellsGwl.size(); tnow = cellsGwlData[2+icount*4]*day; }
	if (icount>0) {icount -= cellsHw.size();};
	for (int i=0; i<cellsHw.size(); i++) {condGwl[i] = cellsGwlData[2+(icount+i)*4+2];}
	}
// get fixed gas mixture 
if (cellsGfx.size()>0)
	{
	float tnow = cellsGfxData[2]*day; //time in files is in days
	while (time>=tnow) {icount += cellsGfx.size(); tnow = cellsGfxData[2+icount*4]*day; }
	if (icount>0) {icount -= cellsGfx.size();};
	for (int i=0; i<cellsGfx.size(); i++) {condGfx[i] = cellsGfxData[2+(icount+i)*4+2];Info<<" gfx "<<i<<" "<<condGfx[i]<<endl;} // discharge must be in s-1
	}

for (ig=0; ig<ph_gcomp;ig++)
    {
		//Info << "in CgiEqn "<< i << endl;
		volScalarField& Cgi = Cg[ig]();
		word nameCi ("Cg" + name(ig));
		
		// --------------------- internal conditions ------------------------
		if (cellsGwl.size()>0)  // gas in well
		{
			//int nlay = int(cellsHwData[0]);int ncell_lay = int(cellsHwData[1]);
			float q=0;
			// set the values in the Su Sp matrix
			CiSu*=0; CiSp*=0; //previous cells with values must be destroyed
			for (int i=0; i<cellsGwl.size(); i++) 
				{ 
				//q = condHw[i]/mesh.V()[cellsHw[i]]; // dataAdd is the rate (from hwel here without dimension)
				conc = solu_conc[(int)condGwl[i]*ph_gcomp+ig];
				if (q<0) {CiSp[cellsHw[i]] = -qHw[i];} //in Sp term as discharge prop to local concentration
				else {CiSu[cellsHw[i]] = -qHw[i]*conc;} 
				Info<<"cell "<<cellsHw[i]<<" q "<<q<<" Su "<<CiSu[cellsCw[i]]<<" Sp "<<CiSp[cellsCw[i]]<<endl;
				}
		}

		fvScalarMatrix CgiEqn
		(		  
            thetag * fvm::ddt(Cgi) +  Cgi * fvc::ddt(thetag)
			+ fvm::div(phig, Cgi, "div(phig,Cgi)") // used only in multiphase
            - fvm::laplacian(Deffg, Cgi, "laplacian(Deffg,Cgi)")
			+ CiSu + fvm::SuSp(CiSp, Cgi)
		);  //			== fvOptions(sw,Cgi)

		//fvOptions.constrain(CgiEqn);
		for (int i=0; i<cellsGfx.size(); i++) {CgiEqn.setReference(cellsGfx[i],gas_conc[(int)condGfx[i]*ph_gcomp+ig]);}
		CgiEqn.solve();
    }
