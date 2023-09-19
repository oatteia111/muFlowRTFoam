thetag = max(eps*(1 - sw),SMALL);
if (diffusionEqn==0) {Deffg = Dg0 * tensor::I;}
if (diffusionEqn==1) {Deffg = Dg0 * pow(thetag,3.333) / pow(eps,2.) * tensor::I;}
if (diffusionEqn==2) {Deffg = Dg0 * pow(thetag,2.333) / pow(eps,2.) * tensor::I;}
//fDe_T = runTime.controlDict().lookupOrDefault("fDe_T",1);
std::vector<double> resu ;
resu.resize(nxyz);
resu = myFunc(fDe_T.type,T,fDe_T.parms);//Info<<"ftype "<<fDe_T.type<<" De parms "<<fDe_T.parms[0]<<" "<<fDe_T.parms[1]<<endl;
for (i=0;i<nxyz;i++) {Deffg[i] *= resu[i];}
Info<<"Deff 0 "<<Deffg[0]<<" resu 0 "<<resu[0]<<endl;
int ig;
// -----------------  get internal conditions  ----------------------
//get the gas mixture well data
icount = 0;
if (cellsGwel.size()>0)
	{
	if (time>=tGwel) {
		output = getCbuffer("gwel",itGwel,cellsGwel.size()); 
		tGwel = output.t; cellsGwelData=output.d; itGwel+=1;
		}
	for (int i=0; i<cellsGwel.size(); i++) {condGwel[i] = cellsGwelData[i*4+2];}
	}
// get fixed gas mixture 
icount = 0;
if (cellsGfix.size()>0)
	{
	if (time>=tGfix) {
		output = getCbuffer("gfix",itGfix,cellsGfix.size()); 
		tGfix = output.t; cellsGfixData=output.d; itGfix+=1;
		}
	for (int i=0; i<cellsGfix.size(); i++) {condGfix[i] = cellsGfixData[i*4+2];Info<<" gfx "<<i<<" "<<condGfix[i]<<endl;} // discharge must be in s-1
	}

for (ig=0; ig<ph_gcomp;ig++)
    {
		//Info << "in CgiEqn "<< i << endl;
		volScalarField& Cgi = Cg[ig]();
		word nameCi ("Cg" + name(ig));
		
		// --------------------- internal conditions ------------------------
		if (cellsGwel.size()>0)  // gas in well
		{
			//int nlay = int(cellsHwData[0]);int ncell_lay = int(cellsHwData[1]);
			float q=0;
			// set the values in the Su Sp matrix
			CiSu*=0; CiSp*=0; //previous cells with values must be destroyed
			for (int i=0; i<cellsGwel.size(); i++) 
				{ 
				//q = condHw[i]/mesh.V()[cellsHw[i]]; // dataAdd is the rate (from hwel here without dimension)
				conc = gas_conc[static_cast<int>(condGwel[i])*ph_gcomp+ig];
				if (q<0) {CiSp[cellsPwel[i]] = -qPwel[i];} //in Sp term as discharge prop to local concentration
				else {CiSu[cellsPwel[i]] = -qPwel[i]*conc;} 
				//Info<<"cell "<<cellsPwel[i]<<" q "<<qPwel[i]<<" Su "<<CiSu[cellsCwel[i]]<<" Sp "<<CiSp[cellsCwel[i]]<<endl;
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
		for (int i=0; i<cellsGfix.size(); i++) {CgiEqn.setReference(cellsGfix[i],gas_conc[static_cast<int>(condGfix[i])*ph_gcomp+ig]);}
		CgiEqn.solve();
    }
