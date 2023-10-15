thetag = max(eps*(1 - sw),SMALL);
if (diffusionEqn==0) {Deffg = Dg0 * tensor::I;}
if (diffusionEqn==1) {Deffg = Dg0 * pow(thetag,3.333) / pow(eps,2.) * tensor::I;}
if (diffusionEqn==2) {Deffg = Dg0 * pow(thetag,2.333) / pow(eps,2.) * tensor::I;}
Info<<"Dg0 "<<Dg0<<" sw "<<sw[0]<<" " <<sw[1]<<" thetag "<<thetag[0]<<" "<<thetag[1]<<endl;
std::vector<double> resu ;
resu.resize(nxyz);
resu = myFunc(fDe_T.type,T,fDe_T.fparms);
//Info<<"ftype "<<fDe_T.type<<" De parms "<<fDe_T.fparms[0]<<" "<<fDe_T.fparms[1]<<endl;
for (int i=0;i<nxyz;i++) {Deffg[i] *= resu[i];}
Info<<"Deff 0 "<<Deffg[0]<<" resu 0 "<<resu[0]<<endl;
int ig;
// -----------------  get internal conditions  ----------------------
//get the gas mixture well data
icount = 0;size_t i;
if (cellsGwel.size()>0)
	{
	if (time>=tGwel) {
		output = getCbuffer("gwel",itGwel,cellsGwel.size()); 
		tGwel = output.t; cellsGwelData=output.d; itGwel +=1;
		}
	for (i=0; i<cellsGwel.size(); i++) {condGwel[i] = static_cast<int>(cellsGwelData[i*4+2]);}
	}
// get fixed gas mixture 
icount = 0;
if (cellsGfix.size()>0)
	{
	if (time>=tGfix) {
		output = getCbuffer("gfix",itGfix,cellsGfix.size()); 
		tGfix = output.t; cellsGfixData=output.d; itGfix +=1;
		}
	for (i=0; i<cellsGfix.size(); i++) {condGfix[i] = static_cast<int>(cellsGfixData[i*4+2]);}
	}

for (ig=0; ig<ph_gcomp;ig++)
    {
		volScalarField& Cgi = Cg[ig]();
		word nameCi ("Cg" + name(ig));
		
		// --------------------- internal conditions ------------------------
		if (cellsGwel.size()>0)  // gas in well
		{
			float q=0;
			// set the values in the Su Sp matrix
			CiSu*=0; CiSp*=0; //previous cells with values must be destroyed
			for (i=0; i<cellsGwel.size(); i++) 
				{ 
				i1=cellsGwel[i];v=mesh.V()[i1];
				if (qPwel[i]<0) {CiSp[i1] = -qPwel[i]/v;} //in Sp term as discharge prop to local concentration
				else { CiSu[i1] = -qPwel[i]/v*gas_conc[condGwel[i]*ph_gcomp+ig]; } 
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
		for (i=0; i<cellsGfix.size(); i++) {CgiEqn.setReference(cellsGfix[i],gas_conc[condGfix[i]*ph_gcomp+ig]);}
		CgiEqn.solve();
    }
