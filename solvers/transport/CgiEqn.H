thetag = max(eps*(sw_max - sw),SMALL);
if (diffusionEqn==0) {Deffg = vDg0 * tensor::I;}
if (diffusionEqn==1) {Deffg = (vDg0*pow(thetag,3.333)+Dw0*pow(thetaw,3.333))/pow(eps,2.) * tensor::I;}
if (diffusionEqn==2) {Deffg = (vDg0*pow(thetag,2.333)+Dw0*pow(thetaw,2.333))/pow(eps,2.) * tensor::I;}
if (diffusionEqn==3) {Deffg = Dg0 * 0.058*max(atan((thetag-0.21)*27)+1.36,SMALL) * tensor::I;} // very rough, a test, valid only for cm2 and min
Deffg *= factDg;
Info<<"Dg0 "<<Dg0<<" sw "<<sw[1]<<" " <<sw[ncell-2]<<" thetag "<<thetag[0]<<" "<<thetag[1]<<" Deffg "<<Deffg[1]<<" "<<Deffg[ncell-2]<<endl;
//Info<<"ftype "<<fDe_T.type<<" De parms "<<fDe_T.fparms[0]<<" "<<fDe_T.fparms[1]<<endl;
//for (int i=0;i<nxyz;i++) {Deffg[i] *= resu[i];}
//Info<<"Deff 0 "<<Deffg[0]<<" resu 0 "<<resu[0]<<endl;
int ig;
// -----------------  get internal conditions  ----------------------
//get the gas mixture well data
size_t i;
if (cellsGwel.size()>0)
	{
	if (time>tGwel) {
		output = getCbuffer("gwel",itGwel,cellsGwel.size()); 
		tGwel = output.t; cellsGwelData=output.d; itGwel +=1;
		}
	for (i=0; i<cellsGwel.size(); i++) {condGwel[i] = static_cast<int>(cellsGwelData[i*4+2]);}
	}
// get fixed gas mixture 
if (cellsGfix.size()>0)
	{
	if (time>tGfix) {
		output = getCbuffer("gfix",itGfix,cellsGfix.size()); 
		tGfix = output.t; cellsGfixData=output.d; itGfix +=1;
		}
	for (i=0; i<cellsGfix.size(); i++) {condGfix[i] = static_cast<int>(cellsGfixData[i*4+2]);}
	}

for (ig=0; ig<ph_gcomp;ig++)
    {
		volScalarField& Cgi = Cg[ig]();
		word nameCi ("Cg" + name(ig));
		GiSu *=0;GiSp *=0;
		
		// --------------------- internal conditions ------------------------
		if (cellsGwel.size()>0)  // gas in well, only for 2phase
		{
			float q=0;
			// set the values in the Su Sp matrix
			GiSu*=0; GiSp*=0; //previous cells with values must be destroyed
			for (i=0; i<cellsGwel.size(); i++) 
				{ 
				i1=cellsGwel[i];v=mesh.V()[i1];
				if (flowType==4) {q=qPwel[i];} else {q=qHwel[i];}
				if (q<0) {GiSp[i1] = -q/v;} //in Sp term as discharge prop to local concentration
				else { GiSu[i1] = -q/v*gas_conc[condGwel[i]*ph_gcomp+ig]/Vmol[i1];} // ; } 
				Info<<"g wel q conc "<<q<<" "<<gas_conc[condGwel[i]*ph_gcomp+ig]<<endl;
				}
		}
		//std::cout<<"before lbda \n";
		lbdag.value() = plugCgi.getLbdag(ig); // this is the plugin for degradation,needs to be here (compg specific)
		Info<<"lbdag "<<ig <<" "<<lbdag<<endl;
		
		if (plugCgi.flag) {
			#include "plugins/plugin_Cgiupdate.H"
			}
			
		fvScalarMatrix CgiEqn
		(		  
            thetag * fvm::ddt(Cgi) +  Cgi * fvc::ddt(thetag)
			+ fvm::div(phig, Cgi, "div(phig,Cgi)") // used only in multiphase
            - fvm::laplacian(Deffg, Cgi, "laplacian(Deffg,Cgi)")
			+ thetag * fvm::Sp(lbdag,Cgi) //decay
			== - GiSu - fvm::SuSp(GiSp, Cgi)
		); 

		for (i=0; i<cellsGfix.size(); i++) {CgiEqn.setReference(cellsGfix[i],gas_conc[condGfix[i]*ph_gcomp+ig]/Vmol[cellsGfix[i]]);}
		CgiEqn.solve();
		Info<<" after solve cell 0 "<<ig<<" " <<Cgi[0]<<endl;
    }
