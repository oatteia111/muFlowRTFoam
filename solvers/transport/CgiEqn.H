thetag = max(eps*(1 - sw),SMALL);
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
//for (j=0;j<5;j++) { for (i=0;i<ph_gcomp;i++) {Info<<Cg[i]()[j]<<" ";} Info<<endl; }
Ggrd.resize(cellsGfix.size()*ph_gcomp,0);
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
//Info<<" gwel "<<endl;
// get fixed gas mixture 
if (cellsGfix.size()>0)
	{
	if (time>tGfix) {
		output = getCbuffer("gfix",itGfix,cellsGfix.size()); 
		tGfix = output.t; cellsGfixData=output.d; itGfix +=1;
		}
	for (i=0; i<cellsGfix.size(); i++) {condGfix[i] = static_cast<int>(cellsGfixData[i*4+2]);}
	}
//Info<<" gfix "<<endl;
//for (j=0;j<nxyz;j++) { for (ig=0; ig<ph_gcomp;ig++) {Info<<Cg[ig]()[j]<<" ";} Info<<endl;}
	
for (ig=0; ig<ph_gcomp;ig++)
    {
		volScalarField& Cgi = Cg[ig]();//Info<<"i gas "<<ig<<endl;
		word nameCi ("Cg" + name(ig));
		GiSu *=0;GiSp *=0;
		
		// --------------------- internal conditions ------------------------
		if (cellsGwel.size()>0)  // gas in well, only for 2phase
		{
			float q=0;
			// set the values in the Su Sp matrix
			for (i=0; i<cellsGwel.size(); i++) 
				{ 
				i1=cellsGwel[i];v=mesh.V()[i1];
				if (flowType==4) {q=qPwel[i];} else {q=qHwel[i];}
				if (q<0) {GiSp[i1] = -q/v;} //in Sp term as discharge prop to local concentration
				else { GiSu[i1] = -q/v*gas_conc[condGwel[i]*ph_gcomp+ig]/Vmol[i1];} // ; } 
				//Info<<"g wel q conc "<<q<<" "<<gas_conc[condGwel[i]*ph_gcomp+ig]<<endl;
				}
		}
		//std::cout<<"before lbda \n";
		lbdag.value() = plugCgi.getLbdag(ig); // this is the plugin for degradation,needs to be here (compg specific)
		Info<<"lbdag "<<ig <<" "<<lbdag<<endl;
		//Info<<"gmax thetag "<<gMax(thetag)<<endl;
		
		if (reactIndicator==1) {
			for (j=0;j<nxyz;j++) {
				GiSu[j] += (Cg[ig]()[j] - freak.gm[ig*nxyz+j]/gvol[j])/float(runTime.deltaTValue())*thetag[j];
				//if (j<5) {Info<<j<< " inGsu Cg, fk.p fk.g/gv gsu "<<Cg[ig]()[j]<<" "<<freak.p[j]<<" "<<freak.gm[ig*nxyz+j]/gvol[j]<<" "<<GiSu[j]<<endl;}
			}
		}
		
		if (plugCgi.flag) {
			#include "plugins/plugin_Cgiupdate.H"
			}
		if (gMax(thetag)>1e-6) //solving this equation for very low thetag may give wrong results
		{
		fvScalarMatrix CgiEqn
		(		  
            thetag * fvm::ddt(Cgi) +  Cgi * fvc::ddt(thetag)
			+ fvc::div(phig, Cgi, "div(phig,Cgi)") // used only in multiphase
            - fvm::laplacian(Deffg, Cgi, "laplacian(Deffg,Cgi)")
			+ thetag * fvm::Sp(lbdag,Cgi) //decay
			+ fvm::SuSp(GiSp, Cgi)
			== - GiSu 
		); 
		
		for (i=0; i<cellsGfix.size(); i++) {CgiEqn.setReference(cellsGfix[i],gas_conc[condGfix[i]*ph_gcomp+ig]/Vmol[cellsGfix[i]]);}
		
		CgiEqn.solve();
		
		// after solve, calc gradient at Gfix cells
		/*
		volVectorField gr3 = fvc::grad(Cgi);
		size_t nc = cellsGfix.size();
		for (i=0; i<nc; i++) {Ggrd[ig*nc+i] = gr3()[cellsGfix[i]].y();}
		*/
	}
		Info<<" after solve, cell 0 "<<ig<<" " <<Cgi[0]<<endl;

}
//for (j=0;j<8;j++) { for (i=0;i<ph_gcomp;i++) {Info<<Cg[i]()[j]<<" ";} Info<<endl; }

reactIndicator = 0.;
//for (j=0;j<nxyz;j++) { for (ig=0; ig<ph_gcomp;ig++) {Info<<Cg[ig]()[j]<<" ";} Info<<endl;}

