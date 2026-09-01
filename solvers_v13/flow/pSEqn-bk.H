// combine p and S equations

// -################# - conditions  ---####################---
	pSu*=0;fSu*=0;size_t i; //previous cells with values must be destroyed
	bcPS.resize(nxyz,0);
	float time = mesh.time().value();float f; 
	icount = 0;	//get the ghb conditions
	//Info<< " p "<<p<<"\n sw "<<sw<<"\n rhow "<<rhow<<"\n Mf "<<Mf<<"\n Lwf "<<Lwf<<"\n phiGr "<<phiGr<<"\n phiw "<<phiw<<endl;
	//------------------------------- SEqn
	// nb as S is already used for solution, we use F instead (fraction of water)
	volVectorField gp = fvc::grad(p);//Info<<"gp after Seq "<<gp<<endl;
	if (cellsFwel.size()>0)
		{
		if (time>tFwel) {
			output = getCbuffer("fwel",itFwel,cellsFwel.size()); 
			tFwel = output.t; cellsFwelData=output.d; itFwel +=1;
			tnext = min(tnext,tFwel);
			}
		for (i=0; i<cellsFwel.size(); i++) 
			{
			i1=cellsFwel[i];qFwel[i] = cellsFwelData[i*4+2] ; fSu[i1] = qFwel[i]/mesh.V()[i1];bcPS[i1]=1;
			//Info<<"f wel "<<i<<" "<<i1<<" "<<cellsFwelData[i*4+2]<<" fSu "<<fSu[i1]<<endl;
			}
		}

	
	/*
	dpcdS == -pc0/(n_vg*m_vg)*pow(Se,-1/m_vg - 1); //a b x^(a - 1);a=-1/m,b=1/n
	dpcdS *= pow(pow(Se,-1/m_vg)-1,1/n_vg-1) ;//(x^a - 1)^(b - 1)
	surfaceScalarField dpcdSf ("dpcdSf",fvc::interpolate(dpcdS,"dpcdS"));//Info<<"krgf "<<krgf<<endl;	
	*/
	fvScalarMatrix SEqn
        (
            eps*fvm::ddt(sw) + fvc::div(phiw) == fSu
            //eps*fvm::ddt(sw) - fvc::laplacian(Mwf,p) + fvm::laplacian(Mwf * dpcdSf,sw) == fSu
        );
		
	if (cellsFfix.size()>0)
		{
		if (time>tFfix) {
			output = getCbuffer("ffix",itFfix,cellsFfix.size()); 
			tFfix = output.t; cellsFfixData=output.d; itFfix +=1;
			tnext= min(tnext,tFfix);
			}
		for (i=0; i<cellsFfix.size(); i++) 
			{SEqn.setReference(cellsFfix[i],cellsFfixData[i*4+2]);bcPS[cellsFfix[i]]=1;
			//Info<<"f fix "<<i<<" "<<cellsFfix[i]<<" "<<cellsFfixData[i*4+2]<<" fSu "<<fSu[cellsFfix[i]]<<endl;
			}
		}
		
	SEqn.solve();	
	//gp = fvc::grad(p);Info<<"gp after Seq "<<gp<<endl;

	//Info<<" after Seq sw 1-Se Mwf p-tmp p pc diffpc phiPc phiGr gp phiP phiw"<<endl;for (i=0;i<ncell_lay-1;i++) {Info<<sw[i]<<" "<<1-Se[i]<<" "<<Mwf[i]<<" "<<p_tmp[i]<<" "<<p[i]<<" "<<pc[i]<<" "<<pc[i+1]-pc[i]<<" "<<phiPc[i]<<" "<<phiGr[i]<<" "<<mag(gp[i])<<" "<<phiP[i]<<" "<<phiw[i]<<endl;}
	//---------------------------- update	
	if (plugPS.flag)
		{
		#include "plugins/plugin_PSupdate.H"
		}
	else
		{
		#include "flow/update2phaseFields.H"
		}
	//gp = fvc::grad(p);Info<<"after upd "<<gp<<endl;
	//Info<<"p 0 end "<<p[0]<<" "<<p[ncell_lay-1]<<" sw 0 end "<<sw[0]<<" "<<sw[ncell_lay-1]<<" rhow "<<rhow[0]<<" muw "<<muw[0]<<endl;
	//Info<<" after upd sw 1-Se Mwf p-tmp p pc diffpc phiPc phiGr gp phiP phiw"<<endl;for (i=0;i<ncell_lay-1;i++) {Info<<sw[i]<<" "<<1-Se[i]<<" "<<Mwf[i]<<" "<<p_tmp[i]<<" "<<p[i]<<" "<<pc[i]<<" "<<pc[i+1]-pc[i]<<" "<<phiPc[i]<<" "<<phiGr[i]<<" "<<mag(gp[i])<<" "<<phiP[i]<<" "<<phiw[i]<<endl;}

	//------------------------- start Peqn
	if (cellsPwel.size()>0)
		{
		if (time>tPwel) {
			output = getCbuffer("pwel",itPwel,cellsPwel.size()); 
			tPwel = output.t; cellsPwelData=output.d; itPwel +=1;
			tnext = min(tnext,tPwel);
			}
		for (i=0; i<cellsPwel.size(); i++) 
			{
			i1=cellsPwel[i];qPwel[i] = cellsPwelData[i*4+2] ; pSu[i1] = qPwel[i]/mesh.V()[i1];bcPS[i1]=1;
			//Info<<"p wel "<<i<<" "<<i1<<" "<<cellsPwelData[i*4+2]<<" pSu "<<pSu[i1]<<endl;
			}
		}
	//Info<<"Pfwel read ";
		if ((reactIndicator==1)&&(activateEbullition==1)) {
			for (j=0;j<nxyz;j++) {
				pSu[j] += (freak.p[j]/(p.oldTime()[j]/atmPa)-1)*eps[j]*(sw_max[j]-sw.oldTime()[j])/float(runTime.deltaTValue());
				if (j<5) {Info<<j<< " inPsu p, fk.p psu "<<p.oldTime()[j]/atmPa<<" "<<freak.p[j]<<" "<<pSu[j]<<endl;}
			}
		}

	for (j=0;j<nxyz;j++) {comp_g[j] = gvol[j]/p[j];}
	fvScalarMatrix pEqn
        (
             //comp_g*fvm::ddt(p) 
			 - fvm::laplacian(Mf, p) + fvc::div(phiGr) 
            + fvc::div(phiPc)*activateCapillarity
			== pSu
        );
		
	if (cellsPfix.size()>0)
		{
		if (time>tPfix) {
			output = getCbuffer("pfix",itPfix,cellsPfix.size()); 
			tPfix = output.t; cellsPfixData=output.d; itPfix +=1;
			tnext= min(tnext,tPfix);
			}
		for (i=0; i<cellsPfix.size(); i++) 
			{pEqn.setReference(cellsPfix[i],cellsPfixData[i*4+2]);bcPS[cellsPfix[i]]=1;
			//Info<<"p fix "<<i<<" "<<cellsPfix[i]<<" "<<cellsPfixData[i*4+2]<<" pSu "<<pSu[cellsPfix[i]]<<endl;
			}
		}
	
	pEqn.solve();
	
	gp = fvc::grad(p);//Info<<gp<<endl;
	//Info<<" after Peq sw 1-Se Mwf p-tmp p pc diffpc phiPc phiGr gp phiP phiw"<<endl;for (i=0;i<ncell_lay-1;i++) {Info<<sw[i]<<" "<<1-Se[i]<<" "<<Mwf[i]<<" "<<p_tmp[i]<<" "<<p[i]<<" "<<pc[i]<<" "<<pc[i+1]-pc[i]<<" "<<phiPc[i]<<" "<<phiGr[i]<<" "<<mag(gp[i])<<" "<<phiP[i]<<" "<<phiw[i]<<endl;}
	//recalculate fluxes (phiw, phig used for transport)
	//phiP = pEqn.flux();
	phiP == - (Mf * fvc::snGrad(p)) * mesh.magSf();
	//phiP ==  ((-Mf * fvc::interpolate(gp)) & mesh.Sf());
	gradpc = fvc::interpolate(fvc::grad(pc,"pc"));									      
	phiPc = (Mwf*gradpc) &  mesh.Sf();
	/*
	Info<<" sw ";for (i=0;i<24;i++) {Info<<sw[i]<<" ";};Info<<"\n";
	Info<<" Mf ";for (i=0;i<24;i++) {Info<<Mf[i]<<" ";};Info<<"\n";
	Info<<" p ";for (i=0;i<24;i++) {Info<<p[i]<<" ";};Info<<"\n";
	Info<<" gradp ";for (i=0;i<24;i++) {Info<<mag(gp[i])<<" ";};Info<<"\n";
	Info<<" phiPc ";for (i=0;i<24;i++) {Info<<phiPc[i]<<" ";};Info<<"\n";
	Info<<" phiP ";for (i=0;i<24;i++) {Info<<phiP[i]<<" ";};Info<<"\n";
	*/				
	phi == phiGr + phiP + phiPc*activateCapillarity;
	//phi = phiGr - ((Mf * fvc::interpolate(gp) ) & mesh.Sf()) + phiPc*activateCapillarity; //not usre about tha tone
    // keep Fwf (not Mwf/Mf, don't know why)
	phiw == (Lwf/Lf)*phiGr + Fwf*phiP  + phiPc*activateCapillarity; // == to transfer also BC //removed (1-Fwf)*
    //phiw == Fwf*phi + (((Lwf-Fwf*Lf)*g) & mesh.Sf()) + (1-Fwf)*phiPc*activateCapillarity; // old: wrong (in fact gravity=0?) but this is the only one that works!
    //phiw == (Lwf/Lf)*phiGr - ((Mwf * fvc::interpolate(gp) ) & mesh.Sf()) + phiPc*activateCapillarity; // does not work
	//phiw = ((Lwf * g) & mesh.Sf()) - ((Mwf * fvc::interpolate(gp) ) & mesh.Sf()) + phiPc*activateCapillarity;
    phig == phi - phiw;
	volScalarField q1 = fvc::div(phiw);
	for (i=0;i<cellsPfix.size();i++) {i1=cellsPfix[i];qPfix[i] = q1[i1]*mesh.V()[i1];} //Info<<" heq hfix "<<i<<" "<<i1<<" "<<qHfix[i]<<endl;
	//for (j=0; j<20;j++) {Info<<"gradp "<<gp[j]<<" phi "<<phi[j]<<" phiw "<<phiw[j]<<" phig "<<phig[j]<<endl;}
	//Info<<" after Peq Mf p-tmp p pc diffpc phiPc "<<endl;for (i=0;i<ncell_lay-1;i++) {Info<<sw[i]<<" "<<Mf[i]<<" "<<p_tmp[i]<<" "<<p[i]<<" "<<pc[i]<<" "<<pc[i+1]-pc[i]<<" "<<phiPc[i]<<" ";Info<<"\n";}
	//Info<<p;
	//Info<<"p solved"<<endl;

    U = fvc::reconstruct(phi);
    U.correctBoundaryConditions();
	//Info<<U;
    Uw = fvc::reconstruct(phiw);
    Ug = U-Uw;
	//Info<<Ug;
    Ug.correctBoundaryConditions();
    Uw.correctBoundaryConditions();  
