// combine p and S equations

// -################# - conditions  ---####################---
	pSu*=0;fSu*=0;size_t i; //previous cells with values must be destroyed
	bcPS.resize(nxyz,0);
	float time = mesh.time().value();float f; 
	icount = 0;	//get the ghb conditions

	// update fields
	Se =  pow(1.0+pow(max(-p/pc0,1.e-16),n_vg),-m_vg);
	pc = pc0*pow(pow(Se,-1/m_vg)-1,1/n_vg); //Se=(1+(p/pc0)^n)^-m; Se^-1/m-1=(p/pc0)^n
	gradpc = fvc::interpolate(fvc::grad(pc,"pc"));//
	phiPc = (Mwf * gradpc) & mesh.Sf();//Info<<"PhiPc(10) "<<phiPc[10]<<endl;
	
	krw = 1 * pow(Se,0.5) * pow(1-pow(1-pow(Se,1/m_vg),m_vg),2);
	krg = pow(1-Se,0.5) * pow(1-pow(Se,1/m_vg),2*m_vg);//Info<<"krg "<<krg<<endl; // OA 12/9 2*m -> 5*m
	krg = max(min(krg,.99999),1e-9); //Info<<"krw "<<krw[0]<<" krg "<<krg[0]<<endl; //krmax = 1
	krwf = fvc::interpolate(krw,"krw");
	krgf = fvc::interpolate(krg,"krg");

	//- mobility computation 
	//volScalarField Mw = K*krw/muw;
	Mwf = Kf *krwf/muwf;
	Mf = Mwf + Kf*krgf/mug;
	Lwf = Kf*rhowf*krwf/muwf;
	Lf = Lwf + Kf*rhog*krgf/mug;
	Fwf = (krwf/muwf) / (krgf/mug + krwf/muwf);
	phiGr = (Lf * g) & mesh.Sf();//Info<<"phiGr "<<"PhiGr(10) "<<phiGr[10]<<endl;
	//
	Info<<"before Peqn sw p pc phiPc "<<endl;for (i=0;i<ncell_lay;i++) {Info<<sw[i]<<" "<<p[i]<<" "<<pc[i]<<" "<<phiPc[i]<<" ";Info<<"\n";}

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

	fvScalarMatrix pEqn
        (
            fvm::laplacian(-Mf, p) + fvc::div(phiGr) 
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
	//------------------------- end Peqn
	// find values of capill curve derivative
	dpcdS = -pc0/(n_vg*m_vg)*pow(Se,-1/m_vg - 1); //a b x^(a - 1);a=-1/m,b=1/n
	dpcdS *= pow(pow(Se,-1/m_vg)-1,1/n_vg-1) ;//(x^a - 1)^(b - 1)
	surfaceScalarField dpcdSf ("dpcdSf",fvc::interpolate(dpcdS,"dpcdS"));//Info<<"krgf "<<krgf<<endl;	
	
	Info<<"after Peqn sw p pc phiPc dpcds "<<endl;for (i=0;i<ncell_lay;i++) {Info<<sw[i]<<" "<<p[i]<<" "<<pc[i]<<" "<<phiPc[i]<<" "<<dpcdS[i]<<" "<<endl;}

	//------------------------------- SEqn
	// nb as S is already used for solution, we use F instead (fraction of water)
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

	fvScalarMatrix SEqn
        (
            eps*fvm::ddt(sw) - fvc::laplacian(Mwf,p) + fvm::laplacian(Mwf * dpcdSf,sw) == fSu
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

	Info<<"after Seqn sw p pc phiPc "<<endl;for (i=0;i<ncell_lay;i++) {Info<<sw[i]<<" "<<p[i]<<" "<<pc[i]<<" "<<phiPc[i]<<" "<<endl;}

	volVectorField gp = fvc::grad(p);
	//recalculate fluxes (phiw, phig used for transport)
	phiP = pEqn.flux();
	phiP = - (Mf * fvc::snGrad(p)) * mesh.magSf();
	gradpc = fvc::interpolate(fvc::grad(pc,"pc"));									      
	phiPc = (Mwf*gradpc) &  mesh.Sf();
							
	phi == phiGr + phiP + phiPc*activateCapillarity;
	//phi = phiGr - ((Mf * fvc::interpolate(gp) ) & mesh.Sf()) + phiPc*activateCapillarity; //not usre about tha tone
    // keep Fwf (not Mwf/Mf, don't know why)
	phiw == (Lwf/Lf)*phiGr + Fwf*phiP  + phiPc*activateCapillarity; // == to transfer also BC //removed (1-Fwf)*
    //phiw == Fwf*phi + (((Lwf-Fwf*Lf)*g) & mesh.Sf()) + (1-Fwf)*phiPc*activateCapillarity; // old: wrong (in fact gravity=0?) but this is the only one that works!
    //phiw = (Lwf/Lf)*phiGr - ((Mwf * fvc::interpolate(gp) ) & mesh.Sf()) + phiPc*activateCapillarity; // does not work
	//phiw = ((Lwf * g) & mesh.Sf()) - ((Mwf * fvc::interpolate(gp) ) & mesh.Sf()) + phiPc*activateCapillarity;
    phig == phi - phiw;
	volScalarField q1 = fvc::div(phiw);
	for (i=0;i<cellsPfix.size();i++) {i1=cellsPfix[i];qPfix[i] = q1[i1]*mesh.V()[i1];} //Info<<" heq hfix "<<i<<" "<<i1<<" "<<qHfix[i]<<endl;
	//for (j=0; j<20;j++) {Info<<"gradp "<<gp[j]<<" phi "<<phi[j]<<" phiw "<<phiw[j]<<" phig "<<phig[j]<<endl;}
	Info<<"p 0 end "<<p[0]<<" "<<p[ncell_lay-1]<<" sw 0 end "<<sw[0]<<" "<<sw[ncell_lay-1]<<" rhow "<<rhow[0]<<" muw "<<muw[0]<<endl;
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
	Info<<" Uw ";for (i=0;i<24;i++) {Info<<mag(Uw[i])<<" ";};Info<<"\n";
	Info<<" Ug ";for (i=0;i<24;i++) {Info<<mag(Ug[i])<<" ";};Info<<"\n";
