hp.storePrevIter();

// -------------- conditions  ----------------------
	hSu*=0; hSp*=0; //previous cells with values must be destroyed
	float time = mesh.time().value(); //Info<<"chw size "<<cellsHw.size()<<" chfx size "<<cellsHfx.size()<<endl;//const Time&
	//get the well discharge
	if (cellsHw.size()>0)
	{
		float tnow = cellsHwData[2]*day; //time in files is in days
		while (time>=tnow) {icount += cellsHw.size(); tnow = cellsHwData[2+icount*4]*day; }
		if (icount>0) {icount -= cellsHw.size();};
		for (int i=0; i<cellsHw.size(); i++) 
			{v=mesh.V()[cellsHw[i]];qHw[i] = cellsHwData[2+(icount+i)*4+2]/day; hSu[cellsHw[i]] = qHw[i]/v;
			//{hSu[cellsHw[i]] = cellsHwData[2+(icount+i)*4+2]/day;} // discharge must be in s-1
			Info<<i<<" ic "<<cellsHw[i]<<" q "<<qHw[i]<<endl; }
	}
	//get the fixed head conditions
	icount = 0;
	if (cellsHfx.size()>0)
	{
		float tnow = cellsHfxData[2]*day; //time in files is in days
		while (time>=tnow) {icount += cellsHfx.size(); tnow = cellsHfxData[2+icount*4]*day; }
		if (icount>0) {icount -= cellsHfx.size();} 
		Info<<"icount "<<icount<<endl;
		for (int i=0; i<cellsHfx.size(); i++) 
			{hSu[cellsHfx[i]] = cellsHfxData[2+(icount+i)*4+2]*1e4; hSp[cellsHfx[i]] = 1e4; }
			//Info<<i<<" ic "<<cellsHfx[i]<<" icdat "<<cellsHfxData[2+(icount+i)*4+1]<<endl; }
	}

fvScalarMatrix hpEqn
	( stor*fvm::ddt(hp) - fvm::laplacian(Mf,hp)  + fvc::div(phiGr) 
		- hSu + fvm::SuSp(hSp,hp) ); //== b*fvOptions(sw,hp) ); //
//fvOptions.constrain(hpEqn);
//Info<< fvOptions(sw,hp) <<endl;
hpEqn.solve(); //Info << "h 500 " << h[500] << " iter "<<runTime.solverInfo()<<endl;//nIterations() //hEaqn.solver() ok //SolverPerformance() not ok
//Info<<"end solve"<<endl;
//phiw =  hpEqn.flux() ;
phiw =  (phiGr - ((Mf * fvc::interpolate(fvc::grad(hp)) ) & mesh.Sf() ) ) ;//fvc::interpolate(sw) *
//phiw ==  phiGr - (Mf * fvc::snGrad(hp)) * mesh.magSf() ;//same as hpeqn
//phiw = phiGr - (Mf * fvc::grad(h))* mesh.magSf();Info<<" phiw calc"<<endl;

//resPicard = gMax( mag(h-h.prevIter())->internalField()/mag(h)->internalField() );	Info<<" res picard "<<	resPicard<<endl;	 //->internalField()	
//resPicard = 0;
//forAll(hp, i) {resPicard = max(resPicard,(hp()[i]-hp.prevIter()[i])/hp()[i]);}
resPicard = gMax((mag(hp-hp.prevIter()))->internalField());
//Info<<"phiGr 90 "<<phiGr[90]<<" phiw 90 "<<phiw[90]<<" res picard "<<	resPicard<<endl;