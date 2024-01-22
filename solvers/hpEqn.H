hp.storePrevIter();

// -------------- conditions  ----------------------
	hSu=hSu*0; hSp=hSp*0;size_t i; //previous cells with values must be destroyed
	time = mesh.time().value()/day; //Info<<"chw size "<<cellsHw.size()<<" chfx size "<<cellsHfx.size()<<endl;//const Time&
	//get the well discharge
	if (cellsHwel.size()>0)
	{
		if (time>=tHwel) {
			output = getCbuffer("hwel",itHwel,cellsHwel.size()); 
			tHwel = output.t; cellsHwelData=output.d; itHwel+=1;
			runTime.setDeltaT(runTime.deltaTValue()/20);tnext= tHwel;
			}
		for (i=0; i<cellsHwel.size(); i++) {
			v=mesh.V()[cellsHwel[i]];qHwel[i] = cellsHwelData[i*4+2]/day; hSu[cellsHwel[i]] = qHwel[i]/v;
			Info<<i<<" ic "<<cellsHwel[i]<<" q "<<qHwel[i]<<endl;
			}
	}

fvScalarMatrix hpEqn
	( stor*fvm::ddt(hp) - fvm::laplacian(Mf,hp)  + fvc::div(phiGr) 
		- hSu + fvm::SuSp(hSp,hp) ); //== b*fvOptions(sw,hp) ); //
//hp fixed
if (cellsHfix.size()>0)
{
	if (time>=tHfix) {
		output = getCbuffer("hfix",itHfix,cellsHfix.size()); 
		tHfix = output.t; cellsHfixData=output.d; itHfix +=1;
		tnext= tHfix;
		}
	for (i=0; i<cellsHfix.size(); i++) {hpEqn.setReference(cellsHfix[i],cellsHfixData[i*4+2]);}
}

hpEqn.solve(); //Info << "h 500 " << h[500] << " iter "<<runTime.solverInfo()<<endl;//nIterations() //hEaqn.solver() ok //SolverPerformance() not ok
//Info<<"end solve"<<endl;
//phiw =  hpEqn.flux() ;
//Info<<"phiGr ";for (int i=0;i<4;i++) {Info<<phiGr[i]<<" ";}
Info<<endl;
phiw =  (phiGr - ((Mf * fvc::interpolate(fvc::grad(hp)) ) & mesh.Sf() ) ) ;//fvc::interpolate(sw) *
//phiw ==  phiGr - (Mf * fvc::snGrad(hp)) * mesh.magSf() ;//same as hpeqn
//phiw = phiGr - (Mf * fvc::grad(h))* mesh.magSf();Info<<" phiw calc"<<endl;

//resPicard = gMax( mag(h-h.prevIter())->internalField()/mag(h)->internalField() );	Info<<" res picard "<<	resPicard<<endl;	 //->internalField()	
//resPicard = 0;
//forAll(hp, i) {resPicard = max(resPicard,(hp()[i]-hp.prevIter()[i])/hp()[i]);}
resPicard = gMax((mag(hp-hp.prevIter()))->internalField());
//Info<<"phiGr 90 "<<phiGr[90]<<" phiw 90 "<<phiw[90]<<" res picard "<<	resPicard<<endl;