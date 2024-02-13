hp.storePrevIter();
sw.storePrevIter();
/*
Info<<"hp ";for (int i=0;i<6;i++) {Info<<hp[i*12]<<" ";};Info<<hp[98]<<" "<<hp[99]<<endl;
Info<<"sw ";for (int i=0;i<8;i++) {Info<<sw[i*12]<<" ";};Info<<sw[98]<<" "<<sw[99]<<endl;
Info<<"stor ";for (int i=0;i<8;i++) {Info<<stor[i*12]<<" ";};Info<<stor[98]<<" "<<stor[99]<<endl;
*/

// -------------- conditions  ----------------------
	hSu=hSu*0; hSp=hSp*0;size_t i; //previous cells with values must be destroyed
	time = mesh.time().value()/day; //Info<<"chw size "<<cellsHw.size()<<" chfx size "<<cellsHfx.size()<<endl;//const Time&
	//get the well discharge
	if (cellsHwel.size()>0)
	{
		if (time>=tHwel) {
			output = getCbuffer("hwel",itHwel,cellsHwel.size()); 
			tHwel = output.t; cellsHwelData=output.d; itHwel+=1;
			runTime.setDeltaT(runTime.deltaTValue()/20);
			}
		for (i=0; i<cellsHwel.size(); i++) {
			v=mesh.V()[cellsHwel[i]];qHwel[i] = cellsHwelData[i*4+2]/day; hSu[cellsHwel[i]] = qHwel[i]/v;
			Info<<i<<" ic "<<cellsHwel[i]<<" q "<<qHwel[i]<<" hSu "<<hSu[cellsHwel[i]]<<endl;
			}
	}
Info<<"dt "<<runTime.deltaT()<<endl;

fvScalarMatrix hpEqn
	(stor*fvm::ddt(hp) - fvm::laplacian(Mf,hp) + fvc::div(phiGr) 
		- hSu + fvm::SuSp(hSp,hp)); //

if(timeStepControl == "Picard")  {hpEqn +=  (stor*(hp.oldTime()-hp.prevIter())+(sw.prevIter() - sw.oldTime()))/ runTime.deltaT();} //

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

//residu = hpEqn.solve().initialResidual(); //initial residual does not seem to include hSu??

hpEqn.solve(); 

//Info<<"hp -relax";for (int i=0;i<6;i++) {Info<<hp[i*12]<<" ";};Info<<hp[98]<<" "<<hp[99]<<endl;
//hp.relax(); // mix with previous solution to stabilise, relax factor in fvsolu
if(timeStepControl == "Picard")  {resPicard=residu;}
//Info<<"end solve"<<endl;
//phiw =  hpEqn.flux() ;
//Info<<"phiGr ";for (int i=0;i<4;i++) {Info<<phiGr[i]<<" ";}
//Info<<"hp ";for (int i=0;i<6;i++) {Info<<hp[i*12]<<" ";};Info<<hp[98]<<" "<<hp[99]<<endl;
//Info<<"sw ";for (int i=0;i<8;i++) {Info<<sw[i*12]<<" ";};Info<<sw[98]<<" "<<sw[99]<<endl;
//Info<<"stor ";for (int i=0;i<8;i++) {Info<<stor[i*12]<<" ";};Info<<stor[98]<<" "<<stor[99]<<endl;
//Info<<"hpprev ";for (int i=0;i<8;i++) {Info<<hp.prevIter()[i*12]<<" ";};Info<<endl;
//Info<<"hpold ";for (int i=0;i<8;i++) {Info<<hp.oldTime()[i*12]<<" ";};Info<<endl;
//Info<<"swprev ";for (int i=0;i<8;i++) {Info<<sw.prevIter()[i*12]<<" ";};Info<<endl;
//Info<<"swold ";for (int i=0;i<8;i++) {Info<<sw.oldTime()[i*12]<<" ";};Info<<endl;
phiw =  (phiGr - ((Mf * fvc::interpolate(fvc::grad(hp)) ) & mesh.Sf() ) ) ;//fvc::interpolate(sw) *
//phiw ==  phiGr - (Mf * fvc::snGrad(hp)) * mesh.magSf() ;//same as hpeqn
//phiw = phiGr - (Mf * fvc::grad(h))* mesh.magSf();Info<<" phiw calc"<<endl;
