hp.storePrevIter();
sw.storePrevIter();
/*
Info<<"hp ";for (int i=0;i<6;i++) {Info<<hp[i*12]<<" ";};Info<<hp[98]<<" "<<hp[99]<<endl;
Info<<"sw ";for (int i=0;i<8;i++) {Info<<sw[i*12]<<" ";};Info<<sw[98]<<" "<<sw[99]<<endl;
Info<<"stor ";for (int i=0;i<8;i++) {Info<<stor[i*12]<<" ";};Info<<stor[98]<<" "<<stor[99]<<endl;
*/

// -------------- conditions  ----------------------
	hSu=hSu*0; hSp=hSp*0;size_t i; //previous cells with values must be destroyed
	time = mesh.time().value() ; //Info<<"chw size "<<cellsHw.size()<<" chfx size "<<cellsHfx.size()<<endl;//const Time&
	//get the well discharge
	if (cellsHwel.size()>0)
	{
		if (time>=tHwel) {
			output = getCbuffer("hwel",itHwel,cellsHwel.size()); 
			tHwel = output.t; cellsHwelData=output.d; itHwel+=1;
			newDeltaT = runTime.deltaTValue()/20;
			}
		for (i=0; i<cellsHwel.size(); i++) {
			i1=cellsHwel[i];qHwel[i] = cellsHwelData[i*4+2] ; hSu[i1] = qHwel[i]/mesh.V()[i1];
			//Info<<i<<" ic "<<cellsHwel[i]<<" q "<<qHwel[i]<<" hSu "<<hSu[cellsHwel[i]]<<endl;
			}
	}
	//---- recharge conditions		
	if (cellsHrch.size()>0)
	{
		if (time>=tHrch) {
			output = getCbuffer("hrch",itHrch,cellsHrch.size()); 
			tHrch = output.t; cellsHrchData=output.d; itHrch +=1;
			}
		float frech;float pi=3.141592654; //frech is used to stop recharge or evapo when hp is too high or too low
		for (i=0; i<cellsHrch.size(); i++) 
			{
			i1 = cellsHrch[i];clayRch[i] = i1;qHrch[i] = cellsHrchData[i*4+2] ; // here the injection is at the exact place of recharge (not like in unconfined)
			//Info<<i<<" ic "<<cellsHrch[i]<<" i lay "<<clayRch[i]<<" rch "<<qHrch[i]<<endl;
			if (hSu[i1]==0) { // apply recharge only if there is nothing else in the cell
				if (hp()[i1]<0) {frech = 0.5+std::atan((hp()[i1]+20.)*10.)/pi;} 
				else {frech = 0.5-std::atan((hp()[i1]-10.)*20.)/pi;}
				hSu[i1] = frech*qHrch[i]/mesh.V()[i1]; 
				}
			} 	
	} //Info<<"Hrch read ";//<<" Hghb "<<hSu[ncell_lay*3+ic0]<<" "<<endl;

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
phiw =  (phiGr - ((Mf * fvc::interpolate(fvc::grad(hp)) ) & mesh.Sf() ) ) ;//fvc::interpolate(sw) *

volScalarField q1 = fvc::div(phiw);
//volScalarField d1 = stor*fvc::ddt(hp);
for (i=0;i<cellsHfix.size();i++) {i1=cellsHfix[i];qHfix[i] = q1[i1]*mesh.V()[i1];Info<<" heq hfix "<<i<<" "<<i1<<" "<<qHfix[i]<<endl;}

//Info<<"hp -relax";for (int i=0;i<6;i++) {Info<<hp[i*12]<<" ";};Info<<hp[98]<<" "<<hp[99]<<endl;
//hp.relax(); // mix with previous solution to stabilise, relax factor in fvsolu
//if(timeStepControl == "Picard")  {resPicard=residu;}
