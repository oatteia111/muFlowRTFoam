hp.storePrevIter();

fvScalarMatrix hpEqn
	( stor*fvm::ddt(hp) - fvm::laplacian(Mf,hp)  + fvc::div(phiGr) == b*fvOptions(sw,hp) ); //
fvOptions.constrain(hpEqn);
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