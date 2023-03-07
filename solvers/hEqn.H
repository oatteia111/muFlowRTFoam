double x;int ilay,nc,nc1;Info<<"nlay "<<nlay<<" ncl "<<ncell_lay<<endl;//double pi = 3.141592654;
double e=5e-4;
//Info << "iter picard n "<< iterPicard << endl;
h.storePrevIter();Info<<"flowtype "<<flowType<<endl;

// -------------- conditions  ----------------------
//get the well discharge 
condHw.resize(cellsHw.size());
float tnow = cellsHwData[2]*day; //time in files is in days
float time = mesh.time().value(); //const Time&
while (time>=tnow)
{
	icount += cellsHw.size(); 
	tnow = cellsHwData[2+icount*4]*day; 
}
if (icount>0) {icount -= cellsHw.size();};
for (int i=0; i<cellsHw.size(); i++) {condHw[i] = cellsHwData[2+(icount+i)*4+2]/day;} // discharge must be in s-1
//get the fixed head conditions
std::vector<double> condHfx(cellsHfx.size(),0.);
while (time>=tnow)
{
	icount += cellsHfx.size(); 
	tnow = cellsHfxData[2+icount*4]*day; 
}
if (icount>0) {icount -= cellsHfx.size();};
for (int i=0; i<cellsHfx.size(); i++) {condHfx[i] = cellsHfxData[2+(icount+i)*4+2];Info<<"hfx "<<condHfx[i]<<endl;} // discharge must be in s-1


// -----------------  set the effective thickness for unconfined flow -------------------
if (flowType == 2) 
	{
	for (ilay=0;ilay<nlay;ilay++) //ilay will be used from top
		{ 
		for (i=0;i<ncell_lay;i++) 
			{
				nc = (nlay-ilay-1)*ncell_lay + i;
				x = (h[nc]-zbot[nc])/thk[nc];
				if (x<e && x>0) {sw[nc]=pow(x,2)/(4*e)+x/2+e/4;} // same eq as in mdf USg (transport manual fg 21)
				else if (x>1-e && x<1) {sw[nc]=1-pow((1-x),2)/(4*e)-(1-x)/2-e/4;}
				else {sw[nc]=std::max(std::min(x,sw_max),sw_min[nc]);}
			}
		}
	forAll( mesh.owner() , iface ) 
		{
		// Get the owner and neighbour cell index
		oCell = mesh.owner()[iface];odCface = mag(oCf[iface]-oC[oCell]);
		nCell = mesh.neighbour()[iface];ndCface = mag(oCf[iface]-oC[nCell]);
		if (abs(mesh.Sf()[iface].z())<mag(mesh.Sf()[iface])*1e-5)  //vertical faces
			{
			scalar swm = (sw[oCell]*sw[nCell])/(odCface*sw[oCell]+ndCface*sw[nCell])*(ndCface+odCface);
			Mf[iface] = M0f[iface]*swm;
			}
		}
		
	}
forAll(stor(), celli) { 
	if (sw[celli]<0.999) {stor[celli] = eps[celli]/mesh.V()[celli];}// //Info << stor[celli]<<endl;}  /thk[celli]
	else {stor[celli] = stor0/mesh.V()[celli] ;} //   /thk[celli]
	//stor[celli] = (stor0+eps[celli]*(0.49695+1.01/pi*std::atan((0.96-sw[celli])*1000)))/h[celli];// approx formula tomicmic a smoothed step
	};  Info << "stor0 " << stor0 <<" stor(4) " << stor[4] << " eps " << eps[4]<<" thk " <<thk[4]<<" vol "<<mesh.V()[4]<<endl;
stor.dimensions().reset(dimless/dimLength);
// ***************** solving
//Info << "Kf 2026 "<< Kf[2026]<< " " <<Kf[2026]<< endl;

fvScalarMatrix hEqn
	( stor*fvm::ddt(h) - fvm::laplacian(Mf,h) == b*fvOptions(sw,h) ); // + fvc::div(phiG)
fvOptions.constrain(hEqn);
hEqn.solve(); //Info << "h 500 " << h[500] << " iter "<<runTime.solverInfo()<<endl;//nIterations() //hEaqn.solver() ok //SolverPerformance() not ok
Info<<"end solve"<<endl;
//phiw= Mf*fvc::flux(fvc::grad(h));
phiw =  hEqn.flux() ;//Info<<"phi "<<phiw<<endl;
//phiw = phiGr - ((Mf * fvc::interpolate(fvc::grad(h)) ) & mesh.Sf() )  ;//problem at boundary take grad wiht the bc face
//phiw = phiGr - ((Mf * fvc::grad(h) ) & mesh.Sf() )  ; doe snot work
//phiw = phiG - (Mf * fvc::grad(h))* mesh.magSf();Info<<" phiw calc"<<endl;

//find the flux at hfix cells
/*for (int ic=0; ic<cellsHfx.size(); ic++) 
	{
	const cell& faces = mesh.cells()[cellsHfx[ic]];
    forAll( faces, i ) {
		Info<<"cell "<<cellsHfx[ic]<<" face "<<faces[i]<<" phi "<<phiw[faces[i]]<<endl;
		if (mesh.owner()[faces[i]] == cellsHfx[ic]) {Info<<"its owner "<<endl;}
		if (mesh.neighbour()[faces[i]] == cellsHfx[ic]) {Info<<"its nghb "<<endl;}
		}
	}*/

