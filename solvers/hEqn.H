double x;int ilay,nc,nc1;//Info<<"nlay "<<nlay<<" ncl "<<ncell_lay<<endl;//double pi = 3.141592654;
double e=1e-5;
//Info << "iter picard n "<< iterPicard << endl;
h.storePrevIter();//Info<<"flowtype "<<flowType<<endl;

//for (i=0;i<150;i++) {Info << "cID " <<i<<" "<< mesh.C()[i] << endl;}

// -----------------  set the effective thickness for unconfined flow -------------------
if (flowType == 2) 
	{
	for (ilay=0;ilay<nlay;ilay++) //find sw
		{ 
		for (i=0;i<ncell_lay;i++) 
			{
				nc = (nlay-ilay-1)*ncell_lay + i; // top layer 1st
				x = (h[nc]-zbot[nc])/thk[nc];
				/*if (x<e && x>0) {sw[nc]=pow(x,2)/(4*e)+x/2+e/4;} // same eq as in mdf USg (transport manual fg 21)
				else if (x>1-e && x<1) {sw[nc]=1-pow((1-x),2)/(4*e)-(1-x)/2-e/4;}
				else {*/
				sw[nc]=std::max(std::min(x,sw_max),sw_min[nc]);//}
			}
		}
		
	forAll( mesh.owner() , iface ) //calculate Mf
		{
		oCell = mesh.owner()[iface];odCface = mag(oCf[iface]-oC[oCell]);
		nCell = mesh.neighbour()[iface];ndCface = mag(oCf[iface]-oC[nCell]);
		if (abs(mesh.Sf()[iface].z())<mag(mesh.Sf()[iface])*1e-5)  //vertical faces
			{
			scalar swm = (sw[oCell]*sw[nCell])/(odCface*sw[oCell]+ndCface*sw[nCell])*(ndCface+odCface);
			Mf[iface] = M0f[iface]*swm;
			}
		else {Mf[iface] = M0f[iface];}
		/*else //other faces if sw=sw_min Mf=Small
			{
			if ((sw[oCell]<sw_min[oCell]*1.1)&&(sw[nCell]<sw_min[nCell]*1.1)) {Mf[iface] = M0f[iface]*sw_min[oCell];}
			}*/
		}
	}
forAll(stor(), celli) 
	{ 
	if (sw[celli]<sw_max) {stor[celli] = eps[celli]/thk[celli];}// //Info << stor[celli]<<endl;}  
	else {stor[celli] = stor0/thk[celli];} //   /thk[celli]
	//stor[celli] = (stor0+eps[celli]*(0.49695+1.01/pi*std::atan((0.96-sw[celli])*1000)))/h[celli];// approx formula tomicmic a smoothed step
	};  //Info << "stor0 " << stor0 <<" stor(4) " << stor[4] << " eps " << eps[4]<<" thk " <<thk[4]<<" vol "<<mesh.V()[4]<<endl;
stor.dimensions().reset(dimless/dimLength);

// -################# - conditions  ---####################---
	hSu=hSu*0; hSp=hSp*0; //previous cells with values must be destroyed
	float time = mesh.time().value();float f; //Info<<"chw size "<<cellsHw.size()<<" chfx size "<<cellsHfx.size()<<endl;//const Time&
	//get the well discharge/recharge
	icount = 0;	//get the ghb conditions
	if (cellsHw.size()>0)
	{
		float tnow = cellsHwData[2]*day;//time in files is in days
		while (time>=tnow) {icount += cellsHw.size(); tnow = cellsHwData[2+icount*4]*day; }
		if (icount>0) {icount -= cellsHw.size();};
		for (int i=0; i<cellsHw.size(); i++)
			{v=mesh.V()[cellsHw[i]];qHw[i] = cellsHwData[2+(icount+i)*4+2]/day; hSu[cellsHw[i]] = qHw[i]/v;}
		//Info<<i<<" ic "<<cellsHw[i]<<" icdat "<<cellsHwData[2+(icount+i)*4+1]<<endl; }
	} Info<<"Hfw read ";
	//----- get the ghb conditions
	icount = 0;	//get the ghb conditions
	int ic0=10;
	if (cellsHghb.size()>0)
	{
		float tnow = cellsHghbData[2]*day;
		while (time>=tnow) {icount += cellsHghb.size(); tnow = cellsHghbData[2+icount*4]*day; }
		if (icount>0) {icount -= cellsHghb.size();} 
		for (int i=0; i<cellsHghb.size(); i++) 
			{i1=cellsHghb[i];v=mesh.V()[i1];hSu[i1] = cellsHghbData[2+(icount+i)*4+2]/v/day; hSp[i1] = cellsHghbData[2+(icount+i)*4+3]/v/day;}
	} Info<<"Hghb read hSu ";//<<hSu[ncell_lay*3+ic0]<<" ";
	//----- get the riv conditions
	icount = 0; 
	if (cellsHriv.size()>0)
	{
		float tnow = cellsHrivData[2]*day;
		while (time>=tnow) {icount += cellsHriv.size(); tnow = cellsHrivData[2+icount*4]*day; }
		if (icount>0) {icount -= cellsHriv.size();} 
		for (int i=0; i<cellsHriv.size(); i++) 
			{i1=cellsHriv[i];v=mesh.V()[i1];hSu[i1] = cellsHrivData[2+(icount+i)*4+2]/v/day; hSp[i1] = cellsHrivData[2+(icount+i)*4+3]/v/day;
			Info<<" i "<<i<<" i1 "<<i1<<" hsu "<<hSu[i1]<<" hsp "<<hSp[i1]<<endl;}
	} //Info<<"Hriv read "<<endl;
	//----- get the drn conditions
	icount = 0; 
	if (cellsHdrn.size()>0)
	{
		float tnow = cellsHdrnData[2]*day; 
		while (time>=tnow) {icount += cellsHdrn.size(); tnow = cellsHdrnData[2+icount*4]*day; }
		if (icount>0) {icount -= cellsHdrn.size();} 
		for (int i=0; i<cellsHdrn.size(); i++) 
			{
			i1 = cellsHdrn[i]; v=mesh.V()[i1];
			hSu[i1] = cellsHdrnData[2+(icount+i)*4+2]/v/day; hSp[i1]=cellsHdrnData[2+(icount+i)*4+3]/v/day;
			//hSp[i1] = cellsHdrnData[2+(icount+i)*4+3]/day;
			//scalar dff = hSu[i1]-hSp[i1]*h[i1];
			//if (dff>0) {hSu[i1]/=1e4;hSp[i1]/=1e4;}
			//if (i<5) {Info<<"drn i1 "<<i1<<" h(i1) "<<h()[i1]<<" data1 "<< cellsHdrnData[2+(icount+i)*4+2]<<" data 2 "<< cellsHdrnData[2+(icount+i)*4+3]<<" hSu "<<hSu[i1]<<" hsp "<<hSp[i1]<<endl;}
			}
	}  Info<<"Hdrn read ";
	//---- recharge conditions		
	icount = 0; //float dv=0;
	if (cellsHrch.size()>0)
	{
		float tnow = cellsHrchData[2]*day;
		while (time>=tnow) {icount += cellsHrch.size(); tnow = cellsHrchData[2+icount*4]*day; }
		if (icount>0) {icount -= cellsHrch.size();} Info<<" Hdrch size "<<cellsHrch.size();
		for (int i=0; i<cellsHrch.size(); i++) 
			{
			i1 = cellsHrch[i];//Info<<" i "<<i<<" i1 "<<i1<<" nla "<<nlay<<" ncl "<<ncell_lay;
			for (ilay=0;ilay<nlay;ilay++) { nc=i1-ilay*ncell_lay; if (sw[nc]>=sw_max) {break;} ;} //search for the first confined layer from the top
			ilay1 = max(0, ilay-1);// ilay1 from top
			nc = i1-ilay1*ncell_lay;clayRch[i] = nc;qHrch[i] = cellsHrchData[2+(icount+i)*4+2]/day;
			//if(hSu[nc]!=0) {Info<<"rch i i1 "<<i<<" "<<i1<<" nc "<<nc<<" h(nc) "<<h()[nc]<<" sw[nc] "<<sw[nc]<<" sw_min "<<sw_min[nc]<<" data "<<cellsHrchData[2+(icount+i)*4+2]<<endl;}
			if (hSu[nc]==0) {
				hSu[nc] = qHrch[i]/mesh.V()[nc]; 
				} // injeciton of recharge in the lowest unconfined layer only if no fix cond here
			} 	
	} Info<<"Hrch read ";//<<" Hghb "<<hSu[ncell_lay*3+ic0]<<" "<<endl;

// ***************** solving

fvScalarMatrix hEqn
	( stor*fvm::ddt(h) - fvm::laplacian(Mf,h) - hSu + fvm::SuSp(hSp,h) ); // + fvc::div(phiG) == b*fvOptions(sw,h) 
//fixed head cells
if (cellsHfx.size()>0)
{
	float tnow = cellsHfxData[2]*day; 
	while (time>=tnow) {icount += cellsHfx.size(); tnow = cellsHfxData[2+icount*4]*day; }
	if (icount>0) {icount -= cellsHfx.size();} 
	for (int i=0; i<cellsHfx.size(); i++) {hEqn.setReference(cellsHfx[i],cellsHfxData[2+(icount+i)*4+2]);}
}

hEqn.solve(); 

//phiw= -((Mf*gradHf) & mesh.Sf());
phiw == (-Mf * fvc::snGrad(h)) * mesh.magSf(); //similar to the above one
Uw = fvc::reconstruct(phiw);
//phiw =  hEqn.flux() ;//Info<<"phi "<<phiw<<endl;
volScalarField q1 = fvc::div(phiw);
volScalarField d1 = stor*fvc::ddt(h);
for (ilay=0;ilay<nlay;ilay++) {
	i1=ilay*ncell_lay+ic0;
	Info<<"h "<<h()[i1]<<" v "<<mesh.V()[i1]<<" sw "<<sw[i1]<<" q1 "<<q1[i1]<<" d1 "<<d1[i1]<<" SuSp "<<hSu[i1]-h()[i1]*hSp[i1]<<endl;
	}
// get the the water discharge/recharge at fixed heads cells
int i1;
for (int i=0; i<cellsHdrn.size(); i++) {i1=cellsHdrn[i];qHdrn[i]=(hSu[i1]-h()[i1]*hSp[i1])*mesh.V()[i1];} //Info<<"qhfx "<<i<<" "<<i1<<" "<<qHfx[i]<<endl;}
for (int i=0; i<cellsHghb.size(); i++) {i1=cellsHghb[i];qHghb[i]=(hSu[i1]-h()[i1]*hSp[i1])*mesh.V()[i1];qGerr[i]=qHghb[i]-(q1[i1]+d1[i1]);}// Info<<"qhghb "<<i<<" "<<i1<<" "<<qHghb[i]<<" q1 "<<q1[i1]*mesh.V()[i1]<<endl;}
//for (int i=0; i<cellsHghb.size(); i++) {i1=cellsHghb[i];qHghb[i]=q1[i1]*mesh.V()[i1];}// Info<<"qhghb "<<i<<" "<<i1<<" "<<qHghb[i]<<" q1 "<<q1[i1]*mesh.V()[i1]<<" qdH "<< (hSu[i1]-h()[i1]*hSp[i1])*mesh.V()[i1]<<endl;}
//for (int i=0; i<cellsHrch.size(); i++) {i1=clayRch[i];qHrch[i]=hSu[i1]*mesh.V()[i1];}// Info<<"qrch "<<i<<" "<<i1<<" "<<qHrch[i]<<endl;}
//for (int i=0; i<cellsHrch.size(); i++) {i1=clayRch[i];qHrch[i]=(q1[i1]-d1[i1])*mesh.V()[i1];}// Info<<"qrch "<<i<<" "<<i1<<" "<<qHrch[i]<<endl;}
//surfaceVectorField gradHf =fvc::interpolate(fvc::grad(h));

int ic;float q2;int iface;

/*for (ilay=1;ilay<nlay-1;ilay++) {
	ic=ilay*ncell_lay+ic0;q2=0;
	forAll(mesh.cells()[ic],i) {
		iface=mesh.cells()[ic][i]; oCell=mesh.owner()[iface];nCell=mesh.neighbour()[iface]; 
		Info<<"o "<<oCell<<" n "<<nCell<<" ifc "<<iface<<" phi "<<phiw[iface]<<endl;
		if (oCell==nc) {q2 += phiw[iface];}
		else {q2 -= phiw[iface];}
		}
	Info<<"h cell "<<ic<<" h "<<h()[ic]<<" sw "<<sw[ic]<<" q1 "<<q1[ic]<<" q2 "<<q2<<" q2/v "<<q2/mesh.V()[ic]<<endl;
	}*/

/* nc = (nlay-2)*ncell_lay + ic0;double sPhi=0;
forAll( mesh.owner() , iface ) {
	oCell=mesh.owner()[iface];nCell=mesh.neighbour()[iface];
	if ((oCell==nc)||(nCell==nc)){ 
		if (oCell==nc) {sPhi+=phiw[iface];}
		else {sPhi-=phiw[iface];}
		Info<<" o "<<oCell<<" n "<<nCell<<" Mf "<<Mf[iface]<<" phi "<<phiw[iface]<<" ho "<<h()[oCell]<<" hn "<<h()[nCell]<<" S "<<mesh.magSf()[iface]<<endl;  
		}
	}
Info<<"sphi "<<sPhi<<endl; */

/*
nc = (nlay-1)*ncell_lay + ic0;
forAll( mesh.owner() , iface ) {
	oCell=mesh.owner()[iface];nCell=mesh.neighbour()[iface];
	if ((oCell==nc)||(nCell==nc)){ Info<<" o "<<oCell<<" n "<<nCell<<" Mf "<<Mf[iface]<<" phi "<<phiw[iface]<<" ho "<<h()[oCell]<<" hn "<<h()[nCell]<<" S "<<mesh.magSf()[iface]<<endl;  }
	}
nc = (nlay-2)*ncell_lay + ic0;
forAll( mesh.owner() , iface ) {
	oCell=mesh.owner()[iface];nCell=mesh.neighbour()[iface];
	if ((oCell==nc)||(nCell==nc)){ Info<<" o "<<oCell<<" n "<<nCell<<" Mf "<<Mf[iface]<<" phi "<<phiw[iface]<<" ho "<<h()[oCell]<<" hn "<<h()[nCell]<<" S "<<mesh.magSf()[iface]<<endl;  }
	}
nc = (nlay-3)*ncell_lay + ic0;
forAll( mesh.owner() , iface ) {
	oCell=mesh.owner()[iface];nCell=mesh.neighbour()[iface];
	if ((oCell==nc)||(nCell==nc)){ Info<<" o "<<oCell<<" n "<<nCell<<" Mf "<<Mf[iface]<<" phi "<<phiw[iface]<<" ho "<<h()[oCell]<<" hn "<<h()[nCell]<<" S "<<mesh.magSf()[iface]<<endl;  }
	}
*/

for (i=0;i<cellsHfx.size();i++) {
	i1=cellsHfx[i];qHfx[i] = (q1[i1]-d1[i1])*mesh.V()[i1];//Info<<"qhfx "<<i<<" "<<qHfx[i]<<endl;//scalar q1 = fvc::surfaceIntegrate(Uw(), mesh.Sf(),0)[i1];
	}
	/*qHfx[i]=0; //ci-dessous, donne la mm chose que div
	forAll(mesh.cells()[i1],j) {
		iface=mesh.cells()[i1][j]; oCell=mesh.owner()[iface];nCell=mesh.neighbour()[iface]; 
		//Info<<"o "<<oCell<<" n "<<nCell<<" ifc "<<iface<<" phi "<<phiw[iface]<<endl;
		if (oCell==i1) {qHfx[i] -= phiw[iface];}
		else {qHfx[i] += phiw[iface];}
		}
	//qHfx[i] *= day; 
	Info<<"i "<<i<<" i1 "<<i1<<" q1 "<<q1[i1]*mesh.V()[i1]<<" qHfx "<<qHfx[i]<<endl;
	}//find the flux at hfix cells */
/*for (int ic=0; ic<cellsHfx.size(); ic++) 
	{
	const cell& faces = mesh.cells()[cellsHfx[ic]];
    forAll( faces, i ) {
		Info<<"cell "<<cellsHfx[ic]<<" face "<<faces[i]<<" phi "<<phiw[faces[i]]<<endl;
		if (mesh.owner()[faces[i]] == cellsHfx[ic]) {Info<<"its owner "<<endl;}
		if (mesh.neighbour()[faces[i]] == cellsHfx[ic]) {Info<<"its nghb "<<endl;}
		}
	}*/

