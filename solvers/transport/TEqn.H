
// diffusion+dispersion
dimensionedScalar USMALL("",dimLength/dimTime,1e-36);
//Deffw = (alphaT * mag(vp))*tensor::I + (alphaL-alphaT) * (vp*vp)/(mag(vp)+USMALL); // Ut rpelaced by vp
thetaw = eps*sw;
Deffw = (alphaT * mag(Uw) + Dw0* pow(thetaw,3.333) /pow(eps,2.))*tensor::I + (alphaL-alphaT) * (Uw*Uw)/(mag(Uw)+USMALL);//Info<<"Deffw "<<Deffw<<endl;

forAll(Deffw.boundaryField(),patchi)
	{ //- Dispersion equal to zero on patch where flux is fixed
	if (C.boundaryField().types()[patchi] == "fixedFlux") Deffw.boundaryFieldRef()[patchi] = tensor(zero());
	}
lbdaT = thetaw*lbdaTw + (1-eps)*lbdaTs;

// prepare the internal conditions
float time = mesh.time().value(); //const Time&
int nlay = int(cellsHwData[0]);int ncell_lay = int(cellsHwData[1]);float q=0;
float tnow = cellsHwData[2]*day; //time in files is in days

//get the temperature well data
while (time>=tnow) // this is done here and in heqn because the timing in hwel ad twel files can be different
{
	icount += cellsTw.size(); 
	tnow = cellsTwData[2+icount*4]*day; 
}
if (icount>0) {icount -= cellsHw.size();};
for (int i=0; i<cellsHw.size(); i++) {a[i] = cellsTwData[2+(icount+i)*4+2];}

// set the values in the matrix
TSu*=0; TSp*=0; //previous cells with values must be destroyed
for (int i=0; i<cellsHw.size(); i++) 
	{ 
	q = condHw[i]/mesh.V()[cellsHw[i]]; // dataAdd is the rate (from hwel here without dimension)
	scalar tloc = mag(T()[cellsHw[i]]);
	if (q<0) {TSp[cellsHw[i]] = -q;} //in Sp term as discharge prop to local concentration
	else {TSu[cellsHw[i]] = -q*a[i];} //the injected heat is the water flow x the temp
	Info<<"cell "<<cellsHw[i]<<" q "<<q<<" tloc "<<tloc<<" Su "<<TSu[cellsHw[i]]<<" Sp "<<TSp[cellsHw[i]]<<endl;
	}
	
// find the place where there is water leaving the domain (e.g. fixed head cells)
//volScalarField dvPhi = fvc::div(phiw, C, "div(phiw,C)"); oes not work, it provides null div at cell with fix h
dffPhi *=0;

//for (int i=0;i<3;i++) {Info<<" Kh "<<Kh[i]<<" Kf "<<Kf[i]<<" Mf "<<Mf[i]<<" grdh "<<fvc::grad(h)()[i]<<" surf "<<mesh.Sf()[i]<<endl;}

//forAll(mesh.neighbour(),iface ) {dffPhi[mesh.neighbour()[iface]] -= phiw[iface];Info<< "i fc "<<iface<<" ng "<<mesh.neighbour()[iface]<<" Kf "<<Kf[iface]<<" Mf "<<Mf[iface]<<" dffP "<<dffPhi[mesh.neighbour()[iface]]<<endl;}
//forAll(mesh.owner(),iface ) {dffPhi[mesh.owner()[iface]] += phiw[iface];Info<< "i fc "<<iface<<" "<<mesh.owner()[iface]<<" dffPh "<< dffPhi[mesh.owner()[iface]]<<endl;}
//Info<<"dffPhi "<<dffPhi<<endl;
/*forAll(mesh.cells(), celli)
    {
		const labelList& cFaces = mesh.cells()[celli];
        forAll(cFaces, cFacei){
			Info<<"cell "<<celli<<" faces "<<cFaces[cFacei]<<" phiw "<<phiw[cFaces[cFacei]]<<" "<<mesh.neighbour()[cFaces[cFacei]]<<endl;
			} 
	}*/
//this is the flux divergence
//forAll(dvPhi,i) {if (dvPhi()[i]>1e-3) {TSp[i] = -dvPhi()[i];Info<<"div pos at "<<i<<" dv "<<dvPhi()[i]<<endl;} } //if the div is too positive there is imposed discharge here

T.storePrevIter();

Info << "before Teq "<<endl;
fvScalarMatrix TEqn
	(
		rhow*cpw * thetaw * fvm::ddt(T) + rhow*cpw * T * fvc::ddt(thetaw) + rhos*cps * (1-eps) * fvm::ddt(T) 
		+ fvm::div(rhow*cpw * phiw, T, "div(phiw,T)") 
		- fvm::laplacian(eps*rhow*cpw * Deffw, T, "laplacian(Deffw,T)")
		- fvm::laplacian(lbdaT, T, "laplacian(lbdaT,T)")
		+ rhow*cpw * TSu + rhow*cpw * fvm::SuSp(TSp, T)
	);

//fvOptions.constrain(TEqn);
TEqn.solve();

//calculating dT
for (i=0;i<ncell;i++) {dT= max(dT,T()[i]-T.prevIter()[i]);}



