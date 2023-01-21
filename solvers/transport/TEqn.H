
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
int icount = 0;float day = 86400.; //Info <<"icount "<< icount << endl;//counter of lines within poption file
int nlay = int(cellsHwData[0]);int ncell_lay = int(cellsHwData[1]);float q=0;
float tnow = cellsHwData[2]*day; //time in files is in days
std::vector<double> a(cellsHw.size(),0.); // the cells in Hw and Tw must be the same!!!
std::vector<double> aAdd(cellsHw.size(),0.); 
while (time>=tnow)
{
	icount += cellsHw.size(); 
	tnow = cellsHwData[2+icount*4]*day; 
}
if (icount>0) {icount -= cellsHw.size();};
for (int i=0; i<cellsHw.size(); i++) {aAdd[i] = cellsHwData[2+(icount+i)*4+2]/day;} // discharge must be in s-1

TSu*=0; TSp*=0; //previous cells with vlaues msut be destroyed
while (time>=tnow) // this is done twice because the timing in hwel ad twel files can be different
{
	icount += cellsTw.size(); 
	tnow = cellsTwData[2+icount*4]*day; 
}
if (icount>0) {icount -= cellsHw.size();};
for (int i=0; i<cellsHw.size(); i++) {a[i] = cellsTwData[2+(icount+i)*4+2];}

for (int i=0; i<cellsHw.size(); i++) 
	{ 
	q = aAdd[i]/mesh.V()[cellsHw[i]]; // dataAdd is the rate (from hwel here without dimension)
	scalar tloc = mag(T()[cellsHw[i]]);
	if (q<0) {TSp[cellsHw[i]] = -q;} //in Sp term as discharge prop to local concentration
	else {TSu[cellsHw[i]] = -q*a[i];} //the injected heat is the water flow x the temp
	Info<<"cell "<<cellsHw[i]<<" q "<<q<<" tloc "<<tloc<<" Su "<<TSu[cellsHw[i]]<<" Sp "<<TSp[cellsHw[i]]<<endl;
	}

Info << "before Teq "<<endl;
fvScalarMatrix TEqn
	(
		rhow*cpw * thetaw * fvm::ddt(T) + rhow*cpw * T * fvc::ddt(thetaw) + rhos*cps * (1-eps) * fvm::ddt(T) 
		+ fvm::div(rhow*cpw * phiw, T, "div(phiw,T)") 
		- fvm::laplacian(rhow*cpw * Deffw, T, "laplacian(Deffw,T)")
		- fvm::laplacian(lbdaT, T, "laplacian(lbdaT,T)")
		+ rhow*cpw * TSu + rhow*cpw * fvm::SuSp(TSp, T)
	);
/*
volScalarField cT = rhow*cpw*T;
volScalarField rCp = rhos*cps/rhow/cpw;
volScalarField lbdaT1 = lbdaT/rhow/cpw;
fvScalarMatrix TEqn
	(
		thetaw * fvm::ddt(cT) + cT * fvc::ddt(thetaw) + rCp *(1-eps) * fvm::ddt(cT) 
		+ fvm::div(phiw, cT, "div(phiw,T)") 
		- fvm::laplacian(Deffw, cT, "laplacian(Deffw,T)")
		- fvm::laplacian(lbdaT1, cT, "laplacian(lbdaT,T)")
		==
		fvOptions(sw,cT) 
	);*/

fvOptions.constrain(TEqn);
TEqn.solve();


