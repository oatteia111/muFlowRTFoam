dimensionedScalar USMALL("",dimLength/dimTime,1e-36);
thetaw = eps*sw;
//Deffw = (alphaT * mag(Ut))*tensor::I + (alphaL-alphaT) * (Ut*Ut)/(mag(Ut)+USMALL);
Deffw = (alphaT * mag(Uw))*tensor::I + (alphaL-alphaT) * (Uw*Uw)/(mag(Uw)+USMALL) + Dw0 * pow(thetaw,3.33)/pow(eps,2.)*tensor::I;
//Deffw = Dw0 * pow(thetaw,3.33)/pow(eps,2.)*tensor::I;

for (i=1; i<ph_ncomp;i++)
    {
		//Info << "in CiEqn "<< i << endl;
		volScalarField& Cwi = Cw[i]();
		word nameCi ("Cw" + name(i));
		
		fvScalarMatrix CwiEqn
		(		  
            thetaw * fvm::ddt(Cwi) + Cwi * fvc::ddt(thetaw)
            + fvm::div(phiw, Cwi, "div(phiw,Cwi)")
            - fvm::laplacian(Deffw, Cwi, "laplacian(Deffw,Cwi)")
			==
			fvOptions(sw,Cwi)
		);

		fvOptions.constrain(CwiEqn);
		CwiEqn.solve();
    }
