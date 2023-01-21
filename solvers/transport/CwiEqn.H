dimensionedScalar USMALL("",dimLength/dimTime,1e-36);
thetaw = eps*sw;
//Deffw = (alphaT * mag(Ut))*tensor::I + (alphaL-alphaT) * (Ut*Ut)/(mag(Ut)+USMALL);
Deffw = (alphaT * mag(Uw))*tensor::I + (alphaL-alphaT) * (Uw*Uw)/(mag(Uw)+USMALL) + Dw0 * pow(thetaw,3.33)/pow(eps,2.)*tensor::I;
//Deffw = Dw0 * pow(thetaw,3.33)/pow(eps,2.)*tensor::I;

double dC=0;double dC0,mxC,sC,sC1;

for (i=0; i<ph_ncomp;i++)
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
		//removing negative and compensating + dC
		if (i>3) //not used for water, O, H,
		{
			dC0=0;sC=0;sC1=0;mxC=0;
			for (j=0;j<ncell;j++) {mxC=max(mxC,Cwi()[j]);sC+=Cwi()[j];}
			Cwi=max(Cwi,Cwi*0); // removing negative
			for (j=0;j<ncell;j++) {dC0= max(dC0,Cwi()[j]-Cwi.prevIter()[j]);sC1+=Cwi()[j];}
			Cwi *= sC/(sC1+SMALL); // removing negative value increase artifially the total
			dC0 = dC0/(mxC+SMALL); // dC is relative
			dC = max(dC,dC0); // we consider the species that leads to the max dC
			//Info <<"comp "<<i<<" dC0 "<<dC0<<" dC "<<dC<<endl;
		}
	}
