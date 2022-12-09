thetag = max(eps*(1 - sw),SMALL);
if (diffusionEqn==0) {Deffg = Dg0 * pow(thetag,3.333) /pow(eps,2.)*tensor::I;}
if (diffusionEqn==1) {Deffg = Dg0 * pow(thetag,2.333) /pow(eps,2.)*tensor::I;}
Info<<"sw "<<sw[2]<<" thetag "<<thetag[2]<<" Deffg 2 "<<Deffg[2]<<endl;

for (i=0; i<ph_gcomp;i++)
    {
		//Info << "in CgiEqn "<< i << endl;
		volScalarField& Cgi = Cg[i]();
		word nameCi ("Cg" + name(i));
		
		fvScalarMatrix CgiEqn
		(		  
            thetag * fvm::ddt(Cgi) +  Cgi * fvc::ddt(thetag)
			+ fvm::div(phig, Cgi, "div(phig,Cgi)") // used only in multiphase
            - thetag * fvm::laplacian(Deffg, Cgi, "laplacian(Deffg,Cgi)")
			==
			fvOptions(sw,Cgi)
		);

		fvOptions.constrain(CgiEqn);
		CgiEqn.solve();
    }
