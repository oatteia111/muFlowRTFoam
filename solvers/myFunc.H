//this is a class to store func properties
class dicFunc
{
	public:
	std::string Ykey; void setYkey(std::string Ykey){this->Ykey = Ykey;}	
	std::string Xvar; void setXvar(std::string Xvar){this->Xvar = Xvar;}
	std::string type; void setType(std::string type){this->type = type;}
	float Xref; void setXref(float Xref){this->Xref = Xref;}
	std::vector<float> fparms; void setParms(std::vector<float> fparms){this->fparms = fparms;}
};

//this is the calculation of the value of the multiplying factor from a function
// it returns a vector of multiplygin factors given the value of fvar
// the first paramaters is the reference
volScalarField calcFunc(fvMesh & mesh, std::string ftype, Foam::volScalarField fvar, float Xref, std::vector<float> fparms)
{
	//volScalarField resu; 
	volScalarField fact = fvar*0.0;
	fact.dimensions().reset(dimless);
	fvar.dimensions().reset(dimless);
		/*int i;
		int n = fvar.size();//Info<<"in func n "<<n<<" type "<<ftype<<endl;
		resu.resize(n);*/
		if (ftype == "linear")
			{
			//for (i=0;i<n;i++) {resu[i] = 1+ fparms[1]*(fvar()[i]-fparms[0]);}
			fact = 1.+ fparms[0]*(fvar-Xref);
			Info<<"in calc func, fvar "<<fvar[4]<<" fparm "<<fparms[0]<<" fact(4) "<<fact[4]<<endl;
			}
		if (ftype == "polynomial")
			{
			//for (i=0;i<n;i++) {resu[i] = 1+ fparms[1]*(fvar()[i]-fparms[0]);}
			fact = 1.+ fparms[0]*(fvar-Xref) + fparms[1]*pow(fvar-Xref,2) ;
			Info<<"in calc func, fvar "<<fvar[0]<<" fparm "<<fparms[0]<<" fact(4) "<<fact[4]<<endl;
			}
		if (ftype == "exponential")
			{
			//for (i=0;i<n;i++) {resu[i] = 1+ fparms[1]*(fvar()[i]-fparms[0]);}
			fact = 1.+ fparms[0]*exp(fparms[1]*(fvar-Xref)) ;
			Info<<"in calc func, fvar "<<fvar[0]<<" fparm "<<fparms[0]<<" fact(4) "<<fact[4]<<endl;
			}
		if (ftype=="kozeny") // kozeny F=eps^3/(1-eps^2), with the initial being 
			{
			volScalarField eps0 = mesh.objectRegistry::lookupObject<volScalarField>("eps0");  //Xvar is the name of the variable
			fact = (pow(fvar,3)/pow(1-fvar,2)) / (pow(eps0,3)/pow(1-eps0,2));
			Info<<"in func koz "<<fvar[0]<<" fact "<<fact[0]<<endl;
			}
		return fact;
};

volScalarField calcFunc2(fvMesh & mesh, std::string ftype, Foam::volScalarField C, Foam::volScalarField U,float Xref, std::vector<float> fparms,float tunits)
{
	//here we use two variables and fixed factors to start; 
	volScalarField fact = C*0.0;
	fact.dimensions().reset(dimless);
	C.dimensions().reset(dimless);
	U.dimensions().reset(dimless);
	if (ftype == "T0_muInf") // antonio
		{
		//C is conc and U is velocity (for reactive case); k 2.2e-11; n 0.3 sqrt = 2.57e-6
		volScalarField gamm = 2.4* max(max(U/tunits,0.),1e-12)/2.57e-6; 
		volScalarField T0 = 0.367*max(C-0.483*Xref,0.);
		volScalarField kC = 4.54e-5+3.26e-3*C/Xref + 5.97e-2*pow(C/Xref,2);
		volScalarField nC = 0.6127-0.1608*log(1+max(C,1e-20)/Xref/6.27);
		volScalarField muInf = 9.76e-4+1.6e-3*C/Xref;
		//fact = pow(max(C/Xref,0.), 1.5)*(0.06+pow(max(gamm/2.5,0.),-0.85))*400;
		fact = 1e3*(T0/gamm+kC*pow(gamm,nC-1)+muInf);//1e3 fro water visco
		fact = max(min(fact,1e4),1.);
		Info<<"in calc func2 Xref T0, kC, nC, muInf "<<Xref<<" "<<T0[4]<<" "<<kC[4]<<" "<<nC[4]<<" "<<muInf[4]<<endl;
		Info<<"C,U(wel) "<<C[4]<<" "<<U[4]<<" gam(wel) "<<gamm[4]<<" fac(wel) "<<fact[4]<<endl;
		}
	if (ftype == "T0_modif") // antonio modified
		{
		//C is conc and U is velocity (for reactive case); k 2.2e-11; n 0.3 sqrt = 2.57e-6
		volScalarField gamm = 2.4* max(max(U/tunits,0.),1e-12)/2.57e-6; 
		volScalarField T0 = fparms[0]*max(C-fparms[1]*Xref,0.);
		volScalarField kC = fparms[2] + fparms[3]*C/Xref + fparms[4]*pow(C/Xref,2);
		volScalarField nC = fparms[5] +C*0;//+ fparms[4] *log(1+max(C,1e-20)/Xref/fparms[5]);
		volScalarField muInf = fparms[6] + fparms[7]*C/Xref;
		//fact = pow(max(C/Xref,0.), 1.5)*(0.06+pow(max(gamm/2.5,0.),-0.85))*400;
		fact = 1e3*(T0/gamm+kC*pow(gamm,nC-1)+muInf);//1e3 fro water visco
		fact = max(min(fact,2e3),1.);
		Info<<"in calc func2 Xref T0, kC, nC, muInf "<<Xref<<" "<<T0[4]<<" "<<kC[4]<<" "<<nC[4]<<" "<<muInf[4]<<endl;
		Info<<"C,U(wel) "<<C[4]<<" "<<U[4]<<" gam(wel) "<<gamm[4]<<" fac(wel) "<<fact[4]<<endl;
		}
	if (ftype == "dble_power") // olivier
		{
		//C is conc and U is velocity (for reactive case); k 2.2e-11; n 0.3 sqrt = 2.57e-6
		volScalarField gamm = 2.4* max(max(U/tunits,0.),1e-12)/2.57e-6; 
		volScalarField fG = min(50.,0.06 + 0.4*pow(gamm/2.5,-0.85));
		fact = fG*pow(max(C/Xref/2,0.), 1.5)/0.001;
		fact = max(min(fact,1e4),1.);
		Info<<"in calc xref, fG "<<Xref<<" "<<fG[4]<<endl;
		Info<<"C,U(wel) "<<C[4]<<" "<<U[4]<<" gam(wel) "<<gamm[4]<<" fac(wel) "<<fact[4]<<endl;
		}
	return fact;
}
