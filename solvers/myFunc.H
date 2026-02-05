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
// it retursn a vector of multiplygin factors given the value of fvar
// the frist paramaters is the reference
volScalarField calcFunc(fvMesh & mesh, std::string ftype, Foam::volScalarField fvar, float Xref, std::vector<float> fparms)
{
	//volScalarField resu; 
	volScalarField fact=fvar*0;
	fact.dimensions().reset(dimless);
	fvar.dimensions().reset(dimless);
		/*int i;
		int n = fvar.size();//Info<<"in func n "<<n<<" type "<<ftype<<endl;
		resu.resize(n);*/
		if (ftype == "linear")
			{
			//for (i=0;i<n;i++) {resu[i] = 1+ fparms[1]*(fvar()[i]-fparms[0]);}
			fact = 1.+ fparms[0]*(fvar-Xref);Info<<"in calc func, fvar "<<fvar[0]<<" fparm "<<fparms[0]<<" fact(5) "<<fact[5]<<endl;
			}
		if (ftype == "polynomial")
			{
			//for (i=0;i<n;i++) {resu[i] = 1+ fparms[1]*(fvar()[i]-fparms[0]);}
			fact = 1.+ fparms[0]*(fvar-Xref) + fparms[1]*pow(fvar-Xref,2) ;
			Info<<"in calc func, fvar "<<fvar[0]<<" fparm "<<fparms[0]<<" fact(5) "<<fact[5]<<endl;
			}
		if (ftype == "exponential")
			{
			//for (i=0;i<n;i++) {resu[i] = 1+ fparms[1]*(fvar()[i]-fparms[0]);}
			fact = 1.+ fparms[0]*exp(fparms[1]*(fvar-Xref)) ;
			Info<<"in calc func, fvar "<<fvar[0]<<" fparm "<<fparms[0]<<" fact(5) "<<fact[5]<<endl;
			}
		if (ftype=="kozeny") // kozeny F=eps^3/(1-eps^2), with the initial being 
			{
			volScalarField eps0 = mesh.objectRegistry::lookupObject<volScalarField>("eps0");  //Xvar is the name of the variable
			fact = (pow(fvar,3)/pow(1-fvar,2)) / (pow(eps0,3)/pow(1-eps0,2));
			Info<<"in func koz "<<fvar[0]<<" fact "<<fact[0]<<endl;
			}
		return fact;
};
