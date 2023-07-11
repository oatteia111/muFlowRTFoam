
class dicFunc
{
	public:
	
	std::string type;
	void setType(std::string type){this->type = type;}
	std::vector<double> parms;
	void setParms(std::vector<double> parms){this->parms = parms;}
};


std::vector<double> myFunc(std::string ftype, Foam::volScalarField fvar, std::vector<double> fparms)
{
	std::vector<double> resu; 
		int i;
		int n = fvar.size();Info<<"in func n "<<n<<" type "<<ftype<<endl;
		resu.resize(n);
		if (ftype == "linear")
		{
			for (i=0;i<n;i++) {resu[i] = fparms[0]+fparms[1]*fvar()[i];}
		}
		return resu;
};
