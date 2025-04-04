// file reader of coupling
std::string fname = cur_dir+"/constant/options/coupling";
//std::ifstream f(fname.c_str());
std::vector<float> fparms {1.,0.};//std::cout<<"fparms "<<fparms<<"\n";
std::vector<dicFunc> listCouples;

if (fexists(fname))
	{
	// open the file
	std::ifstream if0(fname);//std::cout<<"file exist \n";
	// read by keywords
	std::string Ykey,Xvar,ftype;
	int nb;
	float Xref,aa;
	std::vector<float> fparms;
	std::string line;

	while(std::getline(if0,line))
		{
		dicFunc data;
		std::istringstream il0(line);
		il0 >> Ykey;data.setYkey(Ykey);// name of the y var
		il0 >> Xvar;data.setXvar(Xvar); // name of the x variable
		il0 >> Xref;data.setXref(Xref);// the reference value for x
		il0 >> ftype;data.setType(ftype); // type fo funciton
		il0 >> nb; // number of parameters
		std::cout<<"func "<<Ykey<<" "<<Xvar<<" "<<Xref<<" "<<ftype<<" "<<nb<<"\n";
		fparms.resize(nb);
		for (i=0;i<nb;i++) {il0 >> aa;fparms[i] = aa;}	
		data.setParms(fparms);
		listCouples.push_back(data);		
		}
	}
/*
else
	{
	fDe_T.setType("linear");
	fDe_T.setParms(fparms);
	}
	*/