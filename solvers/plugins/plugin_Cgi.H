//This plugin will modify some elements of the Cgi equation
class plugin_Cgi
{
public:
	
	std::vector<double>lbdagi; //this is the deg rate in gas for each species, read from a file
	int flag;
	void init(std::string & cur_dir,IOdictionary & Tprop,fvMesh & mesh,my_phq & freak,std::vector<dicFunc> & listCouples)
	{
		//std::cout<<"in plugin Cgi "<<freak.gcomp.size()<<" ";
		for (size_t i=0;i<listCouples.size();i++) {
			if (listCouples[i].Ykey=="Deffg") {this->flag=1;}
		}

		if (ph_gcomp>0) {this->lbdagi.resize(freak.gcomp.size(),0.);std::cout<<this->lbdagi[0]<<"\n";}
		std::string fname = cur_dir+"/constant/options/gasDecay";std::cout<<"fname "<<fname<<"\n"; // getting Koc
		if (fexists(fname))
			{
			outTable table = readTable(fname);int ncol = table.ncol;std::cout<<"gcomp size "<<freak.gcomp.size()<<" ncol "<<ncol<<'\n';
			for (i=0;i<table.nrow;i++) {
				int ic = std::find(freak.gcomp.begin(),freak.gcomp.end(),table.headers[i])-freak.gcomp.begin(); //find the position of str in compn
				this->lbdagi[ic]=table.data[i*ncol];std::cout<<"head "<<table.headers[i]<<" "<<this->lbdagi[ic]<<"\n";
				}
			}
		std::cout<<"end plug Cg\n";
	}
	float getLbdag(int & ic) //fvMesh & mesh,int & ic, volScalarField & Cgi)
	{
		//std::cout<<"in plugin ic "<<ic<<" lbd "<<lbdagi[ic]<<" \n";
		return this->lbdagi[ic];
	}
};