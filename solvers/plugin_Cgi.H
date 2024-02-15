//This plugin will modify some elements of the Cgi equation
class plugin_Cgi
{
public:
	
	std::vector<double>lbdagi; //this is the deg rate in gas for each species, read from a file
	void init(std::string & cur_dir,IOdictionary & Tprop,fvMesh & mesh,my_phq & freak)
	{
		this->lbdagi.resize(freak.gcomp.size(),0);
		std::string fname = cur_dir+"/constant/options/gasDecay";std::cout<<"fname "<<fname<<"\n"; // getting Koc
		if (fexists(fname))
			{
			outTable table = readTable(fname);int ncol = table.ncol;std::cout<<"gcomp size "<<freak.gcomp.size()<<" ncol "<<ncol<<'\n';
			for (i=0;i<table.nrow;i++) {
				int ic = std::find(freak.gcomp.begin(),freak.gcomp.end(),table.headers[i])-freak.gcomp.begin(); //find the position of str in compn
				this->lbdagi[ic]=table.data[i*ncol];std::cout<<"head "<<table.headers[i]<<" "<<this->lbdagi[ic]<<"\n";
				}
			}
	}
	float setLbdag(fvMesh & mesh,int & ic, volScalarField & Cgi)
	{
		return this->lbdagi[ic];
	}
};