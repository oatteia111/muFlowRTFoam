/*This plugin will modify some elements before the pS equation
including transformation of muw, rhow and k or kr
the equations are read from a unique file called coupling and sends back a multiplying factor 
for the given variable
the file coupling is read to find variables of interest (Y) and the predicting variable (X)
this is stored a series of Y=f(X) fnctions that will be applied when needed
*/
class plugin_PS
{
public:
	//volScalarField factM; don't put that here, it creates pb 
	int flag,flgFoam,isurfact;
	float fmmob,epdry,fmdry,cref,fc; //,wfrac;
	
	void init(std::string & cur_dir,IOdictionary & Tprop,fvMesh & mesh,my_phq & freak,std::vector<dicFunc> & listCouples)
	{
		//this->flag=1;
		std::cout<<"init plugin PS "<<" ";
		// first read the foam parameters in the file constant/foam
		std::vector<float> foamParms;
		std::ifstream inputFoam{cur_dir+"/constant/options/foam" }; 
		foamParms = {std::istream_iterator<float>{inputFoam}, {}};
		this->flag = foamParms.size();
		if (this->flag>0) {
			this->fmmob=foamParms[0];
			this->epdry=foamParms[1];
			this->fmdry=foamParms[2];
			this->cref=foamParms[3];
			this->fc=foamParms[4];
			//this->wfrac=foamParms[5];
			}
		//get the number if the species Surfact
		this->isurfact = std::find(freak.comp.begin(),freak.comp.end(),"Surfact")-freak.comp.begin(); //find the position of str in comp
		std::cout<<"isurfact "<<this->isurfact<<"\n";

	}

};