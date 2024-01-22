//general transport class
class plugin_Cwi
{
public:
	
	std::vector<double>Koci; //this is Koc for each species, read from a file
	std::vector<double>AW_a;std::vector<double>AW_b; // a in umol/cm3 or mol/m3 but data in muflow mol/L
	std::vector<double>lbdw_a;std::vector<double>lbdw_b;
	int sorpAW;
	void init(std::string & cur_dir,IOdictionary & Tprop,fvMesh & mesh,my_phq & freak)
	{
		scalar activateSorptionAW(Tprop.lookupOrDefault<scalar>("activateSorptionAW",0));
		this->sorpAW=activateSorptionAW;std::cout<<"in plugin "<<this->sorpAW<<"\n";
		std::string fname = cur_dir+"/constant/options/sorptionDecay";//std::cout<<"fname "<<fname<<"\n"; // getting Koc
		if (fexists(fname))
			{
			outTable table = readTable(fname);int ncol = table.ncol;
			this->Koci.resize(freak.compn.size(),0);
			this->lbdw_a.resize(freak.compn.size(),0);this->lbdw_b.resize(freak.compn.size(),0);
			this->AW_a.resize(freak.compn.size(),0);this->AW_b.resize(freak.compn.size(),0);
			//std::cout<<"ncol "<<ncol<<" comp4 "<<freak.compn[4]<<"\n";
			for (i=0;i<table.nrow;i++) {
				//std::cout<<"head "<<table.headers[i];
				int ic = std::find(freak.compn.begin(),freak.compn.end(),table.headers[i])-freak.compn.begin(); //find the position of str in compn
				this->Koci[ic]=table.data[i*ncol];this->lbdw_a[ic]=table.data[i*ncol+1];this->lbdw_b[ic]=table.data[i*ncol+2];
				if (activateSorptionAW==1)
				{this->AW_a[ic]=table.data[i*ncol+3];this->AW_b[ic]=table.data[i*ncol+4];		}
				}
			}
	}
	volScalarField run(fvMesh & mesh,int & ic, volScalarField & Cwi)
	{
		//calculate the value of Rf depending on Aw and thus sw
		double R=8.31; //R thermo constant in J.mol-1.K-1
		volScalarField eps = mesh.objectRegistry::lookupObject<volScalarField>("eps");
		volScalarField rhos = mesh.objectRegistry::lookupObject<volScalarField>("rhos");
		volScalarField foc = mesh.objectRegistry::lookupObject<volScalarField>("foc");
		volScalarField thetaw = mesh.objectRegistry::lookupObject<volScalarField>("thetaw");
		volScalarField Rf = mesh.objectRegistry::lookupObject<volScalarField>("Rf");
		if (this->Koci.size()>0) 
			{
			dimensionedScalar Koc1("Koc1",dimVol/dimMass,0);Koc1.value()=this->Koci[ic];
			std::cout<<"in plug ic "<<ic<<" koc "<<Koc1.value()<<"\n";
			Rf = 1 + (1-eps)/eps*rhos*Koc1*foc;
		
			if ((this->sorpAW==1)&&(this->AW_b[ic]>0))
				{
				//std::cout<<" a "<< this->AW_a[ic]<<" b "<<this->AW_b[ic] <<
				//std::cout<<"rhos "<<rhos[0]<<"\n";
				volScalarField dgrain = mesh.objectRegistry::lookupObject<volScalarField>("dgrain");
				volScalarField sw = mesh.objectRegistry::lookupObject<volScalarField>("sw");
				volScalarField T = mesh.objectRegistry::lookupObject<volScalarField>("T");
				Cwi.dimensions().reset(dimless); // *1000 later on from mol/L to mol/m3 (don't do it here it does not wrok!!)
				T.dimensions().reset(dimless);//T no dim but in °C
				// std::cout<<"in plug Rf0 "<<Rf[0]<<"\n";
				//Rf = 600.*(1-eps)/dgrain*(1-sw)*
				Rf = 1e-4*(548*pow(sw,2.) - 1182*sw + 634.); // Aaw from Guo Accusand 1e-4 to pass from cm2 to m2
				Rf *= this->AW_b[ic]/R*(7.572e-2 - 1.515e-4*T)/(T+273.15)/(this->AW_a[ic]+Cwi*1000.+SMALL); //Surf tension in N/m, T in °C, Cw in mol/m3
				Rf += rhos*(1-eps)*Koc1*foc;
				Rf = 1+Rf/thetaw;
				Cwi.dimensions().reset(dimMass/dimVol);
				T.dimensions().reset(dimensionSet(0,0,0,1,0,0,0));
				}
			}
		//std::cout<<"in plug Rf0 "<<Rf[0]<<"\n";
		return Rf;
	}
};