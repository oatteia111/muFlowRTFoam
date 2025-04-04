/*This plugin will modify some elements before the H equation
including transformation of muw, rhow and k
the equations are read from a unique file called coupling and what is sent back is a multiplying factor 
for the given variable
the file coupling is read to find variables of interest (Y) and the predicting variable (X)
this is stored a series of Y=f(X) fnctions that will be applied when needed
*/
class plugin_H
{
public:
	//volScalarField factM; don't put that here, it creates pb 
	std::vector<dicFunc> listCouples;
	
	void init(std::string & cur_dir,IOdictionary & Tprop,fvMesh & mesh,Time & runTime,std::vector<dicFunc> & listCouples)
	{
		std::cout<<"init plugin H "<<" ";
		this->listCouples = listCouples;
		//volScalarField::Internal deltaMvol(IOobject("deltaMvol",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),mesh,0.0);
		//volScalarField reps0koz(IOobject("reps0koz",runTime.constant(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),mesh,0.0);
		//volScalarField eps = mesh.objectRegistry::lookupObject<volScalarField>("eps");  //Xvar is the name of the variable
		//reps0koz = pow(eps,3)/pow(1-eps,2);
	}
	
	volScalarField corrM(fvMesh & mesh,my_phq & freak) //this may modify muw, rhow due to effects of C or T, and K due to minerals (Mf=rhow*mag(g)*Kf/muw)
	{
		volScalarField fact=mesh.objectRegistry::lookupObject<volScalarField>("unity");
		
		for (size_t i=0;i<listCouples.size();i++) {
			volScalarField X = mesh.objectRegistry::lookupObject<volScalarField>(listCouples[i].Xvar);  //Xvar is the name of the variable
			if (listCouples[i].Xvar=="eps") 
				{
				volScalarField eps0 = mesh.objectRegistry::lookupObject<volScalarField>("eps0");  
				volScalarField eps = mesh.objectRegistry::lookupObject<volScalarField>("eps");  
				forAll(eps(),i){eps[i]=eps0[i]-freak.deltaMvol[i]*0.037;} // !! this is valid only for calcite //if (i<5) {Info<<"dm "<<freak.deltaMvol[i]<<endl;}
				X = eps*1;
				} 
			Info<<" X name "<<listCouples[i].Xvar<<" X(0) "<<X[0]<<endl;
			if ((listCouples[i].Ykey=="rhow")||(listCouples[i].Ykey=="K"))
				{fact *= calcFunc(mesh,listCouples[i].type,X,listCouples[i].Xref,listCouples[i].fparms);}
			if (listCouples[i].Ykey=="muw") 
				{fact /= calcFunc(mesh,listCouples[i].type,X,listCouples[i].Xref,listCouples[i].fparms);}
		}
		
		return fact;
	}

};