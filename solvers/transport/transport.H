//general transport class
class my_transport
{
public:
	
	void init(fvMesh & mesh)
	{
	volScalarField eps = mesh.objectRegistry::lookupObject<volScalarField>("eps");
	volScalarField rhos = mesh.objectRegistry::lookupObject<volScalarField>("rhos");
	volScalarField foc = mesh.objectRegistry::lookupObject<volScalarField>("foc");
	//dimensionedScalar Koc = mesh.objectRegistry::lookupObject<dimensionedScalar>("Koc");
	//volScalarField unity = mesh.objectRegistry::lookupObject<volScalarField>("unity");
	volScalarField Rf = mesh.objectRegistry::lookupObject<volScalarField>("Rf");
	//volScalarField Rf=unity;
	//Rf = 1 + (1-eps)/eps*rhos*foc;//Info<< " dim rhos"<<rhos.dimensions()<<" Rf "<<Rf.dimensions()<< endl;
	}
};