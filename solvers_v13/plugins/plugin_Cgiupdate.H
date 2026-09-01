//change the values used for Cgi equation
for (size_t i=0;i<listCouples.size();i++) {
	volScalarField X = mesh.objectRegistry::lookupObject<volScalarField>(listCouples[i].Xvar);  //Xvar is the name of the variable
	Info<<" X name "<<listCouples[i].Xvar<<" X(0) "<<X[0]<<endl;
	if (listCouples[i].Ykey=="Deffg")
		{
		Deffg *= calcFunc(mesh,listCouples[i].type,X,listCouples[i].Xref,listCouples[i].fparms);
		Info<<"Deffg "<<Deffg[0]<<endl;
		}
}
