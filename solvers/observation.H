obsFlow.push_back(mesh.time().value());
obsTrans.push_back(mesh.time().value());
obsSolu.push_back(mesh.time().value());
obsGas.push_back(mesh.time().value());
for (int io=0;io<nobs;io++) 
	{
	for (int iv=0;iv<obsFlowIndx.size();iv++)
		{
		if (obsFlowIndx[iv]==0) {obsFlow.push_back(hp()[icello[io]]);}
		if (obsFlowIndx[iv]==1) {obsFlow.push_back(sw()[icello[io]]);Info<<"in obs store "<<sw()[icello[io]]<<endl;}
		}
	if (activateTransport==1)
		for (int iv=0;iv<obsTransIndx.size();iv++)
			{
			if (obsTransIndx[iv]==0) {obsTrans.push_back(C()[icello[io]]);Info<<"in obs store "<<C()[icello[io]]<<endl;}
			if (obsTransIndx[iv]==1) {obsTrans.push_back(T()[icello[io]]);}
			}
	if (obsChemIndx[0]==1)
		for (int ic=0;ic<ph_ncomp;ic++)
			{obsSolu.push_back(Cw[ic]()[icello[io]]);}
	if ((ph_gcomp>0)&&(obsChemIndx[1]==1))
		for (int ig=0;ig<ph_gcomp;ig++)
			{obsGas.push_back(Cg[ig]()[icello[io]]);}
	}
oindex += 1;std::cout<<"oindex "<<oindex<<std::endl;

std::vector<std::string> obsFlowNames = {"hp","sw"};
std::vector<std::string> obsTransNames = {"C","T"};
int nv;
if (flagW==1) // write the files    
	{
	Info<<"obs size obsFlow "<<obsFlow.size()<<" "<<obsFlow[0]<<" "<<obsFlow[1]<<" "<<obsFlow[2]<<" "<<obsFlow[3]<<endl;
	for (int io=0;io<nobs;io++) //one file for each ob point
		{
		nv=obsFlowIndx.size();
		Info<<"nobs "<<nobs<<" nv "<<nv<<endl;
		for (size_t iv=0;iv<nv;iv++)
			{
			Info <<"in owrite "<<io<<" "<<iv<<" "<<obsFlowIndx[iv]<<" "<<obsFlowNames[obsFlowIndx[iv]]<<endl;
			std::ofstream fich1(cur_dir+"/observation/obs_"+observ.headers[io]+"_"+obsFlowNames[obsFlowIndx[iv]]+".txt", std::ios::app);// Ouvrir le fichier en mode append
			for (i=0;i<oindex;i++) {fich1 << obsFlow[i*(nobs*nv+1)]/day<<" "<< obsFlow[i*(nobs*nv+1)+io*nv+iv+1] << "\n";} //
			fich1.close();
			}
		nv=obsTransIndx.size();
		for (size_t iv=0;iv<nv;iv++)
			{
			std::ofstream fich2(cur_dir+"/observation/obs_"+observ.headers[io]+"_"+obsTransNames[obsTransIndx[iv]]+".txt", std::ios::app);// Ouvrir le fichier en mode append
			for (i=0;i<oindex;i++) {fich2 << obsTrans[i*(nobs*nv+1)]/day<<" "<< obsTrans[i*(nobs*nv+1)+io*nv+iv+1] << "\n";}
			fich2.close();
			}
		int nc=ph_ncomp;
		for (int ic=0;ic<nc;ic++)
			{
			std::ofstream fich3(cur_dir+"/observation/obs_"+observ.headers[io]+"_Cw"+std::to_string(ic)+".txt", std::ios::app);// Ouvrir le fichier en mode append
			for (i=0;i<oindex;i++) {fich3 << obsSolu[i*(nobs*nc+1)]/day<<" "<< obsSolu[i*(nobs*nc+1)+io*nc+ic+1] << std::endl;}
			fich3.close();
			}
		int ng=ph_gcomp;
		for (int ig=0;ig<ng;ig++)
			{
			std::ofstream fich4(cur_dir+"/observation/obs_"+observ.headers[io]+"_Cg"+std::to_string(ig)+".txt", std::ios::app);// Ouvrir le fichier en mode append
			for (i=0;i<oindex;i++) {fich4 << obsGas[i*(nobs*ng+1)]/day<<" "<< obsGas[i*(nobs*ng+1)+io*ng+ig+1] << std::endl;}
			fich4.close();
			}
		}
	oindex=0;obsFlow.resize(0);obsTrans.resize(0);obsSolu.resize(0);obsGas.resize(0);
	}
