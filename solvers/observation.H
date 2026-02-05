obsFlow.push_back(mesh.time().value());
Info<<"nobs "<<nobs<<" len obsIndx "<<obsFlowIndx.size()<<" obTrans "<<obsTransIndx.size()<<" chem "<<obsChemIndx.size()<<endl;
if ((activateTransport+activateThermal>0)&&(obsTransIndx.size()>0)) {obsTrans.push_back(mesh.time().value());}
if ((activateReaction==1)&&(obsChemIndx.size()>0))
	{
	if (obsChemIndx[0]==1) {obsSolu.push_back(mesh.time().value());}
	if ((ph_gcomp>0)&&(obsChemIndx[1]==1)) {obsGas.push_back(mesh.time().value());}
	}

for (int io=0;io<nobs;io++) 
	{
	Info<<"in obs "<<io<<" flow t "<<flowType;
	if (obsFlowIndx.size()>0)
		{
		if (obsFlowIndx[0]==1) 
			{
			if (flowType<=2) {obsFlow.push_back(h()[icello[io]]);} //
			if (flowType==3) {obsFlow.push_back(hp()[icello[io]]);}
			}
		if (obsFlowIndx[1]==1) {obsFlow.push_back(sw()[icello[io]]);}
		if (obsFlowIndx[2]==1) {obsFlow.push_back(mag(Uw()[icello[io]]));}
		}
	if ((activateTransport+activateThermal>0)&&(obsTransIndx.size()>0))
		{
		if (obsTransIndx[0]==1) {obsTrans.push_back(C()[icello[io]]);} //Info<<"in obs store t "<<mesh.time().value()<<" C "<<C()[icello[io]]<<endl;
		if (obsTransIndx[1]==1) {obsTrans.push_back(T()[icello[io]]);} //Info<<"in obs store t "<<mesh.time().value()<<" T "<<T()[icello[io]]<<endl;
		}

	if ((activateReaction==1)&&(obsChemIndx.size()>0))
		{
		if (obsChemIndx[0]==1)
			{
			for (int ic=0;ic<ph_ncomp;ic++) {obsSolu.push_back(Cw[ic]()[icello[io]]);}
			}
		if ((ph_gcomp>0)&&(obsChemIndx[1]==1))
			{
			for (int ig=0;ig<ph_gcomp;ig++) {obsGas.push_back(Cg[ig]()[icello[io]]);}
			}
		}
	}
oindex += 1;std::cout<<"oindex "<<oindex<<std::endl;

std::vector<std::string> obsFlowNames = {"h","sw","Uw"};
std::vector<std::string> obsTransNames = {"C","T"};
int iv,iv1,ic,ig,nv;
if (flagW==1) // write the files    
	{
	Info<<"obs size obsFlow "<<obsFlow.size()<<" "<<obsFlow[0]<<" "<<obsFlow[1]<<" "<<obsFlow[2]<<" "<<obsFlow[3]<<endl;
	for (int io=0;io<nobs;io++) //one file for each ob point
		{
		//fro flow and transp only some variables can be stored
		if (obsFlowIndx.size()>0)
			{
			nv=obsFlowIndx[0]+obsFlowIndx[1]+obsFlowIndx[2];iv1=0;
			//Info<<"nobs flow "<<nobs<<" nv "<<nv<<endl;
			for (iv=0;iv<3;iv++)
				{
				//Info <<"in owrite "<<io<<" "<<iv<<" "<<obsFlowIndx[iv]<<" "<<obsFlowNames[obsFlowIndx[iv]]<<endl;				
				if (obsFlowIndx[iv]==1) // that variable must be printed
					{
					iv1+=1;
					std::ofstream fich1(cur_dir+"/observation/obs_"+observ.headers[io]+"_"+obsFlowNames[iv]+".txt", std::ios::app);// Ouvrir le fichier en mode append
					for (i=0;i<oindex;i++) {fich1 << obsFlow[i*(nobs*nv+1)] <<" "<< obsFlow[i*(nobs*nv+1)+io*nv+iv1] << "\n";} //
					fich1.close();
					}
				}
			}
		if ((activateTransport+activateThermal>0) && (obsTransIndx.size()>0)) //
			{
			nv=obsTransIndx[0]+obsTransIndx[1];iv1=0; // 
			for (iv=0;iv<2;iv++)
				if (obsTransIndx[iv]==1)
					{
					iv1 +=1;
					std::ofstream fich2(cur_dir+"/observation/obs_"+observ.headers[io]+"_"+obsTransNames[iv]+".txt", std::ios::app);// Ouvrir le fichier en mode append
					Info<<"nobs trans "<<nobs<<" nv "<<nv<<" iv "<<iv<<" iv1 "<<iv1<<endl;
					for (i=0;i<oindex;i++) {fich2 << obsTrans[i*(nobs*nv+1)] <<" "<< obsTrans[i*(nobs*nv+1)+io*nv+iv1] << "\n";}
					fich2.close();
					}
			}
		if ((activateReaction==1)&&(obsChemIndx.size()>0))
			{
			int nc=ph_ncomp;
			for (ic=0;ic<nc;ic++)
				{
				std::ofstream fich3(cur_dir+"/observation/obs_"+observ.headers[io]+"_Cw"+std::to_string(ic)+".txt", std::ios::app);// Ouvrir le fichier en mode append
				for (i=0;i<oindex;i++) {fich3 << obsSolu[i*(nobs*nc+1)] <<" "<< obsSolu[i*(nobs*nc+1)+io*nc+ic+1] << std::endl;}
				fich3.close();
				}
			int ng=ph_gcomp;
			for (ig=0;ig<ng;ig++)
				{
				std::ofstream fich4(cur_dir+"/observation/obs_"+observ.headers[io]+"_Cg"+std::to_string(ig)+".txt", std::ios::app);// Ouvrir le fichier en mode append
				for (i=0;i<oindex;i++) {fich4 << obsGas[i*(nobs*ng+1)] <<" "<< obsGas[i*(nobs*ng+1)+io*ng+ig+1] << std::endl;}
				fich4.close();
				}
			}
		}
	oindex=0;obsFlow.resize(0);obsTrans.resize(0);obsSolu.resize(0);obsGas.resize(0);
	}
//used for NN unsaturated
/*
std::ofstream fNNhp(cur_dir+"/observation/nnHp.txt", std::ios::app);// Ouvrir le fichier en mode append
fNNhp << mesh.time().value()<<" "<<qHwel[0]<<" ";
for (i=0;i<ncell;i++) {fNNhp << hp[i] << " " ;}
fNNhp << std::endl;
*/