obs.push_back(mesh.time().value());
for (int io=0;io<nobs;io++) 
	{
	obs.push_back(hp()[icello[io]]);
	obs.push_back(sw()[icello[io]]);
	}
oindex += 1;std::cout<<"oindex "<<oindex<<std::endl;

if (flagW==1) // write the files    
	{
	for (int io=0;io<nobs;io++) //one file for each ob point
		{
		std::ofstream fich1(cur_dir+"/observation/obs_"+observ.headers[io]+"_hp.txt", std::ios::app);// Ouvrir le fichier en mode append
		for (i=0;i<oindex;i++) 
				{fich1 << obs[i*(nobs*2+1)]/day<<" "<< obs[i*(nobs*2+1)+io*2+1] << std::endl;} //+1 for time
		fich1.close();
		std::ofstream fich2(cur_dir+"/observation/obs_"+observ.headers[io]+"_sw.txt", std::ios::app);// Ouvrir le fichier en mode append
		for (i=0;i<oindex;i++) 
				{fich2 << obs[i*(nobs*2+1)]/day<<" "<< obs[i*(nobs*2+1)+io*2+2]<< std::endl;} //+1 for time
		fich2.close();
		}
	oindex=0;obs.resize(0);
	}
