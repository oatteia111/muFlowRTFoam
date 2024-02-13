obs.push_back(mesh.time().value());
for (int io=0;io<nobs;io++) 
	{
	obs.push_back(hp()[icello[io]]);
	obs.push_back(sw()[icello[io]]);
	}
oindex += 1;std::cout<<"oindex "<<oindex<<std::endl;

if (oindex>100) // write the files    
	{
	for (int io=0;io<nobs;io++) //one file for each ob point
		{
		std::ofstream fichier(cur_dir+"/observation/obs"+std::to_string(io)+".txt", std::ios::app);// Ouvrir le fichier en mode append
		if (fichier.is_open()) {
			for (i=0;i<oindex;i++) 
				{fichier << obs[i*(nobs*2+1)]/day<<" "<< obs[i*(nobs*2+1)+io*2+1]<<" "<< obs[i*(nobs*2+1)+io*2+2] << std::endl;} //+1 for time
			fichier.close();
			} 
		else {
			std::cerr << "Impossible d'ouvrir le fichier." << std::endl;
			}
		}
	oindex=0;obs.resize(0);
	}
