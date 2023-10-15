char sep ='/';
int nxyz0=2;ph_ncomp=8;ph_gcomp=0;ph_nsolu=2;
std::vector<int> data={2,8,0,2,1, -1,0,1, 0,-1,0,-1,0,-1,0,-1,0,-1,0,1};

freak.setDB(cur_dir+sep+"phreeqc.dat");
//make a first init to calculate the inital solutions (poro required for the initial solutions)
freak.setChemFile(cur_dir+sep+"initChem.pqi"); //Info << "initCh read " << endl;
freak.setData(data);
std::vector<double> poro(ph_nsolu,0.25);
std::vector<double> wsat(ph_nsolu,0.9999);
freak.setPoro(poro);
freak.setWsat(wsat);
freak.init(); 

//first step make 50 time steps for solution 1
int icell,it,k;
int nt=50;std::vector<double> c(nt*nvar,0.);
for (it=1;it<nt;it++) { //50 time steps
	freak.setTstep(it*86400.);
	freak.run();icell=1;
	for (i=0;i<nvar;i++) { 
		c[it*nvar+i]=freak.c[(i+4)*ph_nsolu+icell];
		Cmax[i] = std::max(Cmax[i],c[it*nvar+i]);
		}
}
//std::vector<double> sbase0(4);
//for (i=0;i<4;i++) {sbase0[i] = freak.c[i*nxyz0+0];}
// then from these 50 solutions, mix them with solution zero by 1/10
std::vector<double> c0(nvar,0.);int nk=26;float stpk=1./(nk-1);
for (i=0;i<nvar;i++) {c0[i]=freak.c[(i+4)*nxyz0+0];} // intial solution
std::vector<double> c1(nt*nvar*nk,0.);
for (it=1;it<nt;it++) {
	for (k=0;k<nk;k++) {
		for (i=0;i<nvar;i++) {
			c1[it*nk*nvar+k*nvar+i]=k*stpk*c0[i]+(nk-k-1)*stpk*c[it*nvar+i];
			}
	}
}

// use these concentrations to make a run for a given time step and get dC/dt
// make conc relative, here we keep c0 as cmin, and store in nndata
nxyz0 = nt*nk;
data={nxyz0,8,0,nxyz0,1, 0,0, 0,-1,0,-1,0,-1,0,-1,0,-1, 0,1};
freak.setData(data);
poro.resize(nxyz0,0.25);
wsat.resize(nxyz0,0.9999);
freak.setPoro(poro);
freak.setWsat(wsat);
freak.init(); //std::cout << "freak H2O_0 " << freak.c[0] << "\n";
double dt0 = 86400.;
freak.setTstep(dt0);
std::vector<double> c_ph(nxyz0*ph_ncomp,0);
for (int j=0;j<nxyz0;j++) {
	for (i=0;i<4;i++) {c_ph[i*nxyz0+j] = freak.c[i*nxyz0+j];}
	for (i=0;i<nvar;i++) {c_ph[(i+4)*nxyz0+j] = c1[j*nvar+i];}
	//if (i==0) {std::cout<<c_ph[0]<<"\n";}
	}
freak.setC(c_ph);
freak.run();//std::cout<<"end run \n";
//shuffle
std::random_device rd;
std::mt19937 shuf(rd());
std::vector<int> idx(nt*nk);
std::iota(idx.begin(), idx.end(), 0);		
std::shuffle(idx.begin(), idx.end(), shuf);
//std::ifstream inputIdx{cur_dir+sep+"idx"};
//std::vector<int> idx{std::istream_iterator<int>{inputIdx}, {}}; 
float grr;
for (int j=1;j<nxyz0;j++) {
	int j1=idx[j];
	for (i=0;i<nvar;i++) {
		grr=(c1[j1*nvar+i]-c0[i])/(Cmax[i]-c0[i]);nndata.push_back(grr);outNNdata<<grr<<" "; //crelative
		grr=(freak.c[(i+4)*nxyz0+j1]-c1[j1*nvar+i])/dt0*5e9;nntarget.push_back(grr);outNNtarget<<grr<<" "; 
	}
	outNNdata<<"\n"; outNNtarget<<"\n";
}
std::flush(outNNdata);std::flush(outNNtarget);outNNdata.close();outNNtarget.close();
// train
Cwgnn.setNd(nxyz0);Cwgnn.setData(nndata);Cwgnn.setRunParms({400.,32.,1e-3,.75});//(epoc,batch,learning rate,%train)
Cwgnn.setTarget(nntarget);
float rmse = Cwgnn.train(); std::cout << "rmse "<<rmse<<std::endl;
outNNrmse <<" " << rmse << "\n";

