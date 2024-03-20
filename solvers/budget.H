// get the the water discharge/recharge at fixed heads cells and sum them for all cells to have the balance
int i1;
for (i=0; i<cellsHwel.size(); i++) {i1=cellsHwel[i];sqHwel+=qHwel[i1];}
for (i=0; i<cellsHdrn.size(); i++) {i1=cellsHdrn[i];qHdrn[i]=(hSu[i1]-h()[i1]*hSp[i1])*mesh.V()[i1];sqHdrn+=qHdrn[i];}
for (i=0; i<cellsHghb.size(); i++) {i1=cellsHghb[i];qHghb[i]=(hSu[i1]-h()[i1]*hSp[i1])*mesh.V()[i1];sqHghb+=qHghb[i];}// Info<<"qhghb "<<i<<" "<<i1<<" "<<qHghb[i]<<" q1 "<<q1[i1]*mesh.V()[i1]<<endl;}
for (i=0; i<cellsHriv.size(); i++) {i1=cellsHriv[i];qHriv[i]=(hSu[i1]-h()[i1]*hSp[i1])*mesh.V()[i1];sqHriv+=qHriv[i];} 
volScalarField q1 = fvc::div(phiw);
for (i=0;i<cellsHfix.size();i++)   {i1=cellsHfix[i];qHfix[i] = (q1[i1])*mesh.V()[i1];sqHfix+=qHfix[i];} //Info<<"qhfx "<<i<<" "<<i1<<" "<<qHfx[i]<<endl;}

/*if (activateReaction==1)
{
	for (int ic=0;ic<ph_ncomp;ic++)
	{
		for (i=0; i<cellsHwel.size(); i++) {i1=cellsHwel[i];sqSwel[ic]+=qHwel[i1]*Cw[ic]()[i1];}
		for (i=0; i<cellsHdrn.size(); i++) {i1=cellsHdrn[i];sqSdrn[ic]+=qHdrn[i1]*Cw[ic];}
		for (i=0; i<cellsHghb.size(); i++) {i1=cellsHghb[i];sqSghb[ic]+=qHghb[i1]*Cw()[ic];}// Info<<"qhghb "<<i<<" "<<i1<<" "<<qHghb[i]<<" q1 "<<q1[i1]*mesh.V()[i1]<<endl;}
		for (i=0; i<cellsHriv.size(); i++) {i1=cellsHriv[i];sqSriv[ic]+=qHriv[i1]*Cw()[ic];} 
//		for (i=0;i<cellsHfix.size();i++)   {i1=cellsHfix[i];sqSfix[ic]+=qHfix[i1]*solu_conc[condSfix[i1]*ph_ncomp+ic];} 
	}
}

if (flagW==1) // write the files    
	{
	float t = mesh.time().value()/day;
	std::ofstream fich1(cur_dir+"/observation/budgetH.txt", std::ios::app);// Ouvrir le fichier en mode append
	fich1 << t<<" "<<qHwel<<" "<<qHdrn<<" "<<qHghb<<" "<<qHriv<< std::endl;
	fich1.close();
	}*/
