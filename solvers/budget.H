// get the the water discharge/recharge at fixed heads cells and sum them for all cells to have the balance
int ib,i1;
std::cout<<"start budget \n";
for (i=0;i<cellsHfix.size();i++)
    {if (qHfix[i]>0) {sqIHfix+=qHfix[i]*dt;} else {sqOHfix+=qHfix[i]*dt;}} 
for (i=0;i<cellsHwel.size(); i++) 
	{if (qHwel[i]>0) {sqIHwel+=qHwel[i]*dt;} else {sqOHwel+=qHwel[i]*dt;}}
for (i=0;i<cellsHdrn.size(); i++)
	{if (qHdrn[i]>0) {sqIHdrn+=qHdrn[i]*dt;} else {sqOHdrn+=qHdrn[i]*dt;}}
for (i=0;i<cellsHghb.size(); i++)
	{if (qHghb[i]>0) {sqIHghb+=qHghb[i]*dt;} else {sqOHghb+=qHghb[i]*dt;}}// Info<<"qhghb "<<i<<" "<<i1<<" "<<qHghb[i]<<" q1 "<<q1[i1]*mesh.V()[i1]<<endl;}
for (i=0;i<cellsHrch.size(); i++) 
	{if (qHrch[i]>0) {sqIHrch+=qHrch[i]*dt;} else {sqOHrch+=qHrch[i]*dt;}}
float sq[11] = {mesh.time().value(),sqIHfix,sqOHfix,sqIHwel,sqOHwel,sqIHdrn,sqOHdrn,sqIHghb,sqOHghb,sqIHrch,sqOHrch};
budFlow.insert(budFlow.end(), sq, sq+11);
//std::cout<<"end budget flow \n";

float a;
//std::cout<<"act trp "<<activateTransport<<"\n";
if ((activateTransport==1)&&(activateReaction==0))
{
	for (i=0;i<cellsHfix.size();i++)
		{i1=cellsHfix[i];a=qHfix[i]*dt;if (qHfix[i]>0) {sqICfix+=a*cellsChfixData[i*4+2];} else {sqOCfix+=a*C[i1];} }
	for (i=0; i<cellsHwel.size(); i++)
		{i1=cellsHwel[i];a=qHwel[i]*dt;if (qHwel[i]>0) {sqICwel+=a*cellsCwelData[i*4+2];} else {sqOCwel+=a*C[i1];} }
	for (i=0; i<cellsHdrn.size(); i++)
		{i1=cellsHdrn[i];a=qHdrn[i]*dt;if (qHdrn[i]<0) {sqOCdrn+=a*C[i1];} } // drn only out
	for (i=0; i<cellsCghb.size(); i++)
		{i1=cellsHghb[i];a=qHghb[i]*dt;if (qHghb[i]>0) {sqICghb+=a*cellsCghbData[i*4+2];} else {sqOCghb+=a*C[i1];} } // in is bc conc, out is local conc
	for (i=0; i<cellsCrch.size(); i++)
	    {i1=cellsHrch[i];a=qHrch[i]*dt;if (qHrch[i]>0) {sqICrch+=a*cellsCrchData[i*4+2];}  } //else {sqOCrch+=a;} no flux or C from an rch BC
float sq[11] = {mesh.time().value(),sqICfix,sqOCfix,sqICwel,sqOCwel,sqICdrn,sqOCdrn,sqICghb,sqOCghb,sqICrch,sqOCrch};
budC.insert(budC.end(), sq, sq+11);
}
if (activateThermal==1)
{
	for (i=0;i<cellsHfix.size();i++)
		{i1=cellsHfix[i];a=qHfix[i]*dt;if (qHfix[i]>0) {sqITfix+=a*cellsThfixData[i*4+2];} else {sqOTfix+=a*T[i1];} }
	for (i=0; i<cellsHwel.size(); i++)
		{i1=cellsHwel[i];a=qHwel[i]*dt;if (qHwel[i]>0) {sqITwel+=a*cellsTwelData[i*4+2];} else {sqOTwel+=a*T[i1];} }
	for (i=0; i<cellsHdrn.size(); i++)
		{i1=cellsHdrn[i];a=qHdrn[i]*dt;if (qHdrn[i]<0) {sqOTdrn+=a*T[i1];} }
	for (i=0; i<cellsHghb.size(); i++)
		{i1=cellsHghb[i];a=qHghb[i]*dt;if (qHghb[i]>0) {sqITghb+=a*cellsTghbData[i*4+2];} else {sqOTghb+=a*T[i1];} }
	for (i=0; i<cellsHrch.size(); i++)
	    {i1=cellsHrch[i];a=qHrch[i]*dt;if (qHrch[i]>0) {sqITrch+=a*cellsTrchData[i*4+2];} } //else {sqOTrch+=a;} 
float sq[11] = {mesh.time().value(),sqICfix,sqOCfix,sqICwel,sqOCwel,sqICdrn,sqOCdrn,sqICghb,sqOCghb,sqICrch,sqOCrch};
budT.insert(budT.end(), sq, sq+11);
}
//std::cout<<"act reac "<<activateReaction<<"\n";
if (activateReaction==1)
{
    for (int ic=4;ic<ph_ncomp;ic++)
    {
        //std::cout<<"bud comp "<<ic<<"\n";
		for (i=0;i<cellsHfix.size();i++)
            {i1=cellsHfix[i];a=qHfix[i]*dt;if (qHfix[i]>0) {sqISfix[ic]+=a*cellsShfixData[i*4+2];} else {sqOSfix[ic]+=a*Cw[ic]()[i1];} }
        for (i=0; i<cellsHwel.size(); i++)
            {i1=cellsHwel[i];a=qHwel[i]*dt;if (qHwel[i]>0) {sqISwel[ic]+=a*cellsSwelData[i*4+2];} else {sqOSwel[ic]+=a*Cw[ic]()[i1];} }
        for (i=0; i<cellsHdrn.size(); i++)
            {i1=cellsHdrn[i];a=qHdrn[i]*dt;if (qHdrn[i]<0) {sqOSdrn[ic]+=a*Cw[ic]()[i1];} }
        for (i=0; i<cellsHghb.size(); i++)
            {i1=cellsHghb[i];a=qHghb[i]*dt;if (qHghb[i]>0) {sqISghb[ic]+=a*cellsSghbData[i*4+2];} else {sqOSghb[ic]+=a*Cw[ic]()[i1];} }
		for (i=0; i<cellsHrch.size(); i++)
            {i1=cellsHrch[i];a=qHrch[i]*dt;if (qHrch[i]>0) {sqISrch[ic]+=a;} }//else {sqOSrch[ic]+=a*Cw[ic]()[i1];} }
	float sq[11] = {mesh.time().value(),sqISfix[ic],sqOSfix[ic],sqISwel[ic],sqOSwel[ic],sqISdrn[ic],sqOSdrn[ic],sqISghb[ic],sqOSghb[ic],sqISrch[ic],sqOSrch[ic]};
	budSolu.insert(budSolu.end(), sq, sq+11);
    }
	size_t nc = cellsGfix.size();
    for (int ig=0;ig<ph_gcomp;ig++)
    {
		//diffusion
		for (i=0; i<nc; i++)
			// the gradC is obtained with 3 cells at a boundray and we suppose conc(calc as gradC not clear)
            {i1=cellsGfix[i];float de3,cg3,th3; // below valid only for 1D
			if (i1>3) {th3=thk()[i1]+thk()[i1-1]+0.5*thk()[i1-2]; cg3 = Cg[ig]()[i1-1];de3 = Deffg[i1-1].xx();}
			else {th3=thk[i1]+thk[i1+1]+0.5*thk[i1+2]; cg3 = Cg[ig]()[i1+1];de3 = Deffg[i1+1].xx();}
			sqOGdiff[ig] = de3*(cg3-0.0)/th3; // xx because it is a scalar * tensor I
			//Info<<"Flux diff "<<ig<<" "<<i1<<" "<<" "<<Deffg[i1-1]<<" "<<de3<<" "<<cg3<<" "<<th3<<" "<<sqOGdiff[ig] <<endl;
			budGas.push_back(sqOGdiff[ig]);
			}//else {sqOSrch[ic]+=a*Cw[ic]()[i1];} }
	}
}
bindex +=1;
//Info<<"bud "<<bindex<<" "<<budFlow.size()<<" "<<budC.size()<<" "<<budT.size()<<" "<<budSolu.size()<<" "<<budGas.size()<<endl;

if (flagW==1) // write the files    
	{
	float t = mesh.time().value() ;
	std::ofstream fich1(cur_dir+"/observation/budgetH.txt", std::ios::app);// Ouvrir le fichier en mode append
	//std::cout<<"bud flow indx "<<bindex<<" size budFlow "<<budFlow.size()<<"\n";
	//fich1 << t<<" "<<sqIHfix<<" "<<sqOHfix<<" "<<sqIHwel<<" "<<sqOHwel<<" "<<sqIHdrn<<" "<<sqOHdrn<<" "<<sqIHghb<<" "<<sqOHghb<<" "<<sqIHrch<<" "<<sqOHrch<< std::endl;
	for (int ib=0;ib<bindex;ib++)
		{
		for (i=0;i<11;i++) {fich1<<budFlow[ib*11+i]<<" ";}
		fich1<<"\n";
		}
	fich1.close();
	sqIHfix=0;sqIHwel=0;sqIHdrn=0;sqIHghb=0;sqIHrch=0;
	sqOHfix=0;sqOHwel=0;sqOHdrn=0;sqOHghb=0;sqOHrch=0;
	//transport
	if ((activateTransport==1)&&(activateReaction==0))
		{
		std::ofstream fich1(cur_dir+"/observation/budgetC.txt", std::ios::app);// Ouvrir le fichier en mode append
		//std::cout<<"bud C indx "<<bindex<<" size budC "<<budC.size()<<"\n";
		for (int ib=0;ib<bindex;ib++)
			{
			for (i=0;i<11;i++) {fich1<<budC[ib*11+i]<<" ";}
			fich1<<"\n";
			}
		fich1.close();
		sqICfix=0;sqICwel=0;sqICdrn=0;sqICghb=0;sqICrch=0;
		sqOCfix=0;sqOCwel=0;sqOCdrn=0;sqOCghb=0;sqOCrch=0;
		}
	//thermal
	if (activateThermal==1)
	{
		std::ofstream fich1(cur_dir+"/observation/budgetT.txt", std::ios::app);// Ouvrir le fichier en mode append
		for (int ib=0;ib<bindex;ib++)
			{
			for (i=0;i<11;i++) {fich1<<budT[ib*11+i]<<" ";}
			fich1<<"\n";
			}
		fich1.close();
		sqITfix=0;sqITwel=0;sqITdrn=0;sqITghb=0;sqITrch=0;
		sqOTfix=0;sqOTwel=0;sqOTdrn=0;sqOTghb=0;sqOTrch=0;
	}
	if (activateReaction==1)
	{
		int nc = ph_ncomp -4;
		for (int ic=4;ic<ph_ncomp;ic++)
			{
			std::string s=std::to_string(ic);
			std::ofstream fich1(cur_dir+"/observation/budgetCw"+s+".txt", std::ios::app);// Ouvrir le fichier en mode append
			for (int ib=0;ib<bindex;ib++)
				{
				for (i=0;i<11;i++) {fich1<<budSolu[ib*nc*11+ic*11+i]<<" ";}
				fich1<<"\n";
				}
			fich1.close();
			sqISfix[ic]=0;sqISwel[ic]=0;sqISdrn[ic]=0;sqISghb[ic]=0,sqISrch[ic]=0;
			sqOSfix[ic]=0;sqOSwel[ic]=0;sqOSdrn[ic]=0;sqOSghb[ic]=0,sqOSrch[ic]=0;
			}
		int ng = ph_gcomp;//std::cout<<"bud gcomp "<<ng<<" size budGas "<<budGas.size()<<"\n";
		for (int ig=0;ig<ph_gcomp;ig++)
			{
			std::string s=std::to_string(ig);
			std::ofstream fich1(cur_dir+"/observation/budgetCg"+s+".txt", std::ios::app);// Ouvrir le fichier en mode append
			//here jus tone cell
			for (int ib=0;ib<bindex;ib++)
				{
				fich1<<budFlow[ib*11]<<" "<<budGas[ib*ng+ig]<<"\n";  // budFlow is time
				}
			fich1.close();
			sqOGdiff[ig]=0;
			}
	}
	bindex=0;budFlow.resize(0);budC.resize(0);budT.resize(0);budSolu.resize(0);budGas.resize(0);
}
