// get the the water discharge/recharge at fixed heads cells and sum them for all cells to have the balance
int i1;
for (i=0; i<cellsHwel.size(); i++) { if (qHwel[i]>0) {sqIHwel+=qHwel[i];} else {sqOHwel+=qHwel[i];} }
for (i=0; i<cellsHdrn.size(); i++)
 {i1=cellsHdrn[i];qHdrn[i]=(hSu[i1]-h()[i1]*hSp[i1])*mesh.V()[i1];if (qHdrn[i]>0) {sqIHdrn+=qHdrn[i];} else {sqOHdrn+=qHdrn[i];}}
for (i=0; i<cellsHghb.size(); i++)
 {i1=cellsHghb[i];qHghb[i]=(hSu[i1]-h()[i1]*hSp[i1])*mesh.V()[i1];if (qHghb[i]>0) {sqIHghb+=qHghb[i];} else {sqOHghb+=qHghb[i];}}// Info<<"qhghb "<<i<<" "<<i1<<" "<<qHghb[i]<<" q1 "<<q1[i1]*mesh.V()[i1]<<endl;}
for (i=0; i<cellsHriv.size(); i++)
 {i1=cellsHriv[i];qHriv[i]=(hSu[i1]-h()[i1]*hSp[i1])*mesh.V()[i1];if (qHriv[i]>0) {sqIHriv+=qHriv[i];} else {sqOHriv+=qHriv[i];}}
volScalarField q1 = fvc::div(phiw);
for (i=0;i<cellsHfix.size();i++)
    {i1=cellsHfix[i];qHfix[i] = q1[i1]*mesh.V()[i1];if (qHfix[i]>0) {sqIHfix+=qHfix[i];} else {sqOHfix+=qHfix[i];}} //Info<<"qhfx "<<i<<" "<<i1<<" "<<qHfx[i]<<endl;}

float a;
if (activateTransport==1)
{
	for (i=0; i<cellsHwel.size(); i++)
		{i1=cellsHwel[i];a=qHwel[i]*C[i1];if (qHwel[i]>0) {sqICwel+=a;} else {sqOCwel+=a;} }
	for (i=0; i<cellsHdrn.size(); i++)
		{i1=cellsHdrn[i];a=qHdrn[i]*C[i1];if (qHdrn[i]>0) {sqICdrn+=a;} else {sqOCdrn+=a;} }
	for (i=0; i<cellsHghb.size(); i++)
		{i1=cellsHghb[i];a=qHghb[i]*C[i1];if (qHghb[i]>0) {sqICghb+=a;} else {sqOCghb+=a;} }
	//for (i=0; i<cellsHriv.size(); i++)
	//    {i1=cellsHriv[i];a=qHriv[i]*Cw()[i1];if (qHriv[i]>0) {sqICriv+=a;} else {sqOCriv+=a;} }
	for (i=0;i<cellsHfix.size();i++)
		{i1=cellsHfix[i];a=qHfix[i]*C[i1];if (qHfix[i]>0) {sqICfix+=a;} else {sqOCfix+=a;} }
}
if (activateThermal==1)
{
	for (i=0; i<cellsHwel.size(); i++)
		{i1=cellsHwel[i];a=qHwel[i]*T[i1];if (qHwel[i]>0) {sqITwel+=a;} else {sqOTwel+=a;} }
	for (i=0; i<cellsHdrn.size(); i++)
		{i1=cellsHdrn[i];a=qHdrn[i]*T[i1];if (qHdrn[i]>0) {sqITdrn+=a;} else {sqOTdrn+=a;} }
	for (i=0; i<cellsHghb.size(); i++)
		{i1=cellsHghb[i];a=qHghb[i]*T[i1];if (qHghb[i]>0) {sqITghb+=a;} else {sqOTghb+=a;} }
	//for (i=0; i<cellsHriv.size(); i++)
	//    {i1=cellsHriv[i];a=qHriv[i]*Tw()[i1];if (qHriv[i]>0) {sqITriv+=a;} else {sqOTriv+=a;} }
	for (i=0;i<cellsHfix.size();i++)
		{i1=cellsHfix[i];a=qHfix[i]*T[i1];if (qHfix[i]>0) {sqITfix+=a;} else {sqOTfix+=a;} }
}
if (activateReaction==1)
{
    for (int ic=0;ic<ph_ncomp;ic++)
    {
        for (i=0; i<cellsHwel.size(); i++)
            {i1=cellsHwel[i];a=qHwel[i]*Cw[ic]()[i1];if (qHwel[i]>0) {sqISwel[ic]+=a;} else {sqOSwel[ic]+=a;} }
        for (i=0; i<cellsHdrn.size(); i++)
            {i1=cellsHdrn[i];a=qHdrn[i]*Cw[ic]()[i1];if (qHdrn[i]>0) {sqISdrn[ic]+=a;} else {sqOSdrn[ic]+=a;} }
        for (i=0; i<cellsHghb.size(); i++)
            {i1=cellsHghb[i];a=qHghb[i]*Cw[ic]()[i1];if (qHghb[i]>0) {sqISghb[ic]+=a;} else {sqOSghb[ic]+=a;} }
		//for (i=0; i<cellsHriv.size(); i++)
        //    {i1=cellsHriv[i];a=qHriv[i]*Cw[ic]()[i1];if (qHriv[i]>0) {sqISriv[ic]+=a;} else {sqOSriv[ic]+=a;} }
        for (i=0;i<cellsHfix.size();i++)
            {i1=cellsHfix[i];a=qHfix[i]*Cw[ic]()[i1];if (qHfix[i]>0) {sqISfix[ic]+=a;} else {sqOSfix[ic]+=a;} }
    }
}


if (flagW==1) // write the files    
	{
	float t = mesh.time().value() ;
	std::ofstream fich1(cur_dir+"/observation/budgetH.txt", std::ios::app);// Ouvrir le fichier en mode append
	fich1 << t<<" "<<sqIHfix<<" "<<sqOHfix<<" "<<sqIHwel<<" "<<sqOHwel<<" "<<sqIHdrn<<" "<<sqOHdrn<<" "<<sqIHghb<<" "<<sqOHghb<< std::endl;
	fich1.close();
	sqIHfix=0;sqIHwel=0;sqIHdrn=0;sqIHghb=0;//sqIHrch=0;
	sqOHfix=0;sqOHwel=0;sqOHdrn=0;sqOHghb=0;//sqOHrch=0;
	//transport
	if (activateTransport==1)
		{
		std::ofstream fich1(cur_dir+"/observation/budgetC.txt", std::ios::app);// Ouvrir le fichier en mode append
		fich1 << t<<" "<<sqICfix<<" "<<sqOCfix<<" "<<sqICwel<<" "<<sqOCwel<<" "<<sqICdrn<<" "<<sqOCdrn<<" "<<sqICghb<<" "<<sqOCghb<< std::endl;
		fich1.close();
		sqICfix=0;sqICwel=0;sqICdrn=0;sqICghb=0;//sqICrch=0;
		sqOCfix=0;sqOCwel=0;sqOCdrn=0;sqOCghb=0;//sqOCrch=0;
		}
	//thermal
	if (activateThermal==1)
	{
		std::ofstream fich1(cur_dir+"/observation/budgetT.txt", std::ios::app);// Ouvrir le fichier en mode append
		fich1 << t<<" "<<sqITfix<<" "<<sqOTfix<<" "<<sqITwel<<" "<<sqOTwel<<" "<<sqITdrn<<" "<<sqOTdrn<<" "<<sqITghb<<" "<<sqOTghb<< std::endl;
		fich1.close();
		sqITfix=0;sqITwel=0;sqITdrn=0;sqITghb=0;//sqITrch=0;
		sqOTfix=0;sqOTwel=0;sqOTdrn=0;sqOTghb=0;//sqOTrch=0;
	}
	if (activateReaction==1)
	{
		for (int ic=0;ic<ph_ncomp;ic++)
			{
			std::string s=std::to_string(ic);
			std::ofstream fich1(cur_dir+"/observation/budgetCw"+s+".txt", std::ios::app);// Ouvrir le fichier en mode append
			fich1 << t<<" "<<sqISfix[ic]<<" "<<sqOSfix[ic]<<" "<<sqISwel[ic]<<" "<<sqOSwel[ic];
			fich1 <<" "<<sqISdrn[ic]<<" "<<sqOSdrn[ic]<<" "<<sqISghb[ic]<<" "<<sqOSghb[ic]<<"\n";//" "<<sqISrch[ic]<<" "<<sqOSrch[ic]<<"\n";
			fich1.close();
			sqISfix[ic]=0;sqISwel[ic]=0;sqISdrn[ic]=0;sqISghb[ic]=0,sqISrch[ic]=0;
			sqOSfix[ic]=0;sqOSwel[ic]=0;sqOSdrn[ic]=0;sqOSghb[ic]=0,sqOSrch[ic]=0;
			}
		}
	}
