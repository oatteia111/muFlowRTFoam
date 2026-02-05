//update 2 phase but with foam
sw = max(sw_min,min(sw,sw_max));
Se = min(max((sw-sw_min)/(sw_max-sw_min),1e-9),1-1e-9);//Info<<"Se "<<Se[0]<<" "<<Se[1]<<endl;
krw = pow(Se,0.5) * pow(1-pow(1-pow(Se,1/m_vg),m_vg),2);//Info<<"krw "<<krw[0]<<endl; //krmax = 1
krg = pow(1-Se,0.5) * pow(1-pow(Se,1/m_vg),2*m_vg);//Info<<"krg "<<krg<<endl; // OA 12/9 2*m -> 5*m
krg = max(min(krg,.99999),1e-5); //Info<< "sw "<<sw[0]<<" krg new "<< krg[0] <<" "<< krg[1] << endl;
//foam
volScalarField Fdry=0.5+atan(plugPS.epdry*(sw-plugPS.fmdry))/pi;
// 
volScalarField Fc=unity;
//use surfactant cocnentration if reaction or transpor is simulated (for flow only, we consider medium saturated with Surfact)
int isurf = plugPS.isurfact; //Info<<" in plug isurf "<<isurf<<endl;
if ((ph_ncomp>0)&&(isurf<ph_ncomp)) {for (j=1;j<ncell;j++) {Fc = 0.5+std::atan(plugPS.fc*(Cw[isurf]()[j]-plugPS.cref))/pi;}}
else if (activateTransport==1) {for (j=1;j<ncell;j++) {Fc = 0.5+std::atan(plugPS.fc*(C()[j]-plugPS.cref))/pi;}} //Cw[isurfact]
Info<<"plug foam ep, fm "<<plugPS.epdry<<" "<<plugPS.fmdry<<" Fdry[0] "<<Fdry[0]<<" Fc[0] "<<Fc[0]<<endl;
//if (ph_ncomp>0) {Info<<"isurf "<<isurf<<" Cw isurf "<<Cw[isurf]()[0]<<" Fc 0,1 "<<Fc[0]<<" "<<Fc[1]<<endl;}
krg /= 1+plugPS.fmmob*Fdry*Fc;Info<< "sw "<<sw[0]<<" krg new "<< krg[0] <<" "<< krg[1] <<"  fmmob "<<plugPS.fmmob<<" Fc0 "<<Fc[0]<<" Fdry0 "<<Fdry[0]<< endl;
krwf = fvc::interpolate(krw,"krw");
krgf = fvc::interpolate(krg,"krg");

	pc = pc0*pow(pow(Se,-1/m_vg)-1,1/n_vg); //Se=(1+(p/pc0)^n)^-m; Se^-1/m-1=(p/pc0)^n
	//Info<<" pc "<<pc;
	dpcdS = -pc0/(n_vg*m_vg)*pow(Se,-1/m_vg - 1); //a b x^(a - 1);a=-1/m,b=1/n
	dpcdS *= pow(pow(Se,-1/m_vg)-1,1/n_vg-1) ;//(x^a - 1)^(b - 1)
	//Info<<" dpcdS "<<dpcdS;

//- mobility computation 
Mwf = Kf *krwf/muwf;
Mf = Kf *(krwf/muwf + krgf/mug);
Info<<"rhog "<<rhog<<" mug "<<mug<<"krw 6 "<<krwf[6]<<" Mf(6) "<<Mf[6]<<endl;
Lwf = Kf*rhowf*krwf/muwf;
Lf = Lwf + Kf*rhog*krgf/mug;
Info<<"Lwf 6 "<<Lwf[6]<<"Lf 6 "<<Lf[6]<<endl;
Fwf = (krwf/muwf) / (krgf/mug + krwf/muwf);
//Info<<"Mf done"<<endl;
//- compute fluxes
phiGr = (Lf * g) & mesh.Sf();//Info<<"phiGr "<<"PhiGr(10) "<<phiGr[10]<<endl;
gradpcf = fvc::interpolate(dpcdS*fvc::grad(sw),"pc");//Info<<"dpcds "<<dpcdS[2]<<" sw "<<sw[1]<<" "<<sw[2]<<endl;
phiPc = (Mwf * gradpcf) & mesh.Sf();//Info<<"PhiPc(10) "<<phiPc[10]<<endl;
