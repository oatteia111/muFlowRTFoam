
sw = min(sw,sw_max*.9999);
//- relative permeability computation
	Se = min(max((sw-sw_min)/(sw_max-sw_min),1e-5),0.99999);//Info<<"Se "<<Se[0]<<" "<<Se[1]<<endl;
	//Info<<"se done"<<endl;
	// water
	krw = 1 * pow(Se,0.5) * pow(1-pow(1-pow(Se,1/m_vg),m_vg),2);Info<<"krw "<<krw[0]<<endl; //krmax = 1
	dkrwdS = 0.5 * (1-pow((1-pow(Se,1/m_vg)),m_vg));
	dkrwdS *= ( 4 * pow(Se,1/m_vg-1/2) * pow( (1-pow(Se,1/m_vg)) , m_vg-1)) - ( pow((1-pow(Se,1/m_vg)),m_vg) -1) / pow(Se,0.5);
	dkrwdS *= 1/(sw_max - sw_min);
	//Info<<"kr done"<<endl;
	//gas;
	krg = pow(1-Se,0.5) * pow(1-pow(Se,1/m_vg),2*m_vg);//Info<<"krg "<<krg<<endl; // OA 12/9 2*m -> 5*m
	//krg = 1.3978 - 3.7694*sw + 12.7090*pow(sw,2) - 20.6420*pow(sw,3) + 10.3090*pow(sw,4); //tough formulation
	//krg = 1-(sw-0.2)*0.5-0.623*exp(-10.5*pow(sw-0.94,2.));
	krg = max(min(krg,.99999),1e-5); //Info<< "sw "<<sw[0]<<" krg new "<< krg[0] <<" "<< krg[1] << endl;
	dkrgdS = - pow((1-pow(Se,1/m_vg)),2*m_vg-1) ;
	dkrgdS *= (-5*pow(Se,1/m_vg+1)+4*pow(Se,1/m_vg)+Se);
	dkrgdS *= 1/(2*pow((1-Se),0.5)*Se);
	dkrgdS *=  1/(sw_max - sw_min);

// capillary pressure
	//VG
	dpcdS = -1/(n_vg*m_vg)*pc0/(sw_max-sw_min); // swx-swn because of Se vs sw
	dpcdS *= pow(pow(Se,-1/m_vg)-1,(1/n_vg)-1); 
	dpcdS *= pow(Se,-(1+m_vg)/m_vg);
	//brooks corey
	//dpcdS = -alpha_vg*pc0*(pow(Se,-alpha_vg-1))/(sw_max-sw_min);
	// same as above to validate
	//dpcdS = - (1/(n_vg*m_vg)) * (pc0/(sw_max-sw_min)); 
	//dpcdS *= pow(pow(Se,-1/m_vg)-1,1/n_vg-1); 
	//dpcdS *= pow(Se,-(1+m_vg)/m_vg);//Info<<"dpc "<<dpcdS<<endl;
//Info<<"dpcpdS done"<<endl;
krwf = fvc::interpolate(krw,"krw");
krgf = fvc::interpolate(krg,"krg");

//- mobility computation 
Mwf = Kf*krwf/muw;
Mf = Mwf + Kf*krgf/mug;
//Info<<"Mf "<<Mf[10]<<endl;
Lwf = Kf*rhow*krwf/muw;
Lf = Lwf + Kf*rhog*krgf/mug;
//Info<<"Lf 10 "<<Lf[10]<<" size "<<Lf.size()<<endl;
Fwf = (krwf/muw) / (krgf/mug + krwf/muw);
//Info<<"Mf done"<<endl;
//- compute fluxes
gradpc = fvc::interpolate(dpcdS*fvc::grad(sw),"pc");//Info<<"dpcds "<<dpcdS[2]<<" sw "<<sw[1]<<" "<<sw[2]<<endl;
phiGr = (Lf * g) & mesh.Sf();//Info<<"phiGr ok "<<"PhiGr "<<phiGr[10]<<endl;
phiPc = (Mwf * gradpc) & mesh.Sf();//Info<<"<<endl;
	//Info<<"PhiG "<<phiG[10]<<endl;
