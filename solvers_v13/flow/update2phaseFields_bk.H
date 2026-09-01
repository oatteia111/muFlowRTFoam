
// derivatives
	//VG
	//Info<<" pc "<<pc;
	//dkrwdS = 0.5 * (1-pow((1-pow(Se,1/m_vg)),m_vg));
	//dkrwdS *= ( 4 * pow(Se,1/m_vg-1/2) * pow( (1-pow(Se,1/m_vg)) , m_vg-1)) - ( pow((1-pow(Se,1/m_vg)),m_vg) -1) / pow(Se,0.5);
	//dkrwdS *= 1/(sw_max - sw_min);
	//dkrgdS = - pow((1-pow(Se,1/m_vg)),2*m_vg-1) ;
	//dkrgdS *= (-5*pow(Se,1/m_vg+1)+4*pow(Se,1/m_vg)+Se);
	//dkrgdS *= 1/(2*pow((1-Se),0.5)*Se);
	//dkrgdS *=  1/(sw_max - sw_min);
	//
	//dpcdS = -pc0/(n_vg*m_vg)*pow(Se,-1/m_vg - 1); //a b x^(a - 1);a=-1/m,b=1/n
	//dpcdS *= pow(pow(Se,-1/m_vg)-1,1/n_vg-1) ;//(x^a - 1)^(b - 1)
	//Info<<" dpcdS "<<dpcdS;
	//dpcdS = -1/(n_vg*m_vg)*pc0/(sw_max-sw_min)*1; // swx-swn because of Se vs sw
	//dpcdS *= pow(pow(Se,-1/m_vg)-1,(1/n_vg)-1); 
	//dpcdS *= pow(Se,-(1+m_vg)/m_vg);
	//brooks corey
	//dpcdS = -alpha_vg*pc0*(pow(Se,-alpha_vg-1))/(sw_max-sw_min);
	// same as above to validate
	//dpcdS = - (1/(n_vg*m_vg)) * (pc0/(sw_max-sw_min)); 
	//dpcdS *= pow(pow(Se,-1/m_vg)-1,1/n_vg-1); 
	//dpcdS *= pow(Se,-(1+m_vg)/m_vg);//Info<<"dpc "<<dpcdS<<endl;
	//Info<<"in upd sw "<<sw[0]<<" "<<sw_min[0]<<" "<<sw_max[0]<<" "<<sw_max[0]-sw[0]<<endl;
	sw = max(sw_min,min(sw,sw_max));//Info<<"in upd sw "<<sw[0]<<" "<<sw_min[0]<<" "<<sw_max[0]<<" "<<sw_max[0]-sw[0]<<endl;
	double mysmall = 1e-18;
	Se = min(max((sw-sw_min)/(sw_max-sw_min),mysmall),1-mysmall);//Info<<"Se0 "<<Se[0]<<" "<<1-Se[0]<<endl;
	//if (activateEbullition) {Seg = min(max((sw-sw_min)/(1-sw_min),mysmall),1-mysmall);} else {Seg = Se;}//Segas can reach a saturation of 1	

	//- mobility computation 
	krw = 1 * pow(Se,0.5) * pow(1-pow(1-pow(Se,1/m_vg),m_vg),2);
	krwf = fvc::interpolate(krw,"krw");
	Mwf = Kf *krwf/muwf;
	krg = pow(1-Se,0.5) * pow(1-pow(Se,1/m_vg),2*m_vg);//Info<<"krg "<<krg<<endl; // OA 12/9 2*m -> 5*m
	krg = max(min(krg,1-mysmall),mysmall); //Info<<"krw "<<krw[0]<<" krg "<<krg[0]<<endl; //krmax = 1
	krgf = fvc::interpolate(krg,"krg");
	Mgf = Kf *krgf/mug;
	Mf = Mwf + Mgf;
	Lwf = Mwf*rhowf;
	Lgf = Mgf*rhog;
	Lf = Lwf + Lgf;
	Fwf = (krwf/muwf) / (krgf/mug + krwf/muwf);
	// fluxes
	phiGr = (Lf * g) & mesh.Sf();Info<<"PhiGr(8) "<<phiGr[8]<<" Lwf/Lf "<<Lwf[8]/Lf[8]<<endl;
	pc == pc0*pow(pow(Se,-1/m_vg)-1,1/n_vg); //Info<<"pc "<<pc<<endl;//Se=(1+(p/pc0)^n)^-m; Se^-1/m-1=(p/pc0)^n
	gradpc = fvc::interpolate(fvc::grad(pc,"pc"));//
	phiPc = (Mwf * gradpc) & mesh.Sf(); //fvc::snGrad(pc)

