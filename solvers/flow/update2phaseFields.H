
sw = max(sw_min,min(sw,sw_max));
//- relative permeability computation
	Se = min(max((sw-sw_min)/(sw_max-sw_min),1e-9),1-VSMALL);//Info<<"Se "<<Se[0]<<" "<<Se[1]<<endl;
	//Info<<"se done"<<endl;
	// water
	krw = 1 * pow(Se,0.5) * pow(1-pow(1-pow(Se,1/m_vg),m_vg),2);
	//dkrwdS = 0.5 * (1-pow((1-pow(Se,1/m_vg)),m_vg));
	//dkrwdS *= ( 4 * pow(Se,1/m_vg-1/2) * pow( (1-pow(Se,1/m_vg)) , m_vg-1)) - ( pow((1-pow(Se,1/m_vg)),m_vg) -1) / pow(Se,0.5);
	//dkrwdS *= 1/(sw_max - sw_min);
	//Info<<"kr done"<<endl;
	//gas;
	krg = pow(1-Se,0.5) * pow(1-pow(Se,1/m_vg),2*m_vg);//Info<<"krg "<<krg<<endl; // OA 12/9 2*m -> 5*m
	//krg = 1.3978 - 3.7694*sw + 12.7090*pow(sw,2) - 20.6420*pow(sw,3) + 10.3090*pow(sw,4); //tough formulation
	//krg = 1-(sw-0.2)*0.5-0.623*exp(-10.5*pow(sw-0.94,2.));
	krg = max(min(krg,.99999),1e-9); //Info<<"krw "<<krw[0]<<" krg "<<krg[0]<<endl; //krmax = 1
	//dkrgdS = - pow((1-pow(Se,1/m_vg)),2*m_vg-1) ;
	//dkrgdS *= (-5*pow(Se,1/m_vg+1)+4*pow(Se,1/m_vg)+Se);
	//dkrgdS *= 1/(2*pow((1-Se),0.5)*Se);
	//dkrgdS *=  1/(sw_max - sw_min);

// capillary pressure
	//VG
	pc == pc0*pow(pow(Se,-1/m_vg)-1,1/n_vg); //Se=(1+(p/pc0)^n)^-m; Se^-1/m-1=(p/pc0)^n
	//Info<<" pc "<<pc;
	//dpcdS = -pc0/(n_vg*m_vg)*pow(Se,-1/m_vg - 1); //a b x^(a - 1);a=-1/m,b=1/n
	//dpcdS *= pow(pow(Se,-1/m_vg)-1,1/n_vg-1) ;//(x^a - 1)^(b - 1)
	Info<<" p ";for (i=0;i<24;i++) {Info<<p[i]<<" ";};Info<<"\n";
	Info<<" 1-Se ";for (i=0;i<24;i++) {Info<<1-Se[i]<<" ";};Info<<"\n";
	Info<<" pc ";for (i=0;i<24;i++) {Info<<pc[i]<<" ";};Info<<"\n";
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
//Info<<"dpcpdS done"<<endl;
krwf = fvc::interpolate(krw,"krw");
krgf = fvc::interpolate(krg,"krg");

//- mobility computation 
Mwf = Kf *krwf/muwf;
Mf = Mwf + Kf*krgf/mug;
Info<<"rhog "<<rhog<<" mug "<<mug<<"krw 0 "<<krwf[0]<<" Mf(0) "<<Mf[0]<<endl;
Lwf = Kf*rhowf*krwf/muwf;
Lf = Lwf + Kf*rhog*krgf/mug;
//Info<<"Lwf 6 "<<Lwf[6]<<"Lf 6 "<<Lf[6]<<endl;
Fwf = (krwf/muwf) / (krgf/mug + krwf/muwf);
//Info<<"Mf done"<<endl;
//- compute fluxes (will be used for p Eqn)
phiGr = (Lf * g) & mesh.Sf();//Info<<"phiGr "<<"PhiGr(10) "<<phiGr[10]<<endl;
gradpc = fvc::interpolate(fvc::grad(pc,"pc"));//
phiPc = (Mwf * gradpc) & mesh.Sf();//Info<<"PhiPc(10) "<<phiPc[10]<<endl;
