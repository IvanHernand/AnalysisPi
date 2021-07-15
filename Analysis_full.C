

// Includes here 
#include "TrReco.h"
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
//////////////////////////////////////////////////////////////////////////////

#define D3_TIGC(name,ecal) if (E.BinaHTrig_VT_N>0){	\
    H.Fill(#name,E.Tpos,TAcc,ecal,1.0);				\
  }


#define D3_PE(name,ecal) if (E.PreScaleTrig_VT_N>0 || E.EarlyTrig_VT_N>0) {		\
      if (E.EarlyTrig_VT_N>0) H.Fill(#name,E.Tpos,TAcc,ecal,1.0);	\
      else H.Fill(#name,E.Tpos,TAcc,ecal,16.0);			\
    }

//////////////////////////////////////////////////////////////////////////////

class Analysis2010_2: public ModuleClass{
	public :
		Analysis2010_2()		{};
		~Analysis2010_2()		{};
		bool	Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog) ;
		bool	Process(EventClass &E, HistogramFactory &H);

	private :
                log4cpp::Category *Log;
		string Name;
                int r1;
                int r2;
		double z_wc1_1, z_wc3_1, z_wc3_2, z_wc3_3;

  double TOF_pi_min1, TOF_pi_min2, TOF_pi_min3 ;
  double TOF_pi_max1,  TOF_pi_max2, TOF_pi_max3;
  
  double Q_pi_min;
  double Q_pi_max;
  double E_B2_min;
  double E_B2_max;
  double min1;

  double QcutB1_max;
  double QcutB1_min;
  double QcutB2_max;
  double QcutB2_min;

  double CA,CB;
  double WC3X;
  double WC3Y;
  double phi;
  TrRecoClass TrReco;
  char file_temp[100];

};

bool Analysis2010_2::Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog)
{
	Log	= TmpLog;
	Log->info( "Register Analysis2010");

	//	 --------- Special list of cut in this class for the global histograms ---------		//
	H.AddCut(Name);
	
	if (not Conf.read<bool>("RecoTest/Do"))
	{
		Log->info( "RecoTest turned OFF");
		return false ;
	}
	r1 = Conf.read<int>("RecoTest/run1");
	r2 = Conf.read<int>("RecoTest/run2");
	
        z_wc1_1 = Conf.read<double>("WC1_1/z");
        z_wc3_1 = Conf.read<double>("WC3_1/z");
	z_wc3_2 = Conf.read<double>("WC3_2/z");
        z_wc3_3 = Conf.read<double>("WC3_3/z");

	TOF_pi_min1		= Conf.read<double>("Analysis2010/TOF_Min1");
	TOF_pi_max1		= Conf.read<double>("Analysis2010/TOF_Max1");
	
	TOF_pi_min2		= Conf.read<double>("Analysis2010/TOF_Min2");
	TOF_pi_max2		= Conf.read<double>("Analysis2010/TOF_Max2");
	
	TOF_pi_min3		= Conf.read<double>("Analysis2010/TOF_Min3");
	TOF_pi_max3		= Conf.read<double>("Analysis2010/TOF_Max3");
	
	
	Q_pi_min		= Conf.read<double>("Analysis2010/Q_pi_Min");
	Q_pi_max		= Conf.read<double>("Analysis2010/Q_pi_Max");
	
	E_B2_min		= Conf.read<double>("Analysis2010/E_B2_Min");
	E_B2_max		= Conf.read<double>("Analysis2010/E_B2_Max");

	QcutB1_min		= Conf.read<double>("Analysis2010/QcutB1_Min");
	QcutB2_min		= Conf.read<double>("Analysis2010/QcutB2_Min");
	
	QcutB1_max		= Conf.read<double>("Analysis2010/QcutB1_Max");
	QcutB2_max		= Conf.read<double>("Analysis2010/QcutB2_Max");

	CA       		= Conf.read<double>("Analysis2010/CA");
	CB		        = Conf.read<double>("Analysis2010/CB");
	TrReco.Init(Conf);

        H.DefineTH2D( "EnergySpectrum","Proton_0","", 500,0.0,100.0,500,0.0,10.0);((TH2D*)gDirectory->Get(file_temp))->Sumw2();
        H.DefineTH2D( "EnergySpectrum","Proton_1","", 500,0.0,100.0,500,0.0,10.0);((TH2D*)gDirectory->Get(file_temp))->Sumw2();
	H.DefineTH2D( "EnergySpectrum","Proton_2","", 500,0.0,100.0,500,0.0,10.0);((TH2D*)gDirectory->Get(file_temp))->Sumw2();

		
		
    return true;
}

bool Analysis2010_2::Process(EventClass &E, HistogramFactory &H)
{

   TrReco.Process(E);
  //////////////////////////////////////////////////////////////////////////////////
  
  PTrack& trkWC12 = E.trks[0].GetTrk(0);
  PTrack& trkS12  = E.trks[1].GetTrk(0); 
  
  //////////////////////////////////////////////////////////////////////////////////

  double zWC1 = -112.548; //z at WC1_2
  double zWC2 = -74.409;  //z at WC2_2
  double zWC3 =  56.755;  //z at WC3_2
  double zB2 = -31.55;
  double B2_X = 999.;
  double B2_Y = 999.;
  double zPro = 0.0;

  double TAcc  = 999;
  double xWC3 = 999;
  double yWC3 = 999;
  double x_low, x_high, y_low, y_high;
  double prompt_low, prompt_high;
  double av[4];
  double muhit1f, muhit2f, muhit3f, muhit4f;
  double trighitj_low, trighit_cutValue1, trighit_cutValue2;
  double BinaScale_Q, BinaScale_PH, B1scale, B2scale;
  double T1scale, T2scale;
  double B1_low, B1_high, B2_low, B2_high;

  // Default factors

  BinaScale_Q = 0.9813*0.9885/0.9786;
  BinaScale_PH = 0.9974*1.004/0.9939;
  T1scale = 0.925;
  T2scale = 0.908;
  B1_low = 3.40;
  B1_high = 5.25;
  B2_low = 1.60;
  B2_high = 2.80;

  x_low = -22; x_high = 16; y_low = -12; y_high = 18; //WC12
  prompt_low = -4395; prompt_high = -4375; //T1prompt, TrCons
  av[0] = 7.9; av[1] = 8.7; av[2] = 9.1; av[3] = 9.1; //T1prompt
  trighit_cutValue1 = 1020; trighit_cutValue2 = 1040;

  muhit1f = 8000.; muhit2f = 16400.; muhit3f = 16600.; muhit4f = 17850.; //MuHit
  x_low = -22; x_high = 16; y_low = -12; y_high = 18; //WC12
  trighitj_low = -3797.0; //trigcut

  if(E.runNo>=31000 && E.runNo<42250)   {
    x_low = -22; x_high = 16; y_low = -12; y_high = 18; //WC12
      prompt_low = -4395; prompt_high = -4375; //T1prompt, TrCons
      av[0] = 7.9; av[1] = 8.7; av[2] = 9.1; av[3] = 9.1; //T1prompt
      trighit_cutValue1 = 1020; trighit_cutValue2 = 1040;
    }
  if(E.runNo>=42250 && E.runNo<45819)   {
    x_low = -22; x_high = 16; y_low = -12; y_high = 18; //WC12
      prompt_low = -4409; prompt_high = -4390; //T1prompt, TrCons
      av[0] = 7.3; av[1] = 8.1; av[2] = 8.4; av[3] = 8.3; //T1prompt
      trighit_cutValue1 = 1030; trighit_cutValue2 = 1050;
    }
  if(E.runNo>=45819 && E.runNo<46816)   {
    x_low = -22; x_high = 16; y_low = -12; y_high = 18; //WC12
      prompt_low = -4410; prompt_high = -4380; //T1prompt, TrCons
      av[0] = 2.9; av[1] = 3.7; av[2] = 3.9; av[3] = 3.8; //T1prompt
      trighit_cutValue1 = 1030; trighit_cutValue2 = 1050;
    }
  if(E.runNo>=46816 && E.runNo<47133)   {
    x_low = -22; x_high = 16; y_low = -12; y_high = 18; //WC12
      prompt_low = -4390; prompt_high = -4370; //T1prompt, TrCons
      av[0] = 3.3; av[1] = 5.2; av[2] = 3.9; av[3] = 5.1; //T1prompt
      trighit_cutValue1 = 1020; trighit_cutValue2 = 1040;
    }
  if(E.runNo>=47133 && E.runNo<49095)   {
    x_low = -22; x_high = 16; y_low = -12; y_high = 18; //WC12
      prompt_low = -4390; prompt_high = -4370; //T1prompt, TrCons
      av[0] = 7.2; av[1] = 9.0; av[2] = 7.8; av[3] = 9.1; //T1prompt
      trighit_cutValue1 = 1020; trighit_cutValue2 = 1040;
    }
  if(E.runNo>=49095 && E.runNo<52006)   {
    x_low = -22; x_high = 16; y_low = -12; y_high = 18; //WC12
      prompt_low = -4399; prompt_high = -4380; //T1prompt, TrCons
      av[0] = 11.2; av[1] = 13.1; av[2] = 11.8; av[3] = 13.1; //T1prompt
      trighit_cutValue1 = 1020; trighit_cutValue2 = 1040;
      B1scale = 0.948;
      B2scale = 0.971;
    }
  //Temporary scaling
   if(E.runNo>=49095 && E.runNo<50952)   {
     BinaScale_Q = 0.9813*0.9885/0.9786;
     BinaScale_PH = 0.9974*1.004/0.9939;
     T1scale = 0.925;
     T2scale = 0.908;
     B1_low = 3.40;
     B1_high = 5.25;
     B2_low = 1.60;
     B2_high = 2.80;
   }
   if(E.runNo>=50952 && E.runNo<52006)   {
     BinaScale_Q = 0.9823*0.9885/0.9786;
     BinaScale_PH = 0.9928*1.004/0.9939;
     T1scale = 0.942;
     T2scale = 0.908;
     B1_low = 3.35;
     B1_high = 5.20;
     B2_low = 1.65;
     B2_high = 2.85;
     }

  //2009 Ivan
  if(E.runNo>=10025 && E.runNo<25750)   {
      x_low = -18;  x_high = 22;
      y_low = -15;  y_high = 19;
      QcutB1_min = 0.78; QcutB1_max = 0.98;
      QcutB2_min = 0.74; QcutB2_max = 0.96;
      
      // prompt_low = -4399; prompt_high = -4380; //T1prompt, TrCons
      // av[0] = 11.9; av[1] = 13.6; av[2] = 12.6; av[3] = 13.6; //T1prompt
      // trighit_cutValue1 = 1010; trighit_cutValue2 = 1040;
      // B1scale = 0.965;
      // B2scale = 0.990;   
      
      //Temporary scaling for Bina
     // BinaScale_Q = 0.9781*0.9885/0.9786;
     // BinaScale_PH = 0.9905*1.004/0.9939;
     // T1scale = 1.000;
     // T2scale = 0.881;
     
     B1_low = 3.55;
     B1_high = 5.45;
     B2_low = 1.9;
     B2_high = 3.0;
    }

  //2011
  if(E.runNo>=57418 && E.runNo<61179)   {
      x_low = -20;  x_high = 18;
      prompt_low = -4399; prompt_high = -4380; //T1prompt, TrCons
      av[0] = 11.9; av[1] = 13.6; av[2] = 12.6; av[3] = 13.6; //T1prompt
      trighit_cutValue1 = 1010; trighit_cutValue2 = 1040;
      B1scale = 0.965;
      B2scale = 0.990;      
      //Temporary scaling for Bina
     BinaScale_Q = 0.9781*0.9885/0.9786;
     BinaScale_PH = 0.9905*1.004/0.9939;
     T1scale = 1.000;
     T2scale = 0.881;     
     B1_low = 3.45;
     B1_high = 5.3;
     B2_low = 1.80;
     B2_high = 3.00;
    }

  //2012
  if(E.runNo>=62491 && E.runNo<81560) {
    x_low = -23; x_high = 19; y_low = -17; y_high = 19; //WC12
    trighit_cutValue1 = 1000;  trighit_cutValue2 = 1040;
    trighitj_low = -3820.0;
    // BinaScale_Q = 0.9786;
    // BinaScale_PH = 0.9939;
    // for nAl2O3
    B1scale = 0.940;
    B2scale = 0.971;
    BinaScale_Q = 0.9885;
    BinaScale_PH = 1.004;
    T1scale = 0.961;
    T2scale = 0.887;
     B1_low = 3.35;
     B1_high = 5.20;
     B2_low = 1.70;
     B2_high = 2.90;
  }
  if(E.runNo>=62491 && E.runNo<70025) {
    prompt_low = -4399;  prompt_high = -4350;
    av[0] = 11.2; av[1] = 13.1; av[2] = 11.8; av[3] = 13.1;
  }
  if(E.runNo>=70025 && E.runNo<81560) {
    prompt_low = -4375;  prompt_high = -4350;
    av[0] = 12.1; av[1] = 13.6; av[2] = 12.8; av[3] = 13.6;
  }
  E.Cal_eB1 = E.Cal_eB1*B1scale;
  E.Cal_eB2 = E.Cal_eB2*B2scale;
  E.Cal_eT1 = E.Cal_eT1*T1scale;
  E.Cal_eT2 = E.Cal_eT2*T2scale;
  E.Cal_eBina = E.Cal_eBina*BinaScale_Q;
  E.Cal_eBinaPH = E.Cal_eBinaPH*BinaScale_PH;
  
  double CsI = E.USIcal + E.USOcal + E.DSIcal + E.DSOcal;
  double DS = E.Cal_eS3 + E.Cal_eT1 + E.Cal_eT2;
  double Ecal = E.Cal_eBina + CsI;
  double Ecal_PH = E.Cal_eBinaPH + CsI;
  double Ecal_DS = Ecal + DS;
  double Ecal_PHDS = Ecal_PH + DS;

  //Minimum radius acceptance
  //xWC3 = -0.65;
  //yWC3 = -2.;
  xWC3 = 0.;
  yWC3 = 0.;
  TAcc = 999;
  WC3X = 999.;
  WC3Y = 999.;
  phi = 999;
  int status = 999;
  int Npts = 999;
  double acctmp = 999;
  double Tg = 999;
  int isave = 0;
  double Zv = 999.;
  double Zv1 = 999.;
  double Zv2 = 999.;
  double S3WC3X0 = 999.;
  int Ntrk = TrReco.trks[2]->GetN();
  if (Ntrk>0){
    for (int i=0;i<Ntrk;i++){
      PTrack& t = TrReco.trks[2]->GetTrk(i);
      acctmp = sqrt((t.tx*zWC3+t.x0+xWC3)*(t.tx*zWC3+t.x0+xWC3)+(t.ty*zWC3+t.y0+yWC3)*(t.ty*zWC3+t.y0+yWC3));
      if (acctmp <= TAcc) {
	TAcc = acctmp;
	isave = i;
	WC3X = t.tx*zWC3+t.x0+xWC3;
	WC3Y = t.ty*zWC3+t.y0+yWC3;
	status = TrReco.trks[2]->status;
	Npts = t.npts;
	Tg = sqrt((t.x0)*(t.x0)+(t.y0)*(t.y0));
      } 
    }
       PTrack& trkDw = TrReco.trks[2]->GetTrk(isave);
       PTrack& trkWCS3  = E.trks[2].GetTrk(isave);
       PTrack& trkWC3  = E.trks[2].GetTrk(0);
       S3WC3X0 = trkDw.x0;
   int Ntrk12 = TrReco.trks[1]->GetN();
   if(Ntrk12>0) {
       PTrack& trkS12 = TrReco.trks[1]->GetTrk(0);
       PTrack& trkSS12  = E.trks[1].GetTrk(0); 
       Zv = -((trkS12.x0-trkDw.x0)*(trkS12.tx-trkDw.tx)+
		 (trkS12.y0-trkDw.y0)*(trkS12.ty-trkDw.ty))/
      ((trkS12.tx-trkDw.tx)*(trkS12.tx-trkDw.tx)+
       (trkS12.ty-trkDw.ty)*(trkS12.ty-trkDw.ty));
	phi = atan2(WC3Y,WC3X);

	// Test Zv for 2009 Ivan

	Zv1 = -((trkSS12.x0-trkWCS3.x0)*(trkSS12.tx-trkWCS3.tx)+
		 (trkSS12.y0-trkWCS3.y0)*(trkSS12.ty-trkWCS3.ty))/
      ((trkSS12.tx-trkWCS3.tx)*(trkSS12.tx-trkWCS3.tx)+
       (trkSS12.ty-trkWCS3.ty)*(trkSS12.ty-trkWCS3.ty));

        Zv2 = -((trkSS12.x0-trkWC3.x0)*(trkSS12.tx-trkWC3.tx)+
		 (trkSS12.y0-trkWC3.y0)*(trkSS12.ty-trkWC3.ty))/
      ((trkSS12.tx-trkWC3.tx)*(trkSS12.tx-trkWC3.tx)+
       (trkSS12.ty-trkWC3.ty)*(trkSS12.ty-trkWC3.ty));
   }
  }
    
  //Blinding cut and data integrity cuts
  bool data_integrity = 
    E.B1_1_WF_E==65535 && //Copper End-code
    E.B1_2_WF_E==65535 && 
    E.B1_3_WF_E==65535 && 
    E.B1_4_WF_E==65535 && 
    E.T1_1_WF_E==65535 && 
    E.T1_2_WF_E==65535 && 
    E.T1_3_WF_E==65535 && 
    E.T1_4_WF_E==65535 &&
    E.CprSamples[0]>=610 && //COPPER # of samples
    E.CprSamples[1]>=610 && 
    E.CprSamples[2]>=610 && 
    E.CprSamples[3]>=610 &&
    E.blindf==1;

  //Pion identification cuts and triggers selection (->no calibration triggers)

  bool pion = 
    E.Cal_eB1>B1_low && 
    E.Cal_eB1<B1_high && 
    E.Cal_eB2>B2_low && 
    E.Cal_eB2<B2_high;
  bool trigger = 
    E.ECalibTrig_VT_N==0 && 
    E.BinaSumTrig_VT_N==0 && 
    E.XeTrig_VT_N==0 && 
    E.PulserTrig_VT_N==0 &&
    E.CsISumTrig_VT_N==0 &&
    E.BeamPreTrig_VT_N==0 &&
    E.PhysTrig_VT_N==1;
  bool pion_trigger = pion && trigger;

  //WC12 acceptance definition
  bool WC12 = 0;
  if (E.trks[0].GetN()>0){
    WC12 = 
      (trkWC12.tx*zWC1+trkWC12.x0)<x_high && (trkWC12.tx*zWC1+trkWC12.x0)>x_low &&
      (trkWC12.ty*zWC1+trkWC12.y0)<y_high && (trkWC12.ty*zWC1+trkWC12.y0)>y_low &&
      (trkWC12.tx*zWC2+trkWC12.x0)<x_high && (trkWC12.tx*zWC2+trkWC12.x0)>x_low &&
      (trkWC12.ty*zWC2+trkWC12.y0)<y_high && (trkWC12.ty*zWC2+trkWC12.y0)>y_low &&
      trkWC12.x0>-20 && trkWC12.x0<20 &&
      trkWC12.y0>-20 && trkWC12.y0<20;

    B2_X = trkWC12.tx*zB2+trkWC12.x0;
    B2_Y = trkWC12.ty*zB2+trkWC12.y0;
   }

  bool WC3 = TAcc < 60.0;

  //Proton rejection cut
  Double_t min1, min2,min3, min4;  
  min1 = min(E.Cal_eT1/(0.318*1.03) ,E.Cal_eT2/(0.25*1.03));
  min4 = min(min1,E.Cal_eS3/(0.057*2.33));
  bool proton = min4<CA*E.Cal_eBina+CB;
  
  //Pileup cut in B1,B2 and Chi2 cut for COPPER fit in B1 and T1

  bool B1B2_pileup = 
    (E.B1_1_WF_N==1 || E.B1_2_WF_N==1 || E.B1_3_WF_N==1 || E.B1_4_WF_N==1) &&
    (E.B2_1_WF_N==1 || E.B2_2_WF_N==1 || E.B2_3_WF_N==1 || E.B2_4_WF_N==1) &&
    QcutB1_min<E.B1_1_WF_Q[0]/E.B1_1_WF_Qw[0] && E.B1_1_WF_Q[0]/E.B1_1_WF_Qw[0]<QcutB1_max &&
    QcutB1_min<E.B1_2_WF_Q[0]/E.B1_2_WF_Qw[0] && E.B1_2_WF_Q[0]/E.B1_2_WF_Qw[0]<QcutB1_max &&
    QcutB1_min<E.B1_3_WF_Q[0]/E.B1_3_WF_Qw[0] && E.B1_3_WF_Q[0]/E.B1_3_WF_Qw[0]<QcutB1_max &&
    QcutB1_min<E.B1_4_WF_Q[0]/E.B1_4_WF_Qw[0] && E.B1_4_WF_Q[0]/E.B1_4_WF_Qw[0]<QcutB1_max &&
    QcutB2_min<E.B2_1_WF_Q[0]/E.B2_1_WF_Qw[0] && E.B2_1_WF_Q[0]/E.B2_1_WF_Qw[0]<QcutB2_max &&
    QcutB2_min<E.B2_2_WF_Q[0]/E.B2_2_WF_Qw[0] && E.B2_2_WF_Q[0]/E.B2_2_WF_Qw[0]<QcutB2_max &&
    QcutB2_min<E.B2_3_WF_Q[0]/E.B2_3_WF_Qw[0] && E.B2_3_WF_Q[0]/E.B2_3_WF_Qw[0]<QcutB2_max &&
    QcutB1_min<E.B2_4_WF_Q[0]/E.B2_4_WF_Qw[0] && E.B2_4_WF_Q[0]/E.B2_4_WF_Qw[0]<QcutB1_max;
  bool T1_chi2 = 
    E.T1_1_WF_Chi2[0]>=0 && E.T1_2_WF_Chi2[0]>=0 && E.T1_3_WF_Chi2[0]>=0 && E.T1_4_WF_Chi2[0]>=0;
  bool B1_chi2 = E.B1_1_WF_Chi2[0]>=0 && E.B1_2_WF_Chi2[0]>=0 && E.B1_3_WF_Chi2[0]>=0 && E.B1_4_WF_Chi2[0]>=0;
  bool Pileup = B1B2_pileup && T1_chi2 && B1_chi2;

  //For T1 pilepup
  double qfullavg1 = (E.T1_2_WF_Qfull + E.T1_3_WF_Qfull + E.T1_4_WF_Qfull)/3.0;
  double Afitavg1 =  (E.T1_2_WF_Afit[0]+E.T1_3_WF_Afit[0]+E.T1_4_WF_Afit[0])/3.0;

  //T1 pileup cut
  bool T1Pileup = (E.T1_1_WF_N==1 || E.T1_2_WF_N==1 || E.T1_3_WF_N==1 || E.T1_4_WF_N==1) || 
    (qfullavg1 / Afitavg1 <= 7.8 + 1.17e-4*Afitavg1*Afitavg1-2.38e-3*Afitavg1);

  double cutthresh;
  if(Afitavg1 < 76)
    cutthresh = 7.6;
  else
    cutthresh = 0.0169∗Afitavg1 + 7.6;
  
  bool T1Pileup_v2 = (E.T1_1_WF_N==1 || E.T1_2_WF_N==1 || E.T1_3_WF_N==1 || E.T1_4_WF_N==1) || 
    (qfullavg1 / Afitavg1 <= cutthresh);

  bool simpleT1Pileup = E.T1_1_WF_N == 1 || E.T1_2_WF_N == 1 || E.T1_3_WF_N == 1 || E.T1_4_WF_N == 1;

    bool B1prompt;

    if(E.runNo<=25751){ //2009

      B1prompt =  // for 2009 Ivan
      E.B1_1_WF_t[0]>-1319.5 &&  E.B1_1_WF_t[0]<-1269.5 
      && E.B1_2_WF_t[0]>-1317.5 &&  E.B1_2_WF_t[0]<-1267.5
      && E.B1_3_WF_t[0]>-1319.5 &&  E.B1_3_WF_t[0]<-1269.5
      && E.B1_4_WF_t[0]>-1321.5 &&  E.B1_4_WF_t[0]<-1275.5;

    }

    else{

      B1prompt =  // 2010,2011,2012
      E.B1_1_WF_t[0]>-1380 &&  E.B1_1_WF_t[0]<-1340 
      && E.B1_2_WF_t[0]>-1380 &&  E.B1_2_WF_t[0]<-1340
      && E.B1_3_WF_t[0]>-1380 &&  E.B1_3_WF_t[0]<-1340
      && E.B1_4_WF_t[0]>-1380 &&  E.B1_4_WF_t[0]<-1340;

    }
   
  //No hits in pre region for B1,B2,Tg,T1,T2 cut
  // bool pion_prepu = 
  //   E.B1_1_WFPre_N==0 && E.B1_2_WFPre_N==0 &&  E.B1_3_WFPre_N==0 && E.B1_4_WFPre_N==0 &&
  //   E.B2_1_WFPre_N==0 && E.B2_2_WFPre_N==0 &&  E.B2_3_WFPre_N==0 && E.B2_4_WFPre_N==0 &&
  //   E.Tg_1_WFPre_N==0 && E.Tg_2_WFPre_N==0 &&  E.Tg_3_WFPre_N==0 && E.Tg_4_WFPre_N==0;
  // bool positron_prepu = 
  //   E.T1_1_WFPre_N==0 && E.T1_2_WFPre_N==0 &&  E.T1_3_WFPre_N==0 && E.T1_4_WFPre_N==0 &&
  //   E.T2_1_WFPre_N==0 && E.T2_2_WFPre_N==0 &&  E.T2_3_WFPre_N==0 && E.T2_4_WFPre_N==0;
  // bool PrePileup = pion_prepu && positron_prepu;
  
  //T1 prompt pulses cut
  bool T1prompt = 1;
  // find B1_1 VT prompt index
  int b1_1_pmpt_idx = -1;
  for(int j=0; j<E.B1_1_VT_N; j++) {
    if(prompt_low<E.B1_1_VT_t[j]-E.Upstr_VT_t[0] && E.B1_1_VT_t[j]-E.Upstr_VT_t[0]<prompt_high) {
      b1_1_pmpt_idx = j;
      break;
    }
  }
  if(b1_1_pmpt_idx>-1) {
     int   nn[4] = {E.T1_1_VT_N,E.T1_2_VT_N,E.T1_3_VT_N,E.T1_4_VT_N};
     Int_t *t[4] = {E.T1_1_VT_t,E.T1_2_VT_t,E.T1_3_VT_t,E.T1_4_VT_t};
     for(int i=0; i<4; i++) {
       for(int j=0; j<nn[i]; j++) {
	 if(fabs(t[i][j]-E.B1_1_VT_t[b1_1_pmpt_idx]-av[i])*0.625<2) { //9
	   T1prompt = 0;
	 }
       }
     }
  } else T1prompt = 0;

  //Trigger Consistency 
  bool TrCons = 1;
  int   nn[4] = {E.B1_1_VT_N,E.B1_2_VT_N,E.B1_3_VT_N,E.B1_4_VT_N};
  Int_t *t[4] = {E.B1_1_VT_t,E.B1_2_VT_t,E.B1_3_VT_t,E.B1_4_VT_t};
  bool ok = true;
  for(int i=0; i<4; i++) {
    bool hit = false;
    for(int j=0; j<nn[i]; j++)
      {
	if(prompt_low<t[i][j]-E.Upstr_VT_t[0]&&t[i][j]-E.Upstr_VT_t[0]<prompt_high) { 
	hit = true;
	break;
	}
    }
    if(!hit) {
      ok = false;
      break;
    } 
  }
  if(!ok) {
    TrCons = 0;
  }

  //hit triggering cut (in T1)
  int trighit1 = -1;
  int trighit2 = -1;
  int trighit3 = -1;
  int trighit4 = -1;
  int rhit = -1;
  int nhit = 0;
  double maxPH = 0;
 
  for(int j = 0; j < E.e_time_VT_N; j++)
    {
      if (E.Dwnstr_VT_N > 0 && E.PhysTrig_VT_N == 1) {
	if(E.e_time_VT_t[j]
	  - E.Dwnstr_VT_t[0] > trighitj_low && E.e_time_VT_t[j] -
	  E.Dwnstr_VT_t[0] < -3757.0) rhit = j;
      }
    }
  for (int i = 0; i < E.T1_1_WF_N; i++)
    {
       if (rhit >= 0 && E.T1_1_WF_t[i] - (E.e_time_VT_t[rhit] -
					 E.Upstr_VT_t[0]) * 0.625 > trighit_cutValue1 && E.T1_1_WF_t[i] -
	  (E.e_time_VT_t[rhit] - E.Upstr_VT_t[0]) * 0.625 < trighit_cutValue2) nhit++; //1050
    }
  if (nhit > 0)
  {
	for (int i = 0; i < nhit; i++)
	{
		if (rhit >= 0 && E.T1_1_WF_t[i] - (E.e_time_VT_t[rhit] - E.Upstr_VT_t[0]) * 0.625 > trighit_cutValue1  &&
		    E.T1_1_WF_t[i] - (E.e_time_VT_t[rhit] - E.Upstr_VT_t[0]) * 0.625 < trighit_cutValue2 &&
		    E.T1_1_WF_PH[i] > maxPH) {trighit1 = i; maxPH = E.T1_1_WF_PH[i];}
	}
  }
  maxPH = 0;
  nhit = 0;
    
  for (int i = 0; i < E.T1_2_WF_N; i++)
    {
      if (rhit >= 0 && E.T1_2_WF_t[i] - (E.e_time_VT_t[rhit] -
					 E.Upstr_VT_t[0]) * 0.625 > trighit_cutValue1 && E.T1_2_WF_t[i] -
	  (E.e_time_VT_t[rhit] - E.Upstr_VT_t[0]) * 0.625 < trighit_cutValue2) nhit++;
    }
   if (nhit > 0)
  {
	for (int i = 0; i < nhit; i++)
	{
		if (rhit >= 0 && E.T1_2_WF_t[i] - (E.e_time_VT_t[rhit] - E.Upstr_VT_t[0]) * 0.625 > trighit_cutValue1 &&
		    E.T1_2_WF_t[i] - (E.e_time_VT_t[rhit] - E.Upstr_VT_t[0]) * 0.625 < trighit_cutValue2 &&
		    E.T1_2_WF_PH[i] > maxPH) {trighit2 = i; maxPH = E.T1_2_WF_PH[i];}
	}
  }
  maxPH = 0;
  nhit = 0;
   
  for (int i = 0; i < E.T1_3_WF_N; i++)
    {
      if (rhit >= 0 && E.T1_3_WF_t[i] - (E.e_time_VT_t[rhit] -
					 E.Upstr_VT_t[0]) * 0.625 > trighit_cutValue1 && E.T1_3_WF_t[i] -
	  (E.e_time_VT_t[rhit] - E.Upstr_VT_t[0]) * 0.625 < trighit_cutValue2) nhit++;
    }
    if (nhit > 0)
  {
	for (int i = 0; i < nhit; i++)
	{
		if (rhit >= 0 && E.T1_3_WF_t[i] - (E.e_time_VT_t[rhit] - E.Upstr_VT_t[0]) * 0.625 > trighit_cutValue1 &&
		    E.T1_3_WF_t[i] - (E.e_time_VT_t[rhit] - E.Upstr_VT_t[0]) * 0.625 < trighit_cutValue2 &&
		    E.T1_3_WF_PH[i] > maxPH) {trighit3 = i; maxPH = E.T1_3_WF_PH[i];}
	}
  }
  maxPH = 0;
  nhit = 0;
  
  for (int i = 0; i < E.T1_4_WF_N; i++)
    {
      if (rhit >= 0 && E.T1_4_WF_t[i] - (E.e_time_VT_t[rhit] -
					 E.Upstr_VT_t[0]) * 0.625 > trighit_cutValue1 && E.T1_4_WF_t[i] -
	  (E.e_time_VT_t[rhit] - E.Upstr_VT_t[0]) * 0.625 < trighit_cutValue2) nhit++;
    }
    if (nhit > 0)
  {
	for (int i = 0; i < nhit; i++)
	{
		if (rhit >= 0 && E.T1_4_WF_t[i] - (E.e_time_VT_t[rhit] - E.Upstr_VT_t[0]) * 0.625 > trighit_cutValue1  &&
		    E.T1_4_WF_t[i] - (E.e_time_VT_t[rhit] - E.Upstr_VT_t[0]) * 0.625 < trighit_cutValue2 &&
		    E.T1_4_WF_PH[i] > maxPH) {trighit4 = i; maxPH = E.T1_4_WF_PH[i];}
	}
  }
  
  bool trigcut = trighit1 == 0 && trighit2 == 0 && trighit3 == 0 && trighit4== 0;
  if(E.runNo < 40969) trigcut = true;
 
  //T1-T2 coincidence cut
  bool T1T2coinc = false;
  vector<double> T1T2_All;
  double T1T2MINIMUM = 10000;
  for (int T1SUM=0;T1SUM<E.T1_Sum_VT_N;T1SUM++){
    for (int T2SUM=0;T2SUM<E.T2_Sum_VT_N;T2SUM++){
      
      double MINI = (E.T1_Sum_VT_t[T1SUM]-E.T2_Sum_VT_t[T2SUM])*0.625+24;
      T1T2_All.push_back(fabs(MINI));
      
      if (fabs(MINI)<fabs(T1T2MINIMUM))	T1T2MINIMUM=MINI;
    }
  }
  if (T1T2MINIMUM>-20 && T1T2MINIMUM<20) T1T2coinc = true;


  // FINAL
  //////////////////////////////////////////////////////////////////////////////////
  // CUT 9a
  //No hits in pre region for B1,B2,Tg
  bool PrePileup_pion = 
    E.B1_1_WFPre_N==0 && E.B1_2_WFPre_N==0 &&  E.B1_3_WFPre_N==0 && E.B1_4_WFPre_N==0 &&
    E.B2_1_WFPre_N==0 && E.B2_2_WFPre_N==0 &&  E.B2_3_WFPre_N==0 && E.B2_4_WFPre_N==0 &&
    E.Tg_1_WFPre_N==0 && E.Tg_2_WFPre_N==0 &&  E.Tg_3_WFPre_N==0 && E.Tg_4_WFPre_N==0 ;

  // FINAL
  //////////////////////////////////////////////////////////////////////////////////
  // CUT 9b
  //No hits in pre region for T1,T2 cut
  bool PrePileup_positron = 
    E.T1_1_WFPre_N==0 && E.T1_2_WFPre_N==0 &&  E.T1_3_WFPre_N==0 && E.T1_4_WFPre_N==0 &&
    E.T2_1_WFPre_N==0 && E.T2_2_WFPre_N==0 &&  E.T2_3_WFPre_N==0 && E.T2_4_WFPre_N==0;

  bool PrePileup = PrePileup_pion && PrePileup_positron;


  //Muhit cut
  bool muhitcut = true;
  if (E.MuHit_VT_N>0){
    for (int mhit=0;mhit<E.MuHit_VT_N;mhit++){
      if ( ( (E.MuHit_VT_t[mhit]*0.625+E.Tpos>muhit1f) && (E.MuHit_VT_t[mhit]*0.625+E.Tpos<muhit2f) ) || ( (E.MuHit_VT_t[mhit]*0.625+E.Tpos>muhit3f) && (E.MuHit_VT_t[mhit]*0.625+E.Tpos<muhit4f) ) ) {
	muhitcut = false;
      }
    }
  }

  bool pihitcut = true;
  if (E.PiHit_VT_N>0){
    for (int i=0;i<E.PiHit_VT_N;i++){

      if ( ( (E.PiHit_VT_t[i]*0.625+E.Tpos>muhit1f) && (E.PiHit_VT_t[i]*0.625+E.Tpos<muhit2f) )|| ( (E.PiHit_VT_t[i]*0.625+E.Tpos>muhit3f) && (E.PiHit_VT_t[i]*0.625+E.Tpos<muhit4f) ) ){
	pihitcut = false;
      }
    }
  }

  //Post-pileup cut
  bool Post = true;
  if(E.T1_Sum_VT_N>0 && E.B1_Sum_VT_N>0 && Ecal>52.){
    for(int b1n=0; b1n<E.B1_Sum_VT_N; b1n++){
      for(int t1n=0; t1n<E.T1_Sum_VT_N; t1n++){
	double t1b1vt=(E.T1_Sum_VT_t[t1n]-E.B1_Sum_VT_t[b1n])*0.625;
	if(t1b1vt>1335 && t1b1vt<1480 && E.Tpos>420 )
	  Post = false;
      }
    }
  }

  bool Post_DS = true;
  if(E.T1_Sum_VT_N>0 && E.B1_Sum_VT_N>0 && Ecal>54.){
    for(int b1n=0; b1n<E.B1_Sum_VT_N; b1n++){
      for(int t1n=0; t1n<E.T1_Sum_VT_N; t1n++){
	double t1b1vt=(E.T1_Sum_VT_t[t1n]-E.B1_Sum_VT_t[b1n])*0.625;
	if(t1b1vt>1335 && t1b1vt<1480 && E.Tpos>420 )
	  Post_DS = false;
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////

   bool nominal_cuts_L1= (

       // group 1
       data_integrity &&  // 1 & 2
       pion_trigger && // 3 & 4
       WC12 &&// 5
       // //WC3 && // 5

       B1B2_pileup &&// 6
       B1_chi2  // 7
       // B1prompt && // 11
       // proton && // 10 // proton or proton_PH
       // TrCons && // 12
      
       // //group 2
       // T1prompt && // 13
       // T1Pileup && // 6b
       // T1_chi2 && // 7
       // trigcut && // 16
       // T1T2coinc //&& // 17

       //group 3
       //PrePileup_pion &&  // 9
       //PrePileup_positron &&  // 9 
       //muhitcut //&& // 14 
       //pihitcut //&& // 14b only for 2012 dataset
       
       
       // group 4
       // TAcc < 60.0 && // 5b        
       //Ftrig //&& // 15      
       //Post // 19

       );

   bool cuts_pion_trigger= (
			    
       data_integrity &&  // 1 & 2
       pion_trigger

       );

   bool cuts_WC12= (

       data_integrity &&  // 1 & 2
       pion_trigger && // 3 & 4
       WC12

       );

   bool cuts_WC3= (

       data_integrity &&  // 1 & 2
       pion_trigger && // 3 & 4
       WC12 &&// 5
       WC3  // 5

       );

   bool cuts_WC3_wPion= (

       data_integrity &&  // 1 & 2
       trigger && // 3 & 4
       WC12 &&// 5
       WC3  // 5

       );

   bool cuts_B1B2_pileup= (

       data_integrity &&  // 1 & 2
       pion_trigger && // 3 & 4
       WC12 &&// 5
       WC3 && // 5
       B1B2_pileup

       );

   bool cuts_PrePileup= (

       data_integrity &&  // 1 & 2
       pion_trigger && // 3 & 4
       WC12 &&// 5
       WC3 && // 5
       B1B2_pileup &&
       T1Pileup &&
       B1_chi2 &&
       T1_chi2 &&
       PrePileup

       );

   bool cuts_PrePileup_2= (

       data_integrity &&  // 1 & 2
       pion_trigger && // 3 & 4
       WC12 &&// 5
       WC3 && // 5
       B1B2_pileup &&
       T1Pileup_v2 &&
       B1_chi2 &&
       T1_chi2 &&
       PrePileup

       );



   bool cuts_proton = (

       data_integrity &&  // 1 & 2
       pion_trigger && // 3 & 4
       WC12 &&// 5
       WC3 && // 5
       B1B2_pileup &&
       T1Pileup &&
       B1_chi2 &&
       T1_chi2 &&
       PrePileup &&
       proton

       );

   bool cuts_B1prompt = (

       data_integrity &&  // 1 & 2
       pion_trigger && // 3 & 4
       WC12 &&// 5
       WC3 && // 5
       B1B2_pileup &&
       T1Pileup &&
       B1_chi2 &&
       T1_chi2 &&
       PrePileup &&
       proton &&
       B1prompt

       );

   H.Fill("Proton_0",E.Cal_eBina,min4,1.0);

   if(cuts_PrePileup)
     H.Fill("Proton_1",E.Cal_eBina,min4,1.0);

   if(cuts_PrePileup_2)
     H.Fill("Proton_2",E.Cal_eBina,min4,1.0);

      
 
   
  return true;
  
}
  
  
