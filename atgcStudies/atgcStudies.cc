#include "/local/cms/user/gsorrent/antgc_analysis/macros/extra_tools.cc"

R__ADD_INCLUDE_PATH(/local/cms/user/wadud/aNTGCmet/xgboost/include/xgboost/)
R__LOAD_LIBRARY(/local/cms/user/wadud/aNTGCmet/xgboost/lib/libxgboost.so)
#include </local/cms/user/wadud/aNTGCmet/xgboost/include/xgboost/c_api.h>


#ifndef TAUVETOSTUDIES
#define TAUVETOSTUDIES

// Barrel-Endcap transition region
#define BETRetaMin 1.4442
#define BETRetaMax 1.566
#define HBetaMax 1.3920                     // ref. Josh H.
#define ZMASS 91.1876
#define pi7 3.1415927
#define REPORT_EVERY 10000
#define CUTFLOWSTEPS 30


const Double_t ECAL_ETA_BINS[57] = {-5., -3., -2.9, -2.7, -2.5, -2.4, -2.3, -2.2, -2.1, -2., -1.9, -1.8, -1.7, -BETRetaMax, -BETRetaMin, -1.3, -1.2, -1.1, -1., -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,  0., 0.1, 0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, BETRetaMin, BETRetaMax, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.3, 2.4, 2.5, 2.7, 2.9, 3., 5.};
const Double_t ECAL_ABS_ETA_BINS[13] = { 0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, BETRetaMin, BETRetaMax, 1.8, 2.0, 2.2, 2.5};
const Double_t ECAL_EB_ETA_BINS[57] = {-BETRetaMin, -1.35, -1.3, -1.25, -1.2, -1.15, -1.1, -1.05, -1., -0.95, -0.9, -0.85, -0.8, -0.75, -0.7, -0.65, -0.6, -0.55, -0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, BETRetaMin};
const Double_t ECAL_FINE_ETA_BINS[85] = {-5., -3., -2.9, -2.7, -2.5, -2.4, -2.3, -2.2, -2.1, -2., -1.9, -1.8, -1.7, -BETRetaMax, -BETRetaMin, -1.35, -1.3, -1.25, -1.2, -1.15, -1.1, -1.05, -1., -0.95, -0.9, -0.85, -0.8, -0.75, -0.7, -0.65, -0.6, -0.55, -0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, BETRetaMin, BETRetaMax, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.7, 2.9, 3., 5.};
const Double_t ECAL_EB_ABS_ETA_BINS[8] = {0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, BETRetaMin};


////////////////////////////////////////////////// Container for categories //////////////////////////////////////////////////////////////////////////////////////
struct eventType{
	TTree *                             tree = nullptr;

	// Cut efficiency tracking
	Float_t                             lastCutStep = 0.;
	TH1F *                              cutFlowCount = nullptr;
	TH1F *                              cutFlowGenWeight = nullptr;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class atgcStudies{
public:
	atgcStudies(std::string FILELIST, std::string OUTFILE, Float_t XSECTION=-1., std::string MCPILEUPHIST="", std::string DATAPILEUPHIST="/hdfs/cms/user/wadud/anTGC/analysis_data/METv5_pileup/pileup_2017_data.root",
		std::string PFECALCLUS_EFFECTIVE_AREAS="/hdfs/cms/user/wadud/anTGC/analysis_data/effAreas/90pc/phoPFClusEcalIso.txt",
		std::string PFHCALCLUS_EFFECTIVE_AREAS="/hdfs/cms/user/wadud/anTGC/analysis_data/effAreas/90pc/phoPFClusHcalIso.txt",
		std::string TKRISO_EFFECTIVE_AREAS="/hdfs/cms/user/wadud/anTGC/analysis_data/effAreas/90pc/phoTrkSumPtHollowConeDR03.txt",
		std::string PFECALCLUS_PTSCALING="/hdfs/cms/user/wadud/anTGC/analysis_data/isoPtScaling/90pc/phoPFClusEcalIso.txt",
		std::string PFHCALCLUS_PTSCALING="/hdfs/cms/user/wadud/anTGC/analysis_data/isoPtScaling/90pc/phoPFClusHcalIso.txt",
		std::string BDT_PATH="/local/cms/user/gsorrent/antgc_analysis/BDTmodel/aNTGC_photon_BDT_2020_12_14_20_11_56.model");

	~atgcStudies(){
		XGBoosterFree(phoBDT_h);
		std::cout<<"END @ "<<getCurrentTime()<<std::endl;
		std::cout<<"*************************************************************************************************************************************************"<<std::endl;
	};

private:
	Bool_t              isMC = false;
	Bool_t              doPUreweight = false;
	Float_t             xSec = -1.;

	effectiveAreaMap    PFHCALClusEffAreas;
	effectiveAreaMap    PFECALClusEffAreas;
	effectiveAreaMap    TkrEffAreas;

	isoPtScalingMap 	PFHCALClusPtScaling;
	isoPtScalingMap 	PFECALClusPtScaling;

	void                analyze();
	Bool_t              selectEvent();

	Bool_t             photonPassingID(UShort_t phoIndex);

        Float_t            deltaR(Float_t eta1, Float_t phi1, Float_t eta2, Float_t phi2);
        Short_t		   photonIsPrompt(Short_t _phoIndex, Float_t _deltaRmax);// Float_t _relDeltaPtMin, Float_t _relDeltaPtMax);
	TFile *             outFile = nullptr;

	/////////////////////////////////////////// Pileup Reweighting /////////////////////////////////////////////////////////
	PileupReWeighting   puReweighter;
	TH1F                pileupPreweight{"pileupUnweighted", "Unweighted Pileup; True # of Interactions", 200, 0., 200.};
	TH1F                pileupPostweight{"pileupWeighted", "Weighted Pileup; True # of Interactions", 200, 0., 200.};
	TH1F                rhoPreweight{"rhoUnweighted", "Unweighted #rho; #rho", 200, 0., 200.};
	TH1F                rhoPostweight{"rhoWeighted", "Weighted #rho; #rho", 200, 0., 200.};
	TH1F                nvtxPreweight{"nvtxUnweighted", "Unweighted # of Vertices; # of Vertices", 200, 0., 200.};
	TH1F                nvtxPostweight{"nvtxWeighted", "Weighted # of Vertices; # of Vertices", 200, 0., 200.};
                             
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////// Input TTree ///////////////////////////////////////////////////////////
	Bool_t                                  initNtuples(std::string FILELIST);
	TTreeReader                             inputTTreeReader;
	TChain *                                inputTree = nullptr;

	TTreeReaderAnyValue<Int_t>              _run;
	TTreeReaderAnyValue<Long64_t>           _event;
	TTreeReaderAnyValue<UShort_t>           _lumis;
	TTreeReaderAnyValue<UChar_t>            _nVtx;
	TTreeReaderAnyValue<Float_t>            _rho;
	TTreeReaderAnyValue<ULong64_t>          _HLTPho;
        TTreeReaderAnyValue<ULong64_t>          _HLTEleMuX;
	TTreeReaderAnyValue<UShort_t>           _beamHaloSummary;
	TTreeReaderAnyValue<UShort_t>           _metFilters;

	TTreeReaderAnyValue<UChar_t>            _puTrue;
	TTreeReaderAnyValue<Float_t>            _genWeight;
	TTreeReaderAnyValue<UShort_t>		_nMC;
        TTreeReaderVectorValue<Int_t>           _mcMomPID;
	TTreeReaderVectorValue<Int_t>           _mcPID;
	TTreeReaderVectorValue<Float_t>		_mcPt;
	TTreeReaderVectorValue<Float_t>         _mcEta;
	TTreeReaderVectorValue<Float_t>         _mcPhi;
        TTreeReaderVectorValue<Float_t>         _mcE;
        TTreeReaderVectorValue<Float_t>         _mcEt;
	TTreeReaderVectorValue<UShort_t>        _mcStatusFlag;
        TTreeReaderVectorValue<Char_t>          _mcPromptStatusType;
	TTreeReaderVectorValue<Short_t>         _mcStatus;
	TTreeReaderVectorValue<Short_t>     	_mcIndex;
        //TTreeReaderVectorValue<Short_t>         _pho_gen_index;

	TTreeReaderAnyValue<UShort_t>           _nPho;
	TTreeReaderVectorValue<Float_t>         _phoCalibEt;
        TTreeReaderVectorValue<Float_t>         _phoCalibE;
	TTreeReaderVectorValue<Float_t>         _phoEt;
	TTreeReaderVectorValue<Float_t>         _phoEta;
	TTreeReaderVectorValue<Float_t>         _phoPhi;	
	TTreeReaderVectorValue<Float_t>         _phoSeedTime;
	TTreeReaderVectorValue<UChar_t>         _phoFiducialRegion;

	TTreeReaderVectorValue<UChar_t>         _phoQualityBits;
        TTreeReaderVectorValue<Float_t>         _phoR9;
	TTreeReaderVectorValue<Float_t>         _phoR9Full5x5;
	TTreeReaderVectorValue<Float_t>         _phoSigmaIEtaIEtaFull5x5;
	TTreeReaderVectorValue<Float_t>         _phoSigmaIEtaIPhiFull5x5;
	TTreeReaderVectorValue<Float_t>         _phoSigmaIPhiIPhiFull5x5;	
	TTreeReaderVectorValue<Float_t>         _phoE2x2Full5x5;
	TTreeReaderVectorValue<Float_t>		_phoE5x5Full5x5;

	TTreeReaderVectorValue<Float_t>		_phoMaxEnergyXtal;
	TTreeReaderVectorValue<Float_t>		_phoE2ndFull5x5;
	TTreeReaderVectorValue<Float_t>		_phoE1x3Full5x5;
	TTreeReaderVectorValue<Float_t>		_phoE1x5Full5x5;
	TTreeReaderVectorValue<Float_t>		_phoE2x5Full5x5;

	TTreeReaderVectorValue<Float_t>		_phoPFClusEcalIso;
	TTreeReaderVectorValue<Float_t>		_phoPFClusHcalIso;
	TTreeReaderVectorValue<Float_t>		_phoTrkSumPtSolidConeDR04;
	TTreeReaderVectorValue<Float_t>		_phoTrkSumPtHollowConeDR04;
	TTreeReaderVectorValue<Float_t>		_phoTrkSumPtSolidConeDR03;
	TTreeReaderVectorValue<Float_t>		_phoTrkSumPtHollowConeDR03;
	TTreeReaderVectorValue<Float_t>		_phoECALIso;
	TTreeReaderVectorValue<Float_t>		_phoHCALIso;

	TTreeReaderVectorValue<Float_t>         _phoHoverE;
	TTreeReaderVectorValue<Float_t>         _phoPFChIso;
	TTreeReaderVectorValue<Float_t>         _phoPFPhoIso;
	TTreeReaderVectorValue<Float_t>         _phoPFNeuIso;
	TTreeReaderVectorValue<Float_t>         _phoPFChWorstIso;
	TTreeReaderVectorValue<Float_t>         _phoIDMVA;
	TTreeReaderVectorValue<UChar_t>         _phoIDbit;
	TTreeReaderVectorValue<Float_t>         _phoMIPTotEnergy;
	
	TTreeReaderVectorValue<Short_t>         _phoDirectEcalSCindex;
	TTreeReaderVectorValue<Float_t>         _ecalSCeta;
	TTreeReaderVectorValue<Float_t>         _ecalSCphi;
	TTreeReaderVectorValue<Float_t>         _ecalSCEn;
	TTreeReaderVectorValue<Float_t>         _ecalSCRawEn;
	TTreeReaderVectorValue<Float_t>         _ecalSCetaWidth;
	TTreeReaderVectorValue<Float_t>         _ecalSCphiWidth;

        TTreeReaderAnyValue<Float_t>            _genMET;
	TTreeReaderAnyValue<Float_t>            _pfMET;
	TTreeReaderAnyValue<Float_t>            _pfMETPhi;
	TTreeReaderAnyValue<Float_t>            _pfMET_metSig;
	TTreeReaderAnyValue<Float_t>            _pfMET_EtSig;

	TTreeReaderAnyValue<UShort_t>           _nMu;
	TTreeReaderVectorValue<Float_t>         _muPt;
	TTreeReaderVectorValue<Int_t>           _muIDbit;

	TTreeReaderAnyValue<UShort_t>           _nAK4CHSJet;
	TTreeReaderVectorValue<Float_t>         _AK4CHSJet_Pt;
	TTreeReaderVectorValue<Float_t>         _AK4CHSJet_Eta;
	TTreeReaderVectorValue<Float_t>         _AK4CHSJet_Phi;
	TTreeReaderVectorValue<Char_t>          _AK4CHSJet_PUFullID;
	TTreeReaderVectorValue<Char_t>          _AK4CHSJet_ID;


	TTreeReaderAnyValue<UChar_t>		_ntrgObjPho;
	TTreeReaderVectorValue<UChar_t>		_trgObjPhoBits;
	TTreeReaderVectorValue<Float_t>		_trgObjPhoPt;
	TTreeReaderVectorValue<Float_t>		_trgObjPhoEta;
	TTreeReaderVectorValue<Float_t>		_trgObjPhoPhi;

        TTreeReaderAnyValue<UShort_t> 		_nTau;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////// Output variables //////////////////////////////////////////////////////////

	Int_t		 run_;
	Long64_t	 event_;
	UShort_t	 lumis_;	
	Float_t	 	 rho_;
	UChar_t	 	 nVtx_;

	UShort_t		metFilters_;
	UShort_t		beamHaloSummary_;

	Float_t         genWeight_ = 1.;
	Float_t         puWeight_ = 1.;
	UChar_t		puTrue_;
	UShort_t	genStatusFlag_;
	Short_t		genStatus_;
	Float_t 	deltaRgenPho_;
	Float_t 	relDeltaPtGenPho_;
	Float_t 	deltaRPt_;
	Int_t		genPDGid_ = 0;

	Float_t 		deltaRtrg_;
	Float_t 		deltaPttrg_;

	Float_t 		phoPt_;
	Float_t 		phoEta_;
	Float_t 		phoPhi_;
        Float_t                 phoE_;
        Float_t                 phoEt_;
	Float_t 		phoSeedTime_;

	UChar_t			phoGenBits_;
	UChar_t			phoQualityBits_;
	Float_t 		phoR9Full5x5_;
	Float_t 		phoS4Full5x5_;

	Float_t 		phoEmaxOESCrFull5x5_;
	Float_t 		phoE2ndOESCrFull5x5_;
	Float_t 		phoE2ndOEmaxFull5x5_;
	Float_t 		phoE1x3OESCrFull5x5_;
	Float_t 		phoE2x5OESCrFull5x5_;
	Float_t 		phoE5x5OESCrFull5x5_;

	Float_t 		phoEmaxOE3x3Full5x5_;
	Float_t 		phoE2ndOE3x3Full5x5_;
	Float_t 		pho2x2OE3x3Full5x5_;
	Float_t 		phoSieieOSipipFull5x5_;
	Float_t			phoEtaWOPhiWFull5x5_;

	Float_t 		phoSigmaIEtaIEta_;
	Float_t 		phoSigmaIPhiIPhi_;
	Float_t 		phoSigmaIEtaIPhi_;

	Float_t 		phoE2x2Full5x5_;
	Float_t 		phoE3x3Full5x5_;
	Float_t 		phoE5x5Full5x5_;
	Float_t			phoMaxEnergyXtal_;
	Float_t			phoE2ndFull5x5_;
	Float_t			phoE1x3Full5x5_;
	Float_t			phoE1x5Full5x5_;
	Float_t			phoE2x5Full5x5_;

	Float_t			phoPFClusEcalIso_;
	Float_t			phoPFClusHcalIso_;
	Float_t			phoTrkSumPtSolidConeDR04_;
	Float_t			phoTrkSumPtHollowConeDR04_;
	Float_t			phoTrkSumPtSolidConeDR03_;
	Float_t			phoTrkSumPtHollowConeDR03_;
	Float_t			phoECALIso_;
	Float_t			phoHCALIso_;
	Float_t 		phoPFECALClusIsoCorr_;
	Float_t 		phoPFHCALClusIsoCorr_;
	Float_t 		phoTkrIsoCorr_;	

	Float_t 		phoHoverE_;
	Float_t 		phoPFChIsoRaw_;
	Float_t 		phoPFPhoIsoRaw_;
	Float_t 		phoPFNeuIsoRaw_;
	Float_t 		phoPFChWorstIsoRaw_;
	Float_t 		phoIDMVA_;
	UChar_t                 phoIDbit_;
	Float_t 		phoBDTpred_;
	UChar_t 		phoPFClusIDbits_;

	Bool_t 			pass95_ = 0;

	Float_t 		phoMIP_;	

	Float_t 		phoSCet_;
	Float_t 		phoSCrawet_;
	Float_t 		phoSCeta_;
	Float_t 		phoSCphi_;
	Float_t 		phoSCEn_;
	Float_t 		phoSCRawEn_;
	Float_t 		phoEtaWidth_;
	Float_t 		phoPhiWidth_;

	UChar_t 		lepVeto_;

        Float_t                 genMET_;
	Float_t 		met_;
	Float_t			metPhi_;
	Float_t 		metSig_;
	Float_t 		EtSig_;

	Float_t 		mT_;
	Float_t 		deltaPhiMetPho_;
	Float_t 		phoPtOverMet_;


	UChar_t 		nJet30_;
	Float_t 		Ht30_;
	Float_t 		minDeltaPhiMetJet30_;
	Float_t 		minDeltaRPhoJet30_;
        Float_t                 mt_;
        Int_t                   mcTrue_ = 0;
        Int_t                   nGenPho_ = 0;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////Buffer variables///////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////// XGBoost ////////////////////////////////////////////////////////////
	DMatrixHandle 		dTest;
	BoosterHandle 		phoBDT_h;
	Bool_t 			predictBDT = 0;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////Categories////////////////////////////////////////////////////////////
	Bool_t          	initEventTypes();
	void            	initEventType(eventType & evType, std::string typeName, std::string typeTitle);
	void            	fillEventType(eventType & evType);
	void            	registerCutFlow(eventType & evType);
	void            	registerAllCutFlow();
	eventType		fullECAL;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
atgcStudies::atgcStudies(std::string FILELIST, std::string OUTFILE, Float_t XSECTION, std::string MCPILEUPHIST, std::string DATAPILEUPHIST, //std::string HoE_EFFECTIVE_AREAS, 
std::string PFECALCLUS_EFFECTIVE_AREAS, std::string PFHCALCLUS_EFFECTIVE_AREAS, std::string TKRISO_EFFECTIVE_AREAS,
	std::string PFECALCLUS_PTSCALING,	std::string PFHCALCLUS_PTSCALING,
	std::string BDT_PATH){

	std::cout<<"*************************************************************************************************************************************************"<<std::endl<<
	getCurrentTime()<<std::endl<<
	"Running atgcStudies"<<std::endl<<
	"\n\nInput parameters:"<<std::endl<<
	"\t\tFile list = "<<FILELIST<<std::endl<<
	"\t\tOutput file = "<<OUTFILE<<std::endl<<
	"\t\tCross section = "<<XSECTION<<std::endl<<
	"\t\tMC pileup histogram = "<<MCPILEUPHIST<<std::endl<<
	"\t\tData pileup histogram = "<<DATAPILEUPHIST<<std::endl<<
	"\t\tECAL PFCLuster Isolation Effective Areas = "<<PFECALCLUS_EFFECTIVE_AREAS<<std::endl<<
	"\t\tHCAL PFCLuster Isolation Effective Areas = "<<PFHCALCLUS_EFFECTIVE_AREAS<<std::endl<<
	"\t\tTracker Isolation Effective Areas = "<<TKRISO_EFFECTIVE_AREAS<<std::endl<<
	"\t\tECAL PFCLuster Isolation pT Scaling = "<<PFECALCLUS_PTSCALING<<std::endl<<
	"\t\tHCAL PFCLuster Isolation pT Scaling = "<<PFHCALCLUS_PTSCALING<<std::endl<<
	"\t\tBDT model file = "<<BDT_PATH<<std::endl;

	xSec = XSECTION;
	if(XSECTION > 0.) isMC = true;
	if(isMC && file_exists(MCPILEUPHIST) && file_exists(DATAPILEUPHIST)) doPUreweight = true;
        std::cout << "####################################" << std::endl;
        std::cout << file_exists(MCPILEUPHIST) << "   " << file_exists(DATAPILEUPHIST) << std::endl;

	std::cout<<"\t\tSample is simulation = "<<std::boolalpha<<isMC<<std::endl<<
	"\t\tDo pileup reweight = "<<doPUreweight<<"\n\n"<<std::endl;

	// doPUreweight = 0;
	if(doPUreweight) {
		std::cout<<"Pileup reweighting:"<<std::endl;
		puReweighter.init(MCPILEUPHIST, DATAPILEUPHIST, "hPUTruew", "pileup");
	}

	if(file_exists(PFECALCLUS_EFFECTIVE_AREAS)) {
		std::cout<<"PF ECAL Cluster effective areas:"<<std::endl;
		PFECALClusEffAreas.init(PFECALCLUS_EFFECTIVE_AREAS, 1, ",", 0);
	}

	if(file_exists(PFHCALCLUS_EFFECTIVE_AREAS)) {
		std::cout<<"PF HCAL Cluster effective areas:"<<std::endl;
		PFHCALClusEffAreas.init(PFHCALCLUS_EFFECTIVE_AREAS, 1, ",", 0);
	}
	
	if(file_exists(TKRISO_EFFECTIVE_AREAS)){
		std::cout<<"Tracker isolation effective areas:"<<std::endl;
		TkrEffAreas.init(TKRISO_EFFECTIVE_AREAS, 1, ",", 0);
	}

	if(file_exists(PFECALCLUS_PTSCALING)){
		std::cout<<"ECAL pT scaling:"<<std::endl;
		PFECALClusPtScaling.init(PFECALCLUS_PTSCALING, 0, 1, ",", 0);
	}

	if(file_exists(PFHCALCLUS_PTSCALING)){
		std::cout<<"HCAL pT scaling:"<<std::endl;
		PFHCALClusPtScaling.init(PFHCALCLUS_PTSCALING, 1, 1, ",", 0);
	}

	if(file_exists(BDT_PATH)){
		std::cout<<"\nLoading BDT model from "<<BDT_PATH <<std::endl;
		XGBoosterCreate(NULL, 0, &phoBDT_h);
		XGBoosterSetParam(phoBDT_h, "seed", "0");
		Int_t mLdSuccess = XGBoosterLoadModel(phoBDT_h, BDT_PATH.c_str());
		if(mLdSuccess == 0) predictBDT=1;
		else{
			std::cout<<"Failed to load BDT model!"<<std::endl;
		}
	}

	std::cout<<"\nCreating TChain... " <<std::endl;
	initNtuples(FILELIST);

	outFile = new TFile(OUTFILE.c_str(), "RECREATE");

	initEventTypes();

	analyze();

	outFile->Write();
	outFile->Close();

	closeTChain(inputTree);

	std::cout<<"\n\nOutput written to file\t"<<OUTFILE <<std::endl<<"Complete!"<<std::endl<<getCurrentTime()<<std::endl;
};

Bool_t atgcStudies::photonPassingID(UShort_t phoIndex){

        Float_t rho                    = _rho;

        Short_t phoSCindex      = _phoDirectEcalSCindex[phoIndex];
        Float_t phoSCeta        = _ecalSCeta[phoSCindex];
        Float_t phoAbsSCEta     = std::abs(phoSCeta);

        Float_t phoPt                  = _phoCalibEt[phoIndex];

        Float_t phoR9Full5x5           = _phoR9Full5x5[phoIndex];
        Float_t phoS4Full5x5           = _phoE2x2Full5x5[phoIndex]/_ecalSCRawEn[phoSCindex];
        Float_t phoEmaxOESCrFull5x5    = _phoMaxEnergyXtal[phoIndex]/_ecalSCRawEn[phoSCindex];
        Float_t phoE2ndOESCrFull5x5    = _phoE2ndFull5x5[phoIndex]/_ecalSCRawEn[phoSCindex];
        Float_t phoE2ndOEmaxFull5x5    = _phoE2ndFull5x5[phoIndex]/_phoMaxEnergyXtal[phoIndex];
        Float_t phoE1x3OESCrFull5x5    = _phoE1x3Full5x5[phoIndex]/_ecalSCRawEn[phoSCindex];
        Float_t phoE2x5OESCrFull5x5    = _phoE2x5Full5x5[phoIndex]/_ecalSCRawEn[phoSCindex];
        Float_t phoE5x5OESCrFull5x5    = _phoE5x5Full5x5[phoIndex]/_ecalSCRawEn[phoSCindex];

        Float_t phoSigmaIEtaIEta       = _phoSigmaIEtaIEtaFull5x5[phoIndex];
        Float_t phoSigmaIEtaIPhi       = _phoSigmaIEtaIPhiFull5x5[phoIndex];
        Float_t phoSigmaIPhiIPhi       = _phoSigmaIPhiIPhiFull5x5[phoIndex];

        Float_t phoE2x2Full5x5         = _phoE2x2Full5x5[phoIndex];
        Float_t phoE3x3Full5x5         = phoR9Full5x5 * _ecalSCRawEn[phoSCindex];
        Float_t phoE5x5Full5x5         = _phoE5x5Full5x5[phoIndex];
        Float_t phoMaxEnergyXtal       = _phoMaxEnergyXtal[phoIndex];
        Float_t phoE2ndFull5x5         = _phoE2ndFull5x5[phoIndex];
        Float_t phoE1x3Full5x5         = _phoE1x3Full5x5[phoIndex];
        Float_t phoE1x5Full5x5         = _phoE1x5Full5x5[phoIndex];
        Float_t phoE2x5Full5x5         = _phoE2x2Full5x5[phoIndex];

        Float_t phoEmaxOE3x3Full5x5     = phoMaxEnergyXtal/phoE3x3Full5x5;
        Float_t phoE2ndOE3x3Full5x5     = phoE2ndFull5x5/phoE3x3Full5x5;
        Float_t pho2x2OE3x3Full5x5      = phoE2x2Full5x5/phoE3x3Full5x5;
        Float_t phoSieieOSipipFull5x5   = phoSigmaIEtaIEta/phoSigmaIPhiIPhi;
        Float_t phoEtaWidth             = _ecalSCetaWidth[phoSCindex];
        Float_t phoPhiWidth             = _ecalSCphiWidth[phoSCindex];

        Float_t phoPFClusEcalIso               = _phoPFClusEcalIso[phoIndex];
        Float_t phoPFClusHcalIso               = _phoPFClusHcalIso[phoIndex];
        Float_t phoTrkSumPtHollowConeDR03      = _phoTrkSumPtHollowConeDR03[phoIndex];
        Float_t phoECALIso                     = _phoECALIso[phoIndex];
        Float_t phoHCALIso                     = _phoHCALIso[phoIndex];

        Float_t phoPFECALClusIsoCorr           =       phoPFClusEcalIso - rho * PFECALClusEffAreas.getEffectiveArea(phoAbsSCEta) - PFECALClusPtScaling.getPtScaling(phoAbsSCEta, phoPt);
        Float_t phoPFHCALClusIsoCorr           =       phoPFClusHcalIso - rho * PFHCALClusEffAreas.getEffectiveArea(phoAbsSCEta) - PFHCALClusPtScaling.getPtScaling(phoAbsSCEta, phoPt);
        Float_t phoTkrIsoCorr                  =       phoTrkSumPtHollowConeDR03 - rho * TkrEffAreas.getEffectiveArea(phoAbsSCEta);

        Float_t phoHoverE              = _phoHoverE[phoIndex];
        Float_t phoEtaWOPhiWFull5x5    = phoEtaWidth/phoPhiWidth;

        if(predictBDT){
        std::vector<Float_t> feats{pho2x2OE3x3Full5x5, phoE1x3OESCrFull5x5, phoE2ndOESCrFull5x5, phoE2x5OESCrFull5x5, phoEmaxOESCrFull5x5, phoEtaWOPhiWFull5x5, phoEtaWidth, phoPhiWidth, phoR9Full5x5, phoS4Full5x5, phoSieieOSipipFull5x5, phoSigmaIEtaIEta, phoSigmaIEtaIPhi, phoSigmaIPhiIPhi};

        XGDMatrixCreateFromMat((float*)feats.data(), 1, feats.size(), -9999999999, &dTest);
        bst_ulong out_len;
        const float *prediction;
        XGBoosterPredict(phoBDT_h, dTest, 0, 0, 0, &out_len, &prediction);
        assert(out_len == 1);
        XGDMatrixFree(dTest);
        phoBDTpred_ = prediction[0];

        pass95_ = (phoBDTpred_ >= 4.9013550283797025e-01) && (phoHoverE < 4.2132449148158196e-02) && (phoPFECALClusIsoCorr < 5.1143879813442599e+00) && (phoPFHCALClusIsoCorr < 1.7020221456277284e+01) && (phoTkrIsoCorr < 3.5001887785563284e+00);

        }
        if (pass95_) {
           setBit(phoPFClusIDbits_, 3, 1);
           return 1;
        } else {
           return 0;
        }
};
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Float_t atgcStudies::deltaR(Float_t eta1, Float_t phi1, Float_t eta2, Float_t phi2) {

              Float_t deltaEta = 0;
              Float_t deltaPhi = 0;

              deltaPhi = fabs(phi1 - phi2);
              if (deltaPhi > acos(-1)) {
              deltaPhi = 2*acos(-1) - deltaPhi;
              }
              deltaEta = eta1 - eta2;
              Float_t deltaR_ = sqrt((deltaPhi*deltaPhi) + (deltaEta*deltaEta));

              return deltaR_;

};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Short_t	atgcStudies::photonIsPrompt(Short_t _phoIndex, Float_t _deltaRmax) { //Float_t _relDeltaPtMin, Float_t _relDeltaPtMax){

	Short_t matchedPromptGenPho = -999;
	Float_t minDeltaR = 999.;

	for(UShort_t iGenP=0; iGenP < _nMC; iGenP++){
		if(_mcPID[iGenP] != 22) continue;

		UShort_t iGenPStFl = _mcStatusFlag[iGenP];
		if(!getBit(iGenPStFl,1)) continue;

		Float_t dRiGenPho = deltaR(_phoEta[_phoIndex], _phoPhi[_phoIndex], _mcEta[iGenP], _mcPhi[iGenP]);
		if(dRiGenPho > _deltaRmax) continue;

		//Float_t relDeltaPtiGenPho = std::abs(_mcPt[iGenP] - _phoCalibEt[_phoIndex])/_mcPt[iGenP];
		//if(relDeltaPtiGenPho > _relDeltaPtMax) continue;
		//if(relDeltaPtiGenPho < _relDeltaPtMin) continue;

		if(dRiGenPho < minDeltaR){
			minDeltaR = dRiGenPho;
			matchedPromptGenPho = iGenP;
		}
	}

	return matchedPromptGenPho;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Bool_t atgcStudies::initEventTypes(){

	pileupPreweight.SetDirectory(outFile->GetDirectory(""));
	pileupPostweight.SetDirectory(outFile->GetDirectory(""));
	rhoPreweight.SetDirectory(outFile->GetDirectory(""));
	rhoPostweight.SetDirectory(outFile->GetDirectory(""));
	nvtxPreweight.SetDirectory(outFile->GetDirectory(""));
 	nvtxPostweight.SetDirectory(outFile->GetDirectory(""));
        
	initEventType(fullECAL, "atgc", "Full ECAL");
	return 1;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void atgcStudies::analyze(){
	std::cout<<"--------------------------------------------------------------------------------------------------"<<std::endl<<
	getCurrentTime()<<std::endl<<
	"Analyzing events.."<<std::endl;
	ULong64_t current_entry =0;
	while(inputTTreeReader.Next()){
		if(current_entry % REPORT_EVERY == 0){
			std::cout<<"\t"<< getCurrentTime()<<"\tAnalyzing entry\t"<<current_entry<<
			",\t\tevent\t"<<(_event)<<"\t\tFile " <<inputTree->GetCurrentFile()->GetName() <<std::endl;
		}

                bool ispassHLT = _HLTPho>>9&1;
                if (!ispassHLT) continue;
                //if (_pfMET < 150) continue; //high MET
                if (_genMET < 40) continue;
                  
		rhoPreweight.Fill(_rho, genWeight_);
		nvtxPreweight.Fill(_nVtx, genWeight_);

                Float_t weight = 1;
		if(isMC) {
                        std::string Map = "/local/cms/user/gsorrent/antgc_analysis/atgcStudies/ntuples_2017_METbugfix.csv";
                        std::string filePath = inputTree->GetCurrentFile()->GetName();
                        std::string sample = findAndReplaceAll(getFileName(filePath), ".root", "");
                        std::string sampleName = sample.substr(0, sample.find("_"));
 
                        Double_t xSection = std::stof(vLookup(sampleName, Map, 0, 2));
                        Double_t sumGenWeight = std::stof(vLookup(sampleName, Map, 0, 7));
                        //std::cout << "xsec: " << xSection << std::endl;
                        //std::cout << "sumGenWeight: " << sumGenWeight << std::endl;

                        mcTrue_ = 1; 
                        Double_t lumi = 41.56;

                        Float_t norm = xSection * 1000 * lumi/sumGenWeight;

			if(doPUreweight){
				genWeight_ = _genWeight;
				puWeight_ = puReweighter.weight(_puTrue) * genWeight_ * norm; //puWeight + MC norm + genWeight
				Float_t genPUweight = genWeight_ * puWeight_;
                                weight = genPUweight;
				pileupPostweight.Fill(_puTrue, genPUweight);
				rhoPostweight.Fill(_rho, genPUweight);
				nvtxPostweight.Fill(_nVtx, genPUweight);
			} else {
				genWeight_ = _genWeight;
				puWeight_ = 1.;    
			}
			pileupPreweight.Fill(_puTrue, genWeight_);
		}

               nGenPho_ = 0;
               Short_t highetsPtIndex  = -999;
               Float_t highestPhoPt    = -999.;

               for(UShort_t iGenP=0; iGenP<_nMC; iGenP++){
                  if(_mcPID[iGenP] != 22) continue;
                  UShort_t iGenPStFl = _mcStatusFlag[iGenP];
                  if(!getBit(iGenPStFl,1)) continue;
                  if(_mcPt[iGenP] < 225.) continue;
                  nGenPho_++;

                  if(_mcPt[iGenP] > highestPhoPt){
                      highestPhoPt = _mcPt[iGenP];
                      highetsPtIndex = iGenP;
                  }
               }
               

               bool foundPho = false; 
               UShort_t iPhoton = 999;
               if (nGenPho_ > 0) {
                 foundPho = true;
                 iPhoton = highetsPtIndex;
               }

               if (foundPho) {
                  run_                    = _run;
                  event_                  = _event;
                  lumis_                  = _lumis;
                  rho_                    = _rho;
                  nVtx_                   = _nVtx;

                  //metFilters_             = _metFilters;
                  //beamHaloSummary_        = _beamHaloSummary;
                  genMET_                    = _genMET;

                  phoPt_                  = _mcPt[iPhoton];
                  phoEta_                 = _mcEta[iPhoton];
                  phoPhi_                 = _mcPhi[iPhoton];
                  phoE_                   = _mcE[iPhoton];
                  phoEt_                  = _mcEt[iPhoton];

                  //Fill tree
                  fillEventType(fullECAL);  
               }

	       current_entry++; 
               //cout << "Weight: " << weight << endl;
	};
	std::cout<<"Done analyzing!"<<std::endl<<
	getCurrentTime()<<std::endl<<
	"--------------------------------------------------------------------------------------------------"<<std::endl;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Bool_t atgcStudies::initNtuples(std::string FILELIST){

	inputTree = openTChain(FILELIST, "ggNtuplizer/EventTree");
	inputTTreeReader.SetTree(inputTree);

	std::cout<<"**************************************************************************************************************************************************************"<<std::endl<<
	"Initializing branches in input ntuples..."<<std::endl;

	_run.set(inputTTreeReader, "run");
	_event.set(inputTTreeReader, "event");
	_lumis.set(inputTTreeReader, "lumis");
	_beamHaloSummary.set(inputTTreeReader, "beamHaloSummary");
	_metFilters.set(inputTTreeReader, "metFilters");
	_nVtx.set(inputTTreeReader, "nVtx");
	_rho.set(inputTTreeReader, "rho");
	_HLTPho.set(inputTTreeReader, "HLTPho");
        _HLTEleMuX.set(inputTTreeReader, "HLTEleMuX");
	if(isMC){
		_puTrue.set(inputTTreeReader, "puTrue");
		_genWeight.set(inputTTreeReader, "genWeight");
		_nMC.set(inputTTreeReader, "nMC");
		_mcPID.set(inputTTreeReader, "mcPID");
                _mcMomPID.set(inputTTreeReader, "mcMomPID");
		_mcPt.set(inputTTreeReader, "mcPt");
		_mcEta.set(inputTTreeReader, "mcEta");
		_mcPhi.set(inputTTreeReader, "mcPhi");
                _mcE.set(inputTTreeReader, "mcE");
                _mcEt.set(inputTTreeReader, "mcEt");
		_mcStatusFlag.set(inputTTreeReader, "mcStatusFlag");
                _mcPromptStatusType.set(inputTTreeReader, "mcPromptStatusType");
		_mcStatus.set(inputTTreeReader, "mcStatus");
		_mcIndex.set(inputTTreeReader, "mcIndex");
                _genMET.set(inputTTreeReader, "genMET");
                //_pho_gen_index.set(inputTTreeReader, "pho_gen_index");
	}

	_nPho.set(inputTTreeReader, "nPho");
	_phoEt.set(inputTTreeReader, "phoEt");
	_phoCalibEt.set(inputTTreeReader, "phoCalibEt");
        _phoCalibE.set(inputTTreeReader, "phoCalibE");
	_phoEta.set(inputTTreeReader, "phoEta");
	_phoPhi.set(inputTTreeReader, "phoPhi");
	_phoSeedTime.set(inputTTreeReader, "phoSeedTime");
	_phoFiducialRegion.set(inputTTreeReader, "phoFiducialRegion");

	_phoQualityBits.set(inputTTreeReader, "phoQualityBits");
        _phoR9.set(inputTTreeReader, "phoR9");
	_phoR9Full5x5.set(inputTTreeReader, "phoR9Full5x5");
	_phoSigmaIEtaIEtaFull5x5.set(inputTTreeReader, "phoSigmaIEtaIEtaFull5x5");
	_phoSigmaIEtaIPhiFull5x5.set(inputTTreeReader, "phoSigmaIEtaIPhiFull5x5");
	_phoSigmaIPhiIPhiFull5x5.set(inputTTreeReader, "phoSigmaIPhiIPhiFull5x5");
	_phoE2x2Full5x5.set(inputTTreeReader, "phoE2x2Full5x5");
	_phoE5x5Full5x5.set(inputTTreeReader, "phoE5x5Full5x5");

	_phoMaxEnergyXtal.set(inputTTreeReader, "phoMaxEnergyXtal");
	_phoE2ndFull5x5.set(inputTTreeReader, "phoE2ndFull5x5");
	_phoE1x3Full5x5.set(inputTTreeReader, "phoE1x3Full5x5");
	_phoE1x5Full5x5.set(inputTTreeReader, "phoE1x5Full5x5");
	_phoE2x5Full5x5.set(inputTTreeReader, "phoE2x5Full5x5");

	_phoPFClusEcalIso.set(inputTTreeReader, "phoPFClusEcalIso");
	_phoPFClusHcalIso.set(inputTTreeReader, "phoPFClusHcalIso");
	_phoTrkSumPtSolidConeDR04.set(inputTTreeReader, "phoTrkSumPtSolidConeDR04");
	_phoTrkSumPtHollowConeDR04.set(inputTTreeReader, "phoTrkSumPtHollowConeDR04");
	_phoTrkSumPtSolidConeDR03.set(inputTTreeReader, "phoTrkSumPtSolidConeDR03");
	_phoTrkSumPtHollowConeDR03.set(inputTTreeReader, "phoTrkSumPtHollowConeDR03");
	_phoECALIso.set(inputTTreeReader, "phoECALIso");
	_phoHCALIso.set(inputTTreeReader, "phoHCALIso");
	
	_phoHoverE.set(inputTTreeReader, "phoHoverE");
	_phoPFChIso.set(inputTTreeReader, "phoPFChIso");
	_phoPFPhoIso.set(inputTTreeReader, "phoPFPhoIso");
	_phoPFNeuIso.set(inputTTreeReader, "phoPFNeuIso");
	_phoPFChWorstIso.set(inputTTreeReader, "phoPFChWorstIso");
	_phoIDMVA.set(inputTTreeReader, "phoIDMVA");
	_phoIDbit.set(inputTTreeReader, "phoIDbit");
	_phoMIPTotEnergy.set(inputTTreeReader, "phoMIPTotEnergy");
	
	_phoDirectEcalSCindex.set(inputTTreeReader, "phoDirectEcalSCindex");
	_ecalSCeta.set(inputTTreeReader, "ecalSCeta");
	_ecalSCphi.set(inputTTreeReader, "ecalSCphi");
	_ecalSCEn.set(inputTTreeReader, "ecalSCEn");
	_ecalSCRawEn.set(inputTTreeReader, "ecalSCRawEn");
	_ecalSCetaWidth.set(inputTTreeReader, "ecalSCetaWidth");
	_ecalSCphiWidth.set(inputTTreeReader, "ecalSCphiWidth");

	_metFilters.set(inputTTreeReader, "metFilters");
	_pfMET.set(inputTTreeReader, "pfMET");
	_pfMETPhi.set(inputTTreeReader, "pfMETPhi");
	_pfMET_metSig.set(inputTTreeReader, "pfMET_metSig");
	_pfMET_EtSig.set(inputTTreeReader, "pfMET_EtSig");

	_nMu.set(inputTTreeReader, "nMu");
	_muPt.set(inputTTreeReader, "muPt");
	_muIDbit.set(inputTTreeReader, "muIDbit");

	_nAK4CHSJet.set(inputTTreeReader, "nAK4CHSJet");
	_AK4CHSJet_Pt.set(inputTTreeReader, "AK4CHSJet_Pt");
	_AK4CHSJet_Eta.set(inputTTreeReader, "AK4CHSJet_Eta");
	_AK4CHSJet_Phi.set(inputTTreeReader, "AK4CHSJet_Phi");
	_AK4CHSJet_PUFullID.set(inputTTreeReader, "AK4CHSJet_PUFullID");
	_AK4CHSJet_ID.set(inputTTreeReader, "AK4CHSJet_ID");

	_ntrgObjPho.set(inputTTreeReader, "ntrgObjPho");
	_trgObjPhoBits.set(inputTTreeReader, "trgObjPhoBits");
	_trgObjPhoPt.set(inputTTreeReader, "trgObjPhoPt");
	_trgObjPhoEta.set(inputTTreeReader, "trgObjPhoEta");
	_trgObjPhoPhi.set(inputTTreeReader, "trgObjPhoPhi");

        _nTau.set(inputTTreeReader,"nTau");
	std::cout<<"Branches initialized!"<<std::endl<<
	"**************************************************************************************************************************************************************"<<std::endl;
	return 1;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void atgcStudies::initEventType(eventType & evType, std::string typeName, std::string typeTitle){

	mkTFileDir(outFile, typeName);

	evType.cutFlowCount = new TH1F((typeName+"cutFlowCount").c_str(), "Cut Flow (Unweighted)", CUTFLOWSTEPS, 0., (Float_t)CUTFLOWSTEPS);
	evType.cutFlowGenWeight = new TH1F((typeName+"cutFlowGenWeight").c_str(), "Cut Flow (Gen Weighted)", CUTFLOWSTEPS, 0., (Float_t)CUTFLOWSTEPS);
	evType.cutFlowCount->SetDirectory(outFile->GetDirectory(typeName.c_str()));
	evType.cutFlowGenWeight->SetDirectory(outFile->GetDirectory(typeName.c_str()));

	//////// initialize tree
	evType.tree = new TTree("fitter_tree", typeTitle.c_str());
	evType.tree->SetDirectory(outFile->GetDirectory(typeName.c_str()));

        //changing branch names to match egamma names
	evType.tree->Branch("run", &run_);
	evType.tree->Branch("event", &event_);
	evType.tree->Branch("lumi", &lumis_);
	evType.tree->Branch("event_rho", &rho_);
	evType.tree->Branch("event_nPV", &nVtx_);
	evType.tree->Branch("truePU", &puTrue_);
	evType.tree->Branch("metFilters", &metFilters_); //not in eg tree
	evType.tree->Branch("beamHaloSummary", &beamHaloSummary_); //not in eg tree

	evType.tree->Branch("puWeight", &puWeight_);
	evType.tree->Branch("genWeight", &genWeight_);
	evType.tree->Branch("genStatusFlag", &genStatusFlag_);
	evType.tree->Branch("genStatus", &genStatus_);
	evType.tree->Branch("deltaRgenPho", &deltaRgenPho_);
	evType.tree->Branch("relDeltaPtGenPho", &relDeltaPtGenPho_);
	evType.tree->Branch("deltaRPt", &deltaRPt_);
	evType.tree->Branch("genPDGid", &genPDGid_);	

	evType.tree->Branch("deltaRtrg", &deltaRtrg_);
	evType.tree->Branch("deltaPttrg", &deltaPttrg_);

        //useful branches
	evType.tree->Branch("ph_pt", &phoPt_);
	evType.tree->Branch("ph_eta", &phoEta_);
	evType.tree->Branch("ph_phi", &phoPhi_);
        evType.tree->Branch("ph_e", &phoE_);
        evType.tree->Branch("ph_et", &phoEt_);

	evType.tree->Branch("ph_seedtime", &phoSeedTime_);
	evType.tree->Branch("ph_GenBits", &phoGenBits_);
	evType.tree->Branch("ph_QualityBits", &phoQualityBits_);

	evType.tree->Branch("ph_R9Full5x5", &phoR9Full5x5_);
	evType.tree->Branch("ph_S4Full5x5", &phoS4Full5x5_);
	evType.tree->Branch("ph_EmaxOESCrFull5x5", &phoEmaxOESCrFull5x5_);
	evType.tree->Branch("ph_E2ndOESCrFull5x5", &phoE2ndOESCrFull5x5_);
	evType.tree->Branch("ph_E2ndOEmaxFull5x5", &phoE2ndOEmaxFull5x5_);	
	evType.tree->Branch("ph_E1x3OESCrFull5x5", &phoE1x3OESCrFull5x5_);
	evType.tree->Branch("ph_E2x5OESCrFull5x5", &phoE2x5OESCrFull5x5_);
	evType.tree->Branch("ph_E5x5OESCrFull5x5", &phoE5x5OESCrFull5x5_);
	evType.tree->Branch("ph_EmaxOE3x3Full5x5", &phoEmaxOE3x3Full5x5_);
	evType.tree->Branch("ph_E2ndOE3x3Full5x5", &phoE2ndOE3x3Full5x5_);
	evType.tree->Branch("ph_2x2OE3x3Full5x5", &pho2x2OE3x3Full5x5_);

	evType.tree->Branch("ph_sieie", &phoSigmaIEtaIEta_);
	evType.tree->Branch("ph_sieip", &phoSigmaIEtaIPhi_);
	evType.tree->Branch("ph_sipip", &phoSigmaIPhiIPhi_);
	evType.tree->Branch("ph_sieieOsipip", &phoSieieOSipipFull5x5_);

	evType.tree->Branch("ph_se", &phoEtaWidth_);
	evType.tree->Branch("ph_sp", &phoPhiWidth_);
	evType.tree->Branch("ph_seOsp", &phoEtaWOPhiWFull5x5_);

	evType.tree->Branch("ph_MaxEnergyXtal", &phoMaxEnergyXtal_);
	evType.tree->Branch("ph_E2ndFull5x5", &phoE2ndFull5x5_);
	evType.tree->Branch("ph_E2x2Full5x5", &phoE2x2Full5x5_);
	evType.tree->Branch("ph_E3x3Full5x5", &phoE3x3Full5x5_);
	evType.tree->Branch("ph_E5x5Full5x5", &phoE5x5Full5x5_);
	evType.tree->Branch("ph_E1x3Full5x5", &phoE1x3Full5x5_);
	evType.tree->Branch("ph_E1x5Full5x5", &phoE1x5Full5x5_);
	evType.tree->Branch("ph_E2x5Full5x5", &phoE2x5Full5x5_);

	evType.tree->Branch("ph_PFClusEcalIso", &phoPFClusEcalIso_);
	evType.tree->Branch("ph_PFClusHcalIso", &phoPFClusHcalIso_);
	evType.tree->Branch("ph_TrkSumPtSolidConeDR04", &phoTrkSumPtSolidConeDR04_);
	evType.tree->Branch("ph_TrkSumPtHollowConeDR04", &phoTrkSumPtHollowConeDR04_);
	evType.tree->Branch("ph_TrkSumPtSolidConeDR03", &phoTrkSumPtSolidConeDR03_);
	evType.tree->Branch("ph_TrkSumPtHollowConeDR03", &phoTrkSumPtHollowConeDR03_);
	evType.tree->Branch("ph_ECALIso", &phoECALIso_);
	evType.tree->Branch("ph_HCALIso", &phoHCALIso_);
	evType.tree->Branch("ph_PFECALClusIsoCorr", &phoPFECALClusIsoCorr_);
	evType.tree->Branch("ph_PFHCALClusIsoCorr", &phoPFHCALClusIsoCorr_);
	evType.tree->Branch("ph_TkrIsoCorr", &phoTkrIsoCorr_);

	evType.tree->Branch("ph_hoe", &phoHoverE_);
	evType.tree->Branch("ph_ChIso", &phoPFChIsoRaw_);
	evType.tree->Branch("ph_PhoIso", &phoPFPhoIsoRaw_);
	evType.tree->Branch("ph_NeuIso", &phoPFNeuIsoRaw_);
	evType.tree->Branch("ph_ChWorstIso", &phoPFChWorstIsoRaw_);
	evType.tree->Branch("ph_EGMidMVA", &phoIDMVA_);
	evType.tree->Branch("ph_BDTpred", &phoBDTpred_);

	evType.tree->Branch("ph_PFClusIDbits", &phoPFClusIDbits_);	
	evType.tree->Branch("ph_IDbit", &phoIDbit_);

	evType.tree->Branch("pass95", &pass95_);	

	evType.tree->Branch("ph_MIP", &phoMIP_);
	
	evType.tree->Branch("ph_sc_et", &phoSCet_);
	evType.tree->Branch("ph_sc_rawet", &phoSCrawet_);
	evType.tree->Branch("ph_sc_eta", &phoSCeta_);
	evType.tree->Branch("ph_sc_phi", &phoSCphi_);
	evType.tree->Branch("ph_sc_e", &phoSCEn_);
	evType.tree->Branch("ph_sc_rawe", &phoSCRawEn_);

	evType.tree->Branch("lepVeto", &lepVeto_);

        //met and nGenPho useful branches	
        evType.tree->Branch("genMET", &genMET_);
	evType.tree->Branch("met", &met_);
	evType.tree->Branch("metPhi", &metPhi_);
	evType.tree->Branch("metSig", &metSig_);
	evType.tree->Branch("EtSig", &EtSig_);

	evType.tree->Branch("mT", &mT_);
	evType.tree->Branch("deltaPhiMetPho", &deltaPhiMetPho_);
	evType.tree->Branch("phoPtOverMet", &phoPtOverMet_);

	evType.tree->Branch("nJet30", &nJet30_);
	evType.tree->Branch("Ht30", &Ht30_);
	evType.tree->Branch("minDeltaPhiMetJet30", &minDeltaPhiMetJet30_);
	evType.tree->Branch("minDeltaRPhoJet30", &minDeltaRPhoJet30_);

        evType.tree->Branch("mt", &mt_);
        evType.tree->Branch("mcTrue", &mcTrue_);
        evType.tree->Branch("nGenPho", &nGenPho_);

	std::cout<<"Created output tree:\t"<<typeName<<"\t"<<typeTitle<<std::endl<<std::endl;
	// evType.tree->Print();
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void atgcStudies::fillEventType(eventType & evType){
	evType.tree->Fill();

	registerCutFlow(evType);
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void atgcStudies::registerAllCutFlow(){
	registerCutFlow(fullECAL);
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Reset lastCutStep before each event
void atgcStudies::registerCutFlow(eventType & evType){
	evType.cutFlowCount->Fill(evType.lastCutStep);
	evType.cutFlowGenWeight->Fill(evType.lastCutStep, genWeight_);
	evType.lastCutStep = evType.lastCutStep + 1.;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif
