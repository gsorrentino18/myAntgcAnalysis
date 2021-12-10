//////////////////////////////////////////////////////////////////////
//  Mohammad Abrar Wadud, Univeristy of Minnesota           		//
//	August/04/2020                                            		//
//  Photon ID Study   												//
//////////////////////////////////////////////////////////////////////

#include "/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/macros/extra_tools.cc"

R__ADD_INCLUDE_PATH(/local/cms/user/wadud/aNTGCmet/xgboost/include/xgboost/)
R__LOAD_LIBRARY(/local/cms/user/wadud/aNTGCmet/xgboost/lib/libxgboost.so)
#include </local/cms/user/wadud/aNTGCmet/xgboost/include/xgboost/c_api.h>

#ifndef GENPHOMATCHER
#define GENPHOMATCHER

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
class EBgjetsValidation{
public:
	EBgjetsValidation(std::string FILELIST,
                           std::string OUTFILE,
                           Float_t XSECTION = -1.,
                           std::string MCPILEUPHIST = "",
                           std::string DATAPILEUPHIST = "/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/data/pileupUL17/pileup_2017_data.root",
                           std::string PFECALCLUS_PUCORRECTIONS = "/data/cmszfs1/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/UL17/cutOptimization/bugfix/2017/isoCorrections_woHCAL_hardHoE0p04/phoPFClusEcalIso_RhoCorrections.txt",
                           std::string PFHCALCLUS_PUCORRECTIONS = "/data/cmszfs1/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/UL17/cutOptimization/bugfix/2017/isoCorrections_woHCAL_hardHoE0p04/phoPFClusHcalIso_RhoCorrections.txt",
                           std::string TKRISO_PUCORRECTIONS = "/data/cmszfs1/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/UL17/cutOptimization/bugfix/2017/isoCorrections_woHCAL_hardHoE0p04/phoTrkSumPtHollowConeDR03_RhoCorrections.txt",
                           std::string PFECALCLUS_PTSCALING = "/data/cmszfs1/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/UL17/cutOptimization/bugfix/2017/isoCorrections_woHCAL_hardHoE0p04/phoPFClusEcalIso_PtCorrections.txt",
                           std::string PFHCALCLUS_PTSCALING = "/data/cmszfs1/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/UL17/cutOptimization/bugfix/2017/isoCorrections_woHCAL_hardHoE0p04/phoPFClusHcalIso_PtCorrections.txt",
                           std::string BDT_PATH = "/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/UL17/BDT/training/tuning/v6/EB/aNTGC_photon_BDT_EB_2021_08_26_09_39_52.model");
 
                ~EBgjetsValidation() {
		XGBoosterFree(phoBDT_h);
		std::cout<<"END @ "<<getCurrentTime()<<std::endl;
		std::cout<<"*************************************************************************************************************************************************"<<std::endl;
	};

private:
	Bool_t              isMC = false;
	Bool_t              doPUreweight = false;
	Float_t             xSec = -1.;

        isoCorrMap          ecalIsoRhoCorrMap;
        isoCorrMap          hcalIsoRhoCorrMap;
        isoCorrMap          tkrIsoRhoCorrMap;

        isoCorrMap          ecalIsoPtCorrMap;
        isoCorrMap          hcalIsoPtCorrMap;

	void                analyze();
	Bool_t              selectEvent();

	Short_t         	photonIsPrompt(Short_t _phoIndex, Float_t _deltaRmax, Float_t _relDeltaPtMin, Float_t _relDeltaPtMax);
	Bool_t			photonIsFake(Short_t _phoIndex, Float_t _deltaRmax = 0.3);
	Short_t 		matchWithRecoPho(Short_t _genIndex, Float_t _deltaRmax, Float_t _relDeltaPtMin, Float_t _relDeltaPtMax);
        Short_t                 matchWithRecoEle(Short_t _genIndex, Float_t _deltaRmax, Float_t _relDeltaPtMin, Float_t _relDeltaPtMax);
	Short_t 		matchWithTrigPho(Short_t _phoIndex, Float_t _deltaRmax, Float_t _relDeltaPtMin, Float_t _relDeltaPtMax);

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
        TTreeReaderAnyValue<ULong64_t>          _HLTMuX;
	TTreeReaderAnyValue<UShort_t>           _beamHaloSummary;
	TTreeReaderAnyValue<UShort_t>           _metFilters;

	TTreeReaderAnyValue<UChar_t>            _puTrue;
	TTreeReaderAnyValue<Float_t>            _genWeight;
	TTreeReaderAnyValue<UShort_t>		_nMC;
	TTreeReaderVectorValue<Int_t>           _mcPID;
	TTreeReaderVectorValue<Float_t>		_mcPt;
	TTreeReaderVectorValue<Float_t>         _mcEta;
	TTreeReaderVectorValue<Float_t>         _mcPhi;
	TTreeReaderVectorValue<UShort_t>        _mcStatusFlag;
        TTreeReaderVectorValue<Char_t>          _mcPromptStatusType;
	TTreeReaderVectorValue<Short_t>         _mcStatus;
	TTreeReaderVectorValue<Short_t>     	_mcIndex;

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
        TTreeReaderVectorValue<Float_t>         _phoESEffSigmaRR;

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
	TTreeReaderVectorValue<Float_t>         _ecalSC_eta;
	TTreeReaderVectorValue<Float_t>         _ecalSC_phi;
	TTreeReaderVectorValue<Float_t>         _ecalSC_En;
	TTreeReaderVectorValue<Float_t>         _ecalSC_RawEn;
	TTreeReaderVectorValue<Float_t>         _ecalSC_etaWidth;
	TTreeReaderVectorValue<Float_t>         _ecalSC_phiWidth;

	TTreeReaderAnyValue<Float_t>            _pfMET;
	TTreeReaderAnyValue<Float_t>            _pfMETPhi;
	TTreeReaderAnyValue<Float_t>            _pfMET_metSig;
	TTreeReaderAnyValue<Float_t>            _pfMET_EtSig;

	TTreeReaderAnyValue<UShort_t>           _nEle;
        TTreeReaderVectorValue<Float_t>         _elePt;
        TTreeReaderVectorValue<Float_t>         _eleEta;
        TTreeReaderVectorValue<Float_t>         _elePhi;
	TTreeReaderVectorValue<Float_t>         _eleCalibPt;
        TTreeReaderVectorValue<Float_t>         _eleCalibEn;
	TTreeReaderVectorValue<UInt_t>         _eleIDbit;
        TTreeReaderVectorValue<Float_t>         _eleHoverE;

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
	TTreeReaderVectorValue<UInt_t>		_trgObjPhoBits1;
	TTreeReaderVectorValue<Float_t>		_trgObjPhoPt;
	TTreeReaderVectorValue<Float_t>		_trgObjPhoEta;
	TTreeReaderVectorValue<Float_t>		_trgObjPhoPhi;
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////// Output variables //////////////////////////////////////////////////////////

	Int_t		 run_;
	Long64_t	 event_;
	UShort_t	 lumis_;	
	Float_t	 	 rho_;
	UChar_t	 	 nVtx_;

	UShort_t       	metFilters_;
	UShort_t	beamHaloSummary_;

	Float_t         genWeight_ = 1.;
	Float_t         puWeight_ = 1.;
	UChar_t		puTrue_;
	UShort_t	genStatusFlag_;
	Short_t		genStatus_;
	Float_t 	deltaRgenPho_;
	Float_t 	relDeltaPtGenPho_;
	Float_t 	deltaRPt_;
	Int_t		genPDGid_ = 0;

	Float_t 	deltaRtrg_;
	Float_t 	deltaPttrg_;

	Float_t 	phoPt_;
	Float_t 	phoEta_;
	Float_t 	phoPhi_;
        Float_t         phoEn_;
	Float_t 	phoSeedTime_;

	UChar_t		phoGenBits_;
	UChar_t		phoQualityBits_;
	Float_t 	phoR9Full5x5_;
	Float_t 	phoS4Full5x5_;
        Float_t         phoESEffSigmaRR_;

	Float_t 	phoEmaxOESCrFull5x5_;
	Float_t 	phoE2ndOESCrFull5x5_;
	Float_t 	phoE2ndOEmaxFull5x5_;
	Float_t 	phoE1x3OESCrFull5x5_;
	Float_t 	phoE2x5OESCrFull5x5_;
	Float_t 	phoE5x5OESCrFull5x5_;

	Float_t 	phoEmaxOE3x3Full5x5_;
	Float_t 	phoE2ndOE3x3Full5x5_;
	Float_t 	pho2x2OE3x3Full5x5_;
	Float_t 	phoSieieOSipipFull5x5_;
	Float_t		phoEtaWOPhiWFull5x5_;

	Float_t 	phoSigmaIEtaIEta_;
	Float_t 	phoSigmaIPhiIPhi_;
	Float_t 	phoSigmaIEtaIPhi_;

	Float_t 	phoE2x2Full5x5_;
	Float_t 	phoE3x3Full5x5_;
	Float_t 	phoE5x5Full5x5_;
	Float_t		phoMaxEnergyXtal_;
	Float_t		phoE2ndFull5x5_;
	Float_t		phoE1x3Full5x5_;
	Float_t		phoE1x5Full5x5_;
	Float_t		phoE2x5Full5x5_;

	Float_t		phoPFClusEcalIso_;
	Float_t		phoPFClusHcalIso_;
	Float_t		phoTrkSumPtSolidConeDR04_;
	Float_t		phoTrkSumPtHollowConeDR04_;
	Float_t		phoTrkSumPtSolidConeDR03_;
	Float_t		phoTrkSumPtHollowConeDR03_;
	Float_t		phoECALIso_;
	Float_t		phoHCALIso_;
	Float_t 	phoPFECALClusIsoCorr_;
	Float_t 	phoPFHCALClusIsoCorr_;
	Float_t 	phoTkrIsoCorr_;	

	Float_t 	phoHoverE_;
	Float_t 	phoPFChIsoRaw_;
	Float_t 	phoPFPhoIsoRaw_;
	Float_t 	phoPFNeuIsoRaw_;
	Float_t 	phoPFChWorstIsoRaw_;
	Float_t 	phoIDMVA_;
	UChar_t         phoIDbit_;
	Float_t 	phoBDTpred_;
        Float_t         phoBDTpredPosShift_[20];
        Float_t         phoBDTpredNegShift_[20];
	UChar_t 	phoPFClusIDbits_;

	Bool_t 		pass95_ = 0;

	Float_t 	phoMIP_;	

	Float_t 	phoSCet_;
	Float_t 	phoSCrawet_;
	Float_t 	phoSCeta_;
	Float_t 	phoSCphi_;
	Float_t 	phoSCEn_;
	Float_t 	phoSCRawEn_;
	Float_t 	phoEtaWidth_;
	Float_t 	phoPhiWidth_;

	UChar_t 	lepVeto_;

	Float_t 	met_;
	Float_t		metPhi_;
	Float_t 	metSig_;
	Float_t 	EtSig_;

	Float_t 	mT_;
	Float_t 	deltaPhiMetPho_;
	Float_t 	phoPtOverMet_;


	UChar_t 	nJet30_;
	Float_t 	Ht30_;
	Float_t 	minDeltaPhiMetJet30_;
	Float_t 	minDeltaRPhoJet30_;
        Float_t 	jetphoDeltaPhi_;
        Float_t         jetphoDeltaR_;
        Int_t           nPhoCand_=0;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////Buffer variables///////////////////////////////////////////////////////
	Short_t			phoTrigMatch;
	Short_t			matchedGenPhoIndex;
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
	eventType		fullEB;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
EBgjetsValidation::EBgjetsValidation(std::string FILELIST, std::string OUTFILE, Float_t XSECTION, std::string MCPILEUPHIST, std::string DATAPILEUPHIST,
    std::string PFECALCLUS_PUCORRECTIONS, std::string PFHCALCLUS_PUCORRECTIONS, std::string TKRISO_PUCORRECTIONS,
    std::string PFECALCLUS_PTSCALING, std::string PFHCALCLUS_PTSCALING,
    std::string BDT_PATH){

	std::cout<<"*************************************************************************************************************************************************"<<std::endl<<
	getCurrentTime()<<std::endl<<
	"Running EBgjetsValidation"<<std::endl<<
        "\n\nInput parameters:" << std::endl <<
            "\t\tFile list = " << FILELIST << std::endl <<
            "\t\tOutput file = " << OUTFILE << std::endl <<
            "\t\tCross section = " << XSECTION << std::endl <<
            "\t\tMC pileup histogram = " << MCPILEUPHIST << std::endl <<
            "\t\tData pileup histogram = " << DATAPILEUPHIST << std::endl <<
            "\t\tECAL PFCLuster Isolation Effective Areas = " << PFECALCLUS_PUCORRECTIONS << std::endl <<
            "\t\tHCAL PFCLuster Isolation Effective Areas = " << PFHCALCLUS_PUCORRECTIONS << std::endl <<
            "\t\tTracker Isolation Effective Areas = " << TKRISO_PUCORRECTIONS << std::endl <<
            "\t\tECAL PFCLuster Isolation pT Scaling = " << PFECALCLUS_PTSCALING << std::endl <<
            "\t\tHCAL PFCLuster Isolation pT Scaling = " << PFHCALCLUS_PTSCALING << std::endl <<
            "\t\tBDT model file = " << BDT_PATH << std::endl;

	xSec = XSECTION;
	if(XSECTION > 0.) isMC = true;
	if(isMC && file_exists(MCPILEUPHIST) && file_exists(DATAPILEUPHIST)) doPUreweight = true;

	std::cout<<"\t\tSample is simulation = "<<std::boolalpha<<isMC<<std::endl<<
	"\t\tDo pileup reweight = "<<doPUreweight<<"\n\n"<<std::endl;

        if (doPUreweight) {
          std::cout << "Pileup reweighting:" << std::endl;
          puReweighter.init(MCPILEUPHIST, DATAPILEUPHIST, "hPUTruew", "pileup");
        }
      
      
        if (file_exists(split_string(PFECALCLUS_PUCORRECTIONS)[0])) {
          std::cout << "PF ECAL Cluster Isolation pileup corrections:" << std::endl;
          ecalIsoRhoCorrMap.init(split_string(PFECALCLUS_PUCORRECTIONS)[0], 2);
        }
      
      
        if (file_exists(split_string(PFHCALCLUS_PUCORRECTIONS)[0])) {
          std::cout << "PF HCAL Cluster Isolation pileup correction:" << std::endl;
          hcalIsoRhoCorrMap.init(split_string(PFHCALCLUS_PUCORRECTIONS)[0], 2);
        }
      
        if (file_exists(split_string(TKRISO_PUCORRECTIONS)[0])) {
          std::cout << "Tracker Isolation pileup corrections:" << std::endl;
          tkrIsoRhoCorrMap.init(split_string(TKRISO_PUCORRECTIONS)[0], 2);
        }
      
        if (file_exists(split_string(PFECALCLUS_PTSCALING)[0])) {
          std::cout << "PF ECAL Cluster Isolation pT scaling corrections:" << std::endl;
          ecalIsoPtCorrMap.init(split_string(PFECALCLUS_PTSCALING)[0], 2);
        }
      
        if (file_exists(split_string(PFHCALCLUS_PTSCALING)[0])) {
          std::cout << "PF HCAL Cluster Isolation pT scaling corrections:" << std::endl;
          hcalIsoPtCorrMap.init(split_string(PFHCALCLUS_PTSCALING)[0], 2);
        }
      
        if (file_exists(BDT_PATH)) {
          std::cout << "\nLoading Photon EE BDT model from " << BDT_PATH << std::endl;
          XGBoosterCreate(NULL, 0, &phoBDT_h);
          XGBoosterSetParam(phoBDT_h, "seed", "0");
          Int_t mLdSuccess = XGBoosterLoadModel(phoBDT_h, BDT_PATH.c_str());
          if (mLdSuccess == 0) predictBDT = 1;
          else {
            std::cout << "Failed to load Photon EE BDT model!" << std::endl;
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
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Short_t	EBgjetsValidation::matchWithRecoPho(Short_t _genIndex, Float_t _deltaRmax, Float_t _relDeltaPtMin, Float_t _relDeltaPtMax){

	Short_t matchedRecoPho = -999;
	Float_t minDeltaR = 999.;

	for(UShort_t iPho = 0; iPho < _nPho; iPho++){

		UChar_t iPhoFidReg = _phoFiducialRegion[iPho];
		if(!getBit(iPhoFidReg, 0)) continue; 			// skip if not EB (0 = EB, 1 = EE, 2 = EB-EE gap)

                if(_phoCalibEt[iPho] < 200.) continue;
                if(_phoHoverE[iPho] > 0.05 ) continue;

		Float_t relDeltaPtiGenPho = (_phoCalibEt[iPho] - _mcPt[_genIndex])/_mcPt[_genIndex];
		if(relDeltaPtiGenPho > _relDeltaPtMax) continue;
		if(relDeltaPtiGenPho < _relDeltaPtMin) continue;

		Float_t dRiGenPho = deltaR(_phoEta[iPho], _phoPhi[iPho], _mcEta[_genIndex], _mcPhi[_genIndex]);
		if(dRiGenPho > _deltaRmax) continue;

		if(dRiGenPho < minDeltaR){
			matchedRecoPho = iPho;
			minDeltaR = dRiGenPho;
		}
	}

	return matchedRecoPho;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Short_t EBgjetsValidation::matchWithRecoEle(Short_t _genIndex, Float_t _deltaRmax, Float_t _relDeltaPtMin, Float_t _relDeltaPtMax){

        Short_t matchedRecoEle = -999;
        Float_t minDeltaR = 999.;

        for(UShort_t iEle = 0; iEle < _nEle; iEle++){

                UShort_t iHEEP = _eleIDbit[iEle];
                if(!getBit(iHEEP,4)) continue;

                if (_elePt[iEle] < 60) continue;
                 
                Float_t relDeltaPtiGenEle = (_elePt[iEle] - _mcPt[_genIndex])/_mcPt[_genIndex];
                if(relDeltaPtiGenEle > _relDeltaPtMax) continue;
                if(relDeltaPtiGenEle < _relDeltaPtMin) continue;
                Float_t dRiGenEle = deltaR(_eleEta[iEle], _elePhi[iEle], _mcEta[_genIndex], _mcPhi[_genIndex]);
                if(dRiGenEle > _deltaRmax) continue;

                if(dRiGenEle < minDeltaR){
                        matchedRecoEle = iEle;
                        minDeltaR = dRiGenEle;
                //std::cout << "electron index: " << iEle << "    over all electrons in the event: " << _nEle << endl;
                }
        }

        return matchedRecoEle;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Short_t	EBgjetsValidation::matchWithTrigPho(Short_t _phoIndex, Float_t _deltaRmax, Float_t _relDeltaPtMin, Float_t _relDeltaPtMax){
	Float_t minDeltaR = 999.;
	Short_t matchedTrigPho = -999;

	for(UShort_t iTrgPho = 0; iTrgPho < _ntrgObjPho; iTrgPho++){

		Float_t relDeltaPt = (_phoCalibEt[_phoIndex] - _trgObjPhoPt[iTrgPho])/_trgObjPhoPt[iTrgPho];
		if(relDeltaPt > _relDeltaPtMax) continue;
		if(relDeltaPt < _relDeltaPtMin) continue;

		Float_t dR = deltaR(_phoEta[_phoIndex], _phoPhi[_phoIndex], _trgObjPhoEta[iTrgPho], _trgObjPhoPhi[iTrgPho]);
		if(dR > _deltaRmax) continue;

		if(dR < minDeltaR){
			matchedTrigPho = iTrgPho;
			minDeltaR = dR;
		}
	}

	return matchedTrigPho;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Bool_t EBgjetsValidation::selectEvent(){
	//// reset event cut flow
	fullEB.lastCutStep = 0.;
	registerAllCutFlow();
        return 1;

};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Bool_t EBgjetsValidation::initEventTypes(){
	initEventType(fullEB, "tnpPhoIDs", "Full ECAL Barrel");
	return 1;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EBgjetsValidation::analyze(){
	std::cout<<"--------------------------------------------------------------------------------------------------"<<std::endl<<
	getCurrentTime()<<std::endl<<
	"Analyzing events.."<<std::endl;
	ULong64_t current_entry =0;
        Int_t fp=0;
	while(inputTTreeReader.Next()){

		if(current_entry % REPORT_EVERY == 0){
			std::cout<<"\t"<< getCurrentTime()<<"\tAnalyzing entry\t"<<current_entry<<
			",\t\tevent\t"<<(_event)<<"\t\tFile " <<inputTree->GetCurrentFile()->GetName() <<std::endl;
		}

                selectEvent();

                //to be away from the signal region, following Shilpi code
                if (_pfMET > 160) continue;

                //single electron trigger
                bool ispassHLT = _HLTPho>>22&1; //17 HLT_Ele35_WPTight_Gsf_v, 22 HLT_Photon200_v
                if (!ispassHLT) continue; 

		rhoPreweight.Fill(_rho, genWeight_);
		nvtxPreweight.Fill(_nVtx, genWeight_);

                Float_t weight = 1;

		if(isMC) {

                        std::string Map = "/local/cms/user/gsorrent/antgc_analysis/samples/aNTGC_Samples-ntuples_2017UL.csv";
                        const char* filePathTemp = inputTree->GetDirectory()->GetName();
                        std::string filePath = filePathTemp;
                        std::string sample = findAndReplaceAll(filePath, "/hdfs/cms/user/wadud/anTGC/ntuplesUL/ntuples2017UL/", "");
                        std::string sampleName = sample.substr(0, sample.find("/"));

                        Double_t xSection = std::stof(vLookup(sampleName, Map, 0, 2));
                        Double_t sumGenWeight = 1;
                        
                        Double_t lumi = 41.48; //UL17

                        Float_t norm = xSection * 1000 * lumi/sumGenWeight;

			if(doPUreweight){
				genWeight_ = _genWeight;
				puWeight_ = puReweighter.weight(_puTrue) * genWeight_ * norm; //puWeight + MC normalization + genWeight
				Float_t genPUweight = genWeight_ * puWeight_;
                                weight = genPUweight;
				pileupPostweight.Fill(_puTrue, genPUweight);
				rhoPostweight.Fill(_rho, genPUweight);
				nvtxPostweight.Fill(_nVtx, genPUweight);
			} else{
				genWeight_ = _genWeight;
				puWeight_ = 1.;    
			}
			pileupPreweight.Fill(_puTrue, genWeight_);
		}

               Short_t highestJetIndex = -999;
               Float_t highestPtJet = -999;
               Short_t highestPhoIndex = -999;
               Float_t highestPtPho = -999;
               Short_t goodJet = -999;
               Short_t goodPho = -999;
               Bool_t  PhoJetfound = false;
               Int_t nCand = 0;

               for(UShort_t iPho=0; iPho<_nPho; iPho++) {

                  if (_phoCalibEt[iPho] < 200 )  continue;
                  if ( fabs(_phoEta[iPho]) > BETRetaMin ) continue; //EB selection
                  if (getBit((_phoQualityBits[iPho]),0)) continue;  //pixel seed veto (we want real photons)
                  nCand++;

                  if (_phoCalibEt[iPho] > highestPtPho) {
                        highestPtPho = _phoCalibEt[iPho];
                        highestPhoIndex = iPho;
                     }
 
                  goodPho=highestPhoIndex;
                  //apply later (in separate macro)  the isolation cuts
               }

               for (UShort_t iJet=0; iJet<_nAK4CHSJet; iJet++) {
                  if (_AK4CHSJet_Pt[iJet] < 200 )  continue;
                  if ( fabs(_AK4CHSJet_Eta[iJet]) > 5 ) continue; //EB selection
                  Char_t iJetID = _AK4CHSJet_ID[iJet];
                  if (!getBit(iJetID, 1)) continue;               // tight ID

                  if (_AK4CHSJet_Pt[iJet] > highestPtJet) {
                     highestJetIndex = iJet;
                     highestPtJet = _AK4CHSJet_Pt[iJet];
                  }
                  goodJet=highestJetIndex;
               }
                  if ((goodPho >= 0) && (goodJet >= 0)) PhoJetfound = true;

               if (PhoJetfound) {

                  fp++;

                  run_                            = _run;
                  event_                          = _event;
                  lumis_                          = _lumis;
                  rho_                            = _rho;
                  nVtx_                           = _nVtx;

                  metFilters_                     = _metFilters;
                  beamHaloSummary_                = _beamHaloSummary;

                  Short_t phoSCindex              = _phoDirectEcalSCindex[goodPho];
                  phoSCeta_                       = _ecalSC_eta[phoSCindex];
                  Float_t phoAbsSCEta_            = std::abs(phoSCeta_);

                  phoPt_                          = _phoCalibEt[goodPho];
                  phoEta_                         = _phoEta[goodPho];
                  phoPhi_                         = _phoPhi[goodPho];
                  phoEn_                          = _phoCalibE[goodPho];
                  phoSeedTime_                    = _phoSeedTime[goodPho];

                  phoQualityBits_                 = _phoQualityBits[goodPho];
                  phoR9Full5x5_                   = _phoR9Full5x5[goodPho];
       	          phoS4Full5x5_                   = _phoE2x2Full5x5[goodPho]/_ecalSC_RawEn[phoSCindex];
                  phoEmaxOESCrFull5x5_            = _phoMaxEnergyXtal[goodPho]/_ecalSC_RawEn[phoSCindex];
                  phoE2ndOESCrFull5x5_            = _phoE2ndFull5x5[goodPho]/_ecalSC_RawEn[phoSCindex];
                  phoE2ndOEmaxFull5x5_            = _phoE2ndFull5x5[goodPho]/_phoMaxEnergyXtal[goodPho];
                  phoE1x3OESCrFull5x5_            = _phoE1x3Full5x5[goodPho]/_ecalSC_RawEn[phoSCindex];
                  phoE2x5OESCrFull5x5_            = _phoE2x5Full5x5[goodPho]/_ecalSC_RawEn[phoSCindex];
                  phoE5x5OESCrFull5x5_            = _phoE5x5Full5x5[goodPho]/_ecalSC_RawEn[phoSCindex];
                  phoESEffSigmaRR_                = _phoESEffSigmaRR[goodPho];

                  phoSigmaIEtaIEta_               = _phoSigmaIEtaIEtaFull5x5[goodPho];
                  phoSigmaIEtaIPhi_               = _phoSigmaIEtaIPhiFull5x5[goodPho];
                  phoSigmaIPhiIPhi_               = _phoSigmaIPhiIPhiFull5x5[goodPho];

                  phoE2x2Full5x5_                 = _phoE2x2Full5x5[goodPho];
                  phoE3x3Full5x5_                 = phoR9Full5x5_ * _ecalSC_RawEn[phoSCindex];
                  phoE5x5Full5x5_                 = _phoE5x5Full5x5[goodPho];
                  phoMaxEnergyXtal_               = _phoMaxEnergyXtal[goodPho];
                  phoE2ndFull5x5_                 = _phoE2ndFull5x5[goodPho];
                  phoE1x3Full5x5_                 = _phoE1x3Full5x5[goodPho];
                  phoE1x5Full5x5_                 = _phoE1x5Full5x5[goodPho];
                  phoE2x5Full5x5_                 = _phoE2x2Full5x5[goodPho];

                  phoEmaxOE3x3Full5x5_            = phoMaxEnergyXtal_/phoE3x3Full5x5_;
                  phoE2ndOE3x3Full5x5_            = phoE2ndFull5x5_/phoE3x3Full5x5_;
                  pho2x2OE3x3Full5x5_             = phoE2x2Full5x5_/phoE3x3Full5x5_;
                  phoSieieOSipipFull5x5_          = phoSigmaIEtaIEta_/phoSigmaIPhiIPhi_;
                  phoEtaWidth_                    = _ecalSC_etaWidth[phoSCindex];
                  phoPhiWidth_                    = _ecalSC_phiWidth[phoSCindex];

                  phoPFClusEcalIso_               = _phoPFClusEcalIso[goodPho];
      	          phoPFClusHcalIso_               = _phoPFClusHcalIso[goodPho];
                  phoTrkSumPtSolidConeDR04_       = _phoTrkSumPtSolidConeDR04[goodPho];
                  phoTrkSumPtHollowConeDR04_      = _phoTrkSumPtHollowConeDR04[goodPho];
                  phoTrkSumPtSolidConeDR03_       = _phoTrkSumPtSolidConeDR03[goodPho];
                  phoTrkSumPtHollowConeDR03_      = _phoTrkSumPtHollowConeDR03[goodPho];
                  phoECALIso_                     = _phoECALIso[goodPho];
                  phoHCALIso_                     = _phoHCALIso[goodPho];

                  phoPFECALClusIsoCorr_           = phoPFClusEcalIso_ - ecalIsoRhoCorrMap.getIsoCorr(phoAbsSCEta_, rho_, 0) - ecalIsoPtCorrMap.getIsoCorr(phoAbsSCEta_, phoPt_, 1);
                  phoPFHCALClusIsoCorr_           = phoPFClusHcalIso_ - hcalIsoRhoCorrMap.getIsoCorr(phoAbsSCEta_, rho_, 0) - hcalIsoPtCorrMap.getIsoCorr(phoAbsSCEta_, phoPt_, 1);
                  phoTkrIsoCorr_                  = phoTrkSumPtHollowConeDR03_ - tkrIsoRhoCorrMap.getIsoCorr(phoAbsSCEta_, rho_, 1);

                  phoHoverE_                      = _phoHoverE[goodPho];
                  phoPFChIsoRaw_                  = _phoPFChIso[goodPho];
                  phoPFPhoIsoRaw_                 = _phoPFPhoIso[goodPho];
                  phoPFNeuIsoRaw_                 = _phoPFNeuIso[goodPho];
                  phoPFChWorstIsoRaw_             = _phoPFChWorstIso[goodPho];
                  phoIDMVA_                       = _phoIDMVA[goodPho];
                  phoIDbit_                       = _phoIDbit[goodPho];
                  phoMIP_                         = _phoMIPTotEnergy[goodPho];

                  phoSCet_                        = (_ecalSC_En[phoSCindex]) / std::cosh(phoSCeta_);
                  phoSCrawet_                     = (_ecalSC_RawEn[phoSCindex]) / std::cosh(phoSCeta_);
                  phoSCphi_                       = _ecalSC_phi[phoSCindex];
                  phoSCEn_                        = _ecalSC_En[phoSCindex];
                  phoSCRawEn_                     = _ecalSC_RawEn[phoSCindex];
                  phoEtaWidth_                    = _ecalSC_etaWidth[phoSCindex];
                  phoPhiWidth_                    = _ecalSC_phiWidth[phoSCindex];

                  phoEtaWOPhiWFull5x5_            = phoEtaWidth_/phoPhiWidth_;

                  Float_t 	scE                    = _ecalSC_En[phoSCindex];
                  
               

                  if(predictBDT){


                     const float *prediction;
                     bst_ulong out_len;

                     std::vector<Float_t> feats{_phoE2x2Full5x5[goodPho] / (_phoR9Full5x5[goodPho] * _ecalSC_RawEn[phoSCindex]),
                                  phoAbsSCEta_,
                                  _phoE1x3Full5x5[goodPho] / _ecalSC_RawEn[phoSCindex],
                                  _phoE2ndFull5x5[goodPho] / _ecalSC_RawEn[phoSCindex],
                                  _phoE2x5Full5x5[goodPho] / _ecalSC_RawEn[phoSCindex],
                                  _phoMaxEnergyXtal[goodPho] / _ecalSC_RawEn[phoSCindex],
                                  _ecalSC_etaWidth[phoSCindex] / _ecalSC_phiWidth[phoSCindex],
                                  _ecalSC_etaWidth[phoSCindex],
                                  _ecalSC_phiWidth[phoSCindex],
                                  _phoCalibEt[phoSCindex],
                                  _phoR9Full5x5[goodPho],
                                  _phoE2x2Full5x5[goodPho] / _ecalSC_RawEn[phoSCindex],
                                  _phoSigmaIEtaIEtaFull5x5[goodPho] / _phoSigmaIPhiIPhiFull5x5[goodPho],
                                  _phoSigmaIEtaIEtaFull5x5[goodPho],
                                  _phoSigmaIEtaIPhiFull5x5[goodPho],
                                  _phoSigmaIPhiIPhiFull5x5[goodPho]};

                     XGDMatrixCreateFromMat((float*)feats.data(), 1, feats.size(), -9999999999, &dTest);
                     XGBoosterPredict(phoBDT_h, dTest, 0, 0, 0, &out_len, &prediction);
                     assert(out_len == 1);
                     phoBDTpred_ = prediction[0];
                     XGDMatrixFree(dTest);

                     pass95_ = (phoBDTpred_ > 0.8361 && phoHoverE_ < 0.04012  &&  phoPFECALClusIsoCorr_ < 1.84  && phoTkrIsoCorr_ < 1.63);
                     if(pass95_) setBit(phoPFClusIDbits_, 3, 1);
                  }

                  jetphoDeltaPhi_ = fabs(_AK4CHSJet_Phi[goodJet] - _phoPhi[goodPho]);
                  if (jetphoDeltaPhi_ > acos(-1)) {
                      jetphoDeltaPhi_ = 2*acos(-1) - jetphoDeltaPhi_;
                  }
                  Double_t deltaEta = _AK4CHSJet_Eta[goodJet] - _phoEta[goodPho];
                  Double_t deltaPt = _AK4CHSJet_Pt[goodJet] - _phoCalibEt[goodPho];
                  jetphoDeltaR_ = sqrt((jetphoDeltaPhi_*jetphoDeltaPhi_)  + (deltaEta*deltaEta));
                  nPhoCand_ = nCand;

                  fillEventType(fullEB);
               } //tag and probe > 0 cycle 
     
		current_entry++; 
                //cout << "Weight: " << weight << endl;
	};
        std::cout << "Number of FoundPair: " << fp << std::endl;
	std::cout<<"Done analyzing!"<<std::endl<<
	getCurrentTime()<<std::endl<<
	"--------------------------------------------------------------------------------------------------"<<std::endl;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Bool_t EBgjetsValidation::initNtuples(std::string FILELIST){

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
        _HLTMuX.set(inputTTreeReader, "HLTMuX");
	if(isMC){
		_puTrue.set(inputTTreeReader, "puTrue");
		_genWeight.set(inputTTreeReader, "genWeight");
		_nMC.set(inputTTreeReader, "nMC");
		_mcPID.set(inputTTreeReader, "mcPID");
		_mcPt.set(inputTTreeReader, "mcPt");
		_mcEta.set(inputTTreeReader, "mcEta");
		_mcPhi.set(inputTTreeReader, "mcPhi");
		_mcStatusFlag.set(inputTTreeReader, "mcStatusFlag");
                _mcPromptStatusType.set(inputTTreeReader, "mcPromptStatusType");
		_mcStatus.set(inputTTreeReader, "mcStatus");
		_mcIndex.set(inputTTreeReader, "mcIndex");
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
        _phoESEffSigmaRR.set(inputTTreeReader, "phoESEffSigmaRR");

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
	_ecalSC_eta.set(inputTTreeReader, "ecalSC_eta");
	_ecalSC_phi.set(inputTTreeReader, "ecalSC_phi");
	_ecalSC_En.set(inputTTreeReader, "ecalSC_En");
	_ecalSC_RawEn.set(inputTTreeReader, "ecalSC_RawEn");
	_ecalSC_etaWidth.set(inputTTreeReader, "ecalSC_etaWidth");
	_ecalSC_phiWidth.set(inputTTreeReader, "ecalSC_phiWidth");

	_metFilters.set(inputTTreeReader, "metFilters");
	_pfMET.set(inputTTreeReader, "pfMET");
	_pfMETPhi.set(inputTTreeReader, "pfMETPhi");
	_pfMET_metSig.set(inputTTreeReader, "pfMET_metSig");
	_pfMET_EtSig.set(inputTTreeReader, "pfMET_EtSig");

	_nEle.set(inputTTreeReader, "nEle");
        _elePt.set(inputTTreeReader, "elePt");
        _eleEta.set(inputTTreeReader, "eleEta");
        _elePhi.set(inputTTreeReader, "elePhi");
	_eleCalibPt.set(inputTTreeReader, "eleCalibPt");
        _eleCalibEn.set(inputTTreeReader, "eleCalibEn");
	_eleIDbit.set(inputTTreeReader, "eleIDbit");
        _eleHoverE.set(inputTTreeReader, "eleHoverE");

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
	_trgObjPhoBits1.set(inputTTreeReader, "trgObjPhoBits1");
	_trgObjPhoPt.set(inputTTreeReader, "trgObjPhoPt");
	_trgObjPhoEta.set(inputTTreeReader, "trgObjPhoEta");
	_trgObjPhoPhi.set(inputTTreeReader, "trgObjPhoPhi");

	std::cout<<"Branches initialized!"<<std::endl<<
	"**************************************************************************************************************************************************************"<<std::endl;
	return 1;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EBgjetsValidation::initEventType(eventType & evType, std::string typeName, std::string typeTitle){

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

	evType.tree->Branch("ph_et", &phoPt_);
	evType.tree->Branch("ph_eta", &phoEta_);
	evType.tree->Branch("ph_phi", &phoPhi_);
        evType.tree->Branch("ph_e", &phoEn_);
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
        evType.tree->Branch("phoESEffSigmaRR", &phoESEffSigmaRR_);

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

        evType.tree->Branch("jetphoDeltaPhi", &jetphoDeltaPhi_);
        evType.tree->Branch("jetphoDeltaR", &jetphoDeltaR_);
        evType.tree->Branch("nPhoCand", &nPhoCand_);
        	
	std::cout<<"Created output tree:\t"<<typeName<<"\t"<<typeTitle<<std::endl<<std::endl;
	// evType.tree->Print();
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EBgjetsValidation::fillEventType(eventType & evType){
	evType.tree->Fill();

	//registerCutFlow(evType);
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EBgjetsValidation::registerAllCutFlow(){
	registerCutFlow(fullEB);
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Reset lastCutStep before each event
void EBgjetsValidation::registerCutFlow(eventType & evType){
	evType.cutFlowCount->Fill(evType.lastCutStep);
	evType.cutFlowGenWeight->Fill(evType.lastCutStep, genWeight_);
	evType.lastCutStep = evType.lastCutStep + 1.;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif
