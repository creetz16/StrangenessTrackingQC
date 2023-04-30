/// ROOT
#include "iostream"
#include "ostream"
#include "TStyle.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TList.h"
#include "TLine.h"
#include "TGeoGlobalMagField.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStopwatch.h"
#include "TTree.h"
#include "TSystemDirectory.h"
#include "TMath.h"
#include "TString.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLorentzVector.h"

/// O2
#include "CommonDataFormat/RangeReference.h"
#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/Propagator.h"
#include "ITSBase/GeometryTGeo.h"
#include "ITStracking/IOUtils.h"
#include "ReconstructionDataFormats/Cascade.h"
#include "ReconstructionDataFormats/PrimaryVertex.h"
#include "ReconstructionDataFormats/TrackTPCITS.h"
#include "ReconstructionDataFormats/V0.h"
#include "ReconstructionDataFormats/VtxTrackIndex.h"
#include "ReconstructionDataFormats/PID.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCEventLabel.h"
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "ITSMFTSimulation/Hit.h"
#include "GPUCommonArray.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsCommonDataFormats/DetectorNameConf.h"
#include "StrangenessTracking/StrangenessTracker.h"


using namespace o2;
using namespace o2::framework;

using StrangeTrack = o2::dataformats::StrangeTrack;
using ClusAttachments = o2::strangeness_tracking::ClusAttachments;

void StrangeTrackTableProducer() {

    TFile OutputFile = TFile("StrangeTracks.root", "recreate");
    TTree *Tree = new TTree("StrTTree", "StrTTree");
    float Pt, XKF, YKF, ZKF, PtKF, MassKF, TopoChi2KF, GeoChi2KF;

    Tree->Branch("Pt", &Pt);
    Tree->Branch("XKF", &XKF);
    Tree->Branch("YKF", &YKF);
    Tree->Branch("ZKF", &ZKF);
    Tree->Branch("PtKF", &PtKF);
    Tree->Branch("MassKF", &MassKF);
    Tree->Branch("TopoChi2KF", &TopoChi2KF);
    Tree->Branch("GeoChi2KF", &GeoChi2KF);

    auto fStrangeTracks = TFile::Open("/home/ceres/reetz/sim/StrangenessTracking/GeneratorPYTHIAbottle/sims/000/tf1/o2_strange_tracks.root");
	auto treeStrangeTracks = (TTree*)fStrangeTracks->Get("o2sim");

    std::vector<StrangeTrack> *strangeTrackVec = nullptr;
    std::vector<ClusAttachments> *nAttachments = nullptr;

    treeStrangeTracks->SetBranchAddress("StrangeTracks", &strangeTrackVec);
    treeStrangeTracks->SetBranchAddress("ClusUpdates", &nAttachments);

    // default tree values
    Pt=-1, XKF=-1, YKF=-1, ZKF=-1, PtKF=-1, MassKF=-1, TopoChi2KF=-1, GeoChi2KF=-1; 

    cout << "strangeTrackVec->size() = " << strangeTrackVec->size() << endl;
    cout << "nAttachments->size() = " << nAttachments->size() << endl;


    for (unsigned int iStTr = 0; iStTr < strangeTrackVec->size(); iStTr++) {
        auto &strangeTrack = strangeTrackVec->at(iStTr);
        auto &itsRef = strangeTrack.mITSRef;
        auto &v0Ref = strangeTrack.mDecayRef;
        if (!(strangeTrack.mPartType == o2::dataformats::kStrkV0) && !(strangeTrack.mPartType == o2::dataformats::kStrkCascade)) continue; // only fill if track is either V0 or cascade
        Pt = sqrt(strangeTrack.mDecayMom[0] * strangeTrack.mDecayMom[0] + strangeTrack.mDecayMom[1] * strangeTrack.mDecayMom[1]);
        XKF = strangeTrack.decayVtxXKF;
        YKF = strangeTrack.decayVtxYKF;
        ZKF = strangeTrack.decayVtxZKF;
        MassKF = strangeTrack.mMassKF;
        PtKF = strangeTrack.mPtKF;
        TopoChi2KF = strangeTrack.mTopoChi2KF;
        GeoChi2KF = strangeTrack.mGeoChi2KF;
        Tree->Fill();
    }


    OutputFile.cd();
    Tree->Write();
}