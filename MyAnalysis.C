#define MyAnalysis_cxx
// The class definition in MyAnalysis.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("MyAnalysis.C")
// Root > T->Process("MyAnalysis.C","some options")
// Root > T->Process("MyAnalysis.C+")
//

#include "MyAnalysis.h"
#include <iostream>
#include <TH1F.h>
#include <TLatex.h>

using namespace std;

void MyAnalysis::BuildEvent() {
   
   Muons.clear();
   for (int i = 0; i < NMuon; ++i) {
      MyMuon muon(Muon_Px[i], Muon_Py[i], Muon_Pz[i], Muon_E[i]);
      muon.SetIsolation(Muon_Iso[i]);
      muon.SetCharge(Muon_Charge[i]);
      Muons.push_back(muon);
   }
   
   Electrons.clear();
   for (int i = 0; i < NElectron; ++i) {
      MyElectron electron(Electron_Px[i], Electron_Py[i], Electron_Pz[i], Electron_E[i]);
      electron.SetIsolation(Electron_Iso[i]);
      electron.SetCharge(Electron_Charge[i]);
      Electrons.push_back(electron);
   }
   
   Photons.clear();
   for (int i = 0; i < NPhoton; ++i) {
      MyPhoton photon(Photon_Px[i], Photon_Py[i], Photon_Pz[i], Photon_E[i]);
      photon.SetIsolation(Photon_Iso[i]);
      Photons.push_back(photon);
   }
   
   Jets.clear();
   for (int i = 0; i < NJet; ++i) {
      MyJet jet(Jet_Px[i], Jet_Py[i], Jet_Pz[i], Jet_E[i]);
      jet.SetBTagDiscriminator(Jet_btag[i]);
      jet.SetJetID(Jet_ID[i]);
      Jets.push_back(jet);
   }
   
   hadB.SetXYZM(MChadronicBottom_px, MChadronicBottom_py, MChadronicBottom_pz, 4.8);
   lepB.SetXYZM(MCleptonicBottom_px, MCleptonicBottom_py, MCleptonicBottom_pz, 4.8);
   hadWq.SetXYZM(MChadronicWDecayQuark_px, MChadronicWDecayQuark_py, MChadronicWDecayQuark_pz, 0.0);
   hadWqb.SetXYZM(MChadronicWDecayQuarkBar_px, MChadronicWDecayQuarkBar_py, MChadronicWDecayQuarkBar_pz, 0.0);
   lepWl.SetXYZM(MClepton_px, MClepton_py, MClepton_pz, 0.0);
   lepWn.SetXYZM(MCneutrino_px, MCneutrino_py, MCneutrino_pz, 0.0);
   met.SetXYZM(MET_px, MET_py, 0., 0.);
   
   EventWeight *= weight_factor;
   
}

void MyAnalysis::Begin(TTree * /*tree*/) {
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).
   
   TString option = GetOption();
   
}

void MyAnalysis::SlaveBegin(TTree * /*tree*/) {
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).
   
   TString option = GetOption();
   
   h_Mmumu = new TH1F("Mmumu", "Invariant di-muon mass", 60, 60, 120);
   h_Mmumu->SetXTitle("m_{#mu#mu}");
   h_Mmumu->Sumw2();
   histograms.push_back(h_Mmumu);
   histograms_MC.push_back(h_Mmumu);
   
   h_NMuon = new TH1F("NMuon", "Number of muons", 7, 0, 7);
   h_NMuon->SetXTitle("No. Muons");
   h_NMuon->Sumw2();
   histograms.push_back(h_NMuon);
   histograms_MC.push_back(h_NMuon);

   h_Melel = new TH1F("Melel", "Invariant di-electron mass", 60, 60, 120);
   h_Melel->SetXTitle("m_{ee}");
   h_Melel->Sumw2();
   histograms.push_back(h_Melel);
   histograms_MC.push_back(h_Melel);
   
   h_NElectron = new TH1F("NElectron", "Number of electrons", 7, 0, 7);
   h_NElectron->SetXTitle("No. Electron");
   h_NElectron->Sumw2();
   histograms.push_back(h_NElectron);
   histograms_MC.push_back(h_NElectron);

   h_Mjeje = new TH1F("Mjeje", "Invariant di-jet mass", 60, 60, 120);
   h_Mjeje->SetXTitle("m_{jeje}");
   h_Mjeje->Sumw2();
   histograms.push_back(h_Mjeje);
   histograms_MC.push_back(h_Mjeje);
   
   h_NJet = new TH1F("NJet", "Number of Jets", 7, 0, 7);
   h_NJet->SetXTitle("No. Jet");
   h_NJet->Sumw2();
   histograms.push_back(h_NJet);
   histograms_MC.push_back(h_NJet);

   h_NbJet = new TH1F("NbJet", "Number of bJets", 7, 0, 7);
   h_NbJet->SetXTitle("No. bJet");
   h_NbJet->Sumw2();
   histograms.push_back(h_NbJet);
   histograms_MC.push_back(h_NbJet);

   h_oneLeptonMjeje = new TH1F("oneLeptonMjeje", "Invariant oneLeptondi-jet mass", 60, 60, 120);
   h_oneLeptonMjeje->SetXTitle("m_{jeje}");
   h_oneLeptonMjeje->Sumw2();
   histograms.push_back(h_oneLeptonMjeje);
   histograms_MC.push_back(h_oneLeptonMjeje);

   h_oneLeptonNJet = new TH1F("oneLeptonNJet", "Number of oneLeptonJets", 7, 0, 7);
   h_oneLeptonNJet->SetXTitle("No. oneLeptonJet");
   h_oneLeptonNJet->Sumw2();
   histograms.push_back(h_oneLeptonNJet);
   histograms_MC.push_back(h_oneLeptonNJet);

   h_oneLeptonMbjeje = new TH1F("oneLeptonMbjeje", "Invariant oneLeptondi-bjet mass", 60, 60, 120);
   h_oneLeptonMbjeje->SetXTitle("m_{bjeje}");
   h_oneLeptonMbjeje->Sumw2();
   histograms.push_back(h_oneLeptonMbjeje);
   histograms_MC.push_back(h_oneLeptonMbjeje);

   h_oneLeptonNbJet = new TH1F("oneLeptonNbJet", "Number of oneLeptonbJets", 7, 0, 7);
   h_oneLeptonNbJet->SetXTitle("No. oneLeptonbJet");
   h_oneLeptonNbJet->Sumw2();
   histograms.push_back(h_oneLeptonNbJet);
   histograms_MC.push_back(h_oneLeptonNbJet);

	h_oneLeptonMjeje2 = new TH1F("oneLeptonMjeje2", "Invariant oneLeptondi-jet mass", 60, 60, 120);
   h_oneLeptonMjeje2->SetXTitle("m2_{jeje}");
   h_oneLeptonMjeje2->Sumw2();
   histograms.push_back(h_oneLeptonMjeje2);
   histograms_MC.push_back(h_oneLeptonMjeje2);

   h_oneLeptonNJet2 = new TH1F("oneLeptonNJet2", "Number of oneLeptonJets", 7, 0, 7);
   h_oneLeptonNJet2->SetXTitle("No. oneLeptonJet2");
   h_oneLeptonNJet2->Sumw2();
   histograms.push_back(h_oneLeptonNJet2);
   histograms_MC.push_back(h_oneLeptonNJet2);

   h_oneLeptonMbjeje2 = new TH1F("oneLeptonMbjeje2", "Invariant oneLeptondi-bjet mass", 60, 60, 120);
   h_oneLeptonMbjeje2->SetXTitle("m2_{bjeje}");
   h_oneLeptonMbjeje2->Sumw2();
   histograms.push_back(h_oneLeptonMbjeje2);
   histograms_MC.push_back(h_oneLeptonMbjeje2);

   h_oneLeptonNbJet2 = new TH1F("oneLeptonNbJet2", "Number of oneLeptonbJets", 7, 0, 7);
   h_oneLeptonNbJet2->SetXTitle("No. oneLeptonbJet2");
   h_oneLeptonNbJet2->Sumw2();
   histograms.push_back(h_oneLeptonNbJet2);
   histograms_MC.push_back(h_oneLeptonNbJet2);

	h_oneLeptonMjeje3 = new TH1F("oneLeptonMjeje3", "Invariant oneLeptondi-jet mass", 60, 60, 120);
   h_oneLeptonMjeje3->SetXTitle("m3_{jeje}");
   h_oneLeptonMjeje3->Sumw2();
   histograms.push_back(h_oneLeptonMjeje3);
   histograms_MC.push_back(h_oneLeptonMjeje3);

   h_oneLeptonNJet3 = new TH1F("oneLeptonNJet3", "Number of oneLeptonJets", 7, 0, 7);
   h_oneLeptonNJet3->SetXTitle("No. oneLeptonJet3");
   h_oneLeptonNJet3->Sumw2();
   histograms.push_back(h_oneLeptonNJet3);
   histograms_MC.push_back(h_oneLeptonNJet3);

   h_oneLeptonMbjeje3 = new TH1F("oneLeptonMbjeje3", "Invariant oneLeptondi-bjet mass", 60, 60, 120);
   h_oneLeptonMbjeje3->SetXTitle("m3_{bjeje}");
   h_oneLeptonMbjeje3->Sumw2();
   histograms.push_back(h_oneLeptonMbjeje3);
   histograms_MC.push_back(h_oneLeptonMbjeje3);

   h_oneLeptonNbJet3 = new TH1F("oneLeptonNbJet3", "Number of oneLeptonbJets", 7, 0, 7);
   h_oneLeptonNbJet3->SetXTitle("No. oneLeptonbJet3");
   h_oneLeptonNbJet3->Sumw2();
   histograms.push_back(h_oneLeptonNbJet3);
   histograms_MC.push_back(h_oneLeptonNbJet3);
   
   
}

Bool_t MyAnalysis::Process(Long64_t entry) {
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either MyAnalysis::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.
   
   ++TotalEvents;
   
   GetEntry(entry);
   
   if (TotalEvents % 10000 == 0)
      cout << "Next event -----> " << TotalEvents << endl;
   
   BuildEvent();
   
   double MuonPtCut = 25.;
   double MuonRelIsoCut = 0.10;

   double ElectronPtCut = 25.;
   double ElectronRelIsoCut = 0.10;

   double JetPtCut = 30.;
   double JetEtaCut = 2.5;

   double bTagCut = 2;
   
   //   cout << "Jets: " << endl;
   //   for (vector<MyJet>::iterator it = Jets.begin(); it != Jets.end(); ++it) {
   //      cout << "pt, eta, phi, btag, id: " << it->Pt() << ", " << it->Eta() << ", " << it->Phi() << ", " << it->IsBTagged() << ", " << it->GetJetID()
   //      << endl;
   //   }
   //   cout << "Muons: " << endl;
   //   for (vector<MyMuon>::iterator it = Muons.begin(); it != Muons.end(); ++it) {
   //      cout << "pt, eta, phi, iso, charge: " << it->Pt() << ", " << it->Eta() << ", " << it->Phi() << ", "
   //      << it->GetIsolation() << ", " << it->GetCharge() << endl;
   //   }
   //   cout << "Electrons: " << endl;
   //   for (vector<MyElectron>::iterator it = Electrons.begin(); it != Electrons.end(); ++it) {
   //      cout << "pt, eta, phi, iso, charge: " << it->Pt() << ", " << it->Eta() << ", " << it->Phi() << ", "
   //      << it->GetIsolation() << ", " << it->GetCharge() << endl;
   //   }
   //   cout << "Photons: " << endl;
   //   for (vector<MyPhoton>::iterator it = Photons.begin(); it != Photons.end(); ++it) {
   //      cout << "pt, eta, phi, iso: " << it->Pt() << ", " << it->Eta() << ", " << it->Phi() << ", " << it->GetIsolation()
   //      << endl;
   //   }
   
   
   //////////////////////////////
   // Exercise 1: Invariant Di-Muon mass
   
   int N_IsoMuon = 0;
   MyMuon *muon1, *muon2;
   
   for (vector<MyMuon>::iterator jt = Muons.begin(); jt != Muons.end(); ++jt) {
      if (jt->IsIsolated(MuonRelIsoCut)) {
         ++N_IsoMuon;
         if (N_IsoMuon == 1) muon1 = &(*jt);
         if (N_IsoMuon == 2) muon2 = &(*jt);
      }
   }
   
   h_NMuon->Fill(N_IsoMuon, EventWeight);
   
   if (N_IsoMuon > 1 && triggerIsoMu24) {
      if (muon1->Pt()>MuonPtCut) {
         h_Mmumu->Fill((*muon1 + *muon2).M(), EventWeight);
      }
   }

   int N_IsoElectron = 0;
   MyElectron *Electron1, *Electron2;
   
   for (vector<MyElectron>::iterator jt = Electrons.begin(); jt != Electrons.end(); ++jt) {
      if (jt->IsIsolated(ElectronRelIsoCut)) {
         ++N_IsoElectron;
         if (N_IsoElectron == 1) Electron1 = &(*jt);
         if (N_IsoElectron == 2) Electron2 = &(*jt);
              }
        }
   
   h_NElectron->Fill(N_IsoElectron, EventWeight);
   
   if (N_IsoElectron > 1 && triggerIsoMu24) {
      if (Electron1->Pt()>ElectronPtCut) {
         h_Melel->Fill((*Electron1 + *Electron2).M(), EventWeight);
              }
       }
   
   int N_lepton = N_IsoMuon + N_IsoElectron;

	int N_Jet = 0;
	int N_bJet = 0;
	MyJet *jet1, *jet2, *bjet1, *bjet2;
      
   for (vector<MyJet>::iterator jt = Jets.begin(); jt != Jets.end(); ++jt) {
		if (jt->Pt()>JetPtCut && abs(jt->Eta())<JetEtaCut){
			++N_Jet;
			if (N_Jet == 1) jet1 = &(*jt);
			if (N_Jet == 2) jet2 = &(*jt);
			if (jt->IsBTagged(bTagCut)) {
				++N_bJet;  //(bTag 개수)
				if (N_bJet == 1) bjet1 = &(*jt);
				if (N_bJet == 2) bjet2 = &(*jt);
			} 
		}
	}
	int N_bTagged = N_Jet + N_bJet;
   
	h_NJet->Fill(N_Jet, EventWeight);

    h_NbJet->Fill(N_bJet, EventWeight);

	h_Mjeje->Fill((*jet1 + *jet2).M(), EventWeight); //page 5(순서는 MyAnalysis.h에서 TH1F에 무슨 순서로 쓰느냐에 따라 정해짐)
	if(N_lepton != 1) return kTRUE; //kTRUE->루프 빠짐
	h_oneLeptonNJet->Fill(N_Jet, EventWeight);
	h_oneLeptonMjeje->Fill((*jet1 + *jet2).M(), EventWeight);

	h_oneLeptonNbJet->Fill(N_bJet, EventWeight);
	h_oneLeptonMbjeje->Fill((*bjet1 + *bjet2).M(), EventWeight);

	if(N_bJet < 1) return kTRUE; //이 상태는 N_Jet!=0의 컷만 줬을때의 ttbar개수가 온전히 남지않음
	h_oneLeptonNJet2->Fill(N_Jet, EventWeight);
	h_oneLeptonMjeje2->Fill((*jet1 + *jet2).M(), EventWeight);

	h_oneLeptonNbJet2->Fill(N_bJet, EventWeight);
	h_oneLeptonMbjeje2->Fill((*bjet1 + *bjet2).M(), EventWeight);



	if(N_bJet < 2) return kTRUE; //이 상태는 N_Jet!=0의 컷만 줬을때의 ttbar개수가 온전히 남지않음
	h_oneLeptonNJet3->Fill(N_Jet, EventWeight);
	h_oneLeptonMjeje3->Fill((*jet1 + *jet2).M(), EventWeight);

	h_oneLeptonNbJet3->Fill(N_bJet, EventWeight);
	h_oneLeptonMbjeje3->Fill((*bjet1 + *bjet2).M(), EventWeight);

	

	



   
   /*if (N_IsoJet > 1 && triggerIsoMu24) {
      if (jet1->Pt()>JetPtCut) {
         h_Mjeje->Fill((*jet1 + *jet2).M(), EventWeight);
      }
   }*/
   //////////////////////////////
   
   return kTRUE;
}

void MyAnalysis::SlaveTerminate() {
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
   
}

void MyAnalysis::Terminate() {
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
   
}
