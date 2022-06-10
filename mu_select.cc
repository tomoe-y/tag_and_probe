#include "TChain.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLorentzVector.h"

void mu_select(){
    TChain* chain = new TChain("physics");
    chain->Add("/mnt/susy11/data04/atlas/data16_13TeV/periodL/mu_sample/user.junpei.00311481.physics_Main.merge.NTUP_MCP.1.f758_m1714.00-00-32_L1TGCNtuple_derivated.02-04-00.root");
    //chain->Add("/mnt/susy11/data04/atlas/data16_13TeV/periodL/mu_sample/*.root");

    chain->SetBranchStatus("*", 0);
    chain->SetBranchStatus("mu_m", 1);
    chain->SetBranchStatus("mu_pt", 1);
    chain->SetBranchStatus("mu_eta", 1);
    chain->SetBranchStatus("mu_phi", 1);
    chain->SetBranchStatus("mu_charge", 1);
    chain->SetBranchStatus("mu_author", 1);
    chain->SetBranchStatus("mu_muonType", 1);
    chain->SetBranchStatus("trigger_info_chain", 1);
    chain->SetBranchStatus("trigger_info_isPassed", 1);
    chain->SetBranchStatus("trigger_info_etaVec", 1);
    chain->SetBranchStatus("trigger_info_phiVec", 1);

    std::vector<float> *mu_m = 0;
    std::vector<float> *mu_pt = 0;
    std::vector<float> *mu_eta = 0;
    std::vector<float> *mu_phi = 0;
    std::vector<int> *mu_charge = 0;
    std::vector<int> *mu_author = 0;
    std::vector<int> *mu_muonType = 0;
    std::vector<std::string> *trigger_info_chain = 0;
    std::vector<int> *trigger_info_isPassed = 0;

    chain->SetBranchAddress("mu_m", &mu_m);
    chain->SetBranchAddress("mu_pt", &mu_pt);
    chain->SetBranchAddress("mu_eta", &mu_eta);
    chain->SetBranchAddress("mu_phi", &mu_phi);
    chain->SetBranchAddress("mu_charge", &mu_charge);
    chain->SetBranchAddress("mu_author", &mu_author);
    chain->SetBranchAddress("mu_muonType", &mu_muonType);
    chain->SetBranchAddress("trigger_info_chain", &trigger_info_chain);
    chain->SetBranchAddress("trigger_info_isPassed", &trigger_info_isPassed);

    int entry = chain->GetEntries();
    std::cout << entry << std::endl;


    TH1D *mass_hist = new TH1D("mass_hist", "mass_hist", 1000, 0, 100000);
    TH1D *deltaR_hist = new TH1D("deltaR_hist", "delatR_hist", 100, 0, 5);

    TCanvas *canvas1 = new TCanvas("canvas1", "canvas1");
    TCanvas *canvas2 = new TCanvas("canvas2", "canvas2");

    for(int i = 0; i < entry; i++){
        chain->GetEntry(i);
        int trig_pass = 0;
        for(int j = 0; j < trigger_info_chain->size(); j++){
            if(trigger_info_chain->at(j) == "HLT_mu26_ivarmedium"){
                if(trigger_info_isPassed->at(j) != 1) continue;
                trig_pass = trigger_info_isPassed->at(j);
            }
        }
        if(trig_pass != 1) continue;
        //std::cout << trig_pass << std::endl;

        for(int j = 0; j < mu_m->size(); j++){
            //mu select
            if (mu_author->at(j) != 1) continue;
            if (mu_muonType->at(j) != 0) continue;

            TLorentzVector mu1;
            mu1.SetPtEtaPhiM(mu_pt->at(j), mu_eta->at(j), mu_phi->at(j), mu_m->at(j));

            for(int k = j + 1; k < mu_m->size(); k++){
                //mu select
                if (mu_author->at(j) != 1) continue;
                if (mu_muonType->at(j) != 0) continue;

                TLorentzVector mu2;
                mu2.SetPtEtaPhiM(mu_pt->at(k), mu_eta->at(k), mu_phi->at(k), mu_m->at(k));

                //charge cut
                int charge = mu_charge->at(j) * mu_charge->at(k);
                if(charge != -1) continue;

                TLorentzVector mu_pair = mu1 + mu2;
                float pair_mass = mu_pair.M();
                float pair_deltaR = mu2.DeltaR(mu1);

                // deltaR cut
                //if(pair_deltaR < 0.4) continue;
                
                //mass cut
                //if(pair_mass < 80000 && 100000 < pair_mass) continue;

                //std::cout << "a" << std::endl; 

                mass_hist->Fill(pair_mass);
                deltaR_hist->Fill(pair_deltaR);

            }
        }
    }

    canvas1->cd();
    mass_hist->Draw();
    canvas2->cd();
    deltaR_hist->Draw();
}