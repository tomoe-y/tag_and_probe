#include <iostream>
#include <vector>
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TEfficiency.h"

void mu_pair_select(){
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
    chain->SetBranchStatus("trigger_info_ptVec", 1);

    std::vector<float> *mu_m = 0;
    std::vector<float> *mu_pt = 0;
    std::vector<float> *mu_eta = 0;
    std::vector<float> *mu_phi = 0;
    std::vector<int> *mu_charge = 0;
    std::vector<int> *mu_author = 0;
    std::vector<int> *mu_muonType = 0;
    std::vector<std::string> *trigger_info_chain = 0;
    std::vector<int> *trigger_info_isPassed = 0;
    std::vector<std::vector<float>> *trigger_info_etaVec = 0;
    std::vector<std::vector<float>> *trigger_info_phiVec = 0;
    std::vector<std::vector<float>> *trigger_info_ptVec = 0;

    chain->SetBranchAddress("mu_m", &mu_m);
    chain->SetBranchAddress("mu_pt", &mu_pt);
    chain->SetBranchAddress("mu_eta", &mu_eta);
    chain->SetBranchAddress("mu_phi", &mu_phi);
    chain->SetBranchAddress("mu_charge", &mu_charge);
    chain->SetBranchAddress("mu_author", &mu_author);
    chain->SetBranchAddress("mu_muonType", &mu_muonType);
    chain->SetBranchAddress("trigger_info_chain", &trigger_info_chain);
    chain->SetBranchAddress("trigger_info_isPassed", &trigger_info_isPassed);
    chain->SetBranchAddress("trigger_info_phiVec", &trigger_info_phiVec);
    chain->SetBranchAddress("trigger_info_etaVec", &trigger_info_etaVec);
    chain->SetBranchAddress("trigger_info_ptVec", &trigger_info_ptVec);

    TCanvas *canvas1 = new TCanvas();
    TCanvas *canvas2 = new TCanvas();
    TCanvas *canvas3 = new TCanvas();
    TCanvas *canvas4 = new TCanvas();
    TCanvas *canvas5 = new TCanvas();
    TCanvas *canvas6 = new TCanvas();
    TCanvas *canvas7 = new TCanvas();
    TCanvas *canvas8 = new TCanvas();
    TCanvas *canvas9 = new TCanvas();

    TH1D *mass_hist = new TH1D("mass_hist", "mass_hist", 1000, 0, 100000);
    TH1D *pt_hist = new TH1D("pt_hist", "pt_hist", 1000, 0, 100000);
    TH1D *deltaR_hist = new TH1D("deltaR_hist", "delatR_hist", 100, 0, 5);
    TH1D *hlt_deltaR_hist = new TH1D("hlt_deltaR_hist", "hlt_deltaR_hist", 100, 0, 5);
    TH2D *hlt_deltaR_pt_hist = new TH2D("hlt_deltaR_pt_hist", "hlt_deltaR_pt_hist", 100, 0, 100000, 100, 0, 0.01);
    TH1D *tag_muon_momentum_hist = new TH1D("tag_muon_momentum_hist", "tag_muon_momentum_hist", 100, 0, 100000);
    TH1D *tag_muon_momentum_cut_hist = new TH1D("tag_muon_momentum_cut_hist", "tag_muon_momentum_cut_hist", 100, 0, 100000);
    TH1D *probe_muon_momentum_hist = new TH1D("probe_muon_momentum_hist", "probe_muon_momentum_hist", 100, 0, 100000);
    TH1D *probe_muon_momentum_cut_hist = new TH1D("probe_muon_momentum_cut_hist", "probe_muon_momentum_cut_hist", 100, 0, 100000);
    TH1D *probe_muon_endcap_momentum_hist = new TH1D("probe_muon_endcap_momentum_hist", "probe_muon_endcap_momentum_hist", 100, 0, 100000);
    TH1D *probe_muon_endcap_momentum_cut_hist = new TH1D("probe_muon_endcap_momentum_cut_hist", "probe_muon_endcap_momentum_cut_hist", 100, 0, 100000);
    TH1D *probe_muon_barrel_momentum_hist = new TH1D("probe_muon_barrel_momentum_hist", "probe_muon_barrel_momentum_hist", 100, 0, 100000);
    TH1D *probe_muon_barrel_momentum_cut_hist = new TH1D("probe_muon_barrel_momentum_cut_hist", "probe_muon_barrel_momentum_cut_hist", 100, 0, 100000);
    TH1D *probe_hlt_efficiency_hist = new TH1D("probe_hlt_efficiency_hist", "probe_hlt_efficiency_hist", 100, 0, 100000);
    TH1D *probe_endcap_hlt_efficiency_hist = new TH1D("probe_endcap_hlt_efficiency_hist", "probe_endcap_hlt_efficiency_hist", 100, 0, 100000);
    TH1D *probe_barrel_hlt_efficiency_hist = new TH1D("probe_barrel_hlt_efficiency_hist", "probe_barrel_hlt_efficiency_hist", 100, 0, 100000);

    //int entry = 10;
    int entry = chain->GetEntries();
    std::cout << entry << std::endl;

    std::cout << "start selecting muon..." << std::endl;

    //std::vector<std::vector<std::pair<int, int>>> tag_and_probe_mu_pair;

    int counts = 0;

    for(int i = 0; i < entry; i++){
        std::vector<std::pair<int, int>> mu_pair_number;
        chain->GetEntry(i);
        int trig_chain = 0;
        bool flag_mu_select = 0;
        for(int j = 0; j < trigger_info_chain->size(); j++){
            if(trigger_info_chain->at(j) == "HLT_mu26_ivarmedium"){
                if(trigger_info_isPassed->at(j) == 1){
                    trig_chain = j;
                    flag_mu_select = 1;
                }
            }
        }
        if(flag_mu_select){
            for(int j = 0; j < mu_m->size(); j++){
                //mu select
                bool flag_mu1_author = mu_author->at(j) == 1;
                bool flag_mu1_Type = mu_muonType->at(j) == 0;

                if(flag_mu1_author && flag_mu1_author){
                    TLorentzVector mu1;
                    mu1.SetPtEtaPhiM(mu_pt->at(j), mu_eta->at(j), mu_phi->at(j), mu_m->at(j));
                    for(int k = 0; k < mu_m->size(); k++){
                        if (j == k) continue;
                        //mu select & charge cut
                        int charge = mu_charge->at(j) * mu_charge->at(k);
                        bool flag_mu2_author = mu_author->at(k) == 1;
                        bool flag_mu2_Type = mu_muonType->at(k) == 0;
                        bool flag_charge = charge == -1;
                        
                        if (flag_mu2_author && flag_mu2_Type && flag_charge){
                            TLorentzVector mu2;
                            mu2.SetPtEtaPhiM(mu_pt->at(k), mu_eta->at(k), mu_phi->at(k), mu_m->at(k));

                            TLorentzVector mu_pair = mu1 + mu2;
                            float pair_mass = mu_pair.M();
                            float pair_deltaR = mu2.DeltaR(mu1);

                            // deltaR cut & mass cut
                            bool flag_pair_DeltaR = pair_deltaR > 0.4;
                            bool flag_pair_mass = pair_mass > 80000 && 100000 > pair_mass;

                            if (flag_pair_DeltaR && flag_pair_mass){
                                mu_pair_number.push_back(std::make_pair(j, k));
                                mass_hist->Fill(pair_mass);
                                pt_hist->Fill(mu_pt->at(k));
                            }
                        }
                    }
                }
            }
        }
        //tag
        std::cout << "start tagging" << std::endl;

        for (int j = 0; j < mu_pair_number.size(); j++){
            std::cout << "good" << std::endl;
            int tag_muon_number = mu_pair_number.at(j).first;
            int probe_muon_number = mu_pair_number.at(j).second;

            for (int k = 0; k < trigger_info_ptVec->at(trig_chain).size(); k++){
                std::cout << "good!" << std::endl;
                TVector3 tag_muon;
                tag_muon.SetPtEtaPhi(mu_pt->at(tag_muon_number), mu_eta->at(tag_muon_number), mu_phi->at(tag_muon_number));
                TVector3 hlt_muon;
                hlt_muon.SetPtEtaPhi(trigger_info_ptVec->at(trig_chain).at(k), trigger_info_etaVec->at(trig_chain).at(k), trigger_info_phiVec->at(trig_chain).at(k));
                float tag_DeltaR = tag_muon.DeltaR(hlt_muon);
                hlt_deltaR_hist->Fill(tag_DeltaR);
                hlt_deltaR_pt_hist->Fill(mu_pt->at(tag_muon_number), tag_DeltaR);
                tag_muon_momentum_hist->Fill(mu_pt->at(tag_muon_number));

                if(tag_DeltaR < 0.01){
                    tag_muon_momentum_cut_hist->Fill(mu_pt->at(tag_muon_number));
                    probe_muon_momentum_hist->Fill(mu_pt->at(probe_muon_number));
                    if(mu_pt->at(probe_muon_number) > 1.05 && mu_pt->at(probe_muon_number) < 2.4){
                        probe_muon_barrel_momentum_hist->Fill(mu_pt->at(probe_muon_number));
                    }
                    else{
                        probe_muon_endcap_momentum_hist->Fill(mu_pt->at(probe_muon_number));
                    }
                    //probe
                    for (int m = 0; m < trigger_info_ptVec->at(trig_chain).size(); m++){
                        TVector3 probe_muon;
                        probe_muon.SetPtEtaPhi(mu_pt->at(probe_muon_number), mu_eta->at(probe_muon_number), mu_phi->at(probe_muon_number));
                        TVector3 hlt_probe_muon;
                        hlt_probe_muon.SetPtEtaPhi(trigger_info_ptVec->at(trig_chain).at(m), trigger_info_etaVec->at(trig_chain).at(m), trigger_info_phiVec->at(trig_chain).at(m));
                        float probe_DeltaR = probe_muon.DeltaR(hlt_probe_muon);
                        if(probe_DeltaR < 0.01){
                            probe_muon_momentum_cut_hist->Fill(mu_pt->at(probe_muon_number));
                            if(mu_pt->at(probe_muon_number) > 1.05 && mu_pt->at(probe_muon_number) < 2.4){
                                probe_muon_barrel_momentum_cut_hist->Fill(mu_pt->at(probe_muon_number));
                            }
                            else{
                                probe_muon_endcap_momentum_cut_hist->Fill(mu_pt->at(probe_muon_number));
                            }
                            break;
                        }
                    }
                }
            }
        }
    }

    TEfficiency *pEff = new TEfficiency(*probe_muon_momentum_cut_hist, *probe_muon_momentum_hist);
    pEff->SetTitle("HLT efficiency;pt[MeV];efficiency");

    TEfficiency *pEff_endcap = new TEfficiency(*probe_muon_endcap_momentum_cut_hist, *probe_muon_endcap_momentum_hist);
    pEff_endcap->SetTitle("HLT efficiency in endcap;pt[MeV];efficiency");

    TEfficiency *pEff_barrel = new TEfficiency(*probe_muon_barrel_momentum_cut_hist, *probe_muon_barrel_momentum_hist);
    pEff_barrel->SetTitle("HLT efficiency in barrel;pt[MeV];efficiency");

    canvas1->cd();
    mass_hist->Draw();
    canvas2->cd();
    pt_hist->Draw();
    canvas3->cd();
    probe_muon_momentum_hist->Draw();
    probe_muon_momentum_cut_hist->Draw("same");
    canvas3->SaveAs("img/efficiency_hist.png");
    canvas4->cd();
    probe_hlt_efficiency_hist->Divide(probe_muon_momentum_cut_hist, probe_muon_momentum_hist);
    probe_hlt_efficiency_hist->Draw();
    canvas5->cd();
    pEff->Draw("AP");
    canvas5->SaveAs("img2/efficiency.png");
    canvas6->cd();
    pEff_endcap->Draw("AP");
    canvas6->SaveAs("img2/endcap_efficiency.png");
    canvas7->cd();
    pEff_barrel->Draw("AP");
    canvas7->SaveAs("img2/barrel_efficiency.png");
    canvas8->cd();
    probe_muon_endcap_momentum_cut_hist->Draw();
    probe_muon_endcap_momentum_hist->Draw("same");
    canvas8->SaveAs("img2/endcap_efficiency_hist.png");
    canvas9->cd();
    probe_muon_barrel_momentum_cut_hist->Draw("same");
    probe_muon_barrel_momentum_hist->Draw("same");
    canvas9->SaveAs("img2/barrel_efficiency_hist.png");

}