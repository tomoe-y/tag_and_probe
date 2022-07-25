#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TEfficiency.h"

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

    int entry = chain->GetEntries();
    std::cout << entry << std::endl;


    TH1D *mass_hist = new TH1D("mass_hist", "mass_hist", 1000, 0, 100000);
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

    TCanvas *canvas1 = new TCanvas("canvas1", "canvas1");
    TCanvas *canvas2 = new TCanvas("canvas2", "canvas2");
    TCanvas *canvas3 = new TCanvas("canvas3", "canvas3");
    TCanvas *canvas4 = new TCanvas("canvas4", "canvas4");
    TCanvas *canvas5 = new TCanvas("canvas5", "canvas5");
    TCanvas *canvas6 = new TCanvas("canvas6", "canvas6");
    TCanvas *canvas7 = new TCanvas("canvas7", "canvas7");
    TCanvas *canvas8 = new TCanvas("canvas8", "canvas8");
    TCanvas *canvas9 = new TCanvas("canvas9", "canvas9");

    std::cout << "start selecting muon..." << std::endl;

    for(int i = 0; i < entry; i++){
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
        if(flag_mu_select == 0) continue;
        //std::cout << trig_chain << std::endl;

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
                if(pair_deltaR < 0.4) continue;
                
                //mass cut
                if(pair_mass < 80000 || 100000 < pair_mass) continue;

                //fill hist
                mass_hist->Fill(pair_mass);
                deltaR_hist->Fill(pair_deltaR);

                //j:tag
                for (int l = 0; l < trigger_info_ptVec->at(trig_chain).size(); l++){
                    TVector3 tag_muon;
                    tag_muon.SetPtEtaPhi(mu_pt->at(j), mu_eta->at(j), mu_phi->at(j));
                    TVector3 hlt_muon;
                    hlt_muon.SetPtEtaPhi(trigger_info_ptVec->at(trig_chain).at(l), trigger_info_etaVec->at(trig_chain).at(l), trigger_info_phiVec->at(trig_chain).at(l));
                    float tag_DeltaR = tag_muon.DeltaR(hlt_muon);
                    hlt_deltaR_hist->Fill(tag_DeltaR);
                    hlt_deltaR_pt_hist->Fill(mu_pt->at(j), tag_DeltaR);
                    tag_muon_momentum_hist->Fill(mu_pt->at(j));

                    if(tag_DeltaR < 0.01){
                        tag_muon_momentum_cut_hist->Fill(mu_pt->at(j));
                        probe_muon_momentum_hist->Fill(mu_pt->at(k));
                        if(mu_pt->at(k) < 1.05){
                            probe_muon_barrel_momentum_hist->Fill(mu_pt->at(k));
                        }
                        else{
                            probe_muon_endcap_momentum_hist->Fill(mu_pt->at(k));
                        }
                        //k:probe
                        for (int m = 0; m < trigger_info_ptVec->at(trig_chain).size(); m++){
                            TVector3 probe_muon;
                            probe_muon.SetPtEtaPhi(mu_pt->at(k), mu_eta->at(k), mu_phi->at(k));
                            TVector3 hlt_probe_muon;
                            hlt_probe_muon.SetPtEtaPhi(trigger_info_ptVec->at(trig_chain).at(m), trigger_info_etaVec->at(trig_chain).at(m), trigger_info_phiVec->at(trig_chain).at(m));
                            float probe_DeltaR = probe_muon.DeltaR(hlt_probe_muon);
                            if(probe_DeltaR < 0.01){
                                probe_muon_momentum_cut_hist->Fill(mu_pt->at(k));
                                if(mu_pt->at(k) < 1.05){
                                    probe_muon_barrel_momentum_cut_hist->Fill(mu_pt->at(k));
                                }
                                else{
                                    probe_muon_endcap_momentum_cut_hist->Fill(mu_pt->at(k));
                                }
                                break;
                            }
                            /*
                            //endcap efficiency
                            if(mu_eta->at(k) < 1.05){
                                probe_muon_endcap_momentum_hist->Fill(mu_pt->at(j));
                                if(probe_DeltaR < 0.01){
                                    probe_muon_endcap_momentum_cut_hist->Fill(mu_pt->at(j));
                                }
                            }
                            //barrel efficiency
                            else{
                                probe_muon_barrel_momentum_hist->Fill(mu_pt->at(j));
                                if(probe_DeltaR < 0.01){
                                    probe_muon_barrel_momentum_cut_hist->Fill(mu_pt->at(j));
                                }
                            }
                            */
                        }
                    }
                }

                //k:tag
                for (int l = 0; l < trigger_info_ptVec->at(trig_chain).size(); l++){
                    TVector3 tag_muon;
                    tag_muon.SetPtEtaPhi(mu_pt->at(k), mu_eta->at(k), mu_phi->at(k));
                    TVector3 hlt_muon;
                    hlt_muon.SetPtEtaPhi(trigger_info_ptVec->at(trig_chain).at(l), trigger_info_etaVec->at(trig_chain).at(l), trigger_info_phiVec->at(trig_chain).at(l));
                    float tag_DeltaR = tag_muon.DeltaR(hlt_muon);
                    hlt_deltaR_hist->Fill(tag_DeltaR);
                    hlt_deltaR_pt_hist->Fill(mu_pt->at(j),tag_DeltaR);
                    tag_muon_momentum_hist->Fill(mu_pt->at(k));
                    if(tag_DeltaR < 0.01){
                        tag_muon_momentum_cut_hist->Fill(mu_pt->at(k));
                        probe_muon_momentum_hist->Fill(mu_pt->at(j));
                        if(mu_eta->at(j) < 1.05){
                            probe_muon_barrel_momentum_hist->Fill(mu_pt->at(j));
                        }
                        else{
                            probe_muon_endcap_momentum_hist->Fill(mu_pt->at(j));
                        }
                        //j:probe
                        for (int m = 0; m < trigger_info_ptVec->at(trig_chain).size(); m++){
                            TVector3 probe_muon;
                            probe_muon.SetPtEtaPhi(mu_pt->at(j), mu_eta->at(j), mu_phi->at(j));
                            TVector3 hlt_probe_muon;
                            hlt_probe_muon.SetPtEtaPhi(trigger_info_ptVec->at(trig_chain).at(m), trigger_info_etaVec->at(trig_chain).at(m), trigger_info_phiVec->at(trig_chain).at(m));
                            float probe_DeltaR = probe_muon.DeltaR(hlt_probe_muon);
                            if(probe_DeltaR < 0.01){
                                probe_muon_momentum_cut_hist->Fill(mu_pt->at(j));
                                if(mu_eta->at(j) < 1.05){
                                    probe_muon_barrel_momentum_cut_hist->Fill(mu_pt->at(j));
                                }
                                else{
                                    probe_muon_endcap_momentum_cut_hist->Fill(mu_pt->at(j));
                                }
                                break;
                            }
                            /*
                            //endcap efficiency
                            if(mu_eta->at(j) < 1.05){
                                probe_muon_endcap_momentum_hist->Fill(mu_pt->at(j));
                                if(probe_DeltaR < 0.01){
                                    probe_muon_endcap_momentum_cut_hist->Fill(mu_pt->at(j));
                                }
                            }
                            //barrel efficiency
                            else{
                                probe_muon_barrel_momentum_hist->Fill(mu_pt->at(j));
                                if(probe_DeltaR < 0.01){
                                    probe_muon_barrel_momentum_cut_hist->Fill(mu_pt->at(j));
                                }
                            }
                            */
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

    //canvas1->cd();
    //mass_hist->Draw();
    //canvas2->cd();
    //deltaR_hist->Draw();
    canvas1->cd();
    hlt_deltaR_hist->Draw();
    canvas2->cd();
    tag_muon_momentum_cut_hist->Draw();
    canvas3->cd();
    probe_muon_momentum_hist->Draw();
    probe_muon_momentum_cut_hist->Draw("same");
    canvas3->SaveAs("img/efficiency_hist.png");
    canvas4->cd();
    probe_hlt_efficiency_hist->Divide(probe_muon_momentum_cut_hist, probe_muon_momentum_hist);
    probe_hlt_efficiency_hist->Draw();
    //canvas4->SaveAs("img/efficiency_hist.png");
    canvas5->cd();
    pEff->Draw("AP");
    canvas5->SaveAs("img/efficiency.png");
    canvas6->cd();
    pEff_endcap->Draw("AP");
    canvas6->SaveAs("img/endcap_efficiency.png");
    canvas7->cd();
    pEff_barrel->Draw("AP");
    canvas7->SaveAs("img/barrel_efficiency.png");
    canvas8->cd();
    probe_muon_endcap_momentum_cut_hist->Draw();
    probe_muon_endcap_momentum_hist->Draw("same");
    canvas8->SaveAs("img/endcap_efficiency_hist.png");
    canvas9->cd();
    probe_muon_barrel_momentum_cut_hist->Draw("same");
    probe_muon_barrel_momentum_hist->Draw("same");
    canvas9->SaveAs("img/barrel_efficiency_hist.png");
}