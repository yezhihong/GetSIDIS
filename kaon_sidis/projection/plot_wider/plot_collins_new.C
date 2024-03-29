#include "plot_collins3.h"
void plot_new(){
  
	gStyle->SetLabelSize(0.05,"Y");
  gStyle->SetTitleSize(0.05,"Y");
  //gStyle->SetTitleOffset(1.2,"Y");
  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetTitleSize(0.05,"X");
	
  TCanvas *c1 = new TCanvas("c1","c1",1000,600);
  c1->SetFillColor(10);
  c1->cd(1);
/*void plot(particle_flag,Int_t Q2_flag =1, Int_t z_flag=1, Int_t flag_t=0, Int_t flag=0)*/
  //flag_t: 1->collins, 2->sivers
  //flag: 1->ploting Asym axis
  //plot(1,1,1,1,1);
  plot(1,2,3,1,1);
   
//  plot_feng(1,1,1);
  plot_e06010(1,1,3);
//  plot_ase(1,1);
//  plot_pretz(1,1);
//  plot_pretz(1,2);

//  TLegend *l1 = new TLegend(0.60,0.70,0.90,0.95);
  TLegend *l1 = new TLegend(0.5,0.81,0.89,0.94);
  Float_t x=0.,y=0.;
  TGraphErrors *g1 = new TGraphErrors(1,&x,&y,0,0);
  TGraphErrors *g2 = new TGraphErrors(1,&x,&y,0,0);
  TGraphErrors *g3 = new TGraphErrors(1,&x,&y,0,0);
  TGraphErrors *g4 = new TGraphErrors(1,&x,&y,0,0);
  TGraphErrors *g5 = new TGraphErrors(1,&x,&y,0,0);
  TGraphErrors *g6 = new TGraphErrors(1,&x,&y,0,0);
  g1->SetMarkerColor(1);
  g1->SetMarkerSize(1.5);  
  g1->SetMarkerStyle(21);
  g2->SetLineColor(3);
  g3->SetLineColor(6);
  g3->SetLineWidth(2.);
  g3->SetLineStyle(1);
  g4->SetLineColor(4);
  g4->SetLineWidth(2.5);
  g5->SetLineColor(7);
  g5->SetLineWidth(2.5);
  g6->SetMarkerColor(2);
  g6->SetMarkerStyle(20);
  g6->SetMarkerSize(1.2);
  l1->SetTextSize(0.04);
  l1->SetTextFont(22);
  l1->AddEntry(g1,"E06010, PRL 107 072003 (2011)","P");
  l1->Draw();
//  l1->AddEntry(g2,"Vogelsang and Yuan (Collins)","L");
//  l1->AddEntry(g3,"Anselmino et al.(Collins)","L");
//  l1->AddEntry(g4,"Pasquini (Collins)","L");
//  l1->AddEntry(g5,"Pasquini (Pretzelosity)","L");
  l1->AddEntry(g6,"90 days SoLID","P");
  l1->Draw();
  
  TCanvas *c2 = new TCanvas("c2","c2",1000,600);
  c2->SetFillColor(10);
  c2->cd(1);
  /*void plot(particle_flag,Int_t Q2_flag =1, Int_t z_flag=1, Int_t flag_t=0, Int_t flag=0)*/
  //flag_t: 1->collins, 2->sivers
  //flag: 1->ploting Asym axis
 // plot(2,1,1,1,1);
 plot(2,2,3,1,1);
   
//  plot_feng(2,1,1);
  plot_e06010(2,1,3);
//  plot_ase(2,1);
//  plot_pretz(2,1);
//  plot_pretz(2,2);
  l1->Draw();
/*
	c1->Print("SoLID_pip_collins.pdf");
	c1->Print("SoLID_pip_collins.png");
	c2->Print("SoLID_pim_collins.pdf");
	c2->Print("SoLID_pim_collins.png");
*/
	c1->Print("SoLID_pip_collins_NoTheo.pdf");
	c1->Print("SoLID_pip_collins_NoTheo.png");
	c2->Print("SoLID_pim_collins_NoTheo.pdf");
	c2->Print("SoLID_pim_collins_NoTheo.png");
}
