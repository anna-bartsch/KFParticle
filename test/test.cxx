/* Copyright (C) 2021-2021 GSI Helmholtzzentrum fuer Schwerionenforschung, Darmstadt
   SPDX-License-Identifier: GPL-3.0-only
   Authors: Sergey Gorbunov [committer] */

/// @file test.cxx
/// @author Sergey Gorbunov
/// @date 2021-08-06

#include "TApplication.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TROOT.h"
#include "Riostream.h"
#include "TDatabasePDG.h"
#include "TRandom.h"
#include "TMath.h"
#include "TFile.h"

#include "KFParticle.h"

using namespace std;

int main(int argc, char** /*argv*/)
{
  gROOT->SetBatch((argc <= 1) ? kTRUE : kFALSE);
  TApplication app("app", 0, nullptr);
  TCanvas* c = new TCanvas();
  TH1F* h = new TH1F("h", "h", 100, 0., 1.);
  TH2F* g = new TH2F("g", "g", 100, 0., 1., 100, 0., 1.);
  //g->Fill(0.5);
  g->Draw();
  h->Fill(0.5);
  h->Draw();
  c->Draw();
  std::cout << "batch mode " << c->IsBatch() << std::endl;

#ifdef HomogeneousField //ALICE
  std::cout << "HomogeneousField option is set" << std::endl;
#endif
#ifdef NonhomogeneousField //cbm
  std::cout << "NonhomogeneousField option is set" << std::endl;
#endif

  // check the field
  KFParticle::SetField(0.);
  {
     float xyz[3] = {0,0,0};
     float B[3];
     KFParticle p;
     p.GetFieldValue(xyz,B);
     std::cout<<"Field is set to " <<B[0]<<" "<<B[1]<<" "<<B[2]<<std::endl;
  }

  gRandom->SetSeed(1);
 
  TFile file("output.root","RECREATE");

  TH1F *x_res = new TH1F("x_res", "X resolution", 1000, -1, 1);   // x, y, z 
  TH1F *x_pull = new TH1F("x_pull", "X pull", 1000, -10, 10);
  TH1F *y_res = new TH1F("y_res", "Y resolution", 1000, -1, 1);
  TH1F *y_pull = new TH1F("y_pull", "Y pull", 1000, -10, 10);
  TH1F *z_res = new TH1F("z_res", "Z resolution", 1000, -0.5, 0.5);
  TH1F *z_pull = new TH1F("z_pull", "Z pull", 1000, -10, 10);

  TH2F *x_res_angle = new TH2F( "x_res_angle", "X Res and angle", 1000, -1, 1, 1000, 0, TMath::TwoPi()  );

  TH1F *Px_res = new TH1F("Px_res", "PX resolution", 1000, -0.1, 0.1);  // Px, Py, Pz
  TH1F *Px_pull = new TH1F("Px_pull", "PX pull", 1000, -10, 10);
  TH1F *Py_res = new TH1F("Py_res", "PY resolution", 1000, -0.1, 0.1);
  TH1F *Py_pull = new TH1F("Py_pull", "PY pull", 1000, -10, 10);
  TH1F *Pz_res = new TH1F("Pz_res", "PZ resolution", 1000, -0.1, 0.1);
  TH1F *Pz_pull = new TH1F("Pz_pull", "PZ pull", 1000, -10, 10);

  TH1F *E_res = new TH1F("e_res", "E resolution", 1000, -0.1, 0.1);  // E
  TH1F *E_pull = new TH1F("e_pull", "E pull", 1000, -10, 10);

  TH1F *Chi2 = new TH1F("chi2/ndf", "Chi2/NDF", 1000, -8, 8);




  for( int iter=0; iter<1000000; iter++) {
    double errors[6] = { .1, .1, .1, .01, .01, .01}; //x, y, z in mm und p_X, p_y, p_z in GeV
 
    double cov[21]; //Kovarianz Matrix
    for( int i=0; i<21; i++ ){ cov[i] = 0;}
    cov[ 0] = errors[0] * errors[0];
    cov[ 2] = errors[1] * errors[1];
    cov[ 5] = errors[2] * errors[2];
    cov[ 9] = errors[3] * errors[3];
    cov[14] = errors[4] * errors[4];
    cov[20] = errors[5] * errors[5];

    double phi1 = gRandom->Uniform(0,TMath::TwoPi());
    double phi2 = gRandom->Uniform(0,TMath::TwoPi());

    double paramReal1[6] = {0,0,0, TMath::Cos(phi1), TMath::Sin(phi1), 0};    
    double paramReal2[6] = {0,0,0, TMath::Cos(phi2), TMath::Sin(phi2), 0};
    
    double paramFit1[6]; // with errors
    double paramFit2[6];
  
    for( int i=0; i<6; i++){
      paramFit1[i] = paramReal1[i] + gRandom->Gaus(0,errors[i]);
      paramFit2[i] = paramReal2[i] + gRandom->Gaus(0,errors[i]);
    }

    KFParticle p1;
    p1.Create( paramFit1, cov, 1, (float) TDatabasePDG::Instance()->GetParticle("pi+")->Mass(), 0, 1);

    KFParticle p2;
    p2.Create( paramFit2, cov, -1, (float) TDatabasePDG::Instance()->GetParticle("pi-")->Mass(), 0, 1);

    KFParticle mother( p1, p2);

    KFParticle pReal1;
    pReal1.Create( paramReal1, cov, 1, (float) TDatabasePDG::Instance()->GetParticle("pi+")->Mass(), 0, 1 );

    KFParticle pReal2;
    pReal2.Create( paramReal2, cov, -1, (float) TDatabasePDG::Instance()->GetParticle("pi-")->Mass(), 0, 1 );

    KFParticle motherReal(pReal1, pReal2);
    

    //ein echtes Mutterteilchen erstellen, damit die Impulse und Energien verglichen werden können? xyz sind real ja null, aber Px und Py nicht.



    //std::cout<<" x " << mother.X()<<" y "<< mother.Y()<<" z "<<mother.Z()<<std::endl; //x, y, z
    x_res->Fill(mother.X() - 0.);
    x_pull->Fill( (mother.X() - 0.)/mother.GetErrX() );
    y_res->Fill(mother.Y() - 0.);
    y_pull->Fill((mother.Y() - 0.) / mother.GetErrY() );
    z_res->Fill(mother.Z() - 0.);
    z_pull->Fill( (mother.Z() - 0.)/mother.GetErrZ() );

    x_res_angle->Fill( (mother.X() -0.), abs( phi2 - phi1) );

    //std::cout<<" Px " << mother.Px()<<" Py "<< mother.Py()<<" Pz "<<mother.Pz()<<std::endl;  //Px, Py, Pz
    Px_res->Fill(mother.Px() - motherReal.Px());
    Px_pull->Fill( (mother.Px() - motherReal.Px())/mother.GetErrPx() );
    Py_res->Fill(mother.Py() - motherReal.Py());
    Py_pull->Fill((mother.Py() - motherReal.Py()) / mother.GetErrPy() );
    Pz_res->Fill(mother.Pz() - motherReal.Pz());
    Pz_pull->Fill( (mother.Pz() - motherReal.Pz())/mother.GetErrPz() );

    E_res->Fill(mother.E() - motherReal.E());
    E_pull->Fill( (mother.E() - motherReal.E())/mother.GetErrE() );

    Chi2->Fill( mother.GetChi2()  ); // mother.GetNDF()

    std::cout << mother.GetNDF() << std::endl;

  }
   file.cd();
  x_res->Write();
  x_pull->Write();
  y_res->Write();
  y_pull->Write();
  z_res->Write();
  z_pull->Write();

  x_res_angle->Write();

  Px_res->Write();
  Px_pull->Write();
  Py_res->Write();
  Py_pull->Write();
  Pz_res->Write();
  Pz_pull->Write();

  E_res->Write();
  E_pull->Write();

  Chi2->Write();

  file.Close();
  return 0;
}