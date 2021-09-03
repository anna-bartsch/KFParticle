/* Copyright (C) 2021-2021 GSI Helmholtzzentrum fuer Schwerionenforschung, Darmstadt
   SPDX-License-Identifier: GPL-3.0-only
   Authors: Sergey Gorbunov [committer] */

/// @file test.cxx
/// @author Sergey Gorbunov
/// @date 2021-08-06

#include "TApplication.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TROOT.h"
#include "Riostream.h"

#include "KFParticle.h"

int main(int argc, char** /*argv*/)
{
  gROOT->SetBatch((argc <= 1) ? kTRUE : kFALSE);
  TApplication app("app", 0, nullptr);
  TCanvas* c = new TCanvas();
  TH1F* h = new TH1F("h", "h", 100, 0., 1.);
  h->Fill(0.5);
  h->Draw();
  c->Draw();
  std::cout << "batch mode " << c->IsBatch() << std::endl;

#ifdef HomogeneousField
  std::cout << "HomogeneousField option is set" << std::endl;
#endif
#ifdef NonhomogeneousField
  std::cout << "NonhomogeneousField option is set" << std::endl;
#endif

  KFParticle p;
  p.Initialize();
  std::cout << p.Chi2() << std::endl;

  return 0;
}