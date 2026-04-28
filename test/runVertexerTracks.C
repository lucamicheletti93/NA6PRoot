#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TParticle.h>
#include "MagneticField.h"
#include "NA6PVertex.h"
#include "NA6PVertexerTracks.h"
#include "NA6PLayoutParam.h"
#include "NA6PVerTelCluster.h"
#include "NA6PVerTelReconstruction.h"
#endif

void runVertexerTracks(
                       const char* dirSimu = "/data/lmichele/datasets/real_config/PYTHIA_PbPb_MB_10000",
                       const char* na6pLayoutFile = "/data/lmichele/datasets/na6pLayout_real.ini",
                       const char* fOutName = "outputs/real_config/PYTHIA_PbPb_MB_newTest/vtx_from_tracks.root")
{
  auto magField = new MagneticField();
  magField->loadField();
  magField->setAsGlobalField();

  TFile* ft = new TFile(Form("%s/TracksVerTel.root", dirSimu));
  if (!ft)
    return;
  TTree* trTree = (TTree*)ft->Get("tracksVerTel");
  std::vector<NA6PTrack>* trArr = nullptr;
  trTree->SetBranchAddress("VerTel", &trArr);

  TFile* fk = new TFile(Form("%s/MCKine.root", dirSimu));
  TTree* mcTree = (TTree*)fk->Get("mckine");
  std::vector<TParticle>* mcArr = nullptr;
  mcTree->SetBranchAddress("tracks", &mcArr);

  int nVtx, targId;
  double genX, genY, genZ;
  int nContrib[20];
  double recsX[20], recsY[20], recsZ[20];

  TFile *fOut = new TFile(fOutName, "RECREATE");
  TTree *tree = new TTree("vertex", "vertex");
  tree -> Branch("nVtx", &nVtx, "nVtx/I");
  tree -> Branch("targId", &targId, "targId/I");
  tree -> Branch("genX", &genX, "genX/D");
  tree -> Branch("genY", &genY, "genY/D");
  tree -> Branch("genZ", &genZ, "genZ/D");
  tree -> Branch("nContrib", nContrib, "nContrib[20]/I");
  tree -> Branch("recsX", recsX, "recsX[20]/D");
  tree -> Branch("recsY", recsY, "recsY[20]/D");
  tree -> Branch("recsZ", recsZ, "recsZ[20]/D");

  TH1F* hncontr = new TH1F("hncontr", ";N_{contributors}", 100, -0.5, 999.5);
  TH1F* hnvert = new TH1F("hnvert", ";Number of reco vertices", 11, -0.5, 10.5);
  TH2F* hxrecgen = new TH2F("hxrecgen", ";x_{gen} (cm);x_{rec} (cm)", 50, -0.05, 0.05, 50, -0.05, 0.05);
  TH2F* hyrecgen = new TH2F("hyrecgen", ";y_{gen} (cm);y_{rec} (cm)", 50, -0.05, 0.05, 50, -0.05, 0.05);
  TH2F* hzrecgen = new TH2F("hzrecgen", ";z_{gen} (cm);z_{rec} (cm)", 50, -4., 2., 50, -4., 2.);
  TH1F* hdx = new TH1F("hdx", ";x_{rec} - x_{gen} (cm)", 100, -0.05, 0.05);
  TH1F* hdy = new TH1F("hdy", ";y_{rec} - y_{gen} (cm)", 100, -0.05, 0.05);
  TH1F* hdz = new TH1F("hdz", ";z_{rec} - z_{gen} (cm)", 100, -0.5, 0.5);

  int indexTarg = -999;

  NA6PVertexerTracks* vertxr = new NA6PVertexerTracks();
  vertxr->setVerbosity(true);

  int nEv = trTree->GetEntries();
  printf("Number of events = %d\n", nEv);
  for (int jEv = 0; jEv < nEv; jEv++) {
    mcTree->GetEvent(jEv);
    trTree->GetEvent(jEv);
    int nPart = mcArr->size();
    int nTracks = trArr->size();
    printf("Event %d particles = %d tracks = %d\n", jEv, nPart, nTracks);

    double xVertGen = 0;
    double yVertGen = 0;
    double zVertGen = 0;
    // get primary vertex position from the Kine Tree
    for (int jp = 0; jp < nPart; jp++) {
      auto curPart = mcArr->at(jp);
      if (curPart.IsPrimary()) {
        xVertGen = curPart.Vx();
        yVertGen = curPart.Vy();
        zVertGen = curPart.Vz();
        genX = xVertGen;
        genY = yVertGen;
        genZ = zVertGen;
        break;
      }
    }


    na6p::conf::ConfigurableParam::updateFromFile(na6pLayoutFile, "", true);
    const auto& layoutPar = NA6PLayoutParam::Instance();
    int nTargs = int(layoutPar.nTargets);


    // Find the position of the primary target
    for (int iTarg = 0;iTarg < nTargs;iTarg++) {
      double minTargZ = layoutPar.posTargetZ[iTarg] - ((layoutPar.thicknessTarget[iTarg])/2.);
      double maxTargZ = layoutPar.posTargetZ[iTarg] + ((layoutPar.thicknessTarget[iTarg])/2.);
      if (zVertGen > minTargZ && zVertGen < maxTargZ) {
        indexTarg = iTarg;
      }
    }
    targId = indexTarg;

    // Check if there are other interactions
    std::vector<int> tmpTargs;
    for (int i = 0;i < 5;i++) {
      if (i != indexTarg) {
        tmpTargs.push_back(i);
      } 
    }

    for (int jp = 0; jp < nPart; jp++) {
      auto curPart = mcArr->at(jp);
      double zGenTmp = curPart.Vz();
      if (!curPart.IsPrimary() && zGenTmp < 1 && zGenTmp > -9) {
        for (int i = tmpTargs.size() - 1; i >= 0; --i) {
          int tmpTarg = tmpTargs[i];
          double minTargZ = layoutPar.posTargetZ[tmpTarg] - ((layoutPar.thicknessTarget[tmpTarg])/2.);
          double maxTargZ = layoutPar.posTargetZ[tmpTarg] + ((layoutPar.thicknessTarget[tmpTarg])/2.);
          if (zGenTmp > minTargZ && zGenTmp < maxTargZ) {
            tmpTargs.erase(tmpTargs.begin() + i);
          }
        }
      }
    }


    vertxr->setBeamX(xVertGen);
    vertxr->setBeamY(yVertGen);
    vertxr->createTracksPool(*trArr);
    std::vector<NA6PVertex> vertices;
    vertxr->findVertices(vertices);
    int nVertices = vertices.size();
    nVtx = nVertices;
    hnvert->Fill(nVertices);
    int jv = 0;
    for (auto vert : vertices) {
      double xRec = vert.getX();
      double yRec = vert.getY();
      double zRec = vert.getZ();
      if (jv == 0) {
        hxrecgen->Fill(xVertGen, xRec);
        hyrecgen->Fill(yVertGen, yRec);
        hzrecgen->Fill(zVertGen, zRec);
        hdx->Fill(xRec - xVertGen);
        hdy->Fill(yRec - yVertGen);
        hdz->Fill(zRec - zVertGen);
        hncontr->Fill(vert.getNContributors());
      }

      nContrib[jv] = vert.getNContributors();
      recsX[jv] = xRec;
      recsY[jv] = yRec;
      recsZ[jv] = zRec;

      printf("Vertex %d, z = %f contrib = %d\n", jv++, zRec, vert.getNContributors());
    }
    tree -> Fill();
  }

  TCanvas* cv = new TCanvas("cv", "", 1000, 500);
  cv->Divide(2, 1);
  cv->cd(1);
  hnvert->Draw();
  cv->cd(2);
  hncontr->Draw();

  TCanvas* coutp = new TCanvas("coutp", "", 1400, 800);
  coutp->Divide(3, 2);
  coutp->cd(1);
  hxrecgen->Draw("colz");
  coutp->cd(2);
  hyrecgen->Draw("colz");
  coutp->cd(3);
  hzrecgen->Draw("colz");
  coutp->cd(4);
  hdx->Draw();
  coutp->cd(5);
  hdy->Draw();
  coutp->cd(6);
  hdz->Draw();

  fOut -> cd();
  cv->Write();
  coutp->Write();
  tree -> Write();
  fOut->Close();
}
