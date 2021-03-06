#include "TH1F.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TMath.h"

TH1F* rebin(TH1F* old) {


  old->SetBinContent(1, old->GetBinContent(1) + old->GetBinContent(0));
  old->SetBinContent(old->GetXaxis()->GetNbins(), old->GetBinContent(old->GetXaxis()->GetNbins()) + old->GetBinContent(old->GetXaxis()->GetNbins()+1));
  old->SetBinError(1, TMath::Sqrt(TMath::Power(old->GetBinError(1), 2) + TMath::Power(old->GetBinError(0), 2)));
  old->SetBinError(old->GetXaxis()->GetNbins(), TMath::Sqrt(TMath::Power(old->GetBinError(old->GetXaxis()->GetNbins()), 2) + TMath::Power(old->GetBinError(old->GetXaxis()->GetNbins()+1), 2)));
  old->SetBinContent(0, 0.);
  old->SetBinContent(old->GetXaxis()->GetNbins()+1, 0.);
  old->SetBinError(0, 0.);
  old->SetBinError(old->GetXaxis()->GetNbins()+1, 0.);

  string name = old->GetName();

  int nb = 0;
  float nbins[100];

  /*if (name.find("_rebin") != string::npos && (name.find("WGG_ele_pho0_pho1_pt") != string::npos || name.find("WGG_muo_pho0_pho1_pt") != string::npos ||
                                              name.find("ZGG_ele_pho0_pho1_pt") != string::npos || name.find("ZGG_muo_pho0_pho1_pt") != string::npos)) {

    nb = 4;
    nbins[0] = 2;
    nbins[1] = 2;
    nbins[2] = 1;
    nbins[3] = 1;

  } else if (name.find("_rebin") != string::npos && (name.find("WGG_ele_pho1_pt") != string::npos || name.find("WGG_muo_pho1_pt") != string::npos ||
                                                     name.find("ZGG_ele_pho1_pt") != string::npos || name.find("ZGG_muo_pho1_pt") != string::npos)) {
    return old;
  } else if (name.find("_rebin") != string::npos) {
    return old;*/

  if (name.find("EB") != string::npos) {
    std::cout << "test" << std::endl;
    nb = 6;
    nbins[0] = 1;
    nbins[1] = 1;
    nbins[2] = 1;
    nbins[3] = 1;
    nbins[4] = 2;
    nbins[5] = 4;

  } else if (name.find("EE") != string::npos) {

    nb = 4;
    nbins[0] = 1;
    nbins[1] = 1;
    nbins[2] = 2;
    nbins[3] = 6;

  /*} if (name.find("phoPt") != string::npos) {

    nb = 19;
    for (int i = 0; i < 15; i++) {
      nbins[i] = 1;
    }
    nbins[15] = 2;
    nbins[16] = 3;
    nbins[17] = 4;
    nbins[18] = 6;

  } if (name.find("lllPt") != string::npos) {

    nb = 19;
    for (int i = 0; i < 14; i++) {
      nbins[i] = 1;
    }
    nbins[14] = 1;
    nbins[15] = 1;
    nbins[16] = 2;
    nbins[17] = 3;
    nbins[18] = 7;*/

  } else {
    return old;
  }

  if (old->GetNbinsX() == nb) {
    cout << "ERROR: " << name << " already rebinned" << endl;
    return old;
  }

  int n = 0;
  for (int i = 0; i < nb; i++) {
    n = n + nbins[i];
  }
  if (n != old->GetNbinsX()) {
    cout << "ERROR: " << name << " - " << n << " " << old->GetNbinsX() << endl;
    return old;
  }

  int s = 0;
  const int nbs = nb + 1;
  double xbins[nbs];
  for (int i = 0; i < nb + 1; i++) {
    if (i > 0) s = s + nbins[i - 1];
    xbins[i] = old->GetXaxis()->GetBinUpEdge(s);
  }

  TH1F* tmp = new TH1F("tmp", old->GetTitle(), nb, xbins);

  for (int i = 0; i < old->GetXaxis()->GetNbins()+2; i++) {
    float c1 = old->GetBinContent(i);
    float e1 = old->GetBinError(i);
    float c1b = tmp->GetBinContent(tmp->GetXaxis()->FindBin(old->GetXaxis()->GetBinCenter(i)));
    float e1b = tmp->GetBinError(tmp->GetXaxis()->FindBin(old->GetXaxis()->GetBinCenter(i)));
    tmp->SetBinContent(tmp->GetXaxis()->FindBin(old->GetXaxis()->GetBinCenter(i)), c1+c1b);
    tmp->SetBinError(tmp->GetXaxis()->FindBin(old->GetXaxis()->GetBinCenter(i)), TMath::Sqrt(e1*e1+e1b*e1b));
  }

  old->Delete();
  tmp->SetName(name.c_str());

  return tmp;
}

/*TH2D* rebin(TH2D* old) {

  for (int i = 1; i < old->GetXaxis()->GetNbins()+1; i++) {
    old->SetBinContent(i, 1, old->GetBinContent(i, 1) + old->GetBinContent(i, 0));
    old->SetBinContent(i, old->GetYaxis()->GetNbins(), old->GetBinContent(i, old->GetYaxis()->GetNbins()) + old->GetBinContent(i, old->GetYaxis()->GetNbins()+1));
    old->SetBinError(i, 1, TMath::Sqrt(TMath::Power(old->GetBinContent(i, 1), 2) + TMath::Power(old->GetBinContent(i, 0), 2)));
    old->SetBinError(i, old->GetYaxis()->GetNbins(), TMath::Sqrt(TMath::Power(old->GetBinContent(i, old->GetYaxis()->GetNbins()), 2) + TMath::Power(old->GetBinContent(i, old->GetYaxis()->GetNbins()+1), 2)));
    old->SetBinContent(i, 0, 0.);
    old->SetBinContent(i, old->GetYaxis()->GetNbins()+1, 0.);
    old->SetBinError(i, 0, 0.);
    old->SetBinError(i, old->GetYaxis()->GetNbins()+1, 0.);
  }
  for (int i = 1; i < old->GetYaxis()->GetNbins()+1; i++) {
    old->SetBinContent(1, i, old->GetBinContent(1, i) + old->GetBinContent(0, i));
    old->SetBinContent(old->GetXaxis()->GetNbins(), i, old->GetBinContent(old->GetXaxis()->GetNbins(), i) + old->GetBinContent(old->GetXaxis()->GetNbins()+1, i));
    old->SetBinError(1, i, TMath::Sqrt(TMath::Power(old->GetBinContent(1, i), 2) + TMath::Power(old->GetBinContent(0, i), 2)));
    old->SetBinError(old->GetXaxis()->GetNbins(), i, TMath::Sqrt(TMath::Power(old->GetBinContent(old->GetXaxis()->GetNbins(), i), 2) + TMath::Power(old->GetBinContent(old->GetXaxis()->GetNbins()+1, i), 2)));
    old->SetBinContent(0, i, 0.);
    old->SetBinContent(old->GetXaxis()->GetNbins()+1, i, 0.);
    old->SetBinError(0, i, 0.);
    old->SetBinError(old->GetXaxis()->GetNbins()+1, i, 0.);
  }

  string name = old->GetName();

  int nxb = 0;
  float nxbins[100];
  int nyb = 0;
  float nybins[100];

  if (name.find("misid") != string::npos) {

    nxb = 8;
    for (int i = 0; i < 8; i++) {
      nxbins[i] = 1;
    }

    nyb = 19;
    for (int i = 0; i < 15; i++) {
      nybins[i] = 2;
    }
    nybins[15] = 6;
    nybins[16] = 6;
    nybins[17] = 6;
    nybins[18] = 12;

  } else {
    return old;
  }

  if (old->GetNbinsX() == nxb && old->GetNbinsY() == nyb) {
    cout << "ERROR: " << name << " already rebinned" << endl;
    return old;
  }

  int nx = 0;
  for (int i = 0; i < nxb; i++) {
    nx = nx + nxbins[i];
  }
  if (nx != old->GetNbinsX()) {
    cout << "ERROR: " << name << " - " << nx << " " << old->GetNbinsX() << endl;
    return old;
  }
  int ny = 0;
  for (int i = 0; i < nyb; i++) {
    ny = ny + nybins[i];
  }
  if (ny != old->GetNbinsY()) {
    cout << "ERROR: " << name << " - " << ny << " " << old->GetNbinsY() << endl;
    return old;
  }

  int sx = 0;
  const int nxbs = nxb + 1;
  double xbins[nxbs];
  for (int i = 0; i < nxb + 1; i++) {
    if (i > 0) sx = sx + nxbins[i-1];
    xbins[i] = old->GetXaxis()->GetBinUpEdge(sx);
  }
  int sy = 0;
  const int nybs = nyb + 1;
  double ybins[nybs];
  for (int i = 0; i < nyb + 1; i++) {
    if (i > 0) sy = sy + nybins[i-1];
    ybins[i] = old->GetYaxis()->GetBinUpEdge(sy);
  }

  TH2D* tmp = new TH2D("tmp", old->GetTitle(), nxb, xbins, nyb, ybins);

  for (int i = 0; i < old->GetXaxis()->GetNbins()+2; i++) {
    for (int j = 0; j < old->GetYaxis()->GetNbins()+2; j++) {
      float c2 = old->GetBinContent(i, j);
      float e2 = old->GetBinError(i, j);
      float c2b = tmp->GetBinContent(tmp->GetXaxis()->FindBin(old->GetXaxis()->GetBinCenter(i)), tmp->GetYaxis()->FindBin(old->GetYaxis()->GetBinCenter(j)));
      float e2b = tmp->GetBinError(tmp->GetXaxis()->FindBin(old->GetXaxis()->GetBinCenter(i)), tmp->GetYaxis()->FindBin(old->GetYaxis()->GetBinCenter(j)));
      tmp->SetBinContent(tmp->GetXaxis()->FindBin(old->GetXaxis()->GetBinCenter(i)), tmp->GetYaxis()->FindBin(old->GetYaxis()->GetBinCenter(j)), c2+c2b);
      tmp->SetBinError(tmp->GetXaxis()->FindBin(old->GetXaxis()->GetBinCenter(i)), tmp->GetYaxis()->FindBin(old->GetYaxis()->GetBinCenter(j)), TMath::Sqrt(e2*e2+e2b*e2b));
    }
  }

  old->Delete();
  tmp->SetName(name.c_str());

  return tmp;
}

TH3D* rebin(TH3D* old) {

  for (int i = 1; i < old->GetXaxis()->GetNbins()+1; i++) {
    old->SetBinContent(i, 1, 1, old->GetBinContent(i, 1, 1) + old->GetBinContent(i, 0, 1));
    old->SetBinContent(i, old->GetYaxis()->GetNbins(), 1, old->GetBinContent(i, old->GetYaxis()->GetNbins(), 1) + old->GetBinContent(i, old->GetYaxis()->GetNbins()+1, 1));
    old->SetBinError(i, 1, 1, TMath::Sqrt(TMath::Power(old->GetBinContent(i, 1, 1), 2) + TMath::Power(old->GetBinContent(i, 0, 1), 2)));
    old->SetBinError(i, old->GetYaxis()->GetNbins(), 1, TMath::Sqrt(TMath::Power(old->GetBinContent(i, old->GetYaxis()->GetNbins(), 1), 2) + TMath::Power(old->GetBinContent(i, old->GetYaxis()->GetNbins()+1, 1), 2)));
    old->SetBinContent(i, 0, 1, 0.);
    old->SetBinContent(i, old->GetYaxis()->GetNbins()+1, 1, 0.);
    old->SetBinError(i, 0, 1, 0.);
    old->SetBinError(i, old->GetYaxis()->GetNbins()+1, 1, 0.);

    old->SetBinContent(i, 1, 1, old->GetBinContent(i, 1, 1) + old->GetBinContent(i, 1, 0));
    old->SetBinContent(i, 1, old->GetZaxis()->GetNbins(), old->GetBinContent(i, 1, old->GetZaxis()->GetNbins()) + old->GetBinContent(i, 1, old->GetZaxis()->GetNbins()+1));
    old->SetBinError(i, 1, 1, TMath::Sqrt(TMath::Power(old->GetBinContent(i, 1, 1), 2) + TMath::Power(old->GetBinContent(i, 1, 0), 2)));
    old->SetBinError(i, 1, old->GetZaxis()->GetNbins(), TMath::Sqrt(TMath::Power(old->GetBinContent(i, 1, old->GetZaxis()->GetNbins()), 2) + TMath::Power(old->GetBinContent(i, 1, old->GetZaxis()->GetNbins()+1), 2)));
    old->SetBinContent(i, 1, 0, 0.);
    old->SetBinContent(i, 1, old->GetZaxis()->GetNbins()+1, 0.);
    old->SetBinError(i, 1, 0, 0.);
    old->SetBinError(i, 1, old->GetZaxis()->GetNbins()+1, 0.);
  }
  for (int i = 1; i < old->GetYaxis()->GetNbins()+1; i++) {
    old->SetBinContent(1, i, 1, old->GetBinContent(1, i, 1) + old->GetBinContent(0, i, 1));
    old->SetBinContent(old->GetXaxis()->GetNbins(), i, 1, old->GetBinContent(old->GetXaxis()->GetNbins(), i, 1) + old->GetBinContent(old->GetXaxis()->GetNbins()+1, i, 1));
    old->SetBinError(1, i, 1, TMath::Sqrt(TMath::Power(old->GetBinContent(1, i, 1), 2) + TMath::Power(old->GetBinContent(0, i, 1), 2)));
    old->SetBinError(old->GetXaxis()->GetNbins(), i, 1, TMath::Sqrt(TMath::Power(old->GetBinContent(old->GetXaxis()->GetNbins(), i, 1), 2) + TMath::Power(old->GetBinContent(old->GetXaxis()->GetNbins()+1, i, 1), 2)));
    old->SetBinContent(0, i, 1, 0.);
    old->SetBinContent(old->GetXaxis()->GetNbins()+1, i, 1, 0.);
    old->SetBinError(0, i, 1, 0.);
    old->SetBinError(old->GetXaxis()->GetNbins()+1, i, 1, 0.);

    old->SetBinContent(1, i, 1, old->GetBinContent(1, i, 1) + old->GetBinContent(1, i, 0));
    old->SetBinContent(1, i, old->GetZaxis()->GetNbins(), old->GetBinContent(1, i, old->GetZaxis()->GetNbins()) + old->GetBinContent(1, i, old->GetZaxis()->GetNbins()+1));
    old->SetBinError(1, i, 1, TMath::Sqrt(TMath::Power(old->GetBinContent(1, i, 1), 2) + TMath::Power(old->GetBinContent(1, i, 0), 2)));
    old->SetBinError(1, i, old->GetZaxis()->GetNbins(), TMath::Sqrt(TMath::Power(old->GetBinContent(1, i, old->GetZaxis()->GetNbins()), 2) + TMath::Power(old->GetBinContent(1, i, old->GetZaxis()->GetNbins()+1), 2)));
    old->SetBinContent(1, i, 0, 0.);
    old->SetBinContent(1, i, old->GetZaxis()->GetNbins()+1, 0.);
    old->SetBinError(1, i, 0, 0.);
    old->SetBinError(1, i, old->GetZaxis()->GetNbins()+1, 0.);
  }
  for (int i = 1; i < old->GetZaxis()->GetNbins()+1; i++) {
    old->SetBinContent(1, 1, i, old->GetBinContent(1, 1, i) + old->GetBinContent(0, 1, i));
    old->SetBinContent(old->GetXaxis()->GetNbins(), 1, i, old->GetBinContent(old->GetXaxis()->GetNbins(), 1, i) + old->GetBinContent(old->GetXaxis()->GetNbins()+1, 1, i));
    old->SetBinError(1, 1, i, TMath::Sqrt(TMath::Power(old->GetBinContent(1, 1, i), 2) + TMath::Power(old->GetBinContent(0, 1, i), 2)));
    old->SetBinError(old->GetXaxis()->GetNbins(), 1, i, TMath::Sqrt(TMath::Power(old->GetBinContent(old->GetXaxis()->GetNbins(), 1, i), 2) + TMath::Power(old->GetBinContent(old->GetXaxis()->GetNbins()+1, 1, i), 2)));
    old->SetBinContent(0, 1, i, 0.);
    old->SetBinContent(old->GetXaxis()->GetNbins()+1, 1, i, 0.);
    old->SetBinError(0, 1, i, 0.);
    old->SetBinError(old->GetXaxis()->GetNbins()+1, 1, i, 0.);

    old->SetBinContent(1, 1, i, old->GetBinContent(1, 1, i) + old->GetBinContent(1, 0, i));
    old->SetBinContent(1, old->GetYaxis()->GetNbins(), i, old->GetBinContent(1, old->GetYaxis()->GetNbins(), i) + old->GetBinContent(1, old->GetYaxis()->GetNbins()+1, i));
    old->SetBinError(1, 1, i, TMath::Sqrt(TMath::Power(old->GetBinContent(1, 1, i), 2) + TMath::Power(old->GetBinContent(1, 0, i), 2)));
    old->SetBinError(1, old->GetYaxis()->GetNbins(), i, TMath::Sqrt(TMath::Power(old->GetBinContent(1, old->GetYaxis()->GetNbins(), i), 2) + TMath::Power(old->GetBinContent(1, old->GetYaxis()->GetNbins()+1, i), 2)));
    old->SetBinContent(1, 0, i, 0.);
    old->SetBinContent(1, old->GetYaxis()->GetNbins()+1, i, 0.);
    old->SetBinError(1, 0, i, 0.);
    old->SetBinError(1, old->GetYaxis()->GetNbins()+1, i, 0.);
  }

  string name = old->GetName();

  int nxb = 0;
  float nxbins[100];
  int nyb = 0;
  float nybins[100];
  int nzb = 0;
  float nzbins[100];

  if (name.find("WGG_ele_pho0_pt_bb") != string::npos || name.find("WGG_muo_pho0_pt_bb") != string::npos ||
      name.find("WGG_ele_pho0_pt_be") != string::npos || name.find("WGG_muo_pho0_pt_be") != string::npos ||
      name.find("WGG_ele_pho0_pt_eb") != string::npos || name.find("WGG_muo_pho0_pt_eb") != string::npos ||
      name.find("WGG_ele_pho0_pt_ee") != string::npos || name.find("WGG_muo_pho0_pt_ee") != string::npos ||
      name.find("ZGG_ele_pho0_pt_bb") != string::npos || name.find("ZGG_muo_pho0_pt_bb") != string::npos ||
      name.find("ZGG_ele_pho0_pt_be") != string::npos || name.find("ZGG_muo_pho0_pt_be") != string::npos ||
      name.find("ZGG_ele_pho0_pt_eb") != string::npos || name.find("ZGG_muo_pho0_pt_eb") != string::npos ||
      name.find("ZGG_ele_pho0_pt_ee") != string::npos || name.find("ZGG_muo_pho0_pt_ee") != string::npos ||
      name.find("WGG_ele_pho1_pt_bb") != string::npos || name.find("WGG_muo_pho1_pt_bb") != string::npos ||
      name.find("WGG_ele_pho1_pt_be") != string::npos || name.find("WGG_muo_pho1_pt_be") != string::npos ||
      name.find("WGG_ele_pho1_pt_eb") != string::npos || name.find("WGG_muo_pho1_pt_eb") != string::npos ||
      name.find("WGG_ele_pho1_pt_ee") != string::npos || name.find("WGG_muo_pho1_pt_ee") != string::npos ||
      name.find("ZGG_ele_pho1_pt_bb") != string::npos || name.find("ZGG_muo_pho1_pt_bb") != string::npos ||
      name.find("ZGG_ele_pho1_pt_be") != string::npos || name.find("ZGG_muo_pho1_pt_be") != string::npos ||
      name.find("ZGG_ele_pho1_pt_eb") != string::npos || name.find("ZGG_muo_pho1_pt_eb") != string::npos ||
      name.find("ZGG_ele_pho1_pt_ee") != string::npos || name.find("ZGG_muo_pho1_pt_ee") != string::npos ||
      name.find("WGG_ele_pho0_pho1_pt_bb") != string::npos || name.find("WGG_muo_pho0_pho1_pt_bb") != string::npos ||
      name.find("WGG_ele_pho0_pho1_pt_be") != string::npos || name.find("WGG_muo_pho0_pho1_pt_be") != string::npos ||
      name.find("WGG_ele_pho0_pho1_pt_eb") != string::npos || name.find("WGG_muo_pho0_pho1_pt_eb") != string::npos ||
      name.find("WGG_ele_pho0_pho1_pt_ee") != string::npos || name.find("WGG_muo_pho0_pho1_pt_ee") != string::npos ||
      name.find("ZGG_ele_pho0_pho1_pt_bb") != string::npos || name.find("ZGG_muo_pho0_pho1_pt_bb") != string::npos ||
      name.find("ZGG_ele_pho0_pho1_pt_be") != string::npos || name.find("ZGG_muo_pho0_pho1_pt_be") != string::npos ||
      name.find("ZGG_ele_pho0_pho1_pt_eb") != string::npos || name.find("ZGG_muo_pho0_pho1_pt_eb") != string::npos ||
      name.find("ZGG_ele_pho0_pho1_pt_ee") != string::npos || name.find("ZGG_muo_pho0_pho1_pt_ee") != string::npos) {

    nxb = 6;
    nxbins[0] = 1;
    nxbins[1] = 1;
    nxbins[2] = 1;
    nxbins[3] = 1;
    nxbins[4] = 5;
    nxbins[5] = 6;

    nyb = 6;
    nybins[0] = 1;
    nybins[1] = 1;
    nybins[2] = 1;
    nybins[3] = 1;
    nybins[4] = 5;
    nybins[5] = 6;

    nzb = 6;
    nzbins[0] = 1;
    nzbins[1] = 1;
    nzbins[2] = 1;
    nzbins[3] = 1;
    nzbins[4] = 5;
    nzbins[5] = 6;

  } else if (name.find("WGG_ele_pho0_pho1_pt") != string::npos || name.find("WGG_muo_pho0_pho1_pt") != string::npos ||
             name.find("ZGG_ele_pho0_pho1_pt") != string::npos || name.find("ZGG_muo_pho0_pho1_pt") != string::npos ||
             name.find("WGG_ele_pho0_pt") != string::npos || name.find("WGG_muo_pho0_pt") != string::npos ||
             name.find("WGG_ele_pho1_pt") != string::npos || name.find("WGG_muo_pho1_pt") != string::npos ||
             name.find("ZGG_ele_pho0_pt") != string::npos || name.find("ZGG_muo_pho0_pt") != string::npos ||
             name.find("ZGG_ele_pho1_pt") != string::npos || name.find("ZGG_muo_pho1_pt") != string::npos) {

    nxb = 6;
    nxbins[0] = 4;
    nxbins[1] = 4;
    nxbins[2] = 4;
    nxbins[3] = 4;
    nxbins[4] = 20;
    nxbins[5] = 24;

    nyb = 6;
    nybins[0] = 4;
    nybins[1] = 4;
    nybins[2] = 4;
    nybins[3] = 4;
    nybins[4] = 20;
    nybins[5] = 24;

    nzb = 6;
    nzbins[0] = 1;
    nzbins[1] = 1;
    nzbins[2] = 1;
    nzbins[3] = 1;
    nzbins[4] = 5;
    nzbins[5] = 6;

  } else if (name.find("WG_ele_t") != string::npos || name.find("WG_muo_t") != string::npos ||
             name.find("ZG_ele_t") != string::npos || name.find("ZG_muo_t") != string::npos) {

    nxb = 6;
    nxbins[0] = 1;
    nxbins[1] = 1;
    nxbins[2] = 1;
    nxbins[3] = 1;
    nxbins[4] = 5;
    nxbins[5] = 6;

    nyb = 2;
    nybins[0] = 1;
    nybins[1] = 1;

    nzb = 2;
    nzbins[0] = 1;
    nzbins[1] = 1;

  } else {
    return old;
  }

  if (old->GetNbinsX() == nxb && old->GetNbinsY() == nyb && old->GetNbinsZ() == nzb) {
    cout << "ERROR: " << name << " already rebinned" << endl;
    return old;
  }

  int nx = 0;
  for (int i = 0; i < nxb; i++) {
    nx = nx + nxbins[i];
  }
  if (nx != old->GetNbinsX()) {
    cout << "ERROR: " << name << " - " << nx << " " << old->GetNbinsX() << endl;
    return old;
  }
  int ny = 0;
  for (int i = 0; i < nyb; i++) {
    ny = ny + nybins[i];
  }
  if (ny != old->GetNbinsY()) {
    cout << "ERROR: " << name << " - "  << ny << " " << old->GetNbinsY() << endl;
    return old;
  }
  int nz = 0;
  for (int i = 0; i < nzb; i++) {
    nz = nz + nzbins[i];
  }
  if (nz != old->GetNbinsZ()) {
    cout << "ERROR: " << name << " - " << nz << " " << old->GetNbinsZ() << endl;
    return old;
  }

  int sx = 0;
  const int nxbs = nxb + 1;
  double xbins[nxbs];
  for (int i = 0; i < nxb + 1; i++) {
    if (i > 0) sx = sx + nxbins[i-1];
    xbins[i] = old->GetXaxis()->GetBinUpEdge(sx);
  }
  int sy = 0;
  const int nybs = nyb + 1;
  double ybins[nybs];
  for (int i = 0; i < nyb + 1; i++) {
    if (i > 0) sy = sy + nybins[i-1];
    ybins[i] = old->GetYaxis()->GetBinUpEdge(sy);
  }
  int sz = 0;
  const int nzbs = nzb + 1;
  double zbins[nzbs];
  for (int i = 0; i < nzb + 1; i++) {
    if (i > 0) sz = sz + nzbins[i-1];
    zbins[i] = old->GetZaxis()->GetBinUpEdge(sz);
  }

  TH3D* tmp = new TH3D("tmp", old->GetTitle(), nxb, xbins, nyb, ybins, nzb, zbins);

  for (int i = 0; i < old->GetXaxis()->GetNbins()+2; i++) {
    for (int j = 0; j < old->GetYaxis()->GetNbins()+2; j++) {
      for (int k = 0; k < old->GetZaxis()->GetNbins()+2; k++) {
        float c2 = old->GetBinContent(i, j, k);
        float e2 = old->GetBinError(i, j, k);
        float c2b = tmp->GetBinContent(tmp->GetXaxis()->FindBin(old->GetXaxis()->GetBinCenter(i)), tmp->GetYaxis()->FindBin(old->GetYaxis()->GetBinCenter(j)), tmp->GetZaxis()->FindBin(old->GetZaxis()->GetBinCenter(k)));
        float e2b = tmp->GetBinError(tmp->GetXaxis()->FindBin(old->GetXaxis()->GetBinCenter(i)), tmp->GetYaxis()->FindBin(old->GetYaxis()->GetBinCenter(j)), tmp->GetZaxis()->FindBin(old->GetZaxis()->GetBinCenter(k)));
        tmp->SetBinContent(tmp->GetXaxis()->FindBin(old->GetXaxis()->GetBinCenter(i)), tmp->GetYaxis()->FindBin(old->GetYaxis()->GetBinCenter(j)), tmp->GetZaxis()->FindBin(old->GetZaxis()->GetBinCenter(k)), c2+c2b);
        tmp->SetBinError(tmp->GetXaxis()->FindBin(old->GetXaxis()->GetBinCenter(i)), tmp->GetYaxis()->FindBin(old->GetYaxis()->GetBinCenter(j)), tmp->GetZaxis()->FindBin(old->GetZaxis()->GetBinCenter(k)), TMath::Sqrt(e2*e2+e2b*e2b));
      }
    }
  }

  old->Delete();
  tmp->SetName(name.c_str());

  return tmp;
}*/

