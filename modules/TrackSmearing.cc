/** \class TrackSmearing
 *
 *  Performs d0, dZ, p, Theta, Phi smearing of tracks.
 *
 *  \authors A. Hart, M. Selvaggi
 *
*/

#include "modules/TrackSmearing.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "TMath.h"
#include "TString.h"
#include "TFormula.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TProfile2D.h"

#include <algorithm> 
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <utility>

using namespace std;

//------------------------------------------------------------------------------

TrackSmearing::TrackSmearing() :
  fD0Formula(0), fDZFormula(0), fPFormula(0), fCtgThetaFormula(0), fPhiFormula(0), fItInputArray(0)
{
  fD0Formula = new DelphesFormula;
  fDZFormula = new DelphesFormula;
  fPFormula = new DelphesFormula;
  fCtgThetaFormula = new DelphesFormula;
  fPhiFormula = new DelphesFormula;
}

//------------------------------------------------------------------------------

TrackSmearing::~TrackSmearing()
{
  if(fD0Formula) delete fD0Formula;
  if(fDZFormula) delete fDZFormula;
  if(fPFormula) delete fPFormula;
  if(fCtgThetaFormula) delete fCtgThetaFormula;
  if(fPhiFormula) delete fPhiFormula;
}

//------------------------------------------------------------------------------

void TrackSmearing::Init()
{
  fBz = GetDouble("Bz", 0.0);

  // read resolution formula

  // !!! IF WE WANT TO KEEP ROOT INPUT !!!
  if (string (GetString("D0ResolutionFormula", "0.0")) != "0.0")
    {
      fD0Formula->Compile(GetString("D0ResolutionFormula", "0.0"));
      fUseD0Formula = true;
    }
  else
    {
      fD0ResolutionFile = GetString("D0ResolutionFile", "errors.root");
      fD0ResolutionHist = GetString("D0ResolutionHist", "d0");
      fUseD0Formula = false;
    }
  if (string (GetString("DZResolutionFormula", "0.0")) != "0.0")
    {
      fDZFormula->Compile(GetString("DZResolutionFormula", "0.0"));
      fUseDZFormula = true;
    }
  else
    {
      fDZResolutionFile = GetString("DZResolutionFile", "errors.root");
      fDZResolutionHist = GetString("DZResolutionHist", "dz");
      fUseDZFormula = false;
    }
  if (string (GetString("PResolutionFormula", "0.0")) != "0.0")
    {
      fPFormula->Compile(GetString("PResolutionFormula", "0.0"));
      fUsePFormula = true;
    }
  else
    {
      fPResolutionFile = GetString("PResolutionFile", "errors.root");
      fPResolutionHist = GetString("PResolutionHist", "p");
      fUsePFormula = false;
    }
  if (string (GetString("CtgThetaResolutionFormula", "0.0")) != "0.0")
    {
      fCtgThetaFormula->Compile(GetString("CtgThetaResolutionFormula", "0.0"));
      fUseCtgThetaFormula = true;
    }
  else
    {
      fCtgThetaResolutionFile = GetString("CtgThetaResolutionFile", "errors.root");
      fCtgThetaResolutionHist = GetString("CtgThetaResolutionHist", "ctgTheta");
      fUseCtgThetaFormula = false;
    }
  if (string (GetString("PhiResolutionFormula", "0.0")) != "0.0")
    {
      fPhiFormula->Compile(GetString("PhiResolutionFormula", "0.0"));
      fUsePhiFormula = true;
    }
  else
    {
      fPhiResolutionFile = GetString("PhiResolutionFile", "errors.root");
      fPhiResolutionHist = GetString("PhiResolutionHist", "phi");
      fUsePhiFormula = false;
    }

  fApplyToPileUp = GetBool("ApplyToPileUp", true);

  // import input array

  fInputArray = ImportArray(GetString("InputArray", "ParticlePropagator/stableParticles"));
  fItInputArray = fInputArray->MakeIterator();
 
  // import beamspot
  try
  {
    fBeamSpotInputArray = ImportArray(GetString("BeamSpotInputArray", "BeamSpotFilter/beamSpotParticle"));
  }
  catch(runtime_error &e)
  {
    fBeamSpotInputArray = 0;
  }  
 
  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "stableParticles"));
}

//------------------------------------------------------------------------------

void TrackSmearing::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void TrackSmearing::Process()
{
  Int_t iCandidate = 0;
  TLorentzVector beamSpotPosition;
  Candidate *candidate, *mother;
  Double_t pt, eta, d0, d0Error, trueD0, dz, dzError, trueDZ, p, pError, trueP, ctgTheta, ctgThetaError, trueCtgTheta, phi, phiError, truePhi;
  Double_t op, ophi, opx, opy, opt;
  Double_t x, y, z, t, px, py, pz, theta;
  Double_t q, r, phi_0;
  Double_t xd, yd, zd, delta;
  Double_t bsx, bsy, bsz;
  Int_t sign;
  const Double_t c_light = 2.99792458E8;
  TProfile2D *d0ErrorHist = NULL,
             *dzErrorHist = NULL,
             *pErrorHist = NULL,
             *ctgThetaErrorHist = NULL,
             *phiErrorHist = NULL;

  if (!fBeamSpotInputArray || fBeamSpotInputArray->GetSize () == 0) 
    beamSpotPosition.SetXYZT(0.0, 0.0, 0.0, 0.0);
  else
  {
    Candidate &beamSpotCandidate = *((Candidate *) fBeamSpotInputArray->At (0));
    beamSpotPosition = beamSpotCandidate.Position;
  }

  if (!fUseD0Formula)
  {
     TFile *fin = TFile::Open (fD0ResolutionFile.c_str ());
     d0ErrorHist = (TProfile2D *) fin->Get (fD0ResolutionHist.c_str ());
     d0ErrorHist->SetDirectory (0);
     fin->Close ();
  }
  if (!fUseDZFormula)
  {
     TFile *fin = TFile::Open (fDZResolutionFile.c_str ());
     dzErrorHist = (TProfile2D *) fin->Get (fDZResolutionHist.c_str ());
     dzErrorHist->SetDirectory (0);
     fin->Close ();
  }
  if (!fUsePFormula)
  {
     TFile *fin = TFile::Open (fPResolutionFile.c_str ());
     pErrorHist = (TProfile2D *) fin->Get (fPResolutionHist.c_str ());
     pErrorHist->SetDirectory (0);
     fin->Close ();
  }
  if (!fUseCtgThetaFormula)
  {
     TFile *fin = TFile::Open (fCtgThetaResolutionFile.c_str ());
     ctgThetaErrorHist = (TProfile2D *) fin->Get (fCtgThetaResolutionHist.c_str ());
     ctgThetaErrorHist->SetDirectory (0);
     fin->Close ();
  }
  if (!fUsePhiFormula)
  {
     TFile *fin = TFile::Open (fPhiResolutionFile.c_str ());
     phiErrorHist = (TProfile2D *) fin->Get (fPhiResolutionHist.c_str ());
     phiErrorHist->SetDirectory (0);
     fin->Close ();
  }

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
   
    const TLorentzVector &momentum = candidate->Momentum;
    const TLorentzVector &position = candidate->InitialPosition;
    const TLorentzVector &omomentum = candidate->OriginalMomentum;
   
    pt = momentum.Pt();
    eta = momentum.Eta();

    d0 = trueD0 = candidate->D0;
    dz = trueDZ = candidate->DZ;
     
    p = trueP = candidate->P;
    ctgTheta = trueCtgTheta = candidate->CtgTheta;
    phi = truePhi = candidate->Phi;
    
    op = omomentum.P();
    opx = omomentum.Px();
    opy = omomentum.Py();
    opt = omomentum.Pt();
    ophi = TMath::ATan2(opy, opx);

    q = candidate->Charge;
    r = opt / (q * fBz) * 1.0E9/c_light;        // in [m]
    phi_0 = TMath::ATan2(opy, opx); // [rad] in [-pi, pi]
    // check which solution for (xd,yd)->(x,y) to use, before smearing
    xd = candidate->Xd*1.0E-3;
    yd = candidate->Yd*1.0E-3;
    x = position.X ()*1.0E-3;
    y = position.Y ()*1.0E-3;
    sign = 1;
    auto sol = solvePos(xd,yd,r,phi_0,sign);
    if(TMath::Abs(x-sol.first)>0.0001) sign = -1;

    if (fUseD0Formula)
      d0Error = fD0Formula->Eval(pt, eta);
    else
    {
       Int_t xbin, ybin;

       xbin = pt < d0ErrorHist->GetXaxis ()->GetXmax () ? d0ErrorHist->GetXaxis ()->FindBin (pt)
            : d0ErrorHist->GetXaxis ()->GetBinCenter (d0ErrorHist->GetXaxis ()->GetNbins ());
       ybin = d0ErrorHist->GetYaxis ()->FindBin (TMath::Abs (eta));
       d0Error = d0ErrorHist->GetBinContent (xbin, ybin);
       if (!d0Error)
          d0Error = -1.0;
    }
    if (d0Error < 0.0)
      continue;

    if (fUseDZFormula)
      dzError = fDZFormula->Eval(pt, eta);
    else
    {
       Int_t xbin, ybin;

       xbin = pt < dzErrorHist->GetXaxis ()->GetXmax () ? dzErrorHist->GetXaxis ()->FindBin (pt)
            : dzErrorHist->GetXaxis ()->GetBinCenter (dzErrorHist->GetXaxis ()->GetNbins ());
       ybin = dzErrorHist->GetYaxis ()->FindBin (TMath::Abs (eta));
       dzError = dzErrorHist->GetBinContent (xbin, ybin);
       if (!dzError)
          dzError = -1.0;
    }
    if (dzError < 0.0)
      continue;

    if (fUsePFormula)
      pError = fPFormula->Eval(pt, eta) * p;
    else
    {
       Int_t xbin, ybin;

       xbin = pt < pErrorHist->GetXaxis ()->GetXmax () ? pErrorHist->GetXaxis ()->FindBin (pt)
            : pErrorHist->GetXaxis ()->GetBinCenter (pErrorHist->GetXaxis ()->GetNbins ());
       ybin = pErrorHist->GetYaxis ()->FindBin (TMath::Abs (eta));
       pError = pErrorHist->GetBinContent (xbin, ybin) * p;
       if (!pError)
          pError = -1.0;
    }
    if (pError < 0.0)
      continue;

    if (fUseCtgThetaFormula)
      ctgThetaError = fCtgThetaFormula->Eval(pt, eta);
    else
    {
       Int_t xbin, ybin;

       xbin = pt < ctgThetaErrorHist->GetXaxis ()->GetXmax () ? ctgThetaErrorHist->GetXaxis ()->FindBin (pt)
            : ctgThetaErrorHist->GetXaxis ()->GetBinCenter (ctgThetaErrorHist->GetXaxis ()->GetNbins ());
       ybin = ctgThetaErrorHist->GetYaxis ()->FindBin (TMath::Abs (eta));
       ctgThetaError = ctgThetaErrorHist->GetBinContent (xbin, ybin);
       if (!ctgThetaError)
         ctgThetaError = -1.0;
    }
    if (ctgThetaError < 0.0)
      continue;

    if (fUsePhiFormula)
      phiError = fPhiFormula->Eval(pt, eta);
    else
    {
       Int_t xbin, ybin;

       xbin = pt < phiErrorHist->GetXaxis ()->GetXmax () ? phiErrorHist->GetXaxis ()->FindBin (pt)
            : phiErrorHist->GetXaxis ()->GetBinCenter (phiErrorHist->GetXaxis ()->GetNbins ());
       ybin = phiErrorHist->GetYaxis ()->FindBin (TMath::Abs (eta));
       phiError = phiErrorHist->GetBinContent (xbin, ybin);
       if (!phiError)
          phiError = -1.0;
    }
    if (phiError < 0.0)
      continue;

    if (fApplyToPileUp || !candidate->IsPU)
    {
       d0 = gRandom->Gaus(d0, d0Error);
       dz = gRandom->Gaus(dz, dzError);
       p = gRandom->Gaus(p, pError);
       ctgTheta = gRandom->Gaus(ctgTheta, ctgThetaError);
       phi = gRandom->Gaus(phi, phiError);
    }

    if(p < 0.0) continue;
    while (phi > TMath::Pi ()) phi -= TMath::TwoPi ();
    while (phi <= -TMath::Pi ()) phi += TMath::TwoPi ();

    mother = candidate;
    candidate = static_cast<Candidate*>(candidate->Clone());
    candidate->D0 = d0;
    candidate->DZ = dz;
    candidate->P = p;
    candidate->CtgTheta = ctgTheta;
    candidate->Phi = phi;

    // propagate smearing to both original and perigee momentum
    theta = TMath::ACos(ctgTheta / TMath::Sqrt (1.0 + ctgTheta * ctgTheta));
    ophi = ophi + (phi - truePhi);
    while (ophi > TMath::Pi ()) ophi -= TMath::TwoPi ();
    while (ophi <= -TMath::Pi ()) ophi += TMath::TwoPi ();
    op = op + (p - trueP);
    candidate->OriginalMomentum.SetPx (op * TMath::Cos (ophi) * TMath::Sin (theta));
    candidate->OriginalMomentum.SetPy (op * TMath::Sin (ophi) * TMath::Sin (theta));
    candidate->OriginalMomentum.SetPz (op * TMath::Cos (theta));
    candidate->OriginalMomentum.SetE (candidate->OriginalMomentum.Pt () * TMath::CosH (eta));
    candidate->Momentum.SetPx (p * TMath::Cos (phi) * TMath::Sin (theta));
    candidate->Momentum.SetPy (p * TMath::Sin (phi) * TMath::Sin (theta));
    candidate->Momentum.SetPz (p * TMath::Cos (theta));
    candidate->Momentum.SetE (candidate->Momentum.Pt () * TMath::CosH (eta));
    candidate->PT = candidate->Momentum.Pt ();

    bsx = beamSpotPosition.X()*1.0E-3;
    bsy = beamSpotPosition.Y()*1.0E-3;
    bsz = beamSpotPosition.Z()*1.0E-3;

    xd = candidate->Xd*1.0E-3;
    yd = candidate->Yd*1.0E-3;
    zd = candidate->Zd*1.0E-3;
    px = candidate->Momentum.Px ();
    py = candidate->Momentum.Py ();
    pz = candidate->Momentum.Pz ();
    pt = candidate->Momentum.Pt ();
    opx = candidate->OriginalMomentum.Px ();
    opy = candidate->OriginalMomentum.Py ();
    opt = candidate->OriginalMomentum.Pt ();
    
    // -- solve for delta: d0' = ( (xd+delta-bsx)*py' - (yd+delta-bsy)*px' )/pt'
    delta = (d0*1.0E-3 * pt - (xd - bsx)*py + (yd - bsy)*px)/(py - px);
    candidate->Xd = (xd + delta)*1.0E3;
    candidate->Yd = (yd + delta)*1.0E3;
    xd = candidate->Xd*1.0E-3;
    yd = candidate->Yd*1.0E-3;
    // -- solve for delta: dz' = zd + delta - (...)
    candidate->Zd = (dz*1.0E-3 + ((xd - bsx) * px + (yd - bsy) * py) / pt * (pz / pt))*1.0E3;
    zd = candidate->Zd*1.0E-3;

    // update initial position correspondingly (use original momentum)
    r = opt / (q * fBz) * 1.0E9/c_light;
    phi_0 = TMath::ATan2(opy, opx);
    sol = solvePos(xd,yd,r,phi_0,sign);
    x = sol.first;
    y = sol.second;
    z = zd - (TMath::Hypot(xd, yd) - TMath::Hypot(x, y))*pz/pt;
    t = position.T ();

    candidate->InitialPosition.SetX (x*1.0E3);
    candidate->InitialPosition.SetY (y*1.0E3);
    candidate->InitialPosition.SetZ (z*1.0E3);
    candidate->InitialPosition.SetT (t);

    if (fApplyToPileUp || !candidate->IsPU)
    {
       candidate->ErrorD0 = d0Error;
       candidate->ErrorDZ = dzError;
       candidate->ErrorP = pError;
       candidate->ErrorCtgTheta = ctgThetaError;
       candidate->ErrorPhi = phiError;
       candidate->ErrorPT = ptError (p, ctgTheta, pError, ctgThetaError);
       candidate->TrackResolution = pError/p;
    }
    
    candidate->AddCandidate(mother);
    fOutputArray->Add(candidate);

    iCandidate++;
  }
}

Double_t TrackSmearing::ptError (const Double_t p, const Double_t ctgTheta, const Double_t dP, const Double_t dCtgTheta) const
{
  Double_t a, b;
  a = (p * p * ctgTheta * ctgTheta * dCtgTheta * dCtgTheta) / ((ctgTheta * ctgTheta + 1) * (ctgTheta * ctgTheta + 1) * (ctgTheta * ctgTheta + 1));
  b = (dP * dP) / (ctgTheta * ctgTheta + 1);
  return sqrt (a + b);
}

std::pair<Double_t,Double_t> TrackSmearing::solvePos(Double_t xd, Double_t yd, Double_t r, Double_t phi_0, Int_t sign) const
{
  sign = TMath::Sign(1.0,sign);
  Double_t rd = TMath::Hypot(xd,yd);
  Double_t x = (xd - r*TMath::Sin(phi_0))+sign*xd*TMath::Abs(r)/rd;
  Double_t y = (yd + r*TMath::Cos(phi_0))+sign*yd*TMath::Abs(r)/rd;
  return std::make_pair(x,y);
}

//------------------------------------------------------------------------------
