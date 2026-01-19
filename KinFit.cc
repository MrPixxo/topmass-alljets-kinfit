#include <TLorentzVector.h>
#include <TMatrixD.h>
#include <TObject.h>
#include <iostream>
#include <vector>

#include "KinFitter/TFitConstraintM.h"
#include "KinFitter/TFitParticlePtEtaPhi.h"
#include "KinFitter/TKinFitter.h"
#include "config.h"
using std::cout;
using std::vector;
// B1,B2,W1Prod1,W1Prod2, W2Prod1,
//  W2Prod2, W1, W2, Top1, Top2, TTBar)
struct Selection {
public:
  double chi2;
  double pgof;
  vector<TLorentzVector *> bestPermutation;
  vector<double> fitJetpt;
  vector<double> fitJeteta;
  vector<double> fitJetphi;
  vector<double> fitJetM;
  int combinationType;
  Selection(size_t size)
      : fitJetpt(size), fitJeteta(size), fitJetphi(size), fitJetM(size) {}
};

TMatrixD setupMatrix(const TLorentzVector *object, bool isB) {
  TMatrixD CovM(3, 3);
  CovM.Zero();
  const double et = object->E() * std::fabs(sin(object->Theta()));
  const double eta = object->Eta();
  CovM(0, 0) = pow(CalcEt(et, eta, isB), 2);
  CovM(0, 0) *= pow(getEtaDependentScaleFactor(*object), 2);
  CovM(1, 1) = pow(CalcEta(et, eta, isB), 2);
  CovM(2, 2) = pow(CalcPhi(et, eta, isB), 2);
  return CovM;
}
// now overloaded to account for adding Jets
TMatrixD setupMatrix(const TLorentzVector *object1,
                     const TLorentzVector *object2, bool isB) {
  TMatrixD CovM(3, 3);
  CovM.Zero();
  const double et1 = object1->E() * std::fabs(sin(object1->Theta()));
  const double eta1 = object1->Eta();
  const double et2 = object2->E() * std::fabs(sin(object2->Theta()));
  const double eta2 = object2->Eta();

  CovM(0, 0) =
      pow(CalcEt(et1, eta1, isB) * getEtaDependentScaleFactor(*object1), 2) +
      pow(CalcEt(et2, eta2, isB) * getEtaDependentScaleFactor(*object2), 2);
  CovM(1, 1) =
      pow(CalcEta(et1, eta1, isB), 2) + pow(CalcEta(et2, eta2, isB), 2);
  CovM(2, 2) =
      pow(CalcPhi(et1, eta1, isB), 2) + pow(CalcPhi(et2, eta2, isB), 2);
  return CovM;
}

void applyKinFit(vector<TLorentzVector *> jetSelection,
                 vector<TMatrixD> decayprodmats,
                 struct Selection &currentSelection) {
  TKinFitter *fitter_ = new TKinFitter("TopKinFitter", "TopKinFitter");
  fitter_->setMaxNbIter(200);
  fitter_->setMaxDeltaS(5e-5);
  fitter_->setMaxF(1e-4);
  fitter_->setVerbosity(0);

  // definiton of covariance matrices
  TMatrixD mB1 = decayprodmats[0];
  TMatrixD mB2 = decayprodmats[1];
  TMatrixD mW1P1 = decayprodmats[2];
  TMatrixD mW1P2 = decayprodmats[3];
  TMatrixD mW2P1 = decayprodmats[4];
  TMatrixD mW2P2 = decayprodmats[5];

  auto B1 = new TFitParticlePtEtaPhi("B1", "B1", jetSelection[0], &mB1);
  auto B2 = new TFitParticlePtEtaPhi("B2", "B2", jetSelection[1], &mB2);
  auto W1Prod1 =
      new TFitParticlePtEtaPhi("W1Prod1", "W1Prod1 ", jetSelection[2], &mW1P1);
  auto W1Prod2 =
      new TFitParticlePtEtaPhi("W1Prod2", "W1Prod2", jetSelection[3], &mW1P2);
  auto W2Prod1 =
      new TFitParticlePtEtaPhi("W2Prod1", "W2Prod1 ", jetSelection[4], &mW2P1);
  auto W2Prod2 =
      new TFitParticlePtEtaPhi("W2Prod2", "W2Prod2", jetSelection[5], &mW2P2);

  for (auto p : {B1, B2, W1Prod1, W1Prod2, W2Prod1, W2Prod2})
    fitter_->addMeasParticle(p);

  // set constants:
  double mW_ = 80.4;

  // set up constraints and add to fitter_
  auto kW1Mass = new TFitConstraintM("W1Mass", "W1Mass", 0, 0, mW_);
  auto kW2Mass = new TFitConstraintM("W2Mass", "W2Mass", 0, 0, mW_);
  auto kEqualTopMasses =
      new TFitConstraintM("EqualTopMasses", "EqualTopMasses", 0, 0, 0);

  // add particles/jetSelection to constraints
  kW1Mass->addParticles1(W1Prod1, W1Prod2);
  kW2Mass->addParticles1(W2Prod1, W2Prod2);
  kEqualTopMasses->addParticles1(B1, W1Prod1, W1Prod2);
  kEqualTopMasses->addParticles2(B2, W2Prod1, W2Prod2);

  for (auto c : {kW1Mass, kW2Mass, kEqualTopMasses})
    fitter_->addConstraint(c);
  // perform fit
  fitter_->fit();

  double thisChi2 = fitter_->getS();
  if (thisChi2 < currentSelection.chi2 && fitter_->getStatus() == 0) {
    // std::cout << " Fit converged !" << std::endl;
    currentSelection.chi2 = thisChi2;
    currentSelection.bestPermutation = jetSelection;
    currentSelection.pgof = TMath::Prob(thisChi2, fitter_->getNDF());
    // add fitted jets to vectors in order (B1,B2,W1Prod1,W1Prod2, W2Prod1,
    // W2Prod2, W1, W2, Top1, Top2, TTBar)
    for (int i = 0; i < 6; i++) {
      currentSelection.fitJetpt[i] = (fitter_->get4Vec(i)->Pt());
      currentSelection.fitJeteta[i] = (fitter_->get4Vec(i)->Eta());
      currentSelection.fitJetphi[i] = (fitter_->get4Vec(i)->Phi());
      currentSelection.fitJetM[i] = (fitter_->get4Vec(i)->M());
    }
  } else if (currentSelection.chi2 == 10000) {
    currentSelection.bestPermutation = jetSelection;
  }

  delete fitter_;
  for (auto p : {B1, B2, W1Prod1, W1Prod2, W2Prod1, W2Prod2})
    delete p;
  for (auto c : {kW1Mass, kW2Mass, kEqualTopMasses})
    delete c;
}

void tryCombinations(TLorentzVector *B1, TLorentzVector *B2, TLorentzVector *J1,
                     TLorentzVector *J2, TLorentzVector *J3, TLorentzVector *J4,
                     TMatrixD mB1, TMatrixD mB2, TMatrixD mJ1, TMatrixD mJ2,
                     TMatrixD mJ3, TMatrixD mJ4,
                     struct Selection &bestSelection) {

  applyKinFit({B1, B2, J1, J2, J3, J4}, {mB1, mB2, mJ1, mJ2, mJ3, mJ4},
              bestSelection);
  applyKinFit({B1, B2, J1, J3, J2, J4}, {mB1, mB2, mJ1, mJ3, mJ2, mJ4},
              bestSelection);
  applyKinFit({B1, B2, J1, J4, J2, J3}, {mB1, mB2, mJ1, mJ4, mJ2, mJ3},
              bestSelection);
  applyKinFit({B1, B2, J3, J4, J1, J2}, {mB1, mB2, mJ3, mJ4, mJ1, mJ2},
              bestSelection);
  applyKinFit({B1, B2, J2, J4, J1, J3}, {mB1, mB2, mJ2, mJ4, mJ1, mJ3},
              bestSelection);
  applyKinFit({B1, B2, J2, J3, J1, J4}, {mB1, mB2, mJ2, mJ3, mJ1, mJ4},
              bestSelection);
}

/*std::tuple<vector<vector<double>>, vector<vector<double>>,
vector<vector<double>>, vector<vector<double>>, vector<vector<double>>>*/
std::tuple<vector<vector<double>>, vector<vector<double>>,
           vector<vector<double>>, vector<vector<double>>, vector<vector<int>>,
           vector<double>, vector<double>>
setBestCombi(
    vector<vector<double>> inputpt, vector<vector<double>> inputeta,
    vector<vector<double>> inputphi,
    vector<vector<double>> inputM /*, vector<vector<double>> gendata*/) {

  /*    vector<double> genpt = gendata[0];
vector<double> geneta = gendata[1];
vector<double> genphi = gendata[2];
vector<double> genM = gendata[3];
*/
  vector<vector<double>> outputpt(inputpt.size());
  vector<vector<double>> outputeta(inputpt.size());
  vector<vector<double>> outputphi(inputpt.size());
  vector<vector<double>> outputm(inputpt.size());
  vector<vector<int>> indexmask(inputpt.size());
  vector<double> outputchi2(inputpt.size());
  vector<double> outputpgof(inputpt.size());
  // vector<vector<TLorentzVector>> outputfitJets(inputpt.size());

  vector<TLorentzVector *> jets;
  for (size_t i = 0; i < inputpt.size(); i++) {
    struct Selection bestSelection(6);
    bestSelection.chi2 = 10000;
    bestSelection.pgof = 0.00
    jets.clear();
    if (inputpt[i].size() > 5) {
      for (size_t j = 0; j < inputpt[i].size(); j++) {
        TLorentzVector *jet = new TLorentzVector();
        jet->SetPtEtaPhiM(inputpt[i][j], inputeta[i][j], inputphi[i][j],
                          inputM[i][j]);
        jets.push_back(jet);
      }
    } else {

      continue;
    }
    TLorentzVector *B1 = jets[0];
    TLorentzVector *B2 = jets[1];
    std::vector<TLorentzVector *> orderedJets(jets.begin() + 2, jets.end());

    // Sort this new vector
    std::sort(
        orderedJets.begin(), orderedJets.end(),
        [](TLorentzVector *a, TLorentzVector *b) { return a->Pt() > b->Pt(); });
    TLorentzVector *J1 = orderedJets[0];
    TLorentzVector *J2 = orderedJets[1];
    TLorentzVector *J3 = orderedJets[2];
    TLorentzVector *J4 = orderedJets[3];

    TMatrixD mB1 = setupMatrix(B1, true);
    TMatrixD mB2 = setupMatrix(B2, true);
    TMatrixD mJ1 = setupMatrix(J1, false);
    TMatrixD mJ2 = setupMatrix(J2, false);
    TMatrixD mJ3 = setupMatrix(J3, false);
    TMatrixD mJ4 = setupMatrix(J4, false);
    vector<TMatrixD> decayprodmats = {mB1, mB2, mJ1, mJ2, mJ3, mJ4};

    tryCombinations(B1, B2, J1, J2, J3, J4, mB1, mB2, mJ1, mJ2, mJ3, mJ4,
                    bestSelection);
    // Create the index vector
    vector<int> indices;
    indices.reserve(bestSelection.bestPermutation.size());

    // For each pointer in bestPermutation, find its index in jets
    for (const auto &ptr : bestSelection.bestPermutation) {
      // Find the iterator position of this pointer in jets
      auto it = std::find(jets.begin(), jets.end(), ptr);

      // Calculate index from iterator position
      if (it != jets.end()) {
        int index = std::distance(jets.begin(), it);
        indices.push_back(index);
      } else {
        // Handle case where pointer is not found in jets
        std::cerr << "Error: Pointer not found in jets vector\n";
      }
    }
    outputpt[i] = bestSelection.fitJetpt;
    outputeta[i] = bestSelection.fitJeteta;
    outputphi[i] = bestSelection.fitJetphi;
    outputm[i] = bestSelection.fitJetM;
    indexmask[i] = indices;
    outputchi2[i] = bestSelection.chi2;
    outputpgof[i] = bestSelection.pgof;
    //  cout << "Fit PT : " << outputpt[i][0] << "Fit eta : " << outputeta[i][0]
    //  << endl;
  }

  /*
vector<TLorentzVector *> genJets;
    for (size_t i = 0; i < gendata[0].size(); i++) {
    TLorentzVector *genjet = new TLorentzVector();
    genjet->SetPtEtaPhiM(genpt[i], geneta[i], genphi[i], genM[i]);
    genJets.push_back(genjet);
}
*/

  for (auto p : jets)
    delete p;
  // for (auto d : genJets)
  //    delete d;
  return std::make_tuple(outputpt, outputeta, outputphi, outputm, indexmask,
                         outputchi2, outputpgof);
}
