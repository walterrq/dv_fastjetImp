#include "SampleAnalyzer/User/Analyzer/MyAnalysis.h"
//Once the fastjet directory has been copied to the Build location, the following paths will be correct.
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include <iostream> // needed for io
#include <cstdio>   // needed for io
using namespace fastjet;
using namespace MA5;
using namespace std;

// -----------------------------------------------------------------------------
// Initialize
// function called one time at the beginning of the analysis
// -----------------------------------------------------------------------------
bool MyAnalysis::Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters)
{

  // Information on the analysis and authors. 
  INFO << "        <><><><><><><><><><><><><><><><><><><><><><><><>" << endmsg;
  INFO << "        <>  Analysis: Jet implementation for MA5      <>" << endmsg;
  INFO << "        <>  Recasted by: W.Rodriguez, D.Zegarra       <>" << endmsg;
  INFO << "        <>  Contact: walter.rodriguez@pucp.edu.pe     <>" << endmsg;
  INFO << "        <>           danilo.zegarra@pucp.edu.pe       <>" << endmsg;
  INFO << "        <>  Based on MadAnalysis 5 v1.8.44            <>" << endmsg;
  INFO << "        <><><><><><><><><><><><><><><><><><><><><><><><>" << endmsg;

  cout << "BEGIN Initialization" << endl;
  // initialize variables, histos
  
  //We define the invisible particles
  PHYSICS->mcConfig().AddInvisibleId(12);  //neutrino e
  PHYSICS->mcConfig().AddInvisibleId(-12); //antineutrino e
  PHYSICS->mcConfig().AddInvisibleId(14);  //neutrino mu
  PHYSICS->mcConfig().AddInvisibleId(-14); //antineutrino mu
  PHYSICS->mcConfig().AddInvisibleId(16);  //neutrino tau
  PHYSICS->mcConfig().AddInvisibleId(-16); //antineutrino tau
  

  //The region analysis is defined
  Manager()->AddRegionSelection("myregion");

  
  Manager()->AddHisto("MET", 100,0.0,1000.0);  

  Manager()->AddHisto("PT leading jet", 100,0.0,400.0);
  Manager()->AddHisto("PT subleading jet", 100,0.0,400.0);

  Manager()->AddHisto("PT leading pion", 100,0.0,400.0);
  Manager()->AddHisto("PT subleading pion", 100,0.0,400.0);

  Manager()->AddHisto("DELTAR for jet and pi", 25,0.0,5.0);

  Manager()->AddHisto("PHI leading pion", 100,-4.0,4.0);

  Manager()->AddHisto("PT leading pion > 20", 100,0.0,400.0);

  Manager()->AddHisto("DELTAR for jet and pi (with PT(pi)>20)", 25,0.0,5.0);



  //Plots for ROOT are defined

  plot_MET = new TH1F("MET", "MET", 100, 0, 1000);

  plot_PTjet = new TH1F("PT leading jet", "PT leading jet", 100, 0, 400); 
  plot_PTsjet = new TH1F("PT subleading jet", "PT subleading jet", 100, 0, 400); 
  

  plot_PTpi = new TH1F("PT leading pion", "PT leading pion", 100, 0, 400); 
  plot_PTspi = new TH1F("PT subleading pion", "PT subleading pion", 100, 0, 400);

  plot_DELTARjpi = new TH1F("DELTAR for jet and pi", "DELTAR for jet and pi", 25, 0.0, 5.0); 

  plot_PHIpi = new TH1F("PHI leading pion", "PHI leading pion", 100, -4.0, 4.0); 

  plot_PTpi20 =new TH1F("PT leading pion > 20", "PT leading pion > 20", 100, 0, 400); 

  plot_DELTARjpi20 = new TH1F("DELTAR for jet and pi (with PT(pi)>20)", "DELTAR for jet and pi (with PT(pi)>20)", 25, 0.0, 5.0); 
  
  
  cout << "END   Initialization" << endl;
  return true;
}

// -----------------------------------------------------------------------------
// Finalize
// function called one time at the end of the analysis
// -----------------------------------------------------------------------------
void MyAnalysis::Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files)
{
  cout << "BEGIN Finalization" << endl;
  //We save hisotgrams in a ROOT file

  TFile* Output1 = new TFile("output1.root", "RECREATE");
  plot_MET -> Write();
  Output1 -> Close();

  TFile* Output2 = new TFile("output2.root", "RECREATE");
  plot_PTjet -> Write();
  Output2 -> Close();

  TFile* Output3 = new TFile("output3.root", "RECREATE");
  plot_PTsjet -> Write();
  Output3 -> Close();

  TFile* Output4 = new TFile("output4.root", "RECREATE");
  plot_PTpi -> Write();
  Output4 -> Close();

  TFile* Output5 = new TFile("output5.root", "RECREATE");
  plot_PTspi -> Write();
  Output5 -> Close();

  TFile* Output6 = new TFile("output6.root", "RECREATE");
  plot_DELTARjpi -> Write();
  Output6 -> Close();

  TFile* Output7 = new TFile("output7.root", "RECREATE");
  plot_PHIpi -> Write();
  Output7 -> Close();

  TFile* Output8 = new TFile("output8.root", "RECREATE");
  plot_PTpi20 -> Write();
  Output8 -> Close();

  TFile* Output9 = new TFile("output9.root", "RECREATE");
  plot_DELTARjpi20 -> Write();
  Output9 -> Close();

  cout << "END   Finalization" << endl;
}

// -----------------------------------------------------------------------------
// Execute
// function called each time one event is read
// -----------------------------------------------------------------------------
bool MyAnalysis::Execute(SampleFormat& sample, const EventFormat& event)
{
    
  if (event.mc()!=0)
 {
  //Initiation of the code for every event
  MAfloat32 __event_weight__ = 1.0;
  if (weighted_events_ && event.mc()!=0) __event_weight__ = event.mc()->weight();
  Manager()->InitializeForNewEvent(__event_weight__);
  
  //Declaration of variables
  double ptjet_leading = 0.0;
  double ptpi_leading = 0.0;

  std::vector<const MCParticleFormat*> basePions;

  //The vector containing the jets is defined
  vector<PseudoJet> jets;

  //A vector of vectors containing the 4-momentum of particles forming a jet is defined
  vector<PseudoJet> parts;

  //Jet is defined under the antikt algorithm with R = 0.4 
  JetDefinition jet_def( antikt_algorithm, .4 ); 

  cout << "---------------NEW EVENT-------------------" << endl;

  //Start of the analysis of every particle
  for (MAuint32 i=0;i<event.mc()->particles().size();i++){
    const MCParticleFormat *iPar = &(event.mc()->particles()[i]);

    //The interest relies only in final state particles, which are the ones passing through the calorimeter
    if (PHYSICS->Id->IsFinalState(iPar)){ 

      //However, this final particles can not be electrons, muons, neutrinos nor photons
      if ((abs(iPar->pdgid()) != 11) && (abs(iPar->pdgid()) != 12) && (abs(iPar->pdgid()) != 13) && (abs(iPar->pdgid()) != 14) && (abs(iPar->pdgid()) != 16) && (abs(iPar->pdgid()) != 18) && (abs(iPar->pdgid()) != 22))
      {
        //A vector that captures the 4-momentum of a particle is defined
        PseudoJet part( iPar->px(), iPar->py(), iPar->pz(), iPar->e() );      

        parts.push_back(part);
      }
    }

    //Jet analysis results are compared to Pions       
    if (abs(iPar->pdgid()) == 211) {
      cout << "id = " << iPar->pdgid() << endl; 
      basePions.push_back(iPar);
    }
  }

  //Once information of every particle from the event that meet the requirements is gathered, the Culstering algorithm proceeds
  ClusterSequence clust(parts, jet_def);

  //jets is given it's value and ordered by pt using fastjet syntax
  jets = sorted_by_pt(clust.inclusive_jets(20));

  //jets are sorted by their PT
  if(jets.size() > 0){ 
    //PT of leading jet is called with a fastjet command, since this object belongs to fastjet
    ptjet_leading = jets[0].perp();
  }

  if(basePions.size() > 0){
  SORTER->sort(basePions, PTordering);
  ptpi_leading = basePions[0]->pt();
  }

  //Start of plots for this event

  //MET for all events is plotted
  Manager()->FillHisto("MET", PHYSICS->Transverse->EventMET(event.mc()) );
  plot_MET -> Fill( PHYSICS->Transverse->EventMET(event.mc()) );

  if (jets.size() > 0){
    Manager()->FillHisto("PT leading jet", jets[0].perp() );
    plot_PTjet -> Fill( jets[0].perp() );

    if(jets.size() > 1){
      Manager()->FillHisto("PT subleading jet", jets[1].perp() );
      plot_PTsjet -> Fill( jets[1].perp() );
    }
  }

  if (basePions.size() > 0){
    Manager()->FillHisto("PT leading pion", basePions[0]->pt() );
    plot_PTpi -> Fill( basePions[0]->pt() );

    if(basePions.size() > 1){
      Manager()->FillHisto("PT subleading pion", basePions[1]->pt() );
      plot_PTspi -> Fill( basePions[1]->pt() );
    }

    if (jets.size() > 0)
    {
      //We define the DELTAR variable between jets and pions because since this data belong to different instances, it can not be conbined in one sintax using both.
      double DELTAR_j_pi = sqrt( pow( (jets[0].eta() - basePions[0]->eta()) ,2) + pow( (jets[0].phi_std() - basePions[0]->phi()) ,2) ) ; //We use phi_std because it has the same domain as Madanalysis' phi (-pi, pi)
      Manager()->FillHisto("DELTAR for jet and pi", DELTAR_j_pi );
      plot_DELTARjpi -> Fill( DELTAR_j_pi );

      if (basePions[0]->pt() > 20.0)
      {
        Manager()->FillHisto("DELTAR for jet and pi (with PT(pi)>20)", DELTAR_j_pi );
        plot_DELTARjpi20 -> Fill( DELTAR_j_pi );
      }
      

    }

    Manager()->FillHisto("PHI leading pion", basePions[0]->phi() );
    plot_PHIpi -> Fill( basePions[0]->phi() );

    //As for jets, a cut in pt of 20 for pions is applied
    if (basePions[0]->pt() > 20.0)
    {
      Manager()->FillHisto("PT leading pion > 20", basePions[0]->pt() );
      plot_PTpi20 -> Fill( basePions[0]->pt() );
    }
  }
   
 }
  return true; 
}

