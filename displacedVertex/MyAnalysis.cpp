#include "SampleAnalyzer/User/Analyzer/MyAnalysis.h"
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
  INFO << "        <>  Analysis: Vertex Production for MA5       <>" << endmsg;
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

  Manager()->AddHisto("Decay vertex time (s)", 200,0,1000);
  Manager()->AddHisto("Decay vertex x (mm)", 200,-1000,1000);
  Manager()->AddHisto("Decay vertex y (mm)", 200,-1000,1000);
  Manager()->AddHisto("Decay vertex z (mm)", 200,-1000,1000);

  //Plots for ROOT are defined
  plot_disvt = new TH1F("Decay vertex time (s)", "Decay vertex time (s)", 100,0,1000);
  plot_disvx = new TH1F("Decay vertex x (mm)", "Decay vertex x (mm)", 200,-1000,1000);
  plot_disvy = new TH1F("Decay vertex y (mm)", "Decay vertex y (mm)", 200,-1000,1000);
  plot_disvz = new TH1F("Decay vertex z (mm)", "Decay vertex z (mm)", 200,-1000,1000);
  
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
  plot_disvt -> Write();
  Output1 -> Close();
  
  TFile* Output2 = new TFile("output2.root", "RECREATE");
  plot_disvx -> Write();
  Output2 -> Close();
  
  TFile* Output3 = new TFile("output3.root", "RECREATE");
  plot_disvy -> Write();
  Output3 -> Close();
  
  TFile* Output4 = new TFile("output4.root", "RECREATE");
  plot_disvz -> Write();
  Output4 -> Close();
  
  
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
  
  std::vector<const MCParticleFormat*> baseMuons;
  
  //Start of the analysis of every particle
  for (MAuint32 i=0;i<event.mc()->particles().size();i++)
    {
      const MCParticleFormat *part = &(event.mc()->particles()[i]);

      //We analyze the case where the muon comes from the decay of the heavy neutrino nu4
      if (part->mothers().size() > 0){
      MCParticleFormat* mother = part->mothers()[0];
        if (mother->pdgid() == 9900012) {
          if (abs(part->pdgid()) == 13) {
            baseMuons.push_back(part);
          }
        }
      }

      if (part->mothers().size() > 1){
      MCParticleFormat* father = part->mothers()[1];
        if (father->pdgid() == 9900012) {
          if (abs(part->pdgid()) == 13) {
            baseMuons.push_back(part);
          }
        }
      }      
    }
     
    
    if(baseMuons.size() == 1){    
      
      SORTER->sort(baseMuons, PTordering);

      //Start of plots for this event
      Manager()->FillHisto("Decay vertex x (mm)", baseMuons[0]->decay_vertex()[0] );
      plot_disvx -> Fill( baseMuons[0]->decay_vertex()[0] );
      
      Manager()->FillHisto("Decay vertex y (mm)", baseMuons[0]->decay_vertex()[1] );
      plot_disvy -> Fill( baseMuons[0]->decay_vertex()[1] );
      
      Manager()->FillHisto("Decay vertex z (mm)", baseMuons[0]->decay_vertex()[2] );
      plot_disvz -> Fill( baseMuons[0]->decay_vertex()[2] );
      
      Manager()->FillHisto("Decay vertex time (s)", baseMuons[0]->decay_vertex()[3] );
      plot_disvt -> Fill( baseMuons[0]->decay_vertex()[3] );
    }

 }
 return true; 
}

