#ifndef analysis_MyAnalysis_h
#define analysis_MyAnalysis_h

#include "SampleAnalyzer/Process/Analyzer/AnalyzerBase.h"
#include "SampleAnalyzer/Interfaces/root/RootMainHeaders.h"

namespace MA5
{
class MyAnalysis : public AnalyzerBase
{
  INIT_ANALYSIS(MyAnalysis,"MyAnalysis")

 public:
  virtual bool Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters);
  virtual void Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files);
  virtual bool Execute(SampleFormat& sample, const EventFormat& event);

 private:
 	TH1F* plot_MET;
 	TH1F* plot_PTjet;
 	TH1F* plot_PTsjet;
 	TH1F* plot_PTpi;
 	TH1F* plot_PTspi;
 	TH1F* plot_DELTARjpi;
 	TH1F* plot_PHIpi;
 	TH1F* plot_PTpi20;
 	TH1F* plot_DELTARjpi20;
       

};
}

#endif
