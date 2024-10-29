#include "AnalysisFramework.hh"
#include "FastY3.hh"
#include "LundEEGenerator.hh"
#include <iomanip>

using namespace std;
using namespace panscales;
using namespace fjcore;
using namespace fjcore::contrib;

//----------------------------------------------------------------------
/// Code to examine a couple of Soft Drop observables in e+e- events.
///
/// example of command line is 
/**
  build-double/example-ee-Lund -shower panglobal -beta 0 -process ee2qq \
      -physical-coupling -rts 91.1876 -match-process -nev 100000 \
      -out example-results/example-lund-ee.dat
*/
/// This should run in a few seconds.
/// 
/// To activate NNLL, one can run
/** 
  build-double/example-ee -shower panglobal -split-dipole-frame -beta 0 -process ee2qq \
      -physical-coupling -nloops 3 -rts 91.1876 -nev 100000 \
      -no-spin-corr -colour CFHalfCA -nnll-sudakov -match-process \
      -out example-results/example-ee-nnll.dat
*/
class ExampleEELund : public  AnalysisFramework {
public:
  JetDefinition _jet_def;
  RecursiveLundEEGenerator _lund_ee_gen; 
  /// ctor
  ExampleEELund(CmdLine * cmdline_in)
    : AnalysisFramework(cmdline_in) {
      // check that we are running with the correct process
      if (!dynamic_cast<ProcessZ2qq*>(f_process.get())){
        throw runtime_error("This analysis is only to be used with a e+e-->qqbar process");
      }
    _jet_def = JetDefinition(ee_genkt_algorithm, 1.0, 0.0);
    }

  //----------------------------------------------------------------------
  /// options needed by the "user" analysis can be set here
  void user_startup() override {
    
    // declare histograms and binning
    Binning zg_binning(0.0, 0.5, 0.01);
    cumul_hists_err["z" ].declare(zg_binning);   
    Binning thg_binning(0.0, 1.0, 0.1); 
    cumul_hists_err["theta"].declare(thg_binning);    

  }

  //----------------------------------------------------------------------
  /// this gets called once for every event and should carry
  /// out the analysis and output histograms
  void user_analyse_event() override {
   // Main loop over the PanScales event
   std::vector<PseudoJet> particles; 
   //std::cout << "Number of particles " << f_event.size() << std::endl;
    for (auto &p:f_event.particles()){
      PseudoJet q(p.px(),p.py(), p.pz(), p.E());
      if(abs(p.pdgid()) != 11) particles.push_back(q); // only final-state particles
    }
    // We cluster the event in C/A ycut=1 jets. 
   // std::cout << "Number of particles " << particles.size() << std::endl;
    if(particles.size()<2) return;
    ClusterSequence cs(particles, _jet_def);

    // Obtain the declusterings 
    std::vector<LundEEDeclustering> declusterings;
    declusterings = _lund_ee_gen.result(cs);
    // it should be only one in this case
    //std::cout << "Number of declusterings " << declusterings.size() << std::endl;
    for (const auto & d : declusterings){
        cumul_hists_err["z"].add_entry(d.z(), event_weight());
        cumul_hists_err["theta"].add_entry(2*atan(exp(-d.eta())), event_weight());
    }
    
  }
};
  
//----------------------------------------------------------------------
/// a minimal main program
int main(int argc, char** argv) {
  CmdLine cmdline(argc,argv,true);
  ExampleEELund driver(&cmdline);
  driver.run();
}
