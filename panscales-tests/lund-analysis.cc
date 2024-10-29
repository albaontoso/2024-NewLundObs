#include "AnalysisFramework.hh"
#include "FastY3.hh"
#include "LundEEGenerator.hh"
#include <iomanip>

using namespace std;
using namespace panscales;
using namespace fjcore;
using namespace fjcore::contrib;

//----------------------------------------------------------------------
/// Code to examine a couple of Lund observables in e+e- events.
///
/// example of command line is 
/**
  ./lund-analysis -shower panglobal -beta 0 -process ee2qq \
      -physical-coupling -rts 1000 -match-process -nev 100000 \
      -out a
*/
/// This should run in a few seconds.

class ExampleEELund : public  AnalysisFramework {
public:
  JetDefinition _jet_def;
  RecursiveLundEEGenerator _lund_ee_gen; 
  double ktmin,ktmax;
  double logthmin, logthmax; // define the collinear window
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
    
    cmdline->start_section("Choice of input parameters"); //---------

    // -- Lund plane windows that define the observable
    ktmin     = cmdline->value("-kt-min", 2.);
    ktmax     = cmdline->value("-kt-max", 200.);
    logthmin  = cmdline->value("-lnth-min", 4);
    logthmax  = cmdline->value("-lnth-max", 5);
    // veto for the subleading jet based on pt and rapidity
    cmdline->end_section("Choice of input parameters"); //-----------

    // declare histograms 
    hists_2d_compact["lnovth_lnkt_wo_coll" ].declare(0.0, 2*logthmax,1.0, log(ktmin), log(ktmax), 0.5);   
    hists_2d_compact["lnovth_lnkt_w_coll"].declare(0.0, 2*logthmax, 1.0,  log(ktmin), log(ktmax), 0.5);    
  }

  //----------------------------------------------------------------------
  /// this gets called once for every event and should carry
  /// out the analysis and output histograms
  void user_analyse_event() override {
    double evwgt = event_weight();
   // Main loop over the PanScales event
   std::vector<PseudoJet> particles; 
    for (auto &p:f_event.particles()){
      PseudoJet q(p.px(),p.py(), p.pz(), p.E());
      if(abs(p.pdgid()) != 11) particles.push_back(q); // only final-state particles
    }
    // We cluster the event in C/A ycut=1 jets. 
    if(particles.size()<2) return;
    ClusterSequence cs(particles, _jet_def);

    // Obtain the declusterings 
    std::vector<LundEEDeclustering> declusterings;
    declusterings = _lund_ee_gen.result(cs);
    // Record the Lund coordinates
    bool coll_emission = false; // bool to check whether there are emissions in the collinear region 
    std::vector<pair<double,double>> lund_coordinates;
    for (const auto & d : declusterings){
      pair<double,double> coords = d.lund_coordinates();
      lund_coordinates.push_back(coords);
      // Check whether there is an emission in the collinear region
      if (coords.first > logthmin && coords.second < logthmax) coll_emission=true;
    }
    // NOTE: I guess there is a smart way of avoiding this second loop
    for (const auto & l : lund_coordinates)
      if (coll_emission){
        hists_2d_compact["lnovth_lnkt_w_coll"].add_entry(l.first, l.second, evwgt); // There was at least one collinear emission
      }
      else hists_2d_compact["lnovth_lnkt_wo_coll"].add_entry(l.first, l.second, evwgt);
  }
};
  
//----------------------------------------------------------------------
/// a minimal main program
int main(int argc, char** argv) {
  CmdLine cmdline(argc,argv,true);
  ExampleEELund driver(&cmdline);
  driver.run();
}
