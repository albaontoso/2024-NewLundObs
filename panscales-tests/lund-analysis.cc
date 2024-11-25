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

class LundAnalysis : public  AnalysisFramework {
public:
  JetDefinition _jet_def;
  RecursiveLundEEGenerator _lund_ee_gen; 
  double ktmin,ktmax;
  double etamin,etamax; // define the collinear window
  double eps_min; // min epsilon value
  double nbins;
  /// ctor
  LundAnalysis(CmdLine * cmdline_in)
    : AnalysisFramework(cmdline_in) {
      // check that we are running with the correct process
      if (!dynamic_cast<ProcessZ2qq*>(f_process.get()) && !dynamic_cast<ProcessH2gg*>(f_process.get())){
        throw runtime_error("This analysis is only to be used with an e+e- process");
      }
    _jet_def = JetDefinition(ee_genkt_algorithm, 1.0, 0.0);
    }

  //----------------------------------------------------------------------
  /// options needed by the "user" analysis can be set here
  void user_startup() override {
    
    cmdline->start_section("Choice of input parameters"); //---------

    // -- Lund plane windows that define the observable
    double rts = cmdline->value<double>("-rts");
    ktmin      = cmdline->value("-kt-min", 1);
    ktmax      = cmdline->value("-kt-max", rts/2.);
    etamin     = cmdline->value("-eta-min", 5);
    etamax     = cmdline->value("-eta-max", 6);

    cmdline->end_section("Choice of input parameters"); //-----------

    // declare histograms 
    hists_2d_compact["eta_lnkt_wo_coll" ].declare(0.0, 2*etamax,10, log(ktmin), log(ktmax), 20);   
    hists_2d_compact["eta_lnkt_w_coll"].declare(0.0, 2*etamax,10, log(ktmin), log(ktmax), 20);
    hists_2d_compact["eta_lnkt_full"].declare(0.0, 2*etamax,10, log(ktmin), log(ktmax), 20); 
    // create the histogram for the fractal story
    const Binning binning(ktmin/rts,1,15);
    eps_min = ktmin/rts;
    nbins = binning.size();  
    cumul_hists_err["eps_area"].declare(binning);  
    // cross sections
    xsections["eta_lnkt_wo_coll"] = AverageAndError();
    xsections["eta_lnkt_w_coll"]  = AverageAndError();
    xsections["eta_lnkt_full"]	  = AverageAndError();
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
    // We cluster the event into C/A ycut=1 jets. 
    if(particles.size()<2) return;
    ClusterSequence cs(particles, _jet_def);

    // Obtain the primary declusterings 
    std::vector<LundEEDeclustering> declusterings;
    declusterings = _lund_ee_gen.result(cs);
    double Q = sqrt(cs.Q2()); 
    // we are gonna first initialize the histogram
    // as if no emissions were resolved 
    // i.e. epsmax = 1
    LundArea(eps_min,Q/Q,evwgt);
    // Record the Lund coordinates
    bool coll_emission = false; // bool to check whether there are emissions in the collinear region 
    std::vector<pair<double,double>> lund_coordinates;
    for (const auto & d : declusterings){
      pair<double,double> coords = d.lund_coordinates();
      lund_coordinates.push_back(coords);
      // Check whether there is an emission in the collinear region
      if (coords.first > etamin && coords.first < etamax) coll_emission=true;
    }
    // NOTE: I guess there is a smart way of avoiding this 2nd loop, but I'm not inspired today..
    for (const auto & l : lund_coordinates){
      // fill up the fractal dimension histo with the declusterings
      double ktovQ = exp(l.second)/Q;
      LundArea(eps_min,ktovQ,evwgt);
      hists_2d_compact["eta_lnkt_full"].add_entry(l.first, l.second, evwgt); // first fill the primary Lund plane
      xsections["eta_lnkt_full"] += evwgt;
      // for our exclusive Lund plane densities, only record emissions in the kt-window   
      if (l.second > log(ktmin) && l.second < log(ktmax)){
    	  if (coll_emission){
        	hists_2d_compact["eta_lnkt_w_coll"].add_entry(l.first, l.second, evwgt); // There was at least one collinear emission
        	xsections["eta_lnkt_w_coll"] += evwgt;
     	  }
        else {
          hists_2d_compact["eta_lnkt_wo_coll"].add_entry(l.first, l.second, evwgt);
          xsections["eta_lnkt_wo_coll"] += evwgt;
      	}
      }
    }
  }
  // This function takes as input a given value of 
  //  epsmax = ktmax/Q (ktmax is the kt of the declustering)
  //  epsmin = ktmin/Q
  // and fills up the histogram corresponding to the 
  // Lund plane area
  void LundArea(double epsmin, double epsmax, double weight){
    double eps_step=(epsmax-epsmin)/to_double(nbins);
    for(unsigned int j=0; j<nbins;j++){
      double eps_val = epsmin+eps_step*to_double(j);
      // note that for the area we need Q/kt 
      // so we invert the value of eps_val
      double area = pow2(1./eps_val)*eps_step;
      cumul_hists_err["eps_area"].add_entry(eps_val,area*weight);
    }
  }
};
  
//----------------------------------------------------------------------
/// a minimal main program
int main(int argc, char** argv) {
  CmdLine cmdline(argc,argv,true);
  LundAnalysis driver(&cmdline);
  driver.run();
}
