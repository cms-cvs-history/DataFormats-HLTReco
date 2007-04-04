#include <TROOT.h>
#include <TFile.h>
#include <THashList.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TAttText.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPDF.h>
#include <TTree.h>
#include <TBranch.h>
#include <TSystem.h>

#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>
#include <string>

#include "DataFormats/Provenance/interface/EventAuxiliary.h"
#include "DataFormats/HLTReco/interface/HLTPerformanceInfo.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

//--- Created by:  
//--- Bryan Dahmes (Bryan.Michael.Dahmes@cern.ch), January 2007

enum histoType {Single, Path, Module, ModuleInPath} ; 

bool useModule(HLTPerformanceInfo::Module module,
               std::vector<std::string> skip) {

    for (unsigned int i=0; i<skip.size(); i++)
        if (module.name() == skip.at(i)) return false ; 
    return true ; 
}

bool useEvent(int evt, std::vector<int> skip) {

    for (unsigned int i=0; i<skip.size(); i++)
        if (evt == skip.at(i)) return false ;
    return true ; 
}

int getNumberOfModules(HLTPerformanceInfo perf,
                       std::vector<std::string> modList) {

    int nMods = perf.numberOfModules() ;
    for (HLTPerformanceInfo::Modules::const_iterator modIter=perf.beginModules();
         modIter!=perf.endModules(); modIter++)
        if (!useModule(*modIter,modList)) nMods-- ;
    return nMods ; 
}

int getNumberOfPathModules(HLTPerformanceInfo::Path path,
                           std::vector<std::string> modList) {

    int nMods = path.numberOfModules() ; 
    for (HLTPerformanceInfo::Path::const_iterator modIter=path.begin();
         modIter!=path.end(); modIter++)
        if (!useModule(*modIter,modList)) nMods-- ;
    return nMods ; 
}

bool firstSighting(HLTPerformanceInfo::Path::const_iterator modIter,
                   HLTPerformanceInfo::Path path) {

    for (HLTPerformanceInfo::Path::const_iterator modIter2=path.begin();
         modIter2!=modIter; modIter2++)
        if (modIter2->name() == modIter->name()) return false ;
    return true ;
}

double modifyPathTime(HLTPerformanceInfo::Path path,
                      std::vector<std::string> modList) {

    double newTime = path.time() ;
    for (HLTPerformanceInfo::Path::const_iterator modIter=path.begin();
         modIter!=path.end(); modIter++) {
        if (!useModule(*modIter,modList)) {
            //--- Check that the module runs in the path ---//
            if (modIter->indexInPath(path) <= int(path.status().index()))
                newTime -= modIter->time() ;
        }
    }
    return newTime ;
}

void plot1D(TH1D* histo, TCanvas* canvas) {

    double defaultSize = 0.04 ; 
    double defaultMargin = 0.10 ;
    canvas->SetBottomMargin(defaultMargin) ; 
    gStyle->SetLabelSize(defaultSize,"x") ; 
    double newSize = 0.5 / double(histo->GetNbinsX()) ;

    THashList* hLabels = histo->GetXaxis()->GetLabels() ; 
    if ((defaultSize > newSize)&&(hLabels!=0)) {
        gStyle->SetLabelSize(newSize,"x") ; 
        histo->LabelsOption("v","x") ;
        //--- Loop through the bins to get the largest bin name size ---//
        double maxSize = defaultMargin ; 
        for (int i=1; i<=histo->GetNbinsX(); i++) {
            std::string label = histo->GetXaxis()->GetBinLabel(i) ;
            double labelSize = 0.33 * label.size() * defaultSize ;
            if (labelSize > maxSize) maxSize = labelSize ; 
        }
        canvas->SetBottomMargin(maxSize) ; 
    }
        
    histo->Draw() ; canvas->Update() ;
}

void plot1D(TH1D* h1, TH1D* h2, TCanvas* canvas) {

    double defaultSize = 0.04 ; 
    double defaultMargin = 0.10 ;
    canvas->SetBottomMargin(defaultMargin) ; 
    gStyle->SetLabelSize(defaultSize,"x") ; 
    double newSize = 0.5 / double(h1->GetNbinsX()) ; 

    THashList* hLabels = h1->GetXaxis()->GetLabels() ; 
    if ((defaultSize > newSize)&&(hLabels!=0)) {
        gStyle->SetLabelSize(newSize,"x") ; 
        h1->LabelsOption("v","x") ;
        h2->LabelsOption("v","x") ;
        //--- Loop through the bins to get the largest bin name size ---//
        double maxSize = defaultMargin ; 
        for (int i=1; i<=h1->GetNbinsX(); i++) {
            std::string label = h1->GetXaxis()->GetBinLabel(i) ;
            double labelSize = 0.35 * label.size() * defaultSize ;
            if (labelSize > maxSize) maxSize = labelSize ; 
        }
        canvas->SetBottomMargin(maxSize) ; 
    }
        
    h1->Draw() ;
    h2->SetFillColor(2) ; h2->Draw("same") ; 
    canvas->Update() ;
}

void plotMany(std::vector<TH1D*> histo, TCanvas* canvas, histoType type,
              HLTPerformanceInfo perf, std::vector<std::string> skip) {
    
    if (type == Single) {
        std::cout << "Error in plotMany.  Trying to plot single histogram" << std::endl ; 
    } else {
        int ctr = 0 ;
        if (type == Path)
            for (HLTPerformanceInfo::PathList::const_iterator pathIter=perf.beginPaths();
                 pathIter!=perf.endPaths(); pathIter++)
                plot1D(histo[ctr++],canvas) ;
        if (type == Module)
            for (HLTPerformanceInfo::Modules::const_iterator modIter=perf.beginModules();
                 modIter!=perf.endModules(); modIter++) {
                if (useModule(*modIter,skip)) plot1D(histo[ctr],canvas) ;
                ctr++ ; 
            }
    }
}

void plotMany(std::vector<TH1D*> h1, std::vector<TH1D*> h2, TCanvas* canvas, histoType type,
              HLTPerformanceInfo perf, std::vector<std::string> skip) {
    
    if (type == Single) {
        std::cout << "Error in plotMany.  Trying to plot single histogram" << std::endl ; 
    } else {
        int ctr = 0 ;
        if (type == Path)
            for (HLTPerformanceInfo::PathList::const_iterator pathIter=perf.beginPaths();
                 pathIter!=perf.endPaths(); pathIter++) {
                plot1D(h1[ctr],h2[ctr],canvas) ;
                ctr++ ;
            }
        if (type == Module)
            for (HLTPerformanceInfo::Modules::const_iterator modIter=perf.beginModules();
                 modIter!=perf.endModules(); modIter++) {
                if (useModule(*modIter,skip)) plot1D(h1[ctr],h2[ctr],canvas) ;
                ctr++ ; 
            }
    }
}

void plotModuleInPath(std::vector< std::vector<TH1D*> > histo, TCanvas* canvas,
                      HLTPerformanceInfo perf, std::vector<std::string> skip) {

    int pCtr = 0 ;
    for (HLTPerformanceInfo::PathList::const_iterator pathIter=perf.beginPaths();
         pathIter!=perf.endPaths(); pathIter++) {
        int mCtr = 0 ; 
        for (HLTPerformanceInfo::Path::const_iterator modIter=pathIter->begin();
             modIter!=pathIter->end(); modIter++) {
            if (useModule(*modIter,skip)) 
                plot1D(histo.at(pCtr).at(mCtr),canvas) ;
            mCtr++ ; 
        }
        pCtr++ ; 
    }
}

void createTOC(TPDF* pdf, std::vector<std::string> TOC) {
    
    double textSize = 0.032 ; 
    double textSpacing = 0.040 ;
    double tocXval = 0.10 ; 
    double tocYval = 0.93 ; 
    double pageXval = 0.90 ; 

    pdf->SetTextSize(textSize) ;
    pdf->SetTextColor(1) ; 

    int nTOCentries = 0 ; 
    std::string page ; std::string TOCentry ; 
    for (unsigned int i=0; i<TOC.size(); i++) {
        std::string TOCline = TOC.at(i) ;
        int locPage = TOCline.find("^^",0) ; 
        int locSub  = TOCline.find("^",0) ; 
        if (locSub == locPage) {
            nTOCentries++ ;
            page = TOCline.substr(locPage+2) ;
            TOCentry = TOCline.substr(0,locPage) ; 
        }
    }

    tocYval -= textSpacing * nTOCentries ;
    pdf->Text(tocXval,tocYval,TOCentry.c_str()) ; 
    pdf->Text(pageXval,tocYval,page.c_str()) ; 
}

int main(int argc, char ** argv) {
    
  //-- Load libraries ---//
  gSystem->Load("libFWCoreFWLite") ;
  AutoLibraryLoader::enable() ;

  //--- Default arguments ---//
  std::string filename = "hlt.root" ;
  std::string outbase  = "hltTimingSummary" ; 
  std::string outname  = outbase + ".root" ;
  std::string pdfname  = outbase + ".pdf" ;
  std::string txtname  = outbase + ".txt" ;
  std::vector<std::string> skipTiming ; skipTiming.clear() ; 
  std::string excludeName ; 

  double userMaxTime = -1. ;
  bool skipFirstEvent = false ; 
  
  //--- Get parameters from command line ---//
  boost::program_options::options_description desc(
      "Available options for hltTimingSummary") ; 
  desc.add_options()
      ("help,h","Print this help message")
      ("infile,i",   boost::program_options::value<std::string>(),
       "Input file name (Default is hlt.root)") 
      ("outbase,o",  boost::program_options::value<std::string>(),
       "Default base name for output .root, .pdf, and .txt files (Default is hltTimingSummary)")
      ("rootfile,r", boost::program_options::value<std::string>(),
       "Output root file name")
      ("bookmark,b", boost::program_options::value<std::string>(),
       "Output bookmark file name")
      ("pdf,p",      boost::program_options::value<std::string>(),
       "Output pdf file name")
      ("time,t",     boost::program_options::value<double>(),
       "All relevant histogram time axes run from 0 to the user-specified value (in msec)")
      ("noFirst,f","Skip ANY event where a module is run for the first time (Default is to run over all events)") 
      ("exclude,e",   boost::program_options::value<std::string>(),
       "Exclude a list of modules from the timing calculation") ;


  std::string usage = "\nSample hltTimingSummary usage::\n" ; 
  usage += "\"hltTimingSummary -t 50\" " ; 
  usage += "inputs hlt.root, outputs hltTimingsummary.root, .txt, and .pdf\n" ; 
  usage += "                         and timing histograms run from 0 to 50 msec\n" ; 
  usage += "\"hltTimingSummary -i input.root -o output\" " ;
  usage += "inputs input.root, outputs output.root, .txt, and .pdf\n" ; 
  usage += "\"hltTimingSummary -o output -b special.txt\" " ;
  usage += "inputs hlt.root, outputs output.root and .pdf, special.txt\n\n" ; 
  usage += "\"hltTimingSummary -f -e exclude.txt\" " ;
  usage += "inputs hlt.root, outputs hltTimingsummary.root, .txt, and .pdf.\n" ;
  usage += "                                     Also skips events where modules are first run, and excludes \n" ;
  usage += "                                     the modules listed in exclude.txt from timing calculations.\n\n" ;
  usage += "NOTE: To exclude files, the following formats are acceptable:\n" ; 
  usage += "      - A comma-separated list of modules (e.g. \"hltTimingSummary -e a,b,c\")\n" ; 
  usage += "      - A text file, where each excluded module appears on its own line\n" ;
  usage += "        (e.g. \"hltTimingSummary -e file.txt\")\n" ; 

  boost::program_options::positional_options_description pos ; 
  boost::program_options::variables_map vmap ;

  try {
      boost::program_options::store(boost::program_options::command_line_parser(argc,argv).
                                    options(desc).positional(pos).run(), vmap) ; 
  } catch (boost::program_options::error const& x) {
      std::cerr << "Unable to parse options:\n"
                << x.what() << "\n\n" ;
      std::cerr << desc << usage << std::endl ;
      return 1 ; 
  }
  
  boost::program_options::notify(vmap) ; 
  if (vmap.count("help")) {
      std::cout << desc << usage <<  std::endl ;
      return 1 ;
  }
  if (vmap.count("infile")) {
      filename = vmap["infile"].as<std::string>() ; 
  }
  if (vmap.count("outbase")) {
      outbase = vmap["outbase"].as<std::string>() ; 
      outname = outbase + ".root" ; 
      pdfname = outbase + ".pdf" ; 
      txtname = outbase + ".txt" ; 
  }
  if (vmap.count("rootfile")) {
      outname = vmap["rootfile"].as<std::string>() ; 
  }
  if (vmap.count("pdf")) {
      pdfname = vmap["pdf"].as<std::string>() ; 
  }
  if (vmap.count("bookmark")) {
      txtname = vmap["bookmark"].as<std::string>() ; 
  }
  if (vmap.count("time")) {
      userMaxTime = vmap["time"].as<double>() ; 
  }
  if (vmap.count("noFirst")) {
      skipFirstEvent = true ; 
  }
  if (vmap.count("exclude")) {
      excludeName = vmap["exclude"].as<std::string>() ; 

      ifstream excludeFile(excludeName.c_str()) ;
      if (excludeFile.is_open()) { //--- Excluded modules listed in a file ---//
          while ( !excludeFile.eof() ) {
              std::string skipped ;
              getline(excludeFile,skipped) ; 
              skipTiming.push_back( skipped ) ;
          }
      } else { //--- Assume the file is a comma-separated list of modules ---//
          unsigned int strStart = 0 ; 
          for (unsigned int itr=excludeName.find(",",0); itr!=std::string::npos;
               itr=excludeName.find(",",itr)) {
              std::string skipped = excludeName.substr(strStart,(itr-strStart)) ; 
              itr++ ; strStart = itr ; 
              skipTiming.push_back( skipped ) ;
          }
          //--- Fill the last entry ---//
          skipTiming.push_back( excludeName.substr(strStart,excludeName.length()) ) ; 
      }
  }

  std::cout << "Opening file " << filename << std::endl ;
  TFile file(filename.c_str());
  if (file.IsZombie()) {
      std::cout << "*** Error opening file: " << filename << " ***" << std::endl;
      std::cout << "\n\n" << desc << usage <<  std::endl ;
      return 1 ;
  }

  TTree * events = dynamic_cast<TTree *>(file.Get("Events") );
  assert(events);

  //--- Find the HLTPerformanceInfo branch ---//
  TIter iter(events->GetListOfBranches()) ;
  std::string hltPerfInfoBranchName ; 
  TBranch* cand ;
  TBranch::ResetCount() ;
  while ( (cand = (TBranch*)iter()) ) {
      std::string branchName = cand->GetName() ;
      unsigned int loc = branchName.find( "HLTPerformanceInfo" ) ;
      if ( loc != std::string::npos ) hltPerfInfoBranchName = branchName + "obj" ;
  }

  TBranch* TBPerfInfo = events->GetBranch( hltPerfInfoBranchName.c_str() );
  assert(TBPerfInfo);
  HLTPerformanceInfo HLTPerformance ;
  TBPerfInfo->SetAddress((void *) & HLTPerformance) ;

  //--- Prepare the output ---//
  TFile* outFile = new TFile(outname.c_str(), "recreate") ;
  std::cout << "Output to files " << outname 
            << ", " << pdfname << ", and " << txtname << std::endl ;
  ofstream txtfile ;
  txtfile.open(txtname.c_str()) ; 
          
  int n_evts = events->GetEntries();
  std::vector<int> pathSuccess ; 
  std::vector<double> pathTimer ; 
  std::vector<double> incPathTimer ; 
  std::vector<double> moduleTimer ; 
  std::vector<double> moduleInitDouble ; 
  std::vector<int> moduleInitInt ; 
  std::vector< std::vector<double> > moduleInPathTimer ; 
  std::vector< std::vector<int> > failures ; 
  std::vector<std::string> tocList ; 
  std::vector<bool> successPerEvent ;
  std::vector<int> uniquePathSuccess ; 
  std::vector< std::vector<int> > successPathVsPath ; 
  std::vector<int> moduleRunStats ; 
  std::vector< std::vector<int> > moduleInPathRunStats ; 

  //--- Events to be skipped ---//
  std::vector<int> skipEvents ; 
  
  //--- Histogram initialization ---//
  std::vector<TH1D*> pathTime ;
  std::vector<TH1D*> moduleTime ; 
  std::vector<TH1D*> moduleScaledTime ; 
  std::vector< std::vector<TH1D*> > moduleInPathScaledTime ;  
  std::vector<TH1D*> incPathTime ; 
  std::vector<TH1D*> moduleInPathTimeSummary ; 
  std::vector<TH1D*> moduleInPathScaledTimeSummary ; 
  std::vector<TH1D*> failedModule ; 
  std::vector<TH1D*> moduleInPathRejection ; 
  std::vector<TH1D*> moduleInPathRejectAll ; 
  std::vector<TH1D*> moduleInPathRejectTime ; 
  std::vector<TH1D*> moduleInitTH1D ; 

  //--- Initial loop to get scale for histograms ---//
  double maxTime = 0. ;
  double xmin = 0. ;
  double xmax = 0. ;
  int nSkippedMods = skipTiming.size() ; 
  int nSkippedEvts = 0 ; 

  //--- Create the list of skipped events ---//
  if (skipFirstEvent) {
      for (int ievt=0; ievt<n_evts; ievt++) {
          TBPerfInfo->GetEntry(ievt) ;
          
          int mCtr = 0 ; 
          for (HLTPerformanceInfo::Modules::const_iterator modIter=HLTPerformance.beginModules();
               modIter!=HLTPerformance.endModules(); modIter++) {
              if (skipEvents.empty())
                  for (unsigned int i=0; i<HLTPerformance.numberOfModules(); i++) skipEvents.push_back(-1) ; 
              if (modIter->time() > 0) {
                  if (skipEvents.at(mCtr) < 0) {
                  skipEvents.at(mCtr) = ievt ; 
                  }
              }
              mCtr++ ; 
          }
      }
      //--- Determine the number of events to be skipped ---//
      for (unsigned int i=0; i<skipEvents.size(); i++) {
          bool duplicateSkip = false ;
          for (unsigned int j=0; j<i; j++) {
              if (skipEvents.at(j) == skipEvents.at(i)) duplicateSkip = true ;
          }
          if ((!duplicateSkip) && (skipEvents.at(i)>=0)) nSkippedEvts++ ; 
      }
      std::cout << "Skipping a total of " << nSkippedEvts
                << " events for module initialization" << std::endl ;
  }  

  //--- We need to account for the excluded modules when we calculated total time ---//
  bool initMessage = false ; bool initWarning = false ; 
  for (int ievt=0; ievt<n_evts; ievt++) {

      if (!useEvent(ievt,skipEvents)) continue ;
      TBPerfInfo->GetEntry(ievt) ; 

      double time = HLTPerformance.totalTime() ;
      if (skipTiming.size() > 0) {  //--- We have skipped modules ---//
          if ( !initMessage ) {
              initMessage = true ; 
              std::cout << "Excluding the following module(s): " << std::endl ;
              for (unsigned int i=0; i<skipTiming.size(); i++)
                  std::cout << skipTiming.at(i) << std::endl ;
          }
          for (unsigned int i=0; i<skipTiming.size(); i++) {
              HLTPerformanceInfo::Modules::const_iterator mod =
                  HLTPerformance.findModule( skipTiming.at(i).c_str() ) ; 
              if ( mod != HLTPerformance.endModules() ) {
                  time -= mod->time() ;
              } else {
                  if (!initWarning) {
                      initWarning = true ;
                      std::cout << "Module \"" << skipTiming.at(i).c_str()
                                << "\" is not found in any path.  Are you sure"
                                << " you have the right exclusion list?" << std::endl ; 
                      nSkippedMods-- ; 
                  }
              }
          }
      }
      if ( time > maxTime ) maxTime = time ; 
  }

  if (userMaxTime > 0) {
      xmax = userMaxTime ;
  } else {
      int xscale = 4 ;
      if (maxTime == 0) {
          xscale = -2 ; 
      } else {
          while (pow(10,xscale--) > maxTime) ; 
      }
      xscale += 2 ; 
      xmax = pow(10,xscale+3) ;
      if ( (xmax/5.) > (maxTime * 1000.) ) xmax /= 5. ; 
      if ( (xmax/2.) > (maxTime * 1000.) ) xmax /= 2. ; 
  }

  std::cout << "For timing histograms, x-axis ranges from 0 to " << xmax << " msec" << std::endl ; 

  TH1D* totalTime = new TH1D("totalTime","Total time for all modules per event",100,xmin,xmax) ;
  TH1D* pathTimeSummary =
      new TH1D("pathTimeSummary","Average time per path",
               HLTPerformance.numberOfPaths(),0.,double(HLTPerformance.numberOfPaths())) ; 
  TH1D* incPathTimeSummary =
      new TH1D("incPathTimeSummary","Average incremental time per path",
               HLTPerformance.numberOfPaths(),0.,double(HLTPerformance.numberOfPaths())) ; 
  TH1D* pathSuccessFraction =
      new TH1D("pathSuccessFraction","Path success rate (%)",
               HLTPerformance.numberOfPaths(),0.,double(HLTPerformance.numberOfPaths())) ; 
  TH1D* uniquePathSuccessFraction =
      new TH1D("uniquePathSuccessFraction","Fraction (%) of events passing due to a single path",
               HLTPerformance.numberOfPaths(),0.,double(HLTPerformance.numberOfPaths())) ; 
  TH1D* pathRejection =
      new TH1D("pathRejection","Rejection for each path",
               HLTPerformance.numberOfPaths(),0.,double(HLTPerformance.numberOfPaths())) ; 
  TH1D* pathRejectAll =
      new TH1D("pathRejectAll","Rejection for each path",
               HLTPerformance.numberOfPaths(),0.,double(HLTPerformance.numberOfPaths())) ; 

  TH2D* pathVsPathSummary = new TH2D("pathVsPathSummary", "Relative path success",
                               HLTPerformance.numberOfPaths(),
                               0.,double(HLTPerformance.numberOfPaths()), 
                               HLTPerformance.numberOfPaths(),
                               0.,double(HLTPerformance.numberOfPaths())) ; 

  totalTime->GetXaxis()->SetTitle("msec") ; 
  pathTimeSummary->GetYaxis()->SetTitle("msec") ;
  incPathTimeSummary->GetYaxis()->SetTitle("msec") ;

  char name[100] ; 
  char title[200] ; 
  int pCtr = 0 ; 
  int mCtr = 0 ; 
  for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
       pathIter!=HLTPerformance.endPaths(); pathIter++) {

      TH1D* histo ; 
              
      pathTimer.push_back(0.) ;
      incPathTimer.push_back(0.) ;
      pathSuccess.push_back(0) ;
      uniquePathSuccess.push_back(0) ;
      successPerEvent.push_back(false) ; 

      pathTimeSummary->GetXaxis()->SetBinLabel(pCtr+1,pathIter->name().c_str()) ;
      incPathTimeSummary->GetXaxis()->SetBinLabel(pCtr+1,pathIter->name().data()) ; 
      pathSuccessFraction->GetXaxis()->SetBinLabel(pCtr+1,pathIter->name().data()) ; 
      uniquePathSuccessFraction->GetXaxis()->SetBinLabel(pCtr+1,pathIter->name().data()) ; 
      pathRejection->GetXaxis()->SetBinLabel(pCtr+1,pathIter->name().data()) ;
      pathRejectAll->GetXaxis()->SetBinLabel(pCtr+1,pathIter->name().data()) ;

      pathVsPathSummary->GetXaxis()->SetBinLabel(pCtr+1,pathIter->name().data()) ; 
      pathVsPathSummary->GetYaxis()->SetBinLabel(pCtr+1,pathIter->name().data()) ; 

      sprintf(name,"pathTime_%s",pathIter->name().c_str()) ;
      sprintf(title,"Per event time for path %s",pathIter->name().c_str()) ;
      histo = new TH1D(name,title,100,xmin,xmax) ;
      histo->GetXaxis()->SetTitle("msec") ;
      pathTime.push_back( histo ) ; 

      sprintf(name,"incPathTime_%s",pathIter->name().c_str()) ;
      sprintf(title,"Per event incremental time for path %s",pathIter->name().c_str()) ;
      histo = new TH1D(name,title,100,xmin,xmax) ; 
      histo->GetXaxis()->SetTitle("msec") ; 
      incPathTime.push_back( histo ) ; 

      sprintf(name,"failedModule_%s",pathIter->name().c_str()) ;
      sprintf(title,"Failure fraction (%%) by module for path %s",pathIter->name().c_str()) ;
      histo = new TH1D(name,title,
                       (1+getNumberOfPathModules(*pathIter,skipTiming)),-1.,
                       double(getNumberOfPathModules(*pathIter,skipTiming))) ; 
      histo->GetXaxis()->SetBinLabel(1,"SUCCESS") ; 
      failedModule.push_back( histo ) ; 

      sprintf(name,"moduleInPathTimeSummary_%s",pathIter->name().c_str()) ;
      sprintf(title,"Average module time for path %s",
              pathIter->name().c_str()) ;
      histo = new TH1D(name,title,
                       (getNumberOfPathModules(*pathIter,skipTiming)),0.,
                       double(getNumberOfPathModules(*pathIter,skipTiming))) ; 
      histo->GetYaxis()->SetTitle("msec") ; 
      moduleInPathTimeSummary.push_back( histo ) ; 
      
      sprintf(name,"moduleInPathScaledTimeSummary_%s",pathIter->name().c_str()) ;
      sprintf(title,"Average module RUNNING time for path %s",
              pathIter->name().c_str()) ;
      histo = new TH1D(name,title,
                       (getNumberOfPathModules(*pathIter,skipTiming)),0.,
                       double(getNumberOfPathModules(*pathIter,skipTiming))) ; 
      histo->GetYaxis()->SetTitle("msec") ; 
      moduleInPathScaledTimeSummary.push_back( histo ) ; 
      
      sprintf(name,"moduleInPathRejection_%s",pathIter->name().c_str()) ;
      sprintf(title,"Rejection per module for path %s",pathIter->name().c_str()) ;
      histo = new TH1D(name,title,
                       (getNumberOfPathModules(*pathIter,skipTiming)),0.,
                       double(getNumberOfPathModules(*pathIter,skipTiming))) ; 
      moduleInPathRejection.push_back( histo ) ; 

      sprintf(name,"moduleInPathRejectAll_%s",pathIter->name().c_str()) ;
      sprintf(title,"Full rejection per module for path %s",pathIter->name().c_str()) ;
      histo = new TH1D(name,title,
                       (getNumberOfPathModules(*pathIter,skipTiming)),0.,
                       double(getNumberOfPathModules(*pathIter,skipTiming))) ; 
      moduleInPathRejectAll.push_back( histo ) ; 

      sprintf(name,"moduleInPathRejectTime_%s",pathIter->name().c_str()) ;
      sprintf(title,"(Rejection / average running time) per module for path %s",pathIter->name().c_str()) ;
      histo = new TH1D(name,title,
                       (getNumberOfPathModules(*pathIter,skipTiming)),0.,
                       double(getNumberOfPathModules(*pathIter,skipTiming))) ; 
      histo->GetYaxis()->SetTitle("1/msec") ; 
      moduleInPathRejectTime.push_back( histo ) ; 

      mCtr = 0 ;
      int mIdx = 1 ;

      for (HLTPerformanceInfo::Path::const_iterator modIter=pathIter->begin();
           modIter!=pathIter->end(); modIter++) {

          sprintf(name,"moduleInPathScaledTime_%s_%s",
                  pathIter->name().c_str(),modIter->name().c_str()) ;
          sprintf(title,"Per event RUNNING time for module %s in path %s",
                  modIter->name().c_str(),pathIter->name().c_str()) ;
          histo = new TH1D(name,title,100,xmin,xmax) ;
          histo->GetXaxis()->SetTitle("msec") ;
          moduleInitTH1D.push_back( histo ) ; 
          
          moduleInitDouble.push_back(0.) ; 
          moduleInitInt.push_back(0) ;

          const char* modName = modIter->name().c_str() ;
          if (useModule(*modIter,skipTiming)) {
              moduleInPathTimeSummary[pCtr]->GetXaxis()->SetBinLabel(mIdx,modName) ;
              moduleInPathScaledTimeSummary[pCtr]->GetXaxis()->SetBinLabel(mIdx,modName) ;
              moduleInPathRejection[pCtr]->GetXaxis()->SetBinLabel(mIdx,modName) ;
              moduleInPathRejectAll[pCtr]->GetXaxis()->SetBinLabel(mIdx,modName) ;
              moduleInPathRejectTime[pCtr]->GetXaxis()->SetBinLabel(mIdx,modName) ;
              failedModule[pCtr]->GetXaxis()->SetBinLabel(++mIdx,modName) ;
          }
          mCtr++ ; 
      }

      pCtr++ ;
      moduleInPathRunStats.push_back( moduleInitInt ) ; 
      moduleInPathTimer.push_back( moduleInitDouble ) ;
      failures.push_back( moduleInitInt ) ;
      moduleInPathScaledTime.push_back( moduleInitTH1D ) ; 
      moduleInitDouble.clear() ; 
      moduleInitInt.clear() ;
      moduleInitTH1D.clear() ; 
  }

  for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
       pathIter!=HLTPerformance.endPaths(); pathIter++) {
      successPathVsPath.push_back( pathSuccess ) ; 
  }
  
  TH1D* moduleTimeSummary =
      new TH1D("moduleTimeSummary","Average time per module",
               getNumberOfModules(HLTPerformance,skipTiming),0.,
               double(getNumberOfModules(HLTPerformance,skipTiming))) ;
  moduleTimeSummary->GetYaxis()->SetTitle("msec") ; 
  TH1D* moduleScaledTimeSummary =
      new TH1D("moduleScaledTimeSummary","Average RUNNING time per module",
               getNumberOfModules(HLTPerformance,skipTiming),0.,
               double(getNumberOfModules(HLTPerformance,skipTiming))) ;
  moduleScaledTimeSummary->GetYaxis()->SetTitle("msec") ; 

  mCtr = 0 ;
  int mIdx = 1 ; 
  for (HLTPerformanceInfo::Modules::const_iterator modIter=HLTPerformance.beginModules();
       modIter!=HLTPerformance.endModules(); modIter++) {
      moduleTimer.push_back(0.) ;
      moduleRunStats.push_back(0) ; 
      if (useModule(*modIter,skipTiming)) { 
          moduleTimeSummary->GetXaxis()->SetBinLabel(mIdx,modIter->name().c_str()) ; 
          moduleScaledTimeSummary->GetXaxis()->SetBinLabel(mIdx++,modIter->name().c_str()) ; 
      }

      TH1D* histo ; 
      
      sprintf(name,"moduleTime_%s",modIter->name().c_str()) ;
      sprintf(title,"Time per event for module %s",modIter->name().c_str()) ;
      histo = new TH1D(name,title,100,xmin,xmax) ;
      histo->GetXaxis()->SetTitle("msec") ;
      moduleTime.push_back( histo ) ; 

      sprintf(name,"moduleScaledTime_%s",modIter->name().c_str()) ;
      sprintf(title,"RUNNING time per event for module %s",modIter->name().c_str()) ;
      histo = new TH1D(name,title,100,xmin,xmax) ;
      histo->GetXaxis()->SetTitle("msec") ;
      moduleScaledTime.push_back( histo ) ; 

  }

  for (int ievt=0; ievt<n_evts; ievt++) {
      
      if (!useEvent(ievt,skipEvents)) continue ;
      TBPerfInfo->GetEntry(ievt) ;

      pCtr = 0 ;
      for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
           pathIter!=HLTPerformance.endPaths(); pathIter++) {

          pathTimer.at(pCtr) += modifyPathTime(*pathIter,skipTiming) ;

          double addedPathTime = 0 ; mCtr = 0 ; 
          for (HLTPerformanceInfo::Path::const_iterator modIter=pathIter->begin();
               modIter!=pathIter->end(); modIter++) {

              if ( modIter->indexInPath(*pathIter) <= int(pathIter->status().index()) ) {
                  moduleInPathRunStats.at(pCtr).at(mCtr)++ ;
                  moduleInPathTimer.at(pCtr).at(mCtr) += modIter->time() ;  
                  moduleInPathScaledTime.at(pCtr).at(mCtr)->Fill(1000. * modIter->time()) ; 
              }
              if ( useModule(*modIter,skipTiming) &&
                   (modIter->indexInPath(*pathIter) <= int(pathIter->status().index())) &&
                   (HLTPerformance.uniqueModule(modIter->name().c_str())) ) {
                  addedPathTime += modIter->time() ; 
              }
              mCtr++ ; 
          }

          incPathTimer.at(pCtr) += addedPathTime ;
          pathTime.at(pCtr)->Fill( 1000. * modifyPathTime(*pathIter,skipTiming) ) ; 
          incPathTime.at(pCtr)->Fill( 1000. * addedPathTime ) ; 

          if (pathIter->status().accept()) {
              pathSuccess.at(pCtr)++ ;
              successPerEvent.at(pCtr) = true ;
          } else { //--- One of the modules caused the path to fail ---//
              successPerEvent.at(pCtr) = false ;
              failures.at(pCtr).at(pathIter->status().index())++ ;
          }
          pCtr++ ;
      }

      for (unsigned int i=0; i<pathTimer.size(); i++) {
          bool uniqueSuccess = successPerEvent.at(i) ; 
          for (unsigned int j=0; j<pathTimer.size(); j++) {
              if ( successPerEvent.at(i) && successPerEvent.at(j) ) {
                  if (i != j) uniqueSuccess = false ;
                  successPathVsPath.at(i).at(j)++ ; 
              }
          }
          if ( successPerEvent.at(i) && uniqueSuccess ) uniquePathSuccess.at(i)++ ; 
      }

      mCtr = 0 ;
      for (HLTPerformanceInfo::Modules::const_iterator modIter=HLTPerformance.beginModules();
           modIter!=HLTPerformance.endModules(); modIter++) {
          moduleTime.at(mCtr)->Fill( 1000. * modIter->time() ) ;

          if (modIter->status().wasrun()) {
              moduleScaledTime.at(mCtr)->Fill( 1000. * modIter->time() ) ;
              moduleRunStats.at(mCtr)++ ;
          }
          moduleTimer.at(mCtr++) += modIter->time() ; 
      }

      double modifiedTotalTime = HLTPerformance.totalTime() ;
      for (unsigned int i=0; i<skipTiming.size(); i++) {
          HLTPerformanceInfo::Modules::const_iterator mod =
              HLTPerformance.findModule( skipTiming.at(i).c_str() ) ; 
          if ( mod != HLTPerformance.endModules() ) modifiedTotalTime -= mod->time() ;
      }
      totalTime->Fill( 1000. * modifiedTotalTime ) ;
  }

  pCtr = 0 ;
  for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
       pathIter!=HLTPerformance.endPaths(); pathIter++) {

      int totalFailures = 0 ;
      int evtsLeft = n_evts - nSkippedEvts ; 
      int mIdx = 0 ; mCtr = 0 ; 
      for (HLTPerformanceInfo::Path::const_iterator modIter=pathIter->begin();
           modIter!=pathIter->end(); modIter++) {

          //--- Times are in msec ---//
          double scaledAvgModuleTime = -1 ; 
          double avgModuleTime = 0. ; 
          if (moduleInPathRunStats.at(pCtr).at(mCtr) > 0)
              scaledAvgModuleTime = 1000. * moduleInPathTimer.at(pCtr).at(mCtr) /
                  double(moduleInPathRunStats.at(pCtr).at(mCtr)) ; 
          avgModuleTime = 1000. * moduleInPathTimer.at(pCtr).at(mCtr) / double(n_evts-nSkippedEvts) ; 

          //--- Calculate the rejection factor for each module in path, ---//
          //--- defined as N(HLT input)/N(HLT accept) ---//
          double modReject = 0 ;
          if (failures.at(pCtr).at(mCtr) < evtsLeft) {
              modReject = double(evtsLeft) / double(evtsLeft - failures.at(pCtr).at(mCtr)) ; 
          } else {
              modReject = double(evtsLeft) ;
              moduleInPathRejectAll.at(pCtr)->Fill( double(mIdx), modReject ) ; 
          }
          if (useModule(*modIter,skipTiming))
              moduleInPathRejection.at(pCtr)->Fill( double(mIdx), modReject ) ; 

          //--- Calculate rejection per unit time for each module in path ---//
          if (useModule(*modIter,skipTiming)) {
              moduleInPathTimeSummary.at(pCtr)->Fill( double(mIdx), avgModuleTime ) ; 
              moduleInPathScaledTimeSummary.at(pCtr)->Fill( double(mIdx), scaledAvgModuleTime ) ; 
              moduleInPathRejectTime.at(pCtr)->Fill( double(mIdx),
                                                     modReject / scaledAvgModuleTime ) ; 
          }
          
          evtsLeft -= failures.at(pCtr).at(mCtr) ; 
          totalFailures += failures.at(pCtr).at(mCtr) ; 

          if (failures.at(pCtr).at(mCtr) >= 0) {
              if (useModule(*modIter,skipTiming))
                  failedModule.at(pCtr)->Fill( double(mIdx),
                                               100. * double(failures.at(pCtr).at(mCtr)) /
                                               double(n_evts - nSkippedEvts) );
          }
          if (useModule(*modIter,skipTiming)) mIdx++ ;
          mCtr++ ;
      }
      failedModule.at(pCtr)->Fill(-1.,100. * (double(n_evts-nSkippedEvts-totalFailures)/double(n_evts-nSkippedEvts))) ; 

      double pathReject = n_evts - nSkippedEvts ;
      if ((n_evts-nSkippedEvts) > totalFailures) {
          pathReject = double(n_evts-nSkippedEvts) / double(n_evts - nSkippedEvts - totalFailures) ;
      } else {
          pathRejectAll->Fill(double(pCtr), pathReject) ; 
      }
      pathRejection->Fill(double(pCtr), pathReject) ; 
      pCtr++ ; 
  }
  
  for (unsigned int i=0; i<pathTimer.size(); i++) {
      pathTimeSummary->Fill(double(i), 1000. * (pathTimer.at(i)/(n_evts-nSkippedEvts))) ;
      incPathTimeSummary->Fill(double(i), 1000. * (incPathTimer.at(i)/(n_evts-nSkippedEvts))) ;
      pathSuccessFraction->Fill(double(i), 100. * double(pathSuccess.at(i))/(n_evts-nSkippedEvts)) ;
      uniquePathSuccessFraction->Fill(double(i),
                                      100. * double(uniquePathSuccess.at(i))/(n_evts-nSkippedEvts)) ; 
      
      for (unsigned int j=0; j<pathTimer.size(); j++) {
          if (pathSuccess.at(i) > 0) {
              pathVsPathSummary->Fill( double(i), double(j),
                                       double(successPathVsPath.at(i).at(j))/double(pathSuccess.at(i)) ) ; 
          } else {
              pathVsPathSummary->Fill( double(i), double(j), 0.) ; 
          }
      }
  }

  mCtr = 0 ; mIdx = 0 ; 
  for (HLTPerformanceInfo::Modules::const_iterator modIter=HLTPerformance.beginModules();
       modIter!=HLTPerformance.endModules(); modIter++) {

      if (useModule(*modIter,skipTiming)) {
          moduleTimeSummary->Fill(double(mIdx), 1000. * (moduleTimer.at(mCtr)/(n_evts-nSkippedEvts))) ;
          if (moduleRunStats.at(mCtr) > 0) {
              moduleScaledTimeSummary->Fill(double(mIdx), 1000. *
                                            (moduleTimer.at(mCtr)/double(moduleRunStats.at(mCtr)))) ;
          } else {
              moduleScaledTimeSummary->Fill(double(mIdx),-1.) ; 
          }
          mIdx++ ; 
      }
      mCtr++ ; 
  }

  //--- Dump results ---//
  {
      TCanvas* c1 = new TCanvas("c1") ;
      pathTimeSummary->SetMinimum(0.) ; 
      incPathTimeSummary->SetMinimum(0.) ; 
      pathSuccessFraction->SetMinimum(0.) ; 
      pathRejection->SetMinimum(0.) ; 
      moduleTimeSummary->SetMinimum(0.) ; 
      moduleScaledTimeSummary->SetMinimum(0.) ; 

      pCtr = 0 ; 
      for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
           pathIter!=HLTPerformance.endPaths(); pathIter++) {
          failedModule.at(pCtr)->SetMinimum(0.) ; 
          moduleInPathTimeSummary.at(pCtr)->SetMinimum(0.) ; 
          moduleInPathScaledTimeSummary.at(pCtr)->SetMinimum(0.) ; 
          moduleInPathRejection.at(pCtr)->SetMinimum(0.) ; 
          moduleInPathRejectAll.at(pCtr)->SetMinimum(0.) ; 
          moduleInPathRejectTime.at(pCtr)->SetMinimum(0.) ; 
          pCtr++ ; 
      }

      TPDF* pdf = new TPDF(pdfname.c_str()) ; 
      int pageNumber = 2 ; 
      double titleSize = 0.050 ; 
      
      gROOT->SetStyle("Plain") ; 
      gStyle->SetPalette(1) ; 
      c1->UseCurrentStyle() ; 
      gROOT->ForceStyle() ;

      pdf->SetTextColor(1) ; 
      pdf->SetTextSize(titleSize) ; 
      pdf->Text(0.33,0.93,"Timing Summary Output") ;


      pdf->SetTextColor(2) ;
      pdf->SetTextSize(0.037) ;
      pdf->Text(0.05,0.02,
                "Documentation at https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideHLTTimingSummary") ;
      
      std::string tocEntry = "Total time per event" ; 
      char tocPage[5] ; sprintf(tocPage,"%d",pageNumber++) ;
      tocList.push_back( tocEntry + "^^" + std::string(tocPage) ) ; 
      createTOC(pdf,tocList) ; 
      tocEntry = "Per event path time" ;
      sprintf(tocPage,"%d",pageNumber) ; 
      tocList.push_back( tocEntry + "^^" + std::string(tocPage) ) ; 
      for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
           pathIter!=HLTPerformance.endPaths(); pathIter++) {
          std::string subEntry = pathIter->name() ;
          sprintf(tocPage,"%d",pageNumber++) ;
          tocList.push_back( tocEntry + "^" + subEntry + "^^" + std::string(tocPage) ) ; 
      }
      createTOC(pdf,tocList) ; 
      if (HLTPerformance.numberOfPaths() > 1) {
          tocEntry = "Per event incremental path time" ;
          sprintf(tocPage,"%d",pageNumber) ; 
          tocList.push_back( tocEntry + "^^" + std::string(tocPage) ) ; 
          for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
               pathIter!=HLTPerformance.endPaths(); pathIter++) {
              std::string subEntry = pathIter->name() ;
              sprintf(tocPage,"%d",pageNumber++) ;
              tocList.push_back( tocEntry + "^" + subEntry + "^^" + std::string(tocPage) ) ; 
          }
          createTOC(pdf,tocList) ; 
      }
      tocEntry = "Per event module time" ;
      sprintf(tocPage,"%d",pageNumber) ; 
      tocList.push_back( tocEntry + "^^" + std::string(tocPage) ) ; 
      for (HLTPerformanceInfo::Modules::const_iterator modIter=HLTPerformance.beginModules();
           modIter!=HLTPerformance.endModules(); modIter++) {
          if (!useModule(*modIter,skipTiming)) continue ;  
          std::string subEntry = modIter->name() ;
          sprintf(tocPage,"%d",pageNumber++) ;
          tocList.push_back( tocEntry + "^" + subEntry + "^^" + std::string(tocPage) ) ; 
      }
      createTOC(pdf,tocList) ; 
      tocEntry = "Per event module running time" ;
      sprintf(tocPage,"%d",pageNumber) ; 
      tocList.push_back( tocEntry + "^^" + std::string(tocPage) ) ; 
      for (HLTPerformanceInfo::Modules::const_iterator modIter=HLTPerformance.beginModules();
           modIter!=HLTPerformance.endModules(); modIter++) {
          if (!useModule(*modIter,skipTiming)) continue ;  
          std::string subEntry = modIter->name() ;
          sprintf(tocPage,"%d",pageNumber++) ;
          tocList.push_back( tocEntry + "^" + subEntry + "^^" + std::string(tocPage) ) ; 
      }
      createTOC(pdf,tocList) ; 
      if (HLTPerformance.numberOfPaths() > 1) {
          tocEntry = "Per event module (in path) running time" ;
          sprintf(tocPage,"%d",pageNumber) ; 
          tocList.push_back( tocEntry + "^^" + std::string(tocPage) ) ; 
          for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
               pathIter!=HLTPerformance.endPaths(); pathIter++) {
              std::string subEntry = pathIter->name() ; 
              sprintf(tocPage,"%d",pageNumber) ;
              tocList.push_back( tocEntry + "^" + subEntry + "^^" + std::string(tocPage) ) ; 
              for (HLTPerformanceInfo::Path::const_iterator modIter=pathIter->begin();
                   modIter!=pathIter->end(); modIter++) {
                  if (!useModule(*modIter,skipTiming)) continue ;  
                  subEntry = pathIter->name() + "^" + modIter->name() ;
                  sprintf(tocPage,"%d",pageNumber++) ;
                  tocList.push_back( tocEntry + "^" + subEntry + "^^" + std::string(tocPage) ) ; 
              }
          }
          createTOC(pdf,tocList) ; 
          tocEntry = "Average path time" ; 
          sprintf(tocPage,"%d",pageNumber++) ;
          tocList.push_back( tocEntry + "^^" + std::string(tocPage) ) ; 
          createTOC(pdf,tocList) ; 

          tocEntry = "Average incremental path time" ; 
          sprintf(tocPage,"%d",pageNumber++) ;
          tocList.push_back( tocEntry + "^^" + std::string(tocPage) ) ; 
          createTOC(pdf,tocList) ;
      }

      tocEntry = "Average module time" ; 
      sprintf(tocPage,"%d",pageNumber++) ;
      tocList.push_back( tocEntry + "^^" + std::string(tocPage) ) ; 
      createTOC(pdf,tocList) ; 
      tocEntry = "Average module running time" ; 
      sprintf(tocPage,"%d",pageNumber++) ;
      tocList.push_back( tocEntry + "^^" + std::string(tocPage) ) ; 
      createTOC(pdf,tocList) ; 
      if (HLTPerformance.numberOfPaths() > 1) {
          tocEntry = "Average module (in path) time" ;
          sprintf(tocPage,"%d",pageNumber) ; 
          tocList.push_back( tocEntry + "^^" + std::string(tocPage) ) ; 
          for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
               pathIter!=HLTPerformance.endPaths(); pathIter++) {
              std::string subEntry = pathIter->name() ;
              sprintf(tocPage,"%d",pageNumber++) ;
              tocList.push_back( tocEntry + "^" + subEntry + "^^" + std::string(tocPage) ) ; 
          }
          createTOC(pdf,tocList) ; 
          tocEntry = "Average module (in path) running time" ;
          sprintf(tocPage,"%d",pageNumber) ; 
          tocList.push_back( tocEntry + "^^" + std::string(tocPage) ) ; 
          for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
               pathIter!=HLTPerformance.endPaths(); pathIter++) {
              std::string subEntry = pathIter->name() ;
              sprintf(tocPage,"%d",pageNumber++) ;
              tocList.push_back( tocEntry + "^" + subEntry + "^^" + std::string(tocPage) ) ; 
          }
          createTOC(pdf,tocList) ; 
      }
      tocEntry = "Path success rate" ;
      sprintf(tocPage,"%d",pageNumber++) ; 
      tocList.push_back( tocEntry + "^^" + std::string(tocPage) ) ; 
      createTOC(pdf,tocList) ; 
      if (HLTPerformance.numberOfPaths() > 1) {
          tocEntry = "Path vs. Path success rate" ;
          sprintf(tocPage,"%d",pageNumber++) ; 
          tocList.push_back( tocEntry + "^^" + std::string(tocPage) ) ; 
          createTOC(pdf,tocList) ; 
          tocEntry = "Fraction of single path success" ;
          sprintf(tocPage,"%d",pageNumber++) ; 
          tocList.push_back( tocEntry + "^^" + std::string(tocPage) ) ; 
          createTOC(pdf,tocList) ;
      }
      tocEntry = "Path rejection factor" ;
      sprintf(tocPage,"%d",pageNumber++) ; 
      tocList.push_back( tocEntry + "^^" + std::string(tocPage) ) ; 
      createTOC(pdf,tocList) ; 
      tocEntry = "Failing module (by path)" ;
      sprintf(tocPage,"%d",pageNumber) ; 
      tocList.push_back( tocEntry + "^^" + std::string(tocPage) ) ; 
      for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
           pathIter!=HLTPerformance.endPaths(); pathIter++) {
          std::string subEntry = pathIter->name() ;
          sprintf(tocPage,"%d",pageNumber++) ;
          tocList.push_back( tocEntry + "^" + subEntry + "^^" + std::string(tocPage) ) ; 
      }
      createTOC(pdf,tocList) ; 
      tocEntry = "Module rejection factor (by path)" ;
      sprintf(tocPage,"%d",pageNumber) ; 
      tocList.push_back( tocEntry + "^^" + std::string(tocPage) ) ; 
      for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
           pathIter!=HLTPerformance.endPaths(); pathIter++) {
          std::string subEntry = pathIter->name() ;
          sprintf(tocPage,"%d",pageNumber++) ;
          tocList.push_back( tocEntry + "^" + subEntry + "^^" + std::string(tocPage) ) ; 
      }
      createTOC(pdf,tocList) ; 
      tocEntry = "Module rejection factor per unit running time (by path)" ;
      sprintf(tocPage,"%d",pageNumber) ; 
      tocList.push_back( tocEntry + "^^" + std::string(tocPage) ) ; 
      for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
           pathIter!=HLTPerformance.endPaths(); pathIter++) {
          std::string subEntry = pathIter->name() ;
          sprintf(tocPage,"%d",pageNumber++) ;
          tocList.push_back( tocEntry + "^" + subEntry + "^^" + std::string(tocPage) ) ; 
      }
      createTOC(pdf,tocList) ; 

      //--- Bookmarks output ---//
      for (unsigned int i=0; i<tocList.size(); i++) {
          std::string entry = tocList.at(i).append("\n") ;
          txtfile << entry ;
      }

      //-----------------------//
      //--- Plot Histograms ---//
      //-----------------------//

      //--- Event timing ---//
      plot1D(totalTime,c1) ;
      plotMany(pathTime,c1,Path,HLTPerformance,skipTiming) ;
      if (HLTPerformance.numberOfPaths() > 1)
          plotMany(incPathTime,c1,Path,HLTPerformance,skipTiming) ;
      plotMany(moduleTime,c1,Module,HLTPerformance,skipTiming) ;
      plotMany(moduleScaledTime,c1,Module,HLTPerformance,skipTiming) ; 
      if (HLTPerformance.numberOfPaths() > 1) {
          plotModuleInPath(moduleInPathScaledTime,c1,HLTPerformance,skipTiming) ; 
      }
      //--- Average time (summary) plots ---//
      if (HLTPerformance.numberOfPaths() > 1) {
          plot1D(pathTimeSummary,c1) ; 
          plot1D(incPathTimeSummary,c1) ; 
      }
      plot1D(moduleTimeSummary,c1) ; 
      plot1D(moduleScaledTimeSummary,c1) ; 
      if (HLTPerformance.numberOfPaths() > 1) {
          plotMany(moduleInPathTimeSummary,c1,Path,HLTPerformance,skipTiming) ;
          plotMany(moduleInPathScaledTimeSummary,c1,Path,HLTPerformance,skipTiming) ;
      }
      //--- Success/Rejection plots ---//
      plot1D(pathSuccessFraction,c1) ; 
      if (HLTPerformance.numberOfPaths() > 1) {
          pathVsPathSummary->SetMaximum(1.) ;
          pathVsPathSummary->SetMinimum(0.) ;
          pathVsPathSummary->SetStats(kFALSE) ; 
          pathVsPathSummary->Draw("colz") ; c1->Update() ;
          plot1D(uniquePathSuccessFraction,c1) ; 
      }
      plot1D(pathRejection,pathRejectAll,c1) ; 
      plotMany(failedModule,c1,Path,HLTPerformance,skipTiming) ;
      plotMany(moduleInPathRejection,moduleInPathRejectAll,c1,Path,HLTPerformance,skipTiming) ;
      plotMany(moduleInPathRejectTime,c1,Path,HLTPerformance,skipTiming) ;

      pdf->Close() ; 
  }
  
  txtfile.close() ; 
  outFile->Write() ;
  outFile->Close() ; 
  file.Close() ; 

  return 0;
}

