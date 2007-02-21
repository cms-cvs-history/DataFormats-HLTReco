#include <TROOT.h>
#include <TFile.h>
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

#include "DataFormats/Common/interface/EventAux.h"
#include "DataFormats/HLTReco/interface/HLTPerformanceInfo.h"
#include "FWCore/FWLite/src/AutoLibraryLoader.h"

//--- Bryan Dahmes, January 2007 ---//

int main(int argc, char ** argv)
{
  //-- Load libraries ---//
  gSystem->Load("libFWCoreFWLite") ;
  AutoLibraryLoader::enable() ;

  //--- Default arguments ---//
  std::string filename = "hlt.root" ;
  std::string outbase  = "hltTimingSummary" ; 
  std::string outname  = outbase + ".root" ;
  std::string pdfname  = outbase + ".pdf" ;
  std::string txtname  = outbase + ".txt" ;

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
      ("time,t",      boost::program_options::value<double>(),
       "All relevant histogram time axes run from 0 to the user-specified value (in msec)")
      ("noFirst,f","Skip the first event (Default is to run over all events)") ; 


  std::string usage = "\nSample hltTimingSummary usage::\n" ; 
  usage += "\"hltTimingSummary -t 50\" " ; 
  usage += "inputs hlt.root, outputs hltTimingsummary.root, .txt, and .pdf\n" ; 
  usage += "                         and timing histograms run from 0 to 50 msec\n" ; 
  usage += "\"hltTimingSummary -i input.root -o output\" " ;
  usage += "inputs input.root, outputs output.root, .txt, and .pdf\n" ; 
  usage += "\"hltTimingSummary -o output -b special.txt\" " ;
  usage += "inputs hlt.root, outputs hltTimingsummary.root and .pdf, special.txt\n\n" ; 

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

  std::cout << "Opening file " << filename << std::endl ;
  TFile file(filename.c_str());
  if (file.IsZombie()) {
      std::cout << "*** Error opening file: " << filename << " ***" << std::endl;
      std::cout << desc << usage <<  std::endl ;
      return 1 ;
  }

  TTree * events = dynamic_cast<TTree *>(file.Get("Events") );
  assert(events);

  TBranch * TBPerfInfo = events->GetBranch("HLTPerformanceInfo_pts__PROD.obj");
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
  std::vector<double> moduleInitTimer ; 
  std::vector<int> moduleInitFailures ; 
  std::vector< std::vector<double> > moduleInPathTimer ; 
  std::vector< std::vector<int> > failures ; 
  std::vector<std::string> tocList ; 
  std::vector<bool> successPerEvent ;
  std::vector<int> uniquePathSuccess ; 
  std::vector< std::vector<int> > successPathVsPath ; 
  
  //--- Histogram initialization ---//
  std::vector<TH1D*> pathTime ; 
  std::vector<TH1D*> incPathTime ; 
  std::vector<TH1D*> moduleInPathTime ; 
  std::vector<TH1D*> moduleInPathScaledTime ; 
  std::vector<TH1D*> failedModule ; 
  std::vector<TH1D*> moduleInPathRejection ; 
  std::vector<TH1D*> moduleInPathRejectTime ; 
  
  //--- Initial loop to get scale for histograms ---//
  double maxTime = 0. ;
  double xmin = 0. ;
  double xmax = 0. ;
  int firstEvent = 0 ;
  if (skipFirstEvent) {
      firstEvent = 1 ;
      std::cout << "Skipping the first event" << std::endl ; 
  }
  for (int ievt=firstEvent; ievt<n_evts; ievt++) {
      TBPerfInfo->GetEntry(ievt) ; 
      if ( HLTPerformance.totalTime() > maxTime ) maxTime = HLTPerformance.totalTime() ; 
  }

  if (userMaxTime > 0) {
      xmax = userMaxTime ;
  } else {
      int xscale = 2 ; 
      while (pow(10,xscale--) > maxTime) ; 
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
      new TH1D("pathSuccessFraction","Success rate (%) for each path",
               HLTPerformance.numberOfPaths(),0.,double(HLTPerformance.numberOfPaths())) ; 
  TH1D* uniquePathSuccessFraction =
      new TH1D("uniquePathSuccessFraction","Fraction (%) of events passing due to a single path",
               HLTPerformance.numberOfPaths(),0.,double(HLTPerformance.numberOfPaths())) ; 
  TH1D* pathRejection =
      new TH1D("pathRejection","Rejection for each path",
               HLTPerformance.numberOfPaths(),0.,double(HLTPerformance.numberOfPaths())) ; 

  TH2D* pathVsPathSummary = new TH2D("pathVsPathSummary", "Relative path success",
                               HLTPerformance.numberOfPaths(),
                               0.,double(HLTPerformance.numberOfPaths()), 
                               HLTPerformance.numberOfPaths(),
                               0.,double(HLTPerformance.numberOfPaths())) ; 

  totalTime->GetXaxis()->SetTitle("msec") ; 
  pathTimeSummary->GetYaxis()->SetTitle("msec") ;
  incPathTimeSummary->GetYaxis()->SetTitle("msec") ;

  char name[50] ; 
  char title[100] ; 
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
              
      pathTimeSummary->GetXaxis()->SetBinLabel(pCtr+1,pathIter->name().data()) ;
      incPathTimeSummary->GetXaxis()->SetBinLabel(pCtr+1,pathIter->name().data()) ; 
      pathSuccessFraction->GetXaxis()->SetBinLabel(pCtr+1,pathIter->name().data()) ; 
      uniquePathSuccessFraction->GetXaxis()->SetBinLabel(pCtr+1,pathIter->name().data()) ; 
      pathRejection->GetXaxis()->SetBinLabel(pCtr+1,pathIter->name().data()) ;

      pathVsPathSummary->GetXaxis()->SetBinLabel(pCtr+1,pathIter->name().data()) ; 
      pathVsPathSummary->GetYaxis()->SetBinLabel(pCtr+1,pathIter->name().data()) ; 
              
      sprintf(name,"pathTime_%s",pathIter->name().c_str()) ;
      sprintf(title,"Total time for path %s",pathIter->name().c_str()) ;
      histo = new TH1D(name,title,100,xmin,xmax) ;
      histo->GetXaxis()->SetTitle("msec") ;
      pathTime.push_back( histo ) ; 

      sprintf(name,"incPathTime_%s",pathIter->name().c_str()) ;
      sprintf(title,"Total incremental time for path %s",pathIter->name().c_str()) ;
      histo = new TH1D(name,title,100,xmin,xmax) ; 
      histo->GetXaxis()->SetTitle("msec") ; 
      incPathTime.push_back( histo ) ; 

      sprintf(name,"failedModule_%s",pathIter->name().c_str()) ;
      sprintf(title,"Failure fraction (%%) by module for path %s",pathIter->name().c_str()) ;
      histo = new TH1D(name,title,
                       (1+pathIter->numberOfModules()),-1.,
                       double(pathIter->numberOfModules())) ; 
      histo->GetXaxis()->SetBinLabel(1,"SUCCESS") ; 
      failedModule.push_back( histo ) ; 

      sprintf(name,"moduleInPathTime_%s",pathIter->name().c_str()) ;
      sprintf(title,"Average time per module for path %s",pathIter->name().c_str()) ;
      histo = new TH1D(name,title,
                       (pathIter->numberOfModules()),0.,
                       double(pathIter->numberOfModules())) ; 
      histo->GetYaxis()->SetTitle("msec") ; 
      moduleInPathTime.push_back( histo ) ; 

      sprintf(name,"moduleInPathScaledTime_%s",pathIter->name().c_str()) ;
      sprintf(title,"Average module RUNNING time for path %s",
              pathIter->name().c_str()) ;
      histo = new TH1D(name,title,
                       (pathIter->numberOfModules()),0.,
                       double(pathIter->numberOfModules())) ; 
      histo->GetYaxis()->SetTitle("msec") ; 
      moduleInPathScaledTime.push_back( histo ) ; 
      
      sprintf(name,"moduleInPathRejection_%s",pathIter->name().c_str()) ;
      sprintf(title,"Rejection per module for path %s",pathIter->name().c_str()) ;
      histo = new TH1D(name,title,
                       (pathIter->numberOfModules()),0.,
                       double(pathIter->numberOfModules())) ; 
      moduleInPathRejection.push_back( histo ) ; 

      sprintf(name,"moduleInPathRejectTime_%s",pathIter->name().c_str()) ;
      sprintf(title,"(Rejection / average running time) per module for path %s",pathIter->name().c_str()) ;
      histo = new TH1D(name,title,
                       (pathIter->numberOfModules()),0.,
                       double(pathIter->numberOfModules())) ; 
      histo->GetYaxis()->SetTitle("1/msec") ; 
      moduleInPathRejectTime.push_back( histo ) ; 

      mCtr = 1 ;
      for (HLTPerformanceInfo::Path::const_iterator modIter=pathIter->begin();
           modIter!=pathIter->end(); modIter++) {
          moduleInitTimer.push_back(0.) ; 
          moduleInitFailures.push_back(0) ; 
          const char* modName = modIter->name().c_str() ;
          moduleInPathTime[pCtr]->GetXaxis()->SetBinLabel(mCtr,modName) ;
          moduleInPathScaledTime[pCtr]->GetXaxis()->SetBinLabel(mCtr,modName) ;
          moduleInPathRejection[pCtr]->GetXaxis()->SetBinLabel(mCtr,modName) ;
          moduleInPathRejectTime[pCtr]->GetXaxis()->SetBinLabel(mCtr,modName) ;
          failedModule[pCtr]->GetXaxis()->SetBinLabel(++mCtr,modName) ;
      }

      pCtr++ ;
      moduleInPathTimer.push_back( moduleInitTimer ) ;
      failures.push_back( moduleInitFailures ) ;
      moduleInitTimer.clear() ; 
      moduleInitFailures.clear() ; 
  }
  
  for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
       pathIter!=HLTPerformance.endPaths(); pathIter++) {
      successPathVsPath.push_back( pathSuccess ) ; 
  }
  
  TH1D* moduleTimeSummary =
      new TH1D("moduleTimeSummary","Average time per module",
               HLTPerformance.numberOfModules(),0.,double(HLTPerformance.numberOfModules())) ;
  moduleTimeSummary->GetYaxis()->SetTitle("msec") ; 
  
  mCtr = 1 ; 
  for (HLTPerformanceInfo::Modules::const_iterator modIter=HLTPerformance.beginModules();
       modIter!=HLTPerformance.endModules(); modIter++) {
      moduleTimer.push_back(0.) ;
      moduleTimeSummary->GetXaxis()->SetBinLabel(mCtr++,modIter->name().c_str()) ; 
  }

  for (int ievt=firstEvent; ievt<n_evts; ievt++) {
      TBPerfInfo->GetEntry(ievt) ;

      pCtr = 0 ;
      for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
           pathIter!=HLTPerformance.endPaths(); pathIter++) {
          pathTimer.at(pCtr) += pathIter->time() ; 

          double addedPathTime = 0 ;
          mCtr = 0 ; 
          for (HLTPerformanceInfo::Path::const_iterator modIter=pathIter->begin();
               modIter!=pathIter->end(); modIter++) {
              moduleInPathTimer.at(pCtr).at(mCtr++) += modIter->time() ; 
              if (HLTPerformance.uniqueModule(modIter->name().c_str())) {
                  addedPathTime += modIter->time() ;
              }
          }

          incPathTimer.at(pCtr) += addedPathTime ;
          pathTime.at(pCtr)->Fill( 1000. * pathIter->time() ) ; 
          incPathTime.at(pCtr)->Fill( 1000. * addedPathTime ) ; 
          if (pathIter->status().accept()) {
              pathSuccess.at(pCtr)++ ;
              successPerEvent.at(pCtr) = true ; 
          } else { //--- One of the modules caused the path to fail ---//
              bool foundFailure = false ; 
              successPerEvent.at(pCtr) = false ; 
              int failMod = 0 ; 
              for (HLTPerformanceInfo::Path::const_iterator modIter=pathIter->begin();
                   modIter!=pathIter->end(); modIter++) {
                  if ((modIter->time() != 0)&&(!foundFailure)) {
                      failMod++ ;
                  } else {
                      foundFailure = true ; 
                  }
              }
              failures.at(pCtr).at(failMod-1)++ ; 
          }
          pCtr++ ; 
      }


      for (unsigned int i=0; i<pathTimer.size(); i++) {
          bool uniqueSuccess = successPerEvent.at(i) ; 
          for (unsigned int j=0; j<pathTimer.size(); j++) {
              if ( successPerEvent.at(i) ) {
                  if ( successPerEvent.at(j) ) {
                      if (i != j) uniqueSuccess = false ;
                      successPathVsPath.at(i).at(j)++ ; 
                  }
              }
          }
          if ( successPerEvent.at(i) && uniqueSuccess ) uniquePathSuccess.at(i)++ ; 
      }
      
      mCtr = 0 ;
      for (HLTPerformanceInfo::Modules::const_iterator modIter=HLTPerformance.beginModules();
           modIter!=HLTPerformance.endModules(); modIter++) {
          moduleTimer.at(mCtr++) += modIter->time() ;
      }

      totalTime->Fill( 1000. * HLTPerformance.totalTime() ) ; 
  }

  for (unsigned int pCtr=0; pCtr<moduleInPathTimer.size(); pCtr++) {
      int totalFailures = 0 ;
      int evtsLeft = n_evts - firstEvent ; 
      for (unsigned int mCtr=0; mCtr<moduleInPathTimer.at(pCtr).size(); mCtr++) { 

          double avgModuleTime = 1000. * moduleInPathTimer.at(pCtr).at(mCtr)/double(n_evts - firstEvent) ; // msec 
          moduleInPathTime.at(pCtr)->Fill( double(mCtr), avgModuleTime ) ; 
          

          //--- Calculate the rejection factor for each module in path, ---//
          //--- defined as N(HLT input)/N(HLT accept) ---//
          double modReject = 0 ;
          if (failures.at(pCtr).at(mCtr) < evtsLeft) {
              modReject = double(evtsLeft) / double(evtsLeft - failures.at(pCtr).at(mCtr)) ; 
          } else {
              modReject = 1000. ;
          }
          moduleInPathRejection.at(pCtr)->Fill( double(mCtr), modReject ) ; 

          //--- Calculate rejection per unit time for each module in path ---//
          double scaledAvgModuleTime = avgModuleTime * (n_evts - firstEvent) / evtsLeft ; 
          moduleInPathScaledTime.at(pCtr)->Fill( double(mCtr), scaledAvgModuleTime ) ; 
          moduleInPathRejectTime.at(pCtr)->Fill( double(mCtr),
                                                 modReject / scaledAvgModuleTime ) ; 
          
          evtsLeft -= failures.at(pCtr).at(mCtr) ; 
          totalFailures += failures.at(pCtr).at(mCtr) ; 

          if (failures.at(pCtr).at(mCtr) >= 0) {
              failedModule.at(pCtr)->Fill( double(mCtr), 100. * double(failures.at(pCtr).at(mCtr))/double(n_evts - firstEvent) );
          }
      }
      failedModule.at(pCtr)->Fill(-1.,100. * (double(n_evts-firstEvent-totalFailures)/double(n_evts-firstEvent))) ; 

      double pathReject = 1000. ;
      if ((n_evts-firstEvent) > totalFailures)
          pathReject = double(n_evts-firstEvent) / double(n_evts - firstEvent - totalFailures) ; 
      pathRejection->Fill(double(pCtr), pathReject) ; 
  }

          

  
  for (unsigned int i=0; i<pathTimer.size(); i++) {
      pathTimeSummary->Fill(double(i), pathTimer.at(i)/(n_evts-firstEvent)) ;
      incPathTimeSummary->Fill(double(i), incPathTimer.at(i)/(n_evts-firstEvent)) ;
      pathSuccessFraction->Fill(double(i), 100. * double(pathSuccess.at(i))/(n_evts-firstEvent)) ;

      uniquePathSuccessFraction->Fill(double(i),
                                      100. * double(uniquePathSuccess.at(i))/(n_evts-firstEvent)) ; 
      
      for (unsigned int j=0; j<pathTimer.size(); j++) {
          pathVsPathSummary->Fill( double(i), double(j),
                                   double(successPathVsPath.at(i).at(j))/(n_evts-firstEvent) ) ; 
      }
  }

  for (unsigned int i=0; i<moduleTimer.size(); i++) {
      moduleTimeSummary->Fill(double(i), 1000. * (moduleTimer.at(i)/(n_evts-firstEvent))) ;
  }

  //--- Dump results ---//
  {
      TCanvas* c1 = new TCanvas("c1") ;
      pathTimeSummary->SetMinimum(0.) ; 
      incPathTimeSummary->SetMinimum(0.) ; 
      pathSuccessFraction->SetMinimum(0.) ; 
      moduleTimeSummary->SetMinimum(0.) ; 
      TPDF* pdf = new TPDF(pdfname.c_str()) ; 
      int pageNumber = 1 ; 
      char number[3] ; 

      double majorMajorSpacing = 0.045 ;
      double majorMinorSpacing = 0.025 ;
      double minorMajorSpacing = 0.035 ;
      double titleSize = 0.050 ; 
      double majorSize = 0.032 ;
      double minorSize = 0.027 ; 

      
      gROOT->SetStyle("Plain") ; 
      gStyle->SetPalette(1) ; 
      c1->UseCurrentStyle() ; 
      gROOT->ForceStyle() ;
      
      pdf->SetTextSize(titleSize) ; 
      pdf->Text(0.33,0.93,"Timing Summary Output") ; 

      double majorXval = 0.10 ;
      double minorXval = 0.15 ;
      double pageXval = 0.90 ; 
      double mmpYval = 0.87 ; 

      //--- Total time.  No subplots ---//
      pdf->SetTextSize(majorSize) ; 
      pdf->Text(majorXval,mmpYval,"Total time") ;
      pageNumber++ ; 
      sprintf(number,"%3i",pageNumber) ; 
      pdf->Text(pageXval,mmpYval,number) ;
      pdf->SetTextSize(minorSize) ; mmpYval -= majorMinorSpacing ; 
      pdf->Text(minorXval,mmpYval,"Cumulative time for all modules in the event") ; 
      std::string tocEntry = "Total time" ; 
      char tocPage[5] ; sprintf(tocPage,"%d",pageNumber) ;
      tocList.push_back( tocEntry + "^^" + std::string(tocPage) ) ; 

      //--- Total time by path.  Plots for each path ---//
      pdf->SetTextSize(majorSize) ; mmpYval -= minorMajorSpacing ; 
      pdf->Text(majorXval,mmpYval,"Total path time") ; 
      pageNumber++ ;
      sprintf(number,"%3i",pageNumber) ; 
      pdf->Text(pageXval,mmpYval,number) ;
      pdf->SetTextSize(minorSize) ; mmpYval -= majorMinorSpacing ; 
      pdf->Text(minorXval,mmpYval,"Cumulative time for individual paths in the event") ; 
      
      tocEntry = "Total path time" ;
      sprintf(tocPage,"%d",pageNumber) ; 
      tocList.push_back( tocEntry + "^^" + std::string(tocPage) ) ; 
      for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
           pathIter!=HLTPerformance.endPaths(); pathIter++) {
          std::string subEntry = pathIter->name() ;
          sprintf(tocPage,"%d",pageNumber++) ;
          tocList.push_back( tocEntry + "^" + subEntry + "^^" + std::string(tocPage) ) ; 
      }

      //--- Average path time.  No subplots ---//
      pdf->SetTextSize(majorSize) ; mmpYval -= minorMajorSpacing ; 
      pdf->Text(majorXval,mmpYval,"Average path time") ; 
      sprintf(number,"%3i",pageNumber) ; 
      pdf->Text(pageXval,mmpYval,number) ; 

      tocEntry = "Average path time" ; 
      sprintf(tocPage,"%d",pageNumber++) ;
      tocList.push_back( tocEntry + "^^" + std::string(tocPage) ) ; 
      
      //--- Incremental path time.  Plots for each path ---//
      pdf->SetTextSize(majorSize) ; mmpYval -= majorMajorSpacing ; 
      pdf->Text(majorXval,mmpYval,"Incremental path time") ; 
      sprintf(number,"%3i",pageNumber) ; 
      pdf->Text(pageXval,mmpYval,number) ;
      pdf->SetTextSize(minorSize) ; mmpYval -= majorMinorSpacing ; 
      pdf->Text(minorXval,mmpYval,"Incremental time due to unique modules in each path") ; 

      tocEntry = "Incremental path time" ;
      sprintf(tocPage,"%d",pageNumber) ; 
      tocList.push_back( tocEntry + "^^" + std::string(tocPage) ) ; 
      for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
           pathIter!=HLTPerformance.endPaths(); pathIter++) {
          std::string subEntry = pathIter->name() ;
          sprintf(tocPage,"%d",pageNumber++) ;
          tocList.push_back( tocEntry + "^" + subEntry + "^^" + std::string(tocPage) ) ; 
      }

      //--- Average incremental path time.  No subplots ---//
      pdf->SetTextSize(majorSize) ; mmpYval -= minorMajorSpacing ; 
      pdf->Text(majorXval,mmpYval,"Average incremental path time") ; 
      sprintf(number,"%3i",pageNumber) ; 
      pdf->Text(pageXval,mmpYval,number) ; 

      tocEntry = "Average incremental path time" ; 
      sprintf(tocPage,"%d",pageNumber++) ;
      tocList.push_back( tocEntry + "^^" + std::string(tocPage) ) ; 

      //--- Average module time.  Plots for each path ---//
      pdf->SetTextSize(majorSize) ; mmpYval -= majorMajorSpacing ; 
      pdf->Text(majorXval,mmpYval,"Average module time (by path)") ; 
      sprintf(number,"%3i",pageNumber) ; 
      pdf->Text(pageXval,mmpYval,number) ;

      tocEntry = "Average module time (by path)" ;
      sprintf(tocPage,"%d",pageNumber) ; 
      tocList.push_back( tocEntry + "^^" + std::string(tocPage) ) ; 
      for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
           pathIter!=HLTPerformance.endPaths(); pathIter++) {
          std::string subEntry = pathIter->name() ;
          sprintf(tocPage,"%d",pageNumber++) ;
          tocList.push_back( tocEntry + "^" + subEntry + "^^" + std::string(tocPage) ) ; 
      }
      
      //--- Average module time.  No subplots ---//
      pdf->SetTextSize(majorSize) ; mmpYval -= majorMajorSpacing ; 
      pdf->Text(majorXval,mmpYval,"Average module time") ; 
      sprintf(number,"%3i",pageNumber) ; 
      pdf->Text(pageXval,mmpYval,number) ; 

      tocEntry = "Average module time" ; 
      sprintf(tocPage,"%d",pageNumber++) ;
      tocList.push_back( tocEntry + "^^" + std::string(tocPage) ) ; 

      //--- Average module running time.  Plot for each path ---//
      pdf->SetTextSize(majorSize) ; mmpYval -= majorMajorSpacing ; 
      pdf->Text(majorXval,mmpYval,"Average module RUNNING time (by path)") ; 
      sprintf(number,"%3i",pageNumber) ; 
      pdf->Text(pageXval,mmpYval,number) ;
      pdf->SetTextSize(minorSize) ; mmpYval -= majorMinorSpacing ; 
      pdf->Text(minorXval,mmpYval,"RUNNING time omits instances when module time is zero") ; 

      tocEntry = "Average module running time (by path)" ;
      sprintf(tocPage,"%d",pageNumber) ; 
      tocList.push_back( tocEntry + "^^" + std::string(tocPage) ) ; 
      for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
           pathIter!=HLTPerformance.endPaths(); pathIter++) {
          std::string subEntry = pathIter->name() ;
          sprintf(tocPage,"%d",pageNumber++) ;
          tocList.push_back( tocEntry + "^" + subEntry + "^^" + std::string(tocPage) ) ; 
      }

      //--- Success rate.  No subplots ---//
      pdf->SetTextSize(majorSize) ; mmpYval -= minorMajorSpacing ; 
      pdf->Text(majorXval,mmpYval,"Path success rate") ; 
      sprintf(number,"%3i",pageNumber) ; 
      pdf->Text(pageXval,mmpYval,number) ; 

      tocEntry = "Path success rate" ;
      sprintf(tocPage,"%d",pageNumber++) ; 
      tocList.push_back( tocEntry + "^^" + std::string(tocPage) ) ; 

      //--- Relative success rate.  No subplots ---//
      pdf->SetTextSize(majorSize) ; mmpYval -= majorMajorSpacing ; 
      pdf->Text(majorXval,mmpYval,"Path vs. Path success rate") ; 
      sprintf(number,"%3i",pageNumber) ; 
      pdf->Text(pageXval,mmpYval,number) ; 
      pdf->SetTextSize(minorSize) ; mmpYval -= majorMinorSpacing ; 
      pdf->Text(minorXval,mmpYval,"Path success rates shown relative to other paths") ; 

      tocEntry = "Path vs. Path success rate" ;
      sprintf(tocPage,"%d",pageNumber++) ; 
      tocList.push_back( tocEntry + "^^" + std::string(tocPage) ) ; 

      //--- Unique success rate.  No subplots ---//
      pdf->SetTextSize(majorSize) ; mmpYval -= minorMajorSpacing ; 
      pdf->Text(majorXval,mmpYval,"Fraction of event success due to a single path") ; 
      sprintf(number,"%3i",pageNumber) ; 
      pdf->Text(pageXval,mmpYval,number) ; 

      tocEntry = "Fraction of single path success" ;
      sprintf(tocPage,"%d",pageNumber++) ; 
      tocList.push_back( tocEntry + "^^" + std::string(tocPage) ) ; 

      //--- Failing module.  One plot per path ---//
      pdf->SetTextSize(majorSize) ; mmpYval -= majorMajorSpacing ; 
      pdf->Text(majorXval,mmpYval,"Failing module (by path)") ; 
      sprintf(number,"%3i",pageNumber) ; 
      pdf->Text(pageXval,mmpYval,number) ;
      pdf->SetTextSize(minorSize) ; mmpYval -= majorMinorSpacing ; 
      pdf->Text(minorXval,mmpYval,"Records failure rate as a function of modules in the path") ; 

      tocEntry = "Failing module (by path)" ;
      sprintf(tocPage,"%d",pageNumber) ; 
      tocList.push_back( tocEntry + "^^" + std::string(tocPage) ) ; 
      for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
           pathIter!=HLTPerformance.endPaths(); pathIter++) {
          std::string subEntry = pathIter->name() ;
          sprintf(tocPage,"%d",pageNumber++) ;
          tocList.push_back( tocEntry + "^" + subEntry + "^^" + std::string(tocPage) ) ; 
      }

      //--- Module rejection factor.  One plot per path ---//
      pdf->SetTextSize(majorSize) ; mmpYval -= minorMajorSpacing ; 
      pdf->Text(majorXval,mmpYval,"Module rejection factor (by path)") ; 
      sprintf(number,"%3i",pageNumber) ; 
      pdf->Text(pageXval,mmpYval,number) ; 
      pdf->SetTextSize(minorSize) ; mmpYval -= majorMinorSpacing ; 
      pdf->Text(minorXval,mmpYval,"Rejection factor = (Number of input events) / (Number of events that pass)") ; 
      
      tocEntry = "Module rejection factor (by path)" ;
      sprintf(tocPage,"%d",pageNumber) ; 
      tocList.push_back( tocEntry + "^^" + std::string(tocPage) ) ; 
      for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
           pathIter!=HLTPerformance.endPaths(); pathIter++) {
          std::string subEntry = pathIter->name() ;
          sprintf(tocPage,"%d",pageNumber++) ;
          tocList.push_back( tocEntry + "^" + subEntry + "^^" + std::string(tocPage) ) ; 
      }

      //--- Path rejection factor.  No subplots ---//
      pdf->SetTextSize(majorSize) ; mmpYval -= minorMajorSpacing ;
      pdf->Text(majorXval,mmpYval,"Path rejection factor") ; 
      sprintf(number,"%3i",pageNumber) ; 
      pdf->Text(pageXval,mmpYval,number) ; 

      tocEntry = "Path rejection factor" ;
      sprintf(tocPage,"%d",pageNumber++) ; 
      tocList.push_back( tocEntry + "^^" + std::string(tocPage) ) ; 

      //--- Module rejection per running time.  One plot per path ---//
      pdf->SetTextSize(majorSize) ; mmpYval -= majorMajorSpacing ;
      pdf->Text(majorXval,mmpYval,"Module rejection factor per unit running time (by path)") ; 
      sprintf(number,"%3i",pageNumber) ; 
      pdf->Text(pageXval,mmpYval,number) ; 
      pdf->SetTextSize(minorSize) ; mmpYval -= majorMinorSpacing ;
      pdf->Text(minorXval,mmpYval,"Use the module RUNNING time to compute Rejection/Time") ; 

      tocEntry = "Module rejection factor per unit running time (by path)" ;
      sprintf(tocPage,"%d",pageNumber) ; 
      tocList.push_back( tocEntry + "^^" + std::string(tocPage) ) ; 
      for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
           pathIter!=HLTPerformance.endPaths(); pathIter++) {
          std::string subEntry = pathIter->name() ;
          sprintf(tocPage,"%d",pageNumber++) ;
          tocList.push_back( tocEntry + "^" + subEntry + "^^" + std::string(tocPage) ) ; 
      }

      //--- Bookmarks output ---//
      for (unsigned int i=0; i<tocList.size(); i++) {
          std::string entry = tocList.at(i).append("\n") ;
          txtfile << entry ;
      }
      
      //--- Plot Histograms ---//
      totalTime->Draw() ; c1->Update() ;
      int ctr = 0 ; 
      for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
           pathIter!=HLTPerformance.endPaths(); pathIter++) {
          pathTime[ctr++]->Draw() ; c1->Update() ;
      }
      pathTimeSummary->Draw() ; c1->Update() ; 
      ctr = 0 ; 
      for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
           pathIter!=HLTPerformance.endPaths(); pathIter++) {
          incPathTime[ctr++]->Draw() ; c1->Update() ;
      }
      incPathTimeSummary->Draw() ; c1->Update() ; 
      
      ctr = 0 ; 
      for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
           pathIter!=HLTPerformance.endPaths(); pathIter++) {
          moduleInPathTime[ctr]->SetMinimum(0.) ; 
          moduleInPathTime[ctr++]->Draw() ; c1->Update() ;
      }
      moduleTimeSummary->Draw() ; c1->Update() ; 
      
      ctr = 0 ; 
      for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
           pathIter!=HLTPerformance.endPaths(); pathIter++) {
          moduleInPathScaledTime[ctr]->SetMinimum(0.) ; 
          moduleInPathScaledTime[ctr++]->Draw() ; c1->Update() ;
      }

      pathSuccessFraction->Draw() ; c1->Update() ;

      pathVsPathSummary->SetMaximum(1.) ;
      pathVsPathSummary->SetMinimum(0.) ;
      pathVsPathSummary->SetStats(kFALSE) ; 
      pathVsPathSummary->Draw("colz") ; c1->Update() ;

      uniquePathSuccessFraction->Draw() ; c1->Update() ;
      
      ctr = 0 ; 
      for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
           pathIter!=HLTPerformance.endPaths(); pathIter++) {
          failedModule[ctr]->SetMinimum(0.) ; 
          failedModule[ctr++]->Draw() ; c1->Update() ;
      }
      
      ctr = 0 ; 
      for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
           pathIter!=HLTPerformance.endPaths(); pathIter++) {
          moduleInPathRejection[ctr]->SetMinimum(0.) ;
          moduleInPathRejection[ctr++]->Draw() ; c1->Update() ;
      }
      pathRejection->SetMinimum(0.) ; pathRejection->Draw() ; c1->Update() ; 
      ctr = 0 ; 
      for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
           pathIter!=HLTPerformance.endPaths(); pathIter++) {
          moduleInPathRejectTime[ctr]->SetMinimum(0.) ;
          moduleInPathRejectTime[ctr++]->Draw() ; c1->Update() ;
      }
      
      pdf->Close() ; 
  }

  txtfile.close() ; 
  outFile->Write() ;
  outFile->Close() ; 
  file.Close() ; 

  return 0;
}
