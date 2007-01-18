#include <TROOT.h>
#include <TFile.h>
#include <TH1D.h>
#include <TAttText.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPDF.h>
#include <TTree.h>
#include <TBranch.h>
#include <TSystem.h>

#include <iostream>
#include <fstream>
#include <string>

#include "DataFormats/Common/interface/EventAux.h"
#include "DataFormats/HLTReco/interface/HLTPerformanceInfo.h"
#include "FWCore/FWLite/src/AutoLibraryLoader.h"

using std::cout; using std::endl; using std::string;

//--- Bryan Dahmes, January 2007 ---//

//--------------------------------------------------------------------------------------//
//--- usage: hltTimingSummary <in.root> <out.root> <out.pdf> <out.txt>               ---//
//---        <in.root> = File with relevant HLT information                          ---//
//---                    (default name: hlt.root)                                    ---//
//---       <out.root> = Output root file                                            ---//
//---                    (default name: hltTimingSummary.root)                       ---//
//---        <out.pdf> = Histogram output in pdf format                              ---//
//---                    (default name: hltTimingSummary.pdf)                        ---//
//---        <out.txt> = With bookmarkPdf.pl, used to add PDF bookmarks to <out.pdf> ---//
//---                    (default name: hltTimingSummary-bookmarks.txt)              ---//
//--------------------------------------------------------------------------------------//


int main(int argc, char ** argv)
{
  //-- Load libraries ---//
  gSystem->Load("libFWCoreFWLite") ;
  AutoLibraryLoader::enable() ;

  //--- Default arguments ---//
  string filename = "hlt.root" ;
  string outname  = "hltTimingSummary.root" ; 
  string pdfname  = "hltTimingSummary.pdf" ; 
  string txtname ; 

  //--- Get parameters from command line ---//
  if (argc >= 2) filename = argv[1];
  if (argc >= 3) outname = argv[2];
  if (argc >= 4) pdfname = argv[3];
  if (argc >= 5) {
      txtname = argv[4];
  } else {
      txtname = pdfname ; 
      string::size_type loc = txtname.find(".pdf",0) ;
      string txtstring = ".txt" ;
      txtname.replace(loc,txtstring.length(),txtstring) ;
  }

  cout << "Opening file " << filename << endl ;
  TFile file(filename.c_str());
  if (file.IsZombie()) {
      cout << "*** Error opening file: " << filename << " ***" << endl;
      exit(-1);
  }
  TTree * events = dynamic_cast<TTree *>(file.Get("Events") );
  assert(events);

  TBranch * TBPerfInfo = events->GetBranch("HLTPerformanceInfo_pts__PROD.obj");
  assert(TBPerfInfo);
  HLTPerformanceInfo HLTPerformance ;
  TBPerfInfo->SetAddress((void *) & HLTPerformance) ;

  //--- Prepare the output ---//
  TFile* outFile = new TFile(outname.c_str(), "recreate") ;
  cout << "Output to files " << outname 
       << ", " << pdfname << ", and " << txtname << endl ;
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
  std::vector<string> tocList ; 
  
  //--- Histogram initialization ---//
  TH1D*              totalTime ;
  std::vector<TH1D*> pathTime ; 
  TH1D*              pathTimeSummary ; 
  std::vector<TH1D*> incPathTime ; 
  TH1D*              incPathTimeSummary ; 
  std::vector<TH1D*> moduleInPathTime ; 
  TH1D*              moduleTimeSummary ; 
  std::vector<TH1D*> moduleInPathScaledTime ; 
  TH1D*              pathSuccessFraction ;
  std::vector<TH1D*> failedModule ; 
  std::vector<TH1D*> moduleInPathRejection ; 
  TH1D*              pathRejection ; 
  std::vector<TH1D*> moduleInPathRejectTime ; 

  //--- Initial loop to get scale for histograms ---//
  double maxTime = 0. ;
  for (int ievt=0; ievt<n_evts; ievt++) {
      TBPerfInfo->GetEntry(ievt) ; 
      if ( HLTPerformance.totalTime() > maxTime ) maxTime = HLTPerformance.totalTime() ; 
  }
  double xmin = 0. ; 
  int xscale = 2 ; 
  while (pow(10,xscale) > maxTime) { xscale-- ; }
  double xmax = pow(10,xscale+3) ; 

  for (int ievt=0; ievt<n_evts; ievt++) {
      TBPerfInfo->GetEntry(ievt) ;

      //--- Initialization ---//
      int pCtr, mCtr ; 
      if (ievt == 0) {

          totalTime = new TH1D("totalTime","Total time for all modules per event",100,xmin,xmax) ;
          pathTimeSummary =
              new TH1D("pathTimeSummary","Average time per path",
                       HLTPerformance.numberOfPaths(),0.,double(HLTPerformance.numberOfPaths())) ; 
          incPathTimeSummary =
              new TH1D("incPathTimeSummary","Average incremental time per path",
                       HLTPerformance.numberOfPaths(),0.,double(HLTPerformance.numberOfPaths())) ; 
          pathSuccessFraction =
              new TH1D("pathSuccessFraction","Success rate (%) for each path",
                       HLTPerformance.numberOfPaths(),0.,double(HLTPerformance.numberOfPaths())) ; 
          pathRejection =
              new TH1D("pathRejection","Rejection for each path",
                       HLTPerformance.numberOfPaths(),0.,double(HLTPerformance.numberOfPaths())) ; 

          totalTime->GetXaxis()->SetTitle("msec") ; 
          pathTimeSummary->GetYaxis()->SetTitle("msec") ;
          incPathTimeSummary->GetYaxis()->SetTitle("msec") ;

          char name[30] ; 
          char title[100] ; 
          int ictr = 0 ; 
          for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
               pathIter!=HLTPerformance.endPaths(); pathIter++) {

              TH1D* histo ; 
              
              pathTimer.push_back(0.) ;
              incPathTimer.push_back(0.) ;
              pathSuccess.push_back(0) ;

              pathTimeSummary->GetXaxis()->SetBinLabel(ictr+1,pathIter->name().data()) ;
              incPathTimeSummary->GetXaxis()->SetBinLabel(ictr+1,pathIter->name().data()) ; 
              pathSuccessFraction->GetXaxis()->SetBinLabel(ictr+1,pathIter->name().data()) ; 
              pathRejection->GetXaxis()->SetBinLabel(ictr+1,pathIter->name().data()) ;
              
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
                  moduleInPathTime[ictr]->GetXaxis()->SetBinLabel(mCtr,modName) ;
                  moduleInPathScaledTime[ictr]->GetXaxis()->SetBinLabel(mCtr,modName) ;
                  moduleInPathRejection[ictr]->GetXaxis()->SetBinLabel(mCtr,modName) ;
                  moduleInPathRejectTime[ictr]->GetXaxis()->SetBinLabel(mCtr,modName) ;
                  failedModule[ictr]->GetXaxis()->SetBinLabel(++mCtr,modName) ;
              }

              ictr++ ;
              moduleInPathTimer.push_back( moduleInitTimer ) ;
              failures.push_back( moduleInitFailures ) ;
              moduleInitTimer.clear() ; 
              moduleInitFailures.clear() ; 
          }

          moduleTimeSummary =
              new TH1D("moduleTimeSummary","Average time per module",
                       HLTPerformance.numberOfModules(),0.,double(HLTPerformance.numberOfModules())) ;
          moduleTimeSummary->GetYaxis()->SetTitle("msec") ; 

          mCtr = 1 ; 
          for (HLTPerformanceInfo::Modules::const_iterator modIter=HLTPerformance.beginModules();
               modIter!=HLTPerformance.endModules(); modIter++) {
              moduleTimer.push_back(0.) ;
              moduleTimeSummary->GetXaxis()->SetBinLabel(mCtr++,modIter->name().c_str()) ; 
          }
      }

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
          } else { //--- One of the modules caused the path to fail ---//
              bool foundFailure = false ; 
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

      mCtr = 0 ;
      for (HLTPerformanceInfo::Modules::const_iterator modIter=HLTPerformance.beginModules();
           modIter!=HLTPerformance.endModules(); modIter++) {
          moduleTimer.at(mCtr++) += modIter->time() ;
      }

      totalTime->Fill( 1000. * HLTPerformance.totalTime() ) ; 
  }

  for (unsigned int pCtr=0; pCtr<moduleInPathTimer.size(); pCtr++) {
      int totalFailures = 0 ;
      int evtsLeft = n_evts ; 
      for (unsigned int mCtr=0; mCtr<moduleInPathTimer.at(pCtr).size(); mCtr++) { 

          double avgModuleTime = 1000. * moduleInPathTimer.at(pCtr).at(mCtr)/double(n_evts) ; // msec 
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
          double scaledAvgModuleTime = avgModuleTime * n_evts / evtsLeft ; 
          moduleInPathScaledTime.at(pCtr)->Fill( double(mCtr), scaledAvgModuleTime ) ; 
          moduleInPathRejectTime.at(pCtr)->Fill( double(mCtr),
                                                 modReject / scaledAvgModuleTime ) ; 
          
          evtsLeft -= failures.at(pCtr).at(mCtr) ; 
          totalFailures += failures.at(pCtr).at(mCtr) ; 

          if (failures.at(pCtr).at(mCtr) >= 0) {
              failedModule.at(pCtr)->Fill( double(mCtr), 100. * double(failures.at(pCtr).at(mCtr))/double(n_evts) );
          }
      }
      failedModule.at(pCtr)->Fill(-1.,100. * (double(n_evts-totalFailures)/double(n_evts))) ; 

      double pathReject = 1000. ;
      if (n_evts > totalFailures)
          pathReject = double(n_evts) / double(n_evts - totalFailures) ; 
      pathRejection->Fill(double(pCtr), pathReject) ; 
  }

          

  
  for (unsigned int i=0; i<pathTimer.size(); i++) {
      pathTimeSummary->Fill(double(i), pathTimer.at(i)/n_evts) ;
      incPathTimeSummary->Fill(double(i), incPathTimer.at(i)/n_evts) ;
      pathSuccessFraction->Fill(double(i), 100. * double(pathSuccess.at(i))/n_evts) ;

  }

  for (unsigned int i=0; i<moduleTimer.size(); i++) {
      moduleTimeSummary->Fill(double(i), 1000. * (moduleTimer.at(i)/n_evts)) ;
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

      gROOT->SetStyle("Plain") ; 
      gROOT->ForceStyle() ;
      c1->UseCurrentStyle() ; 

      pdf->SetTextSize(0.05) ; 
      pdf->Text(0.33,0.93,"Timing Summary Output") ; 
      
      pdf->SetTextSize(0.035) ; 
      pdf->Text(0.10,0.87,"Total time") ;
      pageNumber++ ; 
      sprintf(number,"%3i",pageNumber) ; 
      pdf->Text(0.90,0.87,number) ;
      pdf->SetTextSize(0.030) ;
      pdf->Text(0.15,0.84,"Cumulative time for all modules in the event") ; 
      string tocEntry = "Total time" ; 
      char tocPage[5] ; sprintf(tocPage,"%d",pageNumber) ;
      tocList.push_back( tocEntry + "^^" + string(tocPage) ) ; 

      pdf->SetTextSize(0.035) ; 
      pdf->Text(0.10,0.80,"Total path time") ; 
      pageNumber++ ;
      sprintf(number,"%3i",pageNumber) ; 
      pdf->Text(0.90,0.80,number) ;
      pdf->SetTextSize(0.030) ;
      pdf->Text(0.15,0.77,"Cumulative time for individual paths in the event") ; 
      
      tocEntry = "Total path time" ;
      sprintf(tocPage,"%d",pageNumber) ; 
      tocList.push_back( tocEntry + "^^" + string(tocPage) ) ; 
      for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
           pathIter!=HLTPerformance.endPaths(); pathIter++) {
          string subEntry = pathIter->name() ;
          sprintf(tocPage,"%d",pageNumber++) ;
          tocList.push_back( tocEntry + "^" + subEntry + "^^" + string(tocPage) ) ; 
      }
      
      pdf->SetTextSize(0.035) ; 
      pdf->Text(0.10,0.73,"Average path time") ; 
      sprintf(number,"%3i",pageNumber) ; 
      pdf->Text(0.90,0.73,number) ; 

      tocEntry = "Average path time" ; 
      sprintf(tocPage,"%d",pageNumber++) ;
      tocList.push_back( tocEntry + "^^" + string(tocPage) ) ; 
      

      pdf->Text(0.10,0.68,"Incremental path time") ; 
      sprintf(number,"%3i",pageNumber) ; 
      pdf->Text(0.90,0.68,number) ;
      pdf->SetTextSize(0.030) ;
      pdf->Text(0.15,0.65,"Incremental time due to unique modules in each path") ; 

      tocEntry = "Incremental path time" ;
      sprintf(tocPage,"%d",pageNumber) ; 
      tocList.push_back( tocEntry + "^^" + string(tocPage) ) ; 
      for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
           pathIter!=HLTPerformance.endPaths(); pathIter++) {
          string subEntry = pathIter->name() ;
          sprintf(tocPage,"%d",pageNumber++) ;
          tocList.push_back( tocEntry + "^" + subEntry + "^^" + string(tocPage) ) ; 
      }

      pdf->SetTextSize(0.035) ; 
      pdf->Text(0.10,0.61,"Average incremental path time") ; 
      sprintf(number,"%3i",pageNumber) ; 
      pdf->Text(0.90,0.61,number) ; 

      tocEntry = "Average incremental path time" ; 
      sprintf(tocPage,"%d",pageNumber++) ;
      tocList.push_back( tocEntry + "^^" + string(tocPage) ) ; 

      pdf->Text(0.10,0.56,"Average module time (by path)") ; 
      sprintf(number,"%3i",pageNumber) ; 
      pdf->Text(0.90,0.56,number) ;

      tocEntry = "Average module time (by path)" ;
      sprintf(tocPage,"%d",pageNumber) ; 
      tocList.push_back( tocEntry + "^^" + string(tocPage) ) ; 
      for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
           pathIter!=HLTPerformance.endPaths(); pathIter++) {
          string subEntry = pathIter->name() ;
          sprintf(tocPage,"%d",pageNumber++) ;
          tocList.push_back( tocEntry + "^" + subEntry + "^^" + string(tocPage) ) ; 
      }
      

      pdf->Text(0.10,0.51,"Average module time") ; 
      sprintf(number,"%3i",pageNumber) ; 
      pdf->Text(0.90,0.51,number) ; 

      tocEntry = "Average module time" ; 
      sprintf(tocPage,"%d",pageNumber++) ;
      tocList.push_back( tocEntry + "^^" + string(tocPage) ) ; 

      pdf->Text(0.10,0.46,"Average module RUNNING time (by path)") ; 
      sprintf(number,"%3i",pageNumber) ; 
      pdf->Text(0.90,0.46,number) ;
      pdf->SetTextSize(0.030) ;
      pdf->Text(0.15,0.43,"RUNNING time omits instances when module time is zero") ; 

      tocEntry = "Average module running time (by path)" ;
      sprintf(tocPage,"%d",pageNumber) ; 
      tocList.push_back( tocEntry + "^^" + string(tocPage) ) ; 
      for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
           pathIter!=HLTPerformance.endPaths(); pathIter++) {
          string subEntry = pathIter->name() ;
          sprintf(tocPage,"%d",pageNumber++) ;
          tocList.push_back( tocEntry + "^" + subEntry + "^^" + string(tocPage) ) ; 
      }

      pdf->SetTextSize(0.035) ; 
      pdf->Text(0.10,0.39,"Path success rate") ; 
      sprintf(number,"%3i",pageNumber) ; 
      pdf->Text(0.90,0.39,number) ; 

      tocEntry = "Path success rate" ;
      sprintf(tocPage,"%d",pageNumber++) ; 
      tocList.push_back( tocEntry + "^^" + string(tocPage) ) ; 
      
      pdf->Text(0.10,0.34,"Failing module (by path)") ; 
      sprintf(number,"%3i",pageNumber) ; 
      pdf->Text(0.90,0.34,number) ;
      pdf->SetTextSize(0.030) ; 
      pdf->Text(0.15,0.31,"Records failure rate as a function of modules in the path") ; 

      tocEntry = "Failing module (by path)" ;
      sprintf(tocPage,"%d",pageNumber) ; 
      tocList.push_back( tocEntry + "^^" + string(tocPage) ) ; 
      for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
           pathIter!=HLTPerformance.endPaths(); pathIter++) {
          string subEntry = pathIter->name() ;
          sprintf(tocPage,"%d",pageNumber++) ;
          tocList.push_back( tocEntry + "^" + subEntry + "^^" + string(tocPage) ) ; 
      }

      pdf->SetTextSize(0.035) ; 
      pdf->Text(0.10,0.27,"Module rejection factor (by path)") ; 
      sprintf(number,"%3i",pageNumber) ; 
      pdf->Text(0.90,0.27,number) ; 
      pdf->SetTextSize(0.030) ;
      pdf->Text(0.15,0.24,"Rejection factor = (Number of input events) / (Number of events that pass)") ; 
      
      tocEntry = "Module rejection factor (by path)" ;
      sprintf(tocPage,"%d",pageNumber) ; 
      tocList.push_back( tocEntry + "^^" + string(tocPage) ) ; 
      for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
           pathIter!=HLTPerformance.endPaths(); pathIter++) {
          string subEntry = pathIter->name() ;
          sprintf(tocPage,"%d",pageNumber++) ;
          tocList.push_back( tocEntry + "^" + subEntry + "^^" + string(tocPage) ) ; 
      }

      pdf->SetTextSize(0.035) ; 
      pdf->Text(0.10,0.20,"Path rejection factor") ; 
      sprintf(number,"%3i",pageNumber) ; 
      pdf->Text(0.90,0.20,number) ; 

      tocEntry = "Path rejection factor" ;
      sprintf(tocPage,"%d",pageNumber++) ; 
      tocList.push_back( tocEntry + "^^" + string(tocPage) ) ; 

      pdf->Text(0.10,0.15,"Module rejection factor per unit running time (by path)") ; 
      sprintf(number,"%3i",pageNumber) ; 
      pdf->Text(0.90,0.15,number) ; 
      pdf->SetTextSize(0.030) ;
      pdf->Text(0.15,0.12,"Use the module RUNNING time to compute Rejection/Time") ; 

      tocEntry = "Module rejection factor per unit running time (by path)" ;
      sprintf(tocPage,"%d",pageNumber) ; 
      tocList.push_back( tocEntry + "^^" + string(tocPage) ) ; 
      for (HLTPerformanceInfo::PathList::const_iterator pathIter=HLTPerformance.beginPaths();
           pathIter!=HLTPerformance.endPaths(); pathIter++) {
          string subEntry = pathIter->name() ;
          sprintf(tocPage,"%d",pageNumber++) ;
          tocList.push_back( tocEntry + "^" + subEntry + "^^" + string(tocPage) ) ; 
      }

      //--- Bookmarks output ---//
      for (unsigned int i=0; i<tocList.size(); i++) {
          string entry = tocList.at(i).append("\n") ;
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
