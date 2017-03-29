//--------------------------------------------------------------------
// File:    HistogramCache.cc
// Purpose: Cache histograms from a file and provide access to them by
//          name. The name must be provided in the format:
//                [dir1/[dir2/..]]histogramName
//
// Created: 11-Oct-2004 Harrison B. Prosper
//
// Modifications:
// 11/11/2004, RS+SJ: - put in protection against unknown key types
//                      when reading in histogram files.
// 
// 12/17/2004, S. Jain: - added the ability to Rebin 1D histograms, 
//                        and Normalize 1D, 2D, and 3D histograms
//                        RebinSize = 1 (default = 1) => No Rebinning
//                        Normalization < 0.0 (default = -1.0) => 
//                        No Normalization
// 
// 03/15/2005, RS+SJ: - added option for projecting a 2D histogram onto 
//                      one axis
// 08/08/2010, HBP - overload contents
//$Revision: 1.1.1.1 $
//---------------------------------------------------------------------
#include "TClass.h"
#include "TMath.h"
#ifdef PROJECT_NAME
#include "PhysicsTools/TheNtupleMaker/interface/HistogramCache.h"
#else
#include "HistogramCache.h"
#endif
//---------------------------------------------------------------------
using namespace std;

namespace {
  const int MAXDEPTH=20;
}
 
HistogramCache::HistogramCache(std::string filename, bool verbose)
  : _filename(filename), 
    _verbose(verbose),
    _os(new ostringstream)
{
  _file = new TFile(_filename.c_str());
  if (_file->IsZombie()) {
    cout << "Error opening file \"" <<_filename<<"\"."<<endl;
    exit(-1);
  }  
  _build(_file); // Build name to histogram map
}

HistogramCache::~HistogramCache()
{
  _file->Close();
}

bool
HistogramCache::exists(std::string histname) 
{
  return _hist.find(histname) != _hist.end();
}

vector<double>
HistogramCache::contents(std::string histname, 
                         int RebinSize, 
                         double Normalization,
                         std::string option) 
{
  vector<double> v;
  
  if ( ! exists(histname) ) 
    {
      cout<<"Error in file "<<_filename<<"!"<<endl;
      cout<<"Histogram \""<<histname<<"\" not found."<<endl;
      return v;
    }

  TH1* h = _hist[histname];
  
  if      ( h->IsA()->InheritsFrom("TH3") ) 
    {
      double bincontent;
      for (int binx=1; binx <= h->GetNbinsX(); binx++) 
        {
          for (int biny=1; biny <= h->GetNbinsY(); biny++) 
            {
              for (int binz=1; binz <= h->GetNbinsZ(); binz++)
                {
                  bincontent = h->GetBinContent(binx, biny, binz);
                  if(TMath::IsNaN(bincontent)) 
                    {
                      cout << "ERROR: In HistogramCache:" << endl;
                      cout << "Bin[" << binx<<", "<< biny<<", "
                           << binz<<"] has NaN."<<endl;
                      exit(0); 
                    } // if NaN
                  v.push_back(bincontent);
                } // loop over binz
            } // loop over biny
        } // loop over binx
      
      if (Normalization>=0.0) 
        {
          double total_count = 0.0;
          for (unsigned ibin = 0; ibin < v.size(); ibin++) 
            total_count +=  v[ibin];
          if (total_count > 0.0) 
            {
              for (unsigned ibin = 0; ibin < v.size(); ibin++) 
                v[ibin] *= Normalization / total_count;
            } 
          else 
            cout << "ERROR: Histogram has zero content !" << endl;
        } // Normalize the histogram
    } // if 3D histogram	      
  else if ( h->IsA()->InheritsFrom("TH2") ) 
    {
      // check options:
      bool doProjectionX = option=="ProjectionX";
      bool doProjectionY = option=="ProjectionY";
      
      if(doProjectionX) 
        {
          // use a projection on the X axis
          TH2* h2dtemp=(TH2*) h;
          TH1D *htemp=h2dtemp->ProjectionX("");
          return contents(htemp, RebinSize, Normalization);
        }
      else if(doProjectionY) 
        {
          TH2* h2dtemp=(TH2*) h;
          TH1D *htemp=h2dtemp->ProjectionY("");
          return contents(htemp, RebinSize, Normalization);
        }
      else 
        {
          // no projections, take all the bins!
          double bincontent;
          for (int binx=1; binx <= h->GetNbinsX(); binx++)
            {
              for (int biny=1; biny <= h->GetNbinsY(); biny++)
                {
                  bincontent = h->GetBinContent(binx, biny);
                  if(TMath::IsNaN(bincontent)) 
                    {
                      cout << "ERROR: In HistogramCache:" << endl;
                      cout << "Bin[" << binx<<", "
                           << biny<<"] has NaN."<<endl;
                      exit(0); 
                    } // if NaN
                  v.push_back(bincontent);
                }// loop over biny
            } // loop over binx
          if (Normalization>=0.0) 
            {
              double total_count = 0.0;
              for (unsigned ibin = 0; ibin < v.size(); ibin++) 
                total_count +=  v[ibin];
              if (total_count > 0.0) 
                {
                  for (unsigned ibin = 0; ibin < v.size(); ibin++) 
                    v[ibin] *= Normalization / total_count;
                } 
              else 
                cout << "ERROR: Histogram has zero content !" 
                     << endl;
            } // Normalize the histogram
        }
    } // if 2D histogram
  else 
    {
      // do this for a 1d histogram
      return contents(h, RebinSize, Normalization);
    }     // else if 1D histogram)  
  return v;
}

// helper function to return the contents of a 1d histogram
vector<double>
HistogramCache::contents(TH1 *h, 
                         int RebinSize, 
                         double Normalization) 
{
  vector<double> v;
  double sumRebin; 
  double bincontent;
  for (int binx=1; binx <= h->GetNbinsX();)
    {  
      sumRebin = 0.0 ;
      for (int i=0; i < RebinSize; i++)
        {
          bincontent = h->GetBinContent(binx+i);
          if(TMath::IsNaN(bincontent)) 
            {
              cout << "ERROR: In HistogramCache:" << endl;
              cout << "Bin[" << binx+i<<"] has NaN."<<endl;
              exit(0); 
            } // if NaN
          sumRebin += bincontent;
        } // loop over binx

      v.push_back(sumRebin);
      binx += RebinSize;
    } // loop over bins
  
  if (Normalization>=0.0) 
    {
      double total_count = 0.0;
      for (unsigned ibin = 0; ibin < v.size(); ibin++) 
        total_count +=  v[ibin];
      if (total_count > 0.0) 
        {
          for (unsigned ibin = 0; ibin < v.size(); ibin++) 
            v[ibin] *= Normalization / total_count;
        } 
      else 
        cout << "ERROR: Histogram has zero content !" << endl;
    } // Normalize the histogram
  return v;
}

TH1*
HistogramCache::histogram(std::string histname) 
{
  if ( exists(histname) )
    return _hist[histname];
  else
    return (TH1*)0;
}

string
HistogramCache::type(std::string histname) 
{
  if ( exists(histname) )
    return string(_hist[histname]->ClassName());
  else
    return string("");
}

void 
HistogramCache::ls(ostream& os)
{
  os << _os->str();
}

string
HistogramCache::str()
{
  return _os->str();
}


//------------------------------------------------------
// Private methods
//------------------------------------------------------
void 
HistogramCache::_check(int depth)
{
  if ( depth > MAXDEPTH )
    {
      std::cerr << "**Warning** I'm lost in the trees...goodbye!" 
                << std::endl;
      exit(1);
    }
}

void
HistogramCache::_build(TDirectory* dir, int depth)
{
  depth++;
  _check(depth);
  
  string tab(2*depth, ' ');

  TIter nextkey(dir->GetListOfKeys());
  while ( TKey *key = (TKey*)nextkey() )
    {
      dir->cd();
      
      TObject* o = key->ReadObj();
      
      if ( o->IsA()->InheritsFrom("TDirectory") )
        {
          TDirectory* d = (TDirectory*)o;
          *_os << tab << "BEGIN " << d->GetName() << endl;
          _path.push(d->GetName());
          _build(d, depth);
          *_os << tab << "END " << d->GetName() << endl;
          if ( ! _path.empty() ) _path.pop();
        }
      // Note: All histograms inherit from TH1
      else if ( o->IsA()->InheritsFrom("TH1") )
        {
          TH1* h = (TH1*)o;
          *_os << tab << o->ClassName() << "\t" << h->GetName() << endl;
          string nkey = _name(h);
          _hist[nkey] = h;
          
          if ( _verbose )
            { 
              cout << nkey;
              if ( h->IsA()->InheritsFrom("TH3") )
                cout << "\t3-D";
              else if ( h->IsA()->InheritsFrom("TH2") )
                cout << "\t2-D";
              else
                cout << "\t1-D";
              cout << endl;
            }
        }
    } // end of loop over keys
}

string
HistogramCache::_name(TH1* h)
{
  _path.push(h->GetName());
  stack<string> tmp = _path;
  string delim("");
  string name("");
  while ( ! tmp.empty() )
    {
      name = tmp.top() + delim + name;
      delim= string("/");
      tmp.pop(); //IMPORTANT!
    }
  _path.pop(); 
  return name;
}

#ifdef __TEST__
//-----------------------------------------------------------
// Test program
//-----------------------------------------------------------
int main(int argc, char** argv)
{
  string filename;
  if ( argc > 1 )
    filename = string(argv[1]);
  else
    filename = string("data.root");

  TFile* file = new TFile(filename.c_str());
  if ( ! file ) 
    {
      cout << "**Error** Unable to open file " << filename 
           << endl;
      exit(1);
    }

  HistogramCache cache(file);

  ofstream out("list.txt");
  cache.ls(out);
  out.close();

  vector<double> d1d = cache.contents("01_Create_W_e_VOID/MissingEt_Pt");
  //  cout << "1-D Histogram with " << d1d.size() << " bins" << endl;
  // for(int bin=0; bin < (int)d1d.size(); bin=bin+10)
  //    cout << "\t" << bin+1 << "\t" << d1d[bin] << endl;

  vector<double> d2d = cache.
    contents("04_CommonHistos/08_Histogramming_D0Jet_good_jet/Pt_vs_Ntracks");
  cout << "2-D Histogram with " << d2d.size() << " bins" << endl;
  for(int bin=0; bin < (int)d2d.size(); bin++)
    cout << "\t" << bin+1 << "\t" << d2d[bin] << endl;

  file->Close();
}
#endif

