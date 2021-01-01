#!/usr/bin/env python
#------------------------------------------------------------------------------
#- File: pgammafit.py
#- Description: Example of use of PoissonGammaFit
#               Fit a mixture of QCD (QCD+Z), W+jets, top (tt, tqb, tb, tW)
#               to data. The counts for the sources are specified as
#               2-d histograms as are the data. The example is the expectation,
#               after cuts, for an integrated luminosity of 10/pb. The MC
#               samples were provided by Lukas V.
#- Created: 20-Sep-2010 Harrison B. Prosper
#           11-Feb-2011 HBP - remove dependence of kit.h
#  Updated: 03-Jun-2016 HBP - migrate to github for Sam
#           01-Jan-2020 HBP - update to Python 3
#------------------------------------------------------------------------------
import os, sys, re
import ROOT
from time import sleep
#------------------------------------------------------------------------------
sources = '''
QCD
W
top
'''
# Histogram names
sources = str.split(str.strip(sources))
hname   = ["h%s" % x for x in  sources]
hnameD  = "hdata"
getcount= re.compile(r'[a-zA-Z]+ +[0-9.]+')
vdouble = ROOT.vector("double")
#------------------------------------------------------------------------------
# Some Root utilities
#------------------------------------------------------------------------------
CSCALE = 1
ESCALE = 1

# unfold a 2-D histogram into 1-D
def unfoldContents(h):
    c = vdouble()
    e = vdouble()
    for i in range(h.GetNbinsX()):
        binx = i+1
        for j in range(h.GetNbinsY()):
            biny = j+1
            c.push_back( h.GetBinContent(binx, biny) * CSCALE )
            e.push_back( h.GetBinError(binx, biny) * CSCALE/ESCALE )
    return (c, e)

def setContents(h, c, e):
    for i in range(h.GetNbinsX()):
        bin =  i + 1
        h.SetBinContent(bin, c[i])
        h.SetBinError(bin, e[i])

def histogram(hname, xtitle, ytitle, nbins, xmin, xmax):
    h = ROOT.TH1F(hname, "", nbins, xmin, xmax)

    # Set some reasonable defaults
    h.SetLineColor(ROOT.kBlue)
    h.SetMarkerSize(0.8)
    h.SetMarkerColor(ROOT.kRed+2)
    h.SetMarkerStyle(20)
    
    h.GetXaxis().CenterTitle()
    h.GetXaxis().SetTitle(xtitle)
    h.GetXaxis().SetTitleOffset(1.3)
    h.SetNdivisions(510, "X")
    h.SetMarkerSize(1.0)

    h.GetYaxis().CenterTitle()
    h.GetYaxis().SetTitle(ytitle)
    h.GetYaxis().SetTitleOffset(1.4)
    h.SetNdivisions(504, "Y")
    return h

def setStyle():
    TEXTFONT=42
    style = ROOT.TStyle("SomeStyle", "pgammafit")
    style.SetPalette(1)

    #  For the canvas:
    style.SetCanvasBorderMode(0)
    style.SetCanvasColor(ROOT.kWhite)
    style.SetCanvasDefH(500) # Height of canvas
    style.SetCanvasDefW(500) # Width of canvas
    style.SetCanvasDefX(0)   # Position on screen
    style.SetCanvasDefY(0)

    #  For the Pad:
    style.SetPadBorderMode(0)
    style.SetPadColor(ROOT.kWhite)
    style.SetPadGridX(ROOT.kFALSE)
    style.SetPadGridY(ROOT.kFALSE)
    style.SetGridColor(ROOT.kGreen)
    style.SetGridStyle(3)
    style.SetGridWidth(1)

    #  For the frame:
    style.SetFrameBorderMode(0)
    style.SetFrameBorderSize(1)
    style.SetFrameFillColor(0)
    style.SetFrameFillStyle(0)
    style.SetFrameLineColor(1)
    style.SetFrameLineStyle(1)
    style.SetFrameLineWidth(1)

    #  For the histo:
    style.SetHistLineColor(1)
    style.SetHistLineStyle(0)
    style.SetHistLineWidth(2)

    style.SetEndErrorSize(2)
    style.SetErrorX(0.)

    style.SetMarkerSize(0.5)
    style.SetMarkerStyle(20)

    # For the fit/function:
    style.SetOptFit(1)
    style.SetFitFormat("5.4g")
    style.SetFuncColor(2)
    style.SetFuncStyle(1)
    style.SetFuncWidth(1)

    # For the date:
    style.SetOptDate(0)

    #  For the statistics box:
    style.SetOptFile(0)
    style.SetOptStat("")
    style.SetStatColor(ROOT.kWhite)
    style.SetStatFont(TEXTFONT)
    style.SetStatFontSize(0.03)
    style.SetStatTextColor(1)
    style.SetStatFormat("6.4g")
    style.SetStatBorderSize(1)
    style.SetStatH(0.2)
    style.SetStatW(0.3)

    #  Margins:
    style.SetPadTopMargin(0.05)
    style.SetPadBottomMargin(0.16)
    style.SetPadLeftMargin(0.20)
    style.SetPadRightMargin(0.10)

    #  For the Global title:
    style.SetOptTitle(0) 
    style.SetTitleFont(TEXTFONT)
    style.SetTitleColor(1)
    style.SetTitleTextColor(1)
    style.SetTitleFillColor(10)
    style.SetTitleFontSize(0.05)

    #  For the axis titles:
    style.SetTitleColor(1, "XYZ")
    style.SetTitleFont(TEXTFONT, "XYZ")
    style.SetTitleSize(0.05, "XYZ")
    style.SetTitleXOffset(1.25)    # (0.9)
    style.SetTitleYOffset(1.45)    # (1.25)

    #  For the axis labels:
    style.SetLabelColor(1, "XYZ")
    style.SetLabelFont(TEXTFONT, "XYZ")
    style.SetLabelOffset(0.007, "XYZ")
    style.SetLabelSize(0.05, "XYZ")

    #  For the axis:
    style.SetAxisColor(1, "XYZ")
    style.SetStripDecimals(ROOT.kTRUE)
    style.SetTickLength(0.03, "XYZ")
    style.SetNdivisions(510, "XYZ")
    #  To get tick marks on the opposite side of the frame
    style.SetPadTickX(1)  
    style.SetPadTickY(1)

    #  Change for log plots:
    style.SetOptLogx(0)
    style.SetOptLogy(0)
    style.SetOptLogz(0)

    #  Postscript options:
    style.SetPaperSize(20.,20.)
    style.cd()
    return style
#------------------------------------------------------------------------------
def main():

    # Load PoissonGammaFit class
    ROOT.gSystem.Load("../lib/libpgammafit")
    from ROOT import PoissonGammaFit
    
    # Read expected counts from the counts.txt file and place them in a map
    lines = [str.split(z) for z in getcount.findall(open("counts.txt").read())]
    counts= [(x[0], float(x[1])) for x in lines]
    
    count = {}
    for key, value in counts: count[key] = value
    
    # Open histogram file
    filename = "histograms.root"
    tfile = ROOT.TFile(filename)
    if not tfile.IsOpen(): sys.exit("** can't find %s" % filename)

    # Get histogram of observed counts and..
    hdata = tfile.Get("hdata")
    if not hdata: sys.exit("** can't find %s" % hnameD)

    # Plot it
    # Set some standard style and create a canvas
    setStyle()
    ROOT.gStyle.SetCanvasPreferGL(True)
    
    canvas = ROOT.TCanvas("fig_fitresults",
                     "#font[12]{m_{T}} vs #font[12]{b_{tag}}", 
                     50, 50, 500, 500)        
    canvas.cd()
    hdata.Draw("lego2 gl")
    canvas.Update()
    sleep(1)
    
    # 1. Unfold 2-d histogram into a 1-D vector, D, with uncertainties in dD
    D, dD = unfoldContents(hdata)
       
    # 2. Create 1-d histograms
    h1d = []
    h1d.append( histogram("hdata1d", "bin", "", len(D), 0, len(D)) )
    setContents(h1d[-1], D, dD)

    # 3. Plot it
    canvas.cd()
    h1d[-1].Draw("EP gl")
    canvas.Update()
    sleep(1)

    # Get distribution of sources and unfold into 1-d vectors, A[j] 
    h  = []
    A  = [] # vector<vector<double> > # counts and uncertainties
    
    for i, source in enumerate(sources):
        h.append( tfile.Get(hname[i]) ) # get histogram
        if not h[-1]: sys.exit("** can't find %s" % hname[i])
            
        # get counts and associated uncertainties
        A.append(unfoldContents(h[-1]))

        # plot just for fun
        canvas.cd()
        h[-1].Draw("lego2 gl")
        canvas.Update()
        sleep(1)

        h1d.append( histogram("%s1d" % hname[i],
                              hname[i], "", len(D), 0, len(D)) )
        h1d[-1].SetFillColor(i+1)              
        setContents(h1d[-1], A[i][0], A[i][1])
        
        canvas.cd()
        h1d[-1].Draw('hist')
        h1d[-1].Draw('EP same')
        canvas.Update()
        sleep(1)

    # Construct a PoissonGamma object 
    pgfit = PoissonGammaFit(D)

    # Add sources
    # scale the p_js so that the answers are counts
    scale = True
    for a, da in A: pgfit.add(a, da, scale) # add a source

    # Find mode of posterior density
    # make some reasonable guesses for fit
    total = sum(D)
    guess = vdouble(len(sources), total/3)
    pgfit.execute(guess)

    if not pgfit.good():
        sys.exit('''
        ** :( boo hoo - fit failed!
        ** quit or try better starting guesses for
        ** the source counts in the signal region and/or do some
        ** rebinning, maybe?
        ''')

    # Get mode and width of posterior density, which, because PoissonGammaFit
    # uses a flat prior for the yields, are the same as those of the
    # Poisson/gamma marginal density. 

    mode  = pgfit.mode()
    error = pgfit.width()
    logevidence = pgfit.logEvidence()

    # Print results
    print('')
    print("total observed count: %d" % total)   
    print("%-10s %10s\t%10s    %10s" %  \
          ('source', "true count", " estimated count", ""))
    for index, s in enumerate(sources):
        print("%-10s %10.f\t%10.0f +/-%-10.0f" % \
              (s, count[s], mode[index], error[index]))
#------------------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print("\n\t\tciao!\n")
