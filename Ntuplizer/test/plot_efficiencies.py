import ROOT as R
import numpy as np
from tqdm import tqdm
from argparse import ArgumentParser
import os
R.gROOT.ProcessLine('.L ./tdrstyle.C')
R.gROOT.SetBatch(1)
R.setTDRStyle()
R.gStyle.SetPaintTextFormat("4.2f")

_vars = { 
    "dmu_dsa_dxy":   {"name": r"tag |d_0| (cm)", "nbins": 20,  "xmin": 0.,   "xmax": 200.},
    "dmu_dsa_pt":    {"name": r"tag p_T (GeV)", "nbins": 30,  "xmin": 0.,   "xmax": 90.},
    "dmu_dsa_eta":   {"name": r"tag #eta", "nbins": 20,  "xmin": -1.,   "xmax": 1.},
    "dmu_dgl_dxy":   {"name": r"tag |d0| (cm)", "nbins": 6,  "xbins": np.array([0.,2.,5.,10.,30.,50.,70.],dtype=float)},
    "dmu_dgl_dz":    {"name": r"tag |dz| (cm)", "nbins": 6,  "xbins": np.array([0.,8.,20.,40.,60.,90.,140.],dtype=float)},
}

parser = ArgumentParser()
parser.add_argument("--var", type=str, nargs='*', choices=_vars.keys(), help='Variable(s) to plot')
parser.add_argument("--mcfile",   type=str, help='MC Ntuple file', default="CosmicsMC_Leonardo_Ntuples.root")
parser.add_argument("--datafile", type=str, help='Data Ntuple file', default="CosmicsData_Ntuples.root")
parser.add_argument("--tag", type=str, help='Tag to add in plots name', default="")
args = parser.parse_args()

if __name__=="__main__":

    datafilename = args.datafile
    mcfilename   = args.mcfile
    datafile = R.TFile.Open(datafilename,"READ")
    datatree = datafile.Get("Events")
    mcfile   = R.TFile.Open(mcfilename,"READ")
    mctree   = mcfile.Get("Events")

    if len(args.var)==1: 
        var = args.var[0]

        nbins = _vars[var]["nbins"]
        xmin  = _vars[var]["xmin"]
        xmax  = _vars[var]["xmax"]
        h_dummy = R.TH1F("h_dummy", ";"+_vars[var]["name"]+";Efficiency", nbins, xmin, xmax)
        h_data_total = R.TH1F("h_data_total", ";"+_vars[var]["name"]+";Efficiency", nbins, xmin, xmax)
        h_data_pass  = R.TH1F("h_data_pass",  ";"+_vars[var]["name"]+";Efficiency", nbins, xmin, xmax)
        h_mc_total = R.TH1F("h_mc_total", ";"+_vars[var]["name"]+";Efficiency", nbins, xmin, xmax)
        h_mc_pass  = R.TH1F("h_mc_pass",  ";"+_vars[var]["name"]+";Efficiency", nbins, xmin, xmax)
        h_data_total.Sumw2()
        h_data_pass.Sumw2()
        h_mc_total.Sumw2()
        h_mc_pass.Sumw2()


        total  = "dmu_isDSA && dmu_dsa_passTagID"
        passed = "dmu_isDSA && dmu_dsa_passTagID && dmu_dsa_hasProbe"
        datatree.Draw(f"{var}>>h_data_total", total, "goff")
        datatree.Draw(f"{var}>>h_data_pass", passed, "goff")
        mctree.Draw(f"{var}>>h_mc_total", total, "goff")
        mctree.Draw(f"{var}>>h_mc_pass", passed, "goff")

        h_data_eff = R.TEfficiency(h_data_pass, h_data_total)
        h_mc_eff   = R.TEfficiency(h_mc_pass, h_mc_total)

        h_data_eff.SetMarkerStyle(8)
        h_data_eff.SetMarkerColor(R.kRed)
        h_data_eff.SetLineColor(R.kRed)
        h_data_eff.SetMarkerSize(0.6)
        h_mc_eff.SetMarkerStyle(8)
        h_mc_eff.SetMarkerColor(R.kRed+2)
        h_mc_eff.SetLineColor(R.kRed+2)
        h_mc_eff.SetMarkerSize(0.6)
        
        c = R.TCanvas("c_eff",";"+_vars[var]["name"]+";Efficiency")
        c.cd()

        # PAD 1
        pad1 = R.TPad("pad1", "pad1", 0, 0.3, 1, 1.0) # for the plot
        pad1.SetBottomMargin(0.015)
        pad1.SetTopMargin(0.1)
        pad1.Draw()
        # PAD2
        pad2 = R.TPad("pad2", "pad2", 0, 0.01, 1, 0.3) # for the difference
        pad2.SetTopMargin(0.05);
        pad2.SetBottomMargin(0.4);
        pad2.Draw()

        pad1.cd()
        h_dummy.SetMinimum(0.)
        h_dummy.SetMaximum(1.1)
        h_dummy.GetXaxis().SetLabelSize(0)
        h_dummy.GetYaxis().SetLabelSize(0.04)
        h_dummy.GetYaxis().SetTitleSize(0.04)
        h_dummy.Draw("HIST")
        h_data_eff.Draw("P,SAME")
        h_mc_eff.Draw("P,SAME")

        pad1.Update()
        pad1.RedrawAxis()

        pad2.cd()
        h_ratio = h_data_pass.Clone()
        h_ratio.Sumw2()
        h_ratio.GetYaxis().SetTitle("Data/MC");
        h_ratio.GetYaxis().CenterTitle();
        h_ratio.GetYaxis().SetLabelSize(0.10);
        h_ratio.GetYaxis().SetNdivisions(6);
        h_ratio.GetYaxis().SetTitleSize(0.10);
        h_ratio.GetXaxis().SetLabelSize(0.10);
        h_ratio.GetXaxis().SetTitleSize(0.10);
        h_ratio.GetXaxis().SetLabelOffset(0.03);
        h_ratio.GetXaxis().SetTitleOffset(1.5);
        h_ratio.Multiply(h_mc_total)
        h_ratio.Divide(h_data_total)
        h_ratio.Divide(h_mc_pass)
        h_ratio.SetMarkerStyle(8)
        h_ratio.SetMarkerColor(R.kBlack)
        h_ratio.SetLineColor(R.kBlack)
        h_ratio.SetMarkerSize(0.6)
        h_ratio.SetMaximum(1.2)
        h_ratio.SetMinimum(0.8)
        xmin = h_ratio.GetBinLowEdge(1)
        xmax = h_ratio.GetBinLowEdge(h_ratio.GetNbinsX()+1)
        line = R.TLine(xmin, 1, xmax, 1)
        line.SetLineColor(R.kGray+2)
        line.SetLineWidth(2)
        h_ratio.Draw("PE")
        line.Draw()
        
        pad2.Update()
        pad2.RedrawAxis()
        aux_frame2 = R.TLine()
        aux_frame2.SetLineWidth(2) 
        aux_frame2.DrawLine(pad2.GetUxmax(), pad2.GetUymin(), pad2.GetUxmax(), pad2.GetUymax())

        c.cd()
        leg = R.TLegend(0.14,0.8,0.5,0.86)
        leg.AddEntry(h_data_eff, "Data", "P")
        leg.AddEntry(h_mc_eff, "MC", "P")
        leg.SetFillStyle(0)
        leg.SetLineWidth(0)
        leg.SetTextFont(42)
        leg.Draw()

        latex = R.TLatex()
        latex.SetNDC();
        latex.SetTextAngle(0);
        latex.SetTextColor(R.kBlack);
        latex.SetTextFont(42);
        latex.SetTextAlign(11);
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.12, 0.94, "#bf{CMS} #it{Preliminary}")
        latex.DrawLatex(0.73, 0.94, "13.6 TeV")
        latex.SetTextSize(0.03);
        latex.DrawLatex(0.23, 0.44, datafilename.split('/')[-1])

        plotname = f"ploteff_{var}_{args.tag}"
        c.SaveAs(f"{plotname}.png")
        c.SaveAs(f"{plotname}.pdf")
        c.SaveAs(f"{plotname}.C")

    elif len(args.var)==2:
        R.gStyle.SetPadGridX(False)
        R.gStyle.SetPadGridY(False)
        R.gStyle.SetTitleYOffset(1.5)
        varx = args.var[0]
        vary = args.var[1]

        nbinsx = _vars[varx]["nbins"]
        nbinsy = _vars[vary]["nbins"]
        xbinsx = _vars[varx]["xbins"]
        xbinsy = _vars[vary]["xbins"]
        h_dummy = R.TH2F("h_dummy", r";"+_vars[varx]["name"]+";"+_vars[vary]["name"], nbinsx, xbinsx, nbinsy, xbinsy)
        h_data_total = R.TH2F("h_data_total", r"Cosmics Data Run2023C;"+_vars[varx]["name"]+";"+_vars[vary]["name"], nbinsx, xbinsx, nbinsy, xbinsy)
        h_data_pass  = R.TH2F("h_data_pass",  r"Cosmics Data Run2023C;"+_vars[varx]["name"]+";"+_vars[vary]["name"], nbinsx, xbinsx, nbinsy, xbinsy)
        h_mc_total   = R.TH2F("h_mc_total",   r"Cosmics MC;"+_vars[varx]["name"]+";"+_vars[vary]["name"], nbinsx, xbinsx, nbinsy, xbinsy)
        h_mc_pass    = R.TH2F("h_mc_pass",    r"Cosmics MC;"+_vars[varx]["name"]+";"+_vars[vary]["name"], nbinsx, xbinsx, nbinsy, xbinsy)
        h_data_total.Sumw2()
        h_data_pass.Sumw2()
        h_mc_total.Sumw2()
        h_mc_pass.Sumw2()

        total  = "dmu_isDGL && dmu_dgl_passTagID"
        passed = "dmu_isDGL && dmu_dgl_passTagID && dmu_dgl_hasProbe"
        datatree.Draw(f"{vary}:{varx}>>h_data_total", total, "goff")
        datatree.Draw(f"{vary}:{varx}>>h_data_pass", passed, "goff")
        mctree.Draw(f"{vary}:{varx}>>h_mc_total", total, "goff")
        mctree.Draw(f"{vary}:{varx}>>h_mc_pass", passed, "goff")

        h_data_eff = R.TEfficiency(h_data_pass, h_data_total)
        h_mc_eff   = R.TEfficiency(h_mc_pass, h_mc_total)

        c = R.TCanvas("c","",1200,600)
        c.Divide(2,1)
        c.cd(1)
        h_data_eff.Draw("COLZ,TEXT")
        latex = R.TLatex()
        latex.SetNDC();
        latex.SetTextAngle(0);
        latex.SetTextColor(R.kBlack);
        latex.SetTextFont(42);
        latex.SetTextAlign(11);
        latex.SetTextSize(0.03);
        latex.DrawLatex(0.12, 0.9, "#bf{CMS} #it{Preliminary}"+f" (Cosmics {datafilename.split('_')[1].split('-')[1]})")
        latex.DrawLatex(0.73, 0.9, "13.6 TeV")
        c.Update()
        c.cd(2)
        h_mc_eff.Draw("COLZ,TEXT")
        latex = R.TLatex()
        latex.SetNDC();
        latex.SetTextAngle(0);
        latex.SetTextColor(R.kBlack);
        latex.SetTextFont(42);
        latex.SetTextAlign(11);
        latex.SetTextSize(0.03);
        latex.DrawLatex(0.12, 0.9, "#bf{CMS} #it{Simulation}")
        latex.DrawLatex(0.73, 0.9, "13.6 TeV")
        c.Update()
        plotname = f"ploteff_{varx}_{vary}_{args.tag}"
        c.SaveAs(f"{plotname}.png")
        c.SaveAs(f"{plotname}.pdf")
        c.SaveAs(f"{plotname}.C")
