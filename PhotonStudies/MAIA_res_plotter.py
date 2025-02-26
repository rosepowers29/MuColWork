import numpy as np
from array import array
import matplotlib.pyplot as plt
import mplhep as hep
import uproot
import ROOT
from ROOT import TGraphErrors, TF1
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import pandas as pd

#MAIA format python plots (!!)



#define the fit function
def py_resofit(stoch, noise, const, E):
    return np.sqrt((stoch**2)/E+(noise/E)**2+const**2)

Epoints = np.linspace(0.01,1000.,1000)
resopoints_B=[]
resopoints_nB = []


def py_resofit_simple(a, b, E):
    return np.sqrt((a**2)/E+(b/E)**2)

#define the root fit
def rootFit(nb_arrays, b_arrays):
    e_arr_nB, sigma_arr_nB, sigma_err_arr_nB, e_err_arr_nB = nb_arrays
    e_arr_B, sigma_arr_B, sigma_err_arr_B, e_err_arr_B = b_arrays
    gr_reso_E_B = TGraphErrors(len(e_arr_B), e_arr_B, sigma_arr_B, e_err_arr_B, sigma_err_arr_B)
    gr_reso_E_nB = TGraphErrors(len(e_arr_nB), e_arr_nB, sigma_arr_nB, e_err_arr_nB, sigma_err_arr_nB)

    resoFit_B = TF1(    "resofit", "sqrt([0]*[0]/x+[1]*[1]/(x*x))", 150., 1000., 3) #"sqrt([0]*[0]/x+[1]*[1]/(x*x)+[2]*[2])", 30., 1000., 3)
    resoFit_B.SetParName(0, "Stochastic")
    resoFit_B.SetParName(1, "Noise")
    resoFit_B.SetParameter(0, .5)
    resoFit_B.SetParameter(1, 1.)
    gr_reso_E_B.Fit(resoFit_B)

    resoFit_nB = TF1(    "resofitnB", "sqrt([0]*[0]/x+[1]*[1]/(x*x))", 30., 1000., 3) #"sqrt([0]*[0]/x+[1]*[1]/(x*x)+[2]*[2])", 30., 1000., 3)
    resoFit_nB.SetParName(0, "Stochastic")
    resoFit_nB.SetParName(1, "Noise")
    resoFit_nB.SetParameter(0, 1.)
    resoFit_nB.FixParameter(1, 0.)
    gr_reso_E_nB.Fit(resoFit_nB)

    a_B = resoFit_B.GetParameter(0)
    b_B = resoFit_B.GetParameter(1)

    a_nB = resoFit_nB.GetParameter(0)
    b_nB = resoFit_nB.GetParameter(1)

    nB_params = [a_nB, b_nB]
    B_params = [a_B, b_B]
    
    #clear the arrays
    resopoints_B.clear()
    resopoints_nB.clear()

    for energy in Epoints:
        resopoints_B.append(py_resofit_simple(a_B, b_B, energy))
        resopoints_nB.append(py_resofit_simple(a_nB, b_nB, energy))
    

    return(nB_params, B_params, resopoints_nB, resopoints_B)


#plotting functions to easily switch between regions, variables of interest
def plot_res(nb_arrays, b_arrays, region, nB_params, B_params, resopoints_nB, resopoints_B):
    #region: 1=endcap, 2=barrel, 3=transition, 0 = null (theta res)

    #unpack values
    e_arr_nB, sigma_arr_nB, sigma_err_arr_nB, e_err_arr_nB = nb_arrays
    e_arr_B, sigma_arr_B, sigma_err_arr_B, e_err_arr_B = b_arrays
    a_nB, b_nB = nB_params
    a_B, b_B = B_params
    
    fig, ax = plt.subplots()
    plt.errorbar(e_arr_B, sigma_arr_B, yerr = sigma_err_arr_B, xerr = e_err_arr_B, fmt=".", color='red',
        label = "With BIB overlay")
    plt.errorbar(e_arr_nB, sigma_arr_nB, yerr= sigma_err_arr_nB, xerr = e_err_arr_nB, fmt = ".", color = 'blue',
        label = "Without BIB overlay")
    if region != 0:
        plt.plot(Epoints, resopoints_B, '--', c='red',lw=0.25)#, label="Fit")
        plt.plot(Epoints, resopoints_nB, '--', c='blue',lw=0.25)
        plt.xlim(30., 800.)
        plt.ylim(0.,1.)
        plt.gcf().text(0.65, 0.65, "$\\frac{\sigma}{E} = \sqrt{\\frac{a^2}{E}+\\frac{b^2}{E^2}}$")
        plt.gcf().text(0.65, 0.6, "Fit Results (no BIB):")
        plt.gcf().text(0.65, 0.55, "$a=$"+str(round(abs(a_nB),2)),color='blue')
        plt.gcf().text(0.65, 0.5, "$b=$"+str(round(abs(b_nB),2)), color='blue')
        plt.gcf().text(0.65, 0.45, "Fit Results (with BIB):")
        plt.gcf().text(0.65, 0.4, "$a=$"+str(round(abs(a_B),2)),color='red')
        plt.gcf().text(0.65, 0.35, "$b=$"+str(round(abs(b_B),2)), color='red')
        
    else:
        plt.xlim(.18, 2.96)
        plt.ylim(0., 0.15)

    if region == 0:
        plt.xlabel("True Photon $\\theta$ [rad]", loc='right')
        plt.gcf().text(0.16, 0.7, "True Photon Energy  > 600 GeV")
    else:
        plt.xlabel("True Photon Energy [GeV]", loc = 'right')
    plt.ylabel("Photon Energy Resolution $\sigma((E_{reco}-E_{true})/E_{true}$)", loc='top')
    hep.cms.label(exp = "Muon Collider", data = False, label = 'with BIB (EU24 Lattice) ',
       rlabel='$MAIA$ Detector Concept', loc=2, italic=(1,0,0), pad=(0.0))
    plt.gcf().text(0.16,0.74,"$\sqrt{s}=10$ TeV")
    if region == 3:
        plt.gcf().text(0.16,0.7,"Transition Region ($0.577<\\theta<1.0$ or $2.0<\\theta<2.56$)")
    elif region == 2:
        plt.gcf().text(0.16, 0.7, "Central Barrel Region ($1.0<\\theta<2.0$)")
    elif region == 1:
        plt.gcf().text(0.16,0.7, "Endcap Region ($\\theta < 0.577$ or $\\theta > 2.56$)")
    
    
    ax.tick_params(bottom=True, top=True, left=True, right=True, which='both', direction='in')
    ax.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)
    ax.tick_params(axis= 'y',which='major',length = 10)
    ax.tick_params(axis='y',which='minor', length=5)

    if region == 0:
        ax.xaxis.set_major_locator(MultipleLocator(0.5))
        ax.xaxis.set_major_formatter('{x:.1f}')
        ax.xaxis.set_minor_locator(MultipleLocator(0.1))
        
        ax.yaxis.set_major_locator(MultipleLocator(0.02))
        ax.yaxis.set_major_formatter('{x:.2f}')
        ax.yaxis.set_minor_locator(MultipleLocator(0.005))
    else:
        ax.xaxis.set_major_locator(MultipleLocator(200))
        ax.xaxis.set_major_formatter('{x:.0f}')
        ax.xaxis.set_minor_locator(MultipleLocator(50))

        ax.yaxis.set_major_locator(MultipleLocator(0.1))
        ax.yaxis.set_major_formatter('{x:.1f}')
        ax.yaxis.set_minor_locator(MultipleLocator(0.025))

    plt.grid(axis = 'y', linestyle = ':')
    plt.legend(loc='upper right',frameon=False)
    plt.show()
    var_region = ''
    if region == 0:
        var_region = "th"
    elif region == 1:
        var_region = "E_endcaps_"
    elif region == 2:
        var_region = "E_barrel_"
    elif region == 3:
        var_region = "E_transition_"
    title = "Ereso_v0.8_maia_"+var_region+"v1.pdf"
    plt.savefig(title)
    print(title,"created")
    plt.close()
    return()

#load in data
BIB_csvs = ["BIB_resarrays_th.csv","BIB_resarrays_e.csv", "BIB_resarrays_b.csv", "BIB_resarrays_t.csv"]
noBIB_csvs = ["noBIB_resarrays_th.csv", "noBIB_resarrays_e.csv", "noBIB_resarrays_b.csv", "noBIB_resarrays_t.csv"]
csvs = dict(zip(BIB_csvs, noBIB_csvs))

i = 0
for Bcsv, nBcsv in csvs.items():
    BIBdf = pd.read_csv(Bcsv)
    noBIBdf = pd.read_csv(nBcsv)
    print(i)
    e_arr_B = BIBdf["E"].to_numpy()
    e_err_arr_B = BIBdf["E_err"].to_numpy()
    sigma_arr_B = BIBdf["sigma"].to_numpy()
    sigma_err_arr_B = BIBdf["sigma_err"].to_numpy()
    print(sigma_arr_B)
    e_arr_nB = noBIBdf["E"].to_numpy()
    e_err_arr_nB = noBIBdf["E_err"].to_numpy()
    sigma_arr_nB = noBIBdf["sigma"].to_numpy()
    sigma_err_arr_nB = noBIBdf["sigma_err"].to_numpy()

    #package up the arrays
    B_arrs = [e_arr_B, sigma_arr_B, sigma_err_arr_B, e_err_arr_B]
    nB_arrs = [e_arr_nB, sigma_arr_nB, sigma_err_arr_nB, e_err_arr_nB]

    #call the functions now

    nB_params, B_params, resopoints_nB, resopoints_B = rootFit(nB_arrs, B_arrs)
    print(len(resopoints_B))
    print(nB_arrs[1])
    plot_res(nB_arrs, B_arrs, i, nB_params, B_params, resopoints_nB, resopoints_B)

    #advance the index of regions
    i += 1

print("done")
