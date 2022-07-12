import matplotlib
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)
matplotlib.rc('axes', labelsize=18)
matplotlib.rcParams['lines.linewidth'] = 3

#matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
matplotlib.rc('font',**{'family':'serif','serif':['Palatino']})
matplotlib.rc('text', usetex=True)
legendfontsize = 13

import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import numpy as np
import units
import fluxes
import os
import E2_fluxes_NuTauSnowmass as nutauplot

# YOU NEED TO CHANGE THIS IN E2_fluxes_NuTauSnowmass,
# expdata.py AND modeldata.py 
# TO UPDATE PROPERLY
energyBinsPerDecade = 1.
plotUnitsEnergy = units.eV
plotUnitsEnergyStr = "eV"
plotUnitsFlux = units.GeV * units.cm ** -2 * units.second ** -1 * units.sr ** -1
flavorRatio = 1./3.
DIFFUSE = True
livetime = 10 * units.year
subtitles = ['Embedded in Ice', "Valleys and Mountains", "Balloons and Satellites"]

#subtitles = ['Water and Ice', "Atmosphere, Earth's Limb, Topography", "$\nu_\tau$ Only"]

# identify which subplot to include
show_icecube = [0,1,2]
show_icecube_uhe = [0]
show_anita = [2]
show_pueo = [2]
show_poemma = [2]
show_eusospb =[2]
show_auger = [1]
show_ara = [0]
show_arianna = [0]
show_grand = [1]
show_beacon = [1]
show_taroge = [1]
show_trinity = [1]
show_tambo = [1]
show_radar = [0]
show_rnog = [0]
show_gen2 = [0]

fig, axs = plt.subplots(3, 1, figsize=(12, 12))
second_axs = []     
for i in range(0,3):
    nutauplot.get_E2_limit_figure(fig=fig, ax=axs[i],show_model_legend=True,
                                      diffuse=DIFFUSE,
                                                              show_ice_cube_EHE_limit=i in show_icecube_uhe,
                                                              show_ice_cube_HESE_data=False,
                                                              show_ice_cube_HESE_fit=i in show_icecube,
                                                              show_ice_cube_mu=i in show_icecube,
                                                              nu_mu_show_data_points=False,
                                                              show_ice_cube_mu_extrap=False,
                                                              show_icecube_glashow=False,
                                                              show_anita_I_III_limit=False,
                                                              show_anita_I_IV_limit=i in show_anita,
                                                              show_pueo30=i in show_pueo,
                                                              show_pueo100=i in show_pueo,
                                                              show_poemma=i in show_poemma,
                                                              show_poemma360=i in show_poemma,
                                                              show_poemma_fluor=i in show_poemma,
                                                              show_eusospb=i in show_eusospb,
                                                              show_auger_limit=i in show_auger,
                                                              show_ara=i in show_ara,
                                                              show_ara_2023=False,
                                                              show_ara_2023_TL=False,
                                                              show_arianna=i in show_arianna,
                                                              show_grand_10k=False,
                                                              show_grand_200k=i in show_grand,
                                                              show_beacon=i in show_beacon,
                                                              show_taroge=i in show_taroge,
                                                              show_tambo=i in show_tambo,
                                                              show_trinity=i in show_trinity,
                                                              show_ska=False,
                                                              show_radar=i in show_radar,
                                                              show_RNOG=i in show_rnog,
                                                              show_IceCubeGen2_whitepaper=False,
                                                              show_IceCubeGen2_ICRC2021=i in show_gen2,
                                                              show_IceCubeGen2_proj=i in show_gen2,
                                                              show_ara_1year=False,
                                                              show_prediction_arianna_200=False,
                                                              show_Heinze=False,
                                                              show_Auger_vanvliet=True,
                                                              show_TA=False,
                                                              show_TA_nominal=False,
                                                              show_TA_ICRC2021=False,
                                                              show_neutrino_best_fit=False,
                                                              show_neutrino_best_case=False,
                                                              show_neutrino_worst_case=False,
                                                              show_muf_bestfit=True,
                                                              show_astro=True)
    axs[i].set_yscale('log')
    axs[i].set_xscale('log')
    axs[i].set_ylabel(r'Tau Neutrino $E^2\Phi$ [GeV cm$^{-2}$ s$^{-1}$ sr$^{-1}$]',fontsize=12)
    axs[i].set_title(subtitles[i])
    second_axs.append(axs[i].secondary_yaxis('right', functions=(lambda x: 3. * x, lambda x: x / 3.)))
    second_axs[i].set_ylabel(r"All Flavor $E^2\Phi$ [GeV cm$^{-2}$ s$^{-1}$ sr$^{-1}]$", fontsize=12)
    axs[i].grid(True, which='minor', alpha=0.1)
    axs[i].grid(True, which='major', alpha=0.4)
    #axs[i].legend(loc='lower left', fontsize=12)
    if DIFFUSE:
        axs[i].set_ylim(0.2e-10, 1e-5)
        axs[i].set_xlim(1.2e13 * units.eV / plotUnitsEnergy, 1e21 * units.eV / plotUnitsEnergy)
        #axs[i].set_yticks([1e-12, 1e-11,1e-10,1e-9,1e-8,1e-7, 1e-6, 1e-5])
        axs[i].set_yticks([1e-10,1e-9,1e-8,1e-7, 1e-6, 1e-5])
        axs[i].yaxis.set_minor_locator(tck.AutoMinorLocator())
        axs[i].minorticks_on()
        #axs[i].tick_params(axis='y', which='minor', left=True)
    else:
        axs[i].set_ylim(1e-11, 2e-6)
        axs[i].set_xlim(1e5, 1e11)
    
        
axs[i].set_xlabel(f'Neutrino Energy [{plotUnitsEnergyStr}]', fontsize=12)

fig.suptitle("Diffuse Flux, 1:1:1 Flavor Ratio  ", fontsize=14)
fig.subplots_adjust(top=0.94, hspace=0.18, bottom=0.05, right=0.92, left=0.08)
#fig.tight_layout()
#labels = []
#labels = add_limit(ax, labels, veff[:, 0], veff[:, 1], n_stations=100, livetime=5 * units.year, label=veff_label)
#labels = add_limit(ax, labels, veff[:, 0], veff[:, 1], n_stations=1000, livetime=5 * units.year, label=veff_label)
#plt.legend(handles=labels, loc=2)
if DIFFUSE:
    name_plot = "Limit_diffuse.pdf"
else:
    name_plot = "Limit_sources.pdf"
plt.savefig(name_plot)
plt.show()