# SGLT2i-hyperfiltratiton-model
A Quantitative Systems Pharmacology kidney model of diabetes associated renal hyperfiltration and the effects of SGLT inhibitors.

![gim](https://github.com/PavelBal/SGLT2i-hyperfiltration-model/blob/master/Figure_1_color.png)
_Figure taken from [1], copyright by CPT:PSP._

Within this repository, we distribute a mechanistic model of renal hyperfiltration coupled with physiologically-based pharmacokinetics/pharmacodynamics (PBPK/PD) of sodium-glucose cotransporter 2 inhibitor (SGLT2i) dapagliflozin and integrated into the glucose-insulin-glucagon model (GIM) [[1](#references)]. The GIM (Panel A of the Figure), as available at [the repository](https://github.com/Open-Systems-Pharmacology/Glucose-Insulin-Model), has been extended by a detailed kidney reprersentation (Panel B of the Figure). A PBPK/PD model of dapagliflozin servers as a representative for the class of SGLT2i medication.

## Repository files
The model is provided as two MoBi projects.
**RenalModel_SGLT2i_Balazki_et_al.mbp3** is the PBPK/PD model of dapagliflozin with the extended kidney structure. The model was created with the Open Systems Pharmacology Suite (OSPS) version 7.2. The project includes simulations of dapagliflozin administrations: intravenous (iv) in rat, monkey, dog [[2](#references)], and human [[3](#references)], and various doses oral (po) administration in human [[2,4,5](#references)].

Additionally, a stepped hyperglycemic clamp (SHC) reported in [[6](#references)] is simulated (_Dapagliflozin_2013_SHC_healthy_) and compares simulated urinary glucose excretion with data and assesses the impact of plasma glucose on glomerular filtration rate (GFR). 

**QSP_Diab_RenalModel_Balazki_et_al.mbp3** is the PBPK/PD model of dapagliflozin and the extended kidney model integrated into GIM. The physiology is based on the PBPK model implemented in PK-Sim version 5.6. The project includes simulations of the multiple ascending dose (MAD) study reported in [[3](#references)] and compares observed urinary glucose excretion with simulated values.

The folder **RScript** contains R-script files used to generate figures of the original publication [[1](#references)]. To execute them, an installation of the [OSPS toolbox for R](https://github.com/Open-Systems-Pharmacology/R-Toolbox/releases) is required. The archives **SGLT2i_hyperfiltration_xml.zip** and **ExpData.zip** include the model files and experimentla data extracted from literature, respectively, required for executing the scripts. If you want to run the scripts, adjust the paths to the respective folders in the script files.

## Version information
The physiology of dapagliflozin model is based on the PBPK model implemented in PK-Sim version 7.2.
The physiology of GIM is based on the PBPK model implemented in PK-Sim version 5.6.
The MoBi project files were created in version 7.2.

## Code of conduct
Everyone interacting in the Open Systems Pharmacology community (codebases, issue trackers, chat rooms, mailing lists etc...) is expected to follow the Open Systems Pharmacology [code of conduct](https://gitprint.com/Open-Systems-Pharmacology/Suite/blob/master/CODE_OF_CONDUCT.md).

## License
The model code is distributed under the [GPLv2 License](https://github.com/Open-Systems-Pharmacology/Suite/blob/develop/LICENSE).

## References
[1] [Balazki P, Schaller S, Eissing T, Lehr T. A Quantitative Systems Pharmacology kidney model of diabetes associated renal hyperfiltration and the effects of SGLT inhibitors. CPT: Pharmacometrics & Systems Pharmacology (2018) XXX EDIT 2:e65; doi:10.1038/psp.2013.40.](http://onlinelibrary.wiley.com/doi/10.1038/psp.2013.40/abstract)

[2] [Obermeier M. T. et al. In vitro characterization and pharmacokinetics of dapagliflozin (BMS-512148), a potent sodium-glucose cotransporter type II inhibitor, in animals and humans. Drug Metab. Dispos. 38, 405–414 (2010)]

[3][Boulton D. W. et al. Simultaneous oral therapeutic and intravenous 14 C-microdoses to determine the absolute oral bioavailability of saxagliptin and dapagliflozin. Br. J. Clin. Pharmacol. 75, 763–768 (2013).]

[4][Komoroski B. et al. Dapagliflozin, a novel SGLT2 inhibitor, induces dose-dependent glucosuria in healthy subjects. Clin. Pharmacol. Ther. 85, 520–526 (2009)]

[5][Kasichayanula S. et al. Influence of hepatic impairment on the pharmacokinetics and safety profile of dapagliflozin: An open-label, parallel-group, single-dose study. Clin. Ther. 33, 1798–1808 (2011).]

[6][Defronzo, R. A. Characterization of renal glucose reabsorption in response to dapagliflozin in healthy subjects and subjects with type 2 diabetes. Diabetes Care 36, 3169–3176 (2013).]
