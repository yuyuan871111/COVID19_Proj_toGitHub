# COVID Proj: From Delta to Omicron, Evaluating the Interactions of Multiple Spike Proteins RBD and hACE2 through Systematic Competitive Analyses and Molecular Dynamics

## Analysis  
All analysis are in Jupyter notebook files.  
Root path: `./Analysis/`  
`Gromacs_RMSF.ipynb`: RMSF analysis file.  
`Gromacs_data_wTIME.ipynb`: root mean squre deviation (RMSD)/Radius of gyration (Rg)/Number of hydrogen bonds and pairs(HB-pairs)/correlation between hydrogen bonds and pairs (corr-HB) over time.  
`Gromacs_distribution.ipynb`: distribution of distance and angle of hydrogen bonds.  
`Gromacs_equilibrium.ipynb`: energy/temperature/pressure/density in equilibrium process (energy minimization/NVT/NPT).  
`Gromacs_features_clustering.ipynb`: using molecular dynamics features for clustering into different variants.  
`Gromacs_hbac_fitting.ipynb`: analaysis of autocorrelation of hydrogen bonds (fitting, calculating relaxation time of hydrogen bonds, and so on).  
`VariantsFreqVaccPlot.ipynb`: integrate multiple sources public health data (proportion of vaccinated individuals/variant frequency) into single plot.  

## Data
Data for analysis (from public health resources or molecular dynamics).  
Root path: `./Data/`  

### gromacs parameter files
Path: `./Data/gromacs_mdp`  
`ions.mdp`: parameter file when adding ions.  
`md.mdp`: parameter file when performing molecular dynamics for 100 ns.  
`minim.mdp`: parameter file when performing energy minimization (EM).  
`npt.mdp`: parameter file when perfomring equilibrium in normal pressure temperature (NPT).  
`nvt.mdp`: parameter file when perfomring equilibrium in normal volume temperature (NVT).  

### metadata
Path:  `./Data/metadata`  
`CoVariants_data`: public health data from CoVariants (about variant frequency).  
`Owid_data`: public health data from Our World in Data (about proportion of vaccinated individuals).  
`md_tidydata`: tidy data of molecular dynamics data.  
* `ComparisonOfVariants.xlsx`: All tidy data here.  
* `ComparisonOfVariants_numeric.csv`: Numeric tidy data for further analysis (e.g. clustering).  
* `ComparisonOfVariants_wMutual_numeric.csv`: Numeric tidy data for further analysis (e.g. clustering). Mutual residues among different PDB model were selected in analysis of RMSD/Rg so that the deviation from different PDB can be reduced.  
* `RMSF.csv`: RMSF data for following analysis.  

### raw data
Path: `./Data/rawdata`  

#### Molecular dyanmics (MD) results
`variants_6m0j`: WT/alpha/beta/gamma/delta/delta plus/epsilon/gamma/lambda/mu/omicron.  
`variants_7a91`: WT/alpha/beta/gamma/delta/delta plus/epsilon/gamma/lambda/mu/omicron.  
`variants_7mjn`: alpha.  
`variants_7v80`: beta.  
`variants_7v84`: gamma.  
`variants_7v8b`: delta.  

***Each folder contains the following items:***  
* `6m0j_WT_afterMD_dimplot`: Ligplot+ dimplot with chord plot by python3 script [here](https://github.com/yuyuan871111/dimpyplot_chordplot) (with after MD structure).  
* `6m0j_WT_beforeMD_dimplot`: Ligplot+ dimplot with chord plot by python3 script [here](https://github.com/yuyuan871111/dimpyplot_chordplot) (with before MD structure).  
* `2dproj.xvg`: PCA analysis for trajectory.  
* `6m0j_WT_afterMD.pdb`: initial structure before MD.  
* `6m0j_WT_beforeMD.pdb`: final structure after MD.  
* `density_npt.xvg`: denesity over time in NPT.  
* `gyrate.xvg`: Rg over time.  
* `gyrate_mutual.xvg`: Rg over time with mutual residues among 6m0j, 7a91, 7mjn, 7v80, 7v84, 7v8b.  
* `hbac.xvg`: autocorrelation of hydrogen bonds between hACE2 and S1RBD.  
* `hbang.xvg`: distribution of the angle of hydrogen bonds.  
* `hbdist.xvg`: distribution of the distance of hydrogen bonds.  
* `hbnum.xvg`: number of hydrogen bonds and pairs over time.  
* `md_protein.xtc`: compressed trajectory of 100 ns MD.  
* `potential_em.xvg`: energy over time in EM.  
* `pressure_npt.xvg`: pressure over time in NPT.  
* `rmsd.xvg`: RMSD over time.  
* `rmsd_mutual.xvg`: RMSD over time with mutual residues among 6m0j, 7a91, 7mjn, 7v80, 7v84, 7v8b.  
* `rmsf_S1RBD.xvg`: RMSF in S1RBD regions.  
* `rmsf_hACE2.xvg`: RMSF in hACE2 regions.  
* `temperature_nvt.xvg`: temperature over time in NVT.  

#### pdb input
All pdb input for whole molecular dynamics procedure were backuped in `./Data/rawdata/pdb_input`.  

## Structural Analysis Tools 
A python 3 self-defined packaged for this study, which might be used in analysis.  
Root path: `./StructuralAnalysisTools`

## Visualization
Figures were generated from analysis including graphical abstract, structrual view, MD results over time, fitting curves, ...  
Root path: `./Visualization/`  
### illustration
* `./illu/COVID_variants_v*.png`: mutations information in different PDBs (* is the version number).  
* `./illu/MyCOVID_Project_Workflow_v*.png`: graphical abstract of my analysis (* is the version number).  
### other plots





