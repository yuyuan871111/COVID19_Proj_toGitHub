# 0: Prerequirement 
## (1) Remote connect information
* Taiwania 1: `ssh uxxxxxxx@140.110.148.11`  
* Taiwania 2: `ssh uxxxxxxx@ln01.twcc.ai`, gromacs version 2018.6  

## (2) Parameters settings 
Note: You can change parameters here.  
`pdbname="6m0j_omicron"`  
`INPUT_FILEPATH="../${pdbname}.pdb"`  

## (3) Make working directory
>`mkdir box`  
`mkdir g`  
`mkdir sol`  
`mkdir ion`  
`mkdir em`  
`mkdir pr`  
`mkdir md`  
`mkdir index`  
`mkdir mdp`  
`mkdir analysis`  

## (4) Upload your .mdp files
Via FileZilla to upload the following files.  
* `ions.mdp`
* `minim.mdp`
* `md.mdp`
* `npt.mdp` 
* `nvt.mdp`

 ## (5) Upload your .sh file
 Via FileZilla to upload the following files to specific folder. Also, make sure the path in the .sh file matched correctly with remote path.   
 * `em.sh` ***->*** /your/path/`em/em.sh` (remote)
 * `pr.sh` ***->*** /your/path/`pr/pr.sh` (remote)
 * `md.sh` ***->*** /your/path/`md/md.sh` (remote)

# 1: make a box
```
cd box
gmx editconf -f [inputfile.pdb] -o [outputfile.pdb] -bt [cubic] -box [6 6 6] -c
```
or 
```
gmx editconf -f [inputfile.pdb] -o [outputfile.pdb] -bt [cubic] -d [1.0] -c
```
> `cd box`  

> `gmx editconf -f $INPUT_FILEPATH -o ${pdbname}_box.pdb -bt cubic -d 1.0 -c `  
Note: Zn ion is not allowed in pdb for molecular simulation.  

# 2: g
```
cd g
gmx pdb2gmx -f [inputfile.pdb] -o [outputfile.pdb] -p [forcefield.top] -ff select -ignh -inter 
```
> `cd g`  

> `gmx pdb2gmx -f ../box/${pdbname}_box.pdb -o ${pdbname}_g.pdb -p LM.top -ff select -ignh -inter`   
Note: The ***Gromos45a3 force field*** and the ***spc water model*** are used.  
Note: Amino acid chrage were selected by users depending on pH.  
> e.g. pH 7.4 (LYSH:+1, ARG:+1, GLN:0, ASP:-1, GLU:-1, HISA:ND1 only, Cys linkage: yes for all, terminus: NH3+, COO-) 

Recording:  
* 6m0j (2500.39 nm^3, 142 mM Na+ = 214 Na+/box)  
  -26 (WT), -26 (alpha), -25 (beta), -25 (delta+), -24 (delta), -25 (gamma), -25 (epsilon), -26 (lambda), -24 (mu), -23 (omicron, BA.1, BA.2, BA.2.12.1, BA.4/5)
* 7a91 (3527.20 nm^3, 142 mM Na+ = 302 Na+/box)  
  -20 (WT), -20 (alpha), -19 (beta), -19 (delta+), -18 (delta), -19 (gamma), -19 (epsilon), -20 (lambda), -18 (mu), -16 (omicron, BA.1, different from 6m0j: T547K), -17 (omicron, BA.2, BA.2.12.1, BA.4/5)
* 7mjn (2307.58 nm^3, 142 mM Na+ = 197 Na+/box)  
-23 (alpha)
* 7v80 (2517.75 nm^3, 142 mM Na+ = 215 Na+/box)  
-22 (beta)
* 7v84 (2473.02 nm^3, 142 mM Na+ = 212 Na+/box)  
-22 (gamma)
* 7v8b (2568.76 nm^3, 142 mM Na+ = 220 Na+/box)  
-21 (delta)

# 3: sol
```
cd sol 
gmx solvate -cp [input.pdb] -o [outputfile.pdb] -p [forcefield.top] -cs 
```
## 3.1: add sol
> `cd sol`  

> `gmx solvate -cp ../g/${pdbname}_g.pdb -o ${pdbname}_sol.pdb -p ../g/LM.top -cs`

## 3.2: make index
```
cd index
gmx make_ndx -f [inputfile.pdb] -o [index.ndx] 
```
> `cd index`  

> `gmx make_ndx -f ../sol/${pdbname}_sol.pdb -o index.ndx`  
Note: select chain. e.g. `chainA` -> chA; `chainE` -> chE; `q` -> to save and quit  

# 4: ion [need ions.mdp]
## 4.1 add ion
```
cd ion 
gmx grompp -c [inputfile.pdb] -n [index.ndx] -f [parameters.mdp] -o [ion.tpr] -p [forcefield.top] -maxwarn [e.g. 5]
gmx genion -s [ion.tpr] -o [outputfile.pdb] -nname CL -nn [e.g. 20] -pname NA -np [e.g. 3]
```
Note: Select ***12*** to replace sol(water) to ions.  
> `cd ion`  

> `gmx grompp -c ../sol/${pdbname}_sol.pdb -f ../../../mdp/ions.mdp -o ion.tpr -p ../g/LM.top -maxwarn 5`  

> `gmx genion -s ion.tpr -o ${pdbname}_ion.pdb -pname NA -np 220 -nname CL -neutral`  
(6m0j: NA -np ***214***; 7a91: NA -np ***302***; 7mjn: NA -np ***197***; 7v80: NA -np ***215***; 7v80: NA -np ***212***; 7v8b: NA -np ***220***)  

## 4.2: revise forcefield file top
> `vi ../g/LM.top`  

Add the number of ions and remove exact number of water. e.g. 
```
vi [forcefield.top]
Then, to the bottom (hotkey: G) of .top file: 
...
SOL 123           ->      SOL 120    
                  ->      NA 1
                  ->      CL 2
```

# 5: em [need em.mdp]
```
cd em 
gmx grompp -c [inputfile.pdb] -f [parameters.mdp] -o [em.tpr] -p [forcefield.top] -n [../index/index.ndx]
gmx mdrun -s [em.tpr] -o [em.trr] -c [outputfile.pdb] -g [em.log] -e [em.edr]
```
> `gmx grompp -c ../ion/${pdbname}_ion.pdb -f ../../../mdp/minim.mdp -o em.tpr -p ../g/LM.top -n ../index/index.ndx`  

> ***Taiwania1***: `qsub em.sh` or ***Taiwania2***: `sbatch em.sh`  

# 6: pr [need pr.mdp]
```
cd pr
gmx grompp -c [inputfile.pdb] -r [inputfile.pdb] -f [parameters.mdp] -o [pr.tpr] -p [forcefield.top] (-t [prev.cpt])
gmx mdrun -s [pr.tpr] -o [pr.trr] -c [outputfile.pdb] -g [pr.log] -e [pr.edr] -cpo [pr.cpt] 
```
## 6.1: nvt 
### 6.1.1: run
> e.g.  
`gmx grompp -c ../em/${pdbname}_em.pdb -r ../em/${pdbname}_em.pdb -f ../../../mdp/nvt.mdp -o pr.tpr -p ../g/LM.top`  

> ***Taiwania1***: `qsub pr.sh` or ***Taiwania2***: `sbatch pr.sh`  

### 6.1.2: move nvt results to nvt folder
> `mkdir nvt`  

> `mv p* nvt/.` \#(all pr related files)  

> `mv *.pdb nvt/.` \#(pdb file)  

> `mv m* nvt/.` \#(mdp file)  

> `cp nvt/pr.sh .`

## 6.2: npt 
### 6.2.1: run
> `gmx grompp -c nvt/${pdbname}_pr.pdb -r nvt/${pdbname}_pr.pdb -f ../../../mdp/npt.mdp -o pr.tpr -p ../g/LM.top -t nvt/pr.cpt`

> ***Taiwania1***: `qsub pr.sh` or ***Taiwania2***: `sbatch pr.sh`  
Note when running:  
NOTE 1 [file LM.top, line 51]:
  You are combining position restraints with Parrinello-Rahman pressure
  coupling, which can lead to instabilities. If you really want to combine
  position restraints with pressure coupling, we suggest to use Berendsen
  pressure coupling instead.  

### 6.2.2: move npt results to npt folder
> `mkdir npt`  

> `mv p* npt/.` \#(all pr related files)  

> `mv *.pdb npt/.` \#(pdb file)  

> `mv m* npt/.` \#(mdp file)  

# 7: md [need md.mdp]
## 7.1: run
```
cd md
gmx grompp -c [inputfile.pdb] -f [parameters.mdp] -o [md.tpr] -p [forcefield.top] -t [prev.cpt]
gmx mdrun -s [md.tpr] -o [md.trr] -c [outputfile.pdb] -g [md.log] -e [md.edr] -cpo [md.cpt]
```
> `cd md`

> `gmx grompp -c ../pr/npt/${pdbname}_pr.pdb -f ../../../mdp/md.mdp -o md.tpr -p ../g/LM.top -t ../pr/npt/pr.cpt `  

> ***Taiwania1***: `qsub md.sh` or ***Taiwania2***: `sbatch md.sh`  

## 7.2: md continue [need .edr, .log, .cpt, .tpr, traj_comp.xtc]
```
gmx convert-tpr -s [previous.tpr] -o [continue.tpr] -until [time(ps), e.g. 10000]
gmx mdrun -s [continue.tpr] -o [md.trr] -c [outputfile.pdb] -g [md.log] -e [md.edr] -cpo [md.cpt] -cpi [md_prev.cpt]
```
> `gmx convert-tpr -s md_prev.tpr -o md.tpr -until 10000`  

> ***Taiwania1***: `qsub md_conti.sh` or ***Taiwania2***: `sbatch md_conti.sh`

Note: You need to prepare md_conti.sh if you want to continue your run.

# 8: analysis
You can change your results path.  
`Results_PATH="/home/u1397281/Project/COVID19_proj/gromacs/hACE2_S1RBD_7a91/WT/analysis"`

## 8.1: potential/temperature/pressure/density equilibrium
`echo 10 0 | gmx energy -f ../em/em.edr -o ${Results_PATH}/potential_em.xvg` \#(10 0)  
`echo 16 0 | gmx energy -f ../pr/nvt/pr.edr -o ${Results_PATH}/temperature_nvt.xvg` \#(16 0)  
`echo 18 0 | gmx energy -f ../pr/npt/pr.edr -o ${Results_PATH}/pressure_npt.xvg` \#(18 0)  
`echo 24 0 | gmx energy -f ../pr/npt/pr.edr -o ${Results_PATH}/density_npt.xvg` \#(24 0)

## 8.2: calibrate .pdb file
### 8.2.1: if multiple chain in .pdb, chain information may miss
* METHOD 1: use text editor to pdb file  
-> add `TER` and `HEADER [name]` in the pdb file -> open pdb by pymol.  

* METHOD 2: use gromacs .xtc to calibrate (same as 8.2.2)  
  `echo 1 1 | gmx trjconv -s ../md/md.tpr -f ../md/${pdbname}_md.pdb -o ${Results_PATH}/${pdbname}_md_noPBC.pdb -pbc nojump -center` #(1 1 : protein protein)

### 8.2.2: structure recalibration (when PBC) 
* METHOD 1: gromacs  
SELCET: (centering group)  
`echo 1 1 | gmx trjconv -s ../md/md.tpr -f ../md/${pdbname}_md.pdb -o ${Results_PATH}/${pdbname}_md_noPBC.pdb -pbc nojump -center` #(1 1 : protein protein)  

* METHOD 2: text editor + pymol  
-> (text editor) pdb file modify: HEADER ..., TER  
-> (pymol)  
`create "S1RBD_WT_copy", "S1RBD_WT_100ns"`  
`alter /S1RBD_WT_copy//"", chain="B"`  
`alter /hACE2_WT//"", chain="A"`  
`reset`  
`translate [0, 0, 136.527], "S1RBD_WT_copy"`  
`export "sele"`

## 8.3: create index to separate hACE2 and S1RBD
`cd index`  
`gmx make_ndx -f ../md/${pdbname}_md.pdb -o index_md.ndx `   
> \> `ri xxx-xxx`  
> 6m0j => 19: `ri 1-597` (hACE2), 20: `ri 598-791` (S1RBD)  
> 7a91 => 19: `ri 1-588` (hACE2), 20: `ri 589-821` (S1RBD)  
> 7mjn, 7v80, 7v84, 7v8b => 19: `ri 201-796` (hACE2), 20: `ri 1-200` (S1RBD)  

### 8.3.1: create index for mutual region     
> \> `ri xxx-xxx | ri xxx-xxx | ri xxx-xxx & 4` (4: backbone)  
> 6m0j => 21: `ri 1-115 | ri 123-595 | ri 598-791 & 4`  
> 7a91 => 21: `ri 1-782 & 4`   
> 7mjn, 7v80, 7v84, 7v8b => 21: `ri 201-315 | ri 323-795 | ri 3-196 & 4`  

## 8.4: traj recalibration (correction PBC)
`echo 1 0 | gmx trjconv -s ../md/md.tpr -f ../md/traj_comp.xtc -o md_noPBC.xtc -pbc nojump -center` \#(1 0, better)  
or  
`gmx trjconv -s md.tpr -f traj_comp.xtc -o md_noPBC.xtc -pbc mol -center` \#(1 0)

## 8.5: rmsd
first: fitting backbone; second: analysis group  
`echo 4 4 | gmx rms -s md.tpr -f md_noPBC.xtc -o rmsd.xvg -tu ns`  
or  
`gmx rms -s md.tpr -f md_noPBC.xtc -o rmsd.xvg -tu ns` \#(4 4 backbone)  
or  
(relative to crystall structure): `gmx rms -s ../em/em.tpr -f md_noPBC.xtc -o rmsd_xtal.xvg -tu ns`  
### 8.5.1: rmsd in mutual region
`echo 4 21 | gmx rms -s md.tpr -f md_noPBC.xtc -o rmsd_mutual.xvg -tu ns -n ../index/index_md.ndx` \#(21 mutual region)   

## 8.6: rmsf
`echo 19 | gmx rmsf -s md.tpr -f md_noPBC.xtc -o rmsf_hACE2.xvg -res -n ../index/index_md.ndx` \#(19 hACE2 ...)  
`echo 20 | gmx rmsf -s md.tpr -f md_noPBC.xtc -o rmsf_S1RBD.xvg -res -n ../index/index_md.ndx` \#(20 S1RBD ...)

## 8.7: radius of gyration
`echo 1 | gmx gyrate -s md.tpr -f md_noPBC.xtc -o gyrate.xvg` \#(1 protein)
### 8.7.1: radius of gyration in mutual region
`echo 21 | gmx gyrate -s md.tpr -f md_noPBC.xtc -o gyrate_mutual.xvg -n ../index/index_md.ndx` \#(21 mutual region)  

## 8.8: h-bond (distance distribution, number, angle distribution and auto-correlation)
`echo 19 20 | gmx hbond -f md_noPBC.xtc -s md.tpr -tu ns -n ../index/index_md.ndx -temp 310.15 -num -ac -dist -ang` 

## 8.9: pdb traj dynamic movie
### 8.9.1: (METHOD 1) compressed traj + vmd
1. create correct pdb (pr-npt)  
`echo 1 1 | gmx trjconv -s pr.tpr -f ${pdbname}_pr_npt.pdb -o ${pdbname}_pr_npt_noPBC.pdb -pbc nojump -center` \#(1 1)  
2. create compressed traj  
`echo 1 | gmx trjconv -s md.tpr -f md_noPBC.xtc -o md_protein.xtc -dt 100` \#(1: protein)  
3. import traj and pdb(pr-npt) to vmd
4. fit different variants (vmd)  
`set sel0 [atomselect 29 "chain A"]`  
`set sel1 [atomselect 32 "chain A"]`  
`set M [measure fit $sel1 $sel0]`  
`set sel1 [atomselect 1 protein]`  
`$sel1 move $M`  
Note: choose correct variant type.

### 8.9.2: (METHOD 2) pymol
1. create correct pdb (pr-npt, same as previous 8.9.1.1) 
2. create compressed traj (same as previous 8.9.1.2)  
3. load traj and pdb into pymol (pymol loading heavy, QQ)  
`load ${pdbname}_pr_npt.pdb, ${pdbname}`  
`load_traj md_protein.xtc, ${pdbname}`  
`smooth protein, 30, 3`  


## 8.10: dssp
```
gmx do_dssp -f .xtc -s .tpr -n .ndx -o ss.xpm -sc .xvg -tu ns
gmx xpm2ps -f ss.xpm -o ss.eps -bx 1 -by 30
```
> `echo 19 | gmx do_dssp -f .xtc -s md.tpr -n ../index/index_md.ndx -o dssp_S1RBD.xpm -sc -tu ns`  
`gmx xpm2ps -f dssp_S1RBD.xpm -o dssp_S1RBD.eps -bx 1 -by 30`  
or  
`echo 20 | gmx do_dssp -f .xtc -s md.tpr -n ../index/index_md.ndx -o dssp_hACE2.xpm -sc -tu ns`  
`gmx xpm2ps -f dssp_hACE2.xpm -o dssp_hACE2.eps -bx 1 -by 30`

## 8.11: PCA
`echo 4 4 | gmx covar -s md.tpr -f md_noPBC.xtc -o -v -l -av ${pdbname}_avg.pdb`  \#( 4 4 : backbone backbone )  
`echo 4 4 | gmx anaeig -v eigenvec.trr -f md_noPBC.xtc -s md.tpr -2d 2dproj.xvg -first 1 -last 2` \#( 4 4 : backbone backbone )

## 8.12: MM/PBSA 
Reference: [g_mmpbsa package](https://rashmikumari.github.io/g_mmpbsa/)  
1. Follow the installation instructions in above reference.
2. Prepare the parameters file (`mmpbsa.mdp`).
3. Run:  
```
echo 19 20 | g_mmpbsa -f md_noPBC.xtc        -s md.tpr \
                      -i mmpbsa.mdp          -n index.ndx -pbsa \
                      -mm energy_MM.xvg      -pol polar.xvg \
                      -apol apolar.xvg       -decomp \
                      -mmcon contrib_MM.dat  -pcon contrib_pol.dat \
                      -apcon contrib_apol.dat 
```

# 9: Others
## 9.1: check space usage (Taiwania 1)
`lfs quota -u uxxxxxxx /home`

## 9.2: pymol - alter chain, segi  
pdb structure: `/name/segi/chain/residues`  
(pymol): `alter /name/segi, segi=""`  
(pymol): `alter /name/segi, chain="E"`  