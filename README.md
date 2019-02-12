# CONPER (continous performance task) project

## Data

[download](https://surfdrive.surf.nl/files/index.php/s/Um8zXm4tIJgGKB6)

### Functional (`fDATA_7T_3DMC_HP_CSF_WD.mat`) 

**Matrices**

- TS: BOLD signal (10 subjects by 4 states by 85 brain regions by 600 timepoints)
- PHI: instantaneous phases (10 subjects by 4 states by 85 brain regions by 600 timepoints)
- R: instantaneous coherence (10 subjects by 4 states by 85 brain regions by 600 timepoints)

**States**

1. Rest
2. Task 1: mental rotation
3. Task 2: orientation discrimination
4. Task 3: numerosity

### Structural (`sDATA_HCP_Q3_SLT_GS.mat`)

- CON: structural connectivity (85 by 85 brain regions)

### Regions of interest (`ROIS.mat`)

- cortex_id: indicates all cortical regions (68 of 85)
- ROIDescriptions: names of all 85 brain regions

Other files were used during pre-processing and are not relevant.

## Tasks

All subjects underwent a [staircase procedure](https://en.wikipedia.org/wiki/Psychophysics#Staircase_procedures) to individually fix the task difficulty
at 80 percent correct prior to undergoing the scan. For mental rotation this procedure was not successful. Subjects continued to improve during the scan and achieved levels of around 95 percent correct.

### Mental rotation

Indicate whether a presented letter is mirrored or normal. Letters are rotated, requiring subjects to mentally rotate them.

### Orientation discrimination

Indicate whether a [Gabor patch](http://neuroanatody.com/2016/05/whats-in-a-gabor-patch/) is rotated clockwise or counter-clockwise with respect to a reference orientation.

### Numerosity

Indicate whether [a cloud of dots](https://www.google.com/search?q=numerosity&client=ubuntu&hs=xbW&channel=fs&source=lnms&tbm=isch&sa=X&ved=0ahUKEwjw6NHC5bXgAhVBJ1AKHUOzAQEQ_AUIDigB&biw=1855&bih=959#imgrc=F0b6xHnXiaQc8M:) contains more or less than 15 dots without counting.