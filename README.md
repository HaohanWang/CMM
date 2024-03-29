![CMM](cmm.PNG "Coupled Mixed Model")

# CMM (Coupled Mixed Model)

Implementation of CMM in this paper:

Wang, Haohan, Fen Pei, Michael M. Vanyukov, Ivet Bahar, Wei Wu, and Eric P. Xing. "Coupled mixed model for joint genetic analysis of complex disorders with two independently collected data sets." BMC bioinformatics 22, no. 1 (2021): 1-14. (_[link](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-03959-2)_)

## Introduction

CMM (Coupled Mixed Model) is used to simulnateously conduct genetic analysis for independently collected data sets of different related phenotypes. CMM aims to achieve this goal by inferring all the information enclosed by dashed lines in the following figure. 

![Introduction](intro.PNG "Introduction")

**Replication:** This repository serves for the purpose to guide others to use our tool, if you are interested in the scripts to replicate our results, please contact us and we will share the repository for replication. Contact information is at the bottom of this page.

## File Structure:

* [models/](https://github.com/HaohanWang/CMM/tree/master/model) main method for CMM
* [utility/](https://github.com/HaohanWang/CMM/tree/master/utility) other helper files
* [cmm.py](https://github.com/HaohanWang/CMM/blob/master/cmm.py) main entry point of using CMM to work with your own data

## An Example Command:

```
python cmm.py --file1 ./data/mice1.plink --file2 ./data/mice2.plink -m --snum 20
```
#### Instructions
```
  Options:
  -h, --help          show this help message and exit

  Data Options:
    --file1 INPUT FILE ONE       name of the first input file
    --file2 INPUT FILE TWO       name of the second input file

  Model Options:
    --lambda=LMBD     the weight of the penalizer. If neither lambda or snum
                      is given, cross validation will be run.
    --snum=SNUM       the number of targeted variables the model selects. If
                      neither lambda or snum is given, cross validation will
                      be run.
    -s                Stability selection
    -q                Run in quiet mode
    -m                Run without missing genotype imputation
```

#### Data Support
* CMM currently supports CSV and binary PLINK files.
* Extensions to other data format can be easily implemented through `FileReader` in `utility/dataLoadear`. Feel free to contact us for the support of other data format.

## Python Users
Proficient python users can directly call the CMM method with python code, see example at [Line 261 to Line 266](https://github.com/HaohanWang/CMM/blob/master/cmm.py#L261)

## Installation (Not Required)
You will need to have numpy, scipy and pysnptool installed on your current system.
You can install CMM using pip by doing the following

```
   pip install git+https://github.com/HaohanWang/CMM
```

You can also clone the repository and do a manual install.
```
   git clone https://github.com/HaohanWang/CMM
   python setup.py install
```

## Contact
[Haohan Wang](http://www.cs.cmu.edu/~haohanw/)
&middot;
[@HaohanWang](https://twitter.com/HaohanWang)
