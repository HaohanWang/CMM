![CMM](cmm.PNG "Coupled Mixed Model")

# CMM (Coupled Mixed Model)

Implementation of CMM in this paper:

Wang H., Pei F., Vanyukov MM., Bahar I., Wu W. and Xing EP, Coupled Mixed Model for Joint Genetic Analysis of Complex Disorders with Two Independently Collected Data Sets (under review)

## Introduction

CMM (Coupled Mixed Model) is used to simulnateously conduct genetic analysis for independently collected data sets of different related phenotypes.

## File Structure:

* models/ main method for CMM
* utility/ other helper files
* cmm.py main entry point of using CMM to work with your own data

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

## Python Users
Proficient python users can directly call the CMM method with python code, see example at [Line 261 to Line 266](https://github.com/HaohanWang/CMM/blob/master/cmm.py#L261)

## Contact
[Haohan Wang](http://www.cs.cmu.edu/~haohanw/)
&middot;
[@HaohanWang](https://twitter.com/HaohanWang)
