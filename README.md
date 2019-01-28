# CMM (Coupled Mixed Model)

Implementation of CMM in this paper:

Wang H., Liu M., Pei F., Vanyukov MM., Bahar I., Wu W. and Xing EP, Coupled Mixed Model for Joint GWAS of Pairs of Complex Disorders with Independently Collected Datasets

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
## Software with GUI
Software with GUI will be avaliable through [GenAMap](http://genamap.org/)


## Contact
[Haohan Wang](http://www.cs.cmu.edu/~haohanw/)
