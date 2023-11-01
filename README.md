# Dashboard for label-free quantification in proteomics


## Project description
This is a simple dashboard developed for exploring results of quantitative analysis using shotgun proteomics data. It allows you to select different quantification (NSAF, SIn) and imputation methods (average, minimum in replicate) and filter the results based on FDR and fold change thresholds. The output includes a volcano plot and a table of proteins with the corresponding p-values and fold change values.

---

## Pre-requisites 
It is recommended to use the dashboard in a <a href="https://conda.io/docs/user-guide/install/">conda</a> environment. To create and activate a new environment with all necessary Python packages (provided in "requirements.txt"), run the following commands:

```
conda create --name <env> --file requirements.txt
conda activate <env>
```

---

## Input data
As of now, the dashboard accepts the output provided by <a href="http://tools.proteomecenter.org/wiki/index.php?title=Software:StPeter">StPeter</a> in prot.xml format. Samples of acceptable input files are included in the folder "examples".

---

## Usage

The dashboard is meant for personal use, so you can easily run it on your local device in the debug mode with the following command:

```
python main_app.py -c interact-control_1.ipro.prot.xml -t interact-treatment_2.ipro.prot.xml
```

*Option -c is followed by the paths of files with control data, option -t is followed by the paths of files with treatment data.*

Once you run the command, open localhost:8050 in your browser.


To see the explanation of all command line arguments:

```
python main_app.py -h
```


