# ResidualPDF
Plot residuals of ground-motion models (GMMs) to visualize GMM performance with period

## Setting up your conda environment: 
```bash
conda env create -f env_respdf.yml
```

Then activate with
```bash
conda activate respdf
```
(the name of your environment is set by the first line ```name: ``` so you can edit this as you wish).  Remember, you will need to run the ```conda activate``` command every time you open a new terminal or tab!

## Residual PDF plots for individual ground motion models
To make a residual PDF for the A10 ground motion model, with residuals stored in file A10_PDF.out and sigmas in A10_sigma.out in current working directory:
```bash
python respdf.py A10
```
this will create a file called A10pdf.png

You can quickly regenerate plots for all of the GMMs in your directory by running the script ```remakeplots.sh``` which will loop over all *_PDF.out files in the directory.

## Summary plots showing means and sigmas for many models
If you have lots of *_PDF.out and *_sigma.out files in your current working directory, you can produce figures summarizing the mean residuals with:
```bash
python resmeans.py 
```
This will create a heat map of the GMM means (resmeans_heat.png) and a plot of residuals and sigmas vs. period (resmeans_lines.png)
