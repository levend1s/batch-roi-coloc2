# Usage

## draw_rois.ijm 
Open Fiji and run macro script draw_rois.ijm. The script will prompt you to select two directories, one after the other. 
- For the first one select the directory containing all your images. 
- For the second, select an empty directory that the drawn ROIs will be saved to.

The script will then loop through each image in the directory and wait for you to draw your region of interest. After clicking OK, it will move to the next image. It repeats this for all images in the directory.

## batch_coloc.ijm
From fiji run macro script batch_coloc.ijm. Coloc2 parameters can be changed in this script. It will ask for three directories.
- First, select directory containing all images.
- Second, select directory containing saved ROIs.
- Third select output directory to write colocalisation CSV file to.

The script will then loop through all ROIs and run Coloc2 with the script defined parameters.

## analyze_coloc.R

Run R script to generate plots. You will need to change variables to point at your `sample_metadata.csv` files and `coloc_results.csv`. You will need to create a `sample_medatada.csv` specific for your project. It tells the R script how to plot the data, and should look like this.

```
filename, condition, replicate
Image 2.czi,control, rep1
Image 30.czi,knocksideways, rep1
```