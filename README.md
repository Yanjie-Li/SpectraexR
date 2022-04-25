# SpectraexR
this is a shiny app for single tree multi-spectral extraction
This is a shiny app specially used to extract spectral information of a single tree, and I also provide a sample to show what the data looks like after extraction.

First of all, you need point cloud data with precise positioning information, which is used to segment each individual plant in a large area of forest land. In addition, you need to have multi-spectral or hyperspectral information of this forest land, such as red, green, blue, red tif data in edge and near-infrared bands. After you upload these two parts of data, the app will run automatically.

To run the app, use "shiny::runGitHub( "SpectraexR", "Yanjie-Li")" to download and run.
