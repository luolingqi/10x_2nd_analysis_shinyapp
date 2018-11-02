# 10x_2nd_analysis_shinyapp
  		The repository contains all the require files to build a docker-based shiny app implementing the 
    Bioconductor workflow of the 2nd analysis on droplet-based scRNASeq count data (e.g. 10X, link: 
    https://bioconductor.org/packages/release/workflows/vignettes/simpleSingleCell/inst/doc/work-3-tenx.html#overview). 
		
    	A docker image can be built either locally or on a EC2 instance where Docker has been installed, by running 'cd Shinyapp; docker build -t shiny-app .'

    The suggested system requirement is 8Gb RAM, 8-16Gb HD. As described in the Dockerfile, what docker builder does include: 
		1) Import a docker image with a well-compiled R3.5.1 on Ubuntu (check Evernote for detail)
		2) Install additional R packages required for run a Shiny app (e.g. shiny, shinyjs, shinythemes) and 
	the Bioconductor workflow (e.g. DropletUtils, scater, scran, AnnotationHub, EnsDb.Hsapiens.v86, Rtsne, pheatmap, irlba, XML and more)
		3) Install Shiny Server, the version 1.5.9.923
		4) Copy over the data and conf file of the app and run server upon docker launch
		
	Docker can be run as "docker run -d -p 3838:3838 --name shinyapp_grch38 shiny-app"
  
