loadOrInstall <- function (packageName, type = "CRAN") {
    
    isPackageInstalled <- packageName %in% rownames(installed.packages())
    
    if (!isPackageInstalled) {
        
        if (type == "CRAN") {
            
            install.packages(packageName)
            
        } else if (type == "bioc") {
            
            BiocManager::install(packageName)
            
        }
        
    }
    
    library(packageName, character.only = TRUE)
}

cranPackages <- c(
    "shiny",
    "shinyjs",
    "shinythemes",
    "shinydashboard",
    "openxlsx",
    "DT",
    "weights",
    "BiocManager"
)

biocPackages <- c(
    "Biostrings",
    "BSgenome"
)

for (package in cranPackages) {
    
    loadOrInstall(package)
    
}

for (package in biocPackages) {
    
    loadOrInstall(package, type = "bioc")
    
}
