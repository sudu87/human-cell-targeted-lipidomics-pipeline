Datasets should be processed and analyzed using the R statistical environment (R Foundation for Statistical Computing). 
Raw lipid abundance tables generated from mass spectrometry–based lipidomics should be imported into R and combined with sample metadata describing experimental conditions, infection status, and inhibitor treatments. 
Lipid species will be retained only if they are detected above background levels in the majority of samples within at least one experimental group. 
Lipids with excessive missing values across all conditions will be excluded prior to downstream analyses.
Missing values arising from low-abundance signals will be handled conservatively and will not be imputed unless explicitly required for a specific statistical model. 
All filtering steps will be applied uniformly across experimental groups to avoid bias.
