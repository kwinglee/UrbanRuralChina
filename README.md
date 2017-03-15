# UrbanRuralChina

This directory contains the code used to analyze the data for Winglee et al., Recent Urbanization in China Is Correlated With a Westernized Microbiome Encoding Increased Virulence and Antibiotic Resistance Genes.

For all code (all code here is written in R but saved as a text file), you will need to modify the working directory. This line is usually the second line of code, and is currently written: setwd(""). Put the directory name in these quotes. The code is currently set up to assume everything is in the same directory; if this is not your setup you may need to add file paths where relevant. In each section, if I refer to a file in a different section I give the path name. For example, 16SrRNA/inputData/RDP would indicate the RDP folder in the inputData folder in the 16SrRNA folder.

##Metadata Analysis##

This section refers to the files in the metadata folder and contains the analysis of the metadata, including diet data, fasting blood, and spot urine results.

The metadata table is metadata_11232015.txt. The values of select columns have been set to "removed" to protect subject privacy. metadataConversion.txt indicates what each variable is measuring.

Model for metadata: metadata ~ urbanRural

formatMetadata_2015-11.R.txt -> clean this data for downstream analyses (metadata_cleaned11-2015.txt) and to produce the merged files with metabolite and 16S OTU data for the classifiers. It takes as input metadata_11232015.txt, 16SrRNA/inputData/abundantOTU/abundantOTULogNormal.txt, and metabolite/MetabolonScaledImpData-transpose.txt.
metadataFig.R.txt -> draws a PCoA for the metadata and generates a table of orthonormal site scores (pcoaCorrected_metadata.txt). Its input is metadata_cleaned11-2015.txt from formatMetadata_2015-11.R.txt.
metadataModel.R.txt -> p-values and boxplots for each metadata (output names start with metadataModel). Its input is metadata_cleaned11-2015.txt from formatMetadata_2015-11.R.txt.

##Metabolite Analysis##

This section refers to the files in the metabolite folder and contains the analysis of the Metabolon results.

The data is in MetabolonScaledImpData-transpose.txt. metabolonPathwayInfo.txt contains pathway details for each metabolite.

Model used for metabolite data: variable ~ urbanRural

mergeData.R.txt -> merge metabolite data with the OTU normalized abundances for the classifiers. Its inputs are MetabolonScaledImpData-transpose.txt and 16SrRNA/inputData/abundantOTU/abundantOTULogNormal.txt.
metabolonModel.R.txt -> model p-values and boxplots for each metabolite. Its input is MetabolonScaledImpData-transpose.txt. Output names start with metabolonModel.
pcoaModel_metabolonMetadata.R.txt -> p-values for the PCoA axes for both the metabolite data and metadata The inputs come from figures/metabolonFig.R and metadata/metadataFig.R. Output names start with pcoaModel.

##16S rRNA Analysis##

This section refers to the files in the 16SrRNA folder and contains the 16S rRNA analysis.

Files in the inputData/RDP folder contain the RDP tables. RDPtaxonomy.txt contains the full taxonomic classification, drawn from the outputs of running RDP on several data sets. Files in the inputData/abundantOTU folder contain the AbundantOTU+ results. Files in the inputData/qiime folder contain the QIIME results. Files in the databaseResults folder contain the results of aligning the OTUs to various databases. Files whose name contains LogNorm or log-norm are log normalized; otherwise files whose name contains taxaAsColumns contain the counts.

Model used for 16S rRNA data: variable ~ time + urbanRural + 1|subjectID

The code for the analyses below is in the analysis folder. All of them use as input the files in inputData unless otherwise indicated.
otuModel_16s.R.txt -> calculates p-values using the mixed linear model for each taxa and produces boxplots of the results. Output names start with otuModel.
pcoaWithDiversity.R.txt -> draws individual PCoA plots with logged diversity measures and produces the table of weighted orthonormal site scores for each sample. Output table names start with pcoaCorrected.
pcoaModel_16s.R.txt -> calculates p-values using the mixed linear model for the PCoA axes and produces boxplots of the results. The input for this is the pcoaCorrected tables from pcoaWithDiversity.R.txt. Output names start with pcoaModel.
rarefyDiversity.R.txt -> rarefies the samples, using a subsample size equal to the sample with the fewest number of classified reads at that taxonomic level, and calculates the richness, Shannon diversity index, Inverse Simpson diversity index, and evenness. This is repeated 10 times, and the average is also calculated, along with p-values from the mixed linear model for the average diversity measures. Output file names start with rarefyDiversity.
relAbun.R.txt -> calculates the relative abundances. Output file names end with relAbunFwd.
formatRDPData.R.txt-> gets the RDP tables into the correct format for the classifiers. Outputs files with names ending in forSVMLight.
classifyPlusROCcforest.R.txt -> uses the cforest function (a random forest classifier) to classify the data and writes the predictions and draws ROC curves. The input data comes from the inputData, metadata and metabolite folders and the output of formatRDPData.R.txt and metadata/formatMetadata_2015-11.R.txt and metabolite/mergeData.R.txt. Outputs files with names that start with cforest.
classifyPlusROCrpart.R -> uses the rpart function to classify the data and writes the predictions and draws ROC curves. The input data comes from the inputData, metadata and metabolite folders and the output of formatRDPData.R.txt and metadata/formatMetadata_2015-11.R.txt and metabolite/mergeData.R.txt. Outputs files with names that start with rpart.
model_otu-metadata-metabolon.R.txt -> calculates associations between 16S rRNA taxonomic abundances, metadata and metabolite data and plots the results. Inputs are metadata_cleaned11-2015.txt from metadata/formatMetadata_2015-11.R.txt, metabolite/MetabolonScaledImpData-transpose.txt, the files in inputData/RDP and inputData/abundantOTU/abundantOTUForwardTaxaAsColumnsLogNormalWithMetadata.txt. Outputs files whose names start with model.

##Minikraken Analysis##

This section refers to the files in the minikraken folder and contains the WGS taxonomic analysis using Kraken with the minikraken database.

Files in inputData contain the raw counts from Kraken for each taxa for each sample.

Model for minikraken data: variable ~ urbanRural

The code for the analyses is in the analysis folder:
formatMinikrakenTable.R.txt -> add metadata (from the output of 16SrRNA/analysis/relAbun.R) and transpose and log normalize the tables (use the files in inputData as input). The output is files with names starting with minikraken_merged_taxaAsCol_ or minikraken_merged_taxaAsCol_logNorm.
minikraken_otuModel_China.R.txt -> calculates p-values for each taxa and plots the results. Its inputs are the log normalized tables from formatMinikrakenTable.R.txt. Its outputs are files whose names start with minikraken_otuModel.
minikrakenPCoA.R.txt -> draws PCoAs for each taxonomic level and outputs tables containing the PCoA site scores for each sample. Its inputs are the log normalized tables from formatMinikrakenTable.R.txt. Its outputs are minikrakenPCoA.pdf and tables whose names start with minikraken_PcoaCorrected.
minikraken_pcoaModel_China.R.txt -> calculates p-values for the PCoA axes and plots the results. Its inputs are the tables from minikrakenPCoA.R.txt. Outputs files whose names start with minikraken_pcoaModel.
minikraken_diversity.R.txt -> calculates the Shannon diversity index, Inverse Simpson diversity index, richness and evenness, and the p-values for each of these measures. It looks at both logged and unlogged data but only unlogged data is used in other analyses. Its inputs are from formatMinikrakenTable.R.txt. Rarefaction is not performed because WGS started with the same number of reads for each sample. Outputs files whose names start with minikraken_diversity.
minikraken_model_otu-metadata-metabolon.R.txt -> calculates the associations between minikraken taxonomic abundances, metadata and metabolite data and plots the results. Its inputs are metadata_cleaned11-2015.txt (from metadata/formatMetadata_2015-11.txt), metabolite/MetabolonScaledImpData-transpose.txt, and the log normalized tables from formatMinikrakenTable.R.txt. Outputs files whose names start with minikraken_model.

##HUMAnN Analysis##

This section refers to the files in the humann folder and contains the WGS gene content analysis using HUMAnN.

The files in the inputData folder contain the HUMAnN output files used.

Model for HUMAnN data: variable ~ urbanRural

The code for the analyses is in the analysis folder:
formatHumannTables.R.txt -> transpose the HUMAnN output and add metadata (from the output of 16SrRNA/analysis/relAbun.R). Note this script will also generate logged abundance tables but these were not used. Its inputs are the files in the inputData folder. Outputs files whose names start with humann_keggAsCol or humann_keggAsRow.
humann_otuModel_China_unlog.R.txt -> calculates p-values for each KEGG and plots the results. Its inputs are the outputs of formatHumannTables.R.txt. Its outputs are files whose name starts with humann_otuModel.
humann_pcoa_unlog.R.txt -> draws PCoAs and outputs tables containing the PCoA site scores for each sample. Its inputs are the tables from formatHumannTables.R.txt. It outputs humann_pcoa_unlog.pdf and tables with names starting with humann_PcoaCorrected_unlog.
humann_pcoaModel_China_unlog.R.txt -> calculates p-values for the PCoA axes and plots the results. Its inputs are the tables from humann_pcoa_unlog.R.txt. Its outputs are files whose name ends with humann_pcoaModel.
humann_vegan_diversity.R.txt -> calculates the Shannon diversity index, Inverse Simpson diversity index, richness and evenness, and the p-values for each of these measures. It looks at both logged and unlogged data but only unlogged data is used in other analyses. Its inputs are from formatHumannTables.R.txt. Rarefaction is not performed because we started with the same number of reads for each sample after whole metagenome sequencing. The names of its output tables start with humann_vegan_diversity.
humann_model_kegg-metadata-metabolon_unlog.R.txt -> calculates the associations between the KEGG modules/pathways and metabolite data or metadata. Its inputs are metadata_cleaned11-2015.txt (from metadata/formatMetadata_2015-11.txt), metabolite/MetabolonScaledImpData-transpose.txt, and the output of formatHumannTables.R.txt. Outputs files whose name starts with humann_model_unlog.

##Align WGS Reads to CARD or MvirDB databases##

This section refers to the files in the CARDandMvirDB folder and contains the results of aligning the whole metagenome reads to the CARD protein homolog database or MvirDB.

The file pro_homolog_results.txt contains the proportion of reads that aligned to the CARD protein homolog database.
The file MvirDB_propReads.txt contains the proportion of reads that aligned to MvirDB. MVirDB_results.txt, VFDBcore_results.txt and VFDBfull_results.txt contain the proportion of reads that aligned to each gene in MvirDB, VFDB core database or VFDB full database, respectively (the first couple of columns have the reads that mapped anywhere in the database).

cardsAssoc.R.txt -> analyzes the associations between the proportion of reads mapped to the CARD protein homolog database and the 16S, Kraken, HUMAnN, metadata or metabolite values. Its inputs are pro_homolog_results.txt, the files in 16SrRNA/inputData/RDP and 16SrRNA/inputData/abundantOTU and the output of 16SrRNA/analysis/relAbun.R.txt, the log normalized tables output from minikraken/analysis/formatMinikrakenTable.R.txt, the output of humann/analysis/formatHumannTables.R.txt, metadata_cleaned11-2015.txt from metadata/formatMetadata_2015_11.R.txt, and metabolite/MetabolonScaledImpData-transpose.txt. It outputs files whose name starts with cardsAssoc_ProHomolog_propMapped_v_ or cardsAssoc_ProHomolog_allGenes_v_, which looks at the association between the proportion of reads mapped anywhere in the database or at the association to all CARD genes individually, respectively. Note that running the genes individually is slow; you can comment out all calls to assocGenes if you wish to skip this analysis and only look at the proportion of reads that mapped anywhere in the database.
virulenceAssoc.R.txt -> analyzes the associations between the proportion of reads mapped to MvirDB and VFDB (full and core databases) and the 16S, Kraken, HUMAnN, metadata or metabolite values. Its inputs are MvirDB_results.txt, VFDBcore_results.txt and VFDBfull_results.txt, the files in 16SrRNA/inputData/RDP and 16SrRNA/inputData/abundantOTU, the log normalized tables output from minikraken/analysis/formatMinikrakenTable.R.txt, the output of humann/analysis/formatHumannTables.R.txt, metadata_cleaned11-2015.txt from metadata/formatMetadata_2015_11.R.txt, and metabolite/MetabolonScaledImpData-transpose.txt. It outputs files whose name starts with virulenceAssoc, which looks at the association between the proportion of reads mapped anywhere in the database or at the association to all genes individually (files with allGenes in the name). Note that running the genes individually is slow; you can comment out all calls to assocGenes if you wish to skip this analysis and only look at the proportion of reads that mapped anywhere in the database. Also note this looks at all three virulence databases, but only MvirDB (which includes VFDB) was used for later analyses.

##Accession Numbers##

The files in the accessionNumbers folder contain the tables from SRA (China_16S_srametadata.txt-processed-ok.tsv for the 16S reads and China_wgs_srametadata.txt-processed-ok.tsv for the WGS reads) and MG-RAST (2016-2-10table_all160_from_mgrast.txt; 16S reads only) containing the accession numbers for the sequence uploads. number16Sreads.txt contains the number of reads for each 16S file from BGI.

##Figures##

The code to generate the main figures can be found in the figures folder.

Figure 1 (Differences in microbial composition between urban and rural Chinese subjects) is generated by microbiomeFig_rarefy.R. This script takes as input the files in 16SrRNA/inputData and the outputs of rarefyDiversity.R, pcoaModel_16s.R, and classifyPlusROCcforest.R (all in 16SrRNA/analysis). The resulting figure is called microbiomeFig.tif.

Figure 2 (The urban Chinese microbiome is significantly more prevalent in an American cohort than the rural Chinese microbiome) is generated by microbiomeHeatmapPlusDBFig.R. This script takes as input the files in 16SrRNA/inputData and the outputs of 16SrRNA/analysis/otuModel_16s.R.txt. The resulting figure is called microbiomeHeatmapPlusDB.tif.

Figure 3 (Differences in metabolites between urban and rural Chinese subjects) is generated by metabolonFig.R. This script takes as input the files in the metabolite folder and the outputs of 16SrRNA/analysis/classifyPlusROCcforest.R and metabolite/metabolonModel.R. The p-values reported in the legend come from metabolite/pcoaModel_metabolonMetadata.R. The resulting figure is called metabolonFig.tif.

Figure 4 (Viruses and Archaea are less abundant in the urban Chinese microbiome while Bacteria are more abundant) is generated by minikrakenFig.R. This script takes as input the outputs of formatMinikrakenTable.R, minikraken_pcoaModel_China.R.txt, minikraken_diversity.R.txt and minikraken_otuModel_China.R.txt (all in minikraken/analysis). The resulting figure is called minikrakenFig.tif.

Figure 5 (The microbiome of urban Chinese subjects encodes more genetic diversity, including more genes conferring antibiotic resistance, than rural Chinese subjects) is generated by humannFigWithAbxAndVirHoriz.R. This script takes as input the outputs of formatHumannTables.R, humann_pcoaModel_China_unlog.R, humann_vegan_diversity.R (all in humann/analysis), as well as the tables in the CARDandMvirDB folder and the output of of 16SrRNA/analysis/relAbun.R to add metadata. The resulting figure is called humannFigWithAbxAndVirHoriz.tif.

##Supplemental Figures##

The code to generate the supplemental figures can be found in the supplementalFigures folder.

Supplemental Figure S1 (Correlation between timepoints) is generated by T1vT2SuppFig.R. This script takes as input the tables from 16SrRNA/analysis/pcoaWithDiversity.R. The resulting figure is called T1vT2SuppFig.tif.

Supplemental Figure S2 (Differences between urban and rural microbial composition is not dependent on bioinformatics pathway) is generated by qiimeSuppFig.R. This script takes as input the files in the 16SrRNA/inputData/qiime folder. The resulting figure is called qiimeSuppFig.tif.

Supplemental Figure S3 (Use of an alternative classifier confirms the ability to predict urban/rural status from the microbiome, the metabolome, or metadata) is generated by rpartROCSuppFig.R. This script takes as input the output of 16SrRNA/analysis/classifyPlusROCrpart.R.txt. The resulting figure is called rpartROCSuppFig.tif.

Supplemental Figure S4 (Differences in microbial diversity based on 16S rRNA sequencing are driven by differences in evenness and are opposite differences in richness) is generated by microbiome16sDivSuppFig_unlogRarefy.R. This script takes as input files from 16SrRNA/inputData/RDP and the output of 16SrRNA/analysis/rarefyDiversity.R.txt. The resulting figure is called microbiome16sDivSuppFig_unlogRarefy.tif.

Supplemental Figure S5 (Prevalence across the data set) is generated by rollingWindowPrevalence.R. This script takes as input 16SrRNA/inputData/databaseResults/rollingWindownPrevelance.txt. The resulting figure is called rollingPrevalenceVp.tif.

Supplemental Figure S7 (Variables measured in both the metadata and metabolite analyses are strongly correlated) is generated by sameVarMetadataMetabolite.R. This script takes as input metadata_cleaned11-2015.txt from metadata/formatMetadata_2015-11.R.txt, metabolite/MetabolonScaledImpData-transpose.txt, and model_metadata_v_metabolon_pValues.txt from 16SrRNA/analysis/model_otu-metadata-metabolon.R.txt. The resulting figure is sameVarMetadataMetaboliteCorr.tif.

Supplemental Figure S8 (Differences in microbial diversity based on whole genome sequencing is not dependent on diversity index used) is generated by minikrakenDivSuppFig_unlog.R. This script takes as input 16SrRNA/inputData/RDP/phylum_taxaAsColumnsLogNorm_withMetadata.txt (for metadata) and the outputs of minikraken/analysis/minikraken_diversity.R.txt. The resulting figure is called minikrakenDivSuppFig_unlog.tif.

Supplemental Figure S9 (Taxa significantly different in abundance between urban and rural subjects in whole genome sequencing have higher abundance in rural subjects) is generated by minikraken_volcanoPlot.R. This script takes as input the output of minikraken/analysis/minikraken_otuModel_China.R.txt. The resulting figure is called minikraken_volcanoPlot.tif.

Supplemental Figure S10 (Urban subjects have higher gene diversity and evenness but equal richness compared to rural subjects) is generated by humannDivSuppFig_unlog.R. This script takes as input the outputs of humann/analysis/humann_vegan_diversity.R.txt. The resulting figure is called humannDivSuppFig_unlog.tif.

Supplemental Figure S11 (KEGG pathways significantly different in abundance between urban and rural subjects tend to have a higher abundance in rural subjects) is generated by humann_volcanoPlot.R. This script takes as input the output of humann/analysis/humann_otuModel_China_unlog.R.txt. The resulting figure is called humann_volcanoPlot.tif.

Supplemental Figure S12 (The proportion of reads that map to genes that confer antibiotic resistance is associated with Escherichia and Shigella) is generated by abx_v_eschshig.R. This script takes as input CARDandMvirDB/pro_homolog_results.txt, the files in 16SrRNA/inputData/RDP, and the output of cardsAssoc.R.txt, minikraken/analysis/formatMinikrakenTable.R.txt. The resulting figure is called abx_v_eschshig.tif.

##Table##

The model results in Table 1 (Characteristics of the populations studied) came from metadata/metadataModel.R.

##Supplemental Tables##

The code to generate the supplemental tables can be found in the supplementalTables folder. Each individual text file output was made into a tab in the final Excel spreadsheet.

Supplemental Table S1 (Sequencing statistics and metadata for each sample) is generated by seqStatsTable.R. This script takes as input 16SrRNA/inputData/RDP/phylum_taxaAsColumnsLogNorm_WithMetadata.txt for metadata and the files in the accessionNumbers folder. Output table is called seqStats.txt.

Supplemental Table S2 (Model results for all microbial taxa identified through 16S rRNA sequencing, metabolites, and metadata) is generated by SuppTableModelResults_MicrobiomeMetabolomeMetadata.R. This script takes as input 16SrRNA/inputData/RDP/RDPtaxonomy.txt, 16SrRNA/inputData/abundantOTU/abundantOTU.chinaForward.taxonomy.txt, the output of 16SrRNA/analysis/otuModel_16s.R.txt, metabolite/metabolonPathwayInfo.txt, the output of metabolite/metabolonModel.R.txt, metadata/metadataConversion.txt, and the output of metadata/metadataModel.R.txt. Output table are SuppTableModel_metabolome and SuppTableModel_metadata or have names starting with SuppTableModel_16S.

Supplemental Table S3 (Significant associations between microbial taxa identified through 16S rRNA sequencing, metabolites, and metadata) is generated by AssocTableSupp_16s_metabolite_metadata.R. This script takes as input 16SrRNA/inputData/abundantOTU/abundantOTU.chinaForward.taxonomy.txt, metabolite/metabolonPathwayInfo.txt, metadata/metadataConversion.txt, and the output of 16SrRNA/analysis/model_otu-metadata-metabolon.R.txt. The output tables are assocTable_metabolite-v-16S.txt, assocTable_metadata-v-16S.txt and assocTable_metadata-v-metabolite.txt.

Supplemental Table S4 (Number of taxa and genes identified in analysis) is generated by numTaxa.R. This script takes as input the files in 16SrRNA/inputData/RDP, 16SrRNA/inputData/abundantOTU/abundantOTUForwardTaxaAsColumnsLogNormalWithMetadata.txt, the log normalized tables from minikraken/analysis/formatMinikrakenTable.R.txt and the output of humann/analysis/formatHumannTables.R.txt. The output table is numTaxa.txt.

Supplemental Table S5 (Model results for whole genome sequencing data) is generated by SuppTableModelResults_minikraken.R and SuppTableModelResults_humann.R. SuppTableModelResults_minikraken.R takes as input the output of minikraken/analysis/minikraken_otuModel_China.R.txt and the files in minikraken/inputData. The output tables have names starting with SuppTableModel_minikraken. SuppTableModelResults_humann.R takes as input the output of humann/analysis/humann_otuModel_China_unlog.R.txt. The output tables have names starting with SuppTableModel_humann.

Supplemental Table S6 (Significant associations between whole genome sequencing results and metabolites or metadata) is generated by AssocTableSupp_wgs.R. This script takes as input metadata/metadataConversion.txt, metabolite/metabolonPathwayInfo.txt, the output of minikraken/analysis/minikraken_model_otu-metadata-metabolon.R.txt, the files in minikraken/inputData the output of humann/anlaysis/humann_model_kegg-metadata-metabolon_unlog.R.txt. It outputs tables whose names start with assocTable.

Supplemental Table S7 (Significant associations with antibiotic resistance) is generated by abxAssocTable.R. This script takes as input 16SrRNA/inputData/abundantOTU/abundantOTU.chinaForward.taxonomy.txt, the output of CARDandMvirDB/cardsAssoc.R.txt, and the files in minikraken/inputData. It outputs tables whose names start with abxAssocTable.

Supplemental Table S8 (Significant associations with virulence genes) is generated by virAssocTable.R. This script takes as input 16SrRNA/inputData/abundantOTU/abundantOTU.chinaForward.taxonomy.txt, the output of CARDandMvirDB/virulenceAssoc.R.txt, and the files in minikraken/inputData. It outputs tables whose names start with virAssocTable.
