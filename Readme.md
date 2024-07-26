## Source codes for Myctophid DVM study

**Software required: R (version 4.4.0), QGIS (3.38.1)**\

The "Myctobase" directory contains data downloaded from the freely available [Myctobase](https://zenodo.org/records/6562776). Three points in "event_edit.csv" have been modified from the original "event.csv" based on literature. Modifications made are stated in the "Notes" column.

The "GEBCO_05_Dec_2023_c97a092c1373" directory contains circumpolar ocean bathymetry data between 40-80°S downloaded from the freely available [GEBCO_2022 Grid](https://www.bodc.ac.uk/data/published_data_library/catalogue/10.5285/e0f0bb80-ab44-2739-e053-6c86abc0289c/) on 05/12/2023.\

It is recommended to run the codes in the following order to ensure all prerequisites are produced:

1.  *Data_wrangling.R*

    -   Carries out data quality control and calculates net effect variable
    -   Produces supplementary figure S1
    -   Exports: "event_edited2.csv", "Net_summary.csv" (supplementary table S3)

2.  *Abundance_count.R*

    -   Identifies the most abundant myctophid species in the data
    -   Exports: "Abundance.csv" (supplementary table S2)

3.  *Data_wrangling_for_QGIS.R*

    -   Produces a csv file (Net type.csv) and a shape file directory (Ocean_extent) for making figure 1 in QGIS
    -   Produces csv files (*species name*\_CPUE.csv) for making supplementary figure S3 in QGIS
    -   Exports: "Net type.csv", "Ocean_extent", "*species_name*\_CPUE.csv" for each species

4.  *Net type.qgz* (for use in QGIS)

    -   Produces figure 1
    -   Prerequisites: [Quantarctica3](https://www.npolar.no/quantarctica/)

5.  *Presence-absence.qgz* (for use in QGIS)

    -   Produces supplementary figure S2
    -   Prerequisites: [Quantarctica3](https://www.npolar.no/quantarctica/)

6.  *Distribution_raincloud_plot.R*

    -   Produces figure 2 (exported as a pdf file)

7.  *Model_comparison.R*

    -   Tests different model structures (log-transformed vs untransformed abundance, Tweedie vs Negative Binomial distribution)
    -   Produces supplementary figure S2
    -   Exports: "Model_comparison.csv" (supplementary table S2)

8.  *DVM_analysis.R*

    -   Runs the species-specific GAMs and all-species GAM
    -   Produces figure 3, figure 4, figure 5
    -   Exports: "DVM_pattern.csv" (table 1), "GAM_result.csv" (supplementary table S5)

9.  *Total_biomass_analysis.R*

    -   Calculates total abundance and biomass of myctophids in the Southern Ocean
    -   Estimates average standard length and weight of each species
    -   Exports: "Total_biomass_estimate.csv" (table2), "Size_summary.csv" (supplementary table S4)

10. *Fish_size_differences.R*

    -   Visualizes differences in the standard length of fish caught against net type, latitude and depth
    -   Produces supplementary figure S4, figure S5, figure S6
