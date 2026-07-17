README for Burnett et al.

This set of scripts and inputs generates the SF estimates reported in the article for 33 Pisonia grandis trunks on Tetiaroa Atoll.

Script inputs can be modified for use with digitized cross-sections for new trees. See the article for a description of how to digitize trunk cross-sections.

### Input files ###
1. smoothed_polygons_aligned_clean.gpkg
	- This geopackage includes vector (polygon) representations of the 33 digitized Pisonia trunks.
	- Each polygon is associated with a unique "tree_ID"
	- UTM zone 1 (EPSG 32601)
2. swd.csv
	- This file includes tree_ID values and a sapwood depth value in cm associated with each individual.
	
### Scripts ###
Scripts should be run in the following order, with input/output filepaths modified to reflect the user's filesystem.

1. trunk-data-generator.R
	- Generates an array of geometric properties for each trunk polygon and saves them to "trunks.gpkg" (this filename can be modified)
2. trunk-erosion.R
	- This script requires radial profiles of sapflux density be described by a beta distribution. Parameters u and K are set by the user (default values are those fit for Pisonia in the article)
	- Creates a new geopackage file for each original polygon within a new "/eroded" directory (this directory can be modified)
	- Each new geopackage file contains the nested annuli of thickness dx (0.5 mm by default) from each polygon and calculates the sapflow (volume/time) contributed by each annulus, according to the beta sapflux distribution entered by the user
	- Mode can be changed from default "sapwood" to "full_radius" if the beta sapflux distribution reflects the entire range of depths from cambium to R_circ rather than the range of depths from cambium to sapwood depth.
	- This script uses packages foreach and doParallel to speed up processing.
	- Sapflux densities are modeled using unscaled beta distributions (i.e. the area under the distribution adds up to 1). Configurable parameter A is a scalar that can be used to introduce real units into the calculation; as a simple example, if the sap flux density is known at a certain depth (J_true), the beta distribution sap flux density value at that depth can be upscaled using A if A=J_true/J_beta. SF values will then be estimated in the units of J_true. The script is written such that the same A must be applied to all trees.
3. r-steiner-r-in-computation.R
	- This script spatially estimates r_steiner and r_in from the eroded geopackages, then writes them to "trunks.gpkg"
4. sapflow-estimates.R
	- Estimates SF_circ, SF_steiner, and SF_corr for all the polygons in trunks.gpkg
	- By default, SF_corr is estimated not with the r_steiner and r_in values measured by the previous script, but instead with modeled values from Eqs. 10-13 in the article. Measured values can be used instead by changing "use_regression_for_r" to FALSE; regression coefficients can also be modified in the "toggle & coefficients" section. The script uses Eqs. 10-13 by default because in the article we were examining the overall uncertainty in SF_corr if the models are used for r_steiner and r_in, as they would be in a new tree. Eqs. 10-13 were fit using the data from the previous script.
	- Mode can be changed from default "sapwood" to "full_radius" if the beta sapflux distribution reflects the entire range of depths from cambium to R_circ rather than the range of depths from cambium to sapwood depth. Stay consistent with the mode used in trunk-erosion.R
	- Scaling parameter A should be the same in this script as in Script #2.

### Outputs ###
1. /eroded/01_eroded_v1.gpkg
	- The "eroded" directory will include a geopackage file for every unique polygon present in the original input. Each of these geopackages will include many layers, each reflecting one eroded annulus used to calculate SF_ref.
	- Each geopackage will include trunk-level data fields (which are the same for every layer):
		1. tree_id
		2. perimeter_cm
		3. hull_perimeter_cm
		4. convexity
		5. area_cm2
		6. hull_area_cm2
		7. area_ratio [ratio of trunk area to its convex hull area]
		8. DBH_cm [DBH as estimated as hull_perimeter divided by pi]
		9. perimeter_to_area
		10. deepest_concave_depth_cm [this is h_max in the article]
		11. concave_depth_over_dbh [this is h_rel in the article]
		12. sapwood_depth_cm [as provided in swd.csv]
		13. r_sapwood [sapwood_depth_cm as a proportion of R_circ=DBH/2]
	- Each geopackage will also include annulus-level data fields (which are different for every layer):
		14. layer [counter for annuli starting at 1 for the outermost]
		15. layer.distance.from.edge [distance between each annulus' outer boundary and the trunk's original outer boundary, in mm]
		16. layer.area.cm2 [area of annulus in square cm]
		17. layer.area.mm2 [area of annulus in square mm]
		18. layer.perimeter.cm [outer perimeter of annulus in cm]
		19. swd.cm [redundant with sapwood_depth_cm]
		20. relative.sfd [sapflux density at the middle-depth of each annulus]
		21. layer.relative.sapflow [the volumetric sapflow contributed by each layer]
		22. relative.depth.to.swd [proportion of layer depth to sapwood depth]
2. trunks.gpkg
	- Each layer represents one trunk polygon with the following data fields:
		1. tree_id
		2. perimeter_cm
		3. hull_perimeter_cm
		4. convexity
		5. area_cm2
		6. hull_area_cm2
		7. area_ratio [ratio of trunk area to its convex hull area]
		8. DBH_cm [DBH as estimated as hull_perimeter divided by pi]
		9. perimeter_to_area
		10. deepest_concave_depth_cm [this is h_max in the article]
		11. concave_depth_over_dbh [this is h_rel in the article]
		12. sapwood_depth_cm [as provided in swd.csv]
		13. r_sapwood [sapwood_depth_cm as a proportion of R_circ=DBH/2]
		14. r_steiner [relative depth of first significant topological change, as calculated by Script #3]
		15. swd_minus_r_steiner [r_sapwood minus r_steiner; useful for examining whether sapflow occurs deeper than r_steiner]
		16. incircle_R_cm [radius of largest possible inscribed circle in cm]
		17. R_cm [R_circ in cm; =DBH/2]
		18. r_in [relative incircle radius]
		19. sapflow_ref [SF_ref estimated from eroded annuli in Script #4]
		20. r_in_pred [modeled r_in from Eq. 13 in the article]
		21. r_steiner_pred [modeled r_steiner from Eq. 11 in the article]
		22. sapflow_circ [SF_circ as modeled in Script #4]
		23. sapflow_steiner [SF_circ as modeled in Script #4]
		24. sapflow_corr [SF_corr as modeled in Script #4]
		25. SF_circ_bias [Bias of SF_circ relative to SF_ref]
		26. SF_steiner_bias [Bias of SF_steiner relative to SF_ref]
		27. SF_corr_bias [Bias of SF_corr relative to SF_ref]
	- Note that the SF units are determined by the A scalars in Scripts #2 and #4. If the default A=1 is used, the SF units are aribitrary (but technically cm^2) because sap flux density is unitless (as described by the Beta PDF). If the beta distribution of sap flux density is scaled such that its units are cm^3 cm^-2 hr^-1, for instance, then final SF units are cm^3/hr.
	
	
	
### Other files ###
1. pisonia-modeled-trunk-data-2025-12-18.csv
	- This spreadsheet includes all the primary and modeled data from the 33 digitized Pisonia trunks.
	- Should mimic the final data from running the above scripts on the Pisonia data.
