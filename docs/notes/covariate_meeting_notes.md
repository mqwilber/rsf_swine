## Notes from covariate meeting on 1-25-2018

### Comments on density

- Definitely include density as a covariate
- Colleen suggested using center of gravity measure [to other pigs?] as a measure of density-dependent interactions.
- SQUIDIS data: potentially determining how long has been present in an area?
	- Could be correlations between this and density measures.
- MIS data from Kim: County level measures of covariates.
	- Could be useful for the cross-study comparisons

## Temperature-dependent covariates

- Ryan has R code to calculate temperature-dependent covariates: EMAIL HIM [DONE]


## Agricultural covariates

The caloric approach seems challenging...people thought it was an interesting 
idea, but it might be tough to implement.

**Challenges**

	- If we look at the caloric value of crops, should we also look at the caloric value of the rest of the landscape?
	- Ryan suggested using bushel per acre as a measure of productivity of an area and then converting that into calories could be a good rough estimate.
	- Kurt suggested using data on intake rate for domestic pigs to help think about the interaction between temperature and resource requirement.

**Data sources**

- NAS or Flaps for agricultural data. Also CropScape has fine scale crop data across the country.
- Agricultural data will vary seasonally as well.  Will we be able to see a shift in how pigs select for agricultural and how that varies with time?
- 

- Question: How are pigs moving about the landscape relative to the resources that are available?

## Foraging covariates

- Hard masting layer: Mikey Tabak is building a masting data layer that could be really, really useful. Definitely going to need to include that.

- This does not include soft masting (not sure what that is...) but that could be important in Florida, for example.
- EVI: A time varying covariate based on LIDAR that could be really useful for measuring forest cover more quantitatively.


## Water-based covariates

- Distance water will be important
- But also, we might want to also include wetland layers as covariates. This would account for type of wetland.
- Colleen suggested using connectivity: Given a pig, how connected are they to the various water ways (IFM connectivity that weights size of water patch as well)
- Presence or absence of cattle/livestock as a predictor of stock tanks and water availability.
- Snow depth

## Home-range measures

- Pigs have dynamic home ranges.  Frankly, I didn't quite follow how we might be able to account for this...ask Colleen! 
- Distance to other collared pigs.
- Center of mass measures of home range...again I didn't totally follow this.
- Colleen also suggested using predictors at different scales and then choosing which scale lead to the best model.

## Capture effects

Use the data to trim out capture effects. Try different time ranges and see if that affects inference.

## Interspecific interactions

- Large mammal diversity layer
- Predator diversity layer
- Both of these seem important.
- Hunting regulations
	- Number of permit issued by state, number of pigs taken by density...all potential covariates for hunting.
- Specify a "hunting" hidden variable that can turn on or off with time.

## Model validation

- Might be hard to validate the model at a national scale at we really only have 480 pigs.
- Could validate the model on the local scale...e.g. test at camp bullis.
- Need to think about what model validation actually looks like.

## Other points

1. Interactions
	- Will need to think carefully about which interactions we will want to include in the model.
	- Snow and temperature
	- Water and livestock
	- Etc.

2. What level of Ecoregion should we use?



