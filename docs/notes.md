## RSF analysis notes and reading

### Goals for the the next two weeks

a. Read up on the RSF literature and identify potential methodological approaches that could be useful
b. Read up on RSF data on pig populations
c. Familiarize myself with data that is currently available

Methodological approaches for RSF analysis
-------------------------------------------

## Methodological questions

1. Do the approaches provide different results when looking at the RSF for pigs? 

2. Will incorporating random effects of individuals improve predictive estimates?
	- Others have included random effects (Northrup)

3. The real key 

3. Using this big data set to validate and PREDICT RSF functions on left out pigs.  Playing with model validation of RSF functions.

## Some potentially interesting approaches

1. Hooten et al. 2014, Statistical Methodology
	- Uses a Continuous time Correlated random walk model for animal movement model animal movement and then conditions on this movement model when estimating the RSF function. 
	- Good for GPS data collected on short intervals. 
	- Can impute over missing time intervals.
2. Nielsen and Sawyer, 2013, Ecology and Evolution
	- This approach seems very straightforward.  Using negative binomial regression to estimate the resource selection function when there are multiple movement points. 
	- Again, this is useful for GPS data where there are many points from the same individual.
3. Hanks et al. 2011 Plos One
	- A velocity-based approach for movement.  The approach relies on modeling how an animal moves over a resource gradient.  To they move to gradients and what is their speed when they hit a certain gradient.
	- I like this approach because it seems be reflecting animal movement a bit better than simply location data.
	- Also this approach makes it easy to impute over unequal time intervals/missing values.
4. Northrup et al. 2016, Ecological Applications
	- Uses the Hooten implementation.  Pretty cool approach I think.
5. Wilson et al. 2018 and Hooten et al. 2015
	- Again, using the CTCRW but also layering in a Markov model with umps between grids that depend on resource and covariates. Can include directional covariates as well as locations based depending on the scale of the data.  I like this approach because the covariates are DIRECTLY affecting the animal movement.  Not just ad-hocly being used to to infer the weighted distribution (not saying this is bad and might be necessary depending on the scale of interest!).  But this approach is very appealing for telemetry data.

Notes
-------------
1. While a Bayesian approach seems like it would be for this type of analysis (random effects, imputation, etc), the size of the data might preclude a bayesian analysis. At the same time, given the size of the data, some assumption can potentially be made about the distributions that various coefficients are going to follow (i.e. MVN).
	- Update. This might not be a problem! 

2. Man, there is a lot of material on resource selection functions!  

3. Northrup et al. 2016 suggests ommited the first week of GPS data to account for potential effects of capture behavior.

4. KEY CHALLENGE: Identifying the availability sample.  Hooten et al. 2014 and Northrup et al. 2016 provide what I think is the least ad-hoc approach so far: Use a dynamic model (CTCRW) to predict for each "used" point where general radius where an animal could go next.  This defines the availability sample FOR each data point.  I like this better than trying to find the availability sample around the hull of points.  

5. THOUGHT: I don't see what advantage this would provide over the CTCRW model, but it might be possible to think about using a Gaussian process regression to help define the availability sample for each observed point.  They would be a bit more flexible...I think that would be the main advantage, but probably way slower to fit.

6. "Impossible movement criteria": Bjornerass et al. 2010.  Useful to eliminate outlying GPS data.

7. Minimum convex polygons:  A method for measuring "home range" and available space. 

8. Boyce et al. 2002:  Provides a cross validation approach for testing the efficacy of RSFs.

9. Where to get habitat layers?
	- Natural Earth
	- USDA cropland data

10. Interesting thought:  There are models that do and do not account for spatial autocorrelation in the animal data.  A few of the papers make a big deal about accounting for autcorrelation, and I tend to agree.  However, if we simulate animal movement with resource selection incorporated (need to think about how to do this...definitely not straightforward) and then ignore the autocorrelation when calculating the RSF, what happens to our resource selection functions? 
	- This could be quite a cool simulation study

11. Revisit your idea of using Gaussian process regression to account for spatial and temporal autocorrelation (or at least temporal autocorrelation), to allow for the implementation of simple cross scale models.  See Hefley et al. 2017 or discussion of using basis functions to control for spatial autocorrelation.  This could be helpful to decorrelate observed time series.


Potential covariates
--------------------

1. Pig density (location covariate)
2. Distance to roads (directional covariates)
3. Distance to crops (directional covariate and location covariates)
4. Distance to water (directional covariate)
5. Forest cover (location covariate) [Lewis et al. 2017 suggest that forest cover is important. References there in.]
6. Climate?  Though Ryan Brook said this might not be as important as we initially thought.  Most studies still include it though.


TODOS
-----

1. Identify the questions that we want to ask using RSFs?  Or, really, the questions that we want to ask in general (perhaps shouldn't be dictated by RSFs).
	- I think the primary question should be about the RSF and how the RSF varies across individuals, populations, time, and space.
	- Methdologically, I think it makes a lot of sense to use the approach of Hooten to define AV

2. Explore the CTCRW approach by Johnson for some of the pig data.  Pull in some environmental covariates and potentially try to do some fitting?

3.  COOL APPROACH:  Individual, population, and metapopulation level!  Think about a random effect of population as well as individual. We should have the data for that. It would be nice in this case to have similar covariates across populations....but we might also have different ones and that could be important too.

4. Also, think about temporal/spatial scale analysis. Similar to what we see in a couple of the analyses that I have read about. 


Questions
---------

(Trying to gear them to be useful for management)

Seems like we have a pretty good handle on what drives pig density/presence at large scales.  But at the scale of an individual, how do pigs select resources? How do pigs move on a landscape?

This could be useful for 

A) Given an introduction by a human, predicting the most likely path a pig will take and likely established home ranges given introduction at a point.

B) Establishing how pig movement differs between endemic and invading populations? Can movement patterns describe WHY they are moving northward 12 km per year?  


CHALLENGES
----------

1. Pigs show strongly daily movement patterns (i.e. they sleep probably, lol).  Need to account for these patterns if I am going to try to use a continuous time model.

2. If I want to compare pig movement across habitats 


## Jan. 16, 2018

So it seems like two potential approaches are emerging: continuous-time animal movement models

Question: How do environmental covariates (i.e. resources) affect pig movement on the scale of the individual? How do these effects vary in a population? How do these effects vary across different populations? 

From a management context: How can this help us predict where pigs are/will be?

Broader ecological implications?: Broader ecological question that hasn't been addressed? 


1. Discrete-time models: Generally are measuring step length and bearing/turning angle.
	a. These are advantageous because they can model switches in states easily using a state-space approach.  This can account for resting, foraging, etc.  There is a also a lot of theory for them as well as a nice package to implement most of the fitting (moveHMM).
	b. Problem: While they allow us to say how various covariates effect movement, the don't to as good of a job projecting a utilization distribution/resource selection function onto the landscape/.  Though Hooten (of course) provides a few different approach for doing this in discrete time.

2. Continuous-time models: They are generally tracking velocity, which is inherently directional and can be integrated to give distance
	- The advantage to this approach (which seems useful in our case) is that inference does not depend on scale as much.  Because this approach is in continuous time, we can implicitly account for different time steps. 
	- So far, the most compelling approach that I have seen is the Hanks et al. UD approach:

		a. Model the continuous path of the animal (using CTCRW or potentially the general framework provided Hooten and Johnson or Buderman).  Buderman might be useful because it will allow us to account for the messy data (unequal sampling points, etc.) a bit better and I think will allow for faster data fitting. Fitting this model allows us to impute a continuous time movement path.
			- Approaches: Johnson et al. 2008 CTCRW, Hooten and Johnson 2016 general approach, Buderman et al. 2016 approach 

		b. Step 2: Impute the continuous path to discrete points (such that the animal stays in a cell for a number of time steps).  Then use a Poisson GLM to fit a continuous-time, discrete state Markov model (see Hanks et al. 2015) to estimate how the transition rate from cell j -> i depends on both "location" and "directional" covariates.  The stationary distribution of this Markov chain could be examined to give a sense of the utilization distribution.

	- Challenges: 

	1. Pigs behavior changes over the course of the day.  They seem to forage and rest, forage and rest, etc.  These latent states are really easy to account for in the discrete time setting, but a bit more challenging to account for in the Continuous time setting.

		a. Possible solutions: The Buderman approach might allow for us to account for these different movement regimes inherently.
		b. Johnson et al. 2008 shows how covariate data can be used to allow the auto-correlation parameter to vary with time.  Though I don't think we have this equivalent covariate for the pig data.
		c. McClintock et al. 2014: Might provide some insight in accounting for two states in the CTCRW by allowing Beta and Sigma to vary. 
		d. Beta and sigma could be really useful as additional comparative parameters between individuals and populations, but they are strictly necessary if we can use an approach like Buderman to impute the continuous time path.  Buderman seems to amount to fitting two separate B-splines to the different location parameters (at a function of time) and then using this to impute the continuous path.  A pretty straightforward approach I think. Will just need to think about error distributions and basis function scales.  Can this actually account for different movement patterns on the same scale? I think so.  Just a different type of spline on the smaller scale. 

	2. Covariates need to be comparable across pig populations
		- Pig density, Forest Cover, Agg type, distance to water, distance to roads, others?

	3. How does the length of the the time series between different bigs affect out ability to do joint inference?
		- Will need to weight individual-based parameters by sample size, unless there is some way it will be computationally feasible to do it. 


Moving forward:

1. Fit some discrete time state-switching model with the moveHMM package. First, ignoring covariates.  And then including covariates.

2. Extract 5 or so pigs from the texas data set.  Extract three covariates to go along with these pigs: 1. Cropland 2. Distance to water 3. Distance to road. Fit a discrete time HMM for the 5 pigs over the same time range to assess
	1. Do they show a clear latent state of resting vs. foraging?
	2. Do movement and/or turn angles depend on the environmental covariates?
	3. Do these parameters look the same across individuals?

Fit a continuous time model (I would say CTCRW, but FinAY that thing is finicky.) and apply the Hanks/Wilson methodology to estimate transition probabilities.

## January 19

Got the Hanks et al. model up and running using an interpolation approach to the continuous trajectory. There are a number of challenges that I see with this approach

1.  Is the interpolated trajectory sufficient or should I use the CTCRW.  Obviously the CTCRW can account for uncertainty, but I don't see much difference in simply using the interpolated trajectory (at continuous time steps) for inference. 
	- Man the CTCRW is not as useful as it is made out to be (in my opinion). When animals have distinct resting phases that are not directly measureable,
	it really struggles.  I think interpolation is ok, recognizing that the uncertainty the parameters will not be exact when you don't have the imputation.

2. The the movement of an individual animal depends on both resource (probably) AND a correlation with its previous direction and speed.  Does this model account for this?  Is it inherent in the trajectory? I.e. are we essentially conditioning on this autocorrelations when we compute the trajectory?
	- NO! We can include the past movement vector as a covariate in the model. Hanks does this and I think this is a great idea.

3. Use the generator matrix to simulate a movement trajectory.  Problem: this simulation assumes that resource is the ONLY thing driving movement. 

4. Implement time varying movement



## February 6, 2018

The movement model is implemented and now it is a matter of identifying questions that are interesting and extracting the covariates that are important. I had a couple thoughts on this.

First, in order to make the analysis comparable across populations I think we could focus on the following question

1. How does pig movement vary in response to cultivated resources/cultivated resource gradients and uncultivated resources/gradients?  Do these movements change across scales populations and do they interact with temperature, pig density, and water availability? 

	- To answer this question, I propose that we think about broad scale "uncultivated" resources in terms of NDVI/EVI, which is a rough proxy for productivity.  There are obviously some issues with this as this isn't going to really get at all of the uncultivated resources available to the pigs.  However, this will give a rough and consistent measure of "productivity" across populations such that we can build an uncultivated resource gradient on the landscape.  We can supplement this layer with a masting layer provided by Mikey. 
		- Sources: VegScape (USDA) for NDVI.  MODIS for EVI. Mikey for masting layer.
		- This could be a composite layer oF NDVI, masting, and other variables that define "non-cultivated" resource

	- Second, we also need to consider cultivated resources (e.g. agriculture). In this case, we want to make agriculture comparable across populations, so we have discussed converting crops into calories per bushel and looking at average yield per county for particular crop types in the US (see https://quickstats.nass.usda.gov/ for the data). In this way we could look at a "calorie landscape" for cultivated crops and, at least coarsely, compare across different crop types.  The biological rational is that if pigs are generalists, than perhaps they are selecting for the caloric potential of a resource patch.  This would also allow us to explore how cultivated calories affects movements (or even if it does), how this interacts with temperature, and how this changes with time.
		- Will need to extract measures of calories per pound for each type of crop.
		- 

	- We will also need to include a proxy for cover: i.e. forest habitat.  Multiple studies have shown that this is important for resting behaviors. I think the easiest way to do this will be to use the CropScape data to classify either forest or non-forest across different populations.

	- The questions

	1. How does pig movement respond to the cultivated and non-cultivated resources? Do movement patterns/selection of these resources change over time?

	2. How do biological factors such as pig density effect the way in which pigs interact with cultivated and non-cultivated resources?

	3. How do abiotic factors such as temperature and availability of water affect the movement patterns/selection of pigs?

	4. Are these effects consistent across populations, suggesting a consistent way in which pigs select and use resource?

	5. Finally, for Sarah's question, how do individual attributes of (weight, sex, reproductive status) affect the selection of resources by pigs?


## Feb. 9, 2018

Goals for today

1. Extract NDVI, Cropscape, and Temperature data for Michigan pigs. 
	- Didn't extract NDVI.  VegScape was being tricky for bulk extractions.
2. Clean and format Michigan pigs. [DONE]
3. Run a comparison between Michigan pigs and Tx Camp pigs. 
4. Account for time varying NDVI and include it in model. [DONE]
5. Compute distance2gradient calculation [DONE]
	- Need to decide on appropriate distance decay function
6. Read at least 2 papers on feral swine (movement/influenza).

## Feb. 12, 2018

1. Finalize data summary and  [Done]
2. Talk with Sarah about covariates.  Check if they are done. [Done]
3. Run comparison between Michigan and txcamp pigs
4. Read papers on feral 

## Feb. 20, 2018

1. Find a suitable EVI covariate and compare results to NDVI.
2. Prepare precipitation covariates.
3. Split crop types into groups Ryan suggested.
4. Have Sarah start preparing distance to road layers.
5. Extract and prepare forest canopy covariate.
6. Figure out how to speed up the distance to water calculation.

**Meeting with Sarah**

1. Distance to road covariates

	- For next meeting, do you think that you could create a raster file in which the bounds are for "txcamp" and each cell has the distance to the nearest road (linear distance)? 
		- Save the ArcGIS code for how you did it.

2. Try to find a drought index covariate.

3. Distance to wetlands covariate.

4. For mark: Figure out EVI layer...a pain

5. Discuss cleaning data and removing pigs
	- Judas Pigs and Michigan.  Should we even include these?
	- Include or drop Judas pigs. 
		- We'd be missing pigs from Indiana...but do you think that is reasonable?
		- Email Justin about collar ID [Sara].
		- Talk to Kurt about the Michigan pigs [see if they are actually usable]. 

Updates:

1. Can we capture season? 
	- Is pressure an important covariate in the Kay et al paper?
2. Temperature, precip, snow depth, drought, and pressure potentially define "season"?
	- Kay et al. used pressure as a variable.
3. Ryan Brook pig weights.

TODOs for Sarah:

1. Distance to road covariates
2. Drought Index
3. Search for an EVI layer that is time varying and fine spatial res.
4. Contact Justin
5. Add literature into Dropbox

TODOs for Mark

1. Format CropScape data according to Ryan's suggestions [Done]
2. Format Forest cover layer
3. Format precipitation data [Done]
3. Contact Kurt about Michigan data
4. Try to put together a population-level analysis for Texas and California.
5. Fix data that Kilgo mentioned.






