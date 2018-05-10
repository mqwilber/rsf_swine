##  Summarize which pigs overlap with crops

library(data.table)
library(yaml)
library(ggplot)
library(lubridate)

files = Sys.glob("~/Repos/rsf_swine/results/glmdata_by_study/*.csv")
anal_params = yaml.load_file("~/Repos/rsf_swine/code/analysis_parameters.yml")
studysum = fread("~/Repos/rsf_swine/data/formatted/study_summary.csv")
pigattrib = fread("~/Repos/rsf_swine/data/formatted/pig_attributes.csv")

datasets  = list()
for(fl in files){

	studynm = strsplit(basename(fl), ".", fixed=T)[[1]][1]
	dat = fread(fl)
	datasets[[studynm]] = dat
}

############### Time-overlap analysis ###################

# Bind ecoregion to each study
for(studynm in names(datasets)){
	datasets[[studynm]][, month:=month(date)]
	datasets[[studynm]][, ecoregion:=studysum[study == studynm, l2ecoregion]]
}

# Calculate the months of each study
unqmonths = do.call(rbind, lapply(datasets, function(dt) dt[, list(unqmonth=unique(month)), 
																								by=list(pigID, ecoregion)]))

# Plot the time span in which pigs were collared across ecoregion
dir.create("../results/time_span_plots")
for(eco in unique(unqmonths$ecoregion)){
	print(eco)
	tmonths = unqmonths[ecoregion == eco]
	tplot = ggplot(tmonths) + 
											geom_point(aes(x=unqmonth, y=pigID, color=pigID)) + 
											facet_wrap(~ecoregion, scale="free")
	ggsave(paste0("../results/time_span_plots/", gsub("/", "-", eco), ".pdf"), tplot)

}

################# Crop-use analysis ####################

# Total pigs
totpigs = sum(sapply(datasets, function(x) length(unique(x$pigID))))

# How many pigs are ever in crops?
used_crops = list()
pinc = list()
for(studynm in names(datasets)){

	tdat = datasets[[studynm]]

	# Figure out which pigs were in crops
	incroppigs = tdat[z == 1][, list(incrop=any(crop_loc == 1)), by=pigID][incrop == TRUE, pigID]
	pinc[[studynm]] = tdat[pigID %in% incroppigs]

}

allpigs = do.call(rbind, pinc)
allpigs_sex = merge(allpigs, pigattrib[, list(sex, pigID)], by="pigID", all.x=T)
allpigs_sex[z == 1 & crop_loc == 1][, list(length(unique(pigID))), by=sex]
# 89 males and 62 females using crops

cropusers = allpigs[, list(totpigs=length(unique(pigID))), by=study]
cropusers[, list(sum(totpigs))]

# Amount of crops used
totused = allpigs[z == 1, list(totused = sum(crop_loc == 1), 
															 timeused=sum(tau[crop_loc == 1]),
															 proptime=sum(tau[crop_loc == 1]) / sum(tau)), 
															 				by=list(pigID, study)][order(timeused, decreasing=T)]
totused[, list(meanprop=mean(proptime)), by=study]


# Of crop users, tx_tyler_w2 are using crops over the highest proportion of time.

used = allpigs[z == 1, lapply(.SD, function(x) sum(tau[x == 1]) / (60 * 60)), by=study, 
										.SDcols=paste0(anal_params$croptypes, "_loc")]

mused = melt(used, measure.vars=paste0(anal_params$croptypes, "_loc"), 
							variable.name="croptype", value.name="hours_used")
mused = mused[order(study, hours_used, decreasing=T), ]
mused = mused[croptype != "crop_loc"]

library(ggplot2)
ggplot(mused) + geom_bar(aes(x=study, y=log10(hours_used + 1), fill=croptype), 
													stat="identity") + ylab("log10(total pig hours used + 1)") +
								theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=4))
ggsave("~/Repos/rsf_swine/results/crop_use_by_study.pdf", width=5, height=4)

ggplot(mused) + geom_bar(aes(x=croptype, y=log10(hours_used + 1)), 
													position = "dodge", stat="identity") + ylab("log10(total pig hours used + 1)") +
								theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8))
ggsave("~/Repos/rsf_swine/results/crop_use_overall.pdf", width=5, height=4)


