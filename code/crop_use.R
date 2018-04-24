# Look at which pigs overlap with crops

library(data.table)
library(yaml)

files = Sys.glob("~/Repos/rsf_swine/results/glmdata_by_study/*.csv")
anal_params = yaml.load_file("~/Repos/rsf_swine/code/analysis_parameters.yml")

datasets  = list()
for(fl in files){

	studynm = strsplit(basename(fl), ".", fixed=T)[[1]][1]
	dat = fread(fl)
	datasets[[studynm]] = dat
}

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
													position = "dodge", stat="identity") + ylab("log10(total pig hours used + 1)") +
								theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8))
ggsave("~/Repos/rsf_swine/results/crop_use_by_study.pdf", width=5, height=4)

ggplot(mused) + geom_bar(aes(x=croptype, y=log10(hours_used + 1)), 
													position = "dodge", stat="identity") + ylab("log10(total pig hours used + 1)") +
								theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8))
ggsave("~/Repos/rsf_swine/results/crop_use_overall.pdf", width=5, height=4)


