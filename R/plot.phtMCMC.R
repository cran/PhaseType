plot.phtMCMC <- function(x, ...) {
	# Get chain info
	StartEndThin <- attr(x$samples, "mcpar")
	
	# Trace plots
	print(
		ggplot(melt(as.data.frame(x$samples[seq(StartEndThin[1], StartEndThin[2], StartEndThin[3]),]), id.vars=NULL)) +
		geom_line(aes(x=1:length(value), y=value)) +
		geom_smooth(aes(x=1:length(value), y=value)) +
		#geom_hline(aes(yintercept=value), data=truth, colour="red", linetype="dashed") +
		facet_wrap(~variable, scale="free") +
		#theme_grey(base_family="serif", base_size=11) +
		opts(title = "Parameter Traces") + #, plot.title = theme_text(size=14, face="bold", family="serif")) +
		xlab("Iteration") + ylab("Parameter Value")
	)
	
	# Marginal posterior densities
	print(
		ggplot() +
		geom_density(aes(x=value), melt(as.data.frame(x$samples[seq(StartEndThin[1], StartEndThin[2], StartEndThin[3]),]), id.vars=NULL)) +
		#geom_vline(aes(xintercept=value), data=truth, colour="red") +
		facet_wrap(~variable, scale="free") +
		#theme_grey(base_family="serif", base_size=11) +
		opts(title = "Marginal Posterior Densities") + #, plot.title = theme_text(size=14, face="bold", family="serif")) +
		xlab("Parameter Value") + ylab("Density")
	)
}
