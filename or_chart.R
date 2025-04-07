library(ggplot2)
library(RColorBrewer)

out_dir = "filepath"

chart_data_any <- read.csv(paste0(out_dir, "\\OR_chart_data.csv"))
chart_data_lim <- read.csv(paste0(out_dir, "\\OR_chart_data_lim.csv"))

chart_data_any$outcome <- "Any severity"
chart_data_lim$outcome <- "Activity-limiting"
chart_data <- rbind(chart_data_any, chart_data_lim)

chart_data$outcome <- as.factor(chart_data$outcome)
chart_data$outcome <- factor(chart_data$outcome, levels=c("Any severity", "Activity-limiting"))
chart_data$model <- factor(chart_data$model, levels=c("Adjusted", "Unadjusted"))
chart_data$vaccine <- factor(chart_data$vaccine, levels=c("mRNA", "Adenovirus vector", "Combined"))

chart_data$label <- paste0(
  format(round(chart_data$or, 2), nsmall=2),
  " (",
  format(round(chart_data$lcl, 2), nsmall=2),
  " to ",
  format(round(chart_data$ucl, 2), nsmall=2),
  ")"
)

chart_data_adj <- chart_data[chart_data$model=="Adjusted",]

dodge <- 0.7

or_chart <- ggplot(chart_data_adj, aes(x=or, y=vaccine, group=outcome, colour=outcome, label=label)) +
  geom_point(position=position_dodge(dodge), size=2.5) +
  geom_errorbarh(aes(xmin=lcl, xmax=ucl), position=position_dodge(dodge), height=0) +
  #geom_text(aes(x=0.97), hjust=1, vjust=1.3, position=position_dodge(dodge), size=3.5) +
  geom_text(aes(x=or), hjust=0.5, vjust=-1, position=position_dodge(dodge), size=3.0) +
  coord_trans(x="log") +
  scale_x_continuous(breaks=seq(0.3,1.0,0.1), labels=sprintf("%.1f", seq(0.3,1.0,0.1)), lim=c(0.3,1.0)) +
  geom_vline(xintercept=1, linetype="dashed", colour="grey50", size=0.5) +
  #scale_colour_manual(values=c("#084594", "#9ECAE1")) +
  scale_colour_manual(values=c("black", "grey70")) +
  xlab("Odds ratio (log scale)") +
  theme(
    axis.title.x=element_text(size=10, colour="black", face="plain"),
    axis.text.x=element_text(size=10, colour="black", face="plain"),
    axis.ticks.x=element_line(size=0.5, colour="black"),
    axis.line.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.y=element_text(size=10, colour="black", face="plain"),
    axis.ticks.y=element_blank(),
    axis.line.y=element_blank(),
    panel.border=element_rect(size=0.5, colour="black", fill=NA),
    panel.background=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.spacing=unit(1.5, "lines"),
    plot.margin=margin(0.2,1,0.2,0.5, unit="lines"),
    strip.background=element_blank(),
    strip.text=element_text(size=10, colour="black", face="bold"),
    legend.title=element_blank(),
    legend.position="bottom",
    legend.direction="horizontal",
    legend.justification="center",
    legend.text=element_text(size=10, colour="black", face="plain")
  )

ggsave(plot=or_chart,
       filename=paste(out_dir, "\\OR_chart.jpg", sep=""),
       width=14,
       height=12,
       units="cm")
