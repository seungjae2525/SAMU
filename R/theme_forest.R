theme_forest <- function(base_size=13,
                         base_line_size=base_size / 22,
                         base_rect_size=base_size / 22){
  theme_minimal(base_size=base_size,
                base_line_size=base_line_size,
                base_rect_size=base_rect_size) +
    theme(rect=element_blank(),
          text=element_text(colour="black")) %+replace%
    theme(plot.title=element_text(face="bold", hjust=0),
          axis.text.y=element_text(colour="black"),
          axis.text.y.right=element_text(hjust=1),
          axis.text.x=element_text(colour="black"),
          panel.border=element_blank(),
          strip.text=element_text(face="bold", hjust=0),
          strip.background=element_rect(colour="white", fill="aliceblue"),
          panel.background=element_rect(colour=NA, fill=NA),
          panel.grid.major.x=element_line(colour="gray50", linewidth=0.25, linetype=2),
          panel.grid.minor.x=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor.y=element_blank())
}
