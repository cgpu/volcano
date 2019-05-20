# Adding roxygen comments
#' A function for generating an interactive plotly plot, and an ! interactive annotated gg volcano plot
#'
#' This function allows you to generate a publication ready volcano plot in png from a limma `topTable` output table
#' @param
#' @keywords interactiveVolcano
#' @export
#' @examples 
#' drawVolcano(toPlot, SAVEDIR = "my/fav/dir/" , FILENAME)



drawVolcano <- function(    toPLot                  =  toPLot , 
                            SAVEDIR                 =  SAVEDIR, 
                            FILENAME                = "interactiveVolcano",     
                            png_plot_width          =  1204   , 
                            png_plot_height         =  800    ,
                            png_plot_res            =  NA     , 
                            
                            FDR_correction          = "BH", 
                            FDR_threshold           = 0.05  ,
                            logFC_signif_cutoff     = 1.3, 

                            x_axis_ticks_count      =  10,

                            transparency            = 0.7, 

                            dashed_line_color       = "#990000",
                            dashed_line_colour      = "#990000",


                            xlab_custom             = 'logFoldChanges', 
                            ylab_custom             = '-log10(adjPValue)', 

                            ggtitle_custom          = paste0("Volcano plot for genes in {CLL} ExpressionSet")
                            ggsubtitle_custom       = paste0( nSignif, " features of ", nTotalFeatures, " found statistically significant, ( FDR =  ", FDR_threshold, ", ", FDR_correction, " adjustment)")
                            ggcaption_custom        = paste0("logFoldChange signif cutoff: "  , logFC_signif_cutoff, "\n",  "- log AdjPvalue signif cutoff: ", round(neg_logPvalue_threshold, 2)), 

                            ggtitle_custom_color    = ggtitle_custom_color   ,
                            ggsubtitle_custom_color = ggsubtitle_custom_color, 
                            ggcaption_custom_color  = ggcaption_custom_color , 

                            scale_color_manual_palette =  c(  "#0198F5"  , "#FA5E5B" ,  "#4A637B" )){


                    # Helper variables for dataset info
                    nSignif        <- dim(toPlot[toPlot$significant == TRUE, ])[1]
                    nTotalFeatures <- dim(toPlot)[1]

                    # Cutoff values for dashed lines
                    neg_logPvalue_threshold   <-  ((-1 * log10(FDR_threshold)))  

                    # Creating helper variables for an airy plot, wide and centered plot
                    x_absolut_extremum        <- max(abs(toPlot$logFC))
                    x_absolut_extremum_wider  <- x_absolut_extremum * 1.65

                    y_absolut_extremum        <- max(abs(toPlot$neg_logPvalues))
                    y_absolut_extremum_wider  <- y_absolut_extremum * 1.15


                    # Start building gglayers of plotness
                    p<-ggplot() +   theme_gray()+
                    
                    # Adding aes for scatter plot points
                    geom_point(data    = toPlot, 
                                mapping = aes(  x       = logFC, 
                                                y       = neg_logPvalues,
                                                shape   = significant,
                                                size    = significant,
                                                color   = expression,
                                                label   = featureID),
                                alpha = transparency)   +
                    
                    # Control legend block that shows first. Here we choose 1, aka top position
                    guides(color = guide_legend(order = 1)) +
                    
                    # Add custom hex triplette for higher, lower, no diff
                    scale_color_manual(values  = scale_color_manual_palette ) +
                    
                    # Add x,y axes custom labels
                    xlab(xlab_cutom)                +
                    ylab(ylab_custom)             +
                    
                    # Set explicitly number of x axis ticks:https://stackoverflow.com/questions/11335836/increase-number-of-axis-ticks
                    scale_x_continuous(breaks = pretty(toPlot$logFC         , n = x_axis_ticks_count)) +
                    
                    # Set ylim based on absolut y extremum for leaving some room to breathe in the top  
                    ylim(0, y_absolut_extremum_wider)    +
                    
                    # Setting annotation for title, subtitle, caption
                    labs(title      = ggtitle_custom,
                        subtitle   = ggsubtitle_custom,
                        caption    = ggcaption_custom)   +
                    
                    # Controlling font face, font colour and font size of labs. Also centering with hjust 
                    theme(
                        plot.title    = element_text(color = ggtitle_custom_color   , size = 16, face = "bold" , hjust = 0.5),
                        plot.subtitle = element_text(color = ggsubtitle_custom_color, size = 12                , hjust = 0.5),
                        plot.caption  = element_text(color = ggsubtitle_custom_color, face = "italic"                       ) ) +
                    
                    # Add on plot annotation 
                    # subset the data, to do so only for stat significant features (1-criterion: FDR)
                    geom_text_repel(data      = subset(toPlot, significant == TRUE ), 
                                    mapping   = (aes(x    = logFC,
                                                    y     = neg_logPvalues,
                                                    label = featureID))) +
                    
                    # Add dashed vertical and horizontal lines to denote doeble criterion significance areas
                    geom_hline(yintercept =  neg_logPvalue_threshold     ,  colour= dashed_line_colour, linetype="dashed") + 
                    geom_vline(xintercept =        logFC_signif_cutoff   ,  colour= dashed_line_colour, linetype="dashed") + 
                    geom_vline(xintercept =  (-1 * logFC_signif_cutoff ) ,  colour= dashed_line_colour, linetype="dashed") +
                    
                    
                    # DOWNregulated, low adjPvalue area
                    # Add shadowing for double-filtering criterion, both stat signif adj.Pvalue and stat signif UPregulated
                    annotate ("rect", 
                                xmin = logFC_signif_cutoff, 
                                xmax = Inf, 
                                ymin = neg_logPvalue_threshold, 
                                ymax = Inf, alpha = .2) +
                    
                    # DOWNregulated, low adjPvalue area
                    # Add another shadowing for double-filtering criterion, both stat signif adj.Pvalue and stat signif DOWNregulated
                    annotate ("rect", 
                                xmin = -logFC_signif_cutoff, 
                                xmax = -Inf, 
                                ymin = neg_logPvalue_threshold, 
                                ymax = Inf, alpha = .2)

                    # Save plotly object as .html file
                    htmlwidgets::saveWidget(ggplotly(p), file = paste0(SAVEDIR,"/",FILENAME, ".html"))

                    # Save ggplot object as .png file
                    png(filename =  paste0(FILENAME, ".png"),
                        width    = png_plot_width, 
                        height   = png_plot_height,
                        res      = res_plot_)
                    plot(p)
                    dev.off()        


                    return(p)
                    
                    # The End

                    }




