
#' Graphically Visualize MTLR Survival Curves
#'
#' Plot the survival curves returned from predict.mtlr. Users must have packages ggplot2 and reshape2 installed in order to use this function.
#' Survival curves for MTLR are smoothed using a monotonic cubic spline using a Hyman filtering between time points. For details regarding this
#' smoothing function see \code{\link[stats]{splinefun}}.
#'
#' @param curves survival curves formatted the same as those from predict.mtlr. Time points must be in the first column of the matrix followed by
#' columns representing survival probabilities for each observation.
#' @param index the index of the observation to plot. Here an index of 1 will refer to the second column of the curves object. If over 15 indices are
#' given the legend will be removed as to not take up plotting space. To avoid this behavior set remove_legend = FALSE.
#' @param color the color of the plotted survival curve. The length of color must match the length of index.
#' @param xlim the limits of the x-axis (must be a 2 length vector).
#' @param remove_legend if TRUE the legend will be removed if over 15 indices are supplied. If FALSE the legend will remain, however be aware that the
#' legend may take up lots of space.
#' @examples
#' #Set up the example:
#' library(survival)
#' mod <- mtlr(Surv(time,status)~., data = lung)
#' curves <- predict(mod, type = "survivalcurve")
#'
#' plotcurves(curves, 1:10)
#' plotcurves(curves, 1:3, color = c("red","blue","purple"))
#' plotcurves(curves, 1:10, xlim = c(0,42))
#'
#' #Note the legend is now gone:
#' plotcurves(curves, 1:20)
#'
#' #and it is back again:
#' plotcurves(curves, 1:20, remove_legend = FALSE)
#' @seealso \code{\link[MTLR]{mtlr}} \code{\link[MTLR]{predict.mtlr}}
#' @export
plotcurves <- function(curves, index = 1, color = c(), xlim = c(), remove_legend = TRUE){
  if (!requireNamespace(c("ggplot2","reshape2"), quietly = TRUE)) {
    stop("Package \"ggplot2\" and \"reshape2\" are needed for this function to work. Please install them.",
         call. = FALSE)
  }
  colorOK <- T
  if(!length(color))
    colorOK <- F
  else if(length(color) != length(index)){
    warning("If you would like to select custom colors please make sure the number of colors
            matches the number of curves.")
    colorOK <- F
  }
  time <- curves[,1]
  curves <- curves[,index +1,drop=F]
  plot_times <- seq(min(time),max(time), length.out = length(time)*100)
  plot_probs <- as.data.frame(sapply(curves,
                                   function(curve){
                                     curve <- ifelse(curve < 1e-20,0,curve)
                                     survivialSpline <- stats::splinefun(time, curve, method = "hyman")
                                     return(pmax(survivialSpline(plot_times),0))
                                   }
  ))
  plot_data <- cbind.data.frame(plot_times,plot_probs)
  plot_data <- reshape2::melt(plot_data,measure.vars = names(plot_data)[-1], variable.name = "Index")
  pl <- ggplot2::ggplot(data = plot_data, ggplot2::aes(x = plot_times,y = plot_data$value, colour = plot_data$Index))+
    ggplot2::geom_line(size = 1.5)
  if(colorOK)
    pl <- pl + ggplot2::scale_color_manual(values = color)
  if(length(xlim)==2){
    pl <- pl+ ggplot2::xlim(c(xlim[1],xlim[2]))
  }
  pl <- pl +
    ggplot2::scale_y_continuous( limits= c(0,1),breaks = c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1))+
    ggplot2::theme_bw() +
    ggplot2::theme(
          text = ggplot2::element_text(size=18, face=  "bold"),
          axis.title = ggplot2::element_text(size = 20),
          axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 15)),
          axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 15)))+
    ggplot2::labs(y = "Survival Probability",x = "Time", color = "Index" )
  #If we have too many survival curves then the legend takes up the whole plot. We have an if to catch this and remove the legend.
  if(length(index) > 15 & remove_legend){
    pl <- pl + ggplot2::theme(legend.position = "None")
  }
  return(pl)
}



