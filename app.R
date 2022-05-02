###############
# Autor: Valentin Koob
# Description: Shows the influence of assumption violations
# on a one-way ANOVA
###############


######### Ab hier nur noch Berechnungen und Plots #########
plotPopulation = function(sampleInput,
                          sdInput1,
                          sdInput2, 
                          sdInput3) {
  # find the min and max values of the population
  if (sampleInput == "norm") {
    xmin = floor(-max(sdInput1, sdInput2, sdInput3) * 4)
    xmax = -xmin
  } else if (sampleInput == "beta") {
    # varianz bei alpha = 2 = beta -> 0.05
    scale_1 = sqrt(sdInput1 ^ 2 / 0.05)
    scale_2 = sqrt(sdInput2 ^ 2 / 0.05)
    scale_3 = sqrt(sdInput3 ^ 2 / 0.05)
    xmax = ceiling(max(scale_1, scale_2, scale_3) / 2) # da spaeter auf 0 zentriert
    xmin = -xmax
  } else if (sampleInput == "unif") {
    # Var(unif) = 1/12 [0,1]
    scale_1 = sqrt(sdInput1 ^ 2 / (1 / 12))
    scale_2 = sqrt(sdInput2 ^ 2 / (1 / 12))
    scale_3 = sqrt(sdInput3 ^ 2 / (1 / 12))
    xmax = ceiling(max(scale_1, scale_2, scale_3) / 2) # da spaeter auf 0 zentriert
    xmin = -xmax
  }
  
  # calculate density distributions
  xs = seq(xmin, xmax, length.out = 10000)
  if (sampleInput == "norm") {
    dichte1 = dnorm(xs, mean = 0, sd = sdInput1)
    dichte2 = dnorm(xs, mean = 0, sd = sdInput2)
    dichte3 = dnorm(xs, mean = 0, sd = sdInput3)
    
  } else if (sampleInput == "beta") {
    dichte1 = dbeta(seq(0, 1, length.out = 10000), 2, 2) / (scale_1)
    dichte2 = dbeta(seq(0, 1, length.out = 10000), 2, 2) / (scale_2)
    dichte3 = dbeta(seq(0, 1, length.out = 10000), 2, 2) / (scale_3)
    
  } else if (sampleInput == "unif") {
    dichte1 = dunif(xs, min = -scale_1 / 2, max = scale_1 / 2)
    dichte2 = dunif(xs, min = -scale_2 / 2, max = scale_2 / 2)
    dichte3 = dunif(xs, min = -scale_3 / 2, max = scale_3 / 2)
    
  }
  
  # plot the distributions
  if (sampleInput == "norm" | sampleInput == "unif") {
    plot(
      c(1, 2) ~ c(1, 1),
      xlim = c(xmin, xmax),
      ty = "l",
      col = "white",
      xlab = "moegliche Werte x",
      ylab = "f(x)",
      ylim = c(0, max(c(
        dichte1, dichte2, dichte3
      )) + max(c(
        dichte1, dichte2, dichte3
      )) * 0.35),
      cex.lab = 1.25,
      cex.axis = 1.25
    )
    polygon(x = xs,
            y = dichte1,
            col = rgb(0, 1, 0, 0.2))
    polygon(x = xs,
            y = dichte2,
            col = rgb(1, 0, 0, 0.2))
    polygon(x = xs,
            y = dichte3,
            col = rgb(0, 0, 1, 0.2))
  } else if (sampleInput == "beta") {
    plot(
      c(1, 2) ~ c(1, 1),
      xlim = c(xmin, xmax),
      ty = "l",
      col = "white",
      xlab = "moegliche Werte x",
      ylab = "f(x)",
      ylim = c(0, max(c(
        dichte1, dichte2, dichte3
      )) + max(c(
        dichte1, dichte2, dichte3
      )) * 0.35)
    )
    polygon(
      x = seq(0 - scale_1 / 2, scale_1 / 2, length.out = 10000),
      y = dichte1,
      col = rgb(0, 1, 0, 0.2)
    )
    polygon(
      x = seq(0 - scale_2 / 2, scale_2 / 2, length.out = 10000),
      y = dichte2,
      col = rgb(1, 0, 0, 0.2)
    )
    polygon(
      x = seq(0 - scale_3 / 2, scale_3 / 2, length.out = 10000),
      y = dichte3,
      col = rgb(0, 0, 1, 0.2)
    )
    
  }
  
  legend(
    "topright",
    legend = c("Population 1", "Population 2", "Population 3"),
    col = c(rgb(0, 1, 0, 0.8), rgb(1, 0, 0, 0.8), rgb(0, 0, 1, 0.8)),
    pch = 15,
    cex = 1.25
  )
  
  
}

plotSimulation = function(sampleInput,
                          sdInput1,
                          sdInput2,
                          sdInput3,
                          n,
                          nSim,
                          buttonState) {
  if (!sampleInput %in% c("norm", "beta", "unif"))
    stop("NOT DEFINED DISTRIBUTION")
  
  
  buttonState = as.numeric(buttonState)
  
  # if buttonState > 0  run the simulation
  if (buttonState > 0) {
    # run the simulation
    res =
      lapply(1:nSim, function(oneRep) {
        # get the samples in one repetition
        if (sampleInput == "norm") {
          sample1 = rnorm(n, mean = 0, sd = sdInput1)
          sample2 = rnorm(n, mean = 0, sd = sdInput2)
          sample3 = rnorm(n, mean = 0, sd = sdInput3)
        } else if (sampleInput == "beta") {
          scale_1 = sqrt(sdInput1 ^ 2 / 0.05) # Streckung der Daten
          scale_2 = sqrt(sdInput2 ^ 2 / 0.05) # Streckung der Daten
          scale_3 = sqrt(sdInput3 ^ 2 / 0.05) # Streckung der Daten
          sample1 = (rbeta(n, 2, 2) - 0.5) * scale_1
          sample2 = (rbeta(n, 2, 2) - 0.5) * scale_2
          sample3 = (rbeta(n, 2, 2) - 0.5) * scale_3
          
        } else if (sampleInput == "unif") {
          scale_1 = sqrt(sdInput1 ^ 2 / (1 / 12)) # obere Grenze finden
          scale_2 = sqrt(sdInput2 ^ 2 / (1 / 12)) # obere Grenze finden
          scale_3 = sqrt(sdInput3 ^ 2 / (1 / 12)) # obere Grenze finden
          
          sample1 = runif(n, min = -scale_1 / 2, max = scale_1 / 2)
          sample2 = runif(n, min = -scale_2 / 2, max = scale_2 / 2)
          sample3 = runif(n, min = -scale_3 / 2, max = scale_3 / 2)
          
        }
        
        # Stichproben als long-format abspeichern
        sample <- c(sample1, sample2, sample3)
        
        indices <- c(rep(1, n), rep(2, n), rep(3, n))
        
        df <- data.frame(sample = sample, indices = indices)
        oneway_ergebnis <- oneway.test(sample ~ indices,
                                       data = df,
                                       var.equal = TRUE)
        F_val <- oneway_ergebnis$statistic # F wert
        p <- oneway_ergebnis$p.value # p wert
        return(matrix(c(F_val, p), nrow = 1, ncol = 2)) # Rueckgabe
      })
    res = do.call("rbind", res)
    
    #### plotten der theoretischen F-Verteilung
    par(mar = c(5.1, 4.1, 4.1, 2.1))
    dfs_1 = 3-1 
    dfs_2 = (3 * n) - 3
    xmin = 0
    xmax = 12
    xs = seq(xmin, xmax, length.out = 10000)
    dichteF = df(xs, df1 = dfs_1, df2 = dfs_2)
    histVal = hist(res[, 1][res[, 1] <= xmax], breaks = 50, plot = FALSE)
    plot(
      c(1, 2) ~ c(1, 1),
      col = "white",
      xlim = c(xmin, xmax),
      xlab = "moegliche F-Werte",
      ylab = "f(F)",
      ylim = c(0, max(dichteF) + max(dichteF) * 0.35),
      main = "",
      cex.lab = 1.25,
      cex.axis = 1.25
    )
    polygon(x = c(0,0, histVal$mids, xmax), y = c(0,  histVal$density[1], histVal$density,0), col = "skyblue")
    points(dichteF ~ xs,
           col = "red",
           ty = "l",
           lwd = 2)
    legend(
      "topright",
      legend = c(
        paste("theoretische F-Verteilung\ndf1 =", dfs_1, "und df2 =", dfs_2)
        ,
        "simulierte F-Werte"
      ),
      col = c("red", "skyblue"),
      lty = 1,
      cex = 1.25
    )
    
    nomAlpha = mean(res[, 2] <= 0.05)
    mtext(
      side = 3,
      paste("% Fehlentscheidungen\n", round(nomAlpha, 4) * 100),
      cex = 1.5,
      line = 0.5
    )
    F_krit = qf(0.95, dfs_1 , dfs_2)
    segments(
      x0 = F_krit,
      y0 = 0,
      x1 = F_krit,
      y1 = 0.1
    )
    text(x = F_krit,
         y = 0.15,
         expression(F[krit]),
         cex = 1.25)
  } else{
    # place holder when no simulation was yet run
    par(mar = c(0, 0, 0, 0))
    plot(
      c(0, 1),
      c(0, 1),
      ann = F,
      bty = 'n',
      type = 'n',
      xaxt = 'n',
      yaxt = 'n'
    )
    text(
      x = 0.5,
      y = 0.5,
      paste("Bitte Simulation starten :)"),
      cex = 1.6,
      col = "black"
    )
    
  }
}



################## Shiny Implementierung


ui <- fluidPage(
  titlePanel(
    "Auswirkung der Annahmensverletzungen auf das nominelle Alpha-Niveau (bei der einfaktoriellen between-subjects ANOVA)"
  ),
  sidebarLayout(
    sidebarPanel(
      sliderInput(
        "SdInput1",
        h4("Standardabweichung der 1. Population"),
        min = 5,
        max = 50,
        value = 5,
        step = 5
      ),
      sliderInput(
        "SdInput2",
        h4("Standardabweichung der 2. Population"),
        min = 5,
        max = 50,
        value = 5,
        step = 5
      ),
      sliderInput(
        "SdInput3",
        h4("Standardabweichung der 3. Population"),
        min = 5,
        max = 50,
        value = 5,
        step = 5
      ),
      sliderInput(
        "nsInput",
        h4("Stichprobengroesse (pro Gruppe)"),
        min = 3,
        max = 15,
        value = 7,
        step = 1
      ),
      selectInput(
        "SampleInput",
        h4("Woraus soll gezogen werden?"),
        list(
          Normalverteilung = "norm",
          `skalierte Beta-Verteilung` = "beta",
          Gleichverteilung = "unif"
        ),
        "Normalverteilung"
      ),
      h4("Simulation starten"),
      actionButton("buttonInput",
                   "go!")
      
      
    ),
    
    mainPanel(fluidPage(
      h4(
        "Veranschaulichung der Populationen aus denen gezogen wird (unter der H0)"
      ),
      fluidRow(plotOutput(outputId = "population")),
      
      h4("Theoretische vs. simulierte Verteilung der F-Werte"),
      fluidRow(plotOutput(outputId = "sim")),
      
    ))
  )
)




server <- function(input, output) {
  globalButton = 0
  output$population <- renderPlot({
    plotPopulation(input$SampleInput,
                   input$SdInput1,
                   input$SdInput2, 
                   input$SdInput3)
  })
  
  output$sim <- renderPlot({
    plotSimulation(
      input$SampleInput,
      input$SdInput1,
      input$SdInput2,
      input$SdInput3,
      input$nsInput,
      5000,
      input$buttonInput
    )
  })
  
}

shinyApp(ui = ui, server = server)
