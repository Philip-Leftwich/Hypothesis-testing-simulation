library(shiny)
library(bslib)
library(ggplot2)
library(purrr)
library(dplyr)
library(tidyr)

# UI definition
ui <- page_fluid(
  theme = bs_theme(bootswatch = "flatly"),
  title = "Understanding P-values, Effect Sizes, and Sample Sizes",
  
  card(
    card_header(
      h1("Statistical Inference Simulator"),
      p("This app demonstrates the relationship between p-values, effect sizes, and sample sizes.")
    ),
    
    card_body(
      p("Explore how p-values are distributed under null and alternative hypotheses, and how sample size affects statistical power.")
    )
  ),
  
  layout_sidebar(
    sidebar = sidebar(
      title = "Control Panel",
      
      # Simulation settings
      card(
        card_header("Simulation Parameters"),
        
        radioButtons("test_type", "Test type",
                     choices = c("Two-sample t-test" = "t_test",
                                 "Correlation test" = "cor_test")),
        
        radioButtons("hypothesis", "Simulate under",
                     choices = c("Null hypothesis (H₀)" = "null", 
                                 "Alternative hypothesis (H₁)" = "alternative")),
        
        conditionalPanel(
          condition = "input.test_type == 't_test' && input.hypothesis == 'alternative'",
          sliderInput("true_effect", "True effect size (Cohen's d)",
                      min = 0, max = 2, value = 0.5, step = 0.1)
        ),
        
        conditionalPanel(
          condition = "input.test_type == 'cor_test' && input.hypothesis == 'alternative'",
         sliderInput("true_cor", "True correlation (r)", 0, min = -1, max = 1, value = 0.5, step = 0.1)
        ),
        
        selectInput("sample_sizes", "Sample sizes to compare",
                    choices = list(
                      "Small to medium" = "small_to_medium",
                      "Medium to large" = "medium_to_large",
                      "Custom" = "custom"
                    )),
        
        conditionalPanel(
          condition = "input.sample_sizes == 'custom'",
          textInput("custom_sizes", "Enter sample sizes (comma-separated)", "10, 30, 50, 100")
        ),
        
        actionButton("run_simulation", "Run Simulation", class = "btn-primary")
      )
    ),
    
    # P-value distribution plot
    card(
      card_header("P-value Distribution"),
      card_body(
        plotOutput("p_value_hist", height = "300px")
      ),
      card_footer(
        "Under the null hypothesis, p-values should be uniformly distributed between 0 and 1."
      )
    ),
    
    # Effect size distribution plot
    card(
      card_header("Effect Size Distribution"),
      card_body(
        plotOutput("effect_size_hist", height = "300px")
      ),
      card_footer(
        "This shows the distribution of estimated effect sizes from the simulations."
      )
    ),
    
    # Power analysis plot
    card(
      card_header("Power Analysis"),
      card_body(
        plotOutput("power_plot", height = "300px")
      ),
      card_footer(
        "Power is the probability of rejecting the null hypothesis when the alternative is true."
      )
    ),
    
    # Summary and insights
    card(
      card_header("Key Insights"),
      card_body(
        h4("What to observe:"),
        tags$ul(
          tags$li("P-values under H₀ should be uniformly distributed"),
          tags$li("Effect sizes under H₀ are larger when the sample size is small - this is called the 'Winner's curse'"),
          tags$li("P-values under H₁ should cluster toward 0"),
          tags$li("Larger sample sizes lead to more power (ability to detect effects) under H₁"),
          tags$li("Larger effect sizes are easier to detect (require smaller samples)")
          
        ),
        textOutput("simulation_summary")
      )
    )
  )
)

# Server logic
server <- function(input, output, session) {
  
  # Helper: parse sample sizes
  get_sample_sizes <- reactive({
    if (input$sample_sizes == "small_to_medium") {
      c(10, 30, 50, 100)
    } else if (input$sample_sizes == "medium_to_large") {
      c(50, 100, 200, 500)
    } else {
      sizes <- as.numeric(strsplit(input$custom_sizes, ",")[[1]])
      sizes[!is.na(sizes)]
    }
  })
  
  n_sim <- 1000
  
  # Run simulations
  sim_results <- eventReactive(input$run_simulation, {
    
    sizes <- get_sample_sizes()
    results <- map_dfr(sizes, function(n) {
      p_vals <- numeric(n_sim)
      effs <- numeric(n_sim)
      for (j in seq_len(n_sim)) {
        if (input$test_type == "t_test") {
          if (input$hypothesis == "null") {
            g1 <- rnorm(n); g2 <- rnorm(n)
          } else {
            g1 <- rnorm(n); g2 <- rnorm(n, mean = input$true_effect)
          }
          tst <- t.test(g1, g2)
          p_vals[j] <- tst$p.value
          sd_pooled <- sqrt(((n-1)*var(g1)+(n-1)*var(g2))/(2*n-2))
          effs[j] <- (mean(g2)-mean(g1))/sd_pooled
        } else {
          if (input$hypothesis == "null") {
            x <- rnorm(n); y <- rnorm(n)
          } else {
            x <- rnorm(n)
            y <- input$true_cor * x + sqrt(1-input$true_cor^2)*rnorm(n)
          }
          tst <- cor.test(x, y)
          p_vals[j] <- tst$p.value
          effs[j] <- unname(tst$estimate)
        }
      }
      tibble(
        sample_size = n,
        p_value = p_vals,
        effect_size = effs,
        power = mean(p_vals < 0.05)
      )
    })
    results
  })
  
  # P-value histogram
  output$p_value_hist <- renderPlot({
    req(sim_results())
    df <- sim_results()
    p <- ggplot(df, aes(x = p_value,
                        y = ..density..)) +
      geom_histogram(bins = 30, fill = "steelblue", alpha = 0.6) +
      facet_wrap(~ sample_size) +
      labs(
        title = paste("P-value distribution under", if (input$hypothesis == "null") "H₀" else "H₁"),
        x = "p-value", y = "Density"
      ) +
      theme_minimal()+
      scale_x_continuous(limits = c(-0.05,1))

      p + geom_vline(xintercept = 0.05, linetype = "dashed", color = "red")
 
  })
  
  # Effect size histogram
  output$effect_size_hist <- renderPlot({
    req(sim_results())
    df <- sim_results()
    label <- if (input$test_type == "t_test") "Cohen's d" else "Correlation (r)"
    p1 <- ggplot(df, aes(x = effect_size,
                         y = ..density..)) +
      geom_histogram(bins = 30, fill = "darkgreen", alpha = 0.6) +
      facet_wrap(~ sample_size) +
      labs(title = paste(label, "distribution"), x = label, y = "Density") +
      theme_minimal()
  
    if (input$hypothesis == "null") p1 + geom_vline(xintercept = 0, linetype = "dashed", color = "red") else   p1 + 
      geom_vline(xintercept = if (input$test_type == "t_test") input$true_effect else input$true_cor,
                 linetype = "dashed", color = "red") 
  })
  
  # Power plot
  output$power_plot <- renderPlot({
    req(sim_results())
    power_df <- sim_results() %>% group_by(sample_size) %>% summarise(power = mean(p_value < 0.05))
    ggplot(power_df, aes(x = sample_size, y = power)) +
      geom_line(size = 1) +
      geom_point(size = 3) +
      geom_hline(yintercept = 0.8, linetype = "dashed") +
      scale_y_continuous(limits = c(0, 1)) +
      labs(
        title = "Power vs. Sample Size",
        subtitle = if (input$test_type == "t_test") paste0("Cohen's d = ", input$true_effect) else paste0("r = ", input$true_cor),
        x = "Sample Size", y = "Power"
      ) +
      theme_minimal()
  })
  
  # Simulation summary text
  output$simulation_summary <- renderText({
    req(sim_results())
    sum_df <- sim_results() %>% group_by(sample_size) %>%
      summarise(power = mean(p_value < 0.05), .groups = 'drop')
    lines <- paste0("n = ", sum_df$sample_size, ": power = ", sprintf("%.2f", sum_df$power))
    paste(lines, collapse = "; ")
  })
}

# Launch the app
shinyApp(ui, server)
