#' Create a single plot of the 32-run test case.

#' Decompose a number into powers of two
#' @details If log is TRUE, gives the log2(x)+1 of each power of two.
powers <- function(num, log = FALSE) {
  pows <- c()
  i <- 1L
  x <- as.integer(num)
  while (i <= x) {
    if (bitwAnd(i, x) > 0) {
      pows <- append(pows, i)
    }
    i <- bitwShiftL(i, 1)
  }
  if (log) {
    pows <- sapply(pows, function(x) log2(x) + 1)
  }
  return(pows)
}

# Packages
library(ggplot2)
library(dplyr)

# List of the three test cases
test_cases <- c(
  "results/N32_m1_r3.xlsx",
  "results/N64_m1_r4.xlsx",
  "results/N64_m2_r4.xlsx"
)

# Run code for each test case
for (case in test_cases) {
  # Find and load the data file
  data <- readxl::read_excel(case)
  
  # Find values of m and N
  test.case <- gsub("results/","",gsub(".xlsx","",case))
  N = as.integer(sub("N","",sub("_.*","",test.case)))
  m = as.integer(sub("m","",sapply(strsplit(test.case, "_"), "[", 2)))
  
  # Remove full-factorial designs
  min_n = log2(N) - 2*m
  tidy_data <- dplyr::filter(data, n > min_n)
  max_n = tidy_data %>% 
    group_by(method) %>% 
    summarise(n = max(n)) %>% 
    summarise(min(n)) %>% 
    unlist() %>%
    as.numeric()
  tidy_data <- dplyr::filter(tidy_data, n < max_n +1)
  
  # Plot time vs. $n$ for each method
  p <- ggplot(tidy_data,
              aes(x = n, y = time, linetype = method, shape = method)) +
    geom_point(size = 2) +
    geom_line() +
    labs(
      x = "Number of two-level factors",
      y = "Time (seconds)"
    ) +
    scale_linetype_discrete("Method") +
    scale_shape_discrete("Method") +
    scale_x_continuous(breaks = seq(min_n+1, max_n, 1)) +
    theme_bw(base_size = 15) +
    theme(
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      legend.position = "none"
	  )
    
  
  # Save the figure in writing folder to access in Latex
  # ggsave(paste0("Writing/figures/",N,"_run_m_",m,"_test_case.pdf"), p, dpi = 300,
  #        height = 5, width = 7)
  ggsave(paste0("figures/",N,"_run_m_",m,"_test_case.pdf"), p, dpi = 300,
         height = 5, width = 7)
}

