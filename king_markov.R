# King Markov and his island kingdom: the Metropolis algorithm
# Statistical Rethinking (2nd ed.), chapter 9, section 9.1 (R code 9.1)
#
# Ten islands sit in a ring. Island i has a population proportional to i.
# Each week the King flips a coin to pick a neighbouring island (the proposal)
# and moves there with probability proposal population / current population.
# In the long run he visits each island in proportion to its population.
#
# Run with:  Rscript king_markov.R
# Outputs:   figures/king_markov_figure.png   (static graph, like Figure 9.2)
#            figures/king_markov.gif          (animation)

library(ggplot2)
library(patchwork)
library(gifski)

set.seed(9)

num_islands <- 10
num_weeks <- 1e5

# ---- the Metropolis algorithm (R code 9.1) --------------------------------

positions <- rep(0, num_weeks)
proposals <- rep(0, num_weeks)
accepted <- rep(FALSE, num_weeks)
current <- 10
for (i in 1:num_weeks) {
  # record current position
  positions[i] <- current

  # flip coin to generate proposal
  proposal <- current + sample(c(-1, 1), size = 1)
  # now make sure he loops around the archipelago
  if (proposal < 1) proposal <- num_islands
  if (proposal > num_islands) proposal <- 1

  # move?
  prob_move <- proposal / current
  accepted[i] <- runif(1) < prob_move
  proposals[i] <- proposal
  current <- ifelse(accepted[i], proposal, current)
}

# ---- shared styling -------------------------------------------------------

ink <- "#0b0b0b"
ink_soft <- "#52514e"
muted <- "#898781"
grid_col <- "#e1e0d9"
surface <- "#fcfcfb"
blue <- "#2a78d6"
blue_light <- "#b7d3f6"
orange <- "#eb6834"

theme_markov <- theme_minimal(base_size = 13) +
  theme(
    plot.background = element_rect(fill = surface, colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = grid_col, linewidth = 0.3),
    axis.text = element_text(colour = muted),
    axis.title = element_text(colour = ink_soft),
    plot.title = element_text(colour = ink, face = "bold", size = 14),
    plot.subtitle = element_text(colour = ink_soft, size = 11),
    plot.title.position = "plot"
  )
theme_set(theme_markov)

target <- data.frame(
  island = 1:num_islands,
  share = (1:num_islands) / sum(1:num_islands)
)

# ---- static graph (Figure 9.2) --------------------------------------------

first_weeks <- data.frame(week = 1:100, island = positions[1:100])

p_trace <- ggplot(first_weeks, aes(week, island)) +
  geom_point(colour = blue, size = 2) +
  scale_y_continuous(breaks = 1:num_islands) +
  labs(
    title = "The King's first 100 weeks",
    subtitle = "It looks like aimless wandering...",
    x = "week", y = "island"
  )

visits <- data.frame(
  island = 1:num_islands,
  weeks = tabulate(positions, nbins = num_islands)
)

p_visits <- ggplot(visits, aes(island, weeks)) +
  geom_col(fill = blue, width = 0.35) +
  geom_point(
    data = target, aes(island, share * num_weeks),
    shape = 95, size = 9, colour = ink
  ) +
  scale_x_continuous(breaks = 1:num_islands) +
  scale_y_continuous(labels = scales::label_comma(), expand = expansion(c(0, 0.05))) +
  labs(
    title = "...but after 100,000 weeks",
    subtitle = "Time on each island matches its population (black ticks)",
    x = "island", y = "number of weeks"
  )

dir.create("figures", showWarnings = FALSE)
ggsave(
  "figures/king_markov_figure.png", p_trace + p_visits,
  width = 12, height = 5, dpi = 150, bg = surface
)

# ---- animation ------------------------------------------------------------

# islands on a ring, island 1 at the top, numbered clockwise
angle <- pi / 2 - 2 * pi * (0:(num_islands - 1)) / num_islands
islands <- data.frame(
  island = 1:num_islands,
  x = cos(angle),
  y = sin(angle)
)

# one frame per week at first, then speed up so the long run becomes visible
frame_weeks <- unique(c(1:60, round(exp(seq(log(62), log(5000), length.out = 90)))))

draw_frame <- function(week) {
  here <- islands[positions[week], ]
  there <- islands[proposals[week], ]
  moved <- accepted[week]
  move_col <- if (moved) blue else orange

  # stop the arrow short of both island circles
  arrow_df <- data.frame(
    x = here$x + 0.25 * (there$x - here$x),
    y = here$y + 0.25 * (there$y - here$y),
    xend = here$x + 0.75 * (there$x - here$x),
    yend = here$y + 0.75 * (there$y - here$y)
  )

  p_map <- ggplot(islands, aes(x, y)) +
    geom_point(aes(size = island), colour = blue_light) +
    geom_point(data = here, aes(size = island), colour = blue) +
    geom_segment(
      data = arrow_df, aes(xend = xend, yend = yend),
      colour = move_col, linewidth = 1.2,
      linetype = if (moved) "solid" else "22",
      arrow = arrow(length = unit(0.12, "inches"), type = "closed")
    ) +
    geom_text(aes(label = island), colour = ink, size = 4.5, fontface = "bold") +
    annotate("text", x = 0, y = 0.16, label = paste("Week", format(week, big.mark = ",")),
             colour = ink, size = 7, fontface = "bold") +
    annotate("text", x = 0, y = -0.06, colour = ink_soft, size = 4,
             label = sprintf("King on island %d, proposes %d", here$island, there$island)) +
    annotate("text", x = 0, y = -0.24, colour = move_col, size = 4.2, fontface = "bold",
             label = if (moved) {
               "accepted: he moves"
             } else {
               sprintf("rejected (p = %d/%d): he stays", there$island, here$island)
             }) +
    scale_size_area(max_size = 24, guide = "none") +
    coord_fixed(xlim = c(-1.3, 1.3), ylim = c(-1.3, 1.3)) +
    labs(
      title = "King Markov's island hopping",
      subtitle = "Circle area = island population. Dark circle = the King."
    ) +
    theme_void(base_size = 13) +
    theme(
      plot.title = element_text(colour = ink, face = "bold", size = 14),
      plot.subtitle = element_text(colour = ink_soft, size = 11),
      plot.title.position = "plot"
    )

  so_far <- data.frame(
    island = 1:num_islands,
    share = tabulate(positions[1:week], nbins = num_islands) / week
  )

  p_share <- ggplot(so_far, aes(island, share)) +
    geom_col(fill = blue, width = 0.35) +
    geom_point(data = target, shape = 95, size = 8, colour = ink) +
    scale_x_continuous(breaks = 1:num_islands) +
    scale_y_continuous(labels = scales::label_percent(), expand = expansion(c(0, 0))) +
    coord_cartesian(ylim = c(0, 0.4)) +
    labs(
      title = "Share of weeks spent on each island",
      subtitle = "Black ticks = share of the total population",
      x = "island", y = NULL
    )

  recent <- data.frame(week = max(1, week - 99):week)
  recent$island <- positions[recent$week]

  p_recent <- ggplot(recent, aes(week, island)) +
    geom_point(colour = blue, size = 1.6) +
    scale_y_continuous(breaks = 1:num_islands, limits = c(1, num_islands)) +
    scale_x_continuous(
      limits = c(max(1, week - 99), max(100, week)),
      labels = scales::label_comma()
    ) +
    labs(title = "Where he has been (last 100 weeks)", x = "week", y = "island")

  p_map + (p_share / p_recent) + plot_layout(widths = c(1, 1.1))
}

frame_dir <- file.path(tempdir(), "king_markov_frames")
dir.create(frame_dir, showWarnings = FALSE)
frame_files <- file.path(frame_dir, sprintf("frame_%04d.png", seq_along(frame_weeks)))

for (f in seq_along(frame_weeks)) {
  ggsave(
    frame_files[f], draw_frame(frame_weeks[f]),
    width = 12, height = 6.5, dpi = 80, bg = surface
  )
}

gifski(
  frame_files,
  gif_file = "figures/king_markov.gif",
  width = 960, height = 520, delay = 0.25
)
