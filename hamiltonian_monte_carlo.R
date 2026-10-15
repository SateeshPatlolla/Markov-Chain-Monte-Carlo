# Hamiltonian Monte Carlo on a 2D Gaussian posterior
# Statistical Rethinking (2nd ed.), chapter 9, section 9.3
#
# The posterior is treated as a frictionless bowl: U is the negative
# log-posterior (the height of the bowl). For each sample we give a particle a
# random flick (momentum), let it glide for L leapfrog steps of size epsilon,
# and take the place where it stops as the proposal.
#
# Run with:  Rscript hamiltonian_monte_carlo.R
# Outputs:   figures/hmc_figure.png   (static graph, like Figure 9.6, top row)
#            figures/hmc.gif          (animation)

library(ggplot2)
library(patchwork)
library(gifski)

# ---- the model: x and y are both Normal(mu, 1), with Normal(0, 1) priors --

# test data
set.seed(7)
y <- rnorm(50)
x <- rnorm(50)
x <- as.numeric(scale(x))
y <- as.numeric(scale(y))

# U needs to return neg-log-probability
# (q is ordered mux, muy here so that q[1] is the horizontal axis)
U <- function(q, a = 0, b = 1, k = 0, d = 1) {
  mux <- q[1]
  muy <- q[2]
  U <- sum(dnorm(y, muy, 1, log = TRUE)) + sum(dnorm(x, mux, 1, log = TRUE)) +
    dnorm(muy, a, b, log = TRUE) + dnorm(mux, k, d, log = TRUE)
  return(-U)
}

# gradient function
# need vector of partial derivatives of U with respect to vector q
U_gradient <- function(q, a = 0, b = 1, k = 0, d = 1) {
  mux <- q[1]
  muy <- q[2]
  G1 <- sum(x - mux) + (k - mux) / d^2 # dU/dmux
  G2 <- sum(y - muy) + (a - muy) / b^2 # dU/dmuy
  return(c(-G1, -G2)) # negative bc energy is neg-log-prob
}

# ---- one HMC transition (the book's HMC2) -----------------------------------

HMC2 <- function(U, grad_U, epsilon, L, current_q) {
  q <- current_q
  p <- rnorm(length(q), 0, 1) # random flick - p is momentum
  current_p <- p
  # Make a half step for momentum at the beginning
  p <- p - epsilon * grad_U(q) / 2
  # initialize bookkeeping - saves trajectory
  qtraj <- matrix(NA, nrow = L + 1, ncol = length(q))
  ptraj <- qtraj
  qtraj[1, ] <- current_q
  ptraj[1, ] <- p

  # Alternate full steps for position and momentum
  for (i in 1:L) {
    q <- q + epsilon * p # Full step for the position
    # Make a full step for the momentum, except at end of trajectory
    if (i != L) {
      p <- p - epsilon * grad_U(q)
      ptraj[i + 1, ] <- p
    }
    qtraj[i + 1, ] <- q
  }

  # Make a half step for momentum at the end
  p <- p - epsilon * grad_U(q) / 2
  ptraj[L + 1, ] <- p
  # Negate momentum at end of trajectory to make the proposal symmetric
  p <- -p
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U <- U(current_q)
  current_K <- sum(current_p^2) / 2
  proposed_U <- U(q)
  proposed_K <- sum(p^2) / 2
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  accept <- 0
  if (runif(1) < exp(current_U - proposed_U + current_K - proposed_K)) {
    new_q <- q # accept
    accept <- 1
  } else {
    new_q <- current_q # reject
  }

  return(list(q = new_q, traj = qtraj, ptraj = ptraj, accept = accept, flick = current_p))
}

# run a chain, keeping every trajectory
run_hmc <- function(n_samples, epsilon, L, start) {
  chain <- vector("list", n_samples)
  q <- start
  for (i in 1:n_samples) {
    chain[[i]] <- HMC2(U, U_gradient, epsilon, L, q)
    q <- chain[[i]]$q
  }
  chain
}

# ---- shared styling -------------------------------------------------------

ink <- "#0b0b0b"
ink_soft <- "#52514e"
muted <- "#898781"
grid_col <- "#e1e0d9"
contour_col <- "#c3c2b7"
surface <- "#fcfcfb"
blue <- "#2a78d6"
blue_light <- "#86b6ef"
orange <- "#eb6834"

theme_set(
  theme_minimal(base_size = 13) +
    theme(
      plot.background = element_rect(fill = surface, colour = NA),
      panel.grid = element_blank(),
      axis.text = element_text(colour = muted),
      axis.title = element_text(colour = ink_soft),
      plot.title = element_text(colour = ink, face = "bold", size = 14),
      plot.subtitle = element_text(colour = ink_soft, size = 11),
      plot.title.position = "plot",
      plot.margin = margin(10, 20, 10, 10)
    )
)

# the posterior bowl, drawn as contours of the log-posterior
pr <- 0.5
bowl <- expand.grid(
  mux = seq(-pr, pr, length.out = 101),
  muy = seq(-pr, pr, length.out = 101)
)
bowl$U <- apply(bowl, 1, U)

base_plot <- function() {
  ggplot() +
    geom_contour(
      data = bowl, aes(mux, muy, z = U),
      colour = contour_col, linewidth = 0.3, bins = 9
    ) +
    coord_fixed(xlim = c(-pr, pr), ylim = c(-pr, pr), expand = FALSE) +
    labs(x = expression(mu[x]), y = expression(mu[y]))
}

# all trajectories of a chain as one data frame
paths_df <- function(chain) {
  do.call(rbind, lapply(seq_along(chain), function(i) {
    data.frame(
      sample = i,
      step = seq_len(nrow(chain[[i]]$traj)) - 1,
      mux = chain[[i]]$traj[, 1],
      muy = chain[[i]]$traj[, 2]
    )
  }))
}

# ---- static graph (Figure 9.6, top row) -----------------------------------

start <- c(-0.1, 0.2)
step <- 0.03

plot_chain <- function(L, title, subtitle) {
  set.seed(11)
  chain <- run_hmc(4, step, L, start)
  paths <- paths_df(chain)
  ends <- paths[paths$step == L, ]

  base_plot() +
    geom_path(data = paths, aes(mux, muy, group = sample), colour = blue, linewidth = 0.7) +
    geom_point(data = paths, aes(mux, muy), colour = blue, size = 1.3) +
    annotate("point", x = start[1], y = start[2], shape = 4, size = 3.5, stroke = 1.2, colour = ink) +
    geom_point(data = ends, aes(mux, muy), colour = ink, size = 6) +
    geom_text(data = ends, aes(mux, muy, label = sample), colour = surface, size = 3.5, fontface = "bold") +
    labs(title = title, subtitle = subtitle)
}

p_good <- plot_chain(
  11, "11 leapfrog steps",
  "Each path ends far from where it began"
)
p_uturn <- plot_chain(
  28, "28 leapfrog steps: the U-turn problem",
  "Paths loop back, so samples land near the start"
)

dir.create("figures", showWarnings = FALSE)
ggsave(
  "figures/hmc_figure.png",
  p_good + p_uturn + plot_annotation(
    caption = "x = starting point. Small dots = leapfrog steps (step size 0.03). Numbered dots = samples, in order.",
    theme = theme(plot.caption = element_text(colour = muted, size = 10, hjust = 0))
  ),
  width = 12, height = 6.4, dpi = 150, bg = surface
)

# ---- animation ------------------------------------------------------------

L <- 11
n_samples <- 1000
n_animated <- 20 # trajectories shown step by step before speeding up

set.seed(11)
chain <- run_hmc(n_samples, step, L, start)
paths <- paths_df(chain)
samples <- as.data.frame(rbind(start, t(sapply(chain, function(s) s$q))))
names(samples) <- c("mux", "muy") # row i + 1 is the position after sample i

# the particle part-way through trajectory i, after s leapfrog steps
draw_trajectory_frame <- function(i, s) {
  old_paths <- paths[paths$sample < i, ]
  path <- paths[paths$sample == i & paths$step <= s, ]
  particle <- path[nrow(path), ]
  origin <- samples[i, ]
  flick <- chain[[i]]$flick
  done <- s == L
  accepted <- chain[[i]]$accept == 1

  p <- base_plot() +
    geom_path(data = old_paths, aes(mux, muy, group = sample), colour = grid_col, linewidth = 0.5) +
    geom_point(data = samples[seq_len(i), ], aes(mux, muy), colour = blue_light, size = 2.5) +
    geom_point(data = path, aes(mux, muy), colour = blue, size = 1.5)

  if (s > 0) {
    p <- p + geom_path(data = path, aes(mux, muy), colour = blue, linewidth = 0.9)
  } else {
    # the random flick that sets the particle going
    p <- p + annotate(
      "segment",
      x = origin$mux, y = origin$muy,
      xend = origin$mux + 0.07 * flick[1], yend = origin$muy + 0.07 * flick[2],
      colour = ink, linewidth = 0.9,
      arrow = arrow(length = unit(0.1, "inches"), type = "closed")
    )
  }

  status <- if (s == 0) {
    "a random flick"
  } else if (!done) {
    sprintf("leapfrog step %d of %d", s, L)
  } else if (accepted) {
    "accepted"
  } else {
    "rejected: stay at the start"
  }

  p +
    geom_point(
      data = particle, aes(mux, muy),
      colour = if (done && !accepted) orange else ink, size = 4.5
    ) +
    labs(
      title = sprintf("Hamiltonian Monte Carlo: sample %d", i),
      subtitle = paste0("Grey rings = the posterior bowl. Now: ", status)
    )
}

# later on: just the samples piling up
draw_samples_frame <- function(n) {
  base_plot() +
    geom_point(
      data = samples[seq_len(n + 1), ], aes(mux, muy),
      colour = blue, alpha = 0.45, size = 2
    ) +
    labs(
      title = sprintf("Hamiltonian Monte Carlo: %s samples", format(n, big.mark = ",")),
      subtitle = "Skipping the paths: the samples fill in the posterior"
    )
}

later <- unique(round(exp(seq(log(n_animated + 2), log(n_samples), length.out = 50))))

frames <- c(
  unlist(lapply(seq_len(n_animated), function(i) {
    lapply(0:L, function(s) function() draw_trajectory_frame(i, s))
  })),
  lapply(later, function(n) function() draw_samples_frame(n)),
  # hold on the final picture
  rep(list(function() draw_samples_frame(n_samples)), 10)
)

frame_dir <- file.path(tempdir(), "hmc_frames")
dir.create(frame_dir, showWarnings = FALSE)
frame_files <- file.path(frame_dir, sprintf("frame_%04d.png", seq_along(frames)))

for (f in seq_along(frames)) {
  ggsave(frame_files[f], frames[[f]](), width = 7, height = 7.4, dpi = 90, bg = surface)
}

gifski(
  frame_files,
  gif_file = "figures/hmc.gif",
  width = 630, height = 666, delay = 0.1, progress = FALSE
)
