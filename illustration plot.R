my_theme_large <- theme_bw(base_size = 18) +
  theme(
    text = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(size = 19),
    plot.subtitle = element_text(size = 15),
    legend.title = element_text(size = 17),
    legend.text = element_text(size = 16)
  )

library(ggplot2)
library(grid)

pts <- data.frame(
  x = c(0, 4, 5, 2, 0, 1, 4),
  y = c(0, 4, 3, 0, 2, 3, 0)
)

seg <- data.frame(
  x = pts$x[-nrow(pts)],
  y = pts$y[-nrow(pts)],
  xend = pts$x[-1],
  yend = pts$y[-1],
  segment = factor(seq_len(nrow(pts) - 1))
)

events <- pts[2:nrow(pts), ]
events$event_type <- "Reflection"

start <- pts[1, ]
start$event_type <- "Initial point"

zigzag <- ggplot() +
  geom_segment(
    data = seg,
    aes(x = x, y = y, xend = xend, yend = yend, colour = segment),
    linewidth = 1.3,
    arrow = arrow(length = unit(0.3, "cm"), type = "closed")
  ) +
  geom_point(
    data = events,
    aes(x, y, shape = event_type),
    fill = "white",
    colour = "black",
    size = 3,
    stroke = 1.1
  ) +
  geom_point(
    data = start,
    aes(x, y, shape = event_type),
    colour = "black",
    fill = "black",
    size = 3,
    stroke = 1.1
  ) +
  scale_shape_manual(
    values = c(
      "Initial point" = 16,
      "Reflection" = 21
    )
  ) +
  coord_equal(xlim = c(-2, 6), ylim = c(-2, 6), expand = FALSE) +
  theme_classic(base_size = 14) +
  labs(
    title = "ZigZag sampler",
    x = expression(theta[1]),
    y = expression(theta[2]),
    colour = "Trajectory segment",
    shape = "Event type"
  ) + 
  my_theme_large
zigzag






library(ggplot2)
library(grid)

pts <- data.frame(
  x = c(0, 4, 0.2, 2.1, 3, 1, 0.5),
  y = c(0, -0.5, 3.5, 4, -1, 2, -0.5)
)

seg <- data.frame(
  x = pts$x[-nrow(pts)],
  y = pts$y[-nrow(pts)],
  xend = pts$x[-1],
  yend = pts$y[-1],
  segment = factor(seq_len(nrow(pts) - 1))
)

events <- pts[2:nrow(pts), ]

set.seed(123)
refresh_id <- sample(seq_len(nrow(events)), size = 3)

events$event_type <- "Reflection"
events$event_type[refresh_id] <- "Refreshment"

start <- pts[1, ]
start$event_type <- "Initial point"

bps <- ggplot() +
  geom_segment(
    data = seg,
    aes(x = x, y = y, xend = xend, yend = yend, colour = segment),
    linewidth = 1.3,
    arrow = arrow(length = unit(0.3, "cm"), type = "closed")
  ) +
  geom_point(
    data = events,
    aes(x, y, shape = event_type),
    fill = "white",
    colour = "black",
    size = 3,
    stroke = 1.1
  ) +
  geom_point(
    data = start,
    aes(x, y, shape = event_type),
    colour = "black",
    fill = "black",
    size = 3,
    stroke = 1.1
  ) +
  scale_shape_manual(
    values = c(
      "Initial point" = 16,
      "Reflection" = 21,
      "Refreshment" = 24
    )
  ) +
  coord_equal(xlim = c(-2, 6), ylim = c(-2, 6), expand = FALSE) +
  theme_classic(base_size = 14) +
  labs(
    title = "BPS sampler illustration",
    x = expression(theta[1]),
    y = expression(theta[2]),
    colour = "Trajectory segment",
    shape = "Event type"
  ) +
  my_theme_large
bps







library(ggplot2)
library(grid)
library(dplyr)

make_orbit_arc <- function(p0, p1, center, a = 1, b = 1.5,
                           direction = "ccw", n = 100) {
  
  u0 <- (p0[1] - center[1]) / a
  v0 <- (p0[2] - center[2]) / b
  u1 <- (p1[1] - center[1]) / a
  v1 <- (p1[2] - center[2]) / b
  
  th0 <- atan2(v0, u0)
  th1 <- atan2(v1, u1)
  
  if (direction == "ccw" && th1 < th0) th1 <- th1 + 2*pi
  if (direction == "cw"  && th1 > th0) th1 <- th1 - 2*pi
  
  r0 <- sqrt(u0^2 + v0^2)
  r1 <- sqrt(u1^2 + v1^2)
  
  s <- seq(0, 1, length.out = n)
  th <- (1 - s) * th0 + s * th1
  r  <- (1 - s) * r0  + s * r1
  
  x <- center[1] + a * r * cos(th)
  y <- center[2] + b * r * sin(th)
  
  data.frame(x = x, y = y)
}

# 设定 event points
pts <- data.frame(
  x = c(0, 4, 0.5, 2, 3, 0.8, 1.9),
  y = c(0, 3, 3.5, 5.5, 3.5, 0.7, 1.9)
)

# 每一段各自的 center
centers <- list(
  c(2.0, -1.8),  # 1: 从(0,0)到(4,3)，向上拱，比较弯
  c(1.8, 5.2),   # 2: 从(4,3)到(-0.5,3.8)，向下弯，中等偏弯
  c(1.5, 2),   # 3: 从(-0.5,3.8)到(2,6)，向下弯，超大弧度
  c(1.5, 1.5),   # 4: 从(2,6)到(3,3.5)，比较平缓
  c(1.9, 2.1),   # 5: 从(3,3.5)到(0.8,0.7)，正常偏弯
  c(1.5, 1.1)    # 6: 从(0.8,0.7)到(3.6,2)，随机
)

# 每一段方向
dirs <- c(
  "cw",  # 1: 向上拱
  "cw",   # 2: 向下弯
  "ccw",   # 3: 向下弯，超大弧度
  "cw",  # 4: 平缓一点
  "cw",  # 5: 正常弯度
  "cw"    # 6: 随机
)

# 每一段弧度参数
a_vals <- c(
  1,  # 1
  1.0,  # 2
  0.5,  # 3
  1,  # 4
  0.6,  # 5
  1   # 6
)

b_vals <- c(
  2.2,  # 1: 比较弯
  1.2,  # 2: 中等弯
  1.8,  # 3: 超级大的弧度
  0.2,  # 4: 比较平缓
  1.3,  # 5: 正常偏弯
  1   # 6: 随机
)

curve_list <- list()

for (i in 1:(nrow(pts) - 1)) {
  curve_list[[i]] <- make_orbit_arc(
    p0 = as.numeric(pts[i, ]),
    p1 = as.numeric(pts[i + 1, ]),
    center = centers[[i]],
    a = a_vals[i],
    b = b_vals[i],
    direction = dirs[i],
    n = 100
  ) %>%
    mutate(segment = factor(i))
}

curve_data <- bind_rows(curve_list)

events <- pts[2:nrow(pts), ]
events$event_type <- c("Refreshment", "Reflection", "Reflection",
                       "Refreshment", "Reflection", "Reflection")

start <- pts[1, ]
start$event_type <- "Initial point"

reference <- data.frame(x = 1, y = 0)

Boomerang <- ggplot() +
  geom_path(
    data = curve_data,
    aes(x = x, y = y, colour = segment, group = segment),
    linewidth = 1.2,
    arrow = arrow(length = unit(0.3, "cm"), type = "closed")
  ) +
  geom_point(
    data = events,
    aes(x = x, y = y, shape = event_type),
    fill = "white",
    colour = "black",
    size = 3,
    stroke = 1.1
  ) +
  geom_point(
    data = start,
    aes(x = x, y = y, shape = event_type),
    colour = "black",
    fill = "black",
    size = 3
  ) +
  scale_shape_manual(
    values = c(
      "Initial point" = 16,
      "Reflection" = 21,
      "Refreshment" = 24
    )
  ) +
  coord_equal(xlim = c(-2, 6), ylim = c(-2, 6), expand = FALSE) +
  theme_classic(base_size = 14) +
  labs(
    title = "Boomerang sampler",
    x = expression(theta[1]),
    y = expression(theta[2]),
    colour = "Trajectory segment",
    shape = "Event type"
  ) + 
  my_theme_large
Boomerang






library(ggplot2)
library(grid)

# ---------------------------
# Trajectory points for a NUTS/HMC illustration
# ---------------------------
traj <- data.frame(
  x = c(0.8, 1.2, 1.7, 2.3, 3.0, 3.8, 4.5, 5.0),
  y = c(0.9, 1.6, 2.4, 3.1, 3.8, 4.2, 4.4, 4.3)
)

# segments between consecutive leapfrog steps
seg <- data.frame(
  x    = traj$x[-nrow(traj)],
  y    = traj$y[-nrow(traj)],
  xend = traj$x[-1],
  yend = traj$y[-1]
)

# initial point
start <- traj[1, ]
start$point_type <- "Initial point"

# proposal point
proposal <- traj[nrow(traj), ]
proposal$point_type <- "Proposal point"

# intermediate leapfrog steps
mid_points <- traj[2:(nrow(traj)-1), ]
mid_points$point_type <- "Leapfrog step"

# combine for legend
point_data <- rbind(start, proposal, mid_points)

nuts <- ggplot() +
  # trajectory made of many small arrows
  geom_segment(
    data = seg,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 1.1,
    colour = "steelblue4",
    arrow = arrow(length = unit(0.3, "cm"), type = "closed")
  ) +
  
  # intermediate points
  geom_point(
    data = mid_points,
    aes(x = x, y = y, shape = point_type),
    colour = "black",
    fill = "white",
    size = 2.8,
    stroke = 1
  ) +
  
  # initial point
  geom_point(
    data = start,
    aes(x = x, y = y, shape = point_type),
    colour = "black",
    fill = "black",
    size = 3,
    stroke = 1
  ) +
  
  # proposal point
  geom_point(
    data = proposal,
    aes(x = x, y = y, shape = point_type),
    colour = "black",
    fill = "white",
    size = 3,
    stroke = 1.2
  ) +
  
  scale_shape_manual(
    values = c(
      "Initial point" = 16,
      "Leapfrog step" = 21,
      "Proposal point" = 24
    )
  ) +
  coord_equal(xlim = c(0, 6), ylim = c(0, 6), expand = FALSE) +
  theme_classic(base_size = 14) +
  labs(
    title = "NUTS sampler illustration",
    x = expression(theta[1]),
    y = expression(theta[2]),
    shape = "Point type"
  ) + 
  my_theme_large
nuts


ggsave(
  filename = "ZigZag.png",
  plot = zigzag,
  width = 7,
  height = 6,
  dpi = 300
)

ggsave(
  filename = "BPS.png",
  plot = bps,
  width = 7,
  height = 6,
  dpi = 300
)

ggsave(
  filename = "Boomerang.png",
  plot = Boomerang,
  width = 7,
  height = 6,
  dpi = 300
)

ggsave(
  filename = "NUTS.png",
  plot = nuts,
  width = 7,
  height = 6,
  dpi = 300
)
