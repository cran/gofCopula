# Internal functions for the derivation of the amh copula for the PIOSRn and
# PIOSTn test statistics. Called by .Rn and .Tn, the internal functions for
# their test statistics.

# 2d density
.amh.12.density <- function(theta, u) { 
  if (is.null(dim(u))) {
    u1 <- u[1]
    u2 <- u[2]
  } else {
    u1 <- u[, 1]
    u2 <- u[, 2]
  }

  -(1 + theta * (-2 + u1 + theta * (-1 + u1) * (-1 + u2) + u2 + u1 * u2)) / 
    ((-1 + theta * (-1 + u1) * (-1 + u2))^3)
}

# 2d first derivative wrt theta of the log-density
.amh.12.V <- function(theta, u1, u2) { 
  (-1 + u1 * (2 - 4 * u2) - (theta^2) * ((-1 + u1)^2) * 
     ((-1 + u2)^2) + 2 * u2 - 2 * theta * (-1 + u1) * (-1 + u2) * 
     (-1 + u1 + u2 + u1 * u2)) / ((-1 + theta * (-1 + u1) * (-1 + u2)) * 
                                    (1 + theta * (-2 + u1 + theta * (-1 + u1) * 
                                                    (-1 + u2) + u2 + u1 * u2)))
}

# 2d second derivative wrt theta of the log-density
.amh.12.S <- function(theta, u1, u2) { 
  (1 + (theta^4) * ((-1 + u1)^4) * ((-1 + u2)^4) + 2 * (-2 + u2) * u2 +
    4 * (theta^3) * ((-1 + u1)^3) * ((-1 + u2)^3) * (-1 + u1 + u2 + u1 * u2) + 
     2 * (u1^2) * (1 + (-4 + u2) * u2) - 4 * u1 * (1 + 2 * (-2 + u2) * u2) +
    2 * (theta^2) * ((-1 + u1)^2) * ((-1 + u2)^2) * 
     (3 + (-6 + u2) * u2 + (u1^2) * ((1 + u2)^2) + 
        2 * u1 * (-1 + u2) * (3 + u2)) + 4 * theta * (-1 + u1) * (-1 + u2) * 
     (-1 - (-3 + u2) * u2 + u1 * (3 + (-7 + u2) * u2) + (u1^2) * 
        (-1 + u2 + 2 * (u2^2)))) /
    (((-1 + theta * (-1 + u1) * (-1 + u2))^2) * 
       ((1 + theta * (-2 + u1 + theta * (-1 + u1) * 
                        (-1 + u2) + u2 + u1 * u2))^2))
}

# 3d density
.amh.123.density <- function(theta, u) { 
  if (is.null(dim(u))) {
    u1 <- u[1]
    u2 <- u[2]
    u3 <- u[3]
  } else {
    u1 <- u[, 1]
    u2 <- u[, 2]
    u3 <- u[, 3]
  }

  (1 + (theta^6) * ((-1 + u1)^2) * ((-1 + u2)^2) * 
      ((-1 + u3)^2) + 2 * theta * (-3 + u1 + u2 + u3 + 2 * u1 * u2 * u3) +
    2 * (theta^5) * (-1 + u1) * (-1 + u2) * (-1 + u3) * 
      (3 + u2 * (-2 + u3) - 2 * u3 + u1 * (-2 + u2 + u3)) +
    (theta^2) * (-10 * u1 + (u1^2) + 4 * u1 * u2 + (u2^2) + 4 * 
    (u1 + u2 + u1 * u2 * (-3 + u1 + u2)) * u3 +
    (1 + u1 * u2 * (4 + u1 * u2)) * (u3^2) - 5 * (-3 + 2 * u2 + 2 * u3)) + 
      (theta^4) * (15 - 20 * u3 + 6 * (u3^2) + (u2^2) * (6 + (-6 + u3) * u3) +
      u2 * (-20 - 6 * (-4 + u3) * u3) + (u1^2) * (6 + (u2^2) * 
      ((1 - 2 * u3)^2) + (-6 + u3) * u3 + u2 * (-6 - 4 * (-2 + u3) * u3)) +
      u1 * (-20 - 6 * (-4 + u3) * u3 + 4 * u2 * (-2 + u3) * (-3 + 2 * u3) + 
      (u2^2) * (-6 - 4 * (-2 + u3) * u3))) + 2 * (theta^3) * ((u2^2) * 
      (-2 + u3) + u2 * (-8 + u3) * u3 - 2 * (u3^2) + 10 * (-1 + u2 + u3) + 
      (u1^2) * (-2 + u2 + u3) * (1 + 2 * u2 * u3) + u1 * (10 + (-8 + u3) * 
      u3 + (u2^2) * (1 + 2 * (-2 + u3) * u3) - 2 * u2 * (4 + u3 * 
      (-5 + 2 * u3))))) / ((-1 - theta * (-2 + u1 + u2) + (theta^2) * 
      (-1 + u1) * (-1 + u2) * (-1 + u3) + theta * (-1 + u1 * u2) * u3)^4)
}

# 3d first derivative wrt theta of the log-density
.amh.123.V <- function(theta, u) { 
  if (is.null(dim(u))) {
    u1 <- u[1]
    u2 <- u[2]
    u3 <- u[3]
  } else {
    u1 <- u[, 1]
    u2 <- u[, 2]
    u3 <- u[, 3]
  }

  (-2 * (1 - u2 + (theta^7) * ((u1 - 1)^3) * ((u2 - 1)^3) * 
           ((u3 - 1)^3) - u3 - (theta^6) * ((u1 - 1)^2) * ((u2 - 1)^2) * 
           ((u3 - 1)^2) * (-7 + u2 * (5 - 3 * u3) + u1 * (5 + u2 * (-3 + u3) - 
          3 * u3) + 5 * u3) + u1 * (-1 + 4 * u2 * u3) + 5 * (theta^4) * 
           (-1 + u1) * (-1 + u2) * (-1 + u3) * (-7 + (u2^2) * (-2 + u3) - 
          2 * (-4 + u3) * u3 + (u1^2) * (-2 + u2 + u3) * (1 + 2 * u2 * u3) +
      u2 * (8 + (-7 + u3) * u3) + u1 * (8 + (-7 + u3) * u3 + u2 * (-7 + 2 * 
      (5 - 2 * u3) * u3) + (u2^2) * (1 + 2 * (-2 + u3) * u3))) + theta * 
        (-7 + 9 * u2 + 9 * u3 - 2 * ((u2^2) + 3 * u2 * u3 + (u3^2)) +
      (u1^2) * (-2 + u2 * u3 * (1 + 7 * u2 * u3)) + u1 * (9 - 6 * u3 + u2 * 
      (-6 + u3 * (-5 + u2 + u3)))) + (theta^5) * (-1 + u1) * (-1 + u2) * 
        (-1 + u3) * (21 + 10 * (-3 + u3) * u3 + u2 * (-30 + (39 - 11 * u3) * 
        u3) + (u2^2) * (10 + u3 * (-11 + 2 * u3)) + u1 * (-30 + 
        (39 - 11 * u3) * u3 + u2 * (-3 + 2 * u3) * (-13 + 9 * u3) + 
        (u2^2) * (-11 - 9 * (-2 + u3) * u3)) + (u1^2) * (10 + u3 * 
        (-11 + 2 * u3) + u2 * (-11 - 9 * (-2 + u3) * u3) + (u2^2) * 
        (2 + u3 * (-9 + 8 * u3)))) + (theta^2) * (21 - u2 * (33 + (-13 + u2) * 
        u2) - 33 * u3 + 3 * (13 - 3 * u2) * u2 * u3 + (13 - 9 * u2) * 
        (u3^2) - (u3^3) + (u1^3) * (-1 + u2 * u3) * (1 + u2 * u3 * 
        (4 + u2 * u3)) - u1 * (33 - 39 * u3 + 3 * (u2^3) * u3 + 9 * (u3^2) + 
        (u2^2) * (9 + u3 * (-19 + 12 * u3)) + u2 * (-39 + u3 * (47 + u3 * 
        (-19 + 3 * u3)))) + (u1^2) * (13 - 9 * u3 + u2 *  (-9 + u3 * (19 - 12 * 
        u3 + u2 * (-12 + u3 * (4 + 3 * u2 + 3 * u3)))))) + (theta^3) * 
        ((u2^3) * (5 - 4 * u3) + (u2^2) * (-35 + 3 * (15 - 4 * u3) * u3) + 
        5 * (-1 + u3) * (7 + (-6 + u3) * u3) + u2 * (65 + u3 * (-105 + 
        (45 - 4 * u3) * u3)) - (u1^2) * (35 + 3 * u3 * (-15 + 4 * u3) + 
        (u2^3) * u3 * (13 + (-11 + u3) * u3) + (u2^2) * (12 + u3 * (-68 + 
        (63 - 11 * u3) * u3)) + u2 * (-45 + u3 * (99 + u3 * 
        (-68 + 13 * u3)))) + u1 * (65 - (u2^3) * (-1 + u3) * (-4 + 13 * u3) + 
        u3 * (-105 + (45 - 4 * u3) * u3) + (u2^2) * (45 + u3 * (-99 + 
        (68 - 13 * u3) * u3)) + u2 * (-105 + u3 * (185 + u3 * 
        (-99 + 17 * u3)))) + (u1^3) * (5 - 4 * u3 + u2 * (-4 + u3 * 
        (17 - 13 * u3 + u2 * (-13 + u3 * (11 - u3 + u2 * (-1 + 3 * u3))))))))) /
        ((-1 - theta * (-2 + u1 + u2) + (theta^2) * (-1 + u1) * (-1 + u2) * 
        (-1 + u3) + theta * (-1 + u1 * u2) * u3) * (1 + (theta^6) * 
        ((u1 - 1)^2) * ((u2 - 1)^2) * ((u3 - 1)^2) + 2 * theta * 
        (-3 + u1 + u2 + u3 + 2 * u1 * u2 * u3) + 2 * (theta^5) * (-1 + u1) * 
        (-1 + u2) * (-1 + u3) * (3 + u2 * (-2 + u3) - 2 * u3 + u1 * 
        (-2 + u2 + u3)) + (theta^2) * (-10 * u1 + (u1^2) + 4 * u1 * u2 + 
        (u2^2) + 4 * (u1 + u2 + u1 * u2 * (-3 + u1 + u2)) * u3 + 
        (1 + u1 * u2 * (4 + u1 * u2)) * (u3^2) - 5 * (-3 + 2 * u2 + 2 * u3)) + 
        (theta^4) * (15 - 20 * u3 + 6 * (u3^2) + (u2^2) * (6 + (-6 + u3) * u3) + 
        u2 * (-20 - 6 * (-4 + u3) * u3) + (u1^2) * (6 + (u2^2) * 
        ((1 - 2 * u3)^2) + (-6 + u3) * u3 + u2 * (-6 - 4 * (-2 + u3) * u3)) + 
        u1 * (-20 - 6 * (-4 + u3) * u3 + 4 * u2 * (-2 + u3) * (-3 + 2 * u3) + 
        (u2^2) * (-6 - 4 * (-2 + u3) * u3))) + 2 * (theta^3) * ((u2^2) * 
        (-2 + u3) + u2 * (-8 + u3) * u3 - 2 * (u3^2) + 10 * (-1 + u2 + u3) +
        (u1^2) * (-2 + u2 + u3) * (1 + 2 * u2 * u3) + u1 * 
        (10 + (-8 + u3) * u3 + (u2^2) * (1 + 2 * (-2 + u3) * u3) - 
        2 * u2 * (4 + u3 * (-5 + 2 * u3))))))
}

# 3d second derivative wrt theta of the log-density
.amh.123.S <- function(theta, u) { 
  if (is.null(dim(u))) {
    u1 <- u[1]
    u2 <- u[2]
    u3 <- u[3]
  } else {
    u1 <- u[, 1]
    u2 <- u[, 2]
    u3 <- u[, 3]
  }

  (2 * (-2 * ((1 - u2 + (theta^7) * ((u1 - 1)^3) * ((u2 - 1)^3) * 
  ((u3 - 1)^3) - u3 - (theta^6) * ((u1 - 1)^2) * ((u2 - 1)^2) * 
  ((u3 - 1)^2) * (-7 + u2 * (5 - 3 * u3) + u1 * (5 + u2 * (-3 + u3) - 3 * u3) + 
  5 * u3) + u1 * (-1 + 4 * u2 * u3) + 5 * (theta^4) * (-1 + u1) * (-1 + u2) * 
  (-1 + u3) * (-7 + (u2^2) * (-2 + u3) - 2 * (-4 + u3) * u3 + (u1^2) * 
  (-2 + u2 + u3) * (1 + 2 * u2 * u3) + u2 * (8 + (-7 + u3) * u3) + u1 * 
  (8 + (-7 + u3) * u3 + u2 * (-7 + 2 * (5 - 2 * u3) * u3) + (u2^2) * 
  (1 + 2 * (-2 + u3) * u3))) + theta * (-7 + 9 * u2 + 9 * u3 - 2 * ((u2^2) + 
  3 * u2 * u3 + (u3^2)) + (u1^2) * (-2 + u2 * u3 * (1 + 7 * u2 * u3)) + u1 * 
  (9 - 6 * u3 + u2 * (-6 + u3 * (-5 + u2 + u3)))) + (theta^5) * (-1 + u1) * 
  (-1 + u2) * (-1 + u3) * (21 + 10 * (-3 + u3) * u3 + u2 * (-30 + 
  (39 - 11 * u3) * u3) + (u2^2) * (10 + u3 * (-11 + 2 * u3)) + u1 * (-30 + 
  (39 - 11 * u3) * u3 + u2 * (-3 + 2 * u3) * (-13 + 9 * u3) + (u2^2) * 
  (-11 - 9 * (-2 + u3) * u3)) + (u1^2) * (10 + u3 * (-11 + 2 * u3) + u2 * 
  (-11 - 9 * (-2 + u3) * u3) + (u2^2) * (2 + u3 * (-9 + 8 * u3)))) +
  (theta^2) * (21 - u2 * (33 + (-13 + u2) * u2) - 33 * u3 + 3 * (13 - 3 * u2) * 
  u2 * u3 + (13 - 9 * u2) * (u3^2) - (u3^3) + (u1^3) * (-1 + u2 * u3) * 
  (1 + u2 * u3 * (4 + u2 * u3)) - u1 * (33 - 39 * u3 + 3 * (u2^3) * u3 + 9 * 
  (u3^2) + (u2^2) * (9 + u3 * (-19 + 12 * u3)) + u2 * (-39 + u3 * (47 + u3 * 
  (-19 + 3 * u3)))) + (u1^2) * (13 - 9 * u3 + u2 * (-9 + u3 * (19 - 12 * 
  u3 + u2 * (-12 + u3 * (4 + 3 * u2 + 3 * u3)))))) + (theta^3) * ((u2^3) * 
  (5 - 4 * u3) + (u2^2) * (-35 + 3 * (15 - 4 * u3) * u3) + 5 * (-1 + u3) * 
  (7 + (-6 + u3) * u3) + u2 * (65 + u3 * (-105 + (45 - 4 * u3) * u3)) -
  (u1^2) * (35 + 3 * u3 * (-15 + 4 * u3) + (u2^3) * u3 * (13 + (-11 + u3) * 
  u3) + (u2^2) * (12 + u3 * (-68 + (63 - 11 * u3) * u3)) + u2 * (-45 + u3 * 
  (99 + u3 * (-68 + 13 * u3)))) + u1 * (65 - (u2^3) * (-1 + u3) * 
  (-4 + 13 * u3) + u3 * (-105 + (45 - 4 * u3) * u3) + (u2^2) * (45 + u3 * 
  (-99 + (68 - 13 * u3) * u3)) + u2 * (-105 + u3 * (185 + u3 * 
  (-99 + 17 * u3)))) + (u1^3) * (5 - 4 * u3 + u2 * (-4 + u3 * 
  (17 - 13 * u3 + u2 * (-13 + u3 * (11 - u3 + u2 * (-1 + 3 * u3))))))))^2) +
  (1 + (theta^6) * ((u1 - 1)^2) * ((u2 - 1)^2) * ((u3 - 1)^2) + 2 * theta * 
  (-3 + u1 + u2 + u3 + 2 * u1 * u2 * u3) + 2 * (theta^5) * (-1 + u1) * 
  (-1 + u2) * (-1 + u3) * (3 + u2 * (-2 + u3) - 2 * u3 + u1 * (-2 + u2 + u3)) +
  (theta^2) * (-10 * u1 + (u1^2) + 4 * u1 * u2 + (u2^2) + 4 * 
  (u1 + u2 + u1 * u2 * (-3 + u1 + u2)) * u3 + (1 + u1 * u2 * (4 + u1 * u2)) * 
  (u3^2) - 5 * (-3 + 2 * u2 + 2 * u3)) + (theta^4) * (15 - 20 * u3 + 6 * 
  (u3^2) + (u2^2) * (6 + (-6 + u3) * u3) + u2 * (-20 - 6 * (-4 + u3) * u3) + 
  (u1^2) * (6 + (u2^2) * ((1 - 2 * u3)^2) + (-6 + u3) * u3 + u2 * (-6 - 4 * 
  (-2 + u3) * u3)) + u1 * (-20 - 6 * (-4 + u3) * u3 + 4 * u2 * (-2 + u3) * 
  (-3 + 2 * u3) + (u2^2) * (-6 - 4 * (-2 + u3) * u3))) + 2 * (theta^3) * 
  ((u2^2) * (-2 + u3) + u2 * (-8 + u3) * u3 - 2 * (u3^2) + 10 * 
  (-1 + u2 + u3) + (u1^2) * (-2 + u2 + u3) * (1 + 2 * u2 * u3) + u1 * (10 + 
  (-8 + u3) * u3 + (u2^2) * (1 + 2 * (-2 + u3) * u3) - 2 * u2 * (4 + u3 * 
  (-5 + 2 * u3))))) * (3 - 6 * u1 - 6 * u2 + 3 * (u2^2) + 3 * (theta^8) * 
  ((u1 - 1)^4) * ((u2 - 1)^4) * ((u3 - 1)^4) - 6 * u3 + 4 * u2 * u3 + 
  3 * (u3^2) - 6 * (theta^7) * ((u1 - 1)^3) * ((u2 - 1)^3) * ((u3 - 1)^3) * 
  (-4 + u2 * (3 - 2 * u3) + u1 * (3 + u2 * (-2 + u3) - 2 * u3) + 3 * u3) + 
  4 * u1 * (u2 + u3 - 2 * u2 * u3 * (-5 + 3 * u2 + 3 * u3)) + 3 * (u1^2) * 
  (1 + u2 * u3 * (-8 + 9 * u2 * u3)) + 6 * (theta^5) * ((u1 - 1)^2) * 
  ((u2 - 1)^2) * ((u3 - 1)^2) * (-28 + 5 * (u2^2) * (-2 + u3) + 5 * 
  (7 - 2 * u3) * u3 + u2 * (-5 + u3) * (-7 + 5 * u3) + 5 * (u1^2) * 
  (-2 + u2 + u3) * (1 + 2 * u2 * u3) + u1 * ((-5 + u3) * (-7 + 5 * u3) + 
  u2 * (-32 + (49 - 20 * u3) * u3) + 5 * (u2^2) * (1 + 2 * (-2 + u3) * u3))) +
  (theta^6) * ((u1 - 1)^2) * ((u2 - 1)^2) * ((u3 - 1)^2) * (84 + 9 * u3 * 
  (-14 + 5 * u3) - 2 * u2 * (63 - 86 * u3 + 26 * (u3^2)) + (u2^2) * 
  (45 + 2 * u3 * (-26 + 5 * u3)) + (u1^2) * (45 + 2 * u3 * (-26 + 5 * u3) + 
  u2 * (-52 + 94 * u3 - 48 * (u3^2)) + (u2^2) * (10 + u3 * (-48 + 41 * u3))) +
  2 * u1 * (-63 + 86 * u3 - 26 * (u3^2) + (u2^2) * (-26 + (47 - 24 * u3) * 
  u3) + u2 * (86 + u3 * (-127 + 47 * u3)))) + 6 * theta * (-4 + 9 * u2 + 9 * 
  u3 + (u2 + u3) * ((u2^2) + 3 * u2 * (-2 + u3) + (-6 + u3) * u3) + (u1^3) * 
  (-1 + u2 * u3) * (-1 + u2 * u3 * (2 + 5 * u2 * u3)) + u1 * 
  (((3 - 2 * u3)^2) - 3 * (u2^3) * u3 + (u2^2) * (4 + 2 * (11 - 8 * u3) * u3) - 
  u2 * (-4 + 3 * u3) * (-3 + (-6 + u3) * u3)) + (u1^2) * (-6 + 4 * u3 + u2 * 
  (4 - u3 * (-22 + 16 * u3 + u2 * (16 + u3 * (-14 + 3 * u2 + 3 * u3)))))) +
  3 * (theta^4) * (-1 + u1) * (-1 + u2) * (-1 + u3) * ((u2^3) * 
  (15 - 11 * u3) + (u2^2) * (-85 + 7 * (15 - 4 * u3) * u3) + 5 * (-1 + u3) * 
  (14 + u3 * (-14 + 3 * u3)) + u2 * (140 + u3 * (-230 + (105 - 11 * u3) * 
  u3)) + (u1^2) * (-85 + 7 * (15 - 4 * u3) * u3 + (u2^3) * u3 * (-32 + u3 * 
  (19 + u3)) + u2 * (105 + u3 * (-221 + 8 * (19 - 4 * u3) * u3)) + (u2^2) * 
  (-28 + u3 * (152 + u3 * (-127 + 19 * u3)))) + u1 * (140 - (u2^3) * 
  (-1 + u3) * (-11 + 32 * u3) + u3 * (-230 + (105 - 11 * u3) * u3) + (u2^2) * 
  (105 + u3 * (-221 + 8 * (19 - 4 * u3) * u3)) + u2 * (-230 + u3 * 
  (400 + u3 * (-221 + 43 * u3)))) + (u1^3) * (15 - 11 * u3 + u2 * (-11 + u3 * 
  (43 - 32 * u3 + u2 * (-32 + u3 * (19 + u2 + u3 + 7 * u2 * u3)))))) +
  2 * (theta^3) * ((u2^4) * (-9 + 8 * u3) - 3 * ((u3 - 1)^2) * (28 + 3 * 
  (-7 + u3) * u3) + 3 * (u2^3) * (27 + 2 * u3 * (-20 + 7 * u3)) + 3 * (u2^2) * 
  (-73 + 2 * u3 * (71 + 7 * (-6 + u3) * u3)) + u2 * (231 + 2 * u3 * 
  (-272 + u3 * (213 + 4 * (-15 + u3) * u3))) + 3 * (u1^2) * (-73 + 2 * u3 * 
  (71 + 7 * (-6 + u3) * u3) + (u2^4) * u3 * (7 + 2 * u3 - 8 * (u3^2)) +
  (u2^2) * (-1 + u3) * (84 + u3 * (-105 + 2 * u3 * (21 + u3))) + u2 * 
  (142 + 7 * u3 * (-40 + u3 * (27 + (-8 + u3) * u3))) + (u2^3) * 
  (14 - 8 * u3 * (7 + u3 * (-5 + (-1 + u3) * u3)))) + (u1^3) * (81 + 6 * u3 * 
  (-20 + 7 * u3) + (u2^4) * (u3^2) * (-24 + (24 - 5 * u3) * u3) + 2 * u2 * 
  (-60 + u3 * (121 + 21 * (-4 + u3) * u3)) + 6 * (u2^3) * u3 * (7 + u3 * 
  (4 + u3 * (-13 + 4 * u3))) - 6 * (u2^2) * (-7 + 4 * u3 * (7 + u3 * (-5 + 
  (-1 + u3) * u3)))) + u1 * (231 + (u2^4) * (8 + 7 * u3 * (-4 + 3 * u3)) + 
  2 * u3 * (-272 + u3 * (213 + 4 * (-15 + u3) * u3)) + 2 * (u2^3) * (-60 + 
  u3 * (121 + 21 * (-4 + u3) * u3)) + 3 * (u2^2) * (142 + 7 * u3 * 
  (-40 + u3 * (27 + (-8 + u3) * u3))) - 2 * u2 * (272 + u3 * (-583 + u3 * 
  (420 + u3 * (-121 + 14 * u3))))) + (u1^4) * (-1 + u2 * u3) * (9 - 8 * u3 + 
  u2 * (-8 + u3 * (37 - 29 * u3 + u2 * (-29 + u3 * (31 - 5 * u3 + u2 * 
  (-5 + 7 * u3))))))) + 3 * (theta^2) * ((u2^4) + 4 * (u2^3) * (-4 + 3 * u3) + 
  ((u3 - 1)^2) * (28 + (-14 + u3) * u3) + (u2^2) * (57 + 4 * u3 * 
  (-19 + 6 * u3)) + (u1^4) * ((-1 + u2 * u3)^2) * (1 + u2 * u3 * 
  (4 + u2 * u3)) + 2 * u2 * (-35 + 2 * u3 * (33 + u3 * (-19 + 3 * u3))) + 
  (u1^2) * (57 - 6 * (u2^4) * (u3^2) + 4 * u3 * (-19 + 6 * u3) - 2 * 
  (u2^3) * u3 * (5 + 4 * u3 * (-7 + 5 * u3)) + u2 * (-76 + 2 * u3 * 
  (11 + (24 - 5 * u3) * u3)) + (u2^2) * (24 + u3 * (48 + u3 * (-111 + 56 * 
  u3 - 6 * (u3^2))))) + 2 * u1 * (-35 + (u2^4) * u3 + (u2^3) * (6 + u3 - 5 * 
  (u3^2)) + 2 * u3 * (33 + u3 * (-19 + 3 * u3)) + (u2^2) * (-38 + u3 * 
  (11 + (24 - 5 * u3) * u3)) + u2 * (66 + u3 * (-73 + u3 * (11 + u3 + 
  (u3^2))))) + 2 * (u1^3) * (-8 + 6 * u3 + u2 * (6 + u3 * (1 - 5 * u3 + u2 * 
  (-5 + u3 * (28 - 20 * u3 + u2 * (-20 + u3 * (15 + u2 + u3))))))))))) /
  (((-1 - theta * (-2 + u1 + u2) + (theta^2) * (-1 + u1) * (-1 + u2) * 
  (-1 + u3) + theta * (-1 + u1 * u2) * u3)^2) * ((1 + (theta^6) * 
  ((u1 - 1)^2) * ((u2 - 1)^2) * ((u3 - 1)^2) + 2 * theta * (-3 + u1 + u2 + 
  u3 + 2 * u1 * u2 * u3) + 2 * (theta^5) * (-1 + u1) * (-1 + u2) * (-1 + u3) * 
  (3 + u2 * (-2 + u3) - 2 * u3 + u1 * (-2 + u2 + u3)) + (theta^2) * (-10 * 
  u1 + (u1^2) + 4 * u1 * u2 + (u2^2) + 4 * (u1 + u2 + u1 * u2 * (-3 + u1 + 
  u2)) * u3 + (1 + u1 * u2 * (4 + u1 * u2)) * (u3^2) - 5 * 
  (-3 + 2 * u2 + 2 * u3)) + (theta^4) * (15 - 20 * u3 + 6 * (u3^2) + (u2^2) * 
  (6 + (-6 + u3) * u3) + u2 * (-20 - 6 * (-4 + u3) * u3) + (u1^2) * 
  (6 + (u2^2) * ((1 - 2 * u3)^2) + (-6 + u3) * u3 + u2 * (-6 - 4 * 
  (-2 + u3) * u3)) + u1 * (-20 - 6 * (-4 + u3) * u3 + 4 * u2 * (-2 + u3) * 
  (-3 + 2 * u3) + (u2^2) * (-6 - 4 * (-2 + u3) * u3))) + 2 * (theta^3) * 
  ((u2^2) * (-2 + u3) + u2 * (-8 + u3) * u3 - 2 * (u3^2) + 10 * 
  (-1 + u2 + u3) + (u1^2) * (-2 + u2 + u3) * (1 + 2 * u2 * u3) + u1 * 
  (10 + (-8 + u3) * u3 + (u2^2) * (1 + 2 * (-2 + u3) * u3) - 2 * u2 * 
  (4 + u3 * (-5 + 2 * u3)))))^2))
}
