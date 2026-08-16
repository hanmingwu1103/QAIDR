# Exceptional-value enumeration machinery for the lambda/nu sensitivity
# analysis (Prompt E 2026-07-17, three-provider agreed contract). Sourced by
# b2_sensitivity.r and by b2_breaks_selftest.r. Requires TOL_ROOT, TOL_DISC,
# TOL_ENDPOINT in the sourcing environment (defaults set here if absent).
if (!exists("TOL_ROOT")) TOL_ROOT <- 1e-10
if (!exists("TOL_DISC")) TOL_DISC <- 128 * .Machine$double.eps
if (!exists("TOL_ENDPOINT")) TOL_ENDPOINT <- 8 * TOL_DISC

## ---- exceptional-value machinery (Prompt E, Codex-audited) ----------------

.root_records <- function() {
  data.frame(space = character(), anchor = integer(), j = integer(),
    jprime = integer(), q2 = numeric(), q1 = numeric(), q0 = numeric(),
    s2 = numeric(), s1 = numeric(), s0 = numeric(), raw = numeric(),
    location = character(), tangency = logical())
}

.summarize_internal <- function(records, tol_root = TOL_ROOT) {
  z <- records[records$location == "internal", , drop = FALSE]
  if (!nrow(z)) return(data.frame(key = numeric(), value = numeric(),
    raw_min = numeric(), raw_max = numeric(), n_source = integer(),
    tangency = logical()))
  z$key <- round(z$raw / tol_root)       # dedup key only; raw kept
  out <- lapply(split(seq_len(nrow(z)), z$key), function(ii) {
    zz <- z[ii, , drop = FALSE]
    candidate <- if (any(zz$tangency)) zz$raw[zz$tangency] else zz$raw
    data.frame(key = zz$key[1], value = sort(candidate)[1],
      raw_min = min(zz$raw), raw_max = max(zz$raw),
      n_source = nrow(zz), tangency = any(zz$tangency))
  })
  ans <- do.call(rbind, out)
  ans[order(ans$value), , drop = FALSE]
}

.finish_breaks <- function(records, persistent, domain,
                           tol_root = TOL_ROOT, tol_disc = TOL_DISC) {
  im <- .summarize_internal(records, tol_root)
  list(internal = im$value,
       boundary = sort(unique(records$raw[records$location == "boundary"])),
       tangency = im$value[im$tangency],
       persistent = persistent, internal_meta = im, records = records,
       domain = domain, tol_root = tol_root, tol_disc = tol_disc)
}

classify_polynomials <- function(coefs, domain, tol_root = TOL_ROOT,
                                 tol_disc = TOL_DISC) {
  need <- c("space", "anchor", "j", "jprime", "q2", "q1", "q0")
  stopifnot(all(need %in% names(coefs)), length(domain) == 2L,
            domain[1] < domain[2])
  if (!"s2" %in% names(coefs)) coefs$s2 <- abs(coefs$q2)
  if (!"s1" %in% names(coefs)) coefs$s1 <- abs(coefs$q1)
  if (!"s0" %in% names(coefs)) coefs$s0 <- abs(coefs$q0)
  records <- .root_records()
  persistent <- coefs[FALSE, , drop = FALSE]

  add <- function(r, x, where, tangent) data.frame(
    space = as.character(r$space), anchor = as.integer(r$anchor),
    j = as.integer(r$j), jprime = as.integer(r$jprime),
    q2 = r$q2, q1 = r$q1, q0 = r$q0, s2 = r$s2, s1 = r$s1, s0 = r$s0,
    raw = x, location = where, tangency = tangent)

  for (ii in seq_len(nrow(coefs))) {
    r <- coefs[ii, , drop = FALSE]
    q2 <- r$q2; q1 <- r$q1; q0 <- r$q0
    if (q2 == 0 && q1 == 0 && q0 == 0) {
      persistent <- rbind(persistent, r)   # identically-zero comparison
      next
    }

    zero_at <- function(x) {
      terms <- c(q2 * x^2, q1 * x, q0)
      source_scale <- r$s2 * x^2 + r$s1 * abs(x) + r$s0
      abs(sum(terms)) <= tol_disc *
        max(sum(abs(terms)), source_scale, .Machine$double.xmin)
    }
    endpoint_zero <- vapply(domain, zero_at, logical(1))
    for (x in domain[endpoint_zero])
      records <- rbind(records, add(r, x, "boundary", FALSE))

    tangent <- FALSE
    if (q2 == 0) {
      roots <- if (q1 == 0) numeric() else -q0 / q1
    } else {
      disc <- q1^2 - 4 * q2 * q0
      disc_scale <- max(q1^2, abs(4 * q2 * q0), .Machine$double.xmin)
      if (abs(disc) <= tol_disc * disc_scale) {
        roots <- -q1 / (2 * q2); tangent <- TRUE      # double root/tangency
      } else if (disc < 0) {
        roots <- numeric()
      } else {
        sq <- sqrt(disc)                              # stable quadratic
        z <- -0.5 * (q1 + if (q1 >= 0) sq else -sq)
        roots <- c(z / q2, q0 / z)
      }
    }

    for (x in roots[is.finite(roots)]) {
      # Suppress only duplicates of an independently verified endpoint root.
      at_endpoint <- endpoint_zero &
        (abs(x - domain) <= TOL_ENDPOINT * pmax(1, abs(domain)) |
         round(x / tol_root) == round(domain / tol_root))
      if (any(at_endpoint)) next
      if (x > domain[1] && x < domain[2])
        records <- rbind(records, add(r, x, "internal", tangent))
    }
  }
  .finish_breaks(records, persistent, domain, tol_root, tol_disc)
}

lambda_breaks <- function(centers, radii, space = "space") {
  n <- nrow(centers)
  d0 <- d1 <- matrix(0, n, n)
  for (i in seq_len(n - 1L)) for (j in (i + 1L):n) {
    dc <- abs(centers[i, ] - centers[j, ]); sr <- radii[i, ] + radii[j, ]
    d0[i, j] <- d0[j, i] <- sqrt(sum(pmax(0, dc - sr)^2))
    d1[i, j] <- d1[j, i] <- sqrt(sum((dc + sr)^2))
  }
  rows <- vector("list", n * choose(n - 1, 2)); rr <- 0L
  for (i in seq_len(n)) {
    js <- setdiff(seq_len(n), i)
    cmb <- utils::combn(seq_along(js), 2L)
    a0 <- d0[i, js]; a1 <- d1[i, js] - a0
    for (k in seq_len(ncol(cmb))) {
      a <- cmb[1, k]; b <- cmb[2, k]; rr <- rr + 1L
      rows[[rr]] <- data.frame(space = space, anchor = i, j = js[a],
        jprime = js[b], q2 = 0, q1 = a1[a] - a1[b], q0 = a0[a] - a0[b],
        s2 = 0, s1 = abs(a1[a]) + abs(a1[b]), s0 = abs(a0[a]) + abs(a0[b]))
    }
  }
  classify_polynomials(do.call(rbind, rows), c(0, 1))
}

nu_breaks <- function(centers, radii, space = "space") {
  n <- nrow(centers); p <- ncol(centers)
  L <- centers - radii; U <- centers + radii
  A <- B <- array(0, c(n, n, p))
  for (i in seq_len(n - 1L)) for (j in (i + 1L):n) for (k in seq_len(p)) {
    meet <- max(0, min(U[i, k], U[j, k]) - max(L[i, k], L[j, k]))
    hull <- max(U[i, k], U[j, k]) - min(L[i, k], L[j, k])
    A[i, j, k] <- A[j, i, k] <- hull - meet
    B[i, j, k] <- B[j, i, k] <- 2 * meet - (U[i, k] - L[i, k]) - (U[j, k] - L[j, k])
  }
  rows <- vector("list", n * choose(n - 1, 2)); rr <- 0L
  for (i in seq_len(n)) {
    js <- setdiff(seq_len(n), i)
    c2 <- vapply(js, function(j) sum(B[i, j, ]^2), numeric(1))
    c1 <- vapply(js, function(j) 2 * sum(A[i, j, ] * B[i, j, ]), numeric(1))
    c0 <- vapply(js, function(j) sum(A[i, j, ]^2), numeric(1))
    cmb <- utils::combn(seq_along(js), 2L)
    for (k in seq_len(ncol(cmb))) {
      a <- cmb[1, k]; b <- cmb[2, k]; rr <- rr + 1L
      rows[[rr]] <- data.frame(space = space, anchor = i, j = js[a],
        jprime = js[b], q2 = c2[a] - c2[b], q1 = c1[a] - c1[b],
        q0 = c0[a] - c0[b],
        s2 = abs(c2[a]) + abs(c2[b]), s1 = abs(c1[a]) + abs(c1[b]),
        s0 = abs(c0[a]) + abs(c0[b]))
    }
  }
  classify_polynomials(do.call(rbind, rows), c(0, 0.5))
}

union_breaks <- function(...) {
  xs <- list(...)
  stopifnot(length(xs) > 0L,
            all(vapply(xs, function(x) identical(x$domain, xs[[1]]$domain),
                       logical(1))))
  records <- do.call(rbind, lapply(xs, `[[`, "records"))
  persistent <- do.call(rbind, lapply(xs, `[[`, "persistent"))
  .finish_breaks(records, persistent, xs[[1]]$domain,
                 xs[[1]]$tol_root, xs[[1]]$tol_disc)
}
