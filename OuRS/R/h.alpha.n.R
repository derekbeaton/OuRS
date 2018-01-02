  ## stolen.
h.alpha.n <- function (alpha, n, p)
{
  n2 <- (n + p + 1)%/%2
  floor(2 * n2 - n + 2 * (n - n2) * alpha)
}
