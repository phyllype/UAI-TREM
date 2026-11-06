# T3 + N1b do E0 (velocidade, 2026-11-02): a fatoracao e a inversa seletiva em paralelo pelas
# SUBARVORES da arvore de eliminacao, e a ordenacao AMD (grau minimo aproximado) ao lado do grau
# minimo exato (UAITREM_ORDEM = "gm" | "amd").
#
# Portoes:
#  (1) com a MESMA ordenacao, 1 e 4 threads dao o mesmo fator e a mesma inversa seletiva BIT A BIT,
#      com a particao em subarvores de fato ligada (o diagnostico UAITREM_T3_TEMPOS diz quantas),
#      com e sem cauda densa;
#  (2) a AMD da uma permutacao valida, resolve o sistema, e o log-determinante e a diagonal da
#      inversa batem com os do grau minimo exato (outra ordem de soma: tolerancia, nao bit a bit);
#  (3) o enchimento da AMD fica perto do exato (pedigree, grade, hubs) e o clique genomico continua
#      no FIM (a cauda densa), que e de onde vem a forma fechada;
#  (4) casos de borda: n = 1 e 2, so diagonal, componentes desconexas;
#  (5) o ajuste inteiro (model, model_mt) com a AMD da o mesmo -2logL e as mesmas solucoes que com o
#      exato, dentro da tolerancia; e a AMD e o exato com 1 e 3 threads, cada um bit a bit consigo.

com_ordem <- function(ordem, expr) {
  antes <- Sys.getenv("UAITREM_ORDEM", NA)
  Sys.setenv(UAITREM_ORDEM = ordem)
  on.exit(if (is.na(antes)) Sys.unsetenv("UAITREM_ORDEM") else Sys.setenv(UAITREM_ORDEM = antes))
  force(expr)
}
com_threads <- function(n, expr) {
  antes <- br_threads()
  on.exit(br_threads(antes$threads, lapack = antes$lapack))
  br_threads(n)
  force(expr)
}
# quantas subarvores a fatoracao (chol) e a inversa seletiva (selinv) usaram
subarvores <- function(expr) {
  antes <- Sys.getenv("UAITREM_T3_TEMPOS", NA)
  Sys.setenv(UAITREM_T3_TEMPOS = "1")
  on.exit(if (is.na(antes)) Sys.unsetenv("UAITREM_T3_TEMPOS") else Sys.setenv(UAITREM_T3_TEMPOS = antes))
  msg <- utils::capture.output(v <- force(expr), type = "message")
  num <- function(re) {
    l <- grep(re, msg, value = TRUE)
    if (!length(l)) return(NA_integer_)
    as.integer(sub(".*subarvores [0-9.]+ s \\(([0-9]+)\\).*", "\\1", l[length(l)]))
  }
  list(v = v, chol = num("^\\[t3 chol\\]"), selinv = num("^\\[t3 selinv\\]"))
}

# grade 2D (5 pontos) com diagonal dominante: arvore larga, muito trabalho na parte esparsa
grade <- function(k) {
  n <- k * k
  id <- function(r, c) (c - 1) * k + r
  i <- seq_len(n); j <- seq_len(n); x <- rep(4.5, n)
  r <- rep(seq_len(k), k); cc <- rep(seq_len(k), each = k)
  v <- r < k; i <- c(i, id(r[v] + 1, cc[v])); j <- c(j, id(r[v], cc[v])); x <- c(x, rep(-1, sum(v)))
  h <- cc < k; i <- c(i, id(r[h], cc[h] + 1)); j <- c(j, id(r[h], cc[h])); x <- c(x, rep(-1, sum(h)))
  list(i = i, j = j, x = x, n = n)
}
# a mesma grade com um clique denso (genomico) ligado a uma parte dos nos
grade_clique <- function(k, T, seed = 3) {
  set.seed(seed)
  g <- grade(k)
  M <- matrix(stats::rnorm(T * (T + 4)), T)
  D <- tcrossprod(M) / T + diag(T)
  lo <- which(lower.tri(D, diag = TRUE), arr.ind = TRUE)
  ligados <- sample(g$n, 3 * T)
  list(i = c(g$i, g$n + lo[, 1], g$n + sample(T, length(ligados), TRUE)),
       j = c(g$j, g$n + lo[, 2], ligados),
       x = c(g$x, D[lo], rep(0.05, length(ligados))), n = g$n + T)
}
ped_sim <- function(n = 6000, por_ger = 1000, machos = 20, femeas = 300, seed = 5) {
  set.seed(seed)
  id <- sprintf("a%06d", seq_len(n)); pa <- ma <- rep("0", n)
  for (g in 2:(n / por_ger)) {
    ant <- ((g - 2) * por_ger + 1):((g - 1) * por_ger)
    agora <- ((g - 1) * por_ger + 1):(g * por_ger)
    pa[agora] <- id[sample(ant[1:machos], por_ger, TRUE)]
    ma[agora] <- id[sample(ant[(machos + 1):(machos + femeas)], por_ger, TRUE)]
  }
  data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
}
mme_ped <- function(ped, ncg = 40, lam = 2, seed = 9) {
  set.seed(seed)
  ai <- a_inverse(ped)
  na <- length(ai$id); cg <- sample(ncg, na, TRUE)
  i <- c(seq_len(ncg), ncg + seq_len(na), ncg + ai$i, ncg + seq_len(na))
  j <- c(seq_len(ncg), cg, ncg + ai$j, ncg + seq_len(na))
  x <- c(tabulate(cg, ncg), rep(1, na), lam * ai$x, rep(1, na))
  list(i = i, j = j, x = x, n = ncg + na)
}
densa <- function(a) {
  A <- matrix(0, a$n, a$n); A[cbind(a$i, a$j)] <- a$x; A[cbind(a$j, a$i)] <- a$x; A
}
diag_si <- function(si, n) { k <- si$i == si$j; v <- numeric(n); v[si$i[k]] <- si$x[k]; v }

test_that("(1) com a mesma ordem, 1 e 4 threads dao o mesmo fator e a mesma inversa, bit a bit, com subarvores", {
  for (ordem in c("gm", "amd")) {
    for (a in list(grade(140), grade_clique(110, 300))) {
      r1 <- com_ordem(ordem, com_threads(1, list(ch = sparse_chol(a), si = selected_inverse(a))))
      r4 <- com_ordem(ordem, com_threads(4, subarvores(list(ch = sparse_chol(a), si = selected_inverse(a)))))
      expect_gt(r4$chol, 1L)                # a particao ligou na fatoracao
      expect_gt(r4$selinv, 1L)              # e na inversa seletiva
      expect_identical(r1$ch$perm, r4$v$ch$perm)
      expect_identical(r1$ch$L$x, r4$v$ch$L$x)
      expect_identical(r1$ch$L$i, r4$v$ch$L$i)
      expect_identical(r1$ch$logdet, r4$v$ch$logdet)
      expect_identical(r1$si$x, r4$v$si$x)
      expect_identical(r1$si$i, r4$v$si$i)
    }
  }
})

test_that("(1b) o pedigree (MME de um modelo animal) tambem, com 1, 3 e 4 threads", {
  a <- mme_ped(ped_sim(12000, 1500, 25, 400))
  r1 <- com_threads(1, list(ch = sparse_chol(a), si = selected_inverse(a)))
  for (nt in c(3, 4)) {
    rn <- com_threads(nt, list(ch = sparse_chol(a), si = selected_inverse(a)))
    expect_identical(r1$ch$L$x, rn$ch$L$x)
    expect_identical(r1$si$x, rn$si$x)
  }
})

test_that("(2) a AMD resolve o sistema e bate com o grau minimo exato", {
  for (a in list(grade(60), grade_clique(40, 150), mme_ped(ped_sim(4000, 500, 10, 150)))) {
    fa <- com_ordem("amd", sparse_chol(a))
    fg <- com_ordem("gm", sparse_chol(a))
    expect_equal(sort(fa$perm), seq_len(a$n))
    expect_false(identical(fa$perm, fg$perm))      # e outra ordem de verdade
    expect_equal(fa$logdet, fg$logdet, tolerance = 1e-12)
    b <- stats::rnorm(a$n)
    xa <- com_ordem("amd", sparse_solve(a, b))
    xg <- com_ordem("gm", sparse_solve(a, b))
    expect_lt(max(abs(xa - xg)) / max(abs(xg)), 1e-10)
    da <- diag_si(com_ordem("amd", selected_inverse(a)), a$n)
    dg <- diag_si(com_ordem("gm", selected_inverse(a)), a$n)
    expect_lt(max(abs(da - dg) / abs(dg)), 1e-10)
  }
  # contra a conta densa do R, num caso pequeno
  a <- grade_clique(20, 60)
  A <- densa(a)
  fa <- com_ordem("amd", sparse_chol(a))
  expect_equal(fa$logdet, as.numeric(determinant(A)$modulus), tolerance = 1e-12)
  da <- diag_si(com_ordem("amd", selected_inverse(a)), a$n)
  expect_lt(max(abs(da - diag(solve(A)))), 1e-10)
})

test_that("(3) o enchimento da AMD fica perto do exato e o clique continua no fim", {
  for (a in list(grade(80), mme_ped(ped_sim(6000)))) {
    na <- length(com_ordem("amd", sparse_chol(a))$L$x)
    ng <- length(com_ordem("gm", sparse_chol(a))$L$x)
    expect_lt(na, 1.15 * ng)
  }
  a <- grade_clique(60, 400)
  fa <- com_ordem("amd", sparse_chol(a)); fg <- com_ordem("gm", sparse_chol(a))
  # medido: o exato poe 374 dos 400 nos do clique numa cauda de 376; a AMD, 300 numa cauda de 404,
  # com enchimento menor (248 mil contra 270 mil). O que se exige e que a cauda densa continue
  # existindo e do tamanho do clique, e que o fator nao fique maior
  expect_gte(fa$dense_block, 300L)
  expect_lte(fa$dense_block, fg$dense_block + 40L)
  expect_lte(length(fa$L$x), 1.05 * length(fg$L$x))
  cauda <- fa$perm[(a$n - fa$dense_block + 1):a$n]
  expect_gte(sum((a$n - 399):a$n %in% cauda), 300L)
})

test_that("(4) casos de borda da AMD", {
  um <- list(i = 1L, j = 1L, x = 2, n = 1L)
  expect_equal(com_ordem("amd", sparse_chol(um))$logdet, log(2), tolerance = 1e-14)
  dois <- list(i = c(1L, 2L, 2L), j = c(1L, 1L, 2L), x = c(2, -1, 2), n = 2L)
  expect_equal(com_ordem("amd", sparse_chol(dois))$logdet, log(3), tolerance = 1e-14)
  so_diag <- list(i = 1:50, j = 1:50, x = seq(1, 5, length.out = 50), n = 50L)
  f <- com_ordem("amd", sparse_chol(so_diag))
  expect_equal(sort(f$perm), 1:50)
  expect_equal(f$logdet, sum(log(so_diag$x)), tolerance = 1e-13)
  # duas grades soltas e nos isolados
  g1 <- grade(15); g2 <- grade(12)
  a <- list(i = c(g1$i, g1$n + g2$i, g1$n + g2$n + 1:7), j = c(g1$j, g1$n + g2$j, g1$n + g2$n + 1:7),
            x = c(g1$x, g2$x, rep(3, 7)), n = g1$n + g2$n + 7L)
  fa <- com_ordem("amd", sparse_chol(a)); fg <- com_ordem("gm", sparse_chol(a))
  expect_equal(sort(fa$perm), seq_len(a$n))
  expect_equal(fa$logdet, fg$logdet, tolerance = 1e-12)
})

test_that("(5) model() e model_mt() com a AMD: o mesmo ajuste que com o grau minimo exato", {
  ped <- ped_sim(3000, 500, 12, 150, seed = 21)
  set.seed(22)
  d <- data.frame(id = ped$id, cg = sprintf("c%02d", sample(30, nrow(ped), TRUE)), stringsAsFactors = FALSE)
  # sinal genetico (efeitos de pai, correlacao genetica longe de 1): sem ele a covariancia fica na fronteira,
  # onde a verossimilhanca e plana e o ponto de parada depende da ordem; aqui o otimo e interior
  pai <- stats::rnorm(length(unique(ped$sire)))[match(ped$sire, unique(ped$sire))]
  pai2 <- stats::rnorm(length(unique(ped$sire)))[match(ped$sire, unique(ped$sire))]
  d$y1 <- 10 + pai + stats::rnorm(nrow(d)); d$y2 <- 5 + 0.5 * pai + pai2 + 0.3 * d$y1 + stats::rnorm(nrow(d))
  fg <- com_ordem("gm", model(y1 ~ cg + animal(id), d, ped, verbose = FALSE))
  fa <- com_ordem("amd", model(y1 ~ cg + animal(id), d, ped, verbose = FALSE))
  expect_true(fg$converged); expect_true(fa$converged)
  expect_equal(fa$neg2logl, fg$neg2logl, tolerance = 1e-10)
  expect_equal(fa$theta, fg$theta, tolerance = 1e-6)
  expect_equal(ebv(fa), ebv(fg), tolerance = 1e-6)
  expect_equal(fa$pev, fg$pev, tolerance = 1e-6)
  # AMD com 1 e 3 threads: bit a bit consigo
  fa3 <- com_ordem("amd", com_threads(3, model(y1 ~ cg + animal(id), d, ped, verbose = FALSE)))
  fa1 <- com_ordem("amd", com_threads(1, model(y1 ~ cg + animal(id), d, ped, verbose = FALSE)))
  expect_identical(fa3$neg2logl, fa1$neg2logl)
  expect_identical(fa3$theta, fa1$theta)
  mg <- com_ordem("gm", model_mt(cbind(y1, y2) ~ cg + animal(id), d, ped, verbose = FALSE))
  ma <- com_ordem("amd", model_mt(cbind(y1, y2) ~ cg + animal(id), d, ped, verbose = FALSE))
  expect_equal(ma$neg2logl, mg$neg2logl, tolerance = 1e-10)
  expect_equal(ma$theta, mg$theta, tolerance = 1e-6)
})
