// Multi-caracteristica (Henderson e Quaas, 1976) com registros completos: R0 cheia t x t, grupos com a caracteristica
// como coeficiente extra.
//
// ## O desenho
//
// A caracteristica NAO ganha um caminho proprio: ela vira mais uma dimensao de coeficiente.
// Um termo com n_coef coeficientes e t caracteristicas tem dim = n_coef * t no grupo, a
// covariancia C_g e (n_coef t) x (n_coef t) cheia, o que inclui as covariancias geneticas
// ENTRE caracteristicas, e a penalidade continua kron(C_g^-1, K^-1), sem caso especial.
// A convencao de coluna continua a mesma: coluna = coef_multi * n_niveis + nivel, com
// coef_multi = tau * n_coef + c (caracteristica major).
//
// ## As unidades
//
// O caminho de uma caracteristica escala as MME por s2e. Aqui nao ha um s2e: ha R0. Tudo e
// montado em unidades ABSOLUTAS:
//
//   C = W' (I (x) R0^-1) W + Pen        rhs = W' (I (x) R0^-1) y
//   -2logL = n log|R0| + soma_g [nl log|C_g| - dim log|K^-1|] + log|C| + (y'R^-1y - b'rhs)
//
// O log|C| em unidades absolutas ja carrega o log|X'V^-1X| do REML, nao ha termo (n-p)
// separado. A forma V densa de referencia computa a MESMA quantidade por um caminho sem
// nada em comum: log|V| + log|X'V^-1X| + y'Py.
//
// ## Registros
//
// Registro com QUALQUER caracteristica ausente sai inteiro, e a contagem e reportada. E a
// limitacao declarada desta primeira versao: padroes de ausencia por caracteristica exigem
// R0 condicionais por padrao, e entrar nisso sem gates proprios seria fingir capacidade.

#include "mme.h"

namespace br {

// vech: indice do elemento (i,j), i >= j, numa matriz d x d guardada por colunas.
static std::size_t vech_idx(std::size_t i, std::size_t j, std::size_t d) {
  if (i < j) std::swap(i, j);
  return (j * (2 * d - j + 1)) / 2 + (i - j);
}

// R0 lida do fim de theta.
static Densa r0_de(const std::vector<double>& theta, std::size_t off, std::size_t t) {
  Densa r(t, t);
  for (std::size_t j = 0; j < t; j++)
    for (std::size_t i = j; i < t; i++) {
      const double v = theta[off + vech_idx(i, j, t)];
      r.at(i, j) = v;
      r.at(j, i) = v;
    }
  return r;
}

// MontadoMT vive em mme.h (o nivel 2 da curva de crescimento, curva.cpp, monta por aqui)

// Colunas empilhadas de um registro: (coluna_global, caracteristica, valor).
struct EntradaLinha {
  std::uint32_t col;
  std::uint32_t trait;
  double val;
};

MontadoMT monta_mme_mt(const DesenhoMT& d, const std::vector<double>& theta,
                       const std::vector<Densa>* fatores) {
  MontadoMT M;
  const std::size_t t = d.t;
  const std::size_t off_r0 = d.modelo.offset_residual;
  M.r0 = r0_de(theta, off_r0, t);
  {
    Densa l = M.r0;
    if (!chol_densa(l)) return M;          // R0 inadmissivel
  }
  M.r0inv = inv_pd(M.r0);

  // mascaras por linha e cache por padrao
  M.mask.assign(d.nlin, 0);
  M.rinv_mask.assign(std::size_t(1) << t, Densa());
  M.ld_mask.assign(std::size_t(1) << t, 0.0);
  M.logdet_r = 0.0;
  for (std::size_t r = 0; r < d.nlin; r++) {
    if (!d.usa[r]) continue;
    std::uint32_t mk = 0;
    for (std::size_t tau = 0; tau < t; tau++)
      if (d.obs.empty() || d.obs[r * t + tau]) mk |= (1u << tau);
    M.mask[r] = mk;
    if (M.rinv_mask[mk].nlin == 0) {
      std::vector<std::size_t> quais;
      for (std::size_t tau = 0; tau < t; tau++)
        if (mk & (1u << tau)) quais.push_back(tau);
      Densa sub(quais.size(), quais.size());
      for (std::size_t a = 0; a < quais.size(); a++)
        for (std::size_t b = 0; b < quais.size(); b++)
          sub.at(a, b) = M.r0.at(quais[a], quais[b]);
      Densa subi = inv_pd(sub);
      Densa emb(t, t);
      for (std::size_t a = 0; a < quais.size(); a++)
        for (std::size_t b = 0; b < quais.size(); b++)
          emb.at(quais[a], quais[b]) = subi.at(a, b);
      M.rinv_mask[mk] = emb;
      M.ld_mask[mk] = logdet_pd(sub);
    }
    M.logdet_r += M.ld_mask[mk];
  }

  M.n_fixo = d.n_fixa;
  M.offset_grupo.resize(d.modelo.grupos.size());
  std::size_t acc = M.n_fixo;
  for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
    M.offset_grupo[g] = acc;
    acc += d.largura(g);
  }
  M.total = acc;

  // slots: (indice em aleatorios, primeira coluna global do termo) por grupo
  std::vector<std::vector<std::pair<std::size_t, std::size_t>>> slots(d.modelo.grupos.size());
  for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
    std::size_t col = M.offset_grupo[g];
    for (std::size_t tm : d.modelo.grupos[g].termos)
      for (std::size_t a = 0; a < d.aleatorios.size(); a++)
        if (d.aleatorios[a].termo == tm) {
          slots[g].push_back({a, col});
          col += d.aleatorios[a].z.ncol * d.aleatorios[a].nt(t);
        }
  }

  std::vector<std::vector<std::vector<std::pair<std::uint32_t, double>>>> zl;
  for (const DesenhoTermo& a : d.aleatorios) {
    std::vector<std::vector<std::pair<std::uint32_t, double>>> lin(a.z.nlin);
    for (std::size_t c = 0; c < a.z.ncol; c++)
      for (std::size_t k = a.z.colptr[c]; k < a.z.colptr[c + 1]; k++)
        lin[a.z.linha[k]].push_back({static_cast<std::uint32_t>(c), a.z.valor[k]});
    zl.push_back(std::move(lin));
  }

  std::vector<std::uint32_t> ti, tj;
  std::vector<double> tv;
  M.rhs.assign(M.total, 0.0);
  std::vector<EntradaLinha> lin;

  for (std::size_t r = 0; r < d.nlin; r++) {
    if (!d.usa[r]) continue;
    lin.clear();
    // X: coluna j da caracteristica tau fica em j*t + tau
    for (std::size_t ex = d.x.linptr[r]; ex < d.x.linptr[r + 1]; ex++) {
      const std::size_t j = d.x.coluna[ex];
      const double v = d.x.valor[ex];
      if (v == 0.0) continue;
      for (std::size_t tau = 0; tau < t; tau++) {
        const int e = d.eq_fixa[j * t + tau];
        if (e < 0) continue;
        lin.push_back({static_cast<std::uint32_t>(e),
                       static_cast<std::uint32_t>(tau), v});
      }
    }
    // Z: coluna c do termo (coef ct, nivel nv) na caracteristica tau vai para
    // col0 + (tau * n_coef + ct) * n_niveis + nv , caracteristica major no coeficiente
    for (std::size_t g = 0; g < slots.size(); g++)
      for (const auto& [a, col0] : slots[g]) {
        const DesenhoTermo& dt = d.aleatorios[a];
        for (const auto& [c, v] : zl[a][r]) {
          const std::size_t ct = c / dt.n_niveis;
          const std::size_t nv = c % dt.n_niveis;
          for (std::size_t pt = 0; pt < dt.nt(t); pt++)
            lin.push_back({static_cast<std::uint32_t>(
                               col0 + (pt * dt.n_coef + ct) * dt.n_niveis + nv),
                           static_cast<std::uint32_t>(dt.traco(pt)), v});
        }
      }

    // contribuicao do registro, com a inversa EMBUTIDA do padrao: caracteristica nao
    // observada tem linha e coluna zero, entao a contribuicao dela anula sozinha
    const Densa& rir = M.rinv_mask[M.mask[r]];
    for (std::size_t p = 0; p < lin.size(); p++) {
      double ry = 0.0;
      for (std::size_t tau = 0; tau < t; tau++)
        ry += rir.at(lin[p].trait, tau) * d.y.at(r, tau);
      M.rhs[lin[p].col] += lin[p].val * ry;
      for (std::size_t q = p; q < lin.size(); q++) {
        std::uint32_t i = lin[p].col, j = lin[q].col;
        double v = lin[p].val * rir.at(lin[p].trait, lin[q].trait) * lin[q].val;
        if (i == j && p != q) v *= 2.0;   // (p,q) e (q,p) caem na mesma posicao da diagonal
        if (i > j) std::swap(i, j);
        ti.push_back(i);
        tj.push_back(j);
        tv.push_back(v);
      }
    }
  }

  // o sistema transformado (BlocoT, mme.h): T' C_d T nos triplos dos dados, antes da penalidade
  if (fatores) {
    for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
      if ((*fatores)[g].dados.empty()) continue;
      std::size_t nl = 0;
      for (const DesenhoTermo& a : d.aleatorios)
        if (a.termo == d.modelo.grupos[g].termos[0]) { nl = a.n_niveis; break; }
      M.blocos_t.push_back({M.offset_grupo[g], d.modelo.grupos[g].dim, nl, (*fatores)[g]});
    }
    transforma_triplos(M.blocos_t, ti, tj, tv, M.rhs);
  }

  // penalidade por grupo, unidades absolutas: kron(C_g^-1, K^-1)
  for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
    const Grupo& gr = d.modelo.grupos[g];
    const std::size_t dim = gr.dim;   // ja inclui o fator t
    Densa cg(dim, dim);
    for (std::size_t j = 0; j < dim; j++)
      for (std::size_t i = j; i < dim; i++) {
        const double v = theta[gr.offset + vech_idx(i, j, dim)];
        cg.at(i, j) = v;
        cg.at(j, i) = v;
      }
    // no sistema transformado a penalidade do grupo e I (x) K^-1, e log|I| = 0
    const bool trans = fatores && !(*fatores)[g].dados.empty();
    Densa cinv(dim, dim);
    double ld = 0.0;
    if (trans) {
      for (std::size_t i = 0; i < dim; i++) cinv.at(i, i) = 1.0;
    } else {
      Densa l = cg;
      if (!chol_densa(l)) return M;
      cinv = inv_pd(cg);
      ld = logdet_pd(cg);
    }

    std::size_t nl = 0;
    for (const DesenhoTermo& a : d.aleatorios)
      if (a.termo == gr.termos[0]) { nl = a.n_niveis; break; }
    M.logdet_g += static_cast<double>(nl) * ld
                - static_cast<double>(dim) * d.kinv_logdet[g] / 1.0;

    const std::size_t off = M.offset_grupo[g];
    const bool com_k = d.kinv[g].ncol > 0;
    auto poe = [&](std::size_t gi, std::size_t gj, double val) {
      if (gi > gj) return;
      ti.push_back(static_cast<std::uint32_t>(gi));
      tj.push_back(static_cast<std::uint32_t>(gj));
      tv.push_back(val);
    };
    for (std::size_t a = 0; a < dim; a++)
      for (std::size_t b = 0; b < dim; b++) {
        const double f = cinv.at(a, b);
        // f == 0 NAO pula mais: a covariancia que comeca em zero ganharia um padrao
        // menor na primeira avaliacao, e o cache da simbolica congelaria esse padrao
        // errado para o ajuste inteiro. Zero explicito ocupa o slot e nao muda numero.
        if (com_k) {
          const Csc& k = d.kinv[g];
          for (std::size_t col = 0; col < k.ncol; col++)
            for (std::size_t p = k.colptr[col]; p < k.colptr[col + 1]; p++) {
              const std::size_t rk = k.linha[p];
              poe(off + a * nl + rk, off + b * nl + col, f * k.valor[p]);
              if (rk != col) poe(off + a * nl + col, off + b * nl + rk, f * k.valor[p]);
            }
        } else {
          for (std::size_t l = 0; l < nl; l++)
            poe(off + a * nl + l, off + b * nl + l, f);
        }
      }
  }

  M.c = de_triplos(M.total, M.total, ti, tj, tv);
  M.ok = true;
  return M;
}

// A nota da diagonal duplicada acima: com i == j vindos de colunas distintas p != q do mesmo
// registro (acontece quando a mesma coluna aparece para duas caracteristicas, nao acontece
// aqui porque a coluna carrega a caracteristica, mas a guarda fica) o par (p,q) representa
// as duas ordens. Mantida por seguranca; o gate da forma V pega qualquer dupla contagem.

// e^ = y - W b por registro e caracteristica
static Densa residuo_mt(const DesenhoMT& d, const MontadoMT& M, const std::vector<double>& sol) {
  const std::size_t t = d.t;
  Densa ehat(d.nlin, t);
  for (std::size_t r = 0; r < d.nlin; r++) {
    if (!d.usa[r]) continue;
    for (std::size_t tau = 0; tau < t; tau++) {
      double wb = 0.0;
      for (std::size_t ex = d.x.linptr[r]; ex < d.x.linptr[r + 1]; ex++) {
        const int e = d.eq_fixa[d.x.coluna[ex] * t + tau];
        if (e >= 0) wb += d.x.valor[ex] * sol[e];
      }
      ehat.at(r, tau) = d.y.at(r, tau) - wb;
    }
  }
  for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
    std::size_t col0 = M.offset_grupo[g];
    for (std::size_t tm : d.modelo.grupos[g].termos)
      for (const DesenhoTermo& a : d.aleatorios)
        if (a.termo == tm) {
          for (std::size_t c = 0; c < a.z.ncol; c++) {
            const std::size_t ct = c / a.n_niveis, nv = c % a.n_niveis;
            for (std::size_t pt = 0; pt < a.nt(t); pt++) {
              const std::size_t tau = a.traco(pt);
              const double bc = sol[col0 + (pt * a.n_coef + ct) * a.n_niveis + nv];
              if (bc == 0.0) continue;
              for (std::size_t p = a.z.colptr[c]; p < a.z.colptr[c + 1]; p++)
                if (d.usa[a.z.linha[p]])
                  ehat.at(a.z.linha[p], tau) -= a.z.valor[p] * bc;
            }
          }
          col0 += a.z.ncol * a.nt(t);
        }
  }
  return ehat;
}

// O -2logL do model_mt() pelo sistema transformado (o racional em BlocoT, src/mme.h), com a
// forma quadratica pelos residuos, como em avalia_mt()
static bool neg2logl_transformada_mt(const DesenhoMT& d, const std::vector<double>& theta,
                                     const std::vector<Densa>& F, double& valor) {
  MontadoMT M = monta_mme_mt(d, theta, &F);
  if (!M.ok) return false;
  const std::vector<std::size_t> perm = grau_minimo(M.c);
  const Csc pc = permuta_sim(M.c, perm);
  const Simbolica sb = simbolica(pc);
  Csc L;
  if (!cholesky(pc, sb, L)) return false;
  std::vector<double> pb(M.total);
  for (std::size_t k = 0; k < M.total; k++) pb[k] = M.rhs[perm[k]];
  const std::vector<double> px = resolve(L, pb);
  std::vector<double> sol(M.total, 0.0);
  for (std::size_t k = 0; k < M.total; k++) sol[perm[k]] = px[k];
  double upen = 0.0;
  if (!upen_transformada(d.modelo, d.aleatorios, d.kinv, theta, M.offset_grupo, F, sol, 1.0, upen))
    return false;
  destransforma_solucao(M.blocos_t, sol);
  const Densa ehat = residuo_mt(d, M, sol);
  double ere = 0.0;
  for (std::size_t r = 0; r < d.nlin; r++) {
    if (!d.usa[r]) continue;
    const Densa& rir = M.rinv_mask[M.mask[r]];
    for (std::size_t a = 0; a < d.t; a++)
      for (std::size_t b = 0; b < d.t; b++) ere += ehat.at(r, a) * rir.at(a, b) * ehat.at(r, b);
  }
  valor = M.logdet_r + M.logdet_g + logdet(L) + ere + upen + d.logdet_escala_x;
  return std::isfinite(valor);
}

AvaliacaoMT avalia_mt(const DesenhoMT& d, const std::vector<double>& theta,
                      CacheSimbolica* cache) {
  AvaliacaoMT A;
  MontadoMT M = monta_mme_mt(d, theta);
  if (!M.ok) { A.motivo = 1; return A; }
  const std::size_t t = d.t;

  CacheSimbolica local;
  CacheSimbolica* cs = cache ? cache : &local;
  if (!cs->pronto) {
    cs->perm = grau_minimo(M.c);
    cs->sb = simbolica(permuta_sim(M.c, cs->perm));
    cs->pronto = true;
  }
  const std::vector<std::size_t>& perm = cs->perm;
  const Simbolica& sb = cs->sb;
  const Csc& pc = permuta_cache(M.c, *cs);
  Csc L;
  if (!cholesky(pc, sb, L)) { A.motivo = 2; return A; }

  std::vector<double> pb(M.total);
  for (std::size_t k = 0; k < M.total; k++) pb[k] = M.rhs[perm[k]];
  std::vector<double> px = resolve(L, pb);
  A.solucao.assign(M.total, 0.0);
  for (std::size_t k = 0; k < M.total; k++) A.solucao[perm[k]] = px[k];

  // -2logL: os termos que nao dependem de y aqui; a forma quadratica sai dos residuos,
  // e'R^-1 e + u'G^-1 u, e nao de y'R^-1 y - b'r, que perde digitos com uma constante
  // somada a y (o racional esta em avalia(), src/aireml.cpp).
  const double parte_fixa = M.logdet_r + M.logdet_g + logdet(L);
  double upen = 0.0;   // soma_g soma_ab C_g^-1[a,b] Q[a,b] = u'G^-1 u

  // residuo por registro e caracteristica: e^ = y - W b
  const Densa ehat = residuo_mt(d, M, A.solucao);

  // inversa seletiva na numeracao original
  SelInv z = inversa_seletiva(L, 0);
  std::vector<std::size_t> inv_perm(M.total);
  for (std::size_t k = 0; k < M.total; k++) inv_perm[perm[k]] = k;
  auto z_orig = [&](std::size_t i, std::size_t j, double& out) {
    return z.get(inv_perm[i], inv_perm[j], out);
  };

  // PEV: a diagonal de C^-1 na numeracao original. Este caminho ja esta em unidades
  // absolutas, entao e a diagonal direta, sem o s2e do caminho uni.
  {
    std::vector<double> dz = z.diagonal();
    A.pev.assign(M.total, std::nan(""));
    for (std::size_t k = 0; k < M.total; k++) A.pev[perm[k]] = dz[k];
  }

  const std::size_t ntheta = d.modelo.ntheta;
  A.score.assign(ntheta, 0.0);
  // a proposta EM comeca IGUAL ao theta corrente: os grupos serao sobrescritos abaixo e o
  // R0 fica onde esta, porque o passo M do residuo multivariado ainda nao existe aqui. A
  // proposta move so o que sabe mover, e quem chama confere se ela melhora.
  A.em_theta = theta;
  A.ai = Densa(ntheta, ntheta);

  // ---------------- score dos grupos: nl C^-1 - C^-1 (Q + T) C^-1, unidades absolutas
  std::vector<Densa> cinv_g(d.modelo.grupos.size());
  std::vector<std::vector<std::vector<double>>> vgrupo(d.modelo.grupos.size());
  for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
    const Grupo& gr = d.modelo.grupos[g];
    const std::size_t dim = gr.dim;
    const std::size_t off = M.offset_grupo[g];
    std::size_t nl = 0;
    for (const DesenhoTermo& a : d.aleatorios)
      if (a.termo == gr.termos[0]) { nl = a.n_niveis; break; }

    Densa cg(dim, dim);
    for (std::size_t j = 0; j < dim; j++)
      for (std::size_t i = j; i < dim; i++) {
        const double v = theta[gr.offset + vech_idx(i, j, dim)];
        cg.at(i, j) = v;
        cg.at(j, i) = v;
      }
    Densa cinv = inv_pd(cg);
    cinv_g[g] = cinv;

    // Q[a,b] = u_a' K^-1 u_b
    Densa Q(dim, dim), T(dim, dim);
    std::vector<std::vector<double>> ku(dim);
    for (std::size_t b = 0; b < dim; b++) {
      ku[b].assign(nl, 0.0);
      const double* u = A.solucao.data() + off + b * nl;
      if (d.kinv[g].ncol > 0) {
        const Csc& k = d.kinv[g];
        for (std::size_t c = 0; c < nl; c++)
          for (std::size_t p = k.colptr[c]; p < k.colptr[c + 1]; p++) {
            const std::size_t r = k.linha[p];
            ku[b][r] += k.valor[p] * u[c];
            if (r != c) ku[b][c] += k.valor[p] * u[r];
          }
      } else {
        for (std::size_t l = 0; l < nl; l++) ku[b][l] = u[l];
      }
    }
    for (std::size_t a = 0; a < dim; a++)
      for (std::size_t b = 0; b < dim; b++) {
        double s = 0.0;
        const double* ua = A.solucao.data() + off + a * nl;
        for (std::size_t l = 0; l < nl; l++) s += ua[l] * ku[b][l];
        Q.at(a, b) = s;
      }
    for (std::size_t a = 0; a < dim; a++)
      for (std::size_t b = 0; b < dim; b++) upen += cinv.at(a, b) * Q.at(a, b);

    // T[a,b] = tr(K^-1 [C^-1]_{a,b}), unidades absolutas, sem fator s2e
    for (std::size_t a = 0; a < dim; a++)
      for (std::size_t b = 0; b < dim; b++) {
        double s = 0.0;
        if (d.kinv[g].ncol > 0) {
          const Csc& k = d.kinv[g];
          for (std::size_t c = 0; c < nl; c++)
            for (std::size_t p = k.colptr[c]; p < k.colptr[c + 1]; p++) {
              const std::size_t r = k.linha[p];
              double zv;
              if (z_orig(off + a * nl + r, off + b * nl + c, zv)) s += k.valor[p] * zv;
              else A.fora_do_padrao++;
              if (r != c) {
                if (z_orig(off + a * nl + c, off + b * nl + r, zv)) s += k.valor[p] * zv;
                else A.fora_do_padrao++;
              }
            }
        } else {
          for (std::size_t l = 0; l < nl; l++) {
            double zv;
            if (z_orig(off + a * nl + l, off + b * nl + l, zv)) s += zv;
            else A.fora_do_padrao++;
          }
        }
        T.at(a, b) = s;
      }

    Densa QT(dim, dim);
    for (std::size_t a = 0; a < dim; a++)
      for (std::size_t b = 0; b < dim; b++) QT.at(a, b) = Q.at(a, b) + T.at(a, b);
    Densa CiQT(dim, dim), Mg(dim, dim);
    for (std::size_t a = 0; a < dim; a++)
      for (std::size_t b = 0; b < dim; b++) {
        double s = 0.0;
        for (std::size_t k = 0; k < dim; k++) s += cinv.at(a, k) * QT.at(k, b);
        CiQT.at(a, b) = s;
      }
    for (std::size_t a = 0; a < dim; a++)
      for (std::size_t b = 0; b < dim; b++) {
        double s = 0.0;
        for (std::size_t k = 0; k < dim; k++) s += CiQT.at(a, k) * cinv.at(k, b);
        Mg.at(a, b) = static_cast<double>(nl) * cinv.at(a, b) - s;
      }
    // EM do grupo: C_g <- (Q + T) / nl. QT ja esta montado logo acima para o score.
    for (std::size_t j = 0; j < dim; j++)
      for (std::size_t i = j; i < dim; i++)
        A.em_theta[gr.offset + vech_idx(i, j, dim)] =
            (QT.at(i, j) + QT.at(j, i)) / (2.0 * static_cast<double>(nl));
    for (std::size_t j = 0; j < dim; j++)
      for (std::size_t i = j; i < dim; i++) {
        const std::size_t k = gr.offset + vech_idx(i, j, dim);
        A.score[k] = (i == j) ? Mg.at(i, i) : Mg.at(i, j) + Mg.at(j, i);
      }

    // v_s = soma_d C^-1[s,d] u_d, para a AI
    vgrupo[g].assign(dim, std::vector<double>(nl, 0.0));
    for (std::size_t s2 = 0; s2 < dim; s2++)
      for (std::size_t dd = 0; dd < dim; dd++) {
        const double cq = cinv.at(s2, dd);
        if (cq == 0.0) continue;
        const double* ud = A.solucao.data() + off + dd * nl;
        for (std::size_t l = 0; l < nl; l++) vgrupo[g][s2][l] += cq * ud[l];
      }
  }

  // ---------------- score de R0[a,b]:
  //   tr(P dV) - y'P dV P y,  dV = I_n (x) E_ab(simetrica)
  //   tr(P dV) = n tr(R0^-1 E) - soma_registros tr(C^-1 W_i' R0^-1 E R0^-1 W_i)
  //   y'P dV P y = soma_i (R0^-1 e_i)' E (R0^-1 e_i)
  //
  // O traco por registro le a inversa seletiva nas posicoes do clique do registro, que
  // estao no padrao do fator porque o clique esta em C. Fora do padrao e CONTADO.
  {
    const std::size_t off_r0 = d.modelo.offset_residual;
    const double n_reg = static_cast<double>(d.n_usadas());

    // R0^-1 e_i por registro
    Densa re(d.nlin, t);
    for (std::size_t r = 0; r < d.nlin; r++) {
      if (!d.usa[r]) continue;
      const Densa& rir = M.rinv_mask[M.mask[r]];
      for (std::size_t a = 0; a < t; a++) {
        double s = 0.0;
        for (std::size_t b = 0; b < t; b++) s += rir.at(a, b) * ehat.at(r, b);
        re.at(r, a) = s;
      }
    }
    // e'R^-1 e: a celula ausente tem linha e coluna zero na inversa embutida e nao conta
    double ere = 0.0;
    for (std::size_t r = 0; r < d.nlin; r++) {
      if (!d.usa[r]) continue;
      for (std::size_t a = 0; a < t; a++) ere += ehat.at(r, a) * re.at(r, a);
    }
    // + o termo da padronizacao das covariaveis (fixos.cpp)
    A.neg2logl = parte_fixa + ere + upen + d.logdet_escala_x;
    // perto da singularidade de um grupo, o VALOR vem do sistema transformado (BlocoT,
    // mme.h); o score e a AI seguem das MME, eles so dao a direcao
    {
      std::vector<Densa> cgs(d.modelo.grupos.size()), F;
      for (std::size_t g = 0; g < cgs.size(); g++) cgs[g] = cov_grupo(d.modelo, theta, g);
      double v = 0.0;
      if (fatores_quase_singulares(cgs, F) && neg2logl_transformada_mt(d, theta, F, v))
        A.neg2logl = v;
    }

    // clique de cada registro (colunas e caracteristica), reconstruido como na montagem
    std::vector<std::vector<std::pair<std::uint32_t, double>>> zl2;
    for (const DesenhoTermo& a : d.aleatorios) {
      (void) a;
      zl2.push_back({});
    }

    (void) n_reg;
    for (std::size_t ja = 0; ja < t; ja++)
      for (std::size_t jb = ja; jb < t; jb++) {
        const std::size_t k = off_r0 + vech_idx(jb, ja, t);
        // tr(R^-1 dV) = soma por registro de tr(R0p^-1 E), so os registros cujo padrao
        // OBSERVA as duas caracteristicas contribuem, e a inversa embutida faz isso sozinha
        double tr_rinv = 0.0;
        double quad = 0.0;
        for (std::size_t r = 0; r < d.nlin; r++) {
          if (!d.usa[r]) continue;
          const Densa& rir = M.rinv_mask[M.mask[r]];
          tr_rinv += (ja == jb) ? rir.at(ja, ja) : rir.at(ja, jb) + rir.at(jb, ja);
          quad += (ja == jb) ? re.at(r, ja) * re.at(r, ja)
                             : 2.0 * re.at(r, ja) * re.at(r, jb);
        }
        A.score[k] = tr_rinv - quad;   // o traco do meio entra abaixo
      }

    // termo do meio: + soma_i tr(C^-1 W_i' R0^-1 E R0^-1 W_i), com sinal NEGATIVO no score
    // (tr(P dV) = tr(R^-1 dV) - esse traco). Percorre os cliques por registro.
    std::vector<EntradaLinha> lin;
    // reconstrucao identica a da montagem
    std::vector<std::vector<std::vector<std::pair<std::uint32_t, double>>>> zl3;
    for (const DesenhoTermo& a : d.aleatorios) {
      std::vector<std::vector<std::pair<std::uint32_t, double>>> l2(a.z.nlin);
      for (std::size_t c = 0; c < a.z.ncol; c++)
        for (std::size_t kk = a.z.colptr[c]; kk < a.z.colptr[c + 1]; kk++)
          l2[a.z.linha[kk]].push_back({static_cast<std::uint32_t>(c), a.z.valor[kk]});
      zl3.push_back(std::move(l2));
    }
    std::vector<std::vector<std::pair<std::size_t, std::size_t>>> slots(d.modelo.grupos.size());
    for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
      std::size_t col = M.offset_grupo[g];
      for (std::size_t tm : d.modelo.grupos[g].termos)
        for (std::size_t a = 0; a < d.aleatorios.size(); a++)
          if (d.aleatorios[a].termo == tm) {
            slots[g].push_back({a, col});
            col += d.aleatorios[a].z.ncol * d.aleatorios[a].nt(t);
          }
    }

    for (std::size_t r = 0; r < d.nlin; r++) {
      if (!d.usa[r]) continue;
      lin.clear();
      for (std::size_t ex = d.x.linptr[r]; ex < d.x.linptr[r + 1]; ex++) {
        const std::size_t j = d.x.coluna[ex];
        const double v = d.x.valor[ex];
        if (v == 0.0) continue;
        for (std::size_t tau = 0; tau < t; tau++) {
          const int e = d.eq_fixa[j * t + tau];
          if (e < 0) continue;
          lin.push_back({static_cast<std::uint32_t>(e),
                         static_cast<std::uint32_t>(tau), v});
        }
      }
      for (std::size_t g = 0; g < slots.size(); g++)
        for (const auto& [a, col0] : slots[g]) {
          const DesenhoTermo& dt = d.aleatorios[a];
          for (const auto& [c, v] : zl3[a][r]) {
            const std::size_t ct = c / dt.n_niveis, nv = c % dt.n_niveis;
            for (std::size_t pt = 0; pt < dt.nt(t); pt++)
              lin.push_back({static_cast<std::uint32_t>(
                                 col0 + (pt * dt.n_coef + ct) * dt.n_niveis + nv),
                             static_cast<std::uint32_t>(dt.traco(pt)), v});
          }
        }

      // para cada par de entradas do clique, o peso R0^-1 E R0^-1 entre as caracteristicas
      for (std::size_t p = 0; p < lin.size(); p++)
        for (std::size_t q = 0; q < lin.size(); q++) {
          double zv;
          if (!z_orig(lin[p].col, lin[q].col, zv)) { A.fora_do_padrao++; continue; }
          const double w = lin[p].val * lin[q].val * zv;
          const Densa& rir2 = M.rinv_mask[M.mask[r]];
          for (std::size_t ja = 0; ja < t; ja++)
            for (std::size_t jb = ja; jb < t; jb++) {
              const std::size_t k = off_r0 + vech_idx(jb, ja, t);
              // [R0p^-1 E R0p^-1]_{trait_p, trait_q}, do PADRAO do registro
              double m2;
              if (ja == jb) m2 = rir2.at(lin[p].trait, ja) * rir2.at(ja, lin[q].trait);
              else m2 = rir2.at(lin[p].trait, ja) * rir2.at(jb, lin[q].trait)
                      + rir2.at(lin[p].trait, jb) * rir2.at(ja, lin[q].trait);
              // SUBTRAIDO: tr(P dV) = tr(R^-1 dV) - tr(C^-1 W'R^-1 dV R^-1 W). A primeira
              // versao somava, e o gate de diferencas finitas apontou os TRES parametros de
              // R0 errados por exatamente 2x este traco, o comentario dizia o sinal certo
              // e o codigo fazia o errado. E o motivo de o gate cobrir TODOS os parametros.
              A.score[k] -= w * m2;
            }
        }
    }
  }

  // ---------------- AI por diferencas do score? NAO: a AI e a informacao media exata,
  // AI[i,j] = f_i' P f_j, com f_k = dV/dtheta_k Py. Py = R^-1 e^ (empilhado).
  // Para grupo: f_k tem a mesma forma do caso uni (Z aplicado a dC C^-1 u), por
  // caracteristica-coeficiente. Para R0: f_k por registro = E R^-1 e_i.
  {
    const std::size_t off_r0 = d.modelo.offset_residual;
    const std::size_t nlin = d.nlin;
    // f por parametro, empilhado registro x caracteristica
    std::vector<Densa> f(ntheta, Densa(nlin, t));

    // grupos
    for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
      const Grupo& gr = d.modelo.grupos[g];
      const std::size_t dim = gr.dim;
      std::size_t nl = 0;
      for (const DesenhoTermo& a : d.aleatorios)
        if (a.termo == gr.termos[0]) { nl = a.n_niveis; break; }

      for (std::size_t j = 0; j < dim; j++)
        for (std::size_t i = j; i < dim; i++) {
          const std::size_t k = gr.offset + vech_idx(i, j, dim);
          std::vector<std::vector<double>> w(dim, std::vector<double>(nl, 0.0));
          if (i == j) w[i] = vgrupo[g][i];
          else { w[i] = vgrupo[g][j]; w[j] = vgrupo[g][i]; }
          // f_k = soma_slots Z_slot w_slot, espalhado por caracteristica do slot
          std::size_t slot = 0;
          for (std::size_t tm : d.modelo.grupos[g].termos)
            for (const DesenhoTermo& a : d.aleatorios)
              if (a.termo == tm) {
                for (std::size_t pt = 0; pt < a.nt(t); pt++)
                  for (std::size_t ct = 0; ct < a.n_coef; ct++) {
                    const std::size_t tau = a.traco(pt);
                    const std::size_t s2 = pt * a.n_coef + ct + slot;
                    const std::vector<double>& ws = w[s2];
                    for (std::size_t nv = 0; nv < nl; nv++) {
                      const double x = ws[nv];
                      if (x == 0.0) continue;
                      const std::size_t colz = ct * a.n_niveis + nv;
                      for (std::size_t p = a.z.colptr[colz]; p < a.z.colptr[colz + 1]; p++)
                        if (d.usa[a.z.linha[p]])
                          f[k].at(a.z.linha[p], tau) += a.z.valor[p] * x;
                    }
                  }
                slot += a.n_coef * a.nt(t);
              }
        }
    }

    // R0: f_k por registro = E (R0^-1 e_i)
    Densa re(nlin, t);
    for (std::size_t r = 0; r < nlin; r++) {
      if (!d.usa[r]) continue;
      const Densa& rir = M.rinv_mask[M.mask[r]];
      for (std::size_t a = 0; a < t; a++) {
        double s = 0.0;
        for (std::size_t b = 0; b < t; b++) s += rir.at(a, b) * ehat.at(r, b);
        re.at(r, a) = s;
      }
    }
    for (std::size_t ja = 0; ja < t; ja++)
      for (std::size_t jb = ja; jb < t; jb++) {
        const std::size_t k = off_r0 + vech_idx(jb, ja, t);
        for (std::size_t r = 0; r < nlin; r++) {
          if (!d.usa[r]) continue;
          if (ja == jb) f[k].at(r, ja) += re.at(r, ja);
          else { f[k].at(r, ja) += re.at(r, jb); f[k].at(r, jb) += re.at(r, ja); }
        }
      }

    // AI[i,j] = f_i' R^-1 f_j - (W'R^-1 f_i)' C^-1 (W'R^-1 f_j)
    // W'R^-1 f: mesma montagem por registro
    std::vector<std::vector<double>> wtf(ntheta, std::vector<double>(M.total, 0.0));
    {
      std::vector<EntradaLinha> lin;
      std::vector<std::vector<std::vector<std::pair<std::uint32_t, double>>>> zl4;
      for (const DesenhoTermo& a : d.aleatorios) {
        std::vector<std::vector<std::pair<std::uint32_t, double>>> l2(a.z.nlin);
        for (std::size_t c = 0; c < a.z.ncol; c++)
          for (std::size_t kk = a.z.colptr[c]; kk < a.z.colptr[c + 1]; kk++)
            l2[a.z.linha[kk]].push_back({static_cast<std::uint32_t>(c), a.z.valor[kk]});
        zl4.push_back(std::move(l2));
      }
      std::vector<std::vector<std::pair<std::size_t, std::size_t>>> slots(d.modelo.grupos.size());
      for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
        std::size_t col = M.offset_grupo[g];
        for (std::size_t tm : d.modelo.grupos[g].termos)
          for (std::size_t a = 0; a < d.aleatorios.size(); a++)
            if (d.aleatorios[a].termo == tm) {
              slots[g].push_back({a, col});
              col += d.aleatorios[a].z.ncol * d.aleatorios[a].nt(t);
            }
      }
      for (std::size_t r = 0; r < nlin; r++) {
        if (!d.usa[r]) continue;
        lin.clear();
        for (std::size_t ex = d.x.linptr[r]; ex < d.x.linptr[r + 1]; ex++) {
          const std::size_t j = d.x.coluna[ex];
          const double v = d.x.valor[ex];
          if (v == 0.0) continue;
          for (std::size_t tau = 0; tau < t; tau++) {
            const int e = d.eq_fixa[j * t + tau];
            if (e < 0) continue;
            lin.push_back({static_cast<std::uint32_t>(e),
                           static_cast<std::uint32_t>(tau), v});
          }
        }
        for (std::size_t g = 0; g < slots.size(); g++)
          for (const auto& [a, col0] : slots[g]) {
            const DesenhoTermo& dt = d.aleatorios[a];
            for (const auto& [c, v] : zl4[a][r]) {
              const std::size_t ct = c / dt.n_niveis, nv = c % dt.n_niveis;
              for (std::size_t pt = 0; pt < dt.nt(t); pt++)
                lin.push_back({static_cast<std::uint32_t>(
                                   col0 + (pt * dt.n_coef + ct) * dt.n_niveis + nv),
                               static_cast<std::uint32_t>(dt.traco(pt)), v});
            }
          }
        const Densa& rir3 = M.rinv_mask[M.mask[r]];
        for (std::size_t k = 0; k < ntheta; k++) {
          for (const EntradaLinha& e : lin) {
            double s = 0.0;
            for (std::size_t tau = 0; tau < t; tau++)
              s += rir3.at(e.trait, tau) * f[k].at(r, tau);
            wtf[k][e.col] += e.val * s;
          }
        }
      }
    }
    std::vector<std::vector<double>> cwf(ntheta);
    for (std::size_t k = 0; k < ntheta; k++) {
      std::vector<double> pb2(M.total);
      for (std::size_t q2 = 0; q2 < M.total; q2++) pb2[q2] = wtf[k][perm[q2]];
      std::vector<double> px2 = resolve(L, pb2);
      cwf[k].assign(M.total, 0.0);
      for (std::size_t q2 = 0; q2 < M.total; q2++) cwf[k][perm[q2]] = px2[q2];
    }
    for (std::size_t a = 0; a < ntheta; a++)
      for (std::size_t b = a; b < ntheta; b++) {
        double frf = 0.0;
        for (std::size_t r = 0; r < nlin; r++) {
          if (!d.usa[r]) continue;
          const Densa& rir4 = M.rinv_mask[M.mask[r]];
          for (std::size_t t1 = 0; t1 < t; t1++)
            for (std::size_t t2 = 0; t2 < t; t2++)
              frf += f[a].at(r, t1) * rir4.at(t1, t2) * f[b].at(r, t2);
        }
        double wcw = 0.0;
        for (std::size_t q2 = 0; q2 < M.total; q2++) wcw += wtf[a][q2] * cwf[b][q2];
        const double v = frf - wcw;
        A.ai.at(a, b) = v;
        A.ai.at(b, a) = v;
      }
  }

  A.ok = true;
  return A;
}

// -2logL pela forma V densa, multi-caracteristica. Caminho independente, so para os gates.
double neg2logl_densa_V_mt(const DesenhoMT& d, const std::vector<double>& theta) {
  const std::size_t t = d.t;
  Densa r0 = r0_de(theta, d.modelo.offset_residual, t);
  {
    Densa l = r0;
    if (!chol_densa(l)) return std::nan("");
  }

  // observacoes empilhadas: so os (registro, caracteristica) OBSERVADOS
  std::vector<std::pair<std::size_t, std::size_t>> obs;
  for (std::size_t i = 0; i < d.nlin; i++) {
    if (!d.usa[i]) continue;
    for (std::size_t tau = 0; tau < t; tau++)
      if (d.obs.empty() || d.obs[i * t + tau]) obs.push_back({i, tau});
  }
  const std::size_t N = obs.size();

  // V = I (x) R0 + soma_g Z_g (C_g (x) K) Z_g'
  Densa V(N, N);
  for (std::size_t p = 0; p < N; p++)
    for (std::size_t q = 0; q < N; q++)
      if (obs[p].first == obs[q].first)
        V.at(p, q) = r0.at(obs[p].second, obs[q].second);

  for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
    const Grupo& gr = d.modelo.grupos[g];
    const std::size_t dim = gr.dim;
    std::size_t nl = 0;
    for (const DesenhoTermo& a : d.aleatorios)
      if (a.termo == gr.termos[0]) { nl = a.n_niveis; break; }

    Densa cg(dim, dim);
    for (std::size_t j = 0; j < dim; j++)
      for (std::size_t i = j; i < dim; i++) {
        const double v = theta[gr.offset + vech_idx(i, j, dim)];
        cg.at(i, j) = v;
        cg.at(j, i) = v;
      }
    Densa K;
    if (d.kinv[g].ncol > 0) K = inv_pd(d.kinv[g].densa_simetrica());
    else { K = Densa(nl, nl); for (std::size_t i = 0; i < nl; i++) K.at(i, i) = 1.0; }

    // Z empilhada do grupo: N x (dim * nl)
    const std::size_t larg = dim * nl;
    Densa Z(N, larg);
    std::size_t slot = 0;
    for (std::size_t tm : gr.termos)
      for (const DesenhoTermo& a : d.aleatorios)
        if (a.termo == tm) {
          Densa zd = a.z.densa();
          for (std::size_t p = 0; p < N; p++) {
            const std::size_t r = obs[p].first, tau = obs[p].second;
            for (std::size_t c = 0; c < a.z.ncol; c++) {
              const double v = zd.at(r, c);
              if (v == 0.0) continue;
              const std::size_t ct = c / a.n_niveis, nv = c % a.n_niveis;
              // a coordenada do caracter tau no bloco do grupo; sem ela o termo nao existe
              // naquele caracter e a celula fica zero
              std::size_t pt = 0;
              while (pt < a.nt(t) && a.traco(pt) != tau) pt++;
              if (pt == a.nt(t)) continue;
              Z.at(p, (slot + pt * a.n_coef + ct) * nl + nv) = v;
            }
          }
          slot += a.n_coef * a.nt(t);
        }

    // V += Z (C (x) K) Z'
    Densa CK(larg, larg);
    for (std::size_t a = 0; a < dim; a++)
      for (std::size_t b = 0; b < dim; b++)
        for (std::size_t i = 0; i < nl; i++)
          for (std::size_t j = 0; j < nl; j++)
            CK.at(a * nl + i, b * nl + j) = cg.at(a, b) * K.at(i, j);
    Densa ZC(N, larg);
    for (std::size_t r = 0; r < N; r++)
      for (std::size_t c = 0; c < larg; c++) {
        double s = 0.0;
        for (std::size_t k = 0; k < larg; k++) s += Z.at(r, k) * CK.at(k, c);
        ZC.at(r, c) = s;
      }
    for (std::size_t r = 0; r < N; r++)
      for (std::size_t c = 0; c < N; c++) {
        double s = 0.0;
        for (std::size_t k = 0; k < larg; k++) s += ZC.at(r, k) * Z.at(c, k);
        V.at(r, c) += s;
      }
  }

  const double logdet_V = logdet_pd(V);
  if (std::isnan(logdet_V)) return std::nan("");
  Densa Vinv = inv_pd(V);

  const std::size_t p = d.n_fixa;
  Densa X(N, p);
  std::vector<double> y(N);
  for (std::size_t q = 0; q < N; q++) {
    const std::size_t r = obs[q].first, tau = obs[q].second;
    y[q] = d.y.at(r, tau);
    for (std::size_t ex = d.x.linptr[r]; ex < d.x.linptr[r + 1]; ex++) {
      const int e = d.eq_fixa[d.x.coluna[ex] * t + tau];
      if (e >= 0) X.at(q, e) = d.x.valor[ex];
    }
  }

  Densa XtVi(p, N);
  for (std::size_t a = 0; a < p; a++)
    for (std::size_t r = 0; r < N; r++) {
      double s = 0.0;
      for (std::size_t k = 0; k < N; k++) s += X.at(k, a) * Vinv.at(k, r);
      XtVi.at(a, r) = s;
    }
  Densa XtViX(p, p);
  std::vector<double> XtViy(p, 0.0);
  for (std::size_t a = 0; a < p; a++) {
    for (std::size_t b = 0; b < p; b++) {
      double s = 0.0;
      for (std::size_t r = 0; r < N; r++) s += XtVi.at(a, r) * X.at(r, b);
      XtViX.at(a, b) = s;
    }
    for (std::size_t r = 0; r < N; r++) XtViy[a] += XtVi.at(a, r) * y[r];
  }
  const double logdet_X = logdet_pd(XtViX);
  if (std::isnan(logdet_X)) return std::nan("");
  Densa XtViXinv = inv_pd(XtViX);

  double yViy = 0.0;
  for (std::size_t r = 0; r < N; r++) {
    double s = 0.0;
    for (std::size_t k = 0; k < N; k++) s += Vinv.at(r, k) * y[k];
    yViy += y[r] * s;
  }
  double quad = 0.0;
  for (std::size_t a = 0; a < p; a++)
    for (std::size_t b = 0; b < p; b++) quad += XtViy[a] * XtViXinv.at(a, b) * XtViy[b];

  return logdet_V + logdet_X + (yViy - quad) + d.logdet_escala_x;
}

// Residuos no theta final: as linhas de W na mesma convencao da montagem (X na coluna
// eq_fixa[j*t + tau], Z caracteristica-major no coeficiente), o residuo condicional por
// celula observada e Var(e^_{r,tau}) = R0[tau,tau] - w_{r,tau}' C^-1 w_{r,tau}, com os pares
// de uma linha lidos da inversa seletiva (estao no padrao: a montagem empurra todos eles).
ResiduosMT residuos_mt(const DesenhoMT& d, const std::vector<double>& theta,
                       const CacheSimbolica* cache) {
  ResiduosMT S;
  const std::size_t t = d.t;
  S.t = t;
  S.nlin = d.nlin;
  if (theta.size() != d.modelo.ntheta) return S;
  MontadoMT M = monta_mme_mt(d, theta);
  if (!M.ok) return S;
  // A ORDEM E O SIMBOLICO DO AJUSTE quando o padrao e o mesmo (a assinatura do ultimo padrao que o
  // caminho fatorou), como no residuos_ar1: refazer o grau minimo aqui era ~70 s de um ajuste de
  // 417 s num model_mt de 100 mil (E0, 2026-11-02). grau_minimo() so le o padrao, entao com o mesmo
  // padrao a ordem e a mesma e o resultado e o mesmo bit a bit; com outro padrao, refaz como antes.
  const bool reusa = cache && cache->pronto && cache->assinatura == assinatura_padrao(M.c);
  std::vector<std::size_t> perm_local;
  Simbolica sb_local;
  if (!reusa) perm_local = grau_minimo(M.c);
  const std::vector<std::size_t>& perm = reusa ? cache->perm : perm_local;
  const Csc pc = permuta_sim(M.c, perm);
  if (!reusa) sb_local = simbolica(pc);
  const Simbolica& sb = reusa ? cache->sb : sb_local;
  Csc L;
  if (!cholesky(pc, sb, L)) return S;
  std::vector<double> pb(M.total);
  for (std::size_t k = 0; k < M.total; k++) pb[k] = M.rhs[perm[k]];
  const std::vector<double> px = resolve(L, pb);
  std::vector<double> s(M.total);
  for (std::size_t k = 0; k < M.total; k++) s[perm[k]] = px[k];
  SelInv z = inversa_seletiva(L, 0);
  std::vector<std::size_t> inv_perm(M.total);
  for (std::size_t k = 0; k < M.total; k++) inv_perm[perm[k]] = k;

  std::vector<std::vector<std::pair<std::size_t, std::size_t>>> slots(d.modelo.grupos.size());
  for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
    std::size_t col = M.offset_grupo[g];
    for (std::size_t tm : d.modelo.grupos[g].termos)
      for (std::size_t a = 0; a < d.aleatorios.size(); a++)
        if (d.aleatorios[a].termo == tm) {
          slots[g].push_back({a, col});
          col += d.aleatorios[a].z.ncol * d.aleatorios[a].nt(t);
        }
  }
  std::vector<std::vector<std::vector<std::pair<std::uint32_t, double>>>> zl;
  for (const DesenhoTermo& a : d.aleatorios) {
    std::vector<std::vector<std::pair<std::uint32_t, double>>> lin(a.z.nlin);
    for (std::size_t c = 0; c < a.z.ncol; c++)
      for (std::size_t k = a.z.colptr[c]; k < a.z.colptr[c + 1]; k++)
        lin[a.z.linha[k]].push_back({static_cast<std::uint32_t>(c), a.z.valor[k]});
    zl.push_back(std::move(lin));
  }

  const double nan = std::numeric_limits<double>::quiet_NaN();
  S.residuo.assign(d.nlin * t, nan);
  S.ajustado.assign(d.nlin * t, nan);
  S.ajustado_fixo.assign(d.nlin * t, nan);
  S.variancia.assign(d.nlin * t, nan);
  std::vector<EntradaLinha> lin;
  for (std::size_t r = 0; r < d.nlin; r++) {
    if (!d.usa[r]) continue;
    lin.clear();
    for (std::size_t ex = d.x.linptr[r]; ex < d.x.linptr[r + 1]; ex++) {
      const std::size_t j = d.x.coluna[ex];
      const double v = d.x.valor[ex];
      if (v == 0.0) continue;
      for (std::size_t tau = 0; tau < t; tau++) {
        const int e = d.eq_fixa[j * t + tau];
        if (e < 0) continue;
        lin.push_back({static_cast<std::uint32_t>(e), static_cast<std::uint32_t>(tau), v});
      }
    }
    const std::size_t n_x = lin.size();
    for (std::size_t g = 0; g < slots.size(); g++)
      for (const auto& [a, col0] : slots[g]) {
        const DesenhoTermo& dt = d.aleatorios[a];
        for (const auto& [c, v] : zl[a][r]) {
          const std::size_t ct = c / dt.n_niveis;
          const std::size_t nv = c % dt.n_niveis;
          for (std::size_t pt = 0; pt < dt.nt(t); pt++)
            lin.push_back({static_cast<std::uint32_t>(
                               col0 + (pt * dt.n_coef + ct) * dt.n_niveis + nv),
                           static_cast<std::uint32_t>(dt.traco(pt)), v});
        }
      }
    for (std::size_t tau = 0; tau < t; tau++) {
      if (!d.obs.empty() && !d.obs[r * t + tau]) continue;
      double fa = 0.0, fx = 0.0, q = 0.0;
      for (std::size_t p = 0; p < lin.size(); p++) {
        if (lin[p].trait != tau) continue;
        const double v = lin[p].val * s[lin[p].col];
        fa += v;
        if (p < n_x) fx += v;
        for (std::size_t p2 = p; p2 < lin.size(); p2++) {
          if (lin[p2].trait != tau) continue;
          double c;
          if (!z.get(inv_perm[lin[p].col], inv_perm[lin[p2].col], c)) { S.fora++; continue; }
          q += (p == p2 ? 1.0 : 2.0) * lin[p].val * lin[p2].val * c;
        }
      }
      S.ajustado[r * t + tau] = fa;
      S.ajustado_fixo[r * t + tau] = fx;
      S.residuo[r * t + tau] = d.y.at(r, tau) - fa;
      S.variancia[r * t + tau] = M.r0.at(tau, tau) - q;
    }
  }
  S.ok = true;
  return S;
}

}  // namespace br
