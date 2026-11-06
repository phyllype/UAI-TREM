// Inversa seletiva de Takahashi (Takahashi, Fagan e Chin, 1973): os elementos de C^-1 no padrao do fator.
//
// E o que a PEV, a acuracia e os tracos do score precisam, e e a unica parte de C^-1 que se
// pode pagar. A inversa cheia da matriz de coeficientes e densa: a 20.400 colunas sao 3,3 GB,
// e cada posicao fora do padrao e uma que ninguem le.
//
// ## Por que o padrao de L basta, e fecha
//
// Com C = L D L' e L unitaria inferior, a recorrencia de Takahashi calcula, de baixo para
// cima,
//
//   Z[i,j] = - soma_{k > j, L[k,j] != 0}  Z[i,k] . L[k,j]      para i > j
//   Z[j,j] = 1/D[j] - soma_{k > j, L[k,j] != 0}  L[k,j] . Z[k,j]
//
// Todo Z que ela le do lado direito esta numa posicao (i,k) com i e k ambos no padrao da
// coluna j. A propriedade de preenchimento do fator de Cholesky diz que esse par tambem
// esta no padrao. Logo a recorrencia nunca pede um elemento que nao calculou, e nunca
// precisa de um fora de L. Esse fechamento e a razao de tudo isto funcionar.
//
// Aqui o fator e L L', entao D[j] = L[j,j]^2 e a unitaria e L~[k,j] = L[k,j] / L[j,j].

#include "mme.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include <chrono>
#include <cstdlib>
#include <R_ext/Print.h>

namespace br {

// Quantas colunas finais do fator sao completamente densas.
//
// Num modelo de passo unico o bloco genomico e um clique denso, e o produto de Kronecker com
// a covariancia do grupo o propaga: as ultimas colunas de L ficam cheias. Detectar esse bloco
// e o que permite a forma fechada, e ele so existe se a ordenacao deixou o clique no FIM,
// que e por que o pacote ordena por grau minimo e nao por Cuthill-McKee reverso.
std::size_t bloco_denso_simbolico(const Simbolica& sb) {
  const std::size_t n = sb.n;
  std::size_t t = n;
  while (t > 0) {
    const std::size_t j = t - 1;
    if (sb.colptr[j + 1] - sb.colptr[j] != n - j) break;
    t--;
  }
  return n - t;
}

std::size_t bloco_denso_final(const Csc& L) {
  const std::size_t n = L.ncol;
  std::size_t t = n;
  while (t > 0) {
    const std::size_t j = t - 1;
    if (L.colptr[j + 1] - L.colptr[j] != n - j) break;
    t--;
  }
  return n - t;
}

// Metodos de SelInv (declarada em mme.h), fora da classe para nao virarem inline de uma
// definicao que os outros modulos nao veem.
bool SelInv::get(std::size_t i, std::size_t j, double& out) const {
  std::size_t r = i, c = j;
  if (r < c) std::swap(r, c);
  if (c >= n) return false;
  const std::size_t lo = colptr[c], hi = colptr[c + 1];
  auto it = std::lower_bound(linha.begin() + lo, linha.begin() + hi,
                             static_cast<std::uint32_t>(r));
  if (it != linha.begin() + hi && *it == static_cast<std::uint32_t>(r)) {
    out = valor[static_cast<std::size_t>(it - linha.begin())];
    return true;
  }
  return false;
}

std::vector<double> SelInv::diagonal() const {
  std::vector<double> d(n);
  for (std::size_t j = 0; j < n; j++) d[j] = valor[colptr[j]];
  return d;
}

// `bloco` de 0 significa detectar; passe um valor explicito so para forcar a recorrencia pura.
// UAITREM_T3_TEMPOS=1 no ambiente: o tempo de cada fase no stderr (diagnostico do T3)
static bool t3_tempos() { const char* e = std::getenv("UAITREM_T3_TEMPOS"); return e && e[0] == '1'; }
static double t3_agora() {
  return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
SelInv inversa_seletiva(const Csc& L, std::size_t bloco = 0) {
  const bool diag_t = t3_tempos();
  const double tt0 = diag_t ? t3_agora() : 0.0;
  const std::size_t n = L.ncol;
  SelInv z;
  z.n = n;
  z.colptr = L.colptr;
  z.linha = L.linha;
  z.valor.assign(L.nnz(), 0.0);

  const std::size_t t = bloco == 0 ? bloco_denso_final(L) : std::min(bloco, n);
  const std::size_t inicio = n - std::min(t, n);

  // Forma fechada no bloco denso final.
  //
  // Para um bloco final denso T o complemento de Schur e exatamente L[T,T] L[T,T]', logo
  // Z[T,T] = inv(L[T,T] L[T,T]'). Como L[T,T] JA E o fator desse produto, formar o produto
  // e fatorar de novo seria fazer a mesma fatoracao duas vezes. O bloco vai ao nucleo
  // compactado em ladrilhos (inversa_empacotada), em paralelo: copiado das entradas de L
  // (zero onde o padrao nao tem entrada, quando `bloco` pede mais que a cauda densa) e
  // devolvido so nas posicoes do padrao.
  const int nth = threads();
  if (t >= 2) {
    std::vector<std::size_t> off(t + 1, 0);
    for (std::size_t j = 0; j < t; j++) off[j + 1] = off[j] + (t - j);
    std::vector<double> lp(off[t], 0.0), zp(off[t]);
    for (std::size_t j = inicio; j < n; j++)
      for (std::size_t p = L.colptr[j]; p < L.colptr[j + 1]; p++) {
        const std::size_t i = L.linha[p];
        if (i >= inicio) lp[off[j - inicio] + (i - j)] = L.valor[p];
      }
    inversa_empacotada(lp.data(), t, zp.data(), nth);
    for (std::size_t j = inicio; j < n; j++)
      for (std::size_t p = L.colptr[j]; p < L.colptr[j + 1]; p++) {
        const std::size_t i = L.linha[p];
        if (i >= inicio) z.valor[p] = zp[off[j - inicio] + (i - j)];
      }
  }

  const double tt1 = diag_t ? t3_agora() : 0.0;
  // Recorrencia para cima nas colunas restantes.
  //
  // A leitura ingenua da identidade pede, para cada par ORDENADO (i,k) de linhas abaixo da
  // diagonal, o elemento Z[max,min], e acha-lo por busca binaria custa m^2 log sondagens
  // aleatorias numa coluna com m elementos. Nao precisa: a coluna c guarda Z[r,c] para r >= c
  // em ordem crescente, e pela simetria esse UNICO elemento serve aos dois pares ordenados.
  // Percorrer as colunas guardadas com um vetor de selo dizendo quais linhas estao no padrao
  // faz a mesma soma dupla sem busca nenhuma e com leitura sequencial.
  //
  // Um par fora do padrao continua contribuindo zero, que e a mesma semantica que a busca
  // binaria tinha ao nao encontrar.
  //
  // Em paralelo POR NIVEL da arvore de eliminacao: a coluna j le so as colunas das linhas do
  // seu padrao, que sao ancestrais dela na arvore, e escreve so a propria coluna. Colunas do
  // mesmo nivel (mesma distancia ate a cauda ou a raiz) sao independentes; cada uma e feita
  // inteira por uma thread, com a mesma ordem de soma do laco serial, entao o resultado e bit
  // a bit o mesmo com qualquer numero de threads.
  //
  // E DENTRO da coluna, quando ela e larga. Os niveis estreitos ficam perto da cauda, onde as
  // colunas tem centenas ou milhares de linhas, e sao eles que pesam: num passo unico com
  // 3000 genotipados so 10 de 194 niveis tinham 64 colunas ou mais. A soma dupla de uma
  // coluna com m >= COLUNA_LARGA linhas e cortada em pedacos FIXOS de PEDACO colunas c; cada
  // pedaco acumula no seu proprio vetor e os vetores sao somados na ordem dos pedacos. O
  // corte depende so de m, nunca do numero de threads, entao a soma tem a mesma ordem com 1
  // ou 8 threads e o resultado continua bit a bit o mesmo.
  const std::size_t primeiro = (t >= 2) ? inicio : n;
  // SUBARVORES (T3 do E0): com mais de uma thread, a parte esparsa e dividida em subarvores
  // disjuntas da arvore de eliminacao (particao_arvore, linalg.cpp). O TOPO (os nos divididos) vem
  // primeiro, pelos niveis como sempre; depois cada subarvore e feita inteira por uma thread, de
  // cima para baixo, sem barreira entre niveis. A coluna j so le colunas de ancestrais, que estao
  // no topo (prontas) ou acima dela na mesma subarvore (feitas antes pela mesma thread); cada coluna
  // tem a mesma ordem de soma do laco serial, entao o resultado nao muda com o numero de threads.
  Particao P;
  if (nth > 1 && primeiro > 0) {
    std::vector<std::int64_t> pai_l(primeiro, -1);
    std::vector<double> peso(primeiro);
    for (std::size_t j = 0; j < primeiro; j++) {
      const std::size_t lo = L.colptr[j], hi = L.colptr[j + 1];
      if (hi - lo > 1) pai_l[j] = static_cast<std::int64_t>(L.linha[lo + 1]);
      const double m = static_cast<double>(hi - lo);
      peso[j] = m * m;
    }
    P = particao_arvore(pai_l, peso, primeiro, nth);
  }
  std::vector<char> no_topo;
  if (!P.vazia()) {
    no_topo.assign(primeiro, 0);
    for (std::uint32_t j : P.topo) if (j < primeiro) no_topo[j] = 1;
  }
  auto conta = [&](std::size_t j) { return P.vazia() || no_topo[j]; };
  std::vector<std::uint32_t> nivel(primeiro, 0);
  std::uint32_t maior = 0;
  for (std::size_t j = primeiro; j-- > 0;) {
    if (!conta(j)) continue;
    const std::size_t lo = L.colptr[j], hi = L.colptr[j + 1];
    if (hi - lo > 1) {
      const std::size_t pai = L.linha[lo + 1];
      nivel[j] = pai >= primeiro ? 0 : nivel[pai] + 1;
    }
    maior = std::max(maior, nivel[j]);
  }
  std::vector<std::size_t> ini_nivel(static_cast<std::size_t>(maior) + 2, 0);
  for (std::size_t j = 0; j < primeiro; j++) if (conta(j)) ini_nivel[nivel[j] + 1]++;
  for (std::size_t v = 0; v + 1 < ini_nivel.size(); v++) ini_nivel[v + 1] += ini_nivel[v];
  std::vector<std::size_t> por_nivel(ini_nivel.back());
  {
    std::vector<std::size_t> pos = ini_nivel;
    for (std::size_t j = primeiro; j-- > 0;) if (conta(j)) por_nivel[pos[nivel[j]]++] = j;
  }

  struct Rascunho {
    std::vector<std::uint32_t> marca, ondice, linhas;
    std::vector<double> acc, vals, pedacos;
    std::uint32_t selo = 0;
  };
  constexpr std::size_t COLUNA_LARGA = 256, PEDACO = 32;
  std::vector<Rascunho> rasc(static_cast<std::size_t>(std::max(1, nth)));
  // as contribuicoes das colunas c = linhas[ic], ic em [i0, i1), somadas em a
  // as colunas c >= denso_de do fator sao CHEIAS (linhas c..n-1): Z[rr, c] esta na posicao
  // colptr[c] + (rr - c), sem varrer a coluna (T3, 2026-11-02). Antes, uma coluna j com poucas
  // linhas na cauda varria cada coluna da cauda inteira (milhares de entradas) para achar as suas.
  // As linhas do padrao de j a partir de c sao exatamente as que a varredura encontrava, na mesma
  // ordem crescente: a soma e a mesma, bit a bit.
  const std::size_t denso_de = n - bloco_denso_final(L);
  auto parcial = [&](const Rascunho& r, std::size_t i0, std::size_t i1, double* a) {
    const std::uint32_t ultima = r.linhas.back();
    const std::size_t m = r.linhas.size();
    for (std::size_t ic = i0; ic < i1; ic++) {
      const std::size_t c = r.linhas[ic];
      const double lc = r.vals[ic];
      if (c >= denso_de) {
        const double* zc = z.valor.data() + z.colptr[c];
        a[ic] += zc[0] * lc;
        for (std::size_t ir = ic + 1; ir < m; ir++) {
          const double zv = zc[r.linhas[ir] - c];
          a[ir] += zv * lc;
          a[ic] += zv * r.vals[ir];
        }
        continue;
      }
      for (std::size_t p = z.colptr[c]; p < z.colptr[c + 1]; p++) {
        const std::uint32_t rr = z.linha[p];
        if (rr > ultima) break;                // a coluna esta ordenada
        if (r.marca[rr] != r.selo) continue;
        const std::size_t ir = r.ondice[rr];
        const double zv = z.valor[p];
        if (ir == ic) a[ic] += zv * lc;
        else { a[ir] += zv * lc; a[ic] += zv * r.vals[ir]; }
      }
    }
  };
  // nt_dentro > 1 so quando a coluna e chamada fora de uma regiao paralela
  auto coluna = [&](std::size_t j, Rascunho& r, int nt_dentro) {
    if (r.marca.empty()) { r.marca.assign(n, 0); r.ondice.assign(n, 0); r.acc.assign(n, 0.0); }
    const std::size_t lo = L.colptr[j], hi = L.colptr[j + 1];
    const double djj = L.valor[lo];
    const double d = djj * djj;
    r.linhas.clear();
    r.vals.clear();
    for (std::size_t p = lo + 1; p < hi; p++) {
      r.linhas.push_back(L.linha[p]);
      r.vals.push_back(L.valor[p] / djj);
    }
    const std::size_t m = r.linhas.size();
    if (m == 0) { z.valor[lo] = 1.0 / d; return; }
    r.selo++;
    for (std::size_t idx = 0; idx < m; idx++) {
      r.marca[r.linhas[idx]] = r.selo;
      r.ondice[r.linhas[idx]] = static_cast<std::uint32_t>(idx);
      r.acc[idx] = 0.0;
    }
    if (m < COLUNA_LARGA) {
      parcial(r, 0, m, r.acc.data());
    } else {
      const std::size_t np = (m + PEDACO - 1) / PEDACO;
      if (r.pedacos.size() < np * m) r.pedacos.resize(np * m);
      double* pd = r.pedacos.data();
#ifdef _OPENMP
#pragma omp parallel for num_threads(nt_dentro) if (nt_dentro > 1) schedule(dynamic, 1)
#endif
      for (long k = 0; k < static_cast<long>(np); k++) {
        const std::size_t kk = static_cast<std::size_t>(k);
        double* a = pd + kk * m;
        std::fill(a, a + m, 0.0);
        parcial(r, kk * PEDACO, std::min(m, (kk + 1) * PEDACO), a);
      }
#ifdef _OPENMP
#pragma omp parallel for num_threads(nt_dentro) if (nt_dentro > 1) schedule(static)
#endif
      for (long ii = 0; ii < static_cast<long>(m); ii++) {
        const std::size_t idx = static_cast<std::size_t>(ii);
        double t = 0.0;
        for (std::size_t k = 0; k < np; k++) t += pd[k * m + idx];
        r.acc[idx] = t;
      }
    }
    for (std::size_t idx = 0; idx < m; idx++) z.valor[lo + 1 + idx] = -r.acc[idx];
    double dg = 0.0;
    for (std::size_t idx = 0; idx < m; idx++) dg += r.vals[idx] * z.valor[lo + 1 + idx];
    z.valor[lo] = 1.0 / d - dg;
  };
  for (std::size_t v = 0; v + 1 < ini_nivel.size(); v++) {
    const std::size_t a = ini_nivel[v], b = ini_nivel[v + 1];
    if (b - a < 64 || nth <= 1) {
      for (std::size_t q = a; q < b; q++) coluna(por_nivel[q], rasc[0], nth);
      continue;
    }
#ifdef _OPENMP
#pragma omp parallel num_threads(nth)
#endif
    {
#ifdef _OPENMP
      Rascunho& r = rasc[static_cast<std::size_t>(omp_get_thread_num())];
#pragma omp for schedule(dynamic, 16)
#else
      Rascunho& r = rasc[0];
#endif
      for (long q = static_cast<long>(a); q < static_cast<long>(b); q++)
        coluna(por_nivel[static_cast<std::size_t>(q)], r, 1);
    }
  }

  const double tt2 = diag_t ? t3_agora() : 0.0;
  if (!P.vazia()) {
    // o rascunho de cada thread alocado fora da regiao paralela
    for (auto& r : rasc)
      if (r.marca.empty()) { r.marca.assign(n, 0); r.ondice.assign(n, 0); r.acc.assign(n, 0.0); }
    const std::size_t nsub = P.ini.size() - 1;
#ifdef _OPENMP
#pragma omp parallel num_threads(nth)
#endif
    {
#ifdef _OPENMP
      Rascunho& r = rasc[static_cast<std::size_t>(omp_get_thread_num())];
#pragma omp for schedule(dynamic, 1)
#else
      Rascunho& r = rasc[0];
#endif
      for (long q = 0; q < static_cast<long>(nsub); q++) {
        const std::size_t sq = static_cast<std::size_t>(q);
        for (std::size_t u = P.ini[sq + 1]; u-- > P.ini[sq];) coluna(P.nos[u], r, 1);
      }
    }
  }

  if (diag_t) {
    double wt = 0.0, ws = 0.0;
    for (std::size_t j = 0; j < primeiro; j++) {
      const double m = static_cast<double>(L.colptr[j + 1] - L.colptr[j]);
      if (conta(j)) wt += m * m; else ws += m * m;
    }
    REprintf("[t3 selinv] n %zu cauda %zu: cauda %.3f s, topo %.3f s (%zu colunas, %.1f%% do peso), "
             "subarvores %.3f s (%zu), threads %d\n", n, t, tt1 - tt0, tt2 - tt1,
             P.vazia() ? primeiro : P.topo.size(), 100.0 * wt / std::max(1.0, wt + ws), t3_agora() - tt2,
             P.vazia() ? std::size_t(0) : P.ini.size() - 1, nth);
  }
  return z;
}

}  // namespace br
