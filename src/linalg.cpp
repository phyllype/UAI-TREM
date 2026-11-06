// Cholesky esparsa: arvore de eliminacao, fatoracao simbolica, fatoracao numerica
// olhando-para-cima, resolucao triangular, permutacao simetrica e grau minimo.
//
// E o que tira o pacote do caminho denso. A Cholesky densa de estruturas.cpp fica como
// REFERENCIA: tudo aqui e conferido contra ela, nunca contra si mesmo.
//
// Duas ideias decidem o formato do codigo (Davis, 2006, Direct Methods for Sparse Linear Systems):
//
//  - a ARVORE DE ELIMINACAO diz, para cada coluna, a que coluna anterior o primeiro
//    elemento fora da diagonal pertence. Subir a arvore da o padrao de nao-zeros de uma
//    LINHA de L sem formar nada.
//  - a fatoracao OLHANDO-PARA-CIMA calcula L uma linha de cada vez, resolvendo um sistema
//    triangular contra a parte ja pronta. Ela precisa do padrao da linha primeiro, que e
//    exatamente o que a arvore da, e custa O(nnz(L)) em vez de algo quadratico.
//
// A entrada e o triangulo SUPERIOR de uma simetrica em CSC. Para uma simetrica isso e a
// transposta do triangulo inferior.

#include "mme.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include <R_ext/Utils.h>
#define R_NO_REMAP
#include <Rinternals.h>
#include <chrono>
#include <cstdlib>
#include <R_ext/Print.h>

namespace br {

// INTERRUPCAO do usuario (Ctrl-C, ou um setTimeLimit()) entre passos sequenciais longos. O
// R_CheckUserInterrupt roda dentro de R_ToplevelExec, e a interrupcao volta como excecao
// C++: o longjmp direto pularia os destrutores dos vetores vivos (o fator, o fecho de
// pares), que ficariam alocados ate o fim do processo. Dentro de regiao paralela nao faz
// nada, porque o R so pode ser chamado da thread principal.
namespace {
void checa_r(void*) { R_CheckUserInterrupt(); }
}  // namespace
void checa_interrupcao() {
#ifdef _OPENMP
  if (omp_in_parallel()) return;
#endif
  if (!R_ToplevelExec(checa_r, nullptr)) throw Interrompido();
}

// Triangulo superior, incluindo a diagonal.
Csc triu(const Csc& a) {
  Csc out(a.nlin, a.ncol);
  out.linha.reserve(a.nnz() / 2 + a.ncol);
  out.valor.reserve(a.nnz() / 2 + a.ncol);
  for (std::size_t c = 0; c < a.ncol; c++) {
    for (std::size_t k = a.colptr[c]; k < a.colptr[c + 1]; k++)
      if (a.linha[k] <= c) { out.linha.push_back(a.linha[k]); out.valor.push_back(a.valor[k]); }
    out.colptr[c + 1] = out.linha.size();
  }
  return out;
}

// Arvore de eliminacao, com compressao de caminho.
std::vector<std::int64_t> arvore(const Csc& au) {
  const std::size_t n = au.ncol;
  std::vector<std::int64_t> pai(n, -1), ancestral(n, -1);
  for (std::size_t k = 0; k < n; k++)
    for (std::size_t p = au.colptr[k]; p < au.colptr[k + 1]; p++) {
      std::size_t i = au.linha[p];
      if (i >= k) continue;
      while (true) {
        const std::int64_t prox = ancestral[i];
        ancestral[i] = static_cast<std::int64_t>(k);
        if (prox < 0) { pai[i] = static_cast<std::int64_t>(k); break; }
        if (static_cast<std::size_t>(prox) == k) break;
        i = static_cast<std::size_t>(prox);
      }
    }
  return pai;
}

// Padrao de nao-zeros da linha k de L, em ordem topologica, devolvido em s[topo..n).
static std::size_t alcance(const Csc& au, std::size_t k, const std::vector<std::int64_t>& pai,
                           std::vector<std::uint32_t>& s, std::vector<std::uint32_t>& marca) {
  const std::size_t n = au.ncol;
  std::size_t topo = n;
  marca[k] = static_cast<std::uint32_t>(k) + 1;
  for (std::size_t p = au.colptr[k]; p < au.colptr[k + 1]; p++) {
    std::size_t i = au.linha[p];
    if (i > k) continue;
    std::size_t len = 0;
    while (marca[i] != static_cast<std::uint32_t>(k) + 1) {
      s[len++] = static_cast<std::uint32_t>(i);
      marca[i] = static_cast<std::uint32_t>(k) + 1;
      if (pai[i] < 0) break;
      i = static_cast<std::size_t>(pai[i]);
    }
    while (len > 0) s[--topo] = s[--len];
  }
  return topo;
}

Simbolica simbolica(const Csc& au) {
  const std::size_t n = au.ncol;
  Simbolica sb;
  sb.n = n;
  sb.pai = arvore(au);
  std::vector<std::size_t> conta(n, 0);
  std::vector<std::uint32_t> s(n), marca(n, 0);
  for (std::size_t k = 0; k < n; k++) {
    const std::size_t topo = alcance(au, k, sb.pai, s, marca);
    for (std::size_t t = topo; t < n; t++) conta[s[t]]++;
    conta[k]++;                      // a diagonal
  }
  sb.colptr.assign(n + 1, 0);
  for (std::size_t j = 0; j < n; j++) sb.colptr[j + 1] = sb.colptr[j] + conta[j];
  return sb;
}


// ---------------------------------------------------------------- paralelismo

static int g_threads = 1;
int threads() { return g_threads; }
static bool g_lapack = false;
bool denso_lapack() { return g_lapack; }
void define_denso_lapack(bool v) { g_lapack = v; }
// o pedido e limitado pelo OMP_THREAD_LIMIT do ambiente (o CRAN usa esse limite nos checks):
// um pacote que passa por cima do limite de threads de quem o chama e rejeitado
void define_threads(int n) {
  if (n < 1) n = 1;
#ifdef _OPENMP
  const int lim = omp_get_thread_limit();
  if (lim >= 1 && n > lim) n = lim;
#endif
  g_threads = n;
}
int threads_disponiveis() {
#ifdef _OPENMP
  return omp_get_num_procs();
#else
  return 1;
#endif
}

// PARTICAO DA ARVORE DE ELIMINACAO (T3 do E0, 2026-11-02). Na fatoracao olhando-para-cima a linha
// k le so as colunas j da propria linha, que sao DESCENDENTES de k na arvore, e escreve so nelas
// (L[k,j] na posicao prox[j]). Duas subarvores disjuntas nao compartilham coluna nenhuma: as linhas
// de uma podem ser feitas por uma thread enquanto outra faz as da outra, cada linha com a MESMA
// ordem de soma do laco serial, e o fator sai bit a bit o mesmo com qualquer numero de threads. A
// ordem das entradas de cada coluna tambem nao muda: quem escreve na coluna j sao os ancestrais de
// j, os da subarvore dela (feitos em ordem crescente pela mesma thread) e os do topo, que tem
// indice maior que a raiz da subarvore e vem depois de todas elas. Na inversa seletiva e o
// contrario: a coluna j le as colunas dos ancestrais, entao o topo vem ANTES e as subarvores
// depois, cada uma de cima para baixo.
//
// A divisao e a gulosa de Geist e Ng (1989): parte das raizes e divide a subarvore mais pesada
// (o no vai para o topo, os filhos entram no lugar) ate a mais pesada caber em 1/(8 nth) do total.
// O peso de um no e o custo da coluna dele (c_j^2). Uma cadeia longa no alto da arvore vai inteira
// para o topo: e serial aqui e continua sendo (a cauda densa e a largura das colunas cuidam dela).
Particao particao_arvore(const std::vector<std::int64_t>& pai, const std::vector<double>& peso,
                         std::size_t limite, int nth) {
  Particao P;
  const std::size_t n = pai.size();
  limite = std::min(limite, n);
  if (nth <= 1 || limite < 1024) return P;
  std::vector<double> sw(n, 0.0);
  for (std::size_t j = 0; j < limite; j++) sw[j] += peso[j];
  for (std::size_t j = 0; j < limite; j++)
    if (pai[j] >= 0 && static_cast<std::size_t>(pai[j]) < limite) sw[static_cast<std::size_t>(pai[j])] += sw[j];
  std::vector<std::size_t> cini(limite + 1, 0);
  for (std::size_t j = 0; j < limite; j++)
    if (pai[j] >= 0 && static_cast<std::size_t>(pai[j]) < limite) cini[static_cast<std::size_t>(pai[j]) + 1]++;
  for (std::size_t j = 0; j < limite; j++) cini[j + 1] += cini[j];
  std::vector<std::uint32_t> fil(cini[limite]);
  {
    std::vector<std::size_t> pos(cini.begin(), cini.end() - 1);
    for (std::size_t j = 0; j < limite; j++)
      if (pai[j] >= 0 && static_cast<std::size_t>(pai[j]) < limite)
        fil[pos[static_cast<std::size_t>(pai[j])]++] = static_cast<std::uint32_t>(j);
  }
  double total = 0.0;
  std::vector<std::pair<double, std::uint32_t>> heap;
  for (std::size_t j = 0; j < limite; j++)
    if (pai[j] < 0 || static_cast<std::size_t>(pai[j]) >= limite) {
      heap.emplace_back(sw[j], static_cast<std::uint32_t>(j));
      total += sw[j];
    }
  // pouco trabalho na parte esparsa: o serial ganha (abrir a regiao paralela custa mais)
  if (total < 2e7) return P;
  std::make_heap(heap.begin(), heap.end());
  const double alvo = total / (8.0 * static_cast<double>(nth));
  std::vector<char> topo(n, 0);
  for (std::size_t j = limite; j < n; j++) topo[j] = 1;
  while (!heap.empty() && heap.front().first > alvo) {
    std::pop_heap(heap.begin(), heap.end());
    const std::uint32_t r = heap.back().second;
    heap.pop_back();
    topo[r] = 1;
    for (std::size_t q = cini[r]; q < cini[r + 1]; q++) {
      heap.emplace_back(sw[fil[q]], fil[q]);
      std::push_heap(heap.begin(), heap.end());
    }
  }
  if (heap.size() < 2) return P;
  std::sort(heap.begin(), heap.end(), [](const std::pair<double, std::uint32_t>& x,
                                         const std::pair<double, std::uint32_t>& y) {
    return x.first != y.first ? x.first > y.first : x.second < y.second;
  });
  std::vector<std::int32_t> rot(n, -1);
  for (std::size_t s2 = 0; s2 < heap.size(); s2++) rot[heap[s2].second] = static_cast<std::int32_t>(s2);
  for (std::size_t j = limite; j-- > 0;) {
    if (topo[j]) continue;
    if (rot[j] < 0) rot[j] = rot[static_cast<std::size_t>(pai[j])];
  }
  P.ini.assign(heap.size() + 1, 0);
  for (std::size_t j = 0; j < limite; j++) if (!topo[j]) P.ini[static_cast<std::size_t>(rot[j]) + 1]++;
  for (std::size_t s2 = 0; s2 < heap.size(); s2++) P.ini[s2 + 1] += P.ini[s2];
  P.nos.resize(P.ini.back());
  {
    std::vector<std::size_t> pos(P.ini.begin(), P.ini.end() - 1);
    for (std::size_t j = 0; j < limite; j++)
      if (!topo[j]) P.nos[pos[static_cast<std::size_t>(rot[j])]++] = static_cast<std::uint32_t>(j);
  }
  for (std::size_t j = 0; j < n; j++) if (topo[j]) P.topo.push_back(static_cast<std::uint32_t>(j));
  return P;
}

// A CAUDA DENSA EM LADRILHOS. O bloco denso final do fator (o dense_block, as ultimas T
// colunas completamente cheias) e onde mora o k^3 de cada fatoracao num passo unico. O laco
// olhando-para-cima trata essa cauda como qualquer coluna esparsa: escalar, com acesso
// indireto, sem reuso de cache. Aqui ela e fatorada por blocos, right-looking, com ladrilhos
// de 64: a coluna j da cauda ja esta guardada em CSC como o triangulo inferior compactado
// POR COLUNA (linhas j..n-1 em sequencia), entao a conta e feita no proprio vetor de valores,
// sem copia T x T, e o laco interno corre contiguo em i nas duas colunas envolvidas.
//
// Cada entrada e escrita por UMA thread em cada passo (dono do ladrilho) e a soma sobre p
// vai em ordem fixa: o fator e bit a bit o mesmo para qualquer numero de threads. Nada de
// reduction(+:), cuja ordem de soma muda com o numero de threads.
namespace {
const std::size_t LADRILHO = 64;
const std::size_t CAUDA_MINIMA = 128;

struct Cauda {
  std::size_t base, T;
  const std::vector<std::size_t>* colptr;
  std::vector<double>* valor;
  double& at(std::size_t i, std::size_t j) {        // local, i >= j
    return (*valor)[(*colptr)[base + j] + (i - j)];
  }
};

// C[i, j] -= soma_p A[i, p] B[j, p], p em [0, kw), para j em [j0, j1) e i em [max(i0, j), i1).
// A e B vem EMPACOTADOS por ladrilho: pa[p * LADRILHO + (i - i0)], um bloco contiguo de
// 64 x kw. Na cauda guardada por coluna, cada passo em p pulava para outra coluna (outra
// pagina de memoria quando a cauda passa de alguns milhares); empacotado, o miolo anda em
// sequencia. O miolo e um bloco 4 x 4 em registradores: 8 leituras e 32 flops por passo.
// A soma em p vai num acumulador e so depois sai de C: ordem fixa, a mesma com qualquer
// numero de threads.
void atualiza(Cauda& c, const double* pa, const double* pb, std::size_t i0, std::size_t i1,
              std::size_t j0, std::size_t j1, std::size_t kw) {
  const std::size_t LD = LADRILHO;
  auto escalar = [&](std::size_t i, std::size_t j) {
    double acc = 0.0;
    for (std::size_t p = 0; p < kw; p++) acc += pa[p * LD + (i - i0)] * pb[p * LD + (j - j0)];
    c.at(i, j) -= acc;
  };
  std::size_t j = j0;
  for (; j + 4 <= j1; j += 4) {
    std::size_t i = std::max(i0, j);
    // o triangulo sobre a diagonal (so num ladrilho da diagonal): so i >= j + s existe
    for (; i < std::min(i1, j + 4); i++)
      for (std::size_t s2 = 0; s2 < 4 && j + s2 <= i; s2++) escalar(i, j + s2);
    for (; i + 4 <= i1; i += 4) {
      double acc[4][4] = {{0.0}};
      const double* a = pa + (i - i0);
      const double* b = pb + (j - j0);
      for (std::size_t p = 0; p < kw; p++, a += LD, b += LD)
        for (std::size_t r = 0; r < 4; r++)
          for (std::size_t s2 = 0; s2 < 4; s2++) acc[r][s2] += a[r] * b[s2];
      for (std::size_t s2 = 0; s2 < 4; s2++) {
        double* cj = &c.at(i, j + s2);
        for (std::size_t r = 0; r < 4; r++) cj[r] -= acc[r][s2];
      }
    }
    for (; i < i1; i++)
      for (std::size_t s2 = 0; s2 < 4; s2++) escalar(i, j + s2);
  }
  for (; j < j1; j++)
    for (std::size_t i = std::max(i0, j); i < i1; i++) escalar(i, j);
}

bool fatora_cauda(Cauda& c, int nth) {
  const std::size_t T = c.T;
  const std::size_t nb = (T + LADRILHO - 1) / LADRILHO;
  // o painel da vez, empacotado por ladrilho de linhas (preenchido no passo 2)
  std::vector<double> painel(nb * LADRILHO * LADRILHO, 0.0);
  for (std::size_t kb = 0; kb < nb; kb++) {
    // entre um painel e o proximo (fora das regioes paralelas): uma cauda de 15 mil colunas
    // fatora em minutos
    checa_interrupcao();
    const std::size_t k0 = kb * LADRILHO, k1 = std::min(T, k0 + LADRILHO);
    const std::size_t kw = k1 - k0;
    // 1. ladrilho da diagonal, serial (left-looking dentro dele)
    for (std::size_t j = k0; j < k1; j++) {
      double d = c.at(j, j);
      for (std::size_t p = k0; p < j; p++) d -= c.at(j, p) * c.at(j, p);
      if (!(d > 0.0) || !std::isfinite(d)) return false;
      d = std::sqrt(d);
      c.at(j, j) = d;
      for (std::size_t i = j + 1; i < k1; i++) {
        double s = c.at(i, j);
        for (std::size_t p = k0; p < j; p++) s -= c.at(i, p) * c.at(j, p);
        c.at(i, j) = s / d;
      }
    }
    if (k1 >= T) break;
    // 2. painel: as linhas abaixo, um ladrilho de linhas por thread, que tambem o empacota
    const std::size_t nlt = (T - k1 + LADRILHO - 1) / LADRILHO;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nth)
#endif
    for (long it = 0; it < static_cast<long>(nlt); it++) {
      const std::size_t i0 = k1 + static_cast<std::size_t>(it) * LADRILHO;
      const std::size_t i1 = std::min(T, i0 + LADRILHO);
      for (std::size_t j = k0; j < k1; j++) {
        for (std::size_t p = k0; p < j; p++) {
          const double ljp = c.at(j, p);
          double* col_j = &c.at(i0, j);
          const double* col_p = &c.at(i0, p);
          for (std::size_t i = 0; i < i1 - i0; i++) col_j[i] -= col_p[i] * ljp;
        }
        const double djj = c.at(j, j);
        double* col_j = &c.at(i0, j);
        for (std::size_t i = 0; i < i1 - i0; i++) col_j[i] /= djj;
      }
      double* pk = &painel[static_cast<std::size_t>(it) * LADRILHO * LADRILHO];
      for (std::size_t p = 0; p < kw; p++) {
        const double* col = &c.at(i0, k0 + p);
        for (std::size_t i = 0; i < i1 - i0; i++) pk[p * LADRILHO + i] = col[i];
      }
    }
    // 3. atualizacao do resto: cada ladrilho (ib, jb), jb <= ib, de UMA thread
    const std::size_t t0 = kb + 1;
    const std::size_t nt = nb - t0;
    const long nlad = static_cast<long>(nt * (nt + 1) / 2);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(nth)
#endif
    for (long t = 0; t < nlad; t++) {
      std::size_t ib = static_cast<std::size_t>((std::sqrt(8.0 * static_cast<double>(t) + 1.0) - 1.0) / 2.0);
      while (ib * (ib + 1) / 2 > static_cast<std::size_t>(t)) ib--;
      while ((ib + 1) * (ib + 2) / 2 <= static_cast<std::size_t>(t)) ib++;
      const std::size_t jb = static_cast<std::size_t>(t) - ib * (ib + 1) / 2;
      const std::size_t i0 = (ib + t0) * LADRILHO, i1 = std::min(T, i0 + LADRILHO);
      const std::size_t j0 = (jb + t0) * LADRILHO, j1 = std::min(T, j0 + LADRILHO);
      atualiza(c, &painel[ib * LADRILHO * LADRILHO], &painel[jb * LADRILHO * LADRILHO],
               i0, i1, j0, j1, kw);
    }
  }
  return true;
}
}  // namespace

// O produto X'Y dos kernels densos (a G de VanRaden, os dois produtos da APY), no lugar do
// BLAS de referencia que o R traz no Windows, que roda numa thread. Por bloco de KB valores
// de k: (1) empacota os paineis de 64 linhas de X e de Y, cada um contiguo [k][64], em
// paralelo por painel; (2) cada ladrilho 64 x 64 de C e de UMA thread, com miolo 4 x 4 em
// registradores (8 leituras e 32 flops por passo de k), e so ao fim do bloco o acumulador
// sai para C. A soma de cada entrada anda em k crescente, bloco a bloco: a mesma com
// qualquer numero de threads. x_por_linha: X guardada com k CONTIGUO, o elemento (k, i) em
// x[i ldx + k] (uma matriz por linhas multiplicando pela direita, o Z W de uma Z n x m);
// so o empacotamento muda.
void produto_ladrilhos(const double* x, std::size_t ldx, const double* y, std::size_t ldy,
                       std::size_t K, std::size_t m, std::size_t n, double* c,
                       std::size_t ldc, bool simetrico, int nth, bool x_por_linha) {
  const std::size_t L = LADRILHO, KB = 256;
  const std::size_t mt = (m + L - 1) / L, nt = (n + L - 1) / L;
  if (mt == 0 || nt == 0 || K == 0) return;
  std::vector<double> px(mt * KB * L), py(simetrico ? 0 : nt * KB * L);
  const long nlad = simetrico ? static_cast<long>(mt * (mt + 1) / 2)
                              : static_cast<long>(mt * nt);
  for (std::size_t k0 = 0; k0 < K; k0 += KB) {
    const std::size_t kw = std::min(KB, K - k0);
    auto empacota = [&](const double* src, std::size_t ld, std::size_t lim, std::size_t t,
                        std::vector<double>& dst, bool por_linha) {
      double* d = &dst[t * KB * L];
      const std::size_t i0 = t * L, w = std::min(L, lim - i0);
      if (por_linha) {
        for (std::size_t r = 0; r < L; r++) {
          if (r < w) {
            const double* s = src + (i0 + r) * ld + k0;
            for (std::size_t kk = 0; kk < kw; kk++) d[kk * L + r] = s[kk];
          } else {
            for (std::size_t kk = 0; kk < kw; kk++) d[kk * L + r] = 0.0;
          }
        }
        return;
      }
      for (std::size_t kk = 0; kk < kw; kk++) {
        const double* s = src + (k0 + kk) * ld + i0;
        double* dk = d + kk * L;
        std::size_t r = 0;
        for (; r < w; r++) dk[r] = s[r];
        for (; r < L; r++) dk[r] = 0.0;
      }
    };
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nth)
#endif
    for (long t = 0; t < static_cast<long>(mt + (simetrico ? 0 : nt)); t++) {
      const std::size_t tt = static_cast<std::size_t>(t);
      if (tt < mt) empacota(x, ldx, m, tt, px, x_por_linha);
      else empacota(y, ldy, n, tt - mt, py, false);
    }
    const std::vector<double>& qy = simetrico ? px : py;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(nth)
#endif
    for (long t = 0; t < nlad; t++) {
      std::size_t it, jt;
      if (simetrico) {
        it = static_cast<std::size_t>((std::sqrt(8.0 * static_cast<double>(t) + 1.0) - 1.0) / 2.0);
        while (it * (it + 1) / 2 > static_cast<std::size_t>(t)) it--;
        while ((it + 1) * (it + 2) / 2 <= static_cast<std::size_t>(t)) it++;
        jt = static_cast<std::size_t>(t) - it * (it + 1) / 2;
      } else {
        it = static_cast<std::size_t>(t) / nt;
        jt = static_cast<std::size_t>(t) % nt;
      }
      const double* pa = &px[it * KB * L];
      const double* pb = &qy[jt * KB * L];
      const std::size_t i0 = it * L, j0 = jt * L;
      const std::size_t iw = std::min(L, m - i0), jw = std::min(L, n - j0);
      for (std::size_t r0 = 0; r0 < iw; r0 += 4)
        for (std::size_t s0 = 0; s0 < jw; s0 += 4) {
          if (simetrico && j0 + s0 > i0 + r0 + 3) continue;   // bloco todo acima da diagonal
          double acc[4][4] = {{0.0}};
          const double* a = pa + r0;
          const double* b = pb + s0;
          for (std::size_t kk = 0; kk < kw; kk++, a += L, b += L)
            for (std::size_t r = 0; r < 4; r++)
              for (std::size_t s2 = 0; s2 < 4; s2++) acc[r][s2] += a[r] * b[s2];
          for (std::size_t r = 0; r < 4 && r0 + r < iw; r++)
            for (std::size_t s2 = 0; s2 < 4 && s0 + s2 < jw; s2++) {
              const std::size_t i = i0 + r0 + r, j = j0 + s0 + s2;
              if (simetrico && j > i) continue;
              c[i * ldc + j] += acc[r][s2];
            }
        }
    }
  }
}

bool cholesky_empacotada(std::vector<double>& v, std::size_t n, int nth) {
  std::vector<std::size_t> colptr(n + 1, 0);
  for (std::size_t j = 0; j < n; j++) colptr[j + 1] = colptr[j] + (n - j);
  if (v.size() != colptr[n]) throw Erro("packed triangle of the wrong size");
  for (const double x : v)
    if (!std::isfinite(x)) return false;
  Cauda c{0, n, &colptr, &v};
  return fatora_cauda(c, nth);
}

// Fatoracao numerica olhando-para-cima. Devolve false se a matriz nao for positiva-definida,
// o que, neste engine, e informacao sobre theta e nao uma falha a reportar como tal.
bool cholesky(const Csc& au, const Simbolica& sb, Csc& L) {
  const std::size_t n = sb.n;
  const std::size_t nz = sb.nnz();
  std::vector<std::size_t> colptr = sb.colptr;
  std::vector<std::uint32_t> linha(nz, 0);
  std::vector<double> valor(nz, 0.0);
  std::vector<std::size_t> prox(colptr.begin(), colptr.begin() + n);

  std::vector<double> x(n, 0.0);
  std::vector<std::uint32_t> s(n), marca(n, 0);

  // A cauda densa vai para os ladrilhos quando e grande o bastante para compensar. As linhas
  // dela calculam aqui so a parte FORA da cauda (W = L[cauda, resto]); o bloco da cauda sai
  // depois, como a fatoracao de A_TT - W W'.
  const std::size_t T_cauda = bloco_denso_simbolico(sb);
  const bool hibrido = T_cauda >= CAUDA_MINIMA;
  const std::size_t base = hibrido ? n - T_cauda : n;

  // A LINHA k do fator, com o rascunho (x, s, marca) de quem a faz. 0: feita; 1: a matriz nao e
  // positiva-definida; 2: a simbolica e menor que o padrao (invariante de quem chama).
  auto faz_linha = [&](std::size_t k, std::vector<double>& x, std::vector<std::uint32_t>& s,
                       std::vector<std::uint32_t>& marca) -> int {
    const std::size_t topo = alcance(au, k, sb.pai, s, marca);

    if (k >= base) {
      for (std::size_t p = au.colptr[k]; p < au.colptr[k + 1]; p++) {
        const std::size_t i = au.linha[p];
        if (i < base) x[i] = au.valor[p];
      }
      for (std::size_t t = topo; t < n; t++) {
        const std::size_t j = s[t];
        if (j >= base) continue;
        const double lkj = x[j] / valor[colptr[j]];
        x[j] = 0.0;
        // as linhas da coluna j estao em ordem crescente: as da cauda vem por ultimo
        for (std::size_t p = colptr[j] + 1; p < prox[j] && linha[p] < base; p++)
          x[linha[p]] -= valor[p] * lkj;
        if (prox[j] >= colptr[j + 1]) return 2;
        linha[prox[j]] = static_cast<std::uint32_t>(k);
        valor[prox[j]] = lkj;
        prox[j]++;
      }
      return 0;
    }

    double d = 0.0;
    for (std::size_t p = au.colptr[k]; p < au.colptr[k + 1]; p++) {
      const std::size_t i = au.linha[p];
      if (i < k) x[i] = au.valor[p];
      else if (i == k) d = au.valor[p];
    }

    for (std::size_t t = topo; t < n; t++) {
      const std::size_t j = s[t];
      const double lkj = x[j] / valor[colptr[j]];   // L[j,j] e o primeiro da coluna
      x[j] = 0.0;
      for (std::size_t p = colptr[j] + 1; p < prox[j]; p++) x[linha[p]] -= valor[p] * lkj;
      d -= lkj * lkj;
      // Uma simbolica menor que o padrao real transbordaria a coluna j para a j+1 e
      // corromperia L sem erro nenhum. E invariante de quem chama, entao para aqui.
      if (prox[j] >= colptr[j + 1]) return 2;
      linha[prox[j]] = static_cast<std::uint32_t>(k);
      valor[prox[j]] = lkj;
      prox[j]++;
    }

    if (!(d > 0.0) || !std::isfinite(d)) return 1;
    linha[prox[k]] = static_cast<std::uint32_t>(k);
    valor[prox[k]] = std::sqrt(d);   // a diagonal vai PRIMEIRO, e o laco acima conta com isso
    prox[k]++;
    return 0;
  };
  auto falha = [](int rc) -> bool {
    if (rc == 2) throw Erro("symbolic factorization too small: the matrix pattern grew afterwards");
    return rc == 1;
  };

  // a interrupcao e conferida pelo RELOGIO, olhado a cada 64 colunas: uma coluna custa de
  // nada (A^-1) a muito (enchimento), e contar colunas daria intervalos de segundos a minutos
  using relogio = std::chrono::steady_clock;
  relogio::time_point prox_checa = relogio::now() + std::chrono::milliseconds(250);
  auto checa = [&](std::size_t q) {
    if ((q & 63) == 0 && relogio::now() >= prox_checa) {
      checa_interrupcao();
      prox_checa = relogio::now() + std::chrono::milliseconds(250);
    }
  };

  // as SUBARVORES em paralelo (particao_arvore), depois o TOPO em ordem crescente; com uma thread,
  // ou pouco trabalho, o laco serial de sempre. As duas formas dao o mesmo fator bit a bit.
  const int nth_part = threads();
  Particao P;
  if (nth_part > 1) {
    std::vector<double> peso(n);
    for (std::size_t j = 0; j < n; j++) {
      const double c = static_cast<double>(sb.colptr[j + 1] - sb.colptr[j]);
      peso[j] = c * c;
    }
    P = particao_arvore(sb.pai, peso, base, nth_part);
  }
  const char* t3e = std::getenv("UAITREM_T3_TEMPOS");
  const bool diag_t = t3e && t3e[0] == '1';
  auto agora = [] {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
  };
  const double tt0 = diag_t ? agora() : 0.0;
  double tt1 = tt0;
  if (P.vazia()) {
    for (std::size_t k = 0; k < n; k++) {
      checa(k);
      if (falha(faz_linha(k, x, s, marca))) return false;
    }
  } else {
    const std::size_t nsub = P.ini.size() - 1;
    int rc_geral = 0;
    std::vector<std::vector<double>> xs(static_cast<std::size_t>(nth_part));
    std::vector<std::vector<std::uint32_t>> ss(static_cast<std::size_t>(nth_part)),
        ms(static_cast<std::size_t>(nth_part));
    // o rascunho das outras threads alocado aqui, fora da regiao paralela (um bad_alloc la dentro
    // derrubaria o processo em vez de voltar ao R)
    for (std::size_t t = 1; t < static_cast<std::size_t>(nth_part); t++) {
      xs[t].assign(n, 0.0); ss[t].assign(n, 0); ms[t].assign(n, 0);
    }
#ifdef _OPENMP
#pragma omp parallel num_threads(nth_part)
#endif
    {
#ifdef _OPENMP
      const std::size_t eu = static_cast<std::size_t>(omp_get_thread_num());
#else
      const std::size_t eu = 0;
#endif
      std::vector<double>& xx = eu == 0 ? x : xs[eu];
      std::vector<std::uint32_t>& s1 = eu == 0 ? s : ss[eu];
      std::vector<std::uint32_t>& m1 = eu == 0 ? marca : ms[eu];
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
      for (long q = 0; q < static_cast<long>(nsub); q++) {
        int rc_visto;
#ifdef _OPENMP
#pragma omp atomic read
#endif
        rc_visto = rc_geral;
        if (rc_visto != 0) continue;
        const std::size_t sq = static_cast<std::size_t>(q);
        for (std::size_t t = P.ini[sq]; t < P.ini[sq + 1]; t++) {
          const int rc = faz_linha(P.nos[t], xx, s1, m1);
          if (rc != 0) {
#ifdef _OPENMP
#pragma omp atomic write
#endif
            rc_geral = rc;
            break;
          }
        }
      }
    }
    if (falha(rc_geral)) return false;
    tt1 = diag_t ? agora() : 0.0;
    // a marca da thread 0 so tem selos das linhas dela (k + 1), que o topo nao repete
    for (std::size_t q = 0; q < P.topo.size(); q++) {
      checa(q);
      if (falha(faz_linha(P.topo[q], x, s, marca))) return false;
    }
  }

  const double tt2 = diag_t ? agora() : 0.0;
  if (hibrido) {
    const std::size_t T = T_cauda;
    // S = A_TT, guardada ja no lugar do fator: a coluna j da cauda recebe as linhas j..n-1
    for (std::size_t j = base; j < n; j++) {
      for (std::size_t r = 0; r < n - j; r++) {
        linha[colptr[j] + r] = static_cast<std::uint32_t>(j + r);
        valor[colptr[j] + r] = 0.0;
      }
      prox[j] = colptr[j + 1];
    }
    Cauda c{base, T, &colptr, &valor};
    for (std::size_t k = base; k < n; k++)
      for (std::size_t p = au.colptr[k]; p < au.colptr[k + 1]; p++) {
        const std::size_t i = au.linha[p];
        if (i >= base && i <= k) c.at(k - base, i - base) = au.valor[p];
      }
    // onde comecam as entradas da cauda em cada coluna de fora dela
    std::vector<std::size_t> ini(base);
    for (std::size_t j = 0; j < base; j++) {
      std::size_t q = prox[j];
      while (q > colptr[j] + 1 && linha[q - 1] >= base) q--;
      ini[j] = q;
    }
    // S -= W W': cada thread e dona de um ladrilho de COLUNAS de S e percorre as colunas de
    // W em ordem fixa, aplicando so o que cai nas suas
    const int nth = threads();
    const std::size_t nct = (T + LADRILHO - 1) / LADRILHO;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(nth)
#endif
    for (long cb = 0; cb < static_cast<long>(nct); cb++) {
      const std::size_t c0 = static_cast<std::size_t>(cb) * LADRILHO;
      const std::size_t c1 = std::min(T, c0 + LADRILHO);
      for (std::size_t j = 0; j < base; j++) {
        for (std::size_t qb = ini[j]; qb < prox[j]; qb++) {
          const std::size_t b = linha[qb] - base;
          if (b < c0) continue;
          if (b >= c1) break;
          const double wb = valor[qb];
          double* col_b = &c.at(b, b);
          for (std::size_t qa = qb; qa < prox[j]; qa++)
            col_b[linha[qa] - base - b] -= valor[qa] * wb;
        }
      }
    }
    if (!fatora_cauda(c, nth)) return false;
  }

  if (diag_t)
    REprintf("[t3 chol] n %zu cauda %zu: subarvores %.3f s (%zu), topo %.3f s (%zu nos), cauda %.3f s, threads %d\n",
             n, T_cauda, tt1 - tt0, P.vazia() ? std::size_t(0) : P.ini.size() - 1, tt2 - tt1,
             P.vazia() ? n : P.topo.size(), agora() - tt2, nth_part);
  L = Csc(n, n);
  L.linha.reserve(nz);
  L.valor.reserve(nz);
  for (std::size_t j = 0; j < n; j++) {
    for (std::size_t p = colptr[j]; p < prox[j]; p++) {
      L.linha.push_back(linha[p]);
      L.valor.push_back(valor[p]);
    }
    L.colptr[j + 1] = L.linha.size();
  }
  return true;
}

double logdet(const Csc& L) {
  double s = 0.0;
  for (std::size_t j = 0; j < L.ncol; j++) s += std::log(L.valor[L.colptr[j]]);
  return 2.0 * s;
}

// Resolve L L' x = b.
std::vector<double> resolve(const Csc& L, const std::vector<double>& b) {
  const std::size_t n = L.ncol;
  std::vector<double> x = b;
  for (std::size_t j = 0; j < n; j++) {          // para a frente
    x[j] /= L.valor[L.colptr[j]];
    for (std::size_t p = L.colptr[j] + 1; p < L.colptr[j + 1]; p++)
      x[L.linha[p]] -= L.valor[p] * x[j];
  }
  for (std::size_t j = n; j-- > 0;) {            // para tras
    for (std::size_t p = L.colptr[j] + 1; p < L.colptr[j + 1]; p++)
      x[j] -= L.valor[p] * x[L.linha[p]];
    x[j] /= L.valor[L.colptr[j]];
  }
  return x;
}

// Permuta simetricamente e devolve o triangulo SUPERIOR de P A P'.
//
// A ENTRADA TEM DE TER CADA PAR GUARDADO UMA VEZ SO, um triangulo, qualquer um deles.
// Duas armadilhas ja apareceram aqui, em direcoes opostas:
//
//  - com uma matriz CHEIA, (i,j) e (j,i) caem na mesma posicao permutada e a montagem por
//    triplos SOMA as duas: o resultado fica com o dobro fora da diagonal e a fatoracao passa
//    a dizer que a matriz nao e positiva-definida.
//  - filtrando por um triangulo fixo, uma entrada guardada no OUTRO triangulo e descartada
//    inteira. Foi o que aconteceu na primeira versao deste arquivo: a entrada vinha no
//    inferior, o filtro mantinha o superior, e sobrava so a diagonal. L saia diagonal, o
//    bloco denso nao aparecia, e a resolucao dava numero errado sem erro nenhum.
//
// Por isso nao ha filtro: cada entrada e mapeada uma vez e colocada no triangulo superior do
// indice novo, venha ela de onde vier.
std::uint64_t assinatura_padrao(const Csc& a) {
  std::uint64_t h = 1469598103934665603ULL;
  auto mistura = [&](std::uint64_t v) { h ^= v; h *= 1099511628211ULL; };
  mistura(a.nlin); mistura(a.ncol);
  for (std::size_t v : a.colptr) mistura(v);
  for (std::uint32_t v : a.linha) mistura(v);
  return h;
}

const Csc& permuta_cache(const Csc& a, CacheSimbolica& cs) {
  const std::uint64_t h = assinatura_padrao(a);
  if (h != cs.assinatura) {
    cs.assinatura = h;
    cs.vistos = 0;
    cs.com_mapa = false;
    cs.mapa.clear();
  }
  cs.vistos++;
  const bool cabe = a.nnz() < static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max());
  if (!cs.com_mapa && cs.vistos >= 2 && cabe) {
    // o mapa sai da propria permutacao aplicada ao INDICE de cada entrada; se a permutacao
    // somasse entradas (duplicatas) o numero de entradas cairia, e ai fica sem mapa
    Csc idx = a;
    for (std::size_t k = 0; k < idx.valor.size(); k++) idx.valor[k] = static_cast<double>(k);
    Csc p = permuta_sim(idx, cs.perm);
    if (p.nnz() == a.nnz()) {
      cs.mapa.assign(a.nnz(), 0);
      for (std::size_t q = 0; q < p.valor.size(); q++)
        cs.mapa[static_cast<std::size_t>(p.valor[q])] = static_cast<std::uint32_t>(q);
      cs.permutada = std::move(p);
      cs.com_mapa = true;
    }
  }
  if (!cs.com_mapa) {
    cs.permutada = permuta_sim(a, cs.perm);
    return cs.permutada;
  }
  for (std::size_t k = 0; k < a.valor.size(); k++) cs.permutada.valor[cs.mapa[k]] = a.valor[k];
  return cs.permutada;
}

Csc permuta_sim(const Csc& a, const std::vector<std::size_t>& perm) {
  const std::size_t n = a.ncol;
  if (perm.size() != n) throw Erro("permutation of the wrong size");
  std::vector<std::size_t> inv(n);
  for (std::size_t k = 0; k < n; k++) {
    if (perm[k] >= n) throw Erro("permutation out of range");
    inv[perm[k]] = k;
  }
  std::vector<std::uint32_t> li, cj;
  std::vector<double> v;
  li.reserve(a.nnz()); cj.reserve(a.nnz()); v.reserve(a.nnz());
  for (std::size_t c = 0; c < n; c++)
    for (std::size_t k = a.colptr[c]; k < a.colptr[c + 1]; k++) {
      const std::size_t r = a.linha[k];
      std::size_t i = inv[r], j = inv[c];
      if (i > j) std::swap(i, j);                // triangulo superior no novo indice
      li.push_back(static_cast<std::uint32_t>(i));
      cj.push_back(static_cast<std::uint32_t>(j));
      v.push_back(a.valor[k]);
    }
  return de_triplos(n, n, li, cj, v);
}

// Grau minimo (Tinney e Walker, 1967) sobre o grafo quociente.
//
// Medido, nao assumido. Num pedigree de 3.000 animais mais um clique genomico o preenchimento
// foi: natural 1.332.217, Cuthill-McKee reverso 746.017, Cuthill-McKee simples 2.506.473,
// grau minimo 201.487. O RCM (Cuthill e McKee, 1969) reduz pela metade mas poe os nos de grau alto PRIMEIRO, que e
// exatamente onde a forma fechada do bloco denso final nao os enxerga.
//
// O grafo quociente e o que torna isto pagavel: uma variavel eliminada vira um ELEMENTO que
// guarda o clique que ela criou, entao o preenchimento nunca e materializado.
std::vector<std::size_t> grau_minimo_exato(const Csc& a, double teto, double teto_trabalho,
                                           double* trabalho) {
  const std::size_t n = a.ncol;
  std::vector<std::size_t> ordem;
  if (n == 0) return ordem;
  ordem.reserve(n);

  std::vector<std::vector<std::uint32_t>> av(n), ev(n), le(n);
  for (std::size_t c = 0; c < n; c++)
    for (std::size_t k = a.colptr[c]; k < a.colptr[c + 1]; k++) {
      const std::size_t r = a.linha[k];
      if (r == c) continue;
      av[c].push_back(static_cast<std::uint32_t>(r));
      av[r].push_back(static_cast<std::uint32_t>(c));
    }
  for (auto& l : av) { std::sort(l.begin(), l.end()); l.erase(std::unique(l.begin(), l.end()), l.end()); }

  std::vector<char> vivo(n, 1);
  std::vector<std::uint32_t> marca(n, 0), no_lp(n, 0);
  std::uint32_t selo = 0, selo_lp = 0;

  // no ADIADO (linha densa, ver abaixo) sai do grafo e nao so da fila: fica vivo para ser
  // posto no fim da ordem, mas nao conta em grau, nao entra em L_p e nao segura absorcao.
  // Medido na APY (3.000 animais, 1.500 genotipados, nucleo de 414): com o nucleo adiado mas
  // ainda nas listas, cada jovem eliminado criava um elemento com os 414 dentro, e todo
  // viz() e toda checagem de absorcao varriam os 414 de novo. Foram 1,3e10 passos de viz()
  // e 5e9 de absorcao contra 7e7 no caminho exato, 17 s por ordenacao em vez de 0,3 s. E o
  // nucleo contado no grau e nos vivos arrastava 341 jovens para o bloco denso junto com
  // ele (752 colunas densas em vez de nucleo + 1). O AMD tira a linha densa do grafo
  // (Amestoy, Davis e Duff, 1996); aqui ela so saia da fila.
  std::vector<char> adiado(n, 0);
  std::size_t n_adiados = 0;
  auto ativo = [&](std::uint32_t j) { return vivo[j] && !adiado[j]; };

  // TRABALHO: entradas de lista lidas (vizinhos, membros de elemento, absorcao, limpeza).
  // E a medida deterministica do custo da ordenacao, que nao e a soma de c_j^2: num grafo de
  // muito enchimento cada pivo rele elementos que se sobrepoem, e medido num pedigree de
  // leitegadas a ordenacao andou a ~3e7 c_j^2 por segundo contra ~7e9 da fatoracao em
  // ladrilhos. Serve ao teto de quem pede (teto_trabalho) e a checagem de interrupcao.
  double trab = 0.0;

  // a lista de um elemento so perde membros (quem morre ou e adiado nao volta), entao ela
  // e compactada ao ser lida: cada entrada morta custa uma leitura so, e nao uma a cada
  // viz() e a cada checagem de absorcao que passar por ela
  auto compacta = [&](std::uint32_t e) {
    auto& l = le[e];
    trab += static_cast<double>(l.size());
    l.erase(std::remove_if(l.begin(), l.end(), [&](std::uint32_t x) { return !ativo(x); }),
            l.end());
  };

  // vizinhanca de uma variavel viva: vizinhos-variavel mais os membros de cada elemento
  auto viz = [&](std::size_t i) {
    selo++;
    std::vector<std::uint32_t> out;
    marca[i] = selo;
    trab += static_cast<double>(av[i].size() + ev[i].size());
    for (std::uint32_t j : av[i]) if (ativo(j) && marca[j] != selo) { marca[j] = selo; out.push_back(j); }
    for (std::uint32_t e : ev[i]) {
      compacta(e);
      trab += static_cast<double>(le[e].size());
      for (std::uint32_t j : le[e]) if (marca[j] != selo) { marca[j] = selo; out.push_back(j); }
    }
    return out;
  };

  std::vector<std::size_t> grau(n);
  for (std::size_t i = 0; i < n; i++) grau[i] = av[i].size();

  // LINHA DENSA SAI DO JOGO DE GRAUS (o "dense row handling" do AMD).
  //
  // Medido em dado real, 5.930 animais com 5.924 genotipados: o H^-1 sai 99,7% denso,
  // porque com genotipagem quase completa o A22^-1 e denso, e grau minimo num clique e o
  // PIOR caso desta estrutura de dados. Cada eliminacao chama viz() para cada membro de
  // L_p, e viz() custa O(k); com |L_p| = k dentro de k eliminacoes o custo e k^3 e nao nnz.
  // Foram 10 minutos parados aqui num ajuste que nem tinha chegado a primeira iteracao do
  // AI-REML, e a proxima chamada (a das MME, em aireml.cpp) e sobre uma matriz 2x maior.
  //
  // O conserto e nao deixar esses nos entrarem no jogo: quem passa do limiar vai direto
  // para o fim da ordem. Num clique todos passam, e o que sobra e exatamente a ordem que a
  // medicao de enchimento ja tinha apontado como otima aqui: elimina a parte esparsa por
  // grau minimo e empilha o clique no fim (dense_block = nucleo+1, enchimento ZERO).
  //
  // ORDENACAO NAO MUDA RESULTADO, so velocidade e enchimento: e semelhanca por permutacao.
  // O pior caso de errar o limiar e um ajuste lento, nunca um numero diferente. Isso e o
  // que torna esta troca barata de gatilhar e segura de errar.
  //
  // O LIMIAR E RELATIVO AO QUE AINDA ESTA VIVO, e nao a sqrt(n). Medido: com 10*sqrt(n), o
  // classico do AMD, a ordenacao ficou 8x a 13x mais rapida E PIOR, porque na APY o animal
  // JOVEM tem grau igual ao tamanho do nucleo (ele se liga a todo o nucleo) e passava do
  // limiar junto com o clique. O bloco denso saltou do nucleo para n_genotipados e o
  // enchimento subiu 7% a 12%. Isso troca custo de ordenacao, que e uma vez por ajuste, por
  // custo de fatoracao, que e TODA iteracao: com nucleo de 2.533 em 5.924 genotipados seriam
  // 12,8x mais flops por iteracao. Troca ruim.
  //
  // A assinatura da patologia nao e "grau alto", e "grau igual ao que sobrou": um clique que
  // abrange quase todos os vivos. Um no com grau 600 entre 1.200 vivos e denso e util; o
  // mesmo grau 600 entre 700 vivos e um clique e nao ha ordem que o salve.
  const double fracao = 0.8;
  auto denso_demais = [&](std::size_t d, std::size_t vivos) {
    return vivos > 64 && static_cast<double>(d) > fracao * static_cast<double>(vivos);
  };

  // NA PARTIDA, o corte relativo aos vivos chega tarde quando so parte dos animais e
  // genotipada. Medido na APY com 1.500 genotipados em 3.000: o nucleo (grau ~1.500, metade
  // dos vivos) so passava de 0,8 no passo 1.019, e ate la cada jovem eliminado atualizava
  // os 414 do nucleo, 2,7 s de ordenacao contra 0,13 s no caminho exato. O 10*sqrt(n) do
  // AMD pega o nucleo, mas sozinho tambem pega o JOVEM quando o nucleo e grande (o grau
  // dele e o tamanho do nucleo; ver acima). Os dois cortes juntos separam os dois casos:
  // o nucleo tem o grau do topo do grafo (liga-se a todo genotipado), o jovem tem o
  // tamanho do nucleo, e nos dois regimes medidos (nucleo 414 de 1.500, nucleo 2.533 de
  // 5.924) o jovem fica abaixo de 0,8 do topo. Num grafo esparso comum o corte absoluto
  // do AMD nunca e atingido e nada muda.
  //
  // O TOPO E O PERCENTIL 99 DO GRAU, e nao o maximo. Medido nas MME do mesmo caso: o
  // intercepto liga todo animal com registro (grau 2.700), virava o maximo, e 0,8 dele
  // passava do grau do nucleo; a ordenacao das MME ficava nos 2,5 s. Um punhado de
  // equacoes de efeito fixo nao pode definir o que e denso para os animais.
  std::size_t grau_topo = 0;
  {
    std::vector<std::size_t> g2(grau);
    const std::size_t pos = n - 1 - std::min(n - 1, n / 100);
    std::nth_element(g2.begin(), g2.begin() + pos, g2.end());
    grau_topo = g2[pos];
  }
  const double corte_abs = std::max(16.0, 10.0 * std::sqrt(static_cast<double>(n)));
  auto denso_na_partida = [&](std::size_t d) {
    return static_cast<double>(d) > corte_abs &&
           static_cast<double>(d) > fracao * static_cast<double>(grau_topo);
  };

  // A REGRA DE PARTIDA SO VALE QUANDO O CONJUNTO DENSO E PEQUENO. Tirar a linha densa do
  // grafo e o que o AMD faz, e e inocuo quando elas sao poucas: na APY o nucleo e 14% dos
  // nos e se liga a todo jovem, entao ignora-lo nao muda a ordem dos jovens. No passo unico
  // EXATO o clique genotipado e metade do grafo, e sem ele o grau minimo ordena os nao
  // genotipados sem saber quantos vizinhos genotipados cada um tem: medido, o ajuste exato
  // ficou ~1,5x mais lento por iteracao. Uma tentativa de devolver essa informacao contando
  // os vizinhos adiados na chave piorou a APY (bloco denso 1.012 -> 1.211) sem consertar o
  // exato, e foi descartada. Acima de 25% o conjunto fica com o corte relativo aos vivos,
  // que ja pega o clique quando os nao genotipados saem.
  {
    std::size_t na_partida = 0;
    for (std::size_t i = 0; i < n; i++) if (denso_na_partida(grau[i])) na_partida++;
    const bool vale_partida = na_partida <= n / 4;
    for (std::size_t i = 0; i < n; i++)
      if (denso_demais(grau[i], n) || (vale_partida && denso_na_partida(grau[i]))) {
        adiado[i] = 1; n_adiados++;
      }
  }
  std::vector<std::vector<std::uint32_t>> baldes(n + 1);
  for (std::size_t i = 0; i < n; i++)
    if (!adiado[i]) baldes[std::min(grau[i], n)].push_back(static_cast<std::uint32_t>(i));
  std::size_t lo = 0;

  // GRAU PREGUICOSO PARA OS HUBS. Medido no multicaracter com pedigree real de estrutura
  // comum, poucos pais para muitos filhos: a ordenacao era 98% do tempo de uma avaliacao e
  // crescia ~n^1.8 (0.06 s com 2 mil animais, 2.36 s com 16 mil), enquanto a Cholesky
  // numerica levava 0.014 s. A causa e o hub: um touro com centenas de filhos, ou a
  // equacao de um grupo contemporaneo ligada a todos os animais dele. A cada filho
  // eliminado, o touro estava em L_p, e o laco limpava a lista dele e recalculava o grau
  // dele com viz(), as duas coisas O(filhos). Sao n eliminacoes vezes O(n / touros).
  //
  // E o grau de um hub nao importa ate o fim: ele nunca e o minimo enquanto tem filhos.
  // Entao o no de lista grande so e MARCADO sujo em L_p, e a limpeza e o viz() acontecem
  // uma vez, quando ele e retirado do balde como candidato. No de lista pequena continua
  // sendo atualizado na hora, exatamente como antes, e onde nao ha hub a ordem sai igual.
  // Ordenacao nao muda resultado, so enchimento; o que se confere e que o enchimento nao
  // piorou.
  const std::size_t lista_grande = 16;
  // e so quando o elemento NOVO e pequeno. A preguica compensa porque cada filho do touro
  // traz um elemento de dois ou tres membros, e adiar a absorcao deles e barato. No nucleo
  // da APY cada jovem eliminado traz um elemento do tamanho do nucleo inteiro; adiado, isso
  // incha a lista ate o proximo recalculo, e medido, deixou a fatoracao da APY 2x MAIS
  // LENTA com enchimento identico. Com L_p grande vale o caminho antigo, exato.
  const std::size_t lp_pequeno = 32;
  std::vector<char> sujo(n, 0);
  double acumulado = 0.0;
  const double passo_checa = 6.4e7;
  double prox_checa = passo_checa;
  auto limpa = [&](std::size_t i) {
    trab += static_cast<double>(av[i].size() + ev[i].size());
    auto& a_i = av[i];
    a_i.erase(std::remove_if(a_i.begin(), a_i.end(),
                             [&](std::uint32_t x) { return !ativo(x); }), a_i.end());
    auto& e_i = ev[i];
    std::sort(e_i.begin(), e_i.end());
    e_i.erase(std::unique(e_i.begin(), e_i.end()), e_i.end());
  };

  while (ordem.size() < n) {
    std::size_t p = static_cast<std::size_t>(-1);
    while (lo <= n) {
      if (baldes[lo].empty()) { lo++; continue; }
      const std::uint32_t cand = baldes[lo].back();
      baldes[lo].pop_back();
      if (!vivo[cand] || adiado[cand] || std::min(grau[cand], n) != lo) continue;
      if (sujo[cand]) {
        // o hub chegou a frente com um grau velho: agora sim ele e limpo e medido, e volta
        // ao balde certo. Se o grau verdadeiro for menor que lo, a busca recomeca dali
        limpa(cand);
        sujo[cand] = 0;
        const std::size_t d = viz(cand).size();
        grau[cand] = d;
        if (denso_demais(d, n - ordem.size() - n_adiados)) { adiado[cand] = 1; n_adiados++; continue; }
        const std::size_t b = std::min(d, n);
        if (b != lo) {
          baldes[b].push_back(cand);
          if (b < lo) lo = b;
          continue;
        }
      }
      p = cand;
      break;
    }
    if (p == static_cast<std::size_t>(-1)) {          // o que sobrou esta isolado ou adiado
      for (std::size_t i = 0; i < n; i++) if (vivo[i]) { vivo[i] = 0; ordem.push_back(i); }
      break;
    }

    // INTERRUPCAO e TETO DE TRABALHO, conferidos a cada pivo: um so pivo de um grafo de muito
    // enchimento ja pode custar segundos, e a ordenacao inteira minutos
    if (trab >= prox_checa) {
      checa_interrupcao();
      prox_checa = trab + passo_checa;
    }
    if (teto_trabalho > 0.0 && trab > teto_trabalho) {
      if (trabalho) *trabalho = trab;
      return std::vector<std::size_t>();
    }
    std::vector<std::uint32_t> lp = viz(p);
    // TETO de c_j^2 (so quem pede, a rota automatica da dominancia): |L_p| + 1 e a contagem
    // da coluna p do fator SEM os adiados, entao a soma parcial de (|L_p| + 1)^2 nunca passa
    // a soma de c_j^2 desta ordem. Passou do teto, a fatoracao custaria mais que o teto e a
    // ordem nao interessa: devolve vazio.
    if (teto > 0.0) {
      const double c = static_cast<double>(lp.size()) + 1.0;
      acumulado += c * c;
      if (acumulado > teto) {
        if (trabalho) *trabalho = trab;
        return std::vector<std::size_t>();
      }
    }
    vivo[p] = 0;
    ordem.push_back(p);
    le[p] = lp;

    // pertencer a L_p por SELO, e nao varrendo L_p. A versao com varredura era quadratica
    // dentro de um laco que ja e sobre L_p, e a ordenacao virava a parte mais cara do ajuste.
    selo_lp++;
    for (std::uint32_t x : lp) no_lp[x] = selo_lp;

    for (std::uint32_t iu : lp) {
      const std::size_t i = iu;
      // O MESMO DESVIO TEM DE ESTAR NOS DOIS LACOS, e a primeira versao so o pos no de
      // baixo. Medido: com a guarda so la, o nucleo era adiado certo na inicializacao, mas
      // cada eliminacao de um NAO-nucleo ainda arrastava os k do nucleo por aqui, e av[i] de
      // um no do clique tem k entradas. Da O(k^2) por eliminacao, 2.533 x 5.923 num caso
      // real, e o ganho caiu para ~2x em vez dos ~1000x que a conta pedia.
      //
      // Pular isto e seguro por um invariante: av[i] e ev[i] so sao LIDOS por viz(i), e
      // viz(i) nunca roda para um no adiado, porque ele nao entra em balde (nao vira pivo) e
      // o laco de baixo o desvia. Manter a adjacencia de quem nunca mais sera consultado e
      // trabalho jogado fora. A absorcao le no_lp[x] para x em le[e], e um le[e] antigo
      // ainda pode trazer um no adiado depois dele; por isso ela filtra por ativo(), e
      // o no adiado nao segura a morte de elemento nenhum.
      if (adiado[i]) continue;
      if (lp.size() <= lp_pequeno && (sujo[i] || av[i].size() + ev[i].size() > lista_grande)) {
        // hub: so registra o elemento novo; limpeza e grau ficam para quando ele sair do
        // balde. viz() ja filtra os mortos, entao a lista suja nao muda o grau medido
        ev[i].push_back(static_cast<std::uint32_t>(p));
        sujo[i] = 1;
        continue;
      }
      if (sujo[i]) { limpa(i); sujo[i] = 0; }
      auto& a_i = av[i];
      a_i.erase(std::remove_if(a_i.begin(), a_i.end(),
                               [&](std::uint32_t x) { return x == p || !ativo(x); }), a_i.end());
      auto& e_i = ev[i];
      e_i.erase(std::remove(e_i.begin(), e_i.end(), static_cast<std::uint32_t>(p)), e_i.end());
      e_i.push_back(static_cast<std::uint32_t>(p));
      // absorcao: um elemento cujos membros vivos ja estao dentro de L_p esta morto
      trab += static_cast<double>(a_i.size() + e_i.size());
      e_i.erase(std::remove_if(e_i.begin(), e_i.end(), [&](std::uint32_t e) {
        if (e == p) return false;
        compacta(e);
        trab += static_cast<double>(le[e].size());
        for (std::uint32_t x : le[e]) if (no_lp[x] != selo_lp) return false;
        return true;
      }), e_i.end());
    }
    for (std::uint32_t iu : lp) {
      const std::size_t i = iu;
      if (adiado[i] || sujo[i]) continue;
      const std::size_t d = viz(i).size();
      grau[i] = d;
      // um no VIRA clique durante a eliminacao: o grau dele nao cresce, o numero de vivos e
      // que encolhe ate a razao passar do corte. E assim que o nucleo da APY e pego, depois
      // que os jovens ja sairam, e nao antes deles
      if (denso_demais(d, n - ordem.size() - n_adiados)) { adiado[i] = 1; n_adiados++; continue; }
      const std::size_t b = std::min(d, n);
      baldes[b].push_back(iu);
      if (b < lo) lo = b;
    }
  }
  if (trabalho) *trabalho = trab;
  return ordem;
}

}  // namespace br
