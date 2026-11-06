// ORDENACAO: grau minimo APROXIMADO (AMD) e a escolha entre ele e o grau minimo exato.
//
// O grau minimo exato (grau_minimo_exato, linalg.cpp) recalcula o grau de cada no tocado por um
// pivo pela uniao das listas (viz()), e e a etapa serial que mais pesou no E0: 190 s numa matriz
// de 1 milhao de colunas, a mesma com 1 ou 16 threads (out/velocidade/E0/NOTAS_E0.md). O AMD
// (Amestoy, Davis e Duff, 1996, SIAM J Matrix Anal Appl 17:886; 2004, ACM TOMS 30:381) troca o
// grau exato por um LIMITE SUPERIOR que sai sem uniao nenhuma:
//
//   d_i <= |A_i \ L_p| + |L_p \ i| + soma_{e em E_i, e != p} |L_e \ L_p|
//
// com |L_e \ L_p| = w(e), calculado de uma vez para todo elemento vizinho de L_p (o primeiro toque
// poe w(e) = |L_e| e cada membro de L_p o desconta). Junto vem o que o faz rapido na pratica:
//  - SUPERVARIAVEIS: duas variaveis com a mesma vizinhanca no grafo quociente sao uma so (nv);
//  - ELIMINACAO EM MASSA: um i de L_p cuja vizinhanca e so o elemento p sai junto com o pivo;
//  - ABSORCAO AGRESSIVA: um elemento com w(e) = 0 esta contido em L_p e morre.
//
// ORDENACAO NAO MUDA O RESULTADO, so o enchimento e o tempo: P C P' e semelhante a C. Mas o fator
// muda nos ultimos bits (outra ordem de soma), entao a troca passa por portao de exatidao com
// tolerancia, nunca bit a bit.
//
// AS LINHAS DENSAS seguem a MESMA regra do grau minimo exato (os mesmos dois cortes na partida e o
// corte relativo aos vivos durante a eliminacao), porque e dela que depende a cauda densa do passo
// unico e da APY: o no adiado sai do grafo e vai para o fim da ordem. Durante a eliminacao o corte
// e conferido com o grau EXATO (o aproximado e um limite superior e adiaria nos demais, e cada no a
// mais na cauda densa e custo de TODA iteracao; ver o comentario do limiar em linalg.cpp).
//
// UAITREM_ORDEM no ambiente escolhe: "gm" (grau minimo exato) ou "amd". Quem pede teto de custo
// (a dominancia) fica sempre no exato, que e quem sabe medir o trabalho.

#include "mme.h"
#include <cstdlib>
#include <cstring>

namespace br {

namespace {

enum : unsigned char { VAR = 0, NAO_PRINC = 1, ELEM = 2, MORTO = 3, ADIADO = 4 };

}  // namespace

std::vector<std::size_t> ordem_amd(const Csc& a) {
  const std::size_t n = a.ncol;
  std::vector<std::size_t> ordem;
  if (n == 0) return ordem;
  ordem.reserve(n);

  // grafo simetrico sem a diagonal, listas ordenadas e sem repeticao (como no exato)
  std::vector<std::vector<std::uint32_t>> A(n), E(n), Le(n), filhos(n);
  for (std::size_t c = 0; c < n; c++)
    for (std::size_t k = a.colptr[c]; k < a.colptr[c + 1]; k++) {
      const std::size_t r = a.linha[k];
      if (r == c) continue;
      A[c].push_back(static_cast<std::uint32_t>(r));
      A[r].push_back(static_cast<std::uint32_t>(c));
    }
  for (auto& l : A) { std::sort(l.begin(), l.end()); l.erase(std::unique(l.begin(), l.end()), l.end()); }

  std::vector<unsigned char> st(n, VAR);
  std::vector<std::int64_t> nv(n, 1), grau(n, 0), eldeg(n, 0), w(n, 0);
  std::vector<std::uint32_t> marca(n, 0), wselo(n, 0);
  std::uint32_t selo = 0, selo_w = 0;

  // ---- linhas densas: os mesmos cortes do grau minimo exato (linalg.cpp)
  const double fracao = 0.8;
  auto denso_demais = [&](double d, double vivos) { return vivos > 64 && d > fracao * vivos; };
  {
    std::vector<std::size_t> g0(n);
    for (std::size_t i = 0; i < n; i++) g0[i] = A[i].size();
    std::size_t grau_topo = 0;
    {
      std::vector<std::size_t> g2(g0);
      const std::size_t pos = n - 1 - std::min(n - 1, n / 100);
      std::nth_element(g2.begin(), g2.begin() + pos, g2.end());
      grau_topo = g2[pos];
    }
    const double corte_abs = std::max(16.0, 10.0 * std::sqrt(static_cast<double>(n)));
    auto denso_na_partida = [&](std::size_t d) {
      return static_cast<double>(d) > corte_abs &&
             static_cast<double>(d) > fracao * static_cast<double>(grau_topo);
    };
    std::size_t na_partida = 0;
    for (std::size_t i = 0; i < n; i++) if (denso_na_partida(g0[i])) na_partida++;
    const bool vale_partida = na_partida <= n / 4;
    for (std::size_t i = 0; i < n; i++)
      if (denso_demais(static_cast<double>(g0[i]), static_cast<double>(n)) ||
          (vale_partida && denso_na_partida(g0[i])))
        st[i] = ADIADO;
  }
  std::int64_t n_adiado = 0;
  for (std::size_t i = 0; i < n; i++) if (st[i] == ADIADO) n_adiado++;
  for (std::size_t i = 0; i < n; i++) {
    if (st[i] != VAR) continue;
    std::int64_t d = 0;
    for (std::uint32_t v : A[i]) if (st[v] == VAR) d++;
    grau[i] = d;
  }

  std::vector<std::vector<std::uint32_t>> baldes(n + 1);
  for (std::size_t i = 0; i < n; i++)
    if (st[i] == VAR) baldes[std::min<std::size_t>(static_cast<std::size_t>(grau[i]), n)].push_back(
        static_cast<std::uint32_t>(i));
  std::size_t lo = 0;
  std::int64_t nel = 0;   // peso ja eliminado

  // a saida: o pivo e, em seguida, quem foi fundido nele (supervariavel ou eliminacao em massa)
  std::vector<std::uint32_t> pilha;
  auto emite = [&](std::uint32_t p) {
    pilha.clear();
    pilha.push_back(p);
    while (!pilha.empty()) {
      const std::uint32_t q = pilha.back();
      pilha.pop_back();
      ordem.push_back(q);
      for (std::size_t k = filhos[q].size(); k-- > 0;) pilha.push_back(filhos[q][k]);
      std::vector<std::uint32_t>().swap(filhos[q]);
    }
  };

  // grau EXATO externo de i (uniao), so para o corte das linhas densas
  auto grau_exato = [&](std::uint32_t i) {
    selo++;
    marca[i] = selo;
    std::int64_t d = 0;
    for (std::uint32_t v : A[i]) if (st[v] == VAR && marca[v] != selo) { marca[v] = selo; d += nv[v]; }
    for (std::uint32_t e : E[i]) {
      if (st[e] != ELEM) continue;
      for (std::uint32_t v : Le[e]) if (st[v] == VAR && marca[v] != selo) { marca[v] = selo; d += nv[v]; }
    }
    return d;
  };

  std::vector<std::uint32_t> lp, novo;
  std::vector<std::uint64_t> hash;
  std::vector<char> massa;
  double trab = 0.0, prox_checa = 6.4e7;

  while (true) {
    std::int64_t p = -1;
    while (lo <= n) {
      if (baldes[lo].empty()) { lo++; continue; }
      const std::uint32_t c = baldes[lo].back();
      baldes[lo].pop_back();
      if (st[c] != VAR || std::min<std::size_t>(static_cast<std::size_t>(grau[c]), n) != lo) continue;
      p = c;
      break;
    }
    if (p < 0) break;
    if (trab >= prox_checa) { checa_interrupcao(); prox_checa = trab + 6.4e7; }
    const std::uint32_t pu = static_cast<std::uint32_t>(p);

    // L_p: os vizinhos variaveis e os membros dos elementos vizinhos; esses elementos morrem em p
    selo++;
    marca[pu] = selo;
    lp.clear();
    for (std::uint32_t e : E[pu]) {
      if (st[e] != ELEM) continue;
      trab += static_cast<double>(Le[e].size());
      for (std::uint32_t v : Le[e])
        if (st[v] == VAR && marca[v] != selo) { marca[v] = selo; lp.push_back(v); }
      st[e] = MORTO;
      std::vector<std::uint32_t>().swap(Le[e]);
    }
    trab += static_cast<double>(A[pu].size());
    for (std::uint32_t v : A[pu])
      if (st[v] == VAR && marca[v] != selo) { marca[v] = selo; lp.push_back(v); }
    std::vector<std::uint32_t>().swap(A[pu]);
    std::vector<std::uint32_t>().swap(E[pu]);
    st[pu] = ELEM;
    nel += nv[pu];
    std::int64_t degme = 0;
    for (std::uint32_t v : lp) degme += nv[v];

    // w(e) = |L_e \ L_p| para todo elemento vizinho de L_p
    selo_w++;
    for (std::uint32_t i : lp)
      for (std::uint32_t e : E[i]) {
        if (st[e] != ELEM) continue;
        if (wselo[e] != selo_w) { wselo[e] = selo_w; w[e] = eldeg[e]; }
        w[e] -= nv[i];
      }

    // o limite de grau de cada i de L_p, a poda das listas e a absorcao agressiva
    const std::int64_t vivos_peso = static_cast<std::int64_t>(n) - nel - n_adiado;
    hash.assign(lp.size(), 0);
    massa.assign(lp.size(), 0);
    std::int64_t em_massa = 0;
    for (std::size_t q = 0; q < lp.size(); q++) {
      const std::uint32_t i = lp[q];
      std::int64_t deg_e = 0;
      std::uint64_t h = 0;
      novo.clear();
      trab += static_cast<double>(E[i].size() + A[i].size());
      for (std::uint32_t e : E[i]) {
        if (st[e] != ELEM) continue;
        const std::int64_t we = std::max<std::int64_t>(0, w[e]);
        if (we == 0) { st[e] = MORTO; std::vector<std::uint32_t>().swap(Le[e]); continue; }
        novo.push_back(e);
        deg_e += we;
        h += e;
      }
      novo.push_back(pu);
      h += pu;
      E[i].swap(novo);
      std::int64_t deg_a = 0;
      novo.clear();
      for (std::uint32_t v : A[i]) {
        if (st[v] != VAR || v == i || marca[v] == selo) continue;
        novo.push_back(v);
        deg_a += nv[v];
        h += v;
      }
      A[i].swap(novo);
      std::int64_t d = deg_a + deg_e + (degme - nv[i]);
      d = std::min(d, grau[i] + (degme - nv[i]));
      d = std::min(d, vivos_peso - nv[i]);
      grau[i] = std::max<std::int64_t>(0, d);
      hash[q] = h;
      if (A[i].empty() && E[i].size() == 1) { massa[q] = 1; em_massa += nv[i]; }
    }

    // eliminacao em massa: i so ve o elemento p, e indistinguivel do pivo
    if (em_massa > 0) {
      std::size_t m = 0;
      for (std::size_t q = 0; q < lp.size(); q++) {
        const std::uint32_t i = lp[q];
        if (massa[q]) {
          st[i] = NAO_PRINC;
          filhos[pu].push_back(i);
          nel += nv[i];
          std::vector<std::uint32_t>().swap(A[i]);
          std::vector<std::uint32_t>().swap(E[i]);
          continue;
        }
        lp[m] = i; hash[m] = hash[q]; m++;
      }
      lp.resize(m); hash.resize(m);
      degme -= em_massa;
      for (std::uint32_t i : lp) grau[i] = std::max<std::int64_t>(0, grau[i] - em_massa);
    }

    // supervariaveis: mesma lista de variaveis e de elementos (ordenadas so quando o hash bate)
    if (lp.size() >= 2) {
      std::vector<std::size_t> ix(lp.size());
      for (std::size_t q = 0; q < ix.size(); q++) ix[q] = q;
      std::sort(ix.begin(), ix.end(), [&](std::size_t x, std::size_t y) {
        return hash[x] != hash[y] ? hash[x] < hash[y] : lp[x] < lp[y];
      });
      for (std::size_t r0 = 0; r0 < ix.size();) {
        std::size_t r1 = r0 + 1;
        while (r1 < ix.size() && hash[ix[r1]] == hash[ix[r0]]) r1++;
        if (r1 - r0 >= 2) {
          for (std::size_t x = r0; x < r1; x++) {
            const std::uint32_t i = lp[ix[x]];
            std::sort(A[i].begin(), A[i].end());
            std::sort(E[i].begin(), E[i].end());
          }
          for (std::size_t x = r0; x < r1; x++) {
            const std::uint32_t i = lp[ix[x]];
            if (st[i] != VAR) continue;
            for (std::size_t y = x + 1; y < r1; y++) {
              const std::uint32_t j = lp[ix[y]];
              if (st[j] != VAR) continue;
              if (A[i] != A[j] || E[i] != E[j]) continue;
              nv[i] += nv[j];
              grau[i] = std::max<std::int64_t>(0, grau[i] - nv[j]);
              nv[j] = 0;
              st[j] = NAO_PRINC;
              filhos[i].push_back(j);
              std::vector<std::uint32_t>().swap(A[j]);
              std::vector<std::uint32_t>().swap(E[j]);
            }
          }
        }
        r0 = r1;
      }
    }

    // o elemento novo e os graus de volta aos baldes (com o corte das linhas densas)
    std::int64_t eld = 0;
    novo.clear();
    for (std::uint32_t i : lp) if (st[i] == VAR) { novo.push_back(i); eld += nv[i]; }
    Le[pu].swap(novo);
    eldeg[pu] = eld;
    const double vivos = static_cast<double>(static_cast<std::int64_t>(n) - nel - n_adiado);
    for (std::uint32_t i : Le[pu]) {
      if (st[i] != VAR) continue;
      if (denso_demais(static_cast<double>(grau[i]), vivos)) {
        const std::int64_t ex = grau_exato(i);
        if (denso_demais(static_cast<double>(ex), vivos)) {
          st[i] = ADIADO;
          n_adiado += nv[i];
          for (std::uint32_t e : E[i]) if (st[e] == ELEM) eldeg[e] -= nv[i];
          continue;
        }
        grau[i] = ex;
      }
      const std::size_t b = std::min<std::size_t>(static_cast<std::size_t>(grau[i]), n);
      baldes[b].push_back(i);
      if (b < lo) lo = b;
    }
    emite(pu);
  }

  // o que sobrou (adiados, e o que nunca entrou em balde) vai para o fim, na ordem dos indices
  for (std::size_t i = 0; i < n; i++)
    if (st[i] == VAR || st[i] == ADIADO) emite(static_cast<std::uint32_t>(i));
  if (ordem.size() != n) throw Erro("internal: the AMD ordering lost columns");
  return ordem;
}

int modo_ordem() {
  const char* e = std::getenv("UAITREM_ORDEM");
  if (e && std::strcmp(e, "amd") == 0) return 1;
  if (e && std::strcmp(e, "gm") == 0) return 0;
  return 0;   // padrao: o exato, ate a medida decidir
}

std::vector<std::size_t> grau_minimo(const Csc& a, double teto, double teto_trabalho,
                                     double* trabalho) {
  if (teto > 0.0 || teto_trabalho > 0.0 || trabalho != nullptr || modo_ordem() == 0)
    return grau_minimo_exato(a, teto, teto_trabalho, trabalho);
  return ordem_amd(a);
}

}  // namespace br
