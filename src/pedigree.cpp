// Pedigree: ordenacao topologica, endogamia de Meuwissen e Luo (1992), e A^-1 de Henderson (1976).
//
// Este e o modulo onde um erro passa despercebido com mais facilidade, porque A^-1 sai
// simetrica, positiva-definida e plausivel mesmo quando esta errada. Dois pontos concentram
// isso:
//
//  1. a VARIANCIA MENDELIANA depende de quantos pais sao conhecidos, 1/2 com os dois,
//     3/4 com um, 1 com nenhum, e a endogamia dos pais entra nela. Usar 1/2 sempre e o erro
//     classico, e ele so aparece num pedigree com endogamia de verdade.
//  2. as contribuicoes ACUMULAM na mesma posicao: o d/4 entre um casal chega por cada filho.
//     Uma montagem em que o ultimo valor vence produz uma matriz que ainda converge, para o
//     lugar errado.

#include "mme.h"
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <cmath>
#include <algorithm>

namespace br {

// Constroi o pedigree em ordem topologica.
//
// A ordenacao e por profundidade ITERATIVA e nao recursiva: um pedigree real tem dezenas de
// milhares de animais e a recursao estoura a pilha muito antes disso.
Pedigree constroi_pedigree(const std::vector<std::string>& id0,
                           const std::vector<std::string>& pa0,
                           const std::vector<std::string>& ma0,
                           const std::vector<std::string>& mf,
                           const std::vector<double>& gama,
                           const std::vector<char>& mgs0) {
  if (!mgs0.empty() && mgs0.size() != id0.size())
    throw Erro("the maternal-grandsire flags must have one entry per pedigree row");
  // gamma = sem metafounders = era ignorado calado, e o ajuste saia na A comum
  if (mf.empty() && !gama.empty())
    throw Erro("gamma = was given without metafounders =: Gamma is the relationship among the "
               "metafounders, so it needs their labels (metafounders = c(...)), and without "
               "them the fit would use the plain A and ignore gamma");
  if (!mf.empty() && !mgs0.empty())
    for (char c : mgs0)
      if (c) throw Erro("a sire / maternal-grandsire pedigree does not take metafounders: "
                        "no reference program combines the two, and the grandsire path of "
                        "Henderson (1975) assumes an unknown, unrelated maternal granddam");
  // Gamma chega de duas formas: comprimento q (a DIAGONAL, compatibilidade) ou q*q (a
  // matriz CHEIA, por linhas). O caso diagonal e so um atalho de escrita: ele e expandido
  // aqui e daqui para baixo existe um caminho unico.
  const std::size_t q = mf.size();
  std::vector<double> G(q * q, 0.0);
  if (q) {
    if (gama.size() == q) {
      for (std::size_t k = 0; k < q; k++) G[k * q + k] = gama[k];
    } else if (gama.size() == q * q) {
      G = gama;
    } else {
      throw Erro("gamma must have one entry per metafounder (the diagonal) or "
                 "n_metafounders^2 entries (the full matrix, by rows)");
    }
    // Simetria: exigida, nao presumida. Aceitar so o triangulo de cima em silencio seria a
    // porta de entrada de um Gamma pela metade.
    double assim = 0.0, esc = 0.0;
    for (std::size_t a = 0; a < q; a++)
      for (std::size_t b = 0; b < q; b++) {
        assim = std::max(assim, std::fabs(G[a * q + b] - G[b * q + a]));
        esc = std::max(esc, std::fabs(G[a * q + b]));
      }
    if (assim > 1e-10 * std::max(esc, 1.0)) throw Erro("gamma is not symmetric");
    for (std::size_t a = 0; a < q; a++)
      for (std::size_t b = 0; b < a; b++) {
        const double m2 = 0.5 * (G[a * q + b] + G[b * q + a]);
        G[a * q + b] = m2; G[b * q + a] = m2;
      }
    // A diagonal e limitada por CIMA, e nao por baixo. d_i = 1 - 0,25(gamma_ss + gamma_tt)
    // e o pior caso e os dois pais serem o MESMO metafundador, entao gamma_ii < 2 basta.
    // gamma_ii = 0 e legitimo (e o limite de grupo de pais desconhecidos) e o fora da
    // diagonal pode ser NEGATIVO (bases divergidas por selecao em sentidos opostos), entao
    // a regra antiga "0 < gamma < 2 elemento a elemento" estava errada nos dois extremos.
    for (std::size_t k = 0; k < q; k++)
      if (!(G[k * q + k] < 2.0) || !(G[k * q + k] >= 0.0))
        throw Erro("a diagonal entry of gamma is outside [0, 2): with gamma_ii >= 2 the "
                   "Mendelian variance of a metafounder's offspring is not positive");
  }

  // metafundadores viram linhas-base VIRTUAIS, prefixadas ao pedigree
  std::vector<std::string> id(mf);
  id.insert(id.end(), id0.begin(), id0.end());
  std::vector<std::string> pa(mf.size(), "0"), ma(mf.size(), "0");
  pa.insert(pa.end(), pa0.begin(), pa0.end());
  ma.insert(ma.end(), ma0.begin(), ma0.end());

  const std::size_t n = id.size();
  if (pa.size() != n || ma.size() != n) throw Erro("pedigree with columns of different lengths");
  if (n == 0) throw Erro("empty pedigree");

  std::unordered_map<std::string, std::size_t> pos;
  pos.reserve(n * 2);
  for (std::size_t i = 0; i < n; i++) {
    if (id[i].empty() || id[i] == "0") throw Erro("there is an animal with an empty or '0' identifier");
    if (!pos.emplace(id[i], i).second)
      throw Erro(i < mf.size()
                     ? "metafounder '" + id[i] + "' repeats or collides with an animal id"
                     : "repeated animal in the pedigree: '" + id[i] + "'");
  }

  auto liga = [&](const std::string& s) -> std::int64_t {
    if (s.empty() || s == "0" || s == "NA") return -1;
    auto it = pos.find(s);
    // Um pai citado e ausente NAO pode virar desconhecido em silencio: isso muda a variancia
    // mendeliana do filho e o parentesco de toda a descendencia.
    if (it == pos.end())
      throw Erro("parent '" + s + "' is cited and has no line of its own in the pedigree" +
                 dica_cientifica(s, id));
    return static_cast<std::int64_t>(it->second);
  };
  // varios pais citados sem linha (um pedigree recortado por data): a mensagem conta todos e
  // cita os primeiros, e nao so o primeiro, que fazia o usuario corrigir um por vez
  {
    std::vector<std::string> faltam;
    std::unordered_set<std::string> vistos;
    for (std::size_t i = 0; i < n; i++)
      for (const std::string* s : {&pa[i], &ma[i]})
        if (!(s->empty() || *s == "0" || *s == "NA") && pos.find(*s) == pos.end() &&
            vistos.insert(*s).second)
          faltam.push_back(*s);
    if (faltam.size() > 1) {
      std::string lista;
      for (std::size_t k = 0; k < faltam.size() && k < 5; k++)
        lista += (k ? ", '" : "'") + faltam[k] + "'";
      if (faltam.size() > 5) lista += ", and " + std::to_string(faltam.size() - 5) + " more";
      throw Erro(std::to_string(faltam.size()) + " parents are cited and have no line of their "
                 "own in the pedigree (" + lista + "): give each one a line, with its own "
                 "parents or 0 for unknown" + dica_cientifica(faltam[0], id));
    }
  }
  // o mapa guarda a COLUNA de Gamma, nao o valor: com Gamma cheia e a coluna que importa
  std::unordered_map<std::string, std::size_t> gmap;
  for (std::size_t k = 0; k < q; k++) gmap[mf[k]] = k;

  std::vector<std::int64_t> p(n), m(n);
  for (std::size_t i = 0; i < n; i++) { p[i] = liga(pa[i]); m[i] = liga(ma[i]); }

  // profundidade por varredura iterativa, com deteccao de ciclo
  std::vector<int> estado(n, 0);         // 0 nao visitado, 1 na pilha, 2 pronto
  std::vector<std::size_t> ordem;
  ordem.reserve(n);
  std::vector<std::pair<std::size_t, int>> pilha;
  for (std::size_t s = 0; s < n; s++) {
    if (estado[s]) continue;
    pilha.push_back({s, 0});
    while (!pilha.empty()) {
      auto& [v, fase] = pilha.back();
      if (fase == 0) {
        if (estado[v] == 1) throw Erro("the pedigree has a cycle at '" + id[v] + "'");
        if (estado[v] == 2) { pilha.pop_back(); continue; }
        estado[v] = 1;
        fase = 1;
        if (p[v] >= 0) { pilha.push_back({static_cast<std::size_t>(p[v]), 0}); continue; }
      }
      if (fase == 1) {
        fase = 2;
        if (m[v] >= 0) { pilha.push_back({static_cast<std::size_t>(m[v]), 0}); continue; }
      }
      estado[v] = 2;
      ordem.push_back(v);
      pilha.pop_back();
    }
  }

  std::vector<std::size_t> novo(n);
  for (std::size_t k = 0; k < ordem.size(); k++) novo[ordem[k]] = k;

  Pedigree out;
  out.ids.resize(n);
  out.pai.resize(n);
  out.mae.resize(n);
  out.eh_mf.assign(n, 0);
  out.gama.assign(n, 0.0);
  out.col_mf.assign(n, -1);
  out.n_mf = q;
  bool tem_mgs = false;
  for (char c : mgs0) if (c) { tem_mgs = true; break; }
  if (tem_mgs) out.mgs.assign(n, 0);
  for (std::size_t k = 0; k < ordem.size(); k++) {
    const std::size_t v = ordem[k];
    out.ids[k] = id[v];
    // a flag vai com a LINHA: v e o indice de entrada (metafundadores na frente, que nao
    // existem aqui quando ha flag), entao a linha original do usuario e v - q
    if (tem_mgs && v >= q && mgs0[v - q]) out.mgs[k] = 1;
    out.pai[k] = p[v] >= 0 ? static_cast<std::int64_t>(novo[p[v]]) : -1;
    out.mae[k] = m[v] >= 0 ? static_cast<std::int64_t>(novo[m[v]]) : -1;
    auto it = gmap.find(id[v]);
    if (it != gmap.end()) {
      out.eh_mf[k] = 1;
      out.col_mf[k] = static_cast<std::int64_t>(it->second);
      out.gama[k] = G[it->second * q + it->second];
    }
  }

  // K (Cholesky inferior) e Gamma^-1, uma vez so. Gamma tem tipicamente 2 a 10 linhas,
  // entao o custo e irrelevante e vale pre-computar em vez de refazer por animal.
  //
  // A Cholesky e tambem o TESTE de admissibilidade, e e o teste certo: uma Gamma com
  // gamma_jk grande demais passa folgado no criterio de variancia mendeliana positiva e
  // ainda assim deixa A(Gamma) indefinida. Testar so d > 0 nao pega isso.
  if (q) {
    // Gamma pode ser SINGULAR e ainda assim significar alguma coisa: gamma_ii = 0 e o
    // limite de grupo de pais desconhecidos, e dois metafundadores para a MESMA populacao
    // dao linhas identicas. O artigo manda usar INVERSA GENERALIZADA nesses casos, e diz
    // que com Gamma = 0 a pseudo-inversa reproduz exatamente a A^-1 com grupos de pais
    // desconhecidos. Entao a rota e a decomposicao espectral, nao a Cholesky: ela aceita
    // autovalor nulo (inverte so o que e positivo, que E a pseudo-inversa de Moore-Penrose)
    // e continua RECUSANDO autovalor negativo, que e a inadmissibilidade de verdade.
    //
    // K, usada na endogamia por ||K' l||^2, sai da mesma decomposicao: K = U sqrt(L). Ela
    // nao precisa ser triangular, so satisfazer K K' = Gamma, e a forma quadratica nao
    // distingue as duas.
    std::vector<double> ev(q), U(q * q);
    {
      // Jacobi ciclico: q e tipicamente 2 a 10, entao a simplicidade vale mais que a
      // velocidade, e ele e exato em simetrica sem depender de LAPACK aqui.
      std::vector<double> Aw = G;
      for (std::size_t a = 0; a < q; a++)
        for (std::size_t b = 0; b < q; b++) U[a * q + b] = (a == b) ? 1.0 : 0.0;
      for (int varr = 0; varr < 100; varr++) {
        double fora = 0.0;
        for (std::size_t a = 0; a < q; a++)
          for (std::size_t b = a + 1; b < q; b++) fora += Aw[a * q + b] * Aw[a * q + b];
        if (fora < 1e-30) break;
        for (std::size_t a = 0; a < q; a++)
          for (std::size_t b = a + 1; b < q; b++) {
            if (std::fabs(Aw[a * q + b]) < 1e-300) continue;
            const double th = 0.5 * (Aw[b * q + b] - Aw[a * q + a]) / Aw[a * q + b];
            const double t = (th >= 0 ? 1.0 : -1.0) /
                             (std::fabs(th) + std::sqrt(th * th + 1.0));
            const double c = 1.0 / std::sqrt(t * t + 1.0), s2 = t * c;
            for (std::size_t k2 = 0; k2 < q; k2++) {
              const double ak = Aw[a * q + k2], bk = Aw[b * q + k2];
              Aw[a * q + k2] = c * ak - s2 * bk;
              Aw[b * q + k2] = s2 * ak + c * bk;
            }
            for (std::size_t k2 = 0; k2 < q; k2++) {
              const double ka = Aw[k2 * q + a], kb = Aw[k2 * q + b];
              Aw[k2 * q + a] = c * ka - s2 * kb;
              Aw[k2 * q + b] = s2 * ka + c * kb;
              const double ua = U[k2 * q + a], ub = U[k2 * q + b];
              U[k2 * q + a] = c * ua - s2 * ub;
              U[k2 * q + b] = s2 * ua + c * ub;
            }
          }
      }
      for (std::size_t a = 0; a < q; a++) ev[a] = Aw[a * q + a];
    }
    double maior = 0.0;
    for (double e : ev) maior = std::max(maior, std::fabs(e));
    const double corte = 1e-12 * std::max(maior, 1.0);
    for (double e : ev)
      if (e < -corte)
        throw Erro("gamma is not positive semi-definite: it has a negative eigenvalue, so "
                   "the relationship matrix it generates is not a covariance matrix at "
                   "all. A base relationship matrix may be singular, which is the "
                   "unknown-parent-group limit and is handled, but it may not be "
                   "indefinite");
    out.gama_chol.assign(q * q, 0.0);
    out.gama_inv.assign(q * q, 0.0);
    for (std::size_t a = 0; a < q; a++)
      for (std::size_t c = 0; c < q; c++) {
        const double raiz = ev[c] > corte ? std::sqrt(ev[c]) : 0.0;
        out.gama_chol[a * q + c] = U[a * q + c] * raiz;
      }
    for (std::size_t a = 0; a < q; a++)
      for (std::size_t b = 0; b < q; b++) {
        double s3 = 0.0;
        for (std::size_t c = 0; c < q; c++)
          if (ev[c] > corte) s3 += U[a * q + c] * U[b * q + c] / ev[c];
        out.gama_inv[a * q + b] = s3;
      }
  }
  return out;
}

// Variancia mendeliana, UMA definicao para a endogamia e a A^-1 (antes a formula vivia
// duplicada nas duas e podia divergir). Pai/mae: 1/2 - (F_s + F_d)/4, 3/4 - F/4 com um so,
// 1 sem nenhum. Pai/avo materno (Henderson 1975, 1976): u_i = u_s/2 + u_k/4 + m_i com a avo
// desconhecida, logo 11/16 - F_s/4 - F_k/16, 3/4 - F_s/4 so com o pai, 15/16 - F_k/16 so
// com o avo, 1 sem nenhum.
double variancia_mendeliana(const Pedigree& p, const std::vector<double>& f, std::size_t i) {
  const std::int64_t s = p.pai[i], t = p.mae[i];
  const bool g = !p.mgs.empty() && p.mgs[i];
  if (g) {
    if (s >= 0 && t >= 0) return 11.0 / 16.0 - 0.25 * f[s] - f[t] / 16.0;
    if (s >= 0) return 0.75 - 0.25 * f[s];
    if (t >= 0) return 15.0 / 16.0 - f[t] / 16.0;
    return 1.0;
  }
  if (s >= 0 && t >= 0) return 0.5 - 0.25 * (f[s] + f[t]);
  if (s >= 0) return 0.75 - 0.25 * f[s];
  if (t >= 0) return 0.75 - 0.25 * f[t];
  return 1.0;
}

// Endogamia por Meuwissen e Luo (1992).
//
// Para cada animal com os dois pais conhecidos, sobe o ramo dos ancestrais comuns em vez de
// formar A. O custo e proporcional ao numero de ancestrais, nao a n^2.
std::vector<double> endogamia(const Pedigree& p) {
  const std::size_t n = p.ids.size();
  std::vector<double> f(n, 0.0);
  std::vector<double> l(n, 0.0);
  std::vector<double> d(n, 0.0);

  // variancia mendeliana de cada animal, ja com a endogamia dos pais
  auto mendel = [&](std::size_t i) { return variancia_mendeliana(p, f, i); };

  std::priority_queue<std::size_t> fila;   // maior indice primeiro: o pedigree esta ordenado
  std::vector<char> na_fila(n, 0);

  for (std::size_t i = 0; i < n; i++) {
    if (!p.eh_mf.empty() && p.eh_mf[i]) {
      // metafundador: linha-base com a_ff = gamma_ii, logo F = gamma_ii - 1 e D = gamma_ii.
      // O F e NEGATIVO quando gamma < 1 e tem de continuar negativo: a formula mendeliana
      // le exatamente esse valor.
      f[i] = p.gama[i] - 1.0;
      d[i] = p.gama[i];
      continue;
    }
    d[i] = mendel(i);
    const std::int64_t s = p.pai[i], t = p.mae[i];
    if (s < 0 || t < 0) { f[i] = 0.0; continue; }

    // l ja esta todo em zero aqui: cada entrada que o animal anterior tocou entrou na fila e foi
    // zerada ao sair dela (l[j] = 0.0), e a fila termina vazia. Zerar o vetor inteiro a cada
    // animal custava 8 n bytes por animal, O(n^2) no pedigree, e era a primeira parede de
    // tempo num pedigree de 1 milhao de animais ou mais.
    l[i] = 1.0;
    fila.push(i);
    na_fila[i] = 1;
    double soma = 0.0;

    // Os l dos metafundadores sao COLHIDOS, nao consumidos um a um. Ler cada metafundador
    // como um ancestral qualquer, somando l_j^2 gamma_jj e zerando l_j, faz os termos
    // cruzados gamma_jk desaparecerem EM SILENCIO e reproduz o comportamento diagonal sem
    // erro visivel nenhum. O termo certo e a forma quadratica inteira, l_mf' Gamma l_mf,
    // computada como ||K' l_mf||^2, que e a forma publicada e a estavel.
    std::vector<double> lmf(p.n_mf, 0.0);

    while (!fila.empty()) {
      const std::size_t j = fila.top();
      fila.pop();
      na_fila[j] = 0;
      const double lj = l[j];
      l[j] = 0.0;
      if (p.col_mf.empty() || p.col_mf[j] < 0) {
        soma += lj * lj * d[j];
      } else {
        lmf[static_cast<std::size_t>(p.col_mf[j])] += lj;
        continue;   // metafundador nao tem pais; o termo dele entra na forma quadratica
      }
      const std::int64_t sj = p.pai[j], tj = p.mae[j];
      if (sj >= 0) {
        l[sj] += 0.5 * lj;
        if (!na_fila[sj]) { fila.push(static_cast<std::size_t>(sj)); na_fila[sj] = 1; }
      }
      if (tj >= 0) {
        // o caminho do avo materno pesa 1/4: u_j = u_s/2 + u_k/4 + m_j
        l[tj] += ((!p.mgs.empty() && p.mgs[j]) ? 0.25 : 0.5) * lj;
        if (!na_fila[tj]) { fila.push(static_cast<std::size_t>(tj)); na_fila[tj] = 1; }
      }
    }
    // ||K' l_mf||^2 = l_mf' K K' l_mf = l_mf' Gamma l_mf. Com Gamma diagonal K e
    // diag(sqrt(gamma)) e isto colapsa em soma_j l_j^2 gamma_jj, o comportamento anterior.
    for (std::size_t c = 0; c < p.n_mf; c++) {
      double s = 0.0;
      // a soma vai de ZERO: K vem de uma decomposicao espectral e NAO e triangular, entao
      // comecar em c descartaria silenciosamente as entradas acima da diagonal. O que K
      // precisa satisfazer e K K' = Gamma, e a forma quadratica nao distingue as duas.
      for (std::size_t a = 0; a < p.n_mf; a++) s += lmf[a] * p.gama_chol[a * p.n_mf + c];
      soma += s * s;
    }
    f[i] = soma - 1.0;
    // O grampo em zero absorve ARREDONDAMENTO, e so isso. Num pedigree classico um F
    // negativo e sempre ruido, porque nao ha como dois pais serem menos que nao
    // aparentados. Com Gamma cheia ha: um gamma_jk NEGATIVO (bases divergidas por selecao
    // em sentidos opostos, que o artigo permite explicitamente) da a_st < 0 e portanto
    // F = 0,5 a_st < 0 de verdade. Medido: com Gamma = [[0,5; -0,2],[-0,2; 0,9]] o filho
    // dos dois metafundadores tem F = -0,10 exato, e grampear isso em zero errava tambem
    // os tres descendentes dele. Entao clampa-se so a escala do arredondamento.
    if (f[i] < 0.0 && f[i] > -1e-10) f[i] = 0.0;
    d[i] = mendel(i);
  }
  return f;
}

// A^-1 de Henderson, no triangulo INFERIOR.
//
//   A^-1 = sum_i d_i^-1 v_i v_i',  v = (1, -1/2, -1/2) sobre (animal, pai, mae)
//
// O produto externo e montado a partir dos coeficientes JA FUNDIDOS por indice, e nao
// termo a termo. A razao e um caso que quebrou a primeira versao deste arquivo: quando pai
// e mae sao o MESMO animal, o vetor colapsa de (1, -1/2, -1/2) para (1, -1), e a
// contribuicao na diagonal do progenitor passa a ser 1/d e nao 3/4 de 1/d.
//
// Guardando as duas metades da matriz, somar (s,d) e (d,s) separadamente acerta esse caso
// por acidente. Guardando so o triangulo inferior, as duas caem na mesma posicao e uma
// delas se perde. Fundir os coeficientes antes acerta os dois casos sem ramo nenhum.
//
// Autofecundacao nao acontece em suino, mas acontece em planta e, sobretudo, acontece como
// ERRO DE DADO. Uma A^-1 errada por causa disso continua simetrica e positiva-definida.
Csc a_inversa(const Pedigree& p, const std::vector<double>& f) {
  const std::size_t n = p.ids.size();
  std::vector<std::uint32_t> li, cj;
  std::vector<double> v;
  li.reserve(n * 6); cj.reserve(n * 6); v.reserve(n * 6);

  // (indice, coeficiente) do animal e dos progenitores, com repetidos somados
  std::vector<std::pair<std::size_t, double>> vi;
  vi.reserve(3);

  // O bloco Gamma^-1 entra SOMADO, junto com tudo o mais. As linhas dos metafundadores
  // tambem recebem as contribuicoes de Henderson vindas dos filhos deles, entao escrever o
  // bloco por cima no fim apagaria essas contribuicoes. Emitir aqui, como mais um conjunto
  // de tripletos, deixa a soma por conta do montador de CSC, que ja soma repetidos.
  //
  // Com Gamma diagonal, Gamma^-1 = diag(1/gamma_ii) e isto reproduz exatamente a entrada
  // 1/gamma_ii que o laco por animal emitia antes.
  if (p.n_mf) {
    std::vector<std::size_t> linha_do_mf(p.n_mf, 0);
    for (std::size_t i = 0; i < n; i++)
      if (!p.col_mf.empty() && p.col_mf[i] >= 0)
        linha_do_mf[static_cast<std::size_t>(p.col_mf[i])] = i;
    for (std::size_t a = 0; a < p.n_mf; a++)
      for (std::size_t b = 0; b <= a; b++) {
        std::size_t r = linha_do_mf[a], c = linha_do_mf[b];
        double val = p.gama_inv[a * p.n_mf + b];
        if (r < c) std::swap(r, c);
        li.push_back(static_cast<std::uint32_t>(r));
        cj.push_back(static_cast<std::uint32_t>(c));
        v.push_back(val);
      }
  }

  for (std::size_t i = 0; i < n; i++) {
    // o metafundador ja entrou pelo bloco Gamma^-1 acima; emitir 1/gamma_ii aqui de novo
    // duplicaria a diagonal dele
    if (!p.eh_mf.empty() && p.eh_mf[i]) continue;
    const std::int64_t s = p.pai[i], t = p.mae[i];
    const double di = variancia_mendeliana(p, f, i);
    if (!(di > 0.0))
      throw Erro("non-positive Mendelian variance at '" + p.ids[i] +
                 "'; the pedigree has impossible inbreeding");
    const double a = 1.0 / di;

    vi.clear();
    auto junta = [&](std::size_t idx, double c) {
      for (auto& e : vi) if (e.first == idx) { e.second += c; return; }
      vi.emplace_back(idx, c);
    };
    junta(i, 1.0);
    if (s >= 0) junta(static_cast<std::size_t>(s), -0.5);
    // (1, -1/2, -1/4) no pedigree pai/avo materno; a fusao por indice acima ja cobre o
    // touro acasalado com a propria filha (pai == avo materno, coeficiente -3/4)
    if (t >= 0) junta(static_cast<std::size_t>(t), (!p.mgs.empty() && p.mgs[i]) ? -0.25 : -0.5);

    for (std::size_t x = 0; x < vi.size(); x++)
      for (std::size_t y = 0; y <= x; y++) {
        std::size_t r = vi[x].first, c = vi[y].first;
        double val = a * vi[x].second * vi[y].second;
        if (r < c) std::swap(r, c);
        li.push_back(static_cast<std::uint32_t>(r));
        cj.push_back(static_cast<std::uint32_t>(c));
        v.push_back(val);
      }
  }
  return de_triplos(n, n, li, cj, v);
}

}  // namespace br
