// Multi-caracteristica, parte 2: layout, desenho e o estimador.

#include "mme.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>

// so para o TRACE opcional e a interrupcao: nenhuma conta usa o R
#include <R_ext/Print.h>
#include <R_ext/Utils.h>

namespace br {

// Reindexa o layout do modelo para t caracteristicas: o termo q de um grupo passa a ocupar
// n_coef * nt(q) coordenadas, com nt(q) os caracteres dele (todos, salvo um traits =), e o
// residual passa a vech(R0) com t(t+1)/2 parametros. Os termos de um grupo podem ter
// caracteres diferentes: rn(id, base, traits = "y1", group = "g") + animal(id, traits = "y2",
// group = "g") e a G 3 x 3 de (a1_0, a1_1, a2) de um traco longitudinal com um de registro
// unico. O bloco continua kron(C_g, K): cada coordenada so precisa dos niveis comuns.
static void expande_layout_mt(Modelo& m, const std::vector<std::string>& alvos) {
  const std::size_t t = alvos.size();
  std::size_t off = 0;
  for (Grupo& g : m.grupos) {
    g.tracos.assign(g.termos.size(), {});
    g.dim = 0;
    for (std::size_t q = 0; q < g.termos.size(); q++) {
      const Termo& te = m.termos[g.termos[q]];
      std::vector<char> tem(t, te.tracos.empty() ? 1 : 0);
      for (const std::string& nm : te.tracos) {
        const auto it = std::find(alvos.begin(), alvos.end(), nm);
        if (it == alvos.end())
          throw Erro("term '" + te.nome + "': traits = names '" + nm +
                     "', which is not a trait of cbind()");
        tem[static_cast<std::size_t>(it - alvos.begin())] = 1;
      }
      std::vector<std::size_t> este;
      for (std::size_t tau = 0; tau < t; tau++) if (tem[tau]) este.push_back(tau);
      if (este.size() < t) g.tracos[q] = este;   // todos: a lista fica vazia
      g.dim += te.n_coef() * g.nt(q, t);
    }
    g.offset = off;
    g.nparam = g.dim * (g.dim + 1) / 2;
    off += g.nparam;
  }
  m.offset_residual = off;
  m.ntheta = off + t * (t + 1) / 2;
}

// Nomes de theta no caso multi: coeficiente rotulado "termo[c]@trait".
//
// A ORDEM DOS ROTULOS TEM DE SER A ORDEM DO LAYOUT, e nao era. A dimensao do grupo e
// numerada TERMO por fora, caracteristica, coeficiente por dentro: e como o Z e montado
// (col0 + (tau * n_coef + ct) * n_niveis por termo) e como os EBV sao lidos. Os rotulos
// andavam caracteristica por fora e termo por dentro. Com UM termo por grupo as duas ordens
// coincidem, que e o caso de quase todo ajuste, e por isso nada acusava. Com dois termos no
// mesmo grupo -- direto e indireto, direto e materno -- num ajuste de dois caracteres os
// componentes saiam com o NOME ERRADO: o que se chamava var(indirect@y) era var(animal@y2),
// e rg(), accuracy() e start= por nome liam o numero de outro componente.
//
// Medido pela identidade que tem de fechar exatamente: com as covariancias entre
// caracteres em zero, o -2logL bivariado e a soma dos dois univariados. Montando o theta
// pela ordem dos rotulos antigos a diferenca era 43.9; pela ordem do layout, 1.7e-10.
static std::vector<std::string> nomes_theta_mt(const Modelo& m,
                                               const std::vector<std::string>& alvos) {
  const std::size_t t = alvos.size();
  std::vector<std::string> out(m.ntheta);
  for (const Grupo& g : m.grupos) {
    std::vector<std::string> rot;
    for (std::size_t q = 0; q < g.termos.size(); q++)
      for (std::size_t pt = 0; pt < g.nt(q, t); pt++) {
        const std::size_t tau = g.traco(q, pt);
        const Termo& te = m.termos[g.termos[q]];
        for (std::size_t c = 0; c < te.n_coef(); c++) {
          std::string b = te.nome;
          if (te.n_coef() > 1) b += "[" + std::to_string(c) + "]";
          rot.push_back(b + "@" + alvos[tau]);
        }
      }
    for (std::size_t j = 0; j < g.dim; j++)
      for (std::size_t i = j; i < g.dim; i++) {
        const std::size_t k = g.offset + (j * (2 * g.dim - j + 1)) / 2 + (i - j);
        out[k] = (i == j) ? "var(" + rot[i] + ")" : "cov(" + rot[i] + "," + rot[j] + ")";
      }
  }
  for (std::size_t j = 0; j < t; j++)
    for (std::size_t i = j; i < t; i++) {
      const std::size_t k = m.offset_residual + (j * (2 * t - j + 1)) / 2 + (i - j);
      out[k] = (i == j) ? "var(res@" + alvos[i] + ")"
                        : "cov(res@" + alvos[i] + ",res@" + alvos[j] + ")";
    }
  return out;
}

DesenhoMT monta_desenho_mt(Modelo m, const std::vector<std::string>& alvos,
                           const Tabela& tab, const Pedigree* ped,
                            const std::vector<KernelDecl>* kernels) {
  DesenhoMT d;
  d.alvos = alvos;
  d.t = alvos.size();
  if (d.t < 2) throw Erro("multi-trait needs at least two; for one, use model()");
  d.nlin = tab.nlin;

  // y n x t. Caracteristica ausente NAO derruba o registro: marca obs = 0, e o registro
  // participa com a submatriz de R0 do seu padrao. y ausente vira zero, com a inversa
  // embutida zerada naquela linha, o valor nunca e lido, e zero e mais seguro que NaN,
  // que contamina qualquer soma que um defeito futuro deixar passar.
  d.y = Densa(d.nlin, d.t);
  d.usa.assign(d.nlin, 1);
  d.obs.assign(d.nlin * d.t, 1);
  for (std::size_t tau = 0; tau < d.t; tau++) {
    std::vector<double> col = tab.numerico(alvos[tau]);
    for (std::size_t i = 0; i < d.nlin; i++) {
      const bool falta = !std::isfinite(col[i]) ||
          (m.tem_ausente && std::fabs(col[i] - m.codigo_ausente) < 1e-9);
      d.obs[i * d.t + tau] = falta ? 0 : 1;
      d.y.at(i, tau) = falta ? 0.0 : col[i];
    }
  }
  // registro sem NENHUMA caracteristica observada sai
  for (std::size_t i = 0; i < d.nlin; i++) {
    bool alguma = false;
    for (std::size_t tau = 0; tau < d.t; tau++)
      if (d.obs[i * d.t + tau]) { alguma = true; break; }
    if (!alguma) d.usa[i] = 0;
  }

  const bool precisa_ped = [&]{
    for (const Termo& tm : m.termos)
      if (tm.estrutura == Estrutura::Parentesco) return true;
    return false;
  }();
  if (precisa_ped && !ped) throw Erro("there is a term with relationship and no pedigree was given");

  // a K declarada passa pela MESMA reducao do univariado antes de qualquer conta
  std::vector<KernelDecl> kern_red;
  std::vector<std::unordered_set<std::string> > kern_nulos;
  reduz_kernels(m, kernels, kern_red, kern_nulos);

  Csc ainv;
  double ld_ainv = 0.0;
  std::vector<std::string> niveis_ped;
  if (precisa_ped) {
    std::vector<double> f = endogamia(*ped);
    ainv = a_inversa(*ped, f);
    niveis_ped = ped->ids;
    std::vector<std::size_t> perm = grau_minimo(ainv);
    Csc pa = permuta_sim(ainv, perm);
    Simbolica sb = simbolica(pa);
    Csc L;
    if (!cholesky(pa, sb, L)) throw Erro("the pedigree A^-1 is not positive-definite");
    ld_ainv = logdet(L);
  }
  for (const Grupo& g : m.grupos) {
    if (g.estrutura == Estrutura::Parentesco) {
      d.kinv.push_back(ainv);
      d.kinv_logdet.push_back(ld_ainv);
      d.kinv_escala.push_back(0.0);
    } else if (g.estrutura == Estrutura::Declarada) {
      kinv_declarada(m, g, kernels, d.kinv, d.kinv_logdet, d.kinv_escala);
    } else {
      d.kinv.push_back(Csc());
      d.kinv_logdet.push_back(0.0);
      d.kinv_escala.push_back(0.0);
    }
  }

  // X por registro, esparso (fixos.cpp); a expansao por caracteristica e na montagem
  FixosCandidatos fx = colunas_fixas(m, tab);
  // AS LINHAS SEM PAR CAEM ANTES DO POSTO DE X, como em model() (mme.cpp). Antes caiam
  // depois: um CG cujos animais nao estao no pedigree ficava com coluna em X, o posto nao
  // a via como dependente, e a coluna chegava vazia na montagem. model() derrubava o nivel
  // e ajustava; model_ar1() devolvia "did NOT converge" com decremento 0,46 e model_mt()
  // "SINGULAR", dois diagnosticos que apontavam para o lugar errado.
  // niveis: do pedigree quando ha parentesco, DA K quando declarada, da uniao das colunas
  // num grupo sem estrutura de varios termos. Todo nivel do conjunto ganha equacao, com ou
  // sem registro, exatamente como em model() (monta_aleatorios, em mme.cpp).
  d.aleatorios = monta_aleatorios(m, tab, niveis_ped, kernels);
  {
    std::vector<DesenhoTermo*> pa;
    for (DesenhoTermo& a : d.aleatorios) pa.push_back(&a);
    casa_niveis_nulos(m, pa, kern_nulos, tab, d.nlin);
  }
  {
    std::vector<char> sem_obs(d.nlin, 0);
    for (std::size_t i = 0; i < d.nlin; i++) sem_obs[i] = d.usa[i] ? 0 : 1;
    for (const DesenhoTermo& a : d.aleatorios)
      for (std::size_t i = 0; i < d.nlin; i++)
        if (!a.casou[i]) d.usa[i] = 0;
    // antes de qualquer posto: sem linha, o posto ficava sem coluna e a recusa saia como
    // "no estimable fixed effect for any trait", que aponta para o lugar errado
    if (d.n_usadas() == 0)
      throw Erro(sem_linhas(m, tab, d.aleatorios, sem_obs,
                            "every trait missing (NA, NaN, Inf or the missing code)"));
  }
  // POSTO DE X SOBRE AS LINHAS QUE ENTRAM, nao sobre a tabela inteira. E o mesmo conserto
  // que o univariado ja tem (mme.cpp, mesma marca). Um nivel fixo cujos registros TODOS
  // sairam continua com coluna nao-nula na tabela, entao um posto medido sobre tudo o
  // mantinha; na montagem, que so anda nas linhas usadas, essa coluna nao recebe nada, a
  // matriz de coeficientes ganha coluna VAZIA, e a Cholesky morre com "nao
  // positiva-definida" tres passos longe da causa. Medindo o posto onde o modelo de fato
  // vive, a coluna sai como dependente e e REPORTADA em dropped_x.
  // A REFERENCIA GLOBAL ENTRE OS NIVEIS COM REGISTRO DE CADA CARACTER (P1-8). O posto e
  // Gram-Schmidt na ordem das colunas e derruba a ULTIMA de cada dependencia: com a ordem
  // de sempre a referencia de um CG era o ultimo nivel, mesmo quando ele nao tinha registro
  // de um dos caracteres. Ai a passada por caracter tinha de derrubar outro nivel daquele
  // caracter, um nivel COM registro, e dropped_x dizia 'CG=x' (sem sufixo, vazio no
  // caracter) e 'CG=y|traco' (nao vazio), sem distinguir vazio de referencia (medido em
  // dado real; posto, verossimilhanca e EBV nao mudam, so o rotulo e a parametrizacao de b).
  // As colunas com registro em TODOS os caracteres vao para o fim da ordem do posto (o
  // intercepto fica na frente), entao a referencia que cai e um nivel completo sempre que
  // existe um; a passada por caracter derruba entao so o que de fato nao tem registro, e
  // cada queda leva o motivo em saiu_motivo. Sem nenhum nivel incompleto a ordem e a de
  // sempre, e nada muda.
  // O posto e o esparso, na ordem acima (posto_referencia_completa, fixos.cpp).
  // AS COVARIAVEIS PADRONIZADAS ANTES DO POSTO (padroniza_covariaveis, fixos.cpp), e a volta de
  // b por caracter. Uma aninhada cujo centro nao volta pela classe em algum caracter (a classe
  // confundida com outro efeito fixo, aninhadas_sem_volta) e refeita so escalada, com os dois
  // postos de novo.
  // UM EFEITO FIXO SO EM ALGUNS CARACTERES (fixed(cg, traits = "y1"), cov(x, traits = "y1")):
  // por termo, os caracteres em que ele entra. Nos outros as colunas dele saem do posto daquele
  // caracter antes da conta e nao vao para dropped_x, porque nao faltou registro nem houve
  // dependencia: o modelo declarado nao as tem. E a coluna 0 do efeito no caracter, no BLUPF90.
  std::vector<std::vector<char>> fixo_em(m.termos.size(), std::vector<char>(d.t, 1));
  for (std::size_t k = 0; k < m.termos.size(); k++) {
    const Termo& te = m.termos[k];
    if (te.aleatorio() || te.tracos.empty()) continue;
    std::fill(fixo_em[k].begin(), fixo_em[k].end(), 0);
    for (const std::string& nm : te.tracos) {
      const auto it = std::find(alvos.begin(), alvos.end(), nm);
      if (it == alvos.end())
        throw Erro("term '" + te.nome + "': traits = names '" + nm +
                   "', which is not a trait of cbind()");
      fixo_em[k][static_cast<std::size_t>(it - alvos.begin())] = 1;
    }
  }
  std::vector<std::size_t> fica, sai, sem_centro;
  for (;;) {
    padroniza_covariaveis(fx, d.usa, sem_centro);
    d.nomes_x.clear();
    d.saiu_x.clear();
    {
      std::vector<std::string> motivo;
      posto_referencia_completa(fx, d.usa, d.obs, d.t, 1e-9, fica, sai, motivo);
      d.saiu_motivo = std::move(motivo);
    }
    d.x = fixas_por_linha(fx.x, fica);
    for (std::size_t j : fica) d.nomes_x.push_back(fx.nomes[j]);
    for (std::size_t j : sai) d.saiu_x.push_back(fx.nomes[j]);

    // POSTO DE X POR CARACTER, e nao so ausencia de registro por caracter.
    //
    // O posto acima e medido sobre as linhas USADAS, e uma linha conta como usada quando
    // QUALQUER caracter foi observado nela. Isso deixa passar o caso que aparece em dado de
    // verdade quando dois caracteres foram medidos em coortes diferentes. Medido no desenho que
    // motivou isto, CG aninhado no caracter: X tem posto 8 de 8 sobre as linhas usadas e 4 de 8
    // DENTRO de cada caracter. As equacoes sao montadas por caracter, entao a deficiencia
    // aparece la, a matriz de coeficientes fica singular e a Cholesky morre tres passos longe
    // da causa.
    //
    // O univariado ja derruba coluna inestimavel e a reporta em dropped_x. Aqui a coluna pode
    // ser estimavel para um caracter e nao para outro, entao a unidade que cai e o PAR
    // (coluna, caracter).
    //
    // NAO BASTA testar se o par tem registro: depois de tirar os pares vazios, o que sobra para
    // um caracter continua podendo ser COLINEAR. No desenho acima y1 fica com intercepto mais
    // CG1..CG4, que somam ao intercepto: posto 4 e nao 5. Por isso aqui roda o mesmo
    // posto_esparso() de cima, uma vez por caracter, sobre as linhas em que aquele caracter foi
    // observado. Ausencia de registro e o caso particular de coluna nula, que o posto ja pega.
    //
    // A numeracao final segue a ordem antiga, j externo e tau interno, entao quando nada cai
    // eq_fixa[k] == k e n_fixa == x.ncol * t: o layout antigo e caso particular deste, e nenhum
    // ajuste que ja andava muda de numero.
    d.eq_fixa.assign(d.nomes_x.size() * d.t, -1);
    {
      const Csc xc = fixas_por_coluna(d.x);
      // a coluna j de d.x e fora do caracter tau quando o termo dela nao o tem (traits =)
      auto fora = [&](std::size_t j, std::size_t tau) {
        const std::size_t k = fx.termo[fica[j]];
        return k < fixo_em.size() && !fixo_em[k][tau];
      };
      for (std::size_t tau = 0; tau < d.t; tau++) {
        std::vector<char> lin_tau(d.nlin, 0);
        for (std::size_t i = 0; i < d.nlin; i++)
          lin_tau[i] = d.usa[i] && d.obs[i * d.t + tau];
        // as colunas fora do caracter entram zeradas: o posto as poe em sai_t, e elas ficam
        // sem equacao e fora do relato
        Csc xt = xc;
        bool tira = false;
        for (std::size_t j = 0; j < xt.ncol; j++)
          if (fora(j, tau)) {
            tira = true;
            for (std::size_t q = xt.colptr[j]; q < xt.colptr[j + 1]; q++) xt.valor[q] = 0.0;
          }
        std::vector<std::size_t> fica_t, sai_t;
        posto_esparso(tira ? xt : xc, lin_tau, 1e-9, fica_t, sai_t, false, nullptr, &d.nomes_x);
        for (std::size_t j : fica_t) d.eq_fixa[j * d.t + tau] = -2;
        for (std::size_t j : sai_t) {
          if (fora(j, tau)) continue;
          d.saiu_x.push_back(d.nomes_x[j] + "|" + d.alvos[tau]);
          bool vazia = true;
          for (std::size_t k = xc.colptr[j]; k < xc.colptr[j + 1]; k++)
            if (lin_tau[xc.linha[k]] && xc.valor[k] != 0.0) { vazia = false; break; }
          d.saiu_motivo.push_back((vazia ? "no record of " : "aliased within ") + d.alvos[tau]);
        }
      }
      int prox = 0;
      for (std::size_t j = 0; j < d.x.ncol; j++)
        for (std::size_t tau = 0; tau < d.t; tau++)
          if (d.eq_fixa[j * d.t + tau] == -2) d.eq_fixa[j * d.t + tau] = prox++;
      d.n_fixa = static_cast<std::size_t>(prox);
    }

    std::vector<std::size_t> falha;
    d.volta_x.clear();
    for (std::size_t tau = 0; tau < d.t; tau++) {
      std::vector<char> lin_tau(d.nlin, 0), mantida(fx.x.ncol, 0);
      for (std::size_t i = 0; i < d.nlin; i++) lin_tau[i] = d.usa[i] && d.obs[i * d.t + tau];
      for (std::size_t jj = 0; jj < fica.size(); jj++)
        mantida[fica[jj]] = d.eq_fixa[jj * d.t + tau] >= 0;
      for (std::size_t k : aninhadas_sem_volta(fx, mantida, lin_tau))
        if (std::find(falha.begin(), falha.end(), k) == falha.end()) falha.push_back(k);
      d.volta_x.push_back(monta_volta(fx, fica, mantida, lin_tau));
    }
    if (falha.empty()) break;
    sem_centro.insert(sem_centro.end(), falha.begin(), falha.end());
    fx = colunas_fixas(m, tab);
  }
  if (d.n_fixa == 0) throw Erro("no estimable fixed effect for any trait");
  // o termo do -2logL da padronizacao: cada par (coluna, caracter) com equacao conta uma vez
  for (std::size_t j = 0; j < d.x.ncol; j++)
    for (std::size_t tau = 0; tau < d.t; tau++)
      if (d.eq_fixa[j * d.t + tau] >= 0 && d.volta_x[tau].escala[j] != 1.0)
        d.logdet_escala_x += 2.0 * std::log(d.volta_x[tau].escala[j]);

  if (d.n_usadas() == 0) throw Erro("no row enters the analysis");

  expande_layout_mt(m, alvos);
  for (const Grupo& g : m.grupos)
    for (std::size_t q = 0; q < g.termos.size(); q++)
      for (DesenhoTermo& a : d.aleatorios)
        if (a.termo == g.termos[q]) a.tracos = g.tracos[q];
  d.modelo = std::move(m);
  return d;
}

std::vector<std::string> nomes_theta_do_mt(const DesenhoMT& d) {
  return nomes_theta_mt(d.modelo, d.alvos);
}

static std::vector<double> partida_mt(const DesenhoMT& d) {
  const std::size_t t = d.t;
  std::vector<double> var_t(t, 1.0);
  for (std::size_t tau = 0; tau < t; tau++) {
    double soma = 0.0, soma2 = 0.0;
    std::size_t n = 0;
    for (std::size_t i = 0; i < d.nlin; i++)
      if (d.usa[i] && d.obs[i * t + tau]) {
        soma += d.y.at(i, tau);
        soma2 += d.y.at(i, tau) * d.y.at(i, tau);
        n++;
      }
    const double media = soma / std::max<std::size_t>(n, 1);
    double v = (soma2 - n * media * media) / std::max(1.0, static_cast<double>(n - 1));
    if (!(v > 0.0) || !std::isfinite(v)) v = 1.0;
    var_t[tau] = v;
  }
  std::vector<double> theta(d.modelo.ntheta, 0.0);
  const double frac = 0.5 / static_cast<double>(d.modelo.grupos.size());
  std::size_t ig = 0;
  for (const Grupo& g : d.modelo.grupos) {
    // A escala da K entra na PARTIDA. Uma K declarada nao tem escala canonica: o usuario
    // pode passar D ou 4D para o mesmo modelo, e a verossimilhanca e identica desde que o
    // componente ande junto. Partir de frac*var(y) nos dois casos e partir de dois pontos
    // DIFERENTES da mesma superficie, e a diferenca sobrevive ate o fim (medido: 204.18
    // contra 215.64 de -2logL no mesmo ajuste, com converged TRUE nos dois). A media
    // geometrica dos autovalores de K sai de graca do log-determinante que o desenho ja
    // guarda, |K^-1|, e dividir por ela deixa a partida EQUIVARIANTE: com cK ela cai por c
    // e o caminho e o mesmo. So a K declarada precisa disso; A tem diagonal ~1 e a
    // diagonal e I, e mexer nelas so deslocaria ajustes que ja estao certos.
    double esc = 1.0;
    // Numa precisao aumentada (dominance_inverse()) a escala vem pronta no KernelDecl: a
    // media geometrica dos autovalores de D, nao a de K_aug com os niveis de subclasse, para
    // que kernel(Kinv =) parta do MESMO ponto de kernel(K = D).
    if (g.estrutura == Estrutura::Declarada && ig < d.kinv.size() && d.kinv[ig].ncol > 0) {
      const double gm = escala_kernel(d.kinv_escala, d.kinv_logdet, d.kinv, ig);
      if (std::isfinite(gm) && gm > 0.0) esc = gm;
    }
    ig++;
    // diagonal por caracteristica, na ordem do layout (termo por fora, caracteristica,
    // coeficiente por dentro; ver nomes_theta_mt). Lia "os primeiros dim/t slots sao da
    // caracteristica 0", que so vale com um termo no grupo: com direto e indireto juntos a
    // var(y) de y2 caia na var(animal@y2) trocada com a de y1 na var(indirect@y)
    const std::vector<std::size_t> tr = traco_das_coordenadas(d.modelo, g, d.t);
    for (std::size_t i = 0; i < g.dim; i++) {
      const std::size_t k = g.offset + (i * (2 * g.dim - i + 1)) / 2;   // vech (i,i)
      theta[k] = frac * var_t[i < tr.size() ? tr[i] : 0] / esc;
    }
  }
  for (std::size_t tau = 0; tau < d.t; tau++) {
    const std::size_t k = d.modelo.offset_residual + (tau * (2 * d.t - tau + 1)) / 2;
    theta[k] = 0.5 * var_t[tau];
  }
  return theta;
}

// O fator de cada componente sob a padronizacao por caracteristica: uma covariancia entre
// as coordenadas a e b do bloco escala por s_a s_b, com s o desvio padrao da caracteristica
// de cada coordenada. Nos grupos a coordenada e (termo, caracteristica, coeficiente), na
// ordem do layout (termo por fora, caracteristica, coeficiente por dentro); no residuo e a
// propria caracteristica. O rho do AR(1) fica com fator 1 (nao tem unidade).
std::vector<std::size_t> traco_das_coordenadas(const Modelo& m, const Grupo& g,
                                               std::size_t t) {
  std::vector<std::size_t> tr;
  for (std::size_t q = 0; q < g.termos.size(); q++) {
    const std::size_t nc = m.termos[g.termos[q]].n_coef();
    for (std::size_t pt = 0; pt < g.nt(q, t); pt++)
      for (std::size_t c = 0; c < nc; c++) tr.push_back(g.traco(q, pt));
  }
  if (tr.size() != g.dim) tr.clear();
  return tr;
}

std::vector<double> fatores_escala(const Modelo& m, std::size_t t, std::size_t off_residuo,
                                   const std::vector<double>& esc) {
  std::vector<double> f(m.ntheta, 1.0);
  for (const Grupo& g : m.grupos) {
    const std::vector<std::size_t> tr = traco_das_coordenadas(m, g, t);
    if (tr.empty()) continue;   // layout desconhecido: fica sem escala
    for (std::size_t j = 0; j < g.dim; j++)
      for (std::size_t i = j; i < g.dim; i++)
        f[g.offset + (j * (2 * g.dim - j + 1)) / 2 + (i - j)] = esc[tr[i]] * esc[tr[j]];
  }
  for (std::size_t j = 0; j < t; j++)
    for (std::size_t i = j; i < t; i++)
      f[off_residuo + (j * (2 * t - j + 1)) / 2 + (i - j)] = esc[i] * esc[j];
  return f;
}

// A unidade da padronizacao de cada caracteristica: a POTENCIA DE 2 logo acima do desvio
// padrao nas celulas observadas dos registros usados (o desvio padronizado fica em [0,5; 1)),
// 1 quando nao ha variacao. Potencia de 2 e nao o proprio desvio porque dividir por ela e
// EXATO em ponto flutuante: y padronizado carrega bit a bit a informacao de y, a volta das
// componentes (vezes s_a s_b) tambem e exata, e y e 2^k y dao o MESMO ajuste bit a bit. Com
// outro fator (100 y) a unica diferenca que sobra e o arredondamento do proprio dado.
std::vector<double> desvios_tracos(const Densa& y, const std::vector<char>& usa,
                                   const std::vector<char>& obs, std::size_t t,
                                   std::vector<double>* exatos) {
  std::vector<double> esc(t, 1.0);
  if (exatos) exatos->assign(t, 1.0);
  for (std::size_t tau = 0; tau < t; tau++) {
    double soma = 0.0;
    std::size_t n = 0;
    for (std::size_t i = 0; i < y.nlin; i++)
      if (usa[i] && obs[i * t + tau]) { soma += y.at(i, tau); n++; }
    if (n < 2) continue;
    const double media = soma / static_cast<double>(n);
    double ss = 0.0;
    for (std::size_t i = 0; i < y.nlin; i++)
      if (usa[i] && obs[i * t + tau]) ss += (y.at(i, tau) - media) * (y.at(i, tau) - media);
    const double v = ss / static_cast<double>(n - 1);
    if (v > 0.0 && std::isfinite(v)) {
      int e = 0;
      std::frexp(std::sqrt(v), &e);   // sd = m 2^e, m em [0,5; 1)
      esc[tau] = std::ldexp(1.0, e);
      if (exatos) (*exatos)[tau] = std::sqrt(v);
    }
  }
  return esc;
}

bool grupos_de_univariados(const Modelo& m, std::size_t t,
                           const std::vector<std::vector<double> >& uni,
                           std::vector<double>& theta) {
  if (uni.size() != t) return false;
  // onde o grupo comeca no univariado de cada caracter: o univariado de tau so tem os termos
  // que existem em tau (traits =), e um grupo sem nenhum deles nao esta la
  std::vector<std::size_t> off_u(t, 0);
  for (const Grupo& g : m.grupos) {
    // primeira coordenada de cada termo no grupo expandido
    std::vector<std::size_t> base, nc;
    std::size_t acc = 0;
    for (std::size_t q = 0; q < g.termos.size(); q++) {
      base.push_back(acc);
      nc.push_back(m.termos[g.termos[q]].n_coef());
      acc += nc.back() * g.nt(q, t);
    }
    if (acc != g.dim) return false;
    for (std::size_t tau = 0; tau < t; tau++) {
      // as coordenadas do grupo no univariado de tau, na ordem dele: (termo, coeficiente) dos
      // termos que existem em tau, com a posicao pt de tau no termo
      std::vector<std::size_t> mt;
      for (std::size_t q = 0; q < g.termos.size(); q++)
        for (std::size_t pt = 0; pt < g.nt(q, t); pt++)
          if (g.traco(q, pt) == tau)
            for (std::size_t c = 0; c < nc[q]; c++) mt.push_back(base[q] + pt * nc[q] + c);
      const std::size_t dim_u = mt.size();
      if (dim_u == 0) continue;
      if (uni[tau].size() < off_u[tau] + dim_u * (dim_u + 1) / 2) return false;
      for (std::size_t ju = 0; ju < dim_u; ju++)
        for (std::size_t iu = ju; iu < dim_u; iu++) {
          const double v = uni[tau][off_u[tau] + (ju * (2 * dim_u - ju + 1)) / 2 + (iu - ju)];
          if (!std::isfinite(v)) return false;
          std::size_t i = mt[iu], j = mt[ju];
          if (i < j) std::swap(i, j);
          theta[g.offset + (j * (2 * g.dim - j + 1)) / 2 + (i - j)] = v;
        }
      off_u[tau] += dim_u * (dim_u + 1) / 2;
    }
  }
  return true;
}

// Os otimos univariados viram partida com as covariancias entre caracteristicas em zero. Uma
// variancia que o univariado deixou NA FRONTEIRA DO ZERO (abaixo de 1e-4 da partida
// automatica, var(y) repartida entre os grupos: o piso absoluto do univariado, 1e-10 s2e, e
// o relativo, 1e-6 da maior variancia do grupo, ficam os dois muito abaixo disso) sobe para
// 1% da partida automatica: senao a coordenada nasceria encostada no piso do passo, de onde
// o caminhante sai so pela reentrada. Uma variancia interior pequena fica como o univariado
// a deixou: o piso antigo, 1% em toda variancia, levantava tambem otimos interiores (medido
// no caso_ige de test-accuracy-kernel.R com y2 = 30 N(0,1) + 0,3 y: var(indirect@y2) 0,117
// no univariado partia de 4,087, 35 vezes). O piso so SOBE a diagonal, entao cada bloco
// segue positivo-semidefinido.
std::vector<double> partida_univariados_mt(const DesenhoMT& d,
                                           const std::vector<std::vector<double> >& uni) {
  const std::size_t t = d.t;
  // o tamanho do vetor univariado de cada caracter: so os grupos que o contem, e o residuo
  std::vector<std::size_t> n_u(t, 1);
  for (const Grupo& g : d.modelo.grupos)
    for (std::size_t tau = 0; tau < t; tau++) {
      std::size_t du = 0;   // as coordenadas do grupo no univariado de tau
      for (std::size_t q = 0; q < g.termos.size(); q++)
        for (std::size_t pt = 0; pt < g.nt(q, t); pt++)
          if (g.traco(q, pt) == tau) du += d.modelo.termos[g.termos[q]].n_coef();
      n_u[tau] += du * (du + 1) / 2;
    }
  if (uni.size() != t) return {};
  for (std::size_t tau = 0; tau < t; tau++)
    if (uni[tau].size() != n_u[tau]) return {};
  std::vector<double> theta(d.modelo.ntheta, 0.0);
  if (!grupos_de_univariados(d.modelo, t, uni, theta)) return {};
  for (std::size_t tau = 0; tau < t; tau++) {
    const double v = uni[tau][n_u[tau] - 1];
    if (!std::isfinite(v)) return {};
    theta[d.modelo.offset_residual + (tau * (2 * t - tau + 1)) / 2] = v;
  }
  const std::vector<double> auto_ = partida_mt(d);
  for (std::size_t k = 0; k < theta.size(); k++)
    if (auto_[k] > 0.0 && theta[k] < 1e-4 * auto_[k]) theta[k] = 1e-2 * auto_[k];
  return theta;
}

// PADRONIZACAO INTERNA POR CARACTERISTICA (P1-8, e o D1 do P1-5 do lado do multicaracter).
// Medido em dado real: com um caracter reescalado por 0,01 a estimativa de fronteira mudava
// 0,66% e o -2logL 0,011; em simulado, com a razao de variancias 3600 o veredito SINGULAR
// aparecia, e com 1000 vezes na unidade o ajuste nem convergia (theta 14% errado). Nada
// disso e o modelo: e o caminhante lendo unidades. O piso relativo do passo (L_ii >= 1e-3
// max L_jj no bloco) compara diagonais de caracteristicas diferentes; o colapso da curvatura
// (|AI_z,ii| <= 1e-10 do maior) compara uma coordenada log, sem unidade, com uma fora da
// diagonal, em unidades da caracteristica; o piso absoluto e o criterio relDelta somam
// variancias de unidades diferentes; as folgas do resgate EM e da estagnacao sao relativas a
// |-2logL|, que muda com a unidade. A caminhada inteira passa a acontecer com cada
// caracteristica dividida pela potencia de 2 logo acima do proprio desvio padrao
// (desvios_tracos), o que torna TODAS essas regras
// invariantes a unidade de uma vez: y1 e 0,01 y1 dao o mesmo caminho, e o resultado volta a
// unidade original no fim. O que se reporta (theta, -2logL, solucoes, PEV, EP) sai de uma
// avaliacao final nos dados ORIGINAIS no theta reescalado, nao de uma conversao: e o mesmo
// numero que eval_internal_mt() daria naquele ponto.
//
// y e reescalado NO LUGAR (sem copia do desenho, que carrega X e Z) e volta bit a bit: a
// copia guardada e de y so.
// (a guarda PadronizaY esta em mme.h, compartilhada com o AR(1))

// GANCHO DE TESTE, so para o portao do D2 (test-unidade-partida-fronteira.R): com a variavel
// de ambiente UAITREM_TESTE_SEM_PASSO definida, todo candidato do passo AI amortecido e todo
// passo EM sao recusados, e a saida "nenhum passo melhora" e exercitada de proposito. Nenhum
// dado conhecido chega nela sozinho desde que a partida alem dos pisos vai a eles (abaixo), e
// sem o gancho a mensagem dela ficava sem portao. Fora do teste a variavel nao existe e nada
// muda.
static bool teste_sem_passo() { return std::getenv("UAITREM_TESTE_SEM_PASSO") != nullptr; }

// A COVARIANCIA RESIDUAL DE DOIS CARACTERES QUE NENHUM REGISTRO OBSERVA JUNTOS nao entra na
// verossimilhanca: cada padrao de observacao usa so a submatriz de R0 dos seus caracteres.
// E o desenho de um caracter de registros repetidos com outro de registro unico, cada um na
// sua linha. Score e linha da AI dessa coordenada sao zero exato, a AI fica singular e a
// inversa dela deixava TODO erro-padrao em NaN, nao so o dessa covariancia (medido: as oito
// componentes de um bivariado animal + pe, todas NaN). Ela fica presa em zero e sai da
// inversa; 1 = fora da verossimilhanca, por posicao em theta.
static std::vector<char> residuos_sem_par(const DesenhoMT& d) {
  const std::size_t t = d.t;
  std::vector<char> junto(t * t, 0), fora(d.modelo.ntheta, 0);
  for (std::size_t i = 0; i < d.nlin; i++) {
    if (!d.usa[i]) continue;
    for (std::size_t a = 0; a < t; a++)
      for (std::size_t b = a + 1; b < t; b++)
        if (d.obs[i * t + a] && d.obs[i * t + b]) junto[a * t + b] = 1;
  }
  for (std::size_t a = 0; a < t; a++)
    for (std::size_t b = a + 1; b < t; b++)
      if (!junto[a * t + b])
        fora[d.modelo.offset_residual + (a * (2 * t - a + 1)) / 2 + (b - a)] = 1;
  return fora;
}

// a submatriz de m nas coordenadas k: a AI restrita as que entram na verossimilhanca
static Densa sub_densa(const Densa& m, const std::vector<std::size_t>& k) {
  Densa s(k.size(), k.size());
  for (std::size_t i = 0; i < k.size(); i++)
    for (std::size_t j = 0; j < k.size(); j++) s.at(i, j) = m.at(k[i], k[j]);
  return s;
}

// A CAMINHADA do multicaracter, de uma partida ate um veredito. ajusta_mt (no fim do arquivo)
// a chama e, quando ela termina com um bloco na parede, chama de novo de dentro do cone
// (explora_fronteira, aireml.cpp).
// `it0` e o numero de iteracoes ja gastas por uma caminhada anterior do mesmo ajuste: o laco
// conta a partir dela e para no mesmo maxiter, entao `iters` e as mensagens falam do total.
static AjusteMT caminha_mt(DesenhoMT& d, const std::vector<double>* theta0,
                           std::size_t maxiter, double tol, bool verboso, std::size_t it0) {
  AjusteMT R;
  R.iters = it0;
  const std::size_t t = d.t;
  const bool sem_passo = teste_sem_passo();
  std::vector<double> esc_exato;
  const std::vector<double> esc = desvios_tracos(d.y, d.usa, d.obs, t, &esc_exato);
  const std::vector<double> fator = fatores_escala(d.modelo, t, d.modelo.offset_residual, esc);
  // as unidades do corte de singularidade, com o desvio padrao exato (aviso_informacao)
  const std::vector<double> fator_corte =
      fatores_escala(d.modelo, t, d.modelo.offset_residual, esc_exato);
  // -2logL original = -2logL padronizado + jac, jac = 2 soma_tau (n_tau - p_tau) log s_tau
  // (n_tau celulas observadas, p_tau equacoes fixas da caracteristica): so para o TRACE,
  // porque o que se reporta vem da avaliacao final nos dados originais
  double jac = 0.0;
  for (std::size_t tau = 0; tau < t; tau++) {
    double n_tau = 0.0, p_tau = 0.0;
    for (std::size_t i = 0; i < d.nlin; i++) if (d.usa[i] && d.obs[i * t + tau]) n_tau++;
    for (std::size_t j = 0; j < d.x.ncol; j++) if (d.eq_fixa[j * t + tau] >= 0) p_tau++;
    jac += 2.0 * (n_tau - p_tau) * std::log(esc[tau]);
  }
  auto bruto = [&](const std::vector<double>& th) {
    std::vector<double> o = th;
    for (std::size_t k = 0; k < o.size() && k < fator.size(); k++) o[k] *= fator[k];
    return o;
  };
  PadronizaY guarda(d.y, d.usa, d.obs, esc, t);
  // start= do usuario, ou a partida automatica. Vale a mesma razao do univariado: a
  // maneira de conferir que o otimo nao depende de onde a busca comecou, e a saida que a
  // mensagem de parada recomenda. O start= chega em unidades originais e e levado para as
  // padronizadas; a partida automatica ja nasce nelas (var(y) = 1).
  std::vector<double> theta;
  if (theta0) {
    theta = *theta0;
    for (std::size_t k = 0; k < theta.size() && k < fator.size(); k++) theta[k] /= fator[k];
  } else {
    theta = partida_mt(d);
  }
  const std::vector<char> sem_par = residuos_sem_par(d);
  for (std::size_t k = 0; k < theta.size(); k++) if (sem_par[k]) theta[k] = 0.0;
  R.partida = bruto(theta);
  // o cache vai junto com o ajuste: residuos_mt() reaproveita a ordem e o simbolico no fim
  // (a ordenacao de grau minimo, serial, era ~17% de um model_mt de 100 mil; E0, 2026-11-02)
  R.cache = std::make_shared<CacheSimbolica>();
  CacheSimbolica& cs = *R.cache;
  AvaliacaoMT cur = avalia_mt(d, theta, &cs);
  if (cs.pronto) { R.bloco_denso = bloco_denso_simbolico(cs.sb); R.colunas_fator = cs.sb.n; }
  if (!cur.ok) {
    // QUAL bloco (lista de correcoes, item 13): a recusa nao dizia, e o ajuste voltava com
    // -2logL 0, que parece um numero
    std::string quais;
    if (cur.motivo != 2) {
      const std::vector<double> th = bruto(theta);
      auto pd = [&](std::size_t off, std::size_t dim) {
        Densa c(dim, dim);
        for (std::size_t j = 0; j < dim; j++)
          for (std::size_t i = j; i < dim; i++) {
            const double v = th[off + (j * (2 * dim - j + 1)) / 2 + (i - j)];
            c.at(i, j) = v;
            c.at(j, i) = v;
          }
        return chol_densa(c);
      };
      for (const Grupo& g : d.modelo.grupos)
        if (!pd(g.offset, g.dim)) quais += (quais.empty() ? "" : ", ") + std::string("'") + g.nome + "'";
      if (!pd(d.modelo.offset_residual, t))
        quais += (quais.empty() ? "" : ", ") + std::string("the residual R0");
    }
    R.mensagem = cur.motivo == 2
      ? "the coefficient matrix is SINGULAR at the starting theta: the design does not identify "
        "some effect, and changing start= does not help"
      : "the starting theta is INADMISSIBLE: the covariance matrix of " +
        (quais.empty() ? std::string("some group") : quais) + " is not positive-definite";
    R.neg2logl = std::nan("");
    return R;
  }

  // O PASSO anda em log-Cholesky, nao em theta cru, e a troca nao e cosmetica. O teste que
  // a mede e a EQUIVARIANCIA: -2logL e invariante sob (K, C) -> (cK, C/c), logo o otimo
  // AJUSTADO tem de ser tambem, e isso e propriedade do caminhante, nao da verossimilhanca.
  // Andando em theta cru, com a mesma K escalada por 25, o ajuste terminava 146 unidades de
  // -2logL longe do de c = 1 (453.55 contra 307.20), com o componente errado por duas
  // ordens de grandeza e converged TRUE nos dois. Em z todo ponto e admissivel, a fronteira
  // det(C_g) = 0 vai para o infinito e a escala de cada direcao acompanha a do proprio
  // componente. Aqui o residuo e um BLOCO de t x t, e nao o escalar do univariado.
  MapaZ mz;
  for (const Grupo& gr : d.modelo.grupos) mz.blocos.push_back(std::make_pair(gr.offset, gr.dim));
  mz.blocos.push_back(std::make_pair(d.modelo.offset_residual, d.t));

  // O certificado final, espelhado do ajustador univariado (o racional completo esta em
  // aireml.cpp): passo relativo pequeno nao prova otimo, AI amortecida pode aceitar um
  // passo minusculo com o gradiente longe de zero. converged exige tambem o decremento
  // de Newton, g' AI^-1 g, ~2x o gap em -2logL perto do otimo. Este laco anda em theta
  // cru, sem pisos e sem log-Cholesky (o limite declarado da etapa 1B), entao as
  // fronteiras que o univariado congela via conjunto ativo aparecem aqui como score que
  // nunca zera. O certificado exclui o que este laco consegue reconhecer: uma variancia
  // diagonal em zero numerico (<= 1e-6 da escala residual) com o score empurrando para
  // baixo, e TODO componente de um grupo cuja C_g esta a 1e-3 (razao das diagonais do
  // Cholesky, o mesmo piso relativo do univariado) da singularidade, este laco nao
  // anda sobre a parede, entao a certificacao e CONDICIONAL a ela e a mensagem diz
  // isso. AI que nao inverte no bloco livre nao certifica nada: conta como recusa.
  const double tol_dec = 2e-4;
  // Quantos componentes o certificado excluiu na ultima avaliacao, para a mensagem. O
  // conjunto sai de PassoZ::na_parede, nas mesmas coordenadas e sobre os mesmos pisos
  // que o passo usa; o teste em theta que havia aqui (lmin < 1e-3 lmax) era a exata
  // desigualdade que o grampo do passo torna falsa, e por isso nunca disparava.
  std::size_t na_fronteira = 0;
  // a escala residual e RECALCULADA a cada ponto, como no univariado: e ela que ancora o
  // piso absoluto, e ela se move durante o ajuste
  auto pecas = [&](const AvaliacaoMT& av) {
    double s = 0.0;
    for (std::size_t tau = 0; tau < d.t; tau++)
      s += theta[d.modelo.offset_residual + (tau * (2 * d.t - tau + 1)) / 2];
    PassoZ P = passo_z(mz, d.modelo.ntheta, theta, av.score, av.ai,
                       std::max(1.0, s / static_cast<double>(d.t)));
    // a covariancia residual sem par (residuos_sem_par) nao e uma fronteira: esta fora da
    // verossimilhanca e presa em zero. Sem isto a curvatura nula a punha na parede, a
    // mensagem dizia "1 component resting at a covariance boundary" e a segunda caminhada
    // de dentro do cone rodava por nada
    if (P.ok)
      for (std::size_t k = 0; k < sem_par.size(); k++)
        if (sem_par[k]) { P.na_parede[k] = 0; P.congelado[k] = 1; }
    return P;
  };
  // RESGATE EM, o que faltava para os espelhos terem a maquinaria completa do univariado.
  // O passo EM e multiplicativo, fica no cone e nao encolhe na fronteira, entao ele anda
  // exatamente onde o passo AI amortecido trava. Aqui ele move so os GRUPOS: o passo M do
  // R0 multivariado, e o do rho no AR(1), nao existem em forma fechada, e por isso a
  // proposta deixa esses componentes onde estao. Isso e um EM GENERALIZADO e continua
  // valido, porque cada candidato so entra se BAIXAR a verossimilhanca.
  //
  // Cada candidato passa pelo grampo em z antes de ser avaliado: um EM que encoste
  // numericamente na singularidade daria um -2logL sem sentido, com det(C_g) em underflow.
  // Ate 15 passos, contra os 50 do univariado, porque aqui cada avaliacao custa mais.
  auto resgate_em = [&](AvaliacaoMT& atual) -> bool {
    bool ganhou = false;
    if (sem_passo) return false;
    for (int e = 0; e < 15; e++) {
      if (atual.em_theta.size() != d.modelo.ntheta) break;
      std::vector<double> ze;
      if (!z_de_theta(mz, d.modelo.ntheta, atual.em_theta, ze)) break;
      const PassoZ Pe = pecas(atual);
      if (Pe.ok)
        for (std::size_t i = 0; i < d.modelo.ntheta; i++)
          if (ze[i] < Pe.piso[i]) ze[i] = Pe.piso[i];
      std::vector<double> cand = theta_de_z(mz, d.modelo.ntheta, ze);
      AvaliacaoMT prox = avalia_mt(d, cand, &cs);
      // ganho RELEVANTE, e nao qualquer migalha. Este EM move so os grupos: quando o
      // que falta esta no R0 ele nao alcanca, e aceitar um ganho de 1e-8 so desvia um
      // caminho saudavel para um ponto de onde o passo AI nao sai (medido: uma celula
      // que certificava com decremento 1.3e-04 passou a parar em 8.4e-04). Com um piso
      // relativo, o EM entra onde ele de fato resolve, a fronteira, onde o passo
      // amortecido trava, e sai da frente onde nao resolve.
      const double ganho_em = atual.neg2logl - prox.neg2logl;
      if (!prox.ok || ganho_em < 1e-6 * std::max(1.0, std::fabs(atual.neg2logl))) break;
      double num = 0.0, den = 0.0;
      for (std::size_t i = 0; i < cand.size(); i++) {
        const double dlt = cand[i] - theta[i];
        num += dlt * dlt;
        den += cand[i] * cand[i];
      }
      R.reldelta = std::sqrt(num / std::max(den, 1e-300));
      theta = cand;
      atual = std::move(prox);
      ganhou = true;
    }
    return ganhou;
  };
  // o conjunto excluido do ultimo certificado, para a mensagem e para a exploracao
  std::vector<char> parede_final;
  // as DIAGONAIS de blocos de varias dimensoes no piso (a regra de 2 vezes o piso do model()),
  // que o EP condicional segura; parede_final e o conjunto inteiro que o certificado excluiu, e
  // o passo congela o bloco inteiro quando a curvatura colapsa
  std::vector<char> piso_final;
  auto conta_parede = [&](const PassoZ& P) {
    parede_final = P.na_parede;
    piso_final.assign(P.zc.size(), 0);
    if (P.ok)
      for (const auto& bl : mz.blocos) {
        if (bl.second == 1) continue;
        for (std::size_t i = 0; i < bl.second; i++) {
          const std::size_t k = bl.first + (i * (2 * bl.second - i + 1)) / 2;
          if (P.zc[k] <= P.piso[k] + std::log(2.0)) piso_final[k] = 1;
        }
      }
    std::size_t n = 0;
    for (char c : P.na_parede) n += c;
    return n;
  };
  // A SONDA DE KKT pela verossimilhanca (o racional e a falha medida estao em
  // candidatos_reentrada, src/aireml.cpp): o certificado deste laco exclui o bloco na parede
  // sem olhar o sinal, e o passo congela a coordenada cuja curvatura em z colapsou. Antes de
  // aceitar qualquer certificado com algum bloco na parede (o do passo pequeno, o da
  // estagnacao e o da crista plana no fim do orcamento), pontos dentro do cone sao
  // avaliados: as direcionais primeiro, as por coordenada so se nenhuma direcional melhorar.
  // Um ponto que baixa o -2logL por mais que tol_reentrada recusa o certificado; com
  // `pode_mover` e abaixo do limite de reentradas o caminhante vai para ele e continua, senao
  // o veredito e `limite` e o ajuste sai com converged = FALSE dizendo por que.
  int n_reentradas = 0;
  double ganho_kkt = 0.0;   // quanto o melhor ponto da ultima sonda recusada baixava
  // A SONDA DA FACE (sonda_face, aireml.cpp), depois da do cone: o certificado tira da conta
  // o bloco inteiro na parede, e um ponto da face que nao e o otimo dela passava junto. Os
  // movimentos sobre a face tem limite proprio, 20, porque cada um e um passo de Newton das
  // diferencas finitas e o caminhante volta a andar depois dele.
  int n_face = 0;
  // a avaliacao em z das sondas da face: o -2logL e o score levado a z pelo Jacobiano
  auto avalia_z = [&](const std::vector<double>& zc, double& v, std::vector<double>& gz) {
    const std::size_t nt = d.modelo.ntheta;
    const AvaliacaoMT av = avalia_mt(d, theta_de_z(mz, nt, zc), &cs);
    if (!av.ok) return false;
    v = av.neg2logl;
    const Densa J = jacobiano_z(mz, nt, zc);
    gz.assign(nt, 0.0);
    for (std::size_t m = 0; m < nt; m++)
      for (std::size_t k = 0; k < nt; k++) gz[m] += J.at(k, m) * av.score[k];
    return true;
  };
  auto sonda_da_face = [&](const PassoZ& P, bool pode_mover) -> Kkt {
    const std::size_t nt = d.modelo.ntheta;
    std::vector<double> zn;
    double fn = 0.0;
    if (!sonda_face(mz, nt, P, cur.neg2logl, avalia_z, zn, fn)) return Kkt::otimo;
    if (!(fn < cur.neg2logl - tol_reentrada)) return Kkt::otimo;
    ganho_kkt = cur.neg2logl - fn;
    if (!pode_mover || n_face >= 20) return Kkt::limite;
    const std::vector<double> cand = theta_de_z(mz, nt, zn);
    AvaliacaoMT av = avalia_mt(d, cand, &cs);
    if (!av.ok) return Kkt::otimo;
    double num = 0.0, den = 0.0;
    for (std::size_t i = 0; i < cand.size(); i++) {
      num += (cand[i] - theta[i]) * (cand[i] - theta[i]);
      den += cand[i] * cand[i];
    }
    R.reldelta = std::sqrt(num / std::max(den, 1e-300));
    if (verboso)
      Rprintf("      along the boundary face: -2logL %.6f -> %.6f (the face was not at its "
              "optimum)\n", cur.neg2logl + jac, av.neg2logl + jac);
    theta = cand;
    cur = std::move(av);
    n_face++;
    return Kkt::moveu;
  };
  // A LAPIDACAO DA FACE (o racional esta em tol_lapida, mme.h), so DEPOIS de um certificado
  // aceito com algum bloco na parede: a cada rodada, primeiro o ponto com as covariancias de
  // uma caracteristica no zero zeradas (descorrelada_na_parede, aireml.cpp), depois o passo de
  // Newton da sonda da face; cada um entra se baixa o -2logL por mais que tol_lapida E o ponto
  // novo ainda certifica (decremento dos livres abaixo de tol_dec, nas pecas do proprio
  // ponto). O primeiro ponto que nao certifica fica de fora e o ajuste fica no ultimo que
  // certificou: a lapidacao nunca troca um TRUE por um FALSE e nunca sobe o -2logL. Tentado
  // antes e desfeito: deixar o caminhante andar depois de cada movimento de lapidacao. Na
  // semente 10 (r12 = -0,25) da varredura da revisao, a partida de var(y) foi de TRUE 3,86e-4
  // acima do melhor para FALSE 1,75e-5 acima: o passo AI levava a diagonal na parede de volta
  // para cima do piso (L_22 / L_11 de 1e-3 a 5e-3), o bloco deixava de estar na parede, o
  // score dele perto da singularidade entrava no decremento (0,008) e o caminhante estagnava.
  auto lapida = [&](PassoZ& Pc) {
    const std::size_t nt = d.modelo.ntheta;
    // um candidato vale se baixa por mais que tol_lapida e o ponto novo ainda certifica
    auto aceita = [&](const std::vector<double>& zn, const char* como) -> bool {
      const std::vector<double> cand = theta_de_z(mz, nt, zn);
      AvaliacaoMT av = avalia_mt(d, cand, &cs);
      if (!av.ok || !(av.neg2logl < cur.neg2logl - tol_lapida)) return false;
      const std::vector<double> antes = theta;
      theta = cand;   // pecas le o theta corrente
      PassoZ Pn = pecas(av);
      const double dn = Pn.ok ? decremento_z(Pn) : tol_dec;
      if (!(dn < tol_dec)) {
        theta = antes;
        return false;
      }
      if (verboso)
        Rprintf("      polishing at the boundary, %s: -2logL %.6f -> %.6f (inside the 1e-4 "
                "slack; the new point still certifies)\n", como, cur.neg2logl + jac,
                av.neg2logl + jac);
      cur = std::move(av);
      Pc = std::move(Pn);
      R.decremento = dn;
      return true;
    };
    for (int k = 0; k < max_lapidacoes; k++) {
      const std::vector<double> zd = descorrelada_na_parede(mz, nt, Pc);
      if (!zd.empty() && aceita(zd, "a trait at zero decorrelated")) continue;
      std::vector<double> zn;
      double fn = 0.0;
      if (!sonda_face(mz, nt, Pc, cur.neg2logl, avalia_z, zn, fn)) return;
      if (!(fn < cur.neg2logl - tol_lapida)) return;
      if (!aceita(zn, "along the face")) return;
    }
  };
  auto sonda_kkt = [&](const PassoZ& P, bool pode_mover) -> Kkt {
    if (!P.ok) return Kkt::otimo;
    bool alguma = false;
    for (char c : P.na_parede) if (c) alguma = true;
    if (!alguma) return Kkt::otimo;
    const std::size_t nt = d.modelo.ntheta;
    double esc_res = 0.0;
    for (std::size_t tau = 0; tau < t; tau++)
      esc_res += theta[d.modelo.offset_residual + (tau * (2 * t - tau + 1)) / 2];
    const CandidatosKKT cands = candidatos_reentrada(
        mz, nt, P.zc, P.na_parede, std::max(1.0, esc_res / static_cast<double>(t)),
        cur.score, cur.ai);
    std::vector<double> z_melhor;
    AvaliacaoMT av_melhor;
    bool achou = false;
    for (const std::vector<std::vector<double> >* lista : {&cands.direcionais, &cands.coordenadas}) {
      for (const std::vector<double>& zc : *lista) {
        AvaliacaoMT av = avalia_mt(d, theta_de_z(mz, nt, zc), &cs);
        if (!av.ok) continue;
        if (!achou || av.neg2logl < av_melhor.neg2logl) {
          achou = true;
          z_melhor = zc;
          av_melhor = std::move(av);
        }
      }
      if (achou && av_melhor.neg2logl < cur.neg2logl - tol_reentrada) break;
    }
    if (!achou || !(av_melhor.neg2logl < cur.neg2logl - tol_reentrada))
      return sonda_da_face(P, pode_mover);
    ganho_kkt = cur.neg2logl - av_melhor.neg2logl;
    if (!pode_mover || n_reentradas >= max_reentradas) return Kkt::limite;
    const std::vector<double> cand = theta_de_z(mz, nt, z_melhor);
    double num = 0.0, den = 0.0;
    for (std::size_t i = 0; i < cand.size(); i++) {
      num += (cand[i] - theta[i]) * (cand[i] - theta[i]);
      den += cand[i] * cand[i];
    }
    R.reldelta = std::sqrt(num / std::max(den, 1e-300));
    if (verboso)
      Rprintf("      re-entry into the cone: -2logL %.6f -> %.6f (a component at the "
              "boundary was not at a constrained optimum)\n", cur.neg2logl + jac,
              av_melhor.neg2logl + jac);
    theta = cand;
    cur = std::move(av_melhor);
    n_reentradas++;
    return Kkt::moveu;
  };
  // a mensagem do certificado recusado pela sonda quando nao ha como seguir
  auto msg_limite = [&](const std::string& onde) {
    char buf[96];
    std::snprintf(buf, sizeof(buf), "%.3g", ganho_kkt);
    return std::string("did NOT converge: ") + onde + " a component resting at a covariance "
        "boundary is not at a constrained optimum (a point inside the cone or along the "
        "boundary face lowers the -2logL by " + buf + "), and the walk had already moved "
        "off it the most times allowed (" + std::to_string(max_reentradas) + " re-entries "
        "into the cone, 20 moves along the face). Try a start= inside the cone, for example "
        "the univariate fits with the covariances at zero";
  };

  // UMA PARTIDA ALEM DOS PISOS DO PASSO vai para eles antes da caminhada (so quando ha
  // caminhada: com maxiter = 0 avalia-se o start= como foi dado). O passo grampeia todo
  // candidato nos pisos (diagonal de Cholesky >= 1e-3 da maior do bloco, variancia >= 1e-10
  // da escala), entao o caminhante nao representa um ponto abaixo deles, e ali o MME perde
  // digitos: partindo do otimo de fronteira com 1 - r^2 = 1e-12 no bloco da ninhada, o
  // -2logL esparso dava 0,0028 unidade abaixo da forma V densa, com 1e-15 51,3 abaixo, e
  // como todo candidato grampeado "piorava" contra esse valor, o ajuste parava na primeira
  // iteracao. Levada aos pisos, a partida e avaliada com os digitos todos e a caminhada
  // segue; `start` devolve o ponto de onde ela de fato saiu.
  if (maxiter > 0) {
    const PassoZ P0 = pecas(cur);
    if (P0.ok) {
      std::vector<double> z0 = P0.zc;
      bool mudou = false;
      for (std::size_t i = 0; i < z0.size(); i++)
        if (z0[i] < P0.piso[i]) { z0[i] = P0.piso[i]; mudou = true; }
      if (mudou) {
        const std::vector<double> th0 = theta_de_z(mz, d.modelo.ntheta, z0);
        AvaliacaoMT a0 = avalia_mt(d, th0, &cs);
        if (a0.ok) {
          if (verboso)
            Rprintf("start beyond the step floors, raised to them: -2logL %.6f -> %.6f\n",
                    cur.neg2logl + jac, a0.neg2logl + jac);
          theta = th0;
          cur = std::move(a0);
          R.partida = bruto(theta);
        }
      }
    }
  }

  double lambda = 1e-2;   // multi comeca mais amortecido: as covariancias partem de zero
  int parado = 0;         // iteracoes consecutivas aceitas SEM progresso real
  bool saiu_do_laco = false;   // o laco parou por um veredito, e nao por esgotar o maxiter
  for (std::size_t it = it0 + 1; it <= maxiter; it++) {
    R.iters = it;
    R_CheckUserInterrupt();
    bool aceitou = false;
    const std::size_t ntz = d.modelo.ntheta;
    // z, score e AI em z, pisos e os dois conjuntos ativos, numa chamada compartilhada
    // pelos tres ajustadores. O certificado la embaixo le as MESMAS pecas.
    PassoZ P = pecas(cur);
    // O PONTO CORRENTE ABAIXO DOS PROPRIOS PISOS vai a eles antes do passo. O piso relativo
    // anda com a maior diagonal do bloco: um passo aceito que a sobe deixa a diagonal na
    // parede (congelada no passo) um fio abaixo do piso novo, e todo candidato da iteracao
    // seguinte carrega o mesmo salto ate ele. Perto da singularidade o -2logL tem ruido maior
    // que a folga de aceite (1e-9), e o salto fazia TODO candidato piorar: medido em
    // test-multitrait.R (APY com nucleo = todos), a caminhada parava na iteracao 12 em
    // "nenhum passo melhora", com a diagonal 5e-6 (em log) abaixo do piso. O salto e aceito
    // sem comparar, porque o caminhante nao representa pontos abaixo dos pisos (a mesma
    // regra da partida, acima). Subir uma diagonal abaixo do piso relativo nao muda a maior
    // do bloco, entao uma vez basta.
    if (P.ok) {
      std::vector<double> z0 = P.zc;
      bool abaixo = false;
      for (std::size_t i = 0; i < z0.size(); i++)
        if (z0[i] < P.piso[i]) { z0[i] = P.piso[i]; abaixo = true; }
      if (abaixo) {
        const std::vector<double> th0 = theta_de_z(mz, ntz, z0);
        AvaliacaoMT a0 = avalia_mt(d, th0, &cs);
        if (a0.ok) {
          theta = th0;
          cur = std::move(a0);
          P = pecas(cur);
        }
      }
    }
    for (int tent = 0; tent < 30 && !aceitou; tent++) {
      Densa m2 = P.ok ? P.az : cur.ai;
      for (std::size_t i = 0; i < m2.nlin; i++) {
        const double di = m2.at(i, i);
        m2.at(i, i) = (di == 0.0) ? lambda : di * (1.0 + lambda);
      }
      if (P.ok)
        for (std::size_t i = 0; i < ntz; i++)
          if (P.congelado[i]) {
            for (std::size_t j = 0; j < ntz; j++) { m2.at(i, j) = 0.0; m2.at(j, i) = 0.0; }
            m2.at(i, i) = 1.0;
          }
      Densa minv;
      try { minv = inv_geral(m2); } catch (const Erro&) { lambda *= 10.0; continue; }
      // a DIRECAO, computada uma vez; o comprimento e que e buscado
      std::vector<double> passo(ntz, 0.0);
      for (std::size_t i = 0; i < ntz; i++)
        for (std::size_t j = 0; j < ntz; j++)
          if (!P.ok || !P.congelado[j]) passo[i] += minv.at(i, j) * (P.ok ? P.sz[j] : cur.score[j]);
      // BUSCA DE COMPRIMENTO na mesma direcao, antes de mexer no amortecimento. Sem ela
      // toda reducao de passo passava por multiplicar lambda por dez, o que encurta E
      // GIRA a direcao: o laco alternava entre dois estados e a verossimilhanca ficava
      // parada (medido: 980 iteracoes com -2logL identico ate a ultima casa).
      for (const double alpha : {1.0, 0.5, 0.25}) {
        std::vector<double> cand;
        if (P.ok) {
          std::vector<double> zn = P.zc;
          for (std::size_t i = 0; i < ntz; i++)
            if (!P.congelado[i]) zn[i] -= alpha * passo[i];
          // PISO na diagonal de cada bloco. Sem ele o passo em z vai ate a fronteira
          // exata (autovalor de G0 medido em 1.5e-14 num ajuste sem sinal genetico
          // nenhum), e ali a inversa do bloco que o MME precisa nao existe mais: os
          // efeitos fixos passam a divergir da mesma GLS montada densa na quinta casa.
          for (std::size_t i = 0; i < ntz; i++)
            if (zn[i] < P.piso[i]) zn[i] = P.piso[i];
          cand = theta_de_z(mz, ntz, zn);
        } else {
          // theta fora do cone (uma C_g ja singular): o passo cru ainda serve de resgate
          cand = theta;
          for (std::size_t i = 0; i < ntz; i++) cand[i] -= alpha * passo[i];
        }
        if (sem_passo) continue;
        AvaliacaoMT prox = avalia_mt(d, cand, &cs);
        if (!prox.ok || prox.neg2logl > cur.neg2logl + 1e-9) continue;
        double num = 0.0, den = 0.0;
        for (std::size_t i = 0; i < cand.size(); i++) {
          const double dlt = cand[i] - theta[i];
          num += dlt * dlt;
          den += cand[i] * cand[i];
        }
        // PROGRESSO nao e o mesmo que aceitacao: o teste de aceitacao tem folga de 1e-9
        // e por isso aceita um passo que nao anda nada. Contar esses separadamente e o
        // que impede o laco de virar ponto fixo ate o maxiter.
        const double ganho = cur.neg2logl - prox.neg2logl;
        parado = (ganho > 1e-8 * std::max(1.0, std::fabs(cur.neg2logl))) ? 0 : parado + 1;
        R.reldelta = std::sqrt(num / std::max(den, 1e-300));
        theta = cand;
        cur = std::move(prox);
        lambda = std::max(lambda / 10.0, 1e-10);
        aceitou = true;
      }
      if (!aceitou) lambda *= 10.0;
    }
    if (!aceitou) {
      // A CAUSA DE VERDADE, em ingles e com a iteracao. Esta saida deixava a string interna
      // "nenhum passo amortecido melhorou a verossimilhanca" na mensagem, e o bloco do
      // maxiter ainda prefixava "stopped at 300 iteration(s)" num ajuste que tinha parado
      // na 13a. Antes de desistir, o resgate EM (o mesmo do univariado) e a sonda de KKT.
      if (resgate_em(cur)) {
        lambda = 1e-2;
        parado = 0;
        if (verboso)
          Rprintf("      EM rescue: -2logL %.6f (no damped step improved)\n", cur.neg2logl + jac);
        continue;
      }
      // A ESTAGNACAO NAO CERTIFICA, nem com o decremento pequeno (o model() certifica ali,
      // mas o decremento dele conta toda coordenada fora da parede com a regra do sinal; o
      // daqui exclui os blocos na parede sem sinal). O aceite tolera 1e-9 de piora, entao
      // num otimo numerico o passo continua sendo aceito sem ganho e e o `parado` que leva
      // ao certificado; esta saida so acontece quando TODO candidato piora. O unico caso
      // medido era uma partida alem dos pisos do passo (cada candidato grampeado de volta a
      // eles "piorava" contra um -2logL sem digitos, ver a partida acima), que agora vai aos
      // pisos antes da caminhada; nenhum outro caso e conhecido, e a suite inteira nao passa
      // por aqui. A versao anterior certificava aqui quando o decremento era pequeno; o base
      // dava sempre FALSE, com a string interna em portugues. Fica FALSE, com a causa e a
      // iteracao.
      const PassoZ Pq = pecas(cur);
      const Kkt kq = sonda_kkt(Pq, true);
      if (kq == Kkt::moveu) {
        lambda = 1e-2;
        parado = 0;
        continue;
      }
      saiu_do_laco = true;
      na_fronteira = conta_parede(Pq);
      R.decremento = decremento_z(Pq);
      if (kq == Kkt::limite) {
        R.mensagem = msg_limite("at iteration " + std::to_string(it) + " no damped step "
                                "improved the likelihood, and");
      } else {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "%.3g", R.decremento);
        R.mensagem = std::string("did NOT converge: at iteration ") + std::to_string(it) +
            " no damped AI step (30 damping levels, three step lengths each) and no EM "
            "step improved the likelihood, and the Newton decrement of the off-boundary "
            "components is " + buf + (R.decremento < tol_dec
              ? std::string(" (inside the 2e-4 tolerance, but a point no step can "
                            "leave is not certified as an optimum). Try a different "
                            "start=, and profile any component the message reports at a "
                            "boundary")
              : std::string(" (tolerance 2e-4 on the -2logL scale), so the point is a "
                            "stall and not a certified optimum. Try a different start=, "
                            "and profile any component the message reports at a "
                            "boundary"));
      }
      break;
    }
    if (verboso) {
      Rprintf("iter %3d  -2logL %.6f  relDelta %.3e\n",
              (int) it, cur.neg2logl + jac, R.reldelta);
      imprime_theta(bruto(theta), nomes_theta_do_mt(d), d.modelo);
    }
    if (R.reldelta < tol || parado >= 2) {
      // O CERTIFICADO, nas mesmas coordenadas e sobre os mesmos pisos que o passo. Ele
      // ficou para tras quando o passo foi portado para log-Cholesky, e a assimetria era
      // um defeito de verdade: o teste que havia aqui, lmin < 1e-3 lmax em theta, e a
      // exata desigualdade que o grampo do passo torna FALSA por construcao. No ponto em
      // que o laco para, o grampo deixa lmin = 1e-3 lmax e o `<` estrito nunca dispara,
      // de modo que a direcao grampeada ficava congelada no passo e cobrada integralmente
      // no certificado. Medido: decremento parado em 3.25e+03 com -2logL identico por 980
      // iteracoes, num ponto que uma descida por coordenada melhorava em 0.200 unidade.
      //
      // converged exige as duas coisas, como no univariado: passo pequeno E decremento de
      // Newton dos componentes LIVRES dentro da tolerancia.
      // o EM primeiro: se ele ainda anda, nao havia convergencia nenhuma, e o
      // amortecimento reinicia para o passo AI voltar a tentar do ponto novo
      if (resgate_em(cur)) {
        lambda = 1e-2;
        parado = 0;
        if (verboso)
          Rprintf("      EM rescue: -2logL %.6f (the AI step had stalled)\n", cur.neg2logl + jac);
        continue;
      }
      PassoZ Pc = pecas(cur);
      na_fronteira = conta_parede(Pc);
      R.decremento = decremento_z(Pc);
      if (R.decremento < tol_dec) {
        // o certificado so vale se nenhum bloco na parede melhora entrando no cone
        const Kkt kc = sonda_kkt(Pc, true);
        if (kc == Kkt::moveu) {
          lambda = 1e-2;
          parado = 0;
          R.decremento = std::numeric_limits<double>::quiet_NaN();
          continue;
        }
        saiu_do_laco = true;
        if (kc == Kkt::limite) {
          R.mensagem = msg_limite("the step had converged, but");
          break;
        }
        lapida(Pc);
        na_fronteira = conta_parede(Pc);
        R.convergiu = true;
        break;
      }
      // A recusa nao encerra o ajuste: o passo pode voltar a andar. Mas ele so volta se o
      // amortecimento for solto, depois de uma sequencia de aceites lambda esta no piso
      // e a direcao e sempre a mesma. `parado` conta os aceites sem ganho real e desiste
      // quando nem soltar o amortecimento adianta.
      // O corte de estagnacao ficou mais paciente depois que o resgate EM entrou, porque
      // "nada anda" passou a ser uma afirmacao mais forte: o passo AI travou E o EM
      // tambem nao move. Um EM bem sucedido reinicia o amortecimento, e o passo AI precisa
      // de espaco para retomar do ponto novo antes de o laco declarar estagnacao.
      const Kkt kp = parado >= 25 ? sonda_kkt(Pc, true) : Kkt::otimo;
      if (kp == Kkt::moveu) {
        lambda = 1e-2;
        parado = 0;
        continue;
      }
      if (parado >= 25) {
        saiu_do_laco = true;
        char buf[64];
        std::snprintf(buf, sizeof(buf), "%.3g", R.decremento);
        R.mensagem = std::string("did NOT converge: the step stopped making progress, ") +
            "but the Newton decrement of the off-boundary components is " + buf +
            " (tolerance 2e-4 on the -2logL scale), so the point is a stall and not a "
            "certified optimum. Try a different start=, and profile any component the "
            "message reports at a boundary";
        break;
      }
      lambda = 1e-2;
    }
  }
  if (std::isnan(R.decremento)) {
    const PassoZ Pf = pecas(cur);
    na_fronteira = conta_parede(Pf);
    R.decremento = decremento_z(Pf);
  }
  // A CRISTA PLANA NO FIM DO ORCAMENTO (o veredito do model(), abaixo) tambem passa pela
  // sonda de KKT: sem ela, um bloco congelado numa fronteira que nao e otimo saia com
  // converged = TRUE so por ter esgotado o maxiter. Aqui a sonda nao move nada, porque o
  // orcamento acabou: um ponto melhor dentro do cone recusa o certificado.
  const bool crista = !R.convergiu && !saiu_do_laco && maxiter > 0 && R.iters >= maxiter &&
                      R.decremento < tol_dec;
  const Kkt k_crista = crista ? sonda_kkt(pecas(cur), false) : Kkt::otimo;
  // DE VOLTA AS UNIDADES ORIGINAIS: y restaurado bit a bit e uma avaliacao final no theta
  // reescalado. O decremento e o conjunto na parede ficam os da caminhada (sao invariantes
  // a reescala); tudo o que se reporta em unidades sai daqui.
  guarda.restaura();
  theta = bruto(theta);
  bool final_ok = true;
  {
    AvaliacaoMT fin = avalia_mt(d, theta, &cs);
    if (fin.ok) cur = std::move(fin);
    else final_ok = false;
  }
  if (na_fronteira > 0)
    R.mensagem += std::string(R.mensagem.empty() ? "" : "; ") +
        std::to_string(na_fronteira) + " component(s) resting at a covariance boundary "
        "excluded from the convergence certificate: a component pinned at a boundary "
        "points out of the cone by construction, so its gradient never vanishes and the "
        "certificate is CONDITIONAL on the pinning; read those components with care: " +
        lista_na_parede(parede_final, nomes_theta_mt(d.modelo, d.alvos));
  // O EP CONDICIONAL (lista de correcoes, item 17), a regra do model() (aireml.cpp): o
  // componente de um bloco ESCALAR (1 x 1 no mapa do passo) parado na parede do certificado
  // nao foi estimado, foi segurado la. Ele fica sem EP, a linha e a coluna dele na vcov sao 0
  // (ali ele entra como constante), e os outros saem da informacao condicional a ele, a AI sem
  // a linha e a coluna dele. Antes a AI inteira dava EP explodidos aos outros. O bloco de
  // varias dimensoes na parede (uma correlacao +-1) usa a informacao condicional a parede, no
  // mapa do passo, como no model() (vcov_condicional, aireml.cpp). Com uma covariancia residual
  // sem registro comum (sem_par, segurada em 0 no theta e nao no z) fica a AI inteira
  R.se.assign(d.modelo.ntheta, std::nan(""));
  const std::size_t ntv = d.modelo.ntheta;
  std::vector<char> dado(ntv, 0);
  bool bloco_na_parede = false;
  if (parede_final.size() == ntv)
    for (const auto& bl : mz.blocos)
      if (bl.second == 1 && parede_final[bl.first]) dado[bl.first] = 1;
  if (piso_final.size() == ntv)
    for (std::size_t k = 0; k < ntv; k++)
      if (piso_final[k]) bloco_na_parede = true;
  bool algum_sem_par = false;
  for (std::size_t k = 0; k < sem_par.size(); k++)
    if (sem_par[k]) algum_sem_par = true;
  std::vector<std::size_t> na_veross;   // as coordenadas que entram na verossimilhanca, livres
  for (std::size_t k = 0; k < d.modelo.ntheta; k++)
    if (!sem_par[k] && !dado[k]) na_veross.push_back(k);
  R.vcov.assign(ntv * ntv, std::nan(""));
  for (std::size_t k = 0; k < ntv; k++)
    if (dado[k])
      for (std::size_t j = 0; j < ntv; j++) R.vcov[k * ntv + j] = R.vcov[j * ntv + k] = 0.0;
  std::vector<double> zf;
  if (final_ok && bloco_na_parede && !algum_sem_par && z_de_theta(mz, ntv, theta, zf)) {
    std::vector<std::size_t> fz;
    for (std::size_t m = 0; m < ntv; m++)
      if (!piso_final[m] && !dado[m]) fz.push_back(m);
    vcov_condicional(cur.ai, jacobiano_z(mz, ntv, zf), fz, na_veross, R.vcov, R.se);
  } else if (final_ok)
    try {
      Densa inv = inv_geral(sub_densa(cur.ai, na_veross));
      for (std::size_t i = 0; i < na_veross.size(); i++) {
        const double v = 2.0 * inv.at(i, i);
        if (v > 0.0) R.se[na_veross[i]] = std::sqrt(v);
        for (std::size_t j = 0; j < na_veross.size(); j++)
          R.vcov[na_veross[i] * ntv + na_veross[j]] = 2.0 * inv.at(i, j);
      }
    } catch (const Erro&) {}
  {
    const std::vector<std::string> nm = nomes_theta_mt(d.modelo, d.alvos);
    std::string quais;
    for (std::size_t k = 0; k < sem_par.size(); k++)
      if (sem_par[k]) quais += (quais.empty() ? "" : ", ") + nm[k];
    if (!quais.empty())
      R.mensagem += std::string(R.mensagem.empty() ? "" : "; ") + "no record has both "
          "traits of " + quais + ", so the data do not estimate it. The fit holds it at 0, "
          "and it has no standard error";
  }

  // o mesmo veredicto do model() (aireml.cpp): esgotado o maxiter com o decremento de
  // Newton ja abaixo da tolerancia, o que se tem e uma crista plana, e nao um ajuste
  // que pede mais iteracoes. Aqui ele ainda dizia "Raise maxiter". So quando o laco andou
  // ate o fim (com maxiter = 0 uma avaliacao no start= nao e convergencia) e a sonda de KKT
  // nao achou ponto melhor dentro do cone.
  if (crista && k_crista == Kkt::limite) {
    char buf[320];
    std::snprintf(buf, sizeof(buf), " (relDelta %.3g, Newton decrement %.3g): a component "
                  "resting at a covariance boundary is not at a constrained optimum, a point "
                  "inside the cone or along the boundary face lowers the -2logL by %.3g", R.reldelta, R.decremento,
                  ganho_kkt);
    R.mensagem = "stopped at " + std::to_string(maxiter) + " iteration(s) without a "
        "certified optimum" + buf + ". Raise maxiter=" +
        (R.mensagem.empty() ? std::string() : "; " + R.mensagem);
  }
  if (crista && k_crista == Kkt::otimo) {
    R.convergiu = true;
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.3g", R.reldelta);
    R.mensagem = std::string(R.mensagem.empty() ? "" : R.mensagem + "; ") +
        "converged by the Newton decrement while the components were still moving "
        "(relDelta " + buf + "): the likelihood is flat along some direction, so the "
        "values along it are about equally supported. Read the standard errors and the "
        "correlation between the components before trusting any one of them";
  }

  // o diagnostico de nao convergencia sai sempre que nenhum outro motivo de falha foi
  // escrito, e vai NA FRENTE das notas informativas (exclusao por lista, fronteira, tempo
  // ausente). Antes so saia com a mensagem vazia, e qualquer nota o escondia.
  if (!R.convergiu && maxiter == 0) {
    // uma avaliacao, nao um ajuste: dito como tal
    char buf[160];
    std::snprintf(buf, sizeof(buf), "maxiter = 0: no iteration was run, so the components "
                  "are the start evaluated, not estimated; the Newton decrement there is "
                  "%.3g", R.decremento);
    R.mensagem = std::string(buf) + (R.mensagem.empty() ? std::string() : "; " + R.mensagem);
  } else if (!R.convergiu && R.iters >= maxiter &&
             R.mensagem.find("did NOT converge") == std::string::npos &&
             R.mensagem.find("stopped at ") != 0) {
    // o teto de verdade: so aqui o ajuste parou POR causa do maxiter
    char buf[128];
    std::snprintf(buf, sizeof(buf), " (relDelta %.3g, Newton decrement %.3g against the "
                  "2e-4 tolerance)", R.reldelta, R.decremento);
    // relDelta 0: a caminhada CONGELOU (nenhum passo mudou os componentes), e o mesmo ponto
    // volta com qualquer maxiter (lista de correcoes, item 8). O conselho de subir o maxiter
    // so vale quando ela ainda anda
    const std::string conselho = R.reldelta == 0.0
        ? ": the walk froze, no step changed the components, so more iterations give the same "
          "point. Fit from another start=, or read the components resting at a boundary"
        : ": this model asks for more iterations. Raise maxiter=";
    R.mensagem = "stopped at " + std::to_string(maxiter) + " iteration(s) without a "
        "certified optimum" + buf + conselho +
        (R.mensagem.empty() ? std::string() : "; " + R.mensagem);
  }
  if (cur.fora_do_padrao > 0)
    R.mensagem += std::string(R.mensagem.empty() ? "" : "; ") +
        std::to_string(cur.fora_do_padrao) + " value(s) read outside the sparsity pattern of the factor: the result is NOT reliable";

  if (verboso)
    Rprintf("%s at iter %d, -2logL %.6f\n", R.convergiu ? "converged" : "STOPPED",
            (int) R.iters, final_ok ? cur.neg2logl : cur.neg2logl + jac);
  R.theta = theta;
  R.fora_do_padrao = cur.fora_do_padrao;
  R.na_parede = parede_final;
  R.parou = saiu_do_laco && !R.convergiu;
  if (final_ok) {
    R.neg2logl = cur.neg2logl;
    R.solucao = cur.solucao;
    R.pev = cur.pev;
    // o veredito de singularidade so sobre as coordenadas da verossimilhanca: a covariancia
    // residual sem par ja foi dita acima, e com ela a AI e singular por construcao
    const std::vector<std::string> nm = nomes_theta_mt(d.modelo, d.alvos);
    std::vector<std::string> nm_v;
    std::vector<double> se_v, fc_v;
    for (std::size_t k : na_veross) {
      nm_v.push_back(nm[k]);
      se_v.push_back(R.se[k]);
      fc_v.push_back(k < fator_corte.size() ? fator_corte[k] : 1.0);
    }
    const std::string av = aviso_informacao(sub_densa(cur.ai, na_veross), nm_v, &se_v, &fc_v);
    for (std::size_t i = 0; i < na_veross.size(); i++) R.se[na_veross[i]] = se_v[i];
    if (!av.empty()) R.mensagem += (R.mensagem.empty() ? "" : "; ") + av;
    // o componente sem erro-padrao (numa crista do veredito SINGULAR) fica sem linha na
    // covariancia: o metodo delta de uma funcao que depende dele tambem nao existe; o dado
    // na parede fica com a linha 0
    for (std::size_t k = 0; k < ntv; k++)
      if (!dado[k] && !std::isfinite(R.se[k]))
        for (std::size_t j = 0; j < ntv; j++) {
          if (dado[j]) continue;
          R.vcov[k * ntv + j] = std::nan("");
          R.vcov[j * ntv + k] = std::nan("");
        }
    std::string dados;
    for (std::size_t k = 0; k < ntv; k++)
      if (dado[k]) dados += (dados.empty() ? "" : ", ") + nm[k];
    if (!dados.empty())
      R.mensagem += (R.mensagem.empty() ? "" : "; ") + std::string("held at the boundary, ") +
          dados + " has no standard error, and the standard errors of the other components "
          "are conditional on it";
  } else {
    // nunca visto: o mesmo theta, so reescalado, deixou de ser avaliavel. O -2logL sai da
    // caminhada mais o jacobiano, e as solucoes ficam NA em vez de sair na unidade errada
    R.neg2logl = cur.neg2logl + jac;
    R.mensagem += std::string(R.mensagem.empty() ? "" : "; ") +
        "the final evaluation in the original units failed, so the solutions are not "
        "reported";
  }
  return R;
}

AjusteMT ajusta_mt(DesenhoMT& d, const std::vector<double>* theta0, std::size_t maxiter,
                   double tol, bool verboso) {
  AjusteMT R = caminha_mt(d, theta0, maxiter, tol, verboso, 0);
  MapaZ mz;
  for (const Grupo& gr : d.modelo.grupos) mz.blocos.push_back(std::make_pair(gr.offset, gr.dim));
  mz.blocos.push_back(std::make_pair(d.modelo.offset_residual, d.t));
  // a partida de var(y) em unidades ORIGINAIS: y ja voltou da padronizacao
  return explora_fronteira(std::move(R), mz, partida_mt(d), maxiter,
                           [&](const std::vector<double>& th0, std::size_t it0) {
                             return caminha_mt(d, &th0, maxiter, tol, verboso, it0);
                           }, verboso);
}

}  // namespace br
