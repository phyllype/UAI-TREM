// Declaracoes compartilhadas entre a montagem (mme.cpp), o estimador (aireml.cpp) e a
// travessia (entrada.cpp).
//
// Um unico lugar, de proposito: a primeira versao redefinia `Desenho` em dois .cpp, o que e
// violacao de ODR, os dois corpos tem de ser identicos token a token, e qualquer edicao
// futura num deles quebraria o outro em silencio, no vinculo ou pior, em tempo de execucao.

#ifndef UAITREM_MME_H
#define UAITREM_MME_H

#include <unordered_set>
#include "modelo.h"

#include <functional>
#include <limits>
#include <memory>

namespace br {

// ---- pedigree.cpp
struct Pedigree {
  std::vector<std::string> ids;
  std::vector<std::int64_t> pai, mae;
  // metafundadores (Legarra et al. 2015). eh_mf marca as linhas-base virtuais, prefixadas
  // ao pedigree, e gama traz a DIAGONAL de Gamma por linha. O truque que faz a variancia
  // mendeliana generalizar sozinha: F(mf) = gamma_ii - 1, NEGATIVO quando gamma < 1 e
  // legitimamente negativo, entao nao pode ser grampeado em zero em lugar nenhum.
  //
  // ATENCAO ao pseudo-codigo IMPRESSO no artigo, que traz "F(i) = 1 - gamma(i)": o sinal
  // esta trocado la, e com ele o proprio Exemplo 1 do artigo nao fecha. O texto corrido diz
  // o certo: "a self-relationship of a11 = gamma and an individual inbreeding coefficient
  // of Fi = a11 - 1 = gamma - 1".
  //
  // Gamma pode ser CHEIA. O gamma_jk fora da diagonal e a ancestralidade compartilhada
  // entre duas bases, que e o parametro que interessa em analise multirraca, e ele entra em
  // TRES lugares, dos quais so o primeiro e obvio:
  //   1. o bloco Gamma^-1 de A^-1 (gama_inv), SOMADO por cima das contribuicoes que os
  //      filhos ja jogam nas linhas dos metafundadores;
  //   2. a endogamia, pelo termo l_mf' Gamma l_mf, computado como ||K' l_mf||^2 com K a
  //      Cholesky inferior (gama_chol), que e a forma publicada e a numericamente estavel;
  //   3. INDIRETAMENTE na variancia mendeliana, sempre que um pai for mestico, porque o F
  //      desse pai ja carrega gamma_jk dentro.
  // Com Gamma diagonal os tres colapsam no comportamento anterior por construcao. Nao ha
  // ramo separado para o caso diagonal de proposito: um ramo duplicado e onde os dois
  // comportamentos divergiriam com o tempo.
  //
  // col_mf existe porque a ordenacao topologica REORDENA as linhas: o metafundador k de
  // Gamma nao fica necessariamente na linha k do pedigree ordenado. col_mf[i] devolve a
  // coluna de Gamma da linha i, ou -1. Sem isso, um desalinhamento com Gamma diagonal so
  // trocaria gamma_ii de lugar; com Gamma cheia ele corrompe todos os cruzados de um jeito
  // dificil de perceber.
  std::vector<char> eh_mf;
  std::vector<double> gama;
  std::vector<std::int64_t> col_mf;
  std::size_t n_mf = 0;
  // K com K K' = Gamma, de decomposicao ESPECTRAL e portanto NAO triangular: e o que
  // permite Gamma singular (gamma = 0, o limite de grupo de pais desconhecidos, e dois
  // metafundadores para a mesma populacao). Quem soma ||K' l||^2 tem de varrer a coluna
  // INTEIRA, nao so o triangulo.
  std::vector<double> gama_chol;   // n_mf x n_mf, por linhas
  std::vector<double> gama_inv;    // Gamma^-1 (pseudo-inversa se singular), por linhas
  // PEDIGREE PAI / AVO MATERNO (Henderson 1975, 1976; Mrode e Pocrnic, secs. 3.6-3.7). Com
  // mgs[i] = 1 a SEGUNDA coluna da linha i (mae) e o avo materno: o caminho dele pesa 1/4 e
  // nao 1/2, e a variancia mendeliana e 11/16 - F_s/4 - F_k/16. A flag e POR LINHA e segue
  // a ordenacao topologica junto com a linha (o mesmo cuidado do col_mf: uma flag fora do
  // lugar da uma A^-1 simetrica, positiva-definida e errada). Vazio = tudo pai/mae.
  std::vector<char> mgs;
};
Pedigree constroi_pedigree(const std::vector<std::string>&, const std::vector<std::string>&,
                           const std::vector<std::string>&,
                           const std::vector<std::string>& = {},
                           const std::vector<double>& = {},
                           const std::vector<char>& = {});
// variancia mendeliana da linha i com os F dos pais ja conhecidos (pai/mae ou pai/MGS)
double variancia_mendeliana(const Pedigree&, const std::vector<double>& f, std::size_t i);
std::vector<double> endogamia(const Pedigree&);
// estado das ESTIMATIVAS no verbose, compartilhado pelos tres ajustadores iterativos
void imprime_theta(const std::vector<double>&, const std::vector<std::string>&,
                   const Modelo&);
Csc a_inversa(const Pedigree&, const std::vector<double>&);

// ---- dominancia.cpp
// A precisao da dominancia aumentada por subclasses pai x mae (Hoeschele e VanRaden 1991, na
// forma geral): niveis 0..n_q-1 sao os animais de Q, n_q.. os pares (subclasses cheias e, na
// rota esparsa, os pares do fecho que sobreviveram a eliminacao). Ver dominancia.cpp.
struct DominanciaInversa {
  std::vector<std::size_t> animal;          // linha no pedigree de cada animal de Q
  std::vector<std::int64_t> par_x, par_y;   // linhas do pedigree de cada par, x <= y
  std::vector<char> par_cheio;              // e a subclasse de algum animal de Q
  std::vector<double> priori_h;             // F_cc / 4, a variancia a priori do par
  Csc q;                                    // triangulo inferior de Q
  double logdet_k = 0.0;                    // log|K_aug| pela forma fechada
  double logdet_qhh = 0.0;                  // log|Q_hh|
  double logdet_d = 0.0;                    // log|D| = log|K_aug| + log|Q_hh|
  int rota = 0;                             // 1 densa, 2 esparsa
  std::size_t n_sub = 0, n_criados = 0, n_processados = 0;
  bool fecho_abortado = false;
  // soma de c_j^2 do fator das MME de referencia de cada rota, e o custo que a rota automatica
  // compara (a parte esparsa pesando mais que a cauda densa; ver dominancia.cpp)
  double flops_densa = std::numeric_limits<double>::quiet_NaN();
  double flops_esparsa = std::numeric_limits<double>::quiet_NaN();
  double custo_densa = std::numeric_limits<double>::quiet_NaN();
  double custo_esparsa = std::numeric_limits<double>::quiet_NaN();
  // segundos de relogio: fecho de pares, decisao da rota automatica e montagem de Q
  double seg_fecho = 0.0, seg_decisao = 0.0, seg_montagem = 0.0;
  // a decisao da rota automatica: o motivo (ver dominancia.cpp, Motivo), o trabalho da
  // ordenacao da candidata esparsa e o teto dele, o pico estimado de montar Q pela densa e os
  // triplos que cada rota monta no bloco de pares de Q (o que a saida M_MEMORIA compara)
  int motivo = 0;
  double trabalho_ordem = std::numeric_limits<double>::quiet_NaN();
  double trabalho_densa = std::numeric_limits<double>::quiet_NaN();
  double teto_trabalho = std::numeric_limits<double>::quiet_NaN();
  double bytes_densa = std::numeric_limits<double>::quiet_NaN();
  double entradas_densa = std::numeric_limits<double>::quiet_NaN();
  double entradas_esparsa = std::numeric_limits<double>::quiet_NaN();
};
// em_q: um por animal do pedigree; rota 0 automatica, 1 densa, 2 esparsa
DominanciaInversa dominancia_inversa(const Pedigree&, const std::vector<char>& em_q, int rota,
                                     double max_pares);

// ---- linalg.cpp
// interrupcao do usuario como excecao C++, so fora de regiao paralela. NAO deriva de Erro de
// proposito: varios lacos tratam um Erro de inversao como theta inadmissivel e seguem
// (catch (const Erro&)), e ali a interrupcao seria engolida. A casca (GUARDA) pega toda
// std::exception e a devolve ao R como erro.
struct Interrompido : std::runtime_error {
  Interrompido() : std::runtime_error("interrupted (Ctrl-C or a time limit)") {}
};
void checa_interrupcao();
struct Simbolica {
  std::vector<std::int64_t> pai;
  std::vector<std::size_t> colptr;
  std::size_t n = 0;
  std::size_t nnz() const { return colptr.empty() ? 0 : colptr.back(); }
};
Csc triu(const Csc&);
Simbolica simbolica(const Csc&);
bool cholesky(const Csc&, const Simbolica&, Csc&);
// numero de threads das regioes paralelas (OpenMP). Todo codigo paralelo do pacote segue a
// regra de DONO UNICO: cada numero de saida e escrito por uma thread so, com as somas em
// ordem fixa, entao o resultado e bit a bit o mesmo com 1 ou com 16 threads.
int threads();
void define_threads(int n);
int threads_disponiveis();
// inversas densas grandes pelo LAPACK do R (true) ou pelos ladrilhos proprios (false, padrao)
bool denso_lapack();
void define_denso_lapack(bool v);
// Cholesky em ladrilhos de uma matriz densa guardada como triangulo inferior compactado por
// coluna (coluna j = linhas j..n-1, em sequencia). false se nao for positiva-definida.
bool cholesky_empacotada(std::vector<double>& v, std::size_t n, int nth);
// c[i*ldc + j] += soma_k x[k*ldx + i] y[k*ldy + j], i < m, j < n (x e y "por k": para cada k
// os i sao contiguos). simetrico (x == y, m == n): so o triangulo j <= i e atualizado, e
// quem chama espelha. Ladrilhos 64 x 64 com dono unico e soma em ordem fixa.
void produto_ladrilhos(const double* x, std::size_t ldx, const double* y, std::size_t ldy,
                       std::size_t K, std::size_t m, std::size_t n, double* c,
                       std::size_t ldc, bool simetrico, int nth,
                       bool x_por_linha = false);

// Cache da fatoracao simbolica para os lacos de ajuste: o PADRAO de C nao depende de
// theta (a montagem empurra sempre os mesmos slots, zeros explicitos inclusive), entao
// o grau minimo e a analise simbolica valem para o ajuste inteiro e sao refeitos hoje a
// cada avaliacao por nada. Ponteiro nulo preserva o comportamento avulso dos gates.
struct CacheSimbolica {
  bool pronto = false;
  std::vector<std::size_t> perm;
  Simbolica sb;
  // O MAPA de valores de C para a C permutada. Refazer a permutacao por triplos a cada
  // avaliacao custa uma ordenacao por contagem com acesso espalhado, da ordem da propria
  // fatoracao da cauda quando ela roda em paralelo; com o mapa, cada avaliacao so espalha os
  // valores. O padrao e conferido por uma assinatura a cada uso e o mapa e refeito se mudar.
  std::uint64_t assinatura = 0;
  int vistos = 0;               // o mapa so e montado na SEGUNDA avaliacao do mesmo padrao
  bool com_mapa = false;
  std::vector<std::uint32_t> mapa;
  Csc permutada;
  // a montagem das MME de um desenho (monta_mme): W'W com zeros explicitos nas posicoes da
  // penalidade e a posicao de cada entrada da penalidade no padrao final (W'y nao: y muda
  // entre chamadas na cadeia probit)
  const void* mont_dono = nullptr;
  Csc mont_base;
  std::vector<std::size_t> mont_pos;
};
double logdet(const Csc&);
std::vector<double> resolve(const Csc&, const std::vector<double>&);
Csc permuta_sim(const Csc&, const std::vector<std::size_t>&);
// a mesma permutacao, pelo mapa guardado em cs (cs.perm ja definido)
const Csc& permuta_cache(const Csc&, CacheSimbolica& cs);
// assinatura (FNV-1a) do padrao de uma CSC: dimensoes, colptr e linhas
std::uint64_t assinatura_padrao(const Csc&);
// teto > 0: desiste (devolve vazio) quando a soma parcial de c_j^2 da ordem passa dele;
// teto_trabalho > 0: desiste quando as entradas de lista lidas passam dele. trabalho, se dado,
// recebe as entradas lidas (a medida deterministica do custo da ordenacao)
std::vector<std::size_t> grau_minimo(const Csc&, double teto = 0.0, double teto_trabalho = 0.0,
                                     double* trabalho = nullptr);
// a ordenacao por grau minimo EXATO (linalg.cpp) e a AMD (ordem.cpp); grau_minimo() escolhe pelo
// UAITREM_ORDEM do ambiente ("gm" ou "amd"); com teto ou trabalho pedido, sempre o exato
std::vector<std::size_t> grau_minimo_exato(const Csc&, double teto = 0.0, double teto_trabalho = 0.0,
                                           double* trabalho = nullptr);
std::vector<std::size_t> ordem_amd(const Csc&);
int modo_ordem();
// PARTICAO DA ARVORE DE ELIMINACAO em subarvores disjuntas para as threads (linalg.cpp). Os nos
// >= limite (a cauda densa) e os que foram divididos ficam no TOPO, feito depois (fatoracao) ou
// antes (inversa seletiva) das subarvores. Vazia (ini.size() < 2) quando nao compensa: quem chama
// segue no caminho serial. As subarvores vem da mais pesada para a mais leve; dentro de cada uma,
// os nos em ordem crescente.
struct Particao {
  std::vector<std::uint32_t> topo;
  std::vector<std::size_t> ini;
  std::vector<std::uint32_t> nos;
  bool vazia() const { return ini.size() < 2; }
};
Particao particao_arvore(const std::vector<std::int64_t>& pai, const std::vector<double>& peso,
                         std::size_t limite, int nth);

// ---- selinv.cpp
struct SelInv {
  std::size_t n = 0;
  std::vector<std::size_t> colptr;
  std::vector<std::uint32_t> linha;
  std::vector<double> valor;
  bool get(std::size_t, std::size_t, double&) const;
  std::vector<double> diagonal() const;
};
std::size_t bloco_denso_final(const Csc&);
// o mesmo numero lido do padrao simbolico, antes de fatorar: quantas colunas finais de L
// sao completamente cheias. E o k do custo k^3 de cada fatoracao.
std::size_t bloco_denso_simbolico(const Simbolica&);
// aviso de matriz de informacao singular, com os componentes que o dado nao separa. `fator`
// (opcional) e a unidade de cada componente (s_a s_b com s o desvio padrao EXATO de cada
// caracteristica; 1 para o rho): o corte e feito na AI dessas unidades, F AI F (aireml.cpp)
std::string aviso_informacao(const Densa& ai, const std::vector<std::string>& nomes,
                             std::vector<double>* se = nullptr,
                             const std::vector<double>* fator = nullptr);
SelInv inversa_seletiva(const Csc&, std::size_t);

// Uma matriz de covariancia DECLARADA pelo usuario (kernel(id, K=)): os ids que nomeiam
// linhas e colunas, e a K densa. O vetor de declaracoes anda PARALELO aos termos do
// modelo, entrada vazia para termo que nao e kernel, para que o desenho ache a K do
// termo sem um mapa a parte.
struct KernelDecl {
  std::vector<std::string> ids;
  Densa k;
  // a PRECISAO declarada (kernel(id, Kinv =)): o triangulo inferior de K^-1, ja esparso. O
  // caso que pede isto e a dominancia por subclasses (dominance_inverse()), em que a K
  // densa nunca existe; a rota densa de cima fica intocada.
  bool esparsa = false;
  Csc kinv;
  // escala de PARTIDA dos espelhos e do gibbs: a media geometrica dos autovalores da matriz
  // que o usuario passaria como K. Na precisao aumentada da dominancia ela e a de D, nao a
  // de K_aug (que tem os niveis de subclasse), e e o que faz kernel(Kinv =) partir do mesmo
  // ponto de kernel(K = D). 0 = tirar do log|K^-1| como sempre.
  double escala = 0.0;
  // log|K^-1| dado pronto (o de dominance_inverse(), pela forma fechada), so aceito quando a
  // assinatura dos triplos confere: poupa a ordenacao e a Cholesky de kinv_declarada. NaN =
  // calcular pela Cholesky esparsa.
  double logdet_kinv = std::numeric_limits<double>::quiet_NaN();
  bool vazia() const { return ids.empty(); }
};

// a escala de partida de um grupo declarado: a do KernelDecl quando veio, senao a media
// geometrica dos autovalores de K, exp(-log|K^-1| / n), tirada do log-determinante guardado
inline double escala_kernel(const std::vector<double>& esc, const std::vector<double>& ld,
                            const std::vector<Csc>& kinv, std::size_t g) {
  if (g < esc.size() && esc[g] > 0.0) return esc[g];
  return std::exp(-ld[g] / static_cast<double>(kinv[g].ncol));
}

// ---- o desenho completo
struct Desenho {
  Modelo modelo;
  // X dos efeitos fixos que ficaram, por linhas, em toda linha da tabela (as nao usadas com o
  // valor sem peso): os mesmos nao-zeros do X denso de antes (fixos.cpp)
  Csr x;
  std::vector<std::string> nomes_x;
  std::vector<std::string> saiu_x;
  // a padronizacao das covariaveis (padroniza_covariaveis, fixos.cpp): x guarda as colunas
  // padronizadas e b sai do motor nessa escala; volta_x leva b a unidade do dado (uma por
  // caracteristica, aqui uma so) e logdet_escala_x devolve o -2logL a escala do dado
  std::vector<VoltaB> volta_x;
  double logdet_escala_x = 0.0;
  std::vector<DesenhoTermo> aleatorios;
  std::vector<Csc> kinv;
  std::vector<double> kinv_logdet;
  std::vector<double> kinv_escala;   // escala de partida por grupo (0 = do log-determinante)
  std::vector<double> y;
  std::vector<char> usa;
  std::size_t nlin = 0;
  // a classe de residuo (0..K-1) de cada linha da tabela, com Modelo::classes_residuo; vazio
  // com um residuo so
  std::vector<std::uint32_t> classe_res;
  // Pesos: um registro de peso w tem residual s2e/w. Entram como ESCALA DE LINHA por
  // sqrt(w) em y, X e Z, as equacoes normais da tabela escalada SAO as equacoes
  // ponderadas, entao score, AI e inversa seletiva continuam valendo sem uma linha de
  // mudanca. So a verossimilhanca precisa do jacobiano, que e esta constante.
  double logdet_peso = 0.0;

  std::size_t n_usadas() const {
    std::size_t s = 0;
    for (char u : usa) s += u;
    return s;
  }
  std::size_t largura(std::size_t g) const {
    std::size_t s = 0;
    for (std::size_t t : modelo.grupos[g].termos)
      for (const DesenhoTermo& d : aleatorios)
        if (d.termo == t) s += d.z.ncol;
    return s;
  }
  std::size_t total_colunas() const {
    std::size_t s = x.ncol;
    for (std::size_t g = 0; g < modelo.grupos.size(); g++) s += largura(g);
    return s;
  }
};

// ---- mme.cpp
// O SISTEMA TRANSFORMADO de um grupo perto da singularidade (lista de correcoes, itens 3 e 41).
// Perto de uma correlacao +-1 a soma log|C_g (x) K| + log|C| das MME cancela dois numeros
// enormes: C carrega C_g^-1, e a eliminacao de Cholesky subtrai termos da ordem de
// 1/lambda_min(C_g). O erro cresce como 1/(1 - |r|) e com o numero de niveis. Medido contra a
// forma V densa com 1 - |r| = 1e-12: 0,002 a 0,03 unidades no model() (480 registros), 0,001
// a 0,003 no model_mt() (280) e 0,0007 a 0,014 no model_ar1() (840). No dado real foram 1470
// unidades num salto que a caminhada aceitou, porque o valor errado era o menor.
// Com C_g = F F' (F pela decomposicao espectral), u = (F (x) I) u* e Var(u*) = I (x) K. As MME
// de u* tem T' C_d T + I (x) K^-1 no bloco do grupo, com C_d a parte dos dados e
// T = diag(I, F (x) I), sem C_g^-1; e log|C_g (x) K| + log|C| = -dim log|K^-1| + log|C*| pela
// identidade de Sylvester. O valor e o mesmo -2logL, e fica certo ate com C_g singular. As
// colunas do grupo seguem a convencao das tres montagens: coordenada a, nivel l na coluna
// off + a * nl + l.
struct BlocoT {
  std::size_t off = 0, dim = 0, nl = 0;
  Densa F;
};
// F de C = F F' pela decomposicao espectral, o autovalor negativo do arredondamento em zero
Densa fator_espectral(const Densa& c);
// F[g] do grupo cuja matriz de correlacao tem autovalor minimo abaixo de 1e-6, vazio nos
// outros; `cgs` e a covariancia de cada grupo na escala da montagem. Devolve se algum entrou
bool fatores_quase_singulares(const std::vector<Densa>& cgs, std::vector<Densa>& F);
// T' C_d T nos triplos (triangulo superior; os duplicados a montagem soma) e T' rhs
void transforma_triplos(const std::vector<BlocoT>& blocos, std::vector<std::uint32_t>& ti,
                        std::vector<std::uint32_t>& tj, std::vector<double>& tv,
                        std::vector<double>& rhs);
// u'(C_g^-1 (x) K^-1)u de todos os grupos, em unidades absolutas (o que avalia_mt e avalia_ar1
// chamam upen), com u*'(I (x) K^-1)u* nos grupos de `F` nao vazio: `sol` e a solucao do
// sistema transformado, onde Var(u*) = escala_t (I (x) K) (s2e no model(), 1 nos outros).
// False quando uma C_g nao se inverte
bool upen_transformada(const Modelo& m, const std::vector<DesenhoTermo>& aleatorios,
                       const std::vector<Csc>& kinv, const std::vector<double>& theta,
                       const std::vector<std::size_t>& offset_grupo, const std::vector<Densa>& F,
                       const std::vector<double>& sol, double escala_t, double& upen);
// u = (F (x) I) u* no bloco de cada grupo transformado
void destransforma_solucao(const std::vector<BlocoT>& blocos, std::vector<double>& sol);
// u' K^-1 v de vetores de nl niveis: K^-1 no triangulo inferior, ou K = I com ncol 0
double bilinear_kinv(const Csc& kinv, const double* u, const double* v, std::size_t nl);
inline double forma_kinv(const Csc& kinv, const double* u, std::size_t nl) {
  return bilinear_kinv(kinv, u, u, nl);
}

struct Montado {
  Csc c;
  std::vector<double> rhs;
  double logdet_g = 0.0;
  double s2e = 0.0;
  std::size_t n_fixo = 0;
  std::size_t total = 0;
  std::vector<std::size_t> offset_grupo;
  std::vector<BlocoT> blocos_t;   // os grupos do sistema transformado (fatores = nao nulo)
  bool ok = false;
};
// a K declarada de um grupo, compartilhada pelos tres ajustadores
void reduz_kernels(const Modelo& m, const std::vector<KernelDecl>*& kernels,
                   std::vector<KernelDecl>& kern_red,
                   std::vector<std::unordered_set<std::string> >& kern_nulos);
void casa_niveis_nulos(const Modelo& m, const std::vector<DesenhoTermo*>& aleatorios,
                       const std::vector<std::unordered_set<std::string> >& kern_nulos,
                       const Tabela& t, std::size_t nlin);
// os termos aleatorios dos tres montadores, com os niveis de cada grupo: o pedigree, os ids
// da K declarada (ja reduzida), ou a uniao dos rotulos num grupo sem estrutura de varios
// termos. Recusa grupo cujos termos nao indexam o mesmo conjunto de niveis.
std::vector<DesenhoTermo> monta_aleatorios(const Modelo& m, const Tabela& t,
                                           const std::vector<std::string>& niveis_ped,
                                           const std::vector<KernelDecl>* kernels);
// A mensagem de quando nenhuma linha entra, a mesma nos tres montadores: quantas linhas ha,
// quantas sairam pela observacao (`sem_obs`, com `oque_obs` dizendo o que faltou) e, por
// termo, quantas sairam por um nivel fora do conjunto, com o primeiro rotulo dos dados e o
// primeiro nivel do conjunto. Antes era "no row enters the analysis" para todas as causas.
std::string sem_linhas(const Modelo& m, const Tabela& t,
                       const std::vector<DesenhoTermo>& aleatorios,
                       const std::vector<char>& sem_obs, const std::string& oque_obs);
// Antes de ESTIMAR: recusa o grupo sem estrutura em que um par de termos nao tem nenhum nivel
// com registro (linha usada) nos dois, porque a covariancia do par nao entra na
// verossimilhanca. Os ajustadores a chamam so quando os componentes andam.
void confere_covariancias_iid(const Modelo& m, const std::vector<DesenhoTermo>& aleatorios,
                              const std::vector<char>& usa);
void kinv_declarada(const Modelo&, const Grupo&, const std::vector<KernelDecl>*,
                    std::vector<Csc>&, std::vector<double>&, std::vector<double>&);
// `fatores`, por grupo (vazio = o grupo como esta), monta o sistema transformado; vale so
// sem cache, porque o padrao muda
Montado monta_mme(const Desenho&, const std::vector<double>&, CacheSimbolica* = nullptr,
                  const std::vector<Densa>* fatores = nullptr);
Densa cov_grupo(const Modelo&, const std::vector<double>&, std::size_t);

// ---- aireml.cpp: Cholesky of a small dense covariance group; false when not PD.
// Coordenadas do PASSO, log-Cholesky por bloco (Pinheiro e Bates, 1996). O racional
// completo, com a falha medida que obrigou a troca, esta em src/aireml.cpp. O mapa e o que
// permite aos tres ajustadores usarem a MESMA maquinaria: o univariado tem um residuo
// escalar, o multicaracter tem um R0 de t x t e o AR(1) tem R0 mais o rho, que anda em
// atanh para que |rho| = 1 fique no infinito em vez de ser uma parede.
struct MapaZ {
  std::vector<std::pair<std::size_t, std::size_t> > blocos;  // (offset, dim) no triangulo
  std::vector<std::size_t> escalares_log;                    // variancia isolada, em log
  std::vector<std::size_t> correlacoes;                      // em atanh, |x| < 1
};
bool z_de_theta(const MapaZ& mp, std::size_t ntheta, const std::vector<double>& theta,
                std::vector<double>& z);
std::vector<double> theta_de_z(const MapaZ& mp, std::size_t ntheta,
                               const std::vector<double>& z);
Densa jacobiano_z(const MapaZ& mp, std::size_t ntheta, const std::vector<double>& z);
// A covariancia das componentes CONDICIONAL a uma parede de bloco: V = 2 J_F (J_F' AI J_F)^-1
// J_F', com J = d theta / d z e F as coordenadas de z livres; escreve vcov (ntheta x ntheta,
// por linha) e se nas linhas `livre` e devolve false se J_F' AI J_F nao inverte. Os tres
// ajustadores chamam quando um bloco de varias dimensoes para na parede (aireml.cpp)
bool vcov_condicional(const Densa& ai, const Densa& J, const std::vector<std::size_t>& fz,
                      const std::vector<std::size_t>& livre, std::vector<double>& vcov,
                      std::vector<double>& se);
// Todas as pecas do passo num ponto: z, score e AI levados para z, os pisos, e os DOIS
// conjuntos ativos sobre os mesmos pisos e a mesma regra de score para fora, o do passo
// (piso + 1e-9) e o do certificado, deliberadamente mais largo (piso + log 2). Ver o
// racional em src/aireml.cpp.
struct PassoZ {
  bool ok = false;                 // false quando o theta corrente esta fora do cone
  std::vector<double> zc, sz, piso;
  Densa az;
  std::vector<char> congelado, na_parede;
};
PassoZ passo_z(const MapaZ& mp, std::size_t ntheta, const std::vector<double>& theta,
               const std::vector<double>& score, const Densa& ai, double escala_residual);
double decremento_z(const PassoZ& P);
// os nomes dos componentes na parede (lista de correcoes, item 21): ate 8, e quantos sobram
std::string lista_na_parede(const std::vector<char>& parede, const std::vector<std::string>& nomes);

void pisos_z(const MapaZ& mp, const std::vector<double>& z, double escala_residual,
             std::vector<double>& piso, std::vector<char>& relativo);
// A SONDA DE KKT pela verossimilhanca (o racional esta em aireml.cpp): pontos dentro do
// cone para cada bloco com coordenada na parede. `direcionais` seguem o autovetor de
// autovalor negativo da matriz do gradiente do bloco (a direcao C + e ww' que baixa o
// -2logL em primeira ordem), com o comprimento do modelo quadratico da AI; `coordenadas`
// levantam cada diagonal de Cholesky na parede, e o bloco inteiro descorrelacionado, a
// fracoes da escala. Tudo em z. Quem chama avalia, primeiro as direcionais, e so aceita o
// que baixa o -2logL por mais que tol_reentrada.
struct CandidatosKKT {
  std::vector<std::vector<double> > direcionais, coordenadas;
};
CandidatosKKT candidatos_reentrada(const MapaZ& mp, std::size_t ntheta,
                                   const std::vector<double>& z,
                                   const std::vector<char>& parede, double escala,
                                   const std::vector<double>& score, const Densa& ai);
// A SONDA DA FACE (o racional esta em aireml.cpp): diferencas finitas do -2logL nas
// coordenadas dos blocos na parede que andam sobre a face, o score do resto nos mesmos
// pontos, e o passo de Newton conjunto. `f(z, valor, gradiente_z)` avalia num z (false se
// inadmissivel); true, com o melhor ponto e o valor dele, quando algum ponto baixa o -2logL
// (por qualquer margem: quem chama compara o ganho com tol_lapida e tol_reentrada).
bool sonda_face(const MapaZ& mp, std::size_t ntheta, const PassoZ& P, double f0,
                const std::function<bool(const std::vector<double>&, double&,
                                         std::vector<double>&)>& f,
                std::vector<double>& z_novo, double& f_novo);
// o candidato descorrelacionado da lapidacao (o racional esta em aireml.cpp): as covariancias
// de cada caracteristica cuja variancia marginal no bloco esta no piso, em zero; em z e nos
// pisos de P, ou vazio quando nao ha caracteristica assim
std::vector<double> descorrelada_na_parede(const MapaZ& mp, std::size_t ntheta,
                                           const PassoZ& P);
// o veredito da sonda: nenhum ponto melhora (o otimo restrito passa), o caminhante foi
// levado para dentro do cone, ou um ponto melhora e nao houve como seguir (o limite de
// reentradas, ou o fim do orcamento de iteracoes): ai o certificado e recusado
enum class Kkt { otimo, moveu, limite };
// A FOLGA DA REENTRADA: ganho minimo, em unidades de -2logL, para um ponto de reentrada
// derrubar o certificado. Metade da tolerancia do decremento (2e-4 ~ duas vezes o gap).
constexpr double tol_reentrada = 1e-4;
// A LAPIDACAO DA FACE: depois de um certificado aceito com algum bloco na parede, o ponto
// descorrelacionado de um caracter no zero (descorrelada_na_parede) e os passos de Newton da
// sonda da face, enquanto cada um baixa o -2logL por mais que tol_lapida e o ponto novo ainda
// certifica (ate max_lapidacoes rodadas; o ponto que nao certifica fica de fora e o ajuste
// fica no ultimo que certificou). Sem isso o ponto certificado sobre uma face plana
// dependia do caminho: na revisao do item ar1real, model_mt com o bloco genetico em r = -1
// (sim_ninhadas, sementes 1 a 12, quatro partidas cada), 4 de 96 ajustes certificavam de
// 1,0e-4 a 3,9e-4 acima do melhor das quatro partidas, e um model_ar1 bivariado com uma
// variancia genetica no zero, 6,7e-4. Na semente 10 (r12 = -0,25, partida de var(y)) a sonda
// achava menos que 1e-4, abaixo da folga da reentrada, e o ajuste certificava 3,86e-4 acima;
// o segmento ate o melhor ponto desce monotonico depois de 1% do caminho. O piso de 1e-5
// fica acima do ruido do -2logL esparso perto da singularidade (3e-7 em r ~ 1, e 2,8e-6
// abaixo da forma V densa no canto mais mal condicionado medido, condicao ~1e11), para que
// a lapidacao nao persiga o arredondamento. Medido com ela, nas mesmas varreduras (368
// ajustes, quatro partidas por dado): nenhum certificado mais que 9,4e-5 acima do melhor, os
// 562 ajustes da suite identicos bit a bit, nenhum veredito trocado; o par AR(1) a 9,3e-6 (com
// o candidato descorrelacionado, descorrelada_na_parede). Sobra um caso conhecido, 1,7e-4: um
// model_mt com o bloco genetico em r ~ +1 e a diagonal condicional congelada pela curvatura
// LONGE do piso, que a sonda da face nao olha; abrir a lapidacao a esses blocos foi tentado e
// nao ajudou (1,72e-4 para 1,80e-4).
constexpr double tol_lapida = 1e-5;
constexpr int max_lapidacoes = 10;
// quantas reentradas um ajuste aceita, para que nenhum vai e volta vire laco; esgotado o
// limite a sonda continua avaliando, e um ponto que ainda melhora recusa o certificado
constexpr int max_reentradas = 10;

bool chol_pequena(const Densa&, Densa&);

struct Neg2LogL {
  bool ok = false;
  double valor = 0.0;
  std::vector<double> solucao;
  double logdet_c = 0.0;
  std::vector<std::size_t> perm;
  Csc L;
};
Neg2LogL neg2logl_esparsa(const Desenho&, const std::vector<double>&);
double neg2logl_densa_V(const Desenho&, const std::vector<double>&);

// ---- aireml.cpp
struct Avaliacao {
  bool ok = false;
  double neg2logl = 0.0;
  std::vector<double> solucao;
  std::vector<double> score;
  Densa ai;
  std::vector<double> em_theta;
  std::size_t fora_do_padrao = 0;
};
Avaliacao avalia(const Desenho&, const std::vector<double>&, CacheSimbolica* = nullptr);

struct Ajuste {
  bool convergiu = false;
  // tamanho do bloco denso final do fator, e o total de colunas: o k de k^3 por iteracao
  std::size_t bloco_denso = 0, colunas_fator = 0;
  std::size_t iters = 0;
  double reldelta = 0.0;
  // Newton decrement of the FREE components at the final point, g' AI^-1 g with the
  // active set excluded: the convergence certificate (~ twice the -2logL gap to the
  // optimum near it). NaN when the fit died before a first evaluation.
  double decremento = std::numeric_limits<double>::quiet_NaN();
  double neg2logl = 0.0;
  std::vector<double> theta, se, solucao;
  std::vector<double> pev;   // diagonal de [C_s^-1] vezes s2e, na numeracao das colunas
  std::vector<double> score; // score no ponto final: a evidencia de otimo
  std::vector<double> vcov;  // 2 AI^-1, a covariancia amostral das componentes
  std::string mensagem;
  std::size_t fora_do_padrao = 0;
  // Residuo condicional e^ = y~ - W~ s^ por linha da tabela, o ajustado condicional W~ s^,
  // o ajustado fixo X~ b^ e a alavanca h_ii = w~_i' C_u^-1 w~_i, todos na escala ESCALADA
  // por sqrt(w) (quem chama divide pelo sqrt(w) da linha) e NaN nas linhas que nao entram.
  // Saem do MESMO fator do bloco final da PEV, entao residuo, alavanca e PEV ficam no mesmo
  // theta. Var(e^_i) = s2e (1 - h_ii) / w_i. alavanca_fora conta pares (j, k) de uma linha
  // lidos fora do padrao do fator (deveria ser zero: os pares de uma linha estao em W'W).
  std::vector<double> residuo, ajustado, ajustado_fixo, alavanca;
  std::size_t alavanca_fora = 0;
};
Desenho monta_desenho(const Modelo&, const Tabela&, const Pedigree*,
                      const std::vector<double>* = nullptr,
                      const std::vector<KernelDecl>* = nullptr);

// ---- residuos.cpp
// As entradas nao nulas de cada linha usada de W~ = [X~ Z~], em CSR, nas colunas globais das
// MME do model() (X nas primeiras x.ncol, cada grupo a partir de offset_grupo, termos na
// ordem do grupo). Linha nao usada fica vazia.
struct LinhasW {
  std::vector<std::size_t> ptr;
  std::vector<std::uint32_t> col;
  std::vector<double> val;
};
LinhasW linhas_w(const Desenho& d, const std::vector<std::size_t>& offset_grupo);
// O sorteio de u*_g ~ N(0, C_g (x) K_g) por grupo: o fator de C_g (Cholesky, ou espectral na
// fronteira) e a Cholesky esparsa de K_g^-1 na ordem de grau minimo (vazia num grupo iid).
struct SorteioAleatorio {
  std::vector<std::size_t> offset, dim, nl;
  std::vector<Densa> chol_c;
  std::vector<std::vector<std::size_t>> perm_k;
  std::vector<Csc> lk;
};
SorteioAleatorio prepara_sorteio(const Modelo& mo, const std::vector<Csc>& kinv,
                                 const std::vector<DesenhoTermo>& aleatorios,
                                 const std::vector<double>& theta,
                                 const std::vector<std::size_t>& offset_grupo);
// escreve os slots de cada grupo de u (o resto de u fica intocado); usa o RNG do R
void sorteia_aleatorios(const SorteioAleatorio& S, std::vector<double>& u);
// Replicas do bootstrap parametrico com theta fixo: por replica o residuo condicional e o
// ajustado fixo das linhas usadas (coluna-major, n_usadas x nsim, na escala ESCALADA por
// sqrt(w) no model()), e, com com_u, os u* sorteados (colunas aleatorias x nsim).
struct NulosSaida {
  std::size_t n_usadas = 0, nsim = 0;
  std::vector<double> residuo, fixo, u;
};
// familia do erro das replicas: 0 normal, < 0 Laplace, > 4 t com esses graus de liberdade
double erro_padrao(double familia);
// O bootstrap PREPARADO uma vez: o desenho, o fator das MME em theta, o sorteio dos
// aleatorios e X b^ ficam vivos entre os lotes, e cada lote so sorteia e resolve. Remontar
// o desenho por lote custava, num passo unico, a G^-1 inteira a cada lote. O preparo
// tambem devolve o que confere o ajuste: o residuo dos DADOS no mesmo fator (y - W s nas
// linhas usadas, na escala escalada por sqrt(w)), que so bate com o do ajuste se o desenho,
// o pedigree e os genotipos forem os dele, e a diagonal de K_g de cada grupo (vazia num grupo
// iid), lida da inversa seletiva do fator de K_g^-1 que o sorteio ja tem: a variancia a
// priori EXATA de cada nivel, inclusive o H_ii de um animal nao genotipado num passo unico.
struct PreparoNulos {
  virtual ~PreparoNulos() {}
  std::vector<char> usa;
  // sqrt(w) de cada linha usada (vazio = 1): quem devolve ao R divide por ele
  std::vector<double> raiz;
  std::vector<double> residuo_dados;
  std::vector<std::vector<double>> diag_k;
  std::vector<std::vector<std::string>> niveis;
  // um lote de k replicas, em sequencia no RNG do R: lotes seguidos dao as mesmas replicas
  // que um lote so do mesmo tamanho total
  virtual NulosSaida lote(std::size_t k, bool com_u, double familia) = 0;
};
std::unique_ptr<PreparoNulos> prepara_nulos(Desenho&& d, const std::vector<double>& theta,
                                            const std::vector<double>& b);
// a diagonal de K = (K^-1)^-1 pela inversa seletiva de um fator ja feito na ordem perm
// (L de perm K^-1 perm'), devolvida na ordem original dos niveis
std::vector<double> diagonal_de_k(const Csc& lk, const std::vector<std::size_t>& perm);
// o mesmo fatorando K_g^-1 aqui, por grupo (vazio num grupo iid): o caminho sem theta
std::vector<std::vector<double>> diagonais_k(const std::vector<Csc>& kinv);
// os niveis de cada grupo (os do primeiro termo, que indexa a K do grupo)
template <class DES>
std::vector<std::vector<std::string>> niveis_grupos(const DES& d) {
  std::vector<std::vector<std::string>> out(d.modelo.grupos.size());
  for (std::size_t g = 0; g < d.modelo.grupos.size(); g++)
    for (const DesenhoTermo& a : d.aleatorios)
      if (a.termo == d.modelo.grupos[g].termos[0]) { out[g] = a.niveis; break; }
  return out;
}

// ---- genomica.cpp
// Os genotipos COMO O R OS GUARDA, sem copia: coluna-major (o animal corre mais rapido), em
// double, integer ou raw. Ausente e NaN no double, NA_INTEGER (o menor int) no integer e 5 no
// raw, o codigo de ausente do BLUPF90. Nada aqui imputa nem centra: quem le um marcador recebe
// os valores crus e aplica a media dele (frequencias_z). A copia double linha-major que o lado
// genomico fazia custava 8 bytes por genotipo alem dos 8 do R; com raw e 1 byte e copia nenhuma.
struct Genotipos {
  const double* d = nullptr;
  const int* i = nullptr;
  const unsigned char* r = nullptr;
  std::size_t n = 0, m = 0;
  // o trecho [a0, a1) da coluna j em out, NaN no ausente
  void trecho(std::size_t j, std::size_t a0, std::size_t a1, double* out) const {
    const std::size_t k0 = j * n;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    if (d) {
      for (std::size_t a = a0; a < a1; a++) out[a - a0] = d[k0 + a];
    } else if (i) {
      for (std::size_t a = a0; a < a1; a++) {
        const int x = i[k0 + a];
        out[a - a0] = x == std::numeric_limits<int>::min() ? nan : static_cast<double>(x);
      }
    } else {
      for (std::size_t a = a0; a < a1; a++) {
        const unsigned char x = r[k0 + a];
        out[a - a0] = x == 5 ? nan : static_cast<double>(x);
      }
    }
  }
  void coluna(std::size_t j, double* out) const { trecho(j, 0, n, out); }
  // um genotipo so, o do animal a no marcador j, NaN no ausente
  double valor(std::size_t j, std::size_t a) const {
    const std::size_t k = j * n + a;
    if (d) return d[k];
    if (i) {
      const int x = i[k];
      return x == std::numeric_limits<int>::min() ? std::numeric_limits<double>::quiet_NaN()
                                                  : static_cast<double>(x);
    }
    const unsigned char x = r[k];
    return x == 5 ? std::numeric_limits<double>::quiet_NaN() : static_cast<double>(x);
  }
};

// tr(K^-1 Z_ab) lendo a inversa seletiva Z so onde K^-1 (triangulo inferior) tem nao-zero,
// com os blocos a e b de Z nas colunas off_a e off_b. As colunas de K^-1 vao em paralelo, cada
// uma com a sua soma parcial, e as parciais sao somadas na ordem das colunas: o mesmo
// resultado com qualquer numero de threads. Um pedido fora do padrao do fator e CONTADO em
// `fora`, nunca pulado em silencio (uma soma que pula termos nao e um traco).
// A regiao so abre com 1024 colunas ou mais: uma coluna de K^-1 tem poucas entradas, e abaixo
// disso abrir o time custa mais que a soma. O traco sai a cada avaliacao para cada bloco de
// cada grupo, e num ajuste de 10 animais em 331 iteracoes 4 threads levavam 0,24 s contra
// 0,05 s com 1 (simulacao, tres rodadas). O resultado e o mesmo nos dois caminhos.
template <class Z>
double traco_kinv(const Csc& k, std::size_t nl, std::size_t off_a, std::size_t off_b,
                  const Z& z, std::size_t& fora) {
  std::vector<double> parc(nl, 0.0);
  std::vector<std::size_t> nf(nl, 0);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64) num_threads(threads()) if (nl >= 1024)
#endif
  for (long cc = 0; cc < static_cast<long>(nl); cc++) {
    const std::size_t c = static_cast<std::size_t>(cc);
    double s = 0.0;
    std::size_t f = 0;
    for (std::size_t p = k.colptr[c]; p < k.colptr[c + 1]; p++) {
      const std::size_t r = k.linha[p];
      double zv;
      if (z(off_a + r, off_b + c, zv)) s += k.valor[p] * zv;
      else f++;
      if (r != c) {
        if (z(off_a + c, off_b + r, zv)) s += k.valor[p] * zv;
        else f++;
      }
    }
    parc[c] = s;
    nf[c] = f;
  }
  double s = 0.0;
  for (std::size_t c = 0; c < nl; c++) { s += parc[c]; fora += nf[c]; }
  return s;
}

struct RelatorioG {
  std::size_t n_imputados = 0;
  std::size_t n_monomorficos = 0;
  // Diagonal de G* e a LINHA do pedigree de cada genotipado, para que accuracy() possa
  // dividir pela variancia a priori CERTA. Sob passo unico a priori de um genotipado e a
  // diagonal de H, que naquele bloco e a de G*, e nao 1 + F do pedigree: medido numa
  // populacao simulada de 510 animais todos genotipados, as duas diferem ate 0,18 e movem
  // uma acuracia individual ate 0,067. A media nao muda (0,6966 contra 0,6964); o que muda
  // e o INDIVIDUO, e com ele o ranqueamento por acuracia.
  std::vector<double> diag_gstar;
  std::vector<std::size_t> linha_ped;
};
Densa vanraden_g(const Genotipos&, RelatorioG&, bool meio = false);
// Gram do lado menor de Z (sem a escala, devolvida em `escala`), para o eigen exato
Densa gram_menor(const Genotipos&, RelatorioG&, double& escala);
// Lanczos em bloco sobre G = ZZ'/k (ou Z'Z/k, o lado menor), para a quadratura de Lanczos
// estocastica do espectro: as sondas (dim x nv, por linhas) vem do R.
struct LanczosG {
  Densa alfa, beta;                  // passos x nv
  std::vector<std::size_t> passos;   // passos de cada sonda (para antes no subespaco invariante)
  double traco = 0.0;                // tr(G), exato
  std::size_t dim = 0;
};
LanczosG lanczos_g(const Genotipos& gt, const Densa& sondas, std::size_t passos, RelatorioG& rel);
Densa a22_inversa(const Csc&, const std::vector<std::size_t>&);
// A entre os animais dados (linhas do pedigree), por Colleau: uma coluna por thread, dono unico
Densa a22_colleau(const Pedigree& p, const std::vector<std::size_t>& geno);
Csc a22_inversa_esparsa(const Csc&, const std::vector<std::size_t>&);
Densa ajusta_g_para_a22(const Densa&, const Densa&, double);
Csc constroi_hinv(const Csc&, const std::vector<std::size_t>&, const Densa&, const Densa&);

// Substitui o K^-1 dos grupos com parentesco pelo H^-1 do passo unico.
Csc h_inversa(const Pedigree& ped, const Csc& ainv, const std::vector<std::string>& geno_ids,
              const Genotipos& gt, double mistura, const std::vector<std::string>& nucleo_apy,
              std::size_t vecchia_k, RelatorioG& rel);
RelatorioG aplica_genomica(Desenho&, const Pedigree&, const std::vector<std::string>&,
                           const Genotipos&, double, const std::vector<std::string>&,
                           std::size_t = 0);
std::vector<double> partida(const Desenho&);

// ---- genomica.cpp: as contas de genotipos que eram feitas em R (F e D genomicos, Gamma dos
// metafundadores, efeitos de SNP), lendo a vista Genotipos sem copia
// a G de VanRaden no buffer de quem chama (n x n), sem a copia densa intermediaria
void vanraden_g_em(const Genotipos&, RelatorioG&, bool meio, double* out);
// heterozigotos e observados por animal
void homozigose(const Genotipos&, double* het, double* obs);
// a D de Vitezica no buffer de quem chama (n x n)
void g_dominancia(const Genotipos&, double* out);
// (X - c)' V, m x k, com o ausente valendo imp
void genotipos_tv(const Genotipos&, const double* centro, const double* imp, const double* v,
                  std::size_t k, double* out);
double conta_ausentes(const Genotipos&);
// contagens de 0, 1 e 2 por marcador nos animais marcados em linha (todos com linha nula)
void contagens_marcadores(const Genotipos&, const char* linha, double* n0, double* n1,
                          double* n2);
// as mesmas contagens so nas linhas dadas (base 0), lendo k x m genotipos
void contagens_nas_linhas(const Genotipos&, const std::vector<std::size_t>& linhas,
                          double* n0, double* n1, double* n2);
// chamadas (genotipos presentes) por animal nos marcadores marcados em usa (todos com nulo)
void chamadas_por_animal(const Genotipos&, const char* usa, double* chamadas);
// pares pai-filho de genotipados: marcadores chamados nos dois e homozigotos opostos (0 contra
// 2), um par por posicao de filho/pai (linhas em base 0), nos marcadores marcados em usa
void conflitos_parentesco(const Genotipos&, const std::vector<std::size_t>& filho,
                          const std::vector<std::size_t>& pai, const char* usa, double* ambos,
                          double* opostos);
// a linha do pedigree de cada genotipado, erro declarado para quem nao tem linha
std::vector<std::size_t> linhas_genotipados(const Pedigree&, const std::vector<std::string>&);
// A22^-1 X (n2 x q coluna-major) pelo Schur, sem formar A22^-1; cache opcional da simbolica
// do bloco nao genotipado (o padrao se repete nas iteracoes do pseudo-EM)
std::vector<double> a22_inversa_vezes(const Csc& ainv, const std::vector<std::size_t>& geno,
                                      const std::vector<double>& X, std::size_t q,
                                      CacheSimbolica& (*cache)(const Csc&) = nullptr);
// os tres resumos do ML de um metafundador: 1'A22^-1 1, soma_j z_j'A22^-1 z_j e
// soma_j (z_j'A22^-1 1)^2, z_j = m_j - 1 (sem dividir por m/2)
struct ResumosMl { double s = 0.0, t = 0.0, r = 0.0; };
ResumosMl gamma_ml_resumos(const Csc& ainv, const std::vector<std::size_t>& geno,
                           const Genotipos&);
// fracoes de genes dos metafundadores, n x q por linhas, pedigree em ordem topologica
std::vector<double> fracoes_genes(const std::vector<std::int64_t>& pai,
                                  const std::vector<std::int64_t>& mae,
                                  const std::vector<std::int64_t>& col_mf, std::size_t q);
// efeitos de SNP retro-calculados com a G* do ajuste
struct EfeitosSnp {
  std::vector<double> efeito, freq;    // por marcador em uso
  std::vector<std::size_t> usados;     // os marcadores em uso (base 0)
  std::vector<double> diag_gstar;      // a diagonal da G* refeita, para conferir com o ajuste
  double denom = 0.0, a = 0.0, b = 1.0;
  int rota = 0;                        // 1 exata, 2 APY, 3 Vecchia
};
// (hp: a diagonal de G* guardada no ajuste, conferida antes do trabalho pesado; vazia = sem)
EfeitosSnp efeitos_snp(const Pedigree&, const std::vector<std::string>& geno_ids,
                       const Genotipos&, const std::vector<double>& u, double mistura,
                       const std::vector<std::string>& nucleo_apy, std::size_t vecchia_k,
                       const std::vector<double>& hp);

// ---- sssnp.cpp: o passo unico sem G, marcadores como equacoes, PCG, A22^-1 livre de
// matriz (Liu et al. 2014; Masuda et al. 2017; Vandenplas et al. 2018, 2019). theta e
// dado: e um resolvedor.
struct SnpBlup {
  bool ok = false, convergiu = false;
  std::size_t iters = 0;
  double residuo = 0.0;              // ||r|| / ||rhs|| ao parar
  std::vector<double> solucao;       // as colunas do desenho de base (fixos + grupos)
  std::vector<double> efeitos;       // m x fatias efeitos de marcador (coluna por fatia), NaN
                                     // nos monomorficos
  std::vector<std::string> fatias;   // um componente de um grupo com parentesco por coluna
  std::vector<char> usa_marcador;
  std::size_t n_imputados = 0, n_monomorficos = 0;
  std::size_t n_fixo = 0;
  std::vector<std::size_t> offset_grupo;
};
// ---- pegs.cpp: SNP-BLUP multicaracter por Gauss-Seidel aleatorizado (Xavier e Habier 2022)
struct ResultadoPegs {
  std::vector<double> mu, h2, ve;
  Densa b, gebv, vb;
  double deflate = 1.0;
  std::size_t iters = 0, n_imputados = 0;
  bool convergiu = false;
};
ResultadoPegs pegs(const Densa& y, const Genotipos& x, std::size_t maxit, double tol,
                   double deflate_min, bool atualiza_vc, const Densa* vb0,
                   const std::vector<double>* ve0, int estrutura, std::size_t nfat);
void estrutura_pegs(Densa& vb, int tipo, std::size_t nfat);

// meio = true: pedigree com metafundadores, Z centrada em 0.5 e escala m/2 (a G05)
SnpBlup snp_blup(const Desenho&, const Genotipos&, const std::vector<std::string>&,
                 const std::vector<double>&, double, std::size_t, double, bool = false,
                 bool meio = false);

// ---- gibbs.cpp: o lado bayesiano, bloco unico de localizacao + condicionais conjugadas
struct GibbsSaida {
  bool ok = false;
  // bloco denso final do fator e total de colunas: o k de k^3 por amostra
  std::size_t bloco_denso = 0, colunas_fator = 0;
  std::size_t n_amostras = 0, ntheta = 0;
  std::vector<double> amostras;    // n_amostras x ntheta, por linha
  std::vector<double> media_loc, var_loc;
  // media e variancia a posteriori de b na unidade do dado: cada amostra de b volta da escala
  // padronizada das covariaveis antes de entrar (o intercepto leva o centro, entao o dp dele
  // pede a covariancia com as inclinacoes, que so a amostra tem)
  std::vector<double> media_b, var_b;
  std::vector<double> partida;     // theta com que a cadeia comeca (a primeira iteracao)
  std::string mensagem;
  // residuo realizado e_s = y - W loc_s (no probit, a liability menos W loc_s): a media a
  // posteriori por linha da tabela (NaN fora das usadas) e, para as discrepancias pivotais
  // (Johnson 2007; Yuan e Johnson 2012), n_pdm sorteios espacados de e_s / sqrt(s2e_s),
  // n_pdm x n_usadas por linha, cada um N(0, I) sob o modelo
  std::vector<double> media_residuo, residuos_pdm;
  std::size_t n_pdm = 0;
  // W E[loc] e X E[b] por linha da tabela (NaN fora das usadas): pela linearidade, a media a
  // posteriori do ajustado condicional e do fixo (no probit, na escala da liability)
  std::vector<double> ajustado, ajustado_fixo;
  // a escala A da meia-t de cada variancia que a usa (por posicao de theta; NaN nas outras)
  std::vector<double> escala_meia_t;
  // o passo ASIS (gibbs.cpp), por grupo, nas iteracoes depois do aquecimento: quantas vezes
  // rodou, quantas propostas foram aceitas e quantas cairam fora (algum L*_ii <= 0). Vazios
  // sem o passo
  std::vector<double> asis_tentativas, asis_aceites, asis_fora;
  // os grupos de covariancia que recuaram para priori_escalar por terem niveis de menos para a
  // priori deles (ExtrasGibbs::recuo_poucos_niveis), pelo nome do grupo, com os niveis (nl) e a
  // dimensao (d) de cada um, para a mensagem do R
  std::vector<std::string> grupos_recuo;
  std::vector<int> grupos_recuo_nl, grupos_recuo_d;
};
// priori das variancias: 0 Jeffreys (1/s2), 1 plana na variancia (o padrao dos grupos de
// covariancia no gaussiano), 2 uniforme no
// desvio-padrao (so grupo escalar), 3 qui-quadrado inversa escalada propria com df e escala
// (OPTION prior do gibbsf90), 4 meia-t hierarquica (Huang e Wand 2013) com gl0 = nu e a
// escala automatica vezes um multiplicador: `escala` e o multiplicador comum, ou, com
// nomes_mult, o de cada variancia pelo nome ("var(animal)"; as nao nomeadas ficam com 1)
struct PrioriGibbs {
  int tipo = 0;
  double gl0 = 0.0, escala = 0.0;   // gl0: graus de liberdade (df colide com macro do Rmath)
  std::vector<std::string> nomes_mult;
  std::vector<double> mult;
};
// o resto das opcoes da cadeia (gibbs.cpp): dp_fixos > 0 poe a priori N(0, dp_fixos^2) em cada
// coeficiente fixo como o motor o codifica (0 = a plana de sempre); px liga o passo PX-DA de
// Haar do probit (Wang e Roy 2018, Algoritmo 2); partida, quando nao vazia, e o theta de
// onde a cadeia sai, na ordem de nomes_theta() (nomes_partida, quando nao vazio, confere a
// ordem pelos nomes); asis > 0 liga o passo de entrelacamento (Yu e Meng 2011) de cada grupo:
// 1 e o passo; 2 (sem o jacobiano) e 3 (aceita sempre, sem Metropolis-Hastings) sao os
// bracos ERRADOS do portao de validation/gibbs_asis.R, que tem de reprova-los; com
// asis_so_grupos o passo roda so nos grupos de dimensao > 1 (o padrao do gibbs() no probit,
// onde o passo foi medido so nos grupos de covariancia). Com escalar_a_parte, cada termo
// escalar (grupo de dimensao 1) usa priori_escalar no lugar da `priori` dos grupos: e o
// padrao gaussiano do gibbs(), a meia-t nos escalares e a plana nos grupos de covariancia (o
// V5 reprovou a meia-t nos grupos, o V3 a de Jeffreys, o V6 aprovou a meia-t nos
// escalares). priori_escalar nao leva multiplicador por nome. Com recuo_poucos_niveis (so com
// escalar_a_parte), um grupo de dimensao d > 1 cuja priori deixa a InvWishart condicional sem
// media finita (graus de liberdade ate d + 1; na plana, nl - d - 1, ou seja nl <= 2d + 2; ate
// d - 1 ela e impropria e, sem o recuo, para a cadeia com erro) usa priori_escalar no lugar
// dela, e o nome vai em GibbsSaida::grupos_recuo. Com nl = 2d + 3 ou 2d + 4 a plana fica: a
// media existe, a variancia da condicional nao (graus de liberdade ate d + 3), e o dp a
// posteriori do grupo tambem nao
struct ExtrasGibbs {
  double dp_fixos = 0.0;
  bool px = false;
  int asis = 0;
  bool asis_so_grupos = false;
  std::vector<double> partida;
  std::vector<std::string> nomes_partida;
  PrioriGibbs priori_escalar;
  bool escalar_a_parte = false;
  bool recuo_poucos_niveis = false;
};
// priori_res: a do residuo, quando vem a parte da dos grupos (nulo = a mesma `priori`)
GibbsSaida gibbs(const Desenho&, std::size_t n_iter, std::size_t burnin, std::size_t thin,
                 bool loc_fixa, const std::vector<double>* theta_fixo, bool verboso = false,
                 const PrioriGibbs& priori = PrioriGibbs(), bool probit = false,
                 std::size_t n_pdm = 0, const PrioriGibbs* priori_res = nullptr,
                 const ExtrasGibbs* extras = nullptr);
// pecas do gibbs() que a curva de crescimento (curva.cpp) reaproveita (gibbs.cpp): x com
// L' x = z (L triangular inferior em CSC), o sorteio InvWishart(nu, S) por Bartlett (com
// w_out, tambem a Wishart sorteada, C^-1) e os graus de liberdade da condicional de um bloco
// dim x dim com nl niveis pela priori escolhida
std::vector<double> resolve_lt(const Csc& L, const std::vector<double>& z);
Densa rinvwishart(double nu, const Densa& S, Densa* w_out = nullptr);
double gl_priori(const PrioriGibbs& pr, std::size_t nl, std::size_t dim);
// o portao V1 da meia-t: so as condicionais da priori, sem dado (gibbs.cpp)
std::vector<double> gibbs_priori_meia_t(std::size_t dim, double nu, const std::vector<double>& A,
                                        std::size_t n, std::size_t thin, std::size_t burnin,
                                        int variante);
Ajuste ajusta(const Desenho&, const std::vector<double>*, std::size_t, std::size_t, double,
              bool = false);


// ---- multitrait.cpp / multitrait2.cpp
struct DesenhoMT {
  Modelo modelo;
  std::vector<std::string> alvos;
  std::size_t t = 0;
  // eq_fixa[j*t + tau] -> indice da equacao fixa, ou -1 quando o par (coluna, caracter)
  // nao tem nenhum registro e portanto nao e estimavel. O bloco fixo era kron CHEIO,
  // todas as colunas vezes todos os caracteres, e um par sem registro nascia como equacao
  // vazia que matava a Cholesky. Quando nada cai, eq_fixa[k] == k e n_fixa == x.ncol * t:
  // o layout antigo e caso particular deste, o que deixa intacto todo ajuste que ja andava.
  std::vector<int> eq_fixa;
  std::size_t n_fixa = 0;
  Csr x;   // por registro, como no Desenho
  std::vector<std::string> nomes_x;
  std::vector<std::string> saiu_x;
  // a padronizacao das covariaveis: a volta de b por caracteristica e o termo do -2logL (ver
  // Desenho)
  std::vector<VoltaB> volta_x;
  double logdet_escala_x = 0.0;
  // o motivo de cada queda, paralelo a saiu_x: "no record" / "aliased" para a coluna global,
  // "no record of <caracter>" / "aliased within <caracter>" para o par (coluna, caracter)
  std::vector<std::string> saiu_motivo;
  std::vector<DesenhoTermo> aleatorios;
  std::vector<Csc> kinv;
  std::vector<double> kinv_logdet;
  std::vector<double> kinv_escala;   // escala de partida por grupo (0 = do log-determinante)
  Densa y;
  std::vector<char> usa;
  // obs[r*t + tau]: a caracteristica tau foi observada no registro r. Um registro pode
  // participar com um subconjunto, e a R0 dele e a submatriz do padrao.
  std::vector<char> obs;
  std::size_t nlin = 0;
  std::size_t n_usadas() const { std::size_t s = 0; for (char u : usa) s += u; return s; }
  std::size_t largura(std::size_t g) const {
    std::size_t s = 0;
    for (std::size_t tm : modelo.grupos[g].termos)
      for (const DesenhoTermo& d : aleatorios)
        if (d.termo == tm) s += d.z.ncol * d.nt(t);
    return s;
  }
  std::size_t total_colunas() const {
    std::size_t s = x.ncol * t;
    for (std::size_t g = 0; g < modelo.grupos.size(); g++) s += largura(g);
    return s;
  }
};
struct AvaliacaoMT {
  bool ok = false;
  // POR QUE falhou, quando ok e false, e as duas causas pedem acoes opostas de quem
  // chama. Antes as duas saiam com a mesma mensagem, a de theta, e quem usava mexia no
  // start para sempre sem chegar perto do problema.
  //   1 = alguma covariancia nao e positiva-definida. O THETA e a causa.
  //   2 = a matriz de coeficientes e singular. O DESENHO e a causa, e start nao resolve.
  int motivo = 0;
  double neg2logl = 0.0;
  std::vector<double> solucao;
  // Proposta EM dos GRUPOS, C_g <- (Q + T)/nl, com os componentes que este laco nao
  // atualiza deixados no valor corrente. E uma proposta, nao um passo: quem chama so a
  // aceita se ela baixar a verossimilhanca, entao um EM parcial e seguro por construcao.
  // O passo EM e multiplicativo, fica no cone e nao encolhe na fronteira, que e onde o
  // passo AI amortecido trava. Sem o fator s2e do univariado: a MME daqui ja e absoluta.
  std::vector<double> em_theta;
  std::vector<double> pev;   // diagonal de C^-1 na numeracao original (unidades absolutas)
  std::vector<double> score;
  Densa ai;
  std::size_t fora_do_padrao = 0;
};
// As MME do multicaracter montadas num theta (multitrait.cpp), em unidades ABSOLUTAS. Fica no
// cabecalho porque o nivel 2 da curva de crescimento (curva.cpp) sorteia as localizacoes
// por esta mesma montagem.
struct MontadoMT {
  Csc c;                       // triangulo superior, unidades absolutas
  std::vector<double> rhs;
  double logdet_g = 0.0;
  double logdet_r = 0.0;       // n log|R0|
  Densa r0, r0inv;
  // Por padrao de observacao (mascara de bits sobre as caracteristicas): a inversa da
  // submatriz de R0 EMBUTIDA numa t x t com zeros nas ausentes, e o log|R0_p|. A inversa
  // embutida e o que mantem todas as formulas com lacos sobre t intactas: contribuicao de
  // caracteristica nao observada anula sozinha.
  std::vector<Densa> rinv_mask;
  std::vector<double> ld_mask;
  std::vector<std::uint32_t> mask;   // por linha dos dados
  std::size_t n_fixo = 0;      // p * t
  std::size_t total = 0;
  std::vector<std::size_t> offset_grupo;
  std::vector<BlocoT> blocos_t;
  bool ok = false;
};
MontadoMT monta_mme_mt(const DesenhoMT& d, const std::vector<double>& theta,
                       const std::vector<Densa>* fatores = nullptr);
AvaliacaoMT avalia_mt(const DesenhoMT&, const std::vector<double>&, CacheSimbolica* = nullptr);
double neg2logl_densa_V_mt(const DesenhoMT&, const std::vector<double>&);
DesenhoMT monta_desenho_mt(Modelo, const std::vector<std::string>&, const Tabela&,
                           const Pedigree*,
                           const std::vector<KernelDecl>* = nullptr);
std::vector<std::string> nomes_theta_do_mt(const DesenhoMT&);
struct AjusteMT {
  bool convergiu = false;
  // o cache simbolico (ordem e fator simbolico) do caminho que deu este ajuste: residuos_mt()
  // reaproveita quando o padrao e o mesmo, em vez de refazer a ordenacao no fim
  std::shared_ptr<CacheSimbolica> cache;
  // bloco denso final do fator e total de colunas: o k de k^3 por iteracao
  std::size_t bloco_denso = 0, colunas_fator = 0;
  std::size_t iters = 0;
  double reldelta = 0.0;
  // Newton decrement at the final point, restricted to the OFF-BOUNDARY components and
  // computed in the log-Cholesky coordinates the step walks in (passo_z / decremento_z),
  // over the same floors and the same active set the step uses. NaN when the fit died
  // before a first evaluation.
  double decremento = std::numeric_limits<double>::quiet_NaN();
  double neg2logl = 0.0;
  std::vector<double> theta, se, solucao;
  // 2 AI^-1, a covariancia amostral das componentes (ntheta x ntheta, por linha), NaN na
  // linha e na coluna de um componente sem erro-padrao (lista de correcoes, item 16: sem ela
  // nao ha metodo delta para rg, h2 ou h2(x) do model_mt() e do model_ar1())
  std::vector<double> vcov;
  std::vector<double> pev;   // diagonal de C^-1 no otimo, ja em unidades absolutas
  // a partida de onde a caminhada saiu, em unidades originais: o start= do usuario, a dos
  // univariados ou a automatica (para conferir, e para repetir o ajuste)
  std::vector<double> partida;
  std::string mensagem;
  std::size_t fora_do_padrao = 0;
  // o conjunto que o ultimo certificado excluiu (PassoZ::na_parede, indexado como theta) e se
  // a caminhada parou por um veredito sem certificar (estagnacao ou limite de reentradas),
  // e nao pelo teto de iteracoes: o que a exploracao da fronteira le
  std::vector<char> na_parede;
  bool parou = false;
  // AR(1): o ajuste devolvido e o da segunda caminhada, que partiu do rho trocado de sinal
  bool espelhado = false;
};
// A EXPLORACAO DA FRONTEIRA (o racional e as medidas estao em aireml.cpp): depois de uma
// caminhada que certificou, ou parou antes do teto, com algum bloco de covariancia na
// parede, uma segunda caminhada parte de dentro do cone, com cada bloco na parede
// descorrelacionado e as variancias dele em pelo menos 10% da partida de var(y) (`auto_`,
// em unidades originais), e volta o melhor dos dois. `caminha(start, maxiter)` e a
// caminhada do ajustador, com o orcamento que sobrou.
AjusteMT explora_fronteira(AjusteMT R, const MapaZ& mp, const std::vector<double>& auto_,
                           std::size_t maxiter,
                           const std::function<AjusteMT(const std::vector<double>&,
                                                        std::size_t)>& caminha,
                           bool verboso);
// O desenho entra por referencia NAO constante: o ajuste padroniza y por caracteristica no
// lugar durante a caminhada e devolve os valores originais, bit a bit, antes de sair (o
// racional esta em multitrait2.cpp).
AjusteMT ajusta_mt(DesenhoMT&, const std::vector<double>*, std::size_t, double,
                   bool = false);
// A partida padrao do multicaracter quando o R manda os otimos univariados de cada
// caracteristica (uni[tau] no layout do model() de uma caracteristica): grupos e residuo na
// diagonal por caracteristica, covariancias entre caracteristicas zero. Vazio quando algum
// vetor nao tem o tamanho do layout univariado ou nao e finito.
std::vector<double> partida_univariados_mt(const DesenhoMT&,
                                           const std::vector<std::vector<double> >& uni);
// o mesmo mapa por grupo, compartilhado pelo multicaracter e pelo AR(1); false quando algum
// vetor nao casa com o layout
bool grupos_de_univariados(const Modelo& m, std::size_t t,
                           const std::vector<std::vector<double> >& uni,
                           std::vector<double>& theta);
// A padronizacao interna por caracteristica dos ajustadores multicaracter e AR(1)
// (multitrait2.cpp): o desvio padrao de cada caracteristica nas celulas observadas, o fator
// de cada componente (s_a s_b; 1 para o rho) e a guarda que divide y no lugar e o devolve
// bit a bit (na saida normal por restaura(), em qualquer outra pelo destrutor).
// `exatos`, quando dado, recebe o desvio padrao exato de cada caracteristica (1 sem
// variacao): o do corte de singularidade, que tem de ser invariante a qualquer fator, e nao
// so a potencias de 2 (o racional esta em aviso_informacao, aireml.cpp)
std::vector<double> desvios_tracos(const Densa& y, const std::vector<char>& usa,
                                   const std::vector<char>& obs, std::size_t t,
                                   std::vector<double>* exatos = nullptr);
std::vector<double> fatores_escala(const Modelo& m, std::size_t t, std::size_t off_residuo,
                                   const std::vector<double>& esc);
// a caracteristica de cada coordenada de um grupo expandido, na ordem do layout (termo por
// fora, caracteristica, coeficiente por dentro); vazio quando a dimensao nao casa
std::vector<std::size_t> traco_das_coordenadas(const Modelo& m, const Grupo& g, std::size_t t);
struct PadronizaY {
  Densa& y;
  Densa original;
  bool ativo = true;
  PadronizaY(Densa& yy, const std::vector<char>& usa, const std::vector<char>& obs,
             const std::vector<double>& esc, std::size_t t) : y(yy), original(yy) {
    for (std::size_t i = 0; i < y.nlin; i++)
      if (usa[i])
        for (std::size_t tau = 0; tau < t; tau++)
          if (obs[i * t + tau]) y.at(i, tau) /= esc[tau];
  }
  void restaura() {
    if (!ativo) return;
    y = original;
    ativo = false;
  }
  ~PadronizaY() { restaura(); }
};
// passo unico na multi: o mesmo nucleo do uni, sobre os mesmos campos
RelatorioG aplica_genomica(DesenhoMT&, const Pedigree&, const std::vector<std::string>&,
                           const Genotipos&, double, const std::vector<std::string>&,
                           std::size_t = 0);


// ---- ar1.cpp / ar1b.cpp
// Entrada de linha do esqueleto AR: coluna global, caracteristica e valor.
struct EntAR {
  std::uint32_t col;
  std::uint32_t trait;
  double val;
};
struct DesenhoAR {
  Modelo modelo;
  std::vector<std::string> alvos;   // t = alvos.size(); t = 1 e o caminho classico
  std::size_t t = 1;
  Csr x;   // por registro, como no Desenho
  std::vector<std::string> nomes_x;
  std::vector<std::string> saiu_x;
  // a padronizacao das covariaveis: a volta de b por caracteristica (todas iguais aqui, toda
  // coluna tem equacao em toda caracteristica) e o termo do -2logL (ver Desenho)
  std::vector<VoltaB> volta_x;
  double logdet_escala_x = 0.0;
  // o motivo de cada queda, paralelo a saiu_x: "no record" ou "aliased" (aqui so ha queda
  // global: um par (coluna, caracteristica) sem registro e recusado na montagem)
  std::vector<std::string> saiu_motivo;
  std::vector<DesenhoTermo> aleatorios;
  std::vector<Csc> kinv;
  std::vector<double> kinv_logdet;
  std::vector<double> kinv_escala;   // escala de partida por grupo (0 = do log-determinante)
  Densa y;                          // nlin x t
  std::vector<char> usa;
  std::size_t nlin = 0;
  std::vector<std::vector<std::size_t>> sujeitos;
  // o tempo de cada linha COMO O MOTOR O USA: na grade (tempo_inteiro), o numero de passos
  // desde o primeiro registro usado do sujeito, inteiro exato; fora dela, o tempo do dado.
  // `tempo_dado` e o tempo como veio, para a saida. Quem decide e eixo_tempo_ar1 (ar1b.cpp).
  std::vector<double> tempo;
  std::vector<double> tempo_dado;
  // Os intervalos dentro de um sujeito caem todos numa GRADE? A pergunta decide se rho
  // negativo e um modelo. Gamma(dt) = rho^dt so e uma funcao de correlacao valida com rho < 0
  // quando dt e um numero inteiro de passos: fora da grade a convencao de sinal quebra a
  // multiplicatividade, s(2.5)|rho|^2.5 ao quadrado da +|rho|^5, e s(5)|rho|^5 e -|rho|^5,
  // e com ela a propriedade de Markov de que a Gamma^-1 tridiagonal depende. A rota densa e
  // a esparsa passam a discordar. Em tempo continuo o modelo e o CAR(1), com rho em (0, 1).
  bool tempo_inteiro = true;
  // o PASSO da grade na unidade do dado (o rho e por passo); NaN em tempo continuo
  double passo_tempo = 1.0;
  // na grade, ha intervalos com um numero PAR e com um numero IMPAR de passos, cada paridade
  // em ao menos 10% deles? Ai rho e -rho dao a mesma correlacao nos pares e a oposta nos
  // impares, e a verossimilhanca pode ter um modo de cada sinal (o espelho de ajusta_ar1)
  bool paridade_mista = false;
  // Quantos registros sairam por terem alguma caracteristica AUSENTE. O multicaracter
  // comum mantem esse registro e o ajusta contra a submatriz de R0 do padrao dele; aqui a
  // separabilidade Gamma (x) R0 nao sobrevive a um padrao parcial e o registro sai inteiro.
  // E exclusao por lista, e exclusao por lista tem de ser DITA: sem isso o unico sinal e um
  // n_used menor que nrow(data), que ninguem confere.
  std::size_t n_incompletos = 0;    // registros sem NENHUMA caracteristica observada
  std::size_t n_sem_tempo = 0;     // registros sem tempo finito, fora da serie
  // CELULAS AUSENTES (restricao #6). obs[r*t + tau] diz se a caracteristica tau do registro
  // r foi observada. Cada celula ausente de um registro usado ganha uma coluna de efeito
  // FIXO propria no fim do vetor de solucoes (off_mv + k), com y = 0 na celula: e o termo
  // mv do ASReml, que mantem o residuo Gamma (x) R0 separavel e cheio. O -2logL, o score e
  // a AI dessa MME aumentada sao EXATAMENTE os da verossimilhanca marginal das celulas
  // observadas (complemento de Schur no bloco mv), sem mudar formula nenhuma do avaliador.
  std::vector<char> obs;
  std::size_t n_mv = 0, off_mv = 0;
  // Esqueleto pre-computado no desenho (nada disto depende de theta): as linhas de W ja
  // em colunas globais COM a caracteristica.
  std::vector<std::vector<EntAR>> lw;
  // offset_s2e = inicio do vech(R0) (t = 1: o proprio s2e); offset_rho = logo depois
  std::size_t offset_s2e = 0, offset_rho = 0;
  std::size_t n_usadas() const { std::size_t s = 0; for (char u : usa) s += u; return s; }
  std::size_t largura(std::size_t g) const {
    std::size_t s = 0;
    for (std::size_t tm : modelo.grupos[g].termos)
      for (const DesenhoTermo& d : aleatorios)
        if (d.termo == tm) s += d.z.ncol * t;
    return s;
  }
  std::size_t total_colunas() const {
    std::size_t s = x.ncol * t;
    for (std::size_t g = 0; g < modelo.grupos.size(); g++) s += largura(g);
    return s + n_mv;
  }
};
struct AvaliacaoAR {
  bool ok = false;
  double neg2logl = 0.0;
  std::vector<double> solucao, score;
  // Proposta EM dos GRUPOS, C_g <- (Q + T)/nl, com os componentes que este laco nao
  // atualiza deixados no valor corrente. E uma proposta, nao um passo: quem chama so a
  // aceita se ela baixar a verossimilhanca, entao um EM parcial e seguro por construcao.
  // O passo EM e multiplicativo, fica no cone e nao encolhe na fronteira, que e onde o
  // passo AI amortecido trava. Sem o fator s2e do univariado: a MME daqui ja e absoluta.
  std::vector<double> em_theta;
  std::vector<double> pev;   // diagonal de C^-1 na numeracao original (unidades absolutas)
  Densa ai;
  std::size_t fora_do_padrao = 0;
};
AvaliacaoAR avalia_ar1(const DesenhoAR&, const std::vector<double>&, CacheSimbolica* = nullptr);
double neg2logl_densa_V_ar1(const DesenhoAR&, const std::vector<double>&);
DesenhoAR monta_desenho_ar1(Modelo, const std::vector<std::string>&, const Tabela&,
                            const Pedigree*, const std::string&, const std::string&,
                            const std::vector<KernelDecl>* = nullptr);
std::vector<std::string> nomes_theta_ar1(const DesenhoAR&);
// por referencia nao constante pelo mesmo motivo de ajusta_mt (y padronizado no lugar)
AjusteMT ajusta_ar1(DesenhoAR&, const std::vector<double>*, std::size_t, double,
                    bool = false);
// a partida a partir dos otimos univariados AR(1) de cada caracteristica (layout do
// model_ar1() de uma caracteristica: grupos, var(residual), rho): R0 diagonal e o rho medio.
// `passos`: o time_step de cada ajuste univariado (NaN no continuo); quando algum difere do
// passo deste desenho, os rho sao de unidades diferentes e o rho parte da automatica
std::vector<double> partida_univariados_ar1(const DesenhoAR&,
                                            const std::vector<std::vector<double> >& uni,
                                            const std::vector<double>& passos);
// passo unico no AR(1): o mesmo nucleo do uni, sobre os mesmos campos
RelatorioG aplica_genomica(DesenhoAR&, const Pedigree&, const std::vector<std::string>&,
                           const Genotipos&, double, const std::vector<std::string>&,
                           std::size_t = 0);

// Residuos no theta final do multicaracter e do AR(1): por linha da tabela e caracteristica
// (indice r * t + tau), o residuo condicional, os ajustados condicional e fixo e
// Var(e^_{r,tau}) = R0[tau,tau] - w' C^-1 w (MME absolutas). NaN em linha nao usada e em
// celula ausente. Uma montagem e um fator a mais no fim do ajuste (o fator da ultima
// avaliacao nao sobrevive a ela).
struct ResiduosMT {
  bool ok = false;
  std::size_t t = 0, nlin = 0, fora = 0;
  std::vector<double> residuo, ajustado, ajustado_fixo, variancia;
};
ResiduosMT residuos_mt(const DesenhoMT& d, const std::vector<double>& theta,
                       const CacheSimbolica* cache = nullptr);
ResiduosMT residuos_ar1(const DesenhoAR& d, const std::vector<double>& theta);
// bootstrap parametrico preparado do AR(1)/CAR(1) univariado (residuos.cpp tem o do model())
std::unique_ptr<PreparoNulos> prepara_nulos_ar1(DesenhoAR&& d, const std::vector<double>& theta,
                                                const std::vector<double>& b);

// ---- curva.cpp: a curva de crescimento nao linear de gibbs_growth() (o racional esta la)
struct CurvaConf {
  double nu = 0.0;          // 1 logistica, 0 Gompertz, -1/3 von Bertalanffy, Richards dado
  bool log_escala = true;   // o modelo em log y (erro multiplicativo)
  bool linear = false;      // a curva de teste phi1 + phi2 t + phi3 t^2 do portao E4
};
// a curva na escala de analise (log f na escala log, f na crua, o polinomio na de teste) e o
// jacobiano dela em phi = (log A, ti, log k); false fora do suporte (1 + nu e <= 0, so com
// nu < 0) ou com transbordamento, e ai f = -inf
bool curva_analise(const CurvaConf& c, double t, const double* phi, double& f, double* jac);
// o peso que a curva da em t (0 fora do suporte): a escala dos derivados
double curva_peso(const CurvaConf& c, double t, const double* phi);

// As pesagens agrupadas por animal (contagem estavel): as do animal i em [ini[i], ini[i+1]),
// y na escala de analise, a classe do residuo, o log do peso (o jacobiano do WAIC na escala
// log) e `ordem`, a posicao de entrada de cada uma
struct PesagensCurva {
  std::size_t n_anim = 0, n_classes = 1;
  std::vector<std::size_t> ini, ordem;
  std::vector<double> t, y, logpeso;
  std::vector<int> classe;
};
PesagensCurva agrupa_pesagens(const std::vector<int>& animal, const std::vector<double>& t,
                              const std::vector<double>& y, const std::vector<double>& logpeso,
                              const std::vector<int>& classe, std::size_t n_anim,
                              std::size_t n_classes);

// O modo condicional de cada animal pelo Gauss-Newton do passo 1, com a media mu (n_anim x
// 3, por linhas; uma linha so vale para todos), a covariancia `sigma` do nivel 2 e o s2e de
// cada classe: as partidas S1 e S2 de gibbs_growth(). m e h (a matriz de Gauss-Newton no
// modo, 3 x 3 por linhas) por animal, os passos dados e se o modo ficou no suporte
struct ModosCurva {
  std::vector<double> m, h;
  std::vector<int> passos;
  std::vector<char> ok;
};
ModosCurva modos_curva(const CurvaConf& c, const PesagensCurva& p, const std::vector<double>& mu,
                       const Densa& sigma, const std::vector<double>& s2e);

// O passo 1 sozinho para um animal, com mu, a covariancia, s2e e nu dados (o portao E2): a
// cadeia guardada (n x 3, por linhas, depois do aquecimento com adaptacao), as aceitacoes de
// cada movimento depois dele, a escala congelada do passeio e o modo com a matriz dele
struct PassoCurvaSaida {
  std::vector<double> cadeia;
  double m[3] = {0.0, 0.0, 0.0};
  double h[9] = {0.0};
  double ac_ind = 0.0, ac_rw = 0.0, escala = 0.0;
  int passos = 0;
  bool ok = false;
};
PassoCurvaSaida passo_curva_isolado(const CurvaConf& c, const PesagensCurva& p,
                                    const double* mu, const Densa& sigma,
                                    const std::vector<double>& s2e, const double* phi0,
                                    std::size_t burnin, std::size_t n, int variante,
                                    int movimentos);

// a priori de um bloco da cadeia: tipo 1 plana, 2 uniforme no desvio (escalar), 3 propria
// InvWishart(gl0, gl0 diag(escala)), uma escala por parametro (por classe no residuo, a mesma)
struct PrioriCurva {
  int tipo = 1;
  double gl0 = 0.0;
  std::vector<double> escala;
};
struct OpcoesCurva {
  CurvaConf curva;
  std::size_t n_iter = 0, burnin = 0, thin = 1;
  std::vector<PrioriCurva> priori_grupos;   // um por grupo do desenho
  PrioriCurva priori_pe, priori_res;
  // o dp da priori N(0, dp^2) dos coeficientes fixos de cada parametro (log A, ti, log k), na
  // unidade dele; 0 = plana naquele parametro
  std::vector<double> dp_fixos = std::vector<double>(3, 0.0);
  std::vector<double> theta0;                // grupos + Sigma_p no layout do motor
  std::vector<double> s2e0;                  // por classe
  std::vector<double> phi0;                  // n_anim x 3, por linhas
  int variante = 1;    // 1 o passo; 2 sem a razao q; 3 sem a priori do nivel 2 na aceitacao;
                       // 4 o ruido do nivel 2 vezes raiz de 2 (os bracos errados dos portoes)
  int movimentos = 0;  // 0 os dois, 1 so o passeio, 2 so o independente
  int segura = 0;      // 1 phi e componentes presos (E3); 2 so phi preso (E8)
  std::size_t n_pdm = 0;
  std::vector<double> idades, grade;         // ages = e a grade das curvas
  std::vector<int> classe_idade;             // a classe do residuo de cada idade
  int grupo_genetico = -1;                   // o grupo dos derivados (-1: sem)
  bool verboso = false;
};
struct SaidaCurva {
  std::size_t n_amostras = 0, ncomp = 0;
  std::vector<double> amostras;              // n_amostras x ncomp, por linhas
  std::vector<double> media_phi, var_phi;    // n_anim x 3
  std::vector<double> media_loc, var_loc;
  std::vector<double> media_b, var_b;        // na unidade do dado, na ordem de eq_fixa
  std::vector<double> ac_ind, ac_rw, escala_rw, tentativas;   // por animal
  std::vector<double> lse, media_ld, m2_ld;  // WAIC, por pesagem (ordem agrupada)
  std::vector<double> media_ajuste;          // a curva media por pesagem (escala de analise)
  std::vector<double> pdm;                   // n_pdm x pesagens
  std::size_t n_pdm = 0;
  std::vector<double> media_der, m2_der;     // niveis x (idades + 2): g(t), W_inf, R_max
  std::vector<double> media_gc, m2_gc;       // niveis x grade: a curva genetica
  std::vector<double> media_ac;              // n_anim x grade: a curva de cada animal
  std::vector<double> pop;                   // n_amostras x grade: a curva da populacao
  std::vector<double> h2_idade;              // n_amostras x idades: o h2 de primeira ordem
  std::vector<double> sol_primeira;          // segura = 1: a media condicional da 1a varredura
  std::vector<double> b_primeira;            // a parte fixa dela, na unidade do dado
  std::vector<double> loglik;                // n_amostras: o log L do dado na escala de analise
  std::size_t modos_fora = 0;                // modos que cairam fora do suporte
  std::string mensagem;
};
SaidaCurva gibbs_curva(DesenhoMT& d, const PesagensCurva& p, const OpcoesCurva& o);

}  // namespace br

#endif
