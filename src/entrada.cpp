// A travessia para o R. Converte objetos R em estruturas do pacote e de volta.
//
// Nenhuma conta acontece aqui. Se um dia aparecer aritmetica neste arquivo, ela estara fora
// dos testes que validam a numerica, e sera por ai que os numeros vao divergir.
//
// Uma excecao de C++ atravessando a fronteira do R deixa o interpretador em estado
// indefinido, entao TODO ponto de entrada embrulha o corpo e converte a excecao em Rf_error.

#include "mme.h"
#include <cstdio>
#include <cstring>
#include <fstream>
#include <unordered_map>

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


// Le triplos do R (base 1, triangulo qualquer) e devolve a Csc do triangulo INFERIOR.
static br::Csc csc_do_R(SEXP i, SEXP j, SEXP x, SEXP n) {
  const int nn = Rf_asInteger(n);
  if (nn <= 0) Rf_error("invalid dimension");
  const R_xlen_t nz = XLENGTH(x);
  if (XLENGTH(i) != nz || XLENGTH(j) != nz) Rf_error("i, j and x with different lengths");
  std::vector<std::uint32_t> li, cj;
  std::vector<double> v;
  li.reserve(static_cast<std::size_t>(nz));
  cj.reserve(static_cast<std::size_t>(nz));
  v.reserve(static_cast<std::size_t>(nz));
  for (R_xlen_t k = 0; k < nz; k++) {
    int a = INTEGER(i)[k] - 1, b = INTEGER(j)[k] - 1;
    if (a < 0 || b < 0 || a >= nn || b >= nn) Rf_error("triplet outside the matrix");
    if (a < b) std::swap(a, b);
    li.push_back(static_cast<std::uint32_t>(a));
    cj.push_back(static_cast<std::uint32_t>(b));
    v.push_back(REAL(x)[k]);
  }
  return br::de_triplos(static_cast<std::size_t>(nn), static_cast<std::size_t>(nn), li, cj, v);
}

static SEXP csc_para_R(const br::Csc& m) {
  const R_xlen_t nz = static_cast<R_xlen_t>(m.nnz());
  SEXP li = PROTECT(Rf_allocVector(INTSXP, nz));
  SEXP cj = PROTECT(Rf_allocVector(INTSXP, nz));
  SEXP vv = PROTECT(Rf_allocVector(REALSXP, nz));
  R_xlen_t k = 0;
  for (std::size_t c = 0; c < m.ncol; c++)
    for (std::size_t t = m.colptr[c]; t < m.colptr[c + 1]; t++, k++) {
      INTEGER(li)[k] = static_cast<int>(m.linha[t]) + 1;
      INTEGER(cj)[k] = static_cast<int>(c) + 1;
      REAL(vv)[k] = m.valor[t];
    }
  const char* campos[] = {"i", "j", "x", "n"};
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, 4));
  for (int q = 0; q < 4; q++) SET_STRING_ELT(nms, q, Rf_mkChar(campos[q]));
  SET_VECTOR_ELT(out, 0, li);
  SET_VECTOR_ELT(out, 1, cj);
  SET_VECTOR_ELT(out, 2, vv);
  SET_VECTOR_ELT(out, 3, Rf_ScalarInteger(static_cast<int>(m.ncol)));
  Rf_setAttrib(out, R_NamesSymbol, nms);
  UNPROTECT(5);
  return out;
}

// ASSINATURA de uma precisao em triplos (FNV-1a sobre palavras de 64 bits: a dimensao, (i, j)
// e os bits de x na ordem dos triplos, e por fim os bits do log|K|). Nao e criptografica: so
// diz se os triplos e o log|K| que voltam num kernel(Kinv =) sao os que dominance_inverse()
// devolveu, para que o log|K| pela forma fechada possa ser usado sem refatorar Q. Qualquer
// mudanca, nos triplos (reordenacao inclusive) ou no log|K|, cai na Cholesky.
struct Assinatura {
  std::uint64_t h = 1469598103934665603ULL;
  void palavra(std::uint64_t w) {
    h ^= w;
    h *= 1099511628211ULL;
  }
  void triplo(int i, int j, double x) {
    std::uint64_t bx;
    std::memcpy(&bx, &x, sizeof bx);
    palavra((static_cast<std::uint64_t>(static_cast<std::uint32_t>(i)) << 32) |
            static_cast<std::uint32_t>(j));
    palavra(bx);
  }
  void real(double x) {
    std::uint64_t bx;
    std::memcpy(&bx, &x, sizeof bx);
    palavra(bx);
  }
  std::string hex() const {
    char b[24];
    std::snprintf(b, sizeof b, "%016llx", static_cast<unsigned long long>(h));
    return b;
  }
};

// Converte a excecao do dominio em condicao do R.
//
// Variadica de proposito: chaves NAO protegem virgulas no preprocessador, so parenteses. Com
// um parametro so, qualquer virgula dentro do corpo, e ha uma em cada Rf_allocVector, vira
// separador de argumento e a macro deixa de existir. Foi exatamente esse o erro na primeira
// versao deste arquivo.
//
// Rf_error faz longjmp, entao nada que precise de destrutor pode estar vivo quando ele
// dispara: a mensagem e copiada para um buffer ANTES de sair do bloco try. 2048 e nao 512:
// com ids de 36 caracteres (UUID) a recusa do Delta da dominancia saia cortada no meio do
// remedio, que e a parte que importa.
#define GUARDA(...)                                                              char msg__[2048];                                                               msg__[0] = 0;                                                                  try { __VA_ARGS__ }                                                            catch (const std::exception& e) {                                                std::snprintf(msg__, sizeof msg__, "%s", e.what());                          }                                                                              catch (...) { std::snprintf(msg__, sizeof msg__, "unknown error"); }       Rf_error("%s", msg__);

// Os ajudantes de leitura abaixo (textos, mgs_de, genotipos_do_R) so rodam dentro do try de
// um GUARDA: a recusa e um br::Erro, e nao um Rf_error, para que o longjmp nao pule os
// destrutores dos objetos C++ que a entrada ja montou (o pedigree, a A^-1).
static std::vector<std::string> textos(SEXP v, const char* quem) {
  if (TYPEOF(v) != STRSXP) throw br::Erro(std::string(quem) + ": expected a text vector");
  const R_xlen_t n = XLENGTH(v);
  std::vector<std::string> out;
  out.reserve(static_cast<std::size_t>(n));
  for (R_xlen_t i = 0; i < n; i++) {
    SEXP e = STRING_ELT(v, i);
    if (e == NA_STRING)
      throw br::Erro(std::string(quem) +
                     ": there is NA. A missing identifier would become an animal named 'NA'");
    out.emplace_back(CHAR(e));
  }
  return out;
}

// Uma lista de vetores numericos (os otimos univariados que model_mt() e model_ar1() ajustam
// para a partida padrao). Vazia quando nao e uma lista de doubles: a partida automatica
// fica, sem erro, porque a partida e uma conveniencia e nao o modelo.
static std::vector<std::vector<double> > vetores_do_R(SEXP v) {
  std::vector<std::vector<double> > out;
  if (TYPEOF(v) != VECSXP) return out;
  for (R_xlen_t k = 0; k < XLENGTH(v); k++) {
    SEXP e = VECTOR_ELT(v, k);
    if (TYPEOF(e) != REALSXP) return {};
    out.emplace_back(REAL(e), REAL(e) + XLENGTH(e));
  }
  return out;
}

// a partida do ajuste multicaracter / AR(1), nomeada como theta; vazia quando o ajuste
// parou antes dela
static SEXP sexp_partida(const std::vector<double>& p, SEXP nomes) {
  SEXP v = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) p.size()));
  for (std::size_t k = 0; k < p.size(); k++) REAL(v)[k] = p[k];
  if ((R_xlen_t) p.size() == XLENGTH(nomes)) Rf_setAttrib(v, R_NamesSymbol, Rf_duplicate(nomes));
  UNPROTECT(1);
  return v;
}

// A flag de pedigree pai / avo materno viaja como ATRIBUTO "mgs" do vetor de maes, posto
// por colunas_pedigree() no R. Assim os onze pontos que constroem o pedigree a recebem sem
// mudar assinatura nenhuma, e um vetor sem o atributo e o pedigree pai/mae de sempre.
static std::vector<char> mgs_de(SEXP mae) {
  SEXP a = Rf_getAttrib(mae, Rf_install("mgs"));
  std::vector<char> out;
  if (a == R_NilValue) return out;
  const R_xlen_t n = XLENGTH(a);
  if (n != XLENGTH(mae)) throw br::Erro("the 'mgs' flags do not match the dam column in length");
  if (TYPEOF(a) != LGLSXP) throw br::Erro("the 'mgs' flags must be logical");
  out.resize((std::size_t) n);
  for (R_xlen_t k = 0; k < n; k++) {
    const int v = LOGICAL(a)[k];
    out[(std::size_t) k] = (v == 1) ? 1 : 0;
  }
  return out;
}


// As portas esparsas do R (sparse_chol, sparse_solve, selected_inverse) sao chamadas dentro
// dos lacos de Newton do limiar e da sobrevivencia, sempre com o MESMO padrao e valores
// novos. Cada chamada refazia o grau minimo e a simbolica. Aqui ficam os ultimos quatro
// padroes vistos, pela assinatura do padrao: a mesma permutacao e a mesma simbolica, logo
// o mesmo resultado bit a bit, sem refazer nada.
namespace {
struct PadraoVisto {
  std::uint64_t assinatura = 0;
  std::size_t n = 0, nnz = 0;
  int modo = 0;                 // a ordenacao que fez a permutacao (UAITREM_ORDEM)
  br::CacheSimbolica cs;
};
std::vector<PadraoVisto> g_padroes;

br::CacheSimbolica& simbolica_do_padrao(const br::Csc& a) {
  const std::uint64_t h = br::assinatura_padrao(a);
  for (std::size_t k = 0; k < g_padroes.size(); k++)
    if (g_padroes[k].assinatura == h && g_padroes[k].n == a.ncol &&
        g_padroes[k].nnz == a.nnz() && g_padroes[k].modo == br::modo_ordem()) {
      std::rotate(g_padroes.begin() + k, g_padroes.begin() + k + 1, g_padroes.end());
      return g_padroes.back().cs;
    }
  PadraoVisto p;
  p.assinatura = h; p.n = a.ncol; p.nnz = a.nnz(); p.modo = br::modo_ordem();
  p.cs.perm = br::grau_minimo(a);
  p.cs.sb = br::simbolica(br::permuta_sim(a, p.cs.perm));
  p.cs.pronto = true;
  if (g_padroes.size() >= 4) g_padroes.erase(g_padroes.begin());
  g_padroes.push_back(std::move(p));
  return g_padroes.back().cs;
}
}  // namespace

// A matriz de genotipos do R como ela esta, double, integer ou raw, coluna-major: nenhuma
// copia. As linhas tem de ser os ids declarados.
static br::Genotipos genotipos_do_R(SEXP gm, std::size_t n_ids) {
  SEXP dim = Rf_getAttrib(gm, R_DimSymbol);
  if (dim == R_NilValue || XLENGTH(dim) != 2) throw br::Erro("the genotype matrix must be a matrix");
  br::Genotipos g;
  g.n = (std::size_t) INTEGER(dim)[0];
  g.m = (std::size_t) INTEGER(dim)[1];
  switch (TYPEOF(gm)) {
    case REALSXP: g.d = REAL(gm); break;
    case INTSXP: g.i = INTEGER(gm); break;
    case RAWSXP: g.r = RAW(gm); break;
    default: throw br::Erro("the genotype matrix must be double, integer or raw");
  }
  if (g.n != n_ids)
    throw br::Erro(std::to_string(n_ids) + " identifiers for " + std::to_string(g.n) +
                   " genotype rows");
  return g;
}

// a media de cada marcador sobre os presentes, NaN sem nenhum (fora do GUARDA: pragma nao
// entra em argumento de macro)
static void medias_genotipos(const br::Genotipos& g, double* med) {
#ifdef _OPENMP
#pragma omp parallel num_threads(br::threads())
#endif
  {
    std::vector<double> col(g.n);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 16)
#endif
    for (long jj = 0; jj < (long) g.m; jj++) {
      g.coluna((std::size_t) jj, col.data());
      double soma = 0.0;
      std::size_t k = 0;
      for (double x : col) if (std::isfinite(x)) { soma += x; k++; }
      med[jj] = k ? soma / (double) k : R_NaN;
    }
  }
}

extern "C" {

// Devolve o pedigree em ORDEM TOPOLOGICA, com a endogamia de cada animal.
//
// A ordem importa e por isso ela sai: o A^-1 e os efeitos geneticos vem nessa ordem, e
// juntar de volta pela posicao original seria trocar animal em silencio.
// a priori de H dos genotipados e a linha do pedigree de cada um (1-based), para accuracy().
// A priori sai NOMEADA pelo id do animal (ped.ids na linha dele): a linha so vale no
// pedigree que o motor montou, e accuracy() remonta o pedigree que o usuario passa, que pode
// vir em outra ordem de linhas (a ordem topologica depende da ordem de entrada) e poria a
// diagonal de G* no animal errado. O valor volta desprotegido: vai direto para
// SET_VECTOR_ELT de uma lista ja protegida, ou para um PROTECT do chamador.
static SEXP sexp_priori_h(const std::vector<double>& v, const std::vector<std::size_t>& linha,
                          const br::Pedigree& ped) {
  SEXP x = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) v.size()));
  for (std::size_t q = 0; q < v.size(); q++) REAL(x)[q] = v[q];
  if (v.empty()) { UNPROTECT(1); return x; }
  SEXP nm = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) v.size()));
  for (std::size_t q = 0; q < v.size(); q++)
    if (q < linha.size() && linha[q] < ped.ids.size())
      SET_STRING_ELT(nm, (R_xlen_t) q, Rf_mkChar(ped.ids[linha[q]].c_str()));
  Rf_setAttrib(x, R_NamesSymbol, nm);
  UNPROTECT(2);
  return x;
}
static SEXP sexp_linha_h(const std::vector<std::size_t>& v) {
  SEXP x = Rf_allocVector(INTSXP, (R_xlen_t) v.size());
  for (std::size_t q = 0; q < v.size(); q++) INTEGER(x)[q] = (int) v[q] + 1;
  return x;
}

// H^-1 do passo unico como triplos, para os motores em R (limiar, sobrevivencia) e para quem
// quiser a matriz: o MESMO nucleo dos ajustadores em C++ (h_inversa em genomica.cpp), entao
// G*, a mistura, a APY e a Vecchia sao as mesmas.
SEXP R_h_inversa(SEXP id, SEXP pai, SEXP mae, SEXP mfx, SEXP gmx, SEXP gid, SEXP gm,
                 SEXP mistura, SEXP anucleo, SEXP vk) {
  GUARDA(
    auto i = textos(id, "id");
    auto p = textos(pai, "sire");
    auto m = textos(mae, "dam");
    br::Pedigree ped = br::constroi_pedigree(i, p, m, textos(mfx, "metafounders"),
        std::vector<double>(REAL(gmx), REAL(gmx) + XLENGTH(gmx)), mgs_de(mae));
    std::vector<double> f = br::endogamia(ped);
    br::Csc ainv = br::a_inversa(ped, f);
    auto gids = textos(gid, "genotypes");
    const br::Genotipos mg = genotipos_do_R(gm, gids.size());
    std::vector<std::string> nuc;
    if (XLENGTH(anucleo) > 0) nuc = textos(anucleo, "APY core");
    const std::size_t kv = (std::size_t) std::max(0, Rf_asInteger(vk));
    br::RelatorioG rel;
    br::Csc a = br::h_inversa(ped, ainv, gids, mg, Rf_asReal(mistura), nuc, kv, rel);

    const R_xlen_t nz = static_cast<R_xlen_t>(a.nnz());
    SEXP li = PROTECT(Rf_allocVector(INTSXP, nz));
    SEXP cj = PROTECT(Rf_allocVector(INTSXP, nz));
    SEXP vv = PROTECT(Rf_allocVector(REALSXP, nz));
    R_xlen_t k = 0;
    for (std::size_t c = 0; c < a.ncol; c++)
      for (std::size_t t = a.colptr[c]; t < a.colptr[c + 1]; t++, k++) {
        INTEGER(li)[k] = static_cast<int>(a.linha[t]) + 1;
        INTEGER(cj)[k] = static_cast<int>(c) + 1;
        REAL(vv)[k] = a.valor[t];
      }
    SEXP ids = PROTECT(Rf_allocVector(STRSXP, static_cast<R_xlen_t>(ped.ids.size())));
    for (std::size_t t = 0; t < ped.ids.size(); t++)
      SET_STRING_ELT(ids, static_cast<R_xlen_t>(t), Rf_mkChar(ped.ids[t].c_str()));

    const char* campos[] = {"i", "j", "x", "n", "id", "h_prior", "h_prior_row",
                            "n_imputed", "n_monomorphic"};
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 9));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 9));
    for (int q = 0; q < 9; q++) SET_STRING_ELT(nms, q, Rf_mkChar(campos[q]));
    SET_VECTOR_ELT(out, 0, li);
    SET_VECTOR_ELT(out, 1, cj);
    SET_VECTOR_ELT(out, 2, vv);
    SET_VECTOR_ELT(out, 3, Rf_ScalarInteger(static_cast<int>(a.ncol)));
    SET_VECTOR_ELT(out, 4, ids);
    SET_VECTOR_ELT(out, 5, sexp_priori_h(rel.diag_gstar, rel.linha_ped, ped));
    SET_VECTOR_ELT(out, 6, sexp_linha_h(rel.linha_ped));
    SET_VECTOR_ELT(out, 7, Rf_ScalarInteger(static_cast<int>(rel.n_imputados)));
    SET_VECTOR_ELT(out, 8, Rf_ScalarInteger(static_cast<int>(rel.n_monomorficos)));
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(6);
    return out;
  )
}

SEXP R_pedigree(SEXP id, SEXP pai, SEXP mae, SEXP mfx, SEXP gmx) {
  GUARDA(
    auto i = textos(id, "id");
    auto p = textos(pai, "sire");
    auto m = textos(mae, "dam");
    br::Pedigree ped = br::constroi_pedigree(i, p, m, textos(mfx, "metafounders"), std::vector<double>(REAL(gmx), REAL(gmx) + XLENGTH(gmx)), mgs_de(mae));
    std::vector<double> f = br::endogamia(ped);

    const R_xlen_t n = static_cast<R_xlen_t>(ped.ids.size());
    SEXP ids = PROTECT(Rf_allocVector(STRSXP, n));
    SEXP vp  = PROTECT(Rf_allocVector(INTSXP, n));
    SEXP vm  = PROTECT(Rf_allocVector(INTSXP, n));
    SEXP vf  = PROTECT(Rf_allocVector(REALSXP, n));
    for (R_xlen_t k = 0; k < n; k++) {
      SET_STRING_ELT(ids, k, Rf_mkChar(ped.ids[k].c_str()));
      // base 1 para o R; NA para progenitor desconhecido, que e o que o R espera ver
      INTEGER(vp)[k] = ped.pai[k] >= 0 ? static_cast<int>(ped.pai[k] + 1) : NA_INTEGER;
      INTEGER(vm)[k] = ped.mae[k] >= 0 ? static_cast<int>(ped.mae[k] + 1) : NA_INTEGER;
      REAL(vf)[k] = f[k];
    }
    const char* campos[] = {"id", "sire", "dam", "F"};
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 4));
    for (int j = 0; j < 4; j++) SET_STRING_ELT(nms, j, Rf_mkChar(campos[j]));
    SET_VECTOR_ELT(out, 0, ids);
    SET_VECTOR_ELT(out, 1, vp);
    SET_VECTOR_ELT(out, 2, vm);
    SET_VECTOR_ELT(out, 3, vf);
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(6);
    return out;
  )
}

// A^-1 no formato de triplos do triangulo inferior, pronta para virar dgCMatrix no R.
SEXP R_a_inversa(SEXP id, SEXP pai, SEXP mae, SEXP mfx, SEXP gmx) {
  GUARDA(
    auto i = textos(id, "id");
    auto p = textos(pai, "sire");
    auto m = textos(mae, "dam");
    br::Pedigree ped = br::constroi_pedigree(i, p, m, textos(mfx, "metafounders"), std::vector<double>(REAL(gmx), REAL(gmx) + XLENGTH(gmx)), mgs_de(mae));
    std::vector<double> f = br::endogamia(ped);
    br::Csc a = br::a_inversa(ped, f);

    const R_xlen_t nz = static_cast<R_xlen_t>(a.nnz());
    SEXP li = PROTECT(Rf_allocVector(INTSXP, nz));
    SEXP cj = PROTECT(Rf_allocVector(INTSXP, nz));
    SEXP vv = PROTECT(Rf_allocVector(REALSXP, nz));
    R_xlen_t k = 0;
    for (std::size_t c = 0; c < a.ncol; c++)
      for (std::size_t t = a.colptr[c]; t < a.colptr[c + 1]; t++, k++) {
        INTEGER(li)[k] = static_cast<int>(a.linha[t]) + 1;
        INTEGER(cj)[k] = static_cast<int>(c) + 1;
        REAL(vv)[k] = a.valor[t];
      }
    SEXP ids = PROTECT(Rf_allocVector(STRSXP, static_cast<R_xlen_t>(ped.ids.size())));
    for (std::size_t t = 0; t < ped.ids.size(); t++)
      SET_STRING_ELT(ids, static_cast<R_xlen_t>(t), Rf_mkChar(ped.ids[t].c_str()));

    const char* campos[] = {"i", "j", "x", "n", "id"};
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 5));
    for (int q = 0; q < 5; q++) SET_STRING_ELT(nms, q, Rf_mkChar(campos[q]));
    SET_VECTOR_ELT(out, 0, li);
    SET_VECTOR_ELT(out, 1, cj);
    SET_VECTOR_ELT(out, 2, vv);
    SET_VECTOR_ELT(out, 3, Rf_ScalarInteger(static_cast<int>(a.ncol)));
    SET_VECTOR_ELT(out, 4, ids);
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(6);
    return out;
  )
}

// A precisao aumentada da dominancia por subclasses (dominancia.cpp). Devolve os triplos do
// triangulo inferior de Q (base 1) e o que o R precisa para nomear e documentar os niveis.
SEXP R_dominancia_inversa(SEXP id, SEXP pai, SEXP mae, SEXP em_q, SEXP rota, SEXP max_pares) {
  GUARDA(
    auto i = textos(id, "id");
    auto p = textos(pai, "sire");
    auto m = textos(mae, "dam");
    br::Pedigree ped = br::constroi_pedigree(i, p, m);
    std::vector<std::string> quem = textos(em_q, "animals");
    std::vector<char> flag(ped.ids.size(), quem.empty() ? 1 : 0);
    if (!quem.empty()) {
      std::unordered_map<std::string, std::size_t> pos;
      for (std::size_t t = 0; t < ped.ids.size(); t++) pos.emplace(ped.ids[t], t);
      for (const std::string& a : quem) {
        auto it = pos.find(a);
        if (it == pos.end())
          throw br::Erro("animal '" + a + "' of animals = is not in the pedigree");
        flag[it->second] = 1;
      }
    }
    br::DominanciaInversa r = br::dominancia_inversa(ped, flag, Rf_asInteger(rota),
                                                     Rf_asReal(max_pares));
    SEXP tq = PROTECT(csc_para_R(r.q));
    const R_xlen_t na = static_cast<R_xlen_t>(r.animal.size());
    const R_xlen_t nh = static_cast<R_xlen_t>(r.par_x.size());
    SEXP ids = PROTECT(Rf_allocVector(STRSXP, static_cast<R_xlen_t>(ped.ids.size())));
    for (std::size_t t = 0; t < ped.ids.size(); t++)
      SET_STRING_ELT(ids, static_cast<R_xlen_t>(t), Rf_mkChar(ped.ids[t].c_str()));
    SEXP an = PROTECT(Rf_allocVector(INTSXP, na));
    for (R_xlen_t t = 0; t < na; t++) INTEGER(an)[t] = static_cast<int>(r.animal[t]) + 1;
    SEXP px = PROTECT(Rf_allocVector(INTSXP, nh));
    SEXP py = PROTECT(Rf_allocVector(INTSXP, nh));
    SEXP pc = PROTECT(Rf_allocVector(LGLSXP, nh));
    SEXP ph = PROTECT(Rf_allocVector(REALSXP, nh));
    for (R_xlen_t t = 0; t < nh; t++) {
      INTEGER(px)[t] = static_cast<int>(r.par_x[t]) + 1;
      INTEGER(py)[t] = static_cast<int>(r.par_y[t]) + 1;
      LOGICAL(pc)[t] = r.par_cheio[t] ? 1 : 0;
      REAL(ph)[t] = r.priori_h[t];
    }
    SEXP num = PROTECT(Rf_allocVector(REALSXP, 21));
    REAL(num)[0] = r.logdet_k;
    REAL(num)[1] = r.logdet_qhh;
    REAL(num)[2] = r.logdet_d;
    REAL(num)[3] = static_cast<double>(r.rota);
    REAL(num)[4] = static_cast<double>(r.n_sub);
    REAL(num)[5] = static_cast<double>(r.n_criados);
    REAL(num)[6] = static_cast<double>(r.n_processados);
    REAL(num)[7] = r.flops_densa;
    REAL(num)[8] = r.flops_esparsa;
    REAL(num)[9] = r.custo_densa;
    REAL(num)[10] = r.custo_esparsa;
    REAL(num)[11] = r.seg_fecho;
    REAL(num)[12] = r.seg_decisao;
    REAL(num)[13] = r.seg_montagem;
    REAL(num)[14] = static_cast<double>(r.motivo);
    REAL(num)[15] = r.trabalho_ordem;
    REAL(num)[16] = r.teto_trabalho;
    REAL(num)[17] = r.bytes_densa;
    REAL(num)[18] = r.entradas_densa;
    REAL(num)[19] = r.entradas_esparsa;
    REAL(num)[20] = r.trabalho_densa;
    Assinatura as;
    {
      const R_xlen_t nz = XLENGTH(VECTOR_ELT(tq, 2));
      const int* qi = INTEGER(VECTOR_ELT(tq, 0));
      const int* qj = INTEGER(VECTOR_ELT(tq, 1));
      const double* qx = REAL(VECTOR_ELT(tq, 2));
      as.palavra(static_cast<std::uint64_t>(r.q.ncol));
      for (R_xlen_t t = 0; t < nz; t++) as.triplo(qi[t], qj[t], qx[t]);
      as.real(r.logdet_k);
    }
    const char* campos[] = {"q", "pedigree_id", "animal_row", "pair_x", "pair_y", "pair_full",
                            "pair_prior", "num", "aborted", "signature"};
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 10));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 10));
    for (int q = 0; q < 10; q++) SET_STRING_ELT(nms, q, Rf_mkChar(campos[q]));
    SET_VECTOR_ELT(out, 0, tq);
    SET_VECTOR_ELT(out, 1, ids);
    SET_VECTOR_ELT(out, 2, an);
    SET_VECTOR_ELT(out, 3, px);
    SET_VECTOR_ELT(out, 4, py);
    SET_VECTOR_ELT(out, 5, pc);
    SET_VECTOR_ELT(out, 6, ph);
    SET_VECTOR_ELT(out, 7, num);
    SET_VECTOR_ELT(out, 8, Rf_ScalarLogical(r.fecho_abortado ? 1 : 0));
    SET_VECTOR_ELT(out, 9, Rf_mkString(as.hex().c_str()));
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(10);
    return out;
  )
}

// Inversa de uma simetrica positiva-definida, pela Cholesky em bloco.
SEXP R_inv_pd(SEXP m) {
  GUARDA(
    SEXP dim = Rf_getAttrib(m, R_DimSymbol);
    if (TYPEOF(m) != REALSXP || dim == R_NilValue || XLENGTH(dim) != 2)
      Rf_error("expected a numeric matrix");
    const int nl = INTEGER(dim)[0], nc = INTEGER(dim)[1];
    if (nl != nc) Rf_error("the matrix must be square");
    br::Densa s(static_cast<std::size_t>(nl), static_cast<std::size_t>(nc));
    // R guarda por COLUNA e a Densa por LINHA
    for (int j = 0; j < nc; j++)
      for (int i = 0; i < nl; i++)
        s.at(static_cast<std::size_t>(i), static_cast<std::size_t>(j)) =
            REAL(m)[static_cast<R_xlen_t>(j) * nl + i];
    br::Densa r = br::inv_pd(s);
    SEXP out = PROTECT(Rf_allocMatrix(REALSXP, nl, nc));
    for (int j = 0; j < nc; j++)
      for (int i = 0; i < nl; i++)
        REAL(out)[static_cast<R_xlen_t>(j) * nl + i] =
            r.at(static_cast<std::size_t>(i), static_cast<std::size_t>(j));
    UNPROTECT(1);
    return out;
  )
}

// Fatoracao esparsa completa: ordena por grau minimo, permuta, fatora, e devolve L com a
// permutacao e o log-determinante.
SEXP R_chol_esparsa(SEXP i, SEXP j, SEXP x, SEXP n, SEXP ordenar) {
  GUARDA(
    br::Csc a = csc_do_R(i, j, x, n);
    std::vector<std::size_t> perm(a.ncol);
    br::Csc L;
    if (Rf_asLogical(ordenar) == TRUE) {
      br::CacheSimbolica& cs = simbolica_do_padrao(a);
      perm = cs.perm;
      if (!br::cholesky(br::permuta_cache(a, cs), cs.sb, L))
        throw br::Erro("the matrix is not positive-definite");
    } else {
      for (std::size_t k = 0; k < a.ncol; k++) perm[k] = k;
      br::Csc pa = br::permuta_sim(a, perm);
      if (!br::cholesky(pa, br::simbolica(pa), L))
        throw br::Erro("the matrix is not positive-definite");
    }

    SEXP vperm = PROTECT(Rf_allocVector(INTSXP, static_cast<R_xlen_t>(perm.size())));
    for (std::size_t k = 0; k < perm.size(); k++)
      INTEGER(vperm)[k] = static_cast<int>(perm[k]) + 1;
    SEXP fac = PROTECT(csc_para_R(L));
    const char* campos[] = {"L", "perm", "logdet", "dense_block"};
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 4));
    for (int q = 0; q < 4; q++) SET_STRING_ELT(nms, q, Rf_mkChar(campos[q]));
    SET_VECTOR_ELT(out, 0, fac);
    SET_VECTOR_ELT(out, 1, vperm);
    SET_VECTOR_ELT(out, 2, Rf_ScalarReal(br::logdet(L)));
    SET_VECTOR_ELT(out, 3, Rf_ScalarInteger(static_cast<int>(br::bloco_denso_final(L))));
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(4);
    return out;
  )
}

// Inversa seletiva. `bloco` de 0 detecta a cauda densa; 1 forca a recorrencia pura.
SEXP R_inv_seletiva(SEXP i, SEXP j, SEXP x, SEXP n, SEXP bloco) {
  GUARDA(
    br::Csc a = csc_do_R(i, j, x, n);
    br::CacheSimbolica& cs = simbolica_do_padrao(a);
    const std::vector<std::size_t> perm = cs.perm;
    br::Csc L;
    if (!br::cholesky(br::permuta_cache(a, cs), cs.sb, L))
      throw br::Erro("the matrix is not positive-definite");
    br::SelInv z = br::inversa_seletiva(L, static_cast<std::size_t>(Rf_asInteger(bloco)));

    // devolve na numeracao ORIGINAL, senao quem le PEV le a diagonal trocada
    const R_xlen_t nz = static_cast<R_xlen_t>(z.valor.size());
    SEXP li = PROTECT(Rf_allocVector(INTSXP, nz));
    SEXP cj = PROTECT(Rf_allocVector(INTSXP, nz));
    SEXP vv = PROTECT(Rf_allocVector(REALSXP, nz));
    R_xlen_t k = 0;
    for (std::size_t c = 0; c < z.n; c++)
      for (std::size_t t = z.colptr[c]; t < z.colptr[c + 1]; t++, k++) {
        INTEGER(li)[k] = static_cast<int>(perm[z.linha[t]]) + 1;
        INTEGER(cj)[k] = static_cast<int>(perm[c]) + 1;
        REAL(vv)[k] = z.valor[t];
      }
    const char* campos[] = {"i", "j", "x", "n"};
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 4));
    for (int q = 0; q < 4; q++) SET_STRING_ELT(nms, q, Rf_mkChar(campos[q]));
    SET_VECTOR_ELT(out, 0, li);
    SET_VECTOR_ELT(out, 1, cj);
    SET_VECTOR_ELT(out, 2, vv);
    SET_VECTOR_ELT(out, 3, Rf_ScalarInteger(static_cast<int>(z.n)));
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(5);
    return out;
  )
}

SEXP R_resolve(SEXP i, SEXP j, SEXP x, SEXP n, SEXP b) {
  GUARDA(
    br::Csc a = csc_do_R(i, j, x, n);
    br::CacheSimbolica& cs = simbolica_do_padrao(a);
    const std::vector<std::size_t> perm = cs.perm;
    br::Csc L;
    if (!br::cholesky(br::permuta_cache(a, cs), cs.sb, L))
      throw br::Erro("the matrix is not positive-definite");
    const std::size_t nn = a.ncol;
    if (static_cast<std::size_t>(XLENGTH(b)) != nn) Rf_error("b with the wrong length");
    std::vector<double> pb(nn);
    for (std::size_t k = 0; k < nn; k++) pb[k] = REAL(b)[perm[k]];
    std::vector<double> px = br::resolve(L, pb);
    SEXP out = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(nn)));
    for (std::size_t k = 0; k < nn; k++) REAL(out)[perm[k]] = px[k];
    UNPROTECT(1);
    return out;
  )
}


// ---------------------------------------------------------------- o ajuste completo

// Monta Modelo + Tabela + Pedigree a partir dos objetos do R.
//
// `tdil` e a diluicao por termo, paralela a tsoc. Todo ajustador a envia (o snp_blup()
// inclusive), entao nao tem padrao: esquecer de passa-la e erro de compilacao, e nao um
// ajuste com d = 0 em silencio.
static br::Modelo modelo_do_R(SEXP alvo, SEXP tnome, SEXP tcol, SEXP tcov, SEXP test,
                              SEXP tgrp, SEXP tnest, SEXP tbase, SEXP tsoc, SEXP ausente,
                              SEXP usa_ausente, SEXP tdil, SEXP tfix = R_NilValue) {
  const R_xlen_t nt = XLENGTH(tnome);
  std::vector<br::Termo> termos;
  std::vector<std::pair<std::string, std::vector<std::string>>> grupos;
  for (R_xlen_t k = 0; k < nt; k++) {
    br::Termo t;
    t.nome = CHAR(STRING_ELT(tnome, k));
    t.coluna = CHAR(STRING_ELT(tcol, k));
    t.efeito = LOGICAL(tcov)[k] ? br::Efeito::Covariavel : br::Efeito::Classe;
    const int e = INTEGER(test)[k];
    t.estrutura = e == 3 ? br::Estrutura::Declarada
                : e == 2 ? br::Estrutura::Parentesco
                : e == 1 ? br::Estrutura::Diagonal : br::Estrutura::Fixo;
    const char* nest = CHAR(STRING_ELT(tnest, k));
    if (*nest) t.aninhado = nest;
    if (t.aninhado.rfind("mgs:", 0) == 0) { t.mgs = true; t.aninhado = t.aninhado.substr(4); }
    t.social = LOGICAL(tsoc)[k] != 0;
    if (TYPEOF(tdil) != REALSXP || XLENGTH(tdil) != nt)
      Rf_error("dilution: expected one numeric value per term");
    t.diluicao = REAL(tdil)[k];
    // base: colunas separadas por virgula. E o que faz um termo virar regressao aleatoria,
    // um termo com m coeficientes e uma covariancia m x m, nao m termos independentes.
    const char* base = CHAR(STRING_ELT(tbase, k));
    if (*base) {
      std::string b = base;
      std::size_t ini = 0;
      while (ini <= b.size()) {
        std::size_t fim = b.find(',', ini);
        if (fim == std::string::npos) fim = b.size();
        std::string col = b.substr(ini, fim - ini);
        if (!col.empty()) t.base.push_back(col);
        ini = fim + 1;
      }
    }
    const char* g = CHAR(STRING_ELT(tgrp, k));
    if (*g && t.estrutura != br::Estrutura::Fixo) {
      bool achou = false;
      for (auto& gn : grupos)
        if (gn.first == g) { gn.second.push_back(t.nome); achou = true; break; }
      if (!achou) grupos.push_back({g, {t.nome}});
    }
    termos.push_back(std::move(t));
  }
  br::Modelo m = br::monta_modelo(CHAR(STRING_ELT(alvo, 0)), std::move(termos), grupos,
                                  Rf_asLogical(usa_ausente) == TRUE, Rf_asReal(ausente));
  // kernel(..., fixed = v): PRENDE a variancia daquele termo em v. tfix chega por TERMO e
  // theta_fixo vive por POSICAO EM THETA, entao o mapeamento e feito aqui, onde os grupos
  // ja existem. Um kernel forma grupo escalar proprio, entao e a diagonal (0,0) dele.
  if (tfix != R_NilValue) {
    if (TYPEOF(tfix) != REALSXP || XLENGTH(tfix) != nt)
      Rf_error("fixed: expected one numeric value per term");
    for (R_xlen_t k = 0; k < nt; k++) {
      if (!R_finite(REAL(tfix)[k])) continue;
      const std::string nome = CHAR(STRING_ELT(tnome, k));
      for (const br::Grupo& gr : m.grupos) {
        bool meu = false;
        for (std::size_t t : gr.termos) if (m.termos[t].nome == nome) meu = true;
        if (!meu) continue;
        if (gr.dim != 1)
          Rf_error("kernel(fixed=) holds a SCALAR component: term '%s' shares a "
                   "covariance group, where holding one entry of the matrix and "
                   "estimating the rest is a different problem", nome.c_str());
        if (m.theta_fixo.empty())
          m.theta_fixo.assign(m.ntheta, std::numeric_limits<double>::quiet_NaN());
        m.theta_fixo[gr.offset] = REAL(tfix)[k];
      }
    }
  }
  return m;
}

// traits = de cada termo, so no multicaracter: um texto por termo, os nomes separados por
// virgula, vazio = todos. Casado pelo NOME do termo, que o R ja garante unico.
static void aplica_tracos(br::Modelo& m, SEXP tnome, SEXP ttra) {
  if (TYPEOF(ttra) != STRSXP || XLENGTH(ttra) != XLENGTH(tnome))
    Rf_error("traits: expected one text per term");
  for (R_xlen_t k = 0; k < XLENGTH(ttra); k++) {
    const std::string s = CHAR(STRING_ELT(ttra, k));
    if (s.empty()) continue;
    const std::string nome = CHAR(STRING_ELT(tnome, k));
    for (br::Termo& te : m.termos) {
      if (te.nome != nome) continue;
      std::size_t ini = 0;
      while (ini <= s.size()) {
        std::size_t fim = s.find(',', ini);
        if (fim == std::string::npos) fim = s.size();
        if (fim > ini) te.tracos.push_back(s.substr(ini, fim - ini));
        ini = fim + 1;
      }
    }
  }
}

// 2 AI^-1 (por linha) como matriz do R com os nomes de theta nas duas dimensoes
static SEXP sexp_vcov(const std::vector<double>& v, std::size_t n, SEXP nomes) {
  SEXP m = PROTECT(Rf_allocMatrix(REALSXP, (int) n, (int) n));
  for (std::size_t k = 0; k < n; k++)
    for (std::size_t j = 0; j < n; j++)
      REAL(m)[j * n + k] = v.size() == n * n ? v[k * n + j] : NA_REAL;
  SEXP dn = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dn, 0, Rf_duplicate(nomes));
  SET_VECTOR_ELT(dn, 1, Rf_duplicate(nomes));
  Rf_setAttrib(m, R_DimNamesSymbol, dn);
  UNPROTECT(2);
  return m;
}

// O CONSELHO DO LIMITE DE REENTRADAS sugere os univariados com as covariancias em zero, que e a
// partida padrao de model_mt() e de model_ar1(): quando o ajuste ja partiu deles, o conselho
// repetia o que foi feito (lista de correcoes, item 14). Aqui ele diz isso e da o que sobra
static void conselho_da_partida(std::string& m, const std::string& partida_de) {
  if (partida_de != "univariate fits") return;
  const std::string velho = "Try a start= inside the cone, for example the univariate fits "
                            "with the covariances at zero";
  const std::size_t i = m.find(velho);
  if (i == std::string::npos) return;
  m.replace(i, velho.size(), "This fit already started from the univariate fits with the "
            "covariances at zero. Try a start= inside the cone with that covariance nearer to "
            "zero than where the walk ended, or compare the -2logL with the fit that holds it "
            "at zero (the terms without the group) by a likelihood-ratio test");
}

// O RESIDUO POR CLASSE (model(residual_by =), lista de correcoes, item 35): cres e o codigo
// 1..K da classe de cada linha dos dados, com o atributo "levels". Vazio = um residuo so. A
// classe entra no Modelo (os componentes var(residual|classe)) antes do desenho, e no Desenho
// depois dele
static void classes_no_modelo(br::Modelo& m, SEXP cres) {
  if (XLENGTH(cres) == 0) return;
  SEXP lv = Rf_getAttrib(cres, R_LevelsSymbol);
  if (TYPEOF(lv) != STRSXP || XLENGTH(lv) < 2)
    Rf_error("residual_by: expected the codes of 2 or more classes with their levels");
  m.classes_residuo.clear();
  for (R_xlen_t k = 0; k < XLENGTH(lv); k++) m.classes_residuo.push_back(CHAR(STRING_ELT(lv, k)));
  m.ntheta = m.offset_residual + m.classes_residuo.size();
}
static void classes_no_desenho(br::Desenho& d, SEXP cres) {
  if (XLENGTH(cres) == 0) return;
  if ((std::size_t) XLENGTH(cres) != d.nlin) Rf_error("residual_by: one class per row of the data");
  const std::size_t K = d.modelo.classes_residuo.size();
  d.classe_res.assign(d.nlin, 0);
  for (std::size_t i = 0; i < d.nlin; i++) {
    const int c = INTEGER(cres)[i];
    if (c == NA_INTEGER || c < 1 || (std::size_t) c > K) {
      if (d.usa[i]) Rf_error("residual_by: row %d has no class", (int) i + 1);
      continue;
    }
    d.classe_res[i] = (std::uint32_t) (c - 1);
  }
}

static br::Tabela tabela_do_R(SEXP dados, SEXP nomes) {
  br::Tabela t;
  const R_xlen_t nc = XLENGTH(dados);
  for (R_xlen_t k = 0; k < nc; k++) {
    SEXP col = VECTOR_ELT(dados, k);
    br::Coluna c;
    if (TYPEOF(col) == STRSXP) {
      c.texto = true;
      const R_xlen_t n = XLENGTH(col);
      c.txt.reserve(static_cast<std::size_t>(n));
      for (R_xlen_t i = 0; i < n; i++) {
        SEXP e = STRING_ELT(col, i);
        // o nome da coluna e a linha: "position k" contava na lista das colunas usadas, que
        // o usuario nao ve, e nao dizia qual registro corrigir
        if (e == NA_STRING)
          Rf_error("NA in the text column '%s' at row %d: drop those rows or fill them",
                   CHAR(STRING_ELT(nomes, k)), (int) i + 1);
        c.txt.emplace_back(CHAR(e));
      }
    } else if (TYPEOF(col) == REALSXP) {
      c.num.assign(REAL(col), REAL(col) + XLENGTH(col));
    } else {
      Rf_error("column of unsupported type at position %d", (int) k + 1);
    }
    t.nomes.emplace_back(CHAR(STRING_ELT(nomes, k)));
    if (t.nlin == 0) t.nlin = c.tamanho();
    else if (t.nlin != c.tamanho())
      Rf_error("column with different length at position %d", (int) k + 1);
    t.colunas.push_back(std::move(c));
  }
  return t;
}

// Kernels declarados (kernel(id, K=)): lista PARALELA aos termos, NULL para termo comum e
// list(ids, K) para termo kernel. A matriz cruza como Densa (o R guarda por COLUNA e a
// Densa por LINHA), e a validacao de forma acontece aqui, antes de qualquer conta.
//
// A PRECISAO declarada (kernel(id, Kinv=)) cruza como list(ids, NULL, i, j, x, escala,
// logdet, assinatura, ...): os triplos base 1 do triangulo INFERIOR de K^-1 (i >= j, a
// convencao de a_inverse()), a escala de partida (0 ou NA = tirar do log|K^-1|), o log|K|
// pela forma fechada (NA = nenhum) e a assinatura dos triplos que o acompanha ("" = nenhuma).
// O log|K| so e aceito se a assinatura dos triplos recebidos for a mesma. Duplicados se
// somam, como no A^-1. Os triplos de dominance_inverse() ja vem ordenados por coluna e linha,
// e entao a Csc e montada direto, sem as copias de trabalho de de_triplos (36 bytes por
// entrada a menos no pico de um ajuste).
static std::vector<br::KernelDecl> kernels_do_R(SEXP kern) {
  std::vector<br::KernelDecl> out;
  if (Rf_isNull(kern)) return out;
  if (TYPEOF(kern) != VECSXP) Rf_error("kernels: expected a list parallel to the terms");
  const R_xlen_t nt = XLENGTH(kern);
  out.resize(static_cast<std::size_t>(nt));
  for (R_xlen_t k = 0; k < nt; k++) {
    SEXP e = VECTOR_ELT(kern, k);
    if (Rf_isNull(e)) continue;
    if (TYPEOF(e) == VECSXP && XLENGTH(e) >= 8 && Rf_isNull(VECTOR_ELT(e, 1))) {
      auto ids = textos(VECTOR_ELT(e, 0), "kernel ids");
      SEXP ti = VECTOR_ELT(e, 2), tj = VECTOR_ELT(e, 3), tx = VECTOR_ELT(e, 4);
      SEXP te = VECTOR_ELT(e, 5), tl = VECTOR_ELT(e, 6), ts = VECTOR_ELT(e, 7);
      const int termo = (int) k + 1;
      if (TYPEOF(ti) != INTSXP || TYPEOF(tj) != INTSXP || TYPEOF(tx) != REALSXP ||
          TYPEOF(te) != REALSXP || XLENGTH(te) != 1 || TYPEOF(tl) != REALSXP ||
          XLENGTH(tl) != 1 || TYPEOF(ts) != STRSXP || XLENGTH(ts) != 1)
        throw br::Erro("kernel of term " + std::to_string(termo) + ": expected list(ids, NULL, "
                       "i, j, x, scale, logdet, signature)");
      const R_xlen_t nz = XLENGTH(tx);
      if (XLENGTH(ti) != nz || XLENGTH(tj) != nz)
        throw br::Erro("kernel of term " + std::to_string(termo) + ": i, j and x of Kinv with "
                       "different lengths");
      const std::size_t nk = ids.size();
      const int* pi = INTEGER(ti);
      const int* pj = INTEGER(tj);
      const double* px = REAL(tx);
      // uma passada: confere, assina e ve se ja esta em ordem de coluna e linha sem repetir
      Assinatura as;
      as.palavra(static_cast<std::uint64_t>(nk));
      bool em_ordem = true;
      int ai = 0, aj = 0;
      for (R_xlen_t t = 0; t < nz; t++) {
        const int a = pi[t], b = pj[t];
        if (a == NA_INTEGER || b == NA_INTEGER || a < 1 || b < 1 ||
            (std::size_t) a > nk || (std::size_t) b > nk)
          throw br::Erro("kernel of term " + std::to_string(termo) + ": a Kinv triplet outside "
                         "the " + std::to_string(nk) + " x " + std::to_string(nk) + " matrix");
        if (a < b)
          throw br::Erro("kernel of term " + std::to_string(termo) + ": Kinv triplets must be "
                         "the lower triangle, i >= j (the convention of a_inverse()); entry (" +
                         std::to_string(a) + ", " + std::to_string(b) + ") is above the diagonal");
        if (!std::isfinite(px[t]))
          throw br::Erro("kernel of term " + std::to_string(termo) + ": Kinv has a non-finite "
                         "value");
        if (em_ordem && t > 0 && !(b > aj || (b == aj && a > ai))) em_ordem = false;
        ai = a;
        aj = b;
        as.triplo(a, b, px[t]);
      }
      br::KernelDecl& kd = out[(std::size_t) k];
      kd.esparsa = true;
      if (em_ordem) {
        br::Csc c(nk, nk);
        c.linha.resize(static_cast<std::size_t>(nz));
        c.valor.resize(static_cast<std::size_t>(nz));
        for (R_xlen_t t = 0; t < nz; t++) {
          c.colptr[static_cast<std::size_t>(pj[t])]++;
          c.linha[static_cast<std::size_t>(t)] = static_cast<std::uint32_t>(pi[t] - 1);
          c.valor[static_cast<std::size_t>(t)] = px[t];
        }
        for (std::size_t q = 0; q < nk; q++) c.colptr[q + 1] += c.colptr[q];
        kd.kinv = std::move(c);
      } else {
        std::vector<std::uint32_t> li, cj;
        std::vector<double> v;
        li.reserve(static_cast<std::size_t>(nz));
        cj.reserve(static_cast<std::size_t>(nz));
        v.reserve(static_cast<std::size_t>(nz));
        for (R_xlen_t t = 0; t < nz; t++) {
          li.push_back(static_cast<std::uint32_t>(pi[t] - 1));
          cj.push_back(static_cast<std::uint32_t>(pj[t] - 1));
          v.push_back(px[t]);
        }
        kd.kinv = br::de_triplos(nk, nk, li, cj, v);
      }
      const double esc = REAL(te)[0];
      kd.escala = (std::isfinite(esc) && esc > 0.0) ? esc : 0.0;
      const double ld = REAL(tl)[0];
      as.real(ld);
      if (std::isfinite(ld) && STRING_ELT(ts, 0) != NA_STRING &&
          as.hex() == CHAR(STRING_ELT(ts, 0)))
        kd.logdet_kinv = -ld;
      kd.ids = std::move(ids);
      continue;
    }
    if (TYPEOF(e) != VECSXP || XLENGTH(e) != 2)
      Rf_error("kernel of term %d: expected list(ids, K)", (int) k + 1);
    auto ids = textos(VECTOR_ELT(e, 0), "kernel ids");
    SEXP km = VECTOR_ELT(e, 1);
    SEXP dim = Rf_getAttrib(km, R_DimSymbol);
    if (TYPEOF(km) != REALSXP || dim == R_NilValue || XLENGTH(dim) != 2)
      Rf_error("kernel of term %d: K must be a numeric matrix", (int) k + 1);
    const int nl = INTEGER(dim)[0], nc = INTEGER(dim)[1];
    if (nl != nc || (std::size_t) nl != ids.size())
      Rf_error("kernel of term %d: K of %d x %d for %d ids", (int) k + 1, nl, nc,
               (int) ids.size());
    br::Densa m((std::size_t) nl, (std::size_t) nc);
    for (int j2 = 0; j2 < nc; j2++)
      for (int i2 = 0; i2 < nl; i2++)
        m.at((std::size_t) i2, (std::size_t) j2) = REAL(km)[(R_xlen_t) j2 * nl + i2];
    out[(std::size_t) k].ids = std::move(ids);
    out[(std::size_t) k].k = std::move(m);
  }
  return out;
}

extern "C++" {

// fit$dense_block = c(dense = k, columns = n): o k do custo k^3 de cada fatoracao. O objeto
// e protegido e desprotegido aqui dentro e devolvido para ser gravado na lista na hora.
static SEXP sexp_bloco_denso(std::size_t k, std::size_t n) {
  SEXP db = PROTECT(Rf_allocVector(INTSXP, 2));
  INTEGER(db)[0] = (int) k;
  INTEGER(db)[1] = (int) n;
  SEXP nm = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(nm, 0, Rf_mkChar("dense"));
  SET_STRING_ELT(nm, 1, Rf_mkChar("columns"));
  Rf_setAttrib(db, R_NamesSymbol, nm);
  UNPROTECT(2);
  return db;
}

// Os residuos do multicaracter e do AR(1) nas posicoes q0..q0+4 de `out` (ja protegido):
// residual, fitted, fitted_fixed e resid_var como matrizes nlin x t do R, com as
// caracteristicas nos nomes das colunas e NA em linha nao usada e celula ausente, e used.
static void poe_residuos(SEXP out, int q0, const br::ResiduosMT& res,
                         const std::vector<char>& usa, const std::vector<std::string>& alvos) {
  const std::size_t nlin = usa.size(), t = alvos.size();
  const std::vector<double>* v[4] = {&res.residuo, &res.ajustado, &res.ajustado_fixo,
                                     &res.variancia};
  for (int q = 0; q < 4; q++) {
    SEXP m = PROTECT(Rf_allocMatrix(REALSXP, (int) nlin, (int) t));
    const bool tem = res.ok && v[q]->size() == nlin * t;
    for (std::size_t r = 0; r < nlin; r++)
      for (std::size_t tau = 0; tau < t; tau++) {
        const double x = tem ? (*v[q])[r * t + tau] : NA_REAL;
        REAL(m)[tau * nlin + r] = std::isfinite(x) ? x : NA_REAL;
      }
    SEXP dn = PROTECT(Rf_allocVector(VECSXP, 2));
    SEXP cn = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) t));
    for (std::size_t tau = 0; tau < t; tau++)
      SET_STRING_ELT(cn, (R_xlen_t) tau, Rf_mkChar(alvos[tau].c_str()));
    SET_VECTOR_ELT(dn, 1, cn);
    Rf_setAttrib(m, R_DimNamesSymbol, dn);
    SET_VECTOR_ELT(out, q0 + q, m);
    UNPROTECT(3);
  }
  SEXP u = PROTECT(Rf_allocVector(LGLSXP, (R_xlen_t) nlin));
  for (std::size_t r = 0; r < nlin; r++) LOGICAL(u)[r] = usa[r] ? TRUE : FALSE;
  SET_VECTOR_ELT(out, q0 + 4, u);
  UNPROTECT(1);
}

template <class DES>
std::string genomica_no_desenho(DES& d, const br::Pedigree* pp, br::Pedigree& ped,
                                SEXP gid, SEXP gm, SEXP mistura, SEXP anucleo,
                                SEXP vk, std::vector<double>* dh = nullptr,
                                std::vector<std::size_t>* dr = nullptr) {
  if (XLENGTH(gid) == 0) return "";
  if (!pp) Rf_error("genotypes without a pedigree: H^-1 needs A^-1");
  // as frequencias alelicas da G vem dos proprios genotipados: com um so a G e nula, e com dois
  // ela forca um parentesco de -1 entre eles (catalogo de estrutura, prioridade media)
  if (XLENGTH(gid) < 3)
    throw br::Erro("the single step has " + std::to_string((long) XLENGTH(gid)) +
                   " genotyped animal(s): the allele frequencies of G come from the genotyped "
                   "themselves, so with fewer than 3 G is not defined (one animal) or forces a "
                   "relationship of -1 between two. Fit on the pedigree alone, or add genotypes");
  auto gids = textos(gid, "genotypes");
  // os genotipos lidos direto do objeto do R, sem copia
  const br::Genotipos mg = genotipos_do_R(gm, gids.size());
  const std::size_t nm2 = mg.m;
  std::vector<std::string> nuc;
  if (XLENGTH(anucleo) > 0) nuc = textos(anucleo, "APY core");
  const std::size_t kv = (std::size_t) std::max(0, Rf_asInteger(vk));
  br::RelatorioG rel = br::aplica_genomica(d, ped, gids, mg, Rf_asReal(mistura), nuc, kv);
  // a priori dos genotipados, para accuracy(): so quem pede recebe
  if (dh) *dh = rel.diag_gstar;
  if (dr) *dr = rel.linha_ped;
  std::string nota = "single-step: " + std::to_string(gids.size()) + " genotyped x " +
                     std::to_string(nm2) + " markers, " +
                     std::to_string(rel.n_imputados) + " missing value(s) imputed by the " +
                     "mean, " + std::to_string(rel.n_monomorficos) +
                     " monomorphic marker(s) left out of G";
  if (!nuc.empty())
    nota += "; G* inverted by APY with a core of " + std::to_string(nuc.size()) +
            " (an approximation; the core size is part of the result)";
  if (kv > 0)
    nota += "; G* inverted by the Vecchia recursion with k = " + std::to_string(kv) +
            " (an approximation; k is part of the result)";
  return nota;
}
}  // extern "C++"

// Avalia -2logL pelas DUAS vias, mais score, AI e EM num theta dado. Existe para os gates:
// a identidade MME <-> forma V, e o score contra diferencas finitas centrais.
SEXP R_avaliar(SEXP dados, SEXP nomes, SEXP alvo, SEXP tnome, SEXP tcol, SEXP tcov,
               SEXP test, SEXP tgrp, SEXP tnest, SEXP tbase, SEXP tsoc, SEXP pid, SEXP ppai,
               SEXP pmae, SEXP ausente, SEXP usa_ausente, SEXP theta, SEXP com_densa, SEXP mfx, SEXP gmx, SEXP pesos, SEXP gid, SEXP gm, SEXP mistura, SEXP anucleo, SEXP vk, SEXP kern, SEXP tdil, SEXP cres) {
  GUARDA(
    br::Modelo m = modelo_do_R(alvo, tnome, tcol, tcov, test, tgrp, tnest, tbase, tsoc, ausente, usa_ausente, tdil);
    classes_no_modelo(m, cres);
    br::Tabela t = tabela_do_R(dados, nomes);
    br::Pedigree ped;
    const br::Pedigree* pp = nullptr;
    if (XLENGTH(pid) > 0) {
      auto i = textos(pid, "id");
      auto p = textos(ppai, "sire");
      auto ma = textos(pmae, "dam");
      ped = br::constroi_pedigree(i, p, ma, textos(mfx, "metafounders"), std::vector<double>(REAL(gmx), REAL(gmx) + XLENGTH(gmx)), mgs_de(pmae));
      pp = &ped;
    }
    std::vector<double> pw;
    if (XLENGTH(pesos) > 0) pw.assign(REAL(pesos), REAL(pesos) + XLENGTH(pesos));
    std::vector<br::KernelDecl> ks = kernels_do_R(kern);
    br::Desenho d = br::monta_desenho(m, t, pp, pw.empty() ? nullptr : &pw,
                                      ks.empty() ? nullptr : &ks);
    classes_no_desenho(d, cres);
    // os kernels ja estao no desenho (a K^-1 de cada grupo vai em d.kinv): a copia de
    // entrada sai aqui, antes das MME e do fator, e nao fica no pico do ajuste
    std::vector<br::KernelDecl>().swap(ks);
    std::vector<double> th(REAL(theta), REAL(theta) + XLENGTH(theta));
    if (th.size() != m.ntheta) Rf_error("theta with %d entries; the layout asks for %d",
                                        (int) th.size(), (int) m.ntheta);

    genomica_no_desenho(d, pp, ped, gid, gm, mistura, anucleo, vk);
    br::Avaliacao a = br::avalia(d, th);
    if (!a.ok) Rf_error("theta INADMISSIBLE: some covariance is not positive-definite");
    const double dv = Rf_asLogical(com_densa) == TRUE ? br::neg2logl_densa_V(d, th)
                                                      : std::nan("");

    const std::size_t nt2 = m.ntheta;
    SEXP score = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) nt2));
    SEXP em    = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) nt2));
    SEXP ai    = PROTECT(Rf_allocMatrix(REALSXP, (int) nt2, (int) nt2));
    for (std::size_t k = 0; k < nt2; k++) {
      REAL(score)[k] = a.score[k];
      REAL(em)[k] = a.em_theta[k];
      for (std::size_t j = 0; j < nt2; j++)
        REAL(ai)[j * nt2 + k] = a.ai.at(k, j);
    }
    const char* campos[] = {"neg2logl", "neg2logl_V", "score", "em", "ai", "off_pattern"};
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 6));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 6));
    for (int q = 0; q < 6; q++) SET_STRING_ELT(nms, q, Rf_mkChar(campos[q]));
    SET_VECTOR_ELT(out, 0, Rf_ScalarReal(a.neg2logl));
    SET_VECTOR_ELT(out, 1, Rf_ScalarReal(dv));
    SET_VECTOR_ELT(out, 2, score);
    SET_VECTOR_ELT(out, 3, em);
    SET_VECTOR_ELT(out, 4, ai);
    SET_VECTOR_ELT(out, 5, Rf_ScalarInteger((int) a.fora_do_padrao));
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(5);
    return out;
  )
}

// Conversao + aplicacao da genomica, comum aos tres ajustadores (o desenho e template
// porque uni, multi e AR carregam os mesmos campos que o passo unico toca). Devolve a
// nota do relatorio, vazia sem genotipos. O extern "C++" existe porque este arquivo vive
// num bloco extern "C", e template nao tem linkage de C.
// SEM GRAUS DE LIBERDADE RESIDUAIS (achado na integracao da leva 3): com n registros (celulas de
// uma caracteristica, no multicaracter) e p colunas de efeito fixo que o desenho mantem (o posto
// de X), o REML tem n - p graus de liberdade. Com n - p <= 0 o -2logL nao e definido (as duas
// formas da forma quadratica divergiam por milhares no mesmo theta) e a caminhada parava com um
// decremento de Newton de 1e14, dizendo "did NOT converge" sem a causa. O ajuste que estima os
// componentes recusa; com os componentes dados (maxiter = 0, n_em = 0) as solucoes seguem, e a
// nota vai para a mensagem.
extern "C++" {
static std::string confere_gl_residuais(std::size_t n, std::size_t p, const std::string& traco,
                                        bool estima) {
  if (n > p) return "";
  const std::string onde = traco.empty() ? std::string() : " of trait '" + traco + "'";
  const std::string txt = "the " + std::to_string(n) + " record(s)" + onde +
      " in the fit are not more than the " + std::to_string(p) +
      " fixed-effect column(s) the design keeps: REML has no residual degree of freedom";
  if (estima)
    throw br::Erro(txt + ", so it has no information on the variance components. Remove fixed "
                   "effects, or give the components (start= with maxiter = 0)");
  return txt + ", and the -2logL is not defined; the solutions use the components given";
}
}  // extern "C++"


SEXP R_ajustar(SEXP dados, SEXP nomes, SEXP alvo, SEXP tnome, SEXP tcol, SEXP tcov,
               SEXP test, SEXP tgrp, SEXP tnest, SEXP tbase, SEXP tsoc, SEXP pid, SEXP ppai,
               SEXP pmae, SEXP ausente, SEXP usa_ausente, SEXP maxiter, SEXP tol, SEXP n_em,
               SEXP gid, SEXP gm, SEXP mistura, SEXP anucleo, SEXP vk, SEXP verb, SEXP mfx, SEXP gmx, SEXP pesos, SEXP inicio, SEXP kern, SEXP tdil, SEXP tfix, SEXP cres) {
  GUARDA(
    br::Modelo m = modelo_do_R(alvo, tnome, tcol, tcov, test, tgrp, tnest, tbase, tsoc, ausente, usa_ausente, tdil, tfix);
    classes_no_modelo(m, cres);
    br::Tabela t = tabela_do_R(dados, nomes);
    // a priori dos genotipados sob passo unico, para accuracy(). Vazio sem genomica, e
    // nesse caso accuracy() segue com o 1 + F do pedigree.
    std::vector<double> diag_h_geno;
    std::vector<std::size_t> linha_h_geno;
    br::Pedigree ped;
    const br::Pedigree* pp = nullptr;
    if (XLENGTH(pid) > 0) {
      auto i = textos(pid, "id");
      auto p = textos(ppai, "sire");
      auto ma = textos(pmae, "dam");
      ped = br::constroi_pedigree(i, p, ma, textos(mfx, "metafounders"), std::vector<double>(REAL(gmx), REAL(gmx) + XLENGTH(gmx)), mgs_de(pmae));
      pp = &ped;
    }
    std::vector<double> pw;
    if (XLENGTH(pesos) > 0) pw.assign(REAL(pesos), REAL(pesos) + XLENGTH(pesos));
    std::vector<br::KernelDecl> ks = kernels_do_R(kern);
    br::Desenho d = br::monta_desenho(m, t, pp, pw.empty() ? nullptr : &pw,
                                      ks.empty() ? nullptr : &ks);
    classes_no_desenho(d, cres);
    // os kernels ja estao no desenho (a K^-1 de cada grupo vai em d.kinv): a copia de
    // entrada sai aqui, antes das MME e do fator, e nao fica no pico do ajuste
    std::vector<br::KernelDecl>().swap(ks);
    // so quando os componentes andam: em maxiter = 0 e n_em = 0 o theta e dado
    if (Rf_asInteger(maxiter) > 0 || Rf_asInteger(n_em) > 0)
      br::confere_covariancias_iid(d.modelo, d.aleatorios, d.usa);

    std::string nota = genomica_no_desenho(d, pp, ped, gid, gm, mistura, anucleo, vk,
                                           &diag_h_geno, &linha_h_geno);
    // sem graus de liberdade residuais: recusa ao estimar, nota com os componentes dados
    {
      const std::string gl = confere_gl_residuais(d.n_usadas(), d.x.ncol, "",
                                                  Rf_asInteger(maxiter) > 0 || Rf_asInteger(n_em) > 0);
      if (!gl.empty()) nota = nota.empty() ? gl : nota + "; " + gl;
    }

    std::vector<double> th0;
    if (XLENGTH(inicio) > 0) {
      th0.assign(REAL(inicio), REAL(inicio) + XLENGTH(inicio));
      if (th0.size() != m.ntheta)
        Rf_error("start with %d entries; the layout asks for %d", (int) th0.size(), (int) m.ntheta);
    }
    br::Ajuste r = br::ajusta(d, th0.empty() ? nullptr : &th0, (std::size_t) Rf_asInteger(n_em),
                              (std::size_t) Rf_asInteger(maxiter), Rf_asReal(tol),
                              Rf_asLogical(verb) == TRUE);
    if (!nota.empty()) r.mensagem = r.mensagem.empty() ? nota : nota + "; " + r.mensagem;

    std::vector<std::string> nomes_th = m.nomes_theta();
    const std::size_t nt2 = m.ntheta;
    SEXP theta = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) nt2));
    SEXP se    = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) nt2));
    SEXP nms_t = PROTECT(Rf_allocVector(STRSXP,  (R_xlen_t) nt2));
    for (std::size_t k = 0; k < nt2; k++) {
      REAL(theta)[k] = r.theta.empty() ? NA_REAL : r.theta[k];
      REAL(se)[k] = r.se.empty() ? NA_REAL : r.se[k];
      SET_STRING_ELT(nms_t, (R_xlen_t) k, Rf_mkChar(nomes_th[k].c_str()));
    }
    Rf_setAttrib(theta, R_NamesSymbol, nms_t);
    Rf_setAttrib(se, R_NamesSymbol, nms_t);

    // EBV e PEV por grupo: fatias da solucao e da inversa seletiva, com os niveis como nomes
    //
    // Um ajuste que parou ANTES da primeira avaliacao volta sem theta: montar as MME com
    // esse vetor vazio lia fora dele. Sem theta nao ha solucao para fatiar, e o objeto ja
    // carrega a mensagem, as fatias saem NA, como no resto do caminho de falha.
    br::Montado M;
    if (!r.theta.empty()) M = br::monta_mme(d, r.theta);
    if (M.offset_grupo.empty()) M.offset_grupo.assign(m.grupos.size(), 0);
    SEXP ebv = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t) m.grupos.size()));
    SEXP pev = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t) m.grupos.size()));
    SEXP ebv_nomes = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) m.grupos.size()));
    for (std::size_t g = 0; g < m.grupos.size(); g++) {
      const std::size_t larg = d.largura(g);
      SEXP v = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) larg));
      SEXP pv = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) larg));
      for (std::size_t c = 0; c < larg; c++) {
        REAL(v)[c] = r.solucao.empty() ? NA_REAL : r.solucao[M.offset_grupo[g] + c];
        REAL(pv)[c] = r.pev.empty() ? NA_REAL : r.pev[M.offset_grupo[g] + c];
      }
      for (const br::DesenhoTermo& a : d.aleatorios)
        if (a.termo == m.grupos[g].termos[0]) {
          SEXP nv = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) larg));
          for (std::size_t c = 0; c < larg; c++)
            SET_STRING_ELT(nv, (R_xlen_t) c, Rf_mkChar(a.niveis[c % a.n_niveis].c_str()));
          Rf_setAttrib(v, R_NamesSymbol, nv);
          Rf_setAttrib(pv, R_NamesSymbol, nv);
          UNPROTECT(1);
          break;
        }
      SET_VECTOR_ELT(ebv, (R_xlen_t) g, v);
      SET_VECTOR_ELT(pev, (R_xlen_t) g, pv);
      SET_STRING_ELT(ebv_nomes, (R_xlen_t) g, Rf_mkChar(m.grupos[g].nome.c_str()));
      UNPROTECT(2);
    }
    Rf_setAttrib(ebv, R_NamesSymbol, ebv_nomes);
    Rf_setAttrib(pev, R_NamesSymbol, Rf_duplicate(ebv_nomes));

    SEXP sc = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) nt2));
    SEXP vc = PROTECT(Rf_allocMatrix(REALSXP, (int) nt2, (int) nt2));
    for (std::size_t k = 0; k < nt2; k++) {
      REAL(sc)[k] = r.score.empty() ? NA_REAL : r.score[k];
      for (std::size_t j = 0; j < nt2; j++)
        REAL(vc)[j * nt2 + k] = r.vcov.empty() ? NA_REAL : r.vcov[k * nt2 + j];
    }
    Rf_setAttrib(sc, R_NamesSymbol, Rf_duplicate(nms_t));

    // Solucoes dos efeitos fixos: as primeiras x.ncol posicoes da solucao, nomeadas pelas
    // colunas de X, os MESMOS nomes termo=nivel que dropped_x usa. Nada e recalculado: o
    // solver ja resolveu o sistema inteiro, isto e uma fatia. Sem theta nao ha solucao e a
    // fatia sai NA, como no resto do caminho de falha.
    SEXP bfix = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) d.x.ncol));
    SEXP bnms = PROTECT(Rf_allocVector(STRSXP,  (R_xlen_t) d.x.ncol));
    for (std::size_t k = 0; k < d.x.ncol; k++) {
      REAL(bfix)[k] = r.solucao.empty() ? NA_REAL : r.solucao[k];
      SET_STRING_ELT(bnms, (R_xlen_t) k, Rf_mkChar(d.nomes_x[k].c_str()));
    }
    // o motor resolve com as covariaveis padronizadas; b volta a unidade do dado
    br::desfaz_escala_b(d.volta_x[0], [](std::size_t j) { return (int) j; },
                        REAL(bfix), true);
    Rf_setAttrib(bfix, R_NamesSymbol, bnms);

    // O residuo condicional, os ajustados e a alavanca por linha da tabela, NA nas linhas
    // que nao entram. O motor os da na escala escalada por sqrt(w); aqui voltam a escala do
    // dado (a alavanca e adimensional e nao muda). Var(e^_i) = s2e (1 - h_ii) / w_i.
    if (r.alavanca_fora > 0)
      r.mensagem += std::string(r.mensagem.empty() ? "" : "; ") +
          std::to_string(r.alavanca_fora) + " leverage pair(s) read outside the sparsity "
          "pattern of the factor: the leverages are NOT reliable";
    SEXP res_s[4];
    const std::vector<double>* res_v[4] = {&r.residuo, &r.ajustado, &r.ajustado_fixo,
                                           &r.alavanca};
    // no residuo por classe as linhas tambem vieram escaladas por sqrt(s2_1 / s2_k), no ponto
    // final
    std::vector<double> w_classe;
    if (!d.classe_res.empty() && r.theta.size() == d.modelo.ntheta) {
      w_classe.assign(d.nlin, 1.0);
      for (std::size_t i = 0; i < d.nlin; i++)
        w_classe[i] = r.theta[d.modelo.offset_residual] /
                      r.theta[d.modelo.offset_residual + d.classe_res[i]];
    }
    for (int q = 0; q < 4; q++) {
      res_s[q] = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) d.nlin));
      for (std::size_t i = 0; i < d.nlin; i++) {
        double v = (res_v[q]->size() == d.nlin && d.usa[i]) ? (*res_v[q])[i] : NA_REAL;
        if (q < 3 && !pw.empty() && R_finite(v)) v /= std::sqrt(pw[i]);
        if (q < 3 && !w_classe.empty() && R_finite(v)) v /= std::sqrt(w_classe[i]);
        REAL(res_s[q])[i] = std::isfinite(v) ? v : NA_REAL;
      }
    }
    const char* campos[] = {"theta", "se", "neg2logl", "converged", "iters", "reldelta",
                            "message", "n_used", "n_columns", "ebv", "dropped_x", "pev",
                            "score", "vcov", "b", "newton_dec",
                            "h_prior", "h_prior_row", "used", "dense_block",
                            "residual", "fitted", "fitted_fixed", "leverage"};
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 24));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 24));
    for (int q = 0; q < 24; q++) SET_STRING_ELT(nms, q, Rf_mkChar(campos[q]));
    for (int q = 0; q < 4; q++) SET_VECTOR_ELT(out, 20 + q, res_s[q]);
    SEXP saiu = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) d.saiu_x.size()));
    for (std::size_t k = 0; k < d.saiu_x.size(); k++)
      SET_STRING_ELT(saiu, (R_xlen_t) k, Rf_mkChar(d.saiu_x[k].c_str()));
    SET_VECTOR_ELT(out, 0, theta);
    SET_VECTOR_ELT(out, 1, se);
    SET_VECTOR_ELT(out, 2, Rf_ScalarReal(r.neg2logl));
    SET_VECTOR_ELT(out, 3, Rf_ScalarLogical(r.convergiu));
    SET_VECTOR_ELT(out, 4, Rf_ScalarInteger((int) r.iters));
    SET_VECTOR_ELT(out, 5, Rf_ScalarReal(r.reldelta));
    SET_VECTOR_ELT(out, 6, Rf_mkString(r.mensagem.c_str()));
    SET_VECTOR_ELT(out, 7, Rf_ScalarInteger((int) d.n_usadas()));
    SET_VECTOR_ELT(out, 8, Rf_ScalarInteger((int) d.total_colunas()));
    SET_VECTOR_ELT(out, 9, ebv);
    SET_VECTOR_ELT(out, 10, saiu);
    SET_VECTOR_ELT(out, 11, pev);
    SET_VECTOR_ELT(out, 12, sc);
    SET_VECTOR_ELT(out, 13, vc);
    SET_VECTOR_ELT(out, 14, bfix);
    SET_VECTOR_ELT(out, 15, Rf_ScalarReal(r.decremento));
    // Sob passo unico a variancia a priori de um genotipado e a diagonal de H, que
    // naquele bloco e a de G*, e nao 1 + F do pedigree. Vao as duas coisas: o valor,
    // nomeado pelo id do animal (accuracy() casa por esse nome), e a LINHA do pedigree do
    // motor a que ele pertence. Vazio quando nao houve genomica, e ai accuracy() segue com o
    // F do pedigree.
    SEXP dh = PROTECT(sexp_priori_h(diag_h_geno, linha_h_geno, ped));
    SEXP dr = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) linha_h_geno.size()));
    for (std::size_t q = 0; q < linha_h_geno.size(); q++)
      INTEGER(dr)[q] = (int) linha_h_geno[q] + 1;
    SET_VECTOR_ELT(out, 16, dh);
    SET_VECTOR_ELT(out, 17, dr);
    // QUAIS linhas entraram, uma por linha da tabela. n_used sozinho conta mas nao diz
    // quais, e quem precisa somar alguma coisa NA MESMA base da verossimilhanca precisa
    // saber quais: o perfil de indirect() corrigia o jacobiano dos pesos somando sobre a
    // tabela inteira enquanto o -2logL que ele corrige e computado so sobre estas linhas.
    SEXP usadas = PROTECT(Rf_allocVector(LGLSXP, (R_xlen_t) d.nlin));
    for (std::size_t q = 0; q < d.nlin; q++) LOGICAL(usadas)[q] = d.usa[q] ? TRUE : FALSE;
    SET_VECTOR_ELT(out, 18, usadas);
    // o bloco denso final do fator e o total de colunas: o k do custo k^3 por iteracao
    SEXP db = PROTECT(Rf_allocVector(INTSXP, 2));
    INTEGER(db)[0] = (int) r.bloco_denso;
    INTEGER(db)[1] = (int) r.colunas_fator;
    SEXP dbn = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(dbn, 0, Rf_mkChar("dense"));
    SET_STRING_ELT(dbn, 1, Rf_mkChar("columns"));
    Rf_setAttrib(db, R_NamesSymbol, dbn);
    SET_VECTOR_ELT(out, 19, db);
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(22);
    return out;
  )
}

// A22^-1 pelo Schur esparso, exposta para os gates. `geno` em base 1 na ordem do pedigree
// TOPOLOGICO (a mesma de pedigree()).
SEXP R_a22_inversa(SEXP pid, SEXP ppai, SEXP pmae, SEXP geno) {
  GUARDA(
    auto i = textos(pid, "id");
    auto p = textos(ppai, "sire");
    auto ma = textos(pmae, "dam");
    br::Pedigree ped = br::constroi_pedigree(i, p, ma, {}, {}, mgs_de(pmae));
    std::vector<double> f = br::endogamia(ped);
    br::Csc ainv = br::a_inversa(ped, f);
    std::vector<std::size_t> idx;
    for (R_xlen_t k = 0; k < XLENGTH(geno); k++) {
      const int v = INTEGER(geno)[k];
      if (v < 1 || (std::size_t) v > ped.ids.size()) Rf_error("genotyped index outside the pedigree");
      idx.push_back((std::size_t) v - 1);
    }
    br::Densa a22i = br::a22_inversa(ainv, idx);
    const std::size_t n2 = idx.size();
    SEXP out = PROTECT(Rf_allocMatrix(REALSXP, (int) n2, (int) n2));
    for (std::size_t j2 = 0; j2 < n2; j2++)
      for (std::size_t i2 = 0; i2 < n2; i2++)
        REAL(out)[j2 * n2 + i2] = a22i.at(i2, j2);
    UNPROTECT(1);
    return out;
  )
}


// ---------------------------------------------------------------- multi-caracteristica

SEXP R_avaliar_mt(SEXP dados, SEXP nomes, SEXP alvos, SEXP tnome, SEXP tcol, SEXP tcov,
                  SEXP test, SEXP tgrp, SEXP tnest, SEXP tbase, SEXP tsoc, SEXP pid,
                  SEXP ppai, SEXP pmae, SEXP ausente, SEXP usa_ausente, SEXP theta,
                  SEXP com_densa, SEXP mfx, SEXP gmx, SEXP kern, SEXP tdil, SEXP ttra) {
  GUARDA(
    SEXP alvo1 = PROTECT(Rf_mkString(CHAR(STRING_ELT(alvos, 0))));
    br::Modelo m = modelo_do_R(alvo1, tnome, tcol, tcov, test, tgrp, tnest, tbase, tsoc,
                               ausente, usa_ausente, tdil);
    UNPROTECT(1);
    aplica_tracos(m, tnome, ttra);
    br::Tabela t = tabela_do_R(dados, nomes);
    br::Pedigree ped;
    const br::Pedigree* pp = nullptr;
    if (XLENGTH(pid) > 0) {
      auto i = textos(pid, "id");
      auto p2 = textos(ppai, "sire");
      auto ma = textos(pmae, "dam");
      ped = br::constroi_pedigree(i, p2, ma, textos(mfx, "metafounders"), std::vector<double>(REAL(gmx), REAL(gmx) + XLENGTH(gmx)), mgs_de(pmae));
      pp = &ped;
    }
    std::vector<std::string> alv = textos(alvos, "traits");
    std::vector<br::KernelDecl> kd = kernels_do_R(kern);
    br::DesenhoMT d = br::monta_desenho_mt(m, alv, t, pp, &kd);
    // os kernels ja estao no desenho (a K^-1 de cada grupo vai em d.kinv): a copia de
    // entrada sai aqui, antes das MME e do fator, e nao fica no pico do ajuste
    std::vector<br::KernelDecl>().swap(kd);
    std::vector<double> th(REAL(theta), REAL(theta) + XLENGTH(theta));
    if (th.size() != d.modelo.ntheta)
      Rf_error("theta with %d entries; the layout asks for %d", (int) th.size(),
               (int) d.modelo.ntheta);

    br::AvaliacaoMT a = br::avalia_mt(d, th);
    if (!a.ok) Rf_error("theta INADMISSIBLE");
    const double dv = Rf_asLogical(com_densa) == TRUE ? br::neg2logl_densa_V_mt(d, th)
                                                      : std::nan("");
    const std::size_t nt2 = d.modelo.ntheta;
    SEXP score = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) nt2));
    SEXP ai    = PROTECT(Rf_allocMatrix(REALSXP, (int) nt2, (int) nt2));
    for (std::size_t k = 0; k < nt2; k++) {
      REAL(score)[k] = a.score[k];
      for (std::size_t j = 0; j < nt2; j++) REAL(ai)[j * nt2 + k] = a.ai.at(k, j);
    }
    const char* campos[] = {"neg2logl", "neg2logl_V", "score", "ai", "off_pattern"};
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 5));
    for (int q = 0; q < 5; q++) SET_STRING_ELT(nms, q, Rf_mkChar(campos[q]));
    SET_VECTOR_ELT(out, 0, Rf_ScalarReal(a.neg2logl));
    SET_VECTOR_ELT(out, 1, Rf_ScalarReal(dv));
    SET_VECTOR_ELT(out, 2, score);
    SET_VECTOR_ELT(out, 3, ai);
    SET_VECTOR_ELT(out, 4, Rf_ScalarInteger((int) a.fora_do_padrao));
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(4);
    return out;
  )
}

SEXP R_ajustar_mt(SEXP dados, SEXP nomes, SEXP alvos, SEXP tnome, SEXP tcol, SEXP tcov,
                  SEXP test, SEXP tgrp, SEXP tnest, SEXP tbase, SEXP tsoc, SEXP pid,
                  SEXP ppai, SEXP pmae, SEXP ausente, SEXP usa_ausente, SEXP maxiter,
                  SEXP tol, SEXP gid, SEXP gm, SEXP mistura, SEXP anucleo, SEXP vk, SEXP verb, SEXP mfx, SEXP gmx, SEXP kern, SEXP start, SEXP tdil, SEXP uni, SEXP ttra) {
  GUARDA(
    SEXP alvo1 = PROTECT(Rf_mkString(CHAR(STRING_ELT(alvos, 0))));
    br::Modelo m = modelo_do_R(alvo1, tnome, tcol, tcov, test, tgrp, tnest, tbase, tsoc,
                               ausente, usa_ausente, tdil);
    UNPROTECT(1);
    aplica_tracos(m, tnome, ttra);
    br::Tabela t = tabela_do_R(dados, nomes);
    br::Pedigree ped;
    const br::Pedigree* pp = nullptr;
    if (XLENGTH(pid) > 0) {
      auto i = textos(pid, "id");
      auto p2 = textos(ppai, "sire");
      auto ma = textos(pmae, "dam");
      ped = br::constroi_pedigree(i, p2, ma, textos(mfx, "metafounders"), std::vector<double>(REAL(gmx), REAL(gmx) + XLENGTH(gmx)), mgs_de(pmae));
      pp = &ped;
    }
    std::vector<std::string> alv = textos(alvos, "traits");
    std::vector<br::KernelDecl> kd = kernels_do_R(kern);
    br::DesenhoMT d = br::monta_desenho_mt(m, alv, t, pp, &kd);
    // os kernels ja estao no desenho (a K^-1 de cada grupo vai em d.kinv): a copia de
    // entrada sai aqui, antes das MME e do fator, e nao fica no pico do ajuste
    std::vector<br::KernelDecl>().swap(kd);
    if (Rf_asInteger(maxiter) > 0)
      br::confere_covariancias_iid(d.modelo, d.aleatorios, d.usa);

    // a priori do genotipado para accuracy(): diag(G*), e nao 1 + F (so o model() a
    // exportava, e o multicaracter seguia com o F do pedigree no passo unico)
    std::vector<double> diag_h_geno;
    std::vector<std::size_t> linha_h_geno;
    std::string nota = genomica_no_desenho(d, pp, ped, gid, gm, mistura, anucleo, vk,
                                           &diag_h_geno, &linha_h_geno);
    // sem graus de liberdade residuais, por caracteristica: as celulas observadas dela contra
    // as equacoes fixas dela
    for (std::size_t tau = 0; tau < d.t; tau++) {
      std::size_t n_tau = 0, p_tau = 0;
      for (std::size_t i = 0; i < d.nlin; i++) if (d.usa[i] && d.obs[i * d.t + tau]) n_tau++;
      for (std::size_t j = 0; j < d.x.ncol; j++) if (d.eq_fixa[j * d.t + tau] >= 0) p_tau++;
      const std::string gl = confere_gl_residuais(n_tau, p_tau, d.alvos[tau],
                                                  Rf_asInteger(maxiter) > 0);
      if (!gl.empty()) nota = nota.empty() ? gl : nota + "; " + gl;
    }

    std::vector<double> th0(REAL(start), REAL(start) + XLENGTH(start));
    if (!th0.empty() && th0.size() != d.modelo.ntheta)
      Rf_error("start with %d entries; the layout asks for %d", (int) th0.size(),
               (int) d.modelo.ntheta);
    // A PARTIDA PADRAO a partir dos otimos univariados que o R ajustou (model_mt(), sem
    // start=): vazia quando algum vetor nao casa com o layout, e ai fica a automatica
    std::string partida_de = th0.empty() ? "var(y)" : "start=";
    if (th0.empty() && XLENGTH(uni) > 0) {
      th0 = br::partida_univariados_mt(d, vetores_do_R(uni));
      if (!th0.empty()) partida_de = "univariate fits";
    }
    br::AjusteMT r = br::ajusta_mt(d, th0.empty() ? nullptr : &th0,
                                   (std::size_t) Rf_asInteger(maxiter), Rf_asReal(tol),
                                  Rf_asLogical(verb) == TRUE);
    conselho_da_partida(r.mensagem, partida_de);
    if (!nota.empty()) r.mensagem = r.mensagem.empty() ? nota : nota + "; " + r.mensagem;

    std::vector<std::string> nomes_th = br::nomes_theta_do_mt(d);
    const std::size_t nt2 = d.modelo.ntheta;
    SEXP theta = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) nt2));
    SEXP se    = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) nt2));
    SEXP nms_t = PROTECT(Rf_allocVector(STRSXP,  (R_xlen_t) nt2));
    for (std::size_t k = 0; k < nt2; k++) {
      REAL(theta)[k] = r.theta.empty() ? NA_REAL : r.theta[k];
      REAL(se)[k] = r.se.empty() ? NA_REAL : r.se[k];
      SET_STRING_ELT(nms_t, (R_xlen_t) k, Rf_mkChar(nomes_th[k].c_str()));
    }
    Rf_setAttrib(theta, R_NamesSymbol, nms_t);
    Rf_setAttrib(se, R_NamesSymbol, nms_t);

    // EBV e PEV por grupo: a fatia do grupo, com nomes "nivel|caracteristica", a coluna
    // (tau * n_coef + ct) * n_niveis + nv e caracteristica-major no coeficiente
    SEXP ebv = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t) d.modelo.grupos.size()));
    SEXP pev = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t) d.modelo.grupos.size()));
    SEXP ebv_nomes = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) d.modelo.grupos.size()));
    {
      // ONDE O BLOCO ALEATORIO COMECA, e isto tem de ser n_fixa e nao x.ncol * t.
      // Desde que o posto passou a ser medido por caracter, o bloco fixo nao e mais o kron
      // cheio: o par (coluna, caracter) sem posto nao tem equacao. Com o tamanho antigo o
      // offset passa do lugar e a leitura dos EBV sai do fim de r.solucao. E corrupcao de
      // heap silenciosa: a sessao sobrevive aos primeiros ajustes e cai num deles adiante,
      // em celula que varia entre execucoes, que foi exatamente como ela apareceu.
      std::size_t off = d.n_fixa;
      for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
        const std::size_t larg = d.largura(g);
        SEXP v = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) larg));
        SEXP pv = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) larg));
        SEXP nv = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) larg));
        for (std::size_t c = 0; c < larg; c++) {
          REAL(v)[c] = r.solucao.empty() ? NA_REAL : r.solucao[off + c];
          REAL(pv)[c] = r.pev.empty() ? NA_REAL : r.pev[off + c];
        }
        std::size_t c0 = 0;
        for (std::size_t tm : d.modelo.grupos[g].termos)
          for (const br::DesenhoTermo& a : d.aleatorios)
            if (a.termo == tm) {
              for (std::size_t pt = 0; pt < a.nt(d.t); pt++)
                for (std::size_t ct = 0; ct < a.n_coef; ct++)
                  for (std::size_t l = 0; l < a.n_niveis; l++) {
                    std::string nm = a.niveis[l] + "|" + d.alvos[a.traco(pt)];
                    if (a.n_coef > 1) nm += "[" + std::to_string(ct) + "]";
                    SET_STRING_ELT(nv, (R_xlen_t) (c0 + (pt * a.n_coef + ct) * a.n_niveis + l),
                                   Rf_mkChar(nm.c_str()));
                  }
              c0 += a.z.ncol * a.nt(d.t);
            }
        Rf_setAttrib(v, R_NamesSymbol, nv);
        Rf_setAttrib(pv, R_NamesSymbol, Rf_duplicate(nv));
        SET_VECTOR_ELT(ebv, (R_xlen_t) g, v);
        SET_VECTOR_ELT(pev, (R_xlen_t) g, pv);
        SET_STRING_ELT(ebv_nomes, (R_xlen_t) g, Rf_mkChar(d.modelo.grupos[g].nome.c_str()));
        UNPROTECT(3);
        off += larg;
      }
    }
    Rf_setAttrib(ebv, R_NamesSymbol, ebv_nomes);
    Rf_setAttrib(pev, R_NamesSymbol, Rf_duplicate(ebv_nomes));

    // Solucoes dos efeitos fixos: a coluna j da caracteristica tau vive em j*t + tau (a
    // convencao do bloco fixo da montagem), nomeada "coluna|caracteristica" como o ebv.
    // b segue o MAPA, nao o produto cheio: o par (coluna, caracter) sem registro nao tem
    // equacao, entao nao tem estimativa, e sai em dropped_x com o nome dele.
    const std::size_t nb = d.n_fixa;
    SEXP bfix = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) nb));
    SEXP bnms = PROTECT(Rf_allocVector(STRSXP,  (R_xlen_t) nb));
    for (std::size_t j = 0; j < d.x.ncol; j++)
      for (std::size_t tau = 0; tau < d.t; tau++) {
        const int e = d.eq_fixa[j * d.t + tau];
        if (e < 0) continue;
        REAL(bfix)[e] = r.solucao.empty() ? NA_REAL : r.solucao[e];
        const std::string nm = d.nomes_x[j] + "|" + d.alvos[tau];
        SET_STRING_ELT(bnms, (R_xlen_t) e, Rf_mkChar(nm.c_str()));
      }
    // b de cada caracteristica na unidade do dado: o intercepto dela absorve o centro
    for (std::size_t tau = 0; tau < d.t; tau++)
      br::desfaz_escala_b(d.volta_x[tau],
                          [&](std::size_t j) { return d.eq_fixa[j * d.t + tau]; },
                          REAL(bfix), true);
    Rf_setAttrib(bfix, R_NamesSymbol, bnms);
    SEXP saiu = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) d.saiu_x.size()));
    for (std::size_t k = 0; k < d.saiu_x.size(); k++)
      SET_STRING_ELT(saiu, (R_xlen_t) k, Rf_mkChar(d.saiu_x[k].c_str()));

    // residuos no theta final: uma montagem e um fator a mais. Uma falha aqui (memoria) nao
    // leva junto o ajuste ja feito: os residuos ficam NA e a mensagem diz por que
    br::ResiduosMT res;
    if (!r.theta.empty()) {
      try {
        res = br::residuos_mt(d, r.theta, r.cache.get());
      } catch (const br::Interrompido&) {
        throw;
      } catch (const std::exception& e) {
        res = br::ResiduosMT();
        r.mensagem += std::string(r.mensagem.empty() ? "" : "; ") + "residuals not computed (" +
            e.what() + "): the fit is kept, residual and resid_var are NA";
      }
    }
    if (res.fora > 0)
      r.mensagem += std::string(r.mensagem.empty() ? "" : "; ") + std::to_string(res.fora) +
          " residual-variance pair(s) read outside the sparsity pattern of the factor: the "
          "residual variances are NOT reliable";
    const char* campos[] = {"theta", "se", "neg2logl", "converged", "iters", "reldelta",
                            "message", "n_used", "n_columns", "ebv", "pev", "b",
                            "dropped_x", "newton_dec", "dense_block", "h_prior",
                            "h_prior_row", "residual", "fitted", "fitted_fixed",
                            "resid_var", "used", "start_from", "dropped_x_reason", "start",
                            "vcov"};
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 26));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 26));
    for (int q = 0; q < 26; q++) SET_STRING_ELT(nms, q, Rf_mkChar(campos[q]));
    SET_VECTOR_ELT(out, 25, sexp_vcov(r.vcov, nt2, nms_t));
    poe_residuos(out, 17, res, d.usa, d.alvos);
    SET_VECTOR_ELT(out, 22, Rf_mkString(partida_de.c_str()));
    SET_VECTOR_ELT(out, 24, sexp_partida(r.partida, nms_t));
    {
      SEXP mot = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) d.saiu_motivo.size()));
      for (std::size_t k = 0; k < d.saiu_motivo.size(); k++)
        SET_STRING_ELT(mot, (R_xlen_t) k, Rf_mkChar(d.saiu_motivo[k].c_str()));
      SET_VECTOR_ELT(out, 23, mot);
      UNPROTECT(1);
    }
    SET_VECTOR_ELT(out, 15, sexp_priori_h(diag_h_geno, linha_h_geno, ped));
    SET_VECTOR_ELT(out, 16, sexp_linha_h(linha_h_geno));
    SET_VECTOR_ELT(out, 14, sexp_bloco_denso(r.bloco_denso, r.colunas_fator));
    SET_VECTOR_ELT(out, 0, theta);
    SET_VECTOR_ELT(out, 1, se);
    SET_VECTOR_ELT(out, 2, Rf_ScalarReal(r.neg2logl));
    SET_VECTOR_ELT(out, 3, Rf_ScalarLogical(r.convergiu));
    SET_VECTOR_ELT(out, 4, Rf_ScalarInteger((int) r.iters));
    SET_VECTOR_ELT(out, 5, Rf_ScalarReal(r.reldelta));
    SET_VECTOR_ELT(out, 6, Rf_mkString(r.mensagem.c_str()));
    SET_VECTOR_ELT(out, 7, Rf_ScalarInteger((int) d.n_usadas()));
    SET_VECTOR_ELT(out, 8, Rf_ScalarInteger((int) d.total_colunas()));
    SET_VECTOR_ELT(out, 9, ebv);
    SET_VECTOR_ELT(out, 10, pev);
    SET_VECTOR_ELT(out, 11, bfix);
    SET_VECTOR_ELT(out, 12, saiu);
    SET_VECTOR_ELT(out, 13, Rf_ScalarReal(r.decremento));
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(11);
    return out;
  )
}


// ---------------------------------------------------------------- residuo AR(1)/CAR(1)

SEXP R_avaliar_ar1(SEXP dados, SEXP nomes, SEXP alvo, SEXP tnome, SEXP tcol, SEXP tcov,
                   SEXP test, SEXP tgrp, SEXP tnest, SEXP tbase, SEXP tsoc, SEXP pid,
                   SEXP ppai, SEXP pmae, SEXP ausente, SEXP usa_ausente, SEXP sujeito,
                   SEXP tempo, SEXP theta, SEXP com_densa, SEXP mfx, SEXP gmx, SEXP kern, SEXP tdil) {
  GUARDA(
    SEXP alvo1 = PROTECT(Rf_mkString(CHAR(STRING_ELT(alvo, 0))));
    br::Modelo m = modelo_do_R(alvo1, tnome, tcol, tcov, test, tgrp, tnest, tbase, tsoc,
                               ausente, usa_ausente, tdil);
    UNPROTECT(1);
    br::Tabela t = tabela_do_R(dados, nomes);
    br::Pedigree ped;
    const br::Pedigree* pp = nullptr;
    if (XLENGTH(pid) > 0) {
      auto i = textos(pid, "id");
      auto p2 = textos(ppai, "sire");
      auto ma = textos(pmae, "dam");
      ped = br::constroi_pedigree(i, p2, ma, textos(mfx, "metafounders"), std::vector<double>(REAL(gmx), REAL(gmx) + XLENGTH(gmx)), mgs_de(pmae));
      pp = &ped;
    }
    std::vector<std::string> alv = textos(alvo, "traits");
    std::vector<br::KernelDecl> kd = kernels_do_R(kern);
    br::DesenhoAR d = br::monta_desenho_ar1(m, alv, t, pp,
                                            CHAR(STRING_ELT(sujeito, 0)),
                                            CHAR(STRING_ELT(tempo, 0)), &kd);
    // os kernels ja estao no desenho (a K^-1 de cada grupo vai em d.kinv): a copia de
    // entrada sai aqui, antes das MME e do fator, e nao fica no pico do ajuste
    std::vector<br::KernelDecl>().swap(kd);
    std::vector<double> th(REAL(theta), REAL(theta) + XLENGTH(theta));
    if (th.size() != d.modelo.ntheta)
      Rf_error("theta with %d entries; the layout asks for %d", (int) th.size(),
               (int) d.modelo.ntheta);
    // A FORMA V DENSA e a referencia dos gates: O(N^3) em tempo e N^2 em memoria, com N as
    // celulas dos registros usados. Sem pedido (NA) ela so sai ate 2000 celulas; pedida
    // (TRUE), ate 6000, e acima disso a recusa diz o tamanho ANTES de alocar. Antes o padrao
    // era calcular sempre: 7500 registros levavam 161 s so nela (medido em f2c0d69), e 20000
    // pediam 3 GB por copia de V.
    const int pede_densa = Rf_asLogical(com_densa);
    const std::size_t n_v = d.n_usadas() * d.t;
    if (pede_densa == TRUE && n_v > 6000) {
      char tam[64];
      std::snprintf(tam, sizeof tam, "%.0f", 8.0 * n_v * n_v / 1048576.0);
      throw br::Erro("with_dense = TRUE builds the dense V form over the " +
                     std::to_string(n_v) + " observed cells (a " + std::to_string(n_v) + " x " +
                     std::to_string(n_v) + " matrix, " + tam + " MB per copy, and O(n^3) "
                     "time): it is the reference of the small gates and is refused above 6000 "
                     "cells; use with_dense = FALSE, the sparse evaluation gives the same "
                     "-2logL");
    }
    br::AvaliacaoAR a = br::avalia_ar1(d, th);
    // um rho negativo no tempo continuo: a mesma recusa da partida em R_ajustar_ar1 (abaixo),
    // porque o CAR(1) tem o parametro em [0, 1) e "|rho| >= 1" nao dizia a causa
    if (!a.ok && !d.tempo_inteiro && d.offset_rho < th.size() && th[d.offset_rho] < 0.0) {
      // a mensagem em std::string: num buffer fixo um nome de coluna longo cortava a causa
      char num[32];
      std::snprintf(num, sizeof num, "%.4g", th[d.offset_rho]);
      throw br::Erro(std::string("theta INADMISSIBLE: rho(residual) = ") + num +
                     " is negative, but the time '" + CHAR(STRING_ELT(tempo, 0)) +
                     "' is continuous here (the gaps between consecutive records of a subject "
                     "are not whole numbers of one step), and on a continuous axis the "
                     "correlation rho^dt (CAR(1): corCAR1 of nlme, EXP of ASReml) exists only "
                     "for rho in [0, 1)");
    }
    if (!a.ok) Rf_error("theta INADMISSIBLE (s2e <= 0 or |rho| >= 1)");
    const bool faz_densa = pede_densa == TRUE || (pede_densa == NA_LOGICAL && n_v <= 2000);
    const double dv = faz_densa ? br::neg2logl_densa_V_ar1(d, th) : std::nan("");
    const std::size_t nt2 = d.modelo.ntheta;
    SEXP score = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) nt2));
    for (std::size_t k = 0; k < nt2; k++) REAL(score)[k] = a.score[k];
    const char* campos[] = {"neg2logl", "neg2logl_V", "score", "off_pattern", "time_step"};
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 5));
    for (int q = 0; q < 5; q++) SET_STRING_ELT(nms, q, Rf_mkChar(campos[q]));
    SET_VECTOR_ELT(out, 0, Rf_ScalarReal(a.neg2logl));
    SET_VECTOR_ELT(out, 1, Rf_ScalarReal(dv));
    SET_VECTOR_ELT(out, 2, score);
    SET_VECTOR_ELT(out, 3, Rf_ScalarInteger((int) a.fora_do_padrao));
    SET_VECTOR_ELT(out, 4, Rf_ScalarReal(std::isfinite(d.passo_tempo) ? d.passo_tempo
                                                                       : NA_REAL));
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(3);
    return out;
  )
}

SEXP R_ajustar_ar1(SEXP dados, SEXP nomes, SEXP alvo, SEXP tnome, SEXP tcol, SEXP tcov,
                   SEXP test, SEXP tgrp, SEXP tnest, SEXP tbase, SEXP tsoc, SEXP pid,
                   SEXP ppai, SEXP pmae, SEXP ausente, SEXP usa_ausente, SEXP sujeito,
                   SEXP tempo, SEXP maxiter, SEXP tol, SEXP gid, SEXP gm, SEXP mistura,
                   SEXP anucleo, SEXP vk, SEXP verb, SEXP mfx, SEXP gmx, SEXP kern, SEXP start, SEXP tdil, SEXP uni) {
  GUARDA(
    SEXP alvo1 = PROTECT(Rf_mkString(CHAR(STRING_ELT(alvo, 0))));
    br::Modelo m = modelo_do_R(alvo1, tnome, tcol, tcov, test, tgrp, tnest, tbase, tsoc,
                               ausente, usa_ausente, tdil);
    UNPROTECT(1);
    br::Tabela t = tabela_do_R(dados, nomes);
    br::Pedigree ped;
    const br::Pedigree* pp = nullptr;
    if (XLENGTH(pid) > 0) {
      auto i = textos(pid, "id");
      auto p2 = textos(ppai, "sire");
      auto ma = textos(pmae, "dam");
      ped = br::constroi_pedigree(i, p2, ma, textos(mfx, "metafounders"), std::vector<double>(REAL(gmx), REAL(gmx) + XLENGTH(gmx)), mgs_de(pmae));
      pp = &ped;
    }
    std::vector<std::string> alv = textos(alvo, "traits");
    std::vector<br::KernelDecl> kd = kernels_do_R(kern);
    br::DesenhoAR d = br::monta_desenho_ar1(m, alv, t, pp,
                                            CHAR(STRING_ELT(sujeito, 0)),
                                            CHAR(STRING_ELT(tempo, 0)), &kd);
    // os kernels ja estao no desenho (a K^-1 de cada grupo vai em d.kinv): a copia de
    // entrada sai aqui, antes das MME e do fator, e nao fica no pico do ajuste
    std::vector<br::KernelDecl>().swap(kd);
    if (Rf_asInteger(maxiter) > 0)
      br::confere_covariancias_iid(d.modelo, d.aleatorios, d.usa);

    std::vector<double> diag_h_geno;
    std::vector<std::size_t> linha_h_geno;
    std::string nota = genomica_no_desenho(d, pp, ped, gid, gm, mistura, anucleo, vk,
                                           &diag_h_geno, &linha_h_geno);
    // sem graus de liberdade residuais, por caracteristica: as celulas observadas dela contra as
    // colunas fixas (toda coluna tem equacao em toda caracteristica; a celula ausente e mv)
    for (std::size_t tau = 0; tau < d.t; tau++) {
      std::size_t n_tau = 0;
      for (std::size_t i = 0; i < d.nlin; i++)
        if (d.usa[i] && (d.obs.empty() || d.obs[i * d.t + tau])) n_tau++;
      const std::string gl = confere_gl_residuais(n_tau, d.x.ncol, d.t > 1 ? d.alvos[tau] : "",
                                                  Rf_asInteger(maxiter) > 0);
      if (!gl.empty()) nota = nota.empty() ? gl : nota + "; " + gl;
    }

    std::vector<double> th0(REAL(start), REAL(start) + XLENGTH(start));
    if (!th0.empty() && th0.size() != d.modelo.ntheta)
      Rf_error("start with %d entries; the layout asks for %d", (int) th0.size(),
               (int) d.modelo.ntheta);
    // UM RHO NEGATIVO EM start= NO TEMPO CONTINUO: rho^dt com dt nao inteiro nao e real, e o
    // ajuste voltava com theta NA, -2logL 0 e "the starting theta is INADMISSIBLE ... or rho is
    // outside (-1, 1)", que nao diz a causa. O CAR(1) (corCAR1 do nlme, EXP do ASReml) tem o
    // parametro em [0, 1): a partida e recusada dizendo isso
    if (!th0.empty() && !d.tempo_inteiro && d.offset_rho < th0.size() &&
        th0[d.offset_rho] < 0.0) {
      char num[32];
      std::snprintf(num, sizeof num, "%.4g", th0[d.offset_rho]);
      throw br::Erro(std::string("start: rho(residual) = ") + num +
                     " is negative, but the time '" + CHAR(STRING_ELT(tempo, 0)) +
                     "' is continuous here (the gaps between consecutive records of a subject "
                     "are not whole numbers of one step), and on a continuous axis the "
                     "correlation rho^dt (CAR(1): corCAR1 of nlme, EXP of ASReml) exists only "
                     "for rho in [0, 1). Give a start with rho in [0, 1), or put the times on a "
                     "grid where every gap is a whole number of one step; fitted without "
                     "start=, the check of the time axis names the records whose gaps are off "
                     "the grid of the others");
    }
    // a partida padrao a partir dos otimos univariados AR(1) que o R ajustou (cbind() sem
    // start=); vazia quando nao casa, e ai fica a automatica
    std::string partida_de = th0.empty() ? "var(y)" : "start=";
    if (th0.empty() && XLENGTH(uni) > 0) {
      // o passo de cada univariado (atributo time_step da lista, posto pelo R): o rho deles e
      // por passo, e so vira partida quando o passo e o deste desenho
      std::vector<double> passos;
      SEXP ps = Rf_getAttrib(uni, Rf_install("time_step"));
      if (TYPEOF(ps) == REALSXP)
        for (R_xlen_t k = 0; k < XLENGTH(ps); k++)
          passos.push_back(ISNA(REAL(ps)[k]) ? std::nan("") : REAL(ps)[k]);
      th0 = br::partida_univariados_ar1(d, vetores_do_R(uni), passos);
      if (!th0.empty()) partida_de = "univariate fits";
    }
    br::AjusteMT r = br::ajusta_ar1(d, th0.empty() ? nullptr : &th0,
                                    (std::size_t) Rf_asInteger(maxiter), Rf_asReal(tol),
                                   Rf_asLogical(verb) == TRUE);
    conselho_da_partida(r.mensagem, partida_de);
    if (r.espelhado) partida_de = "sign-flipped rho";
    if (!nota.empty()) r.mensagem = r.mensagem.empty() ? nota : nota + "; " + r.mensagem;
    std::vector<std::string> nomes_th = br::nomes_theta_ar1(d);
    const std::size_t nt2 = d.modelo.ntheta;
    SEXP theta = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) nt2));
    SEXP se    = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) nt2));
    SEXP nms_t = PROTECT(Rf_allocVector(STRSXP,  (R_xlen_t) nt2));
    for (std::size_t k = 0; k < nt2; k++) {
      REAL(theta)[k] = r.theta.empty() ? NA_REAL : r.theta[k];
      REAL(se)[k] = r.se.empty() ? NA_REAL : r.se[k];
      SET_STRING_ELT(nms_t, (R_xlen_t) k, Rf_mkChar(nomes_th[k].c_str()));
    }
    Rf_setAttrib(theta, R_NamesSymbol, nms_t);
    Rf_setAttrib(se, R_NamesSymbol, Rf_duplicate(nms_t));
    // EBV e PEV por grupo, como no uni: fatias com os niveis como nomes
    SEXP ebv = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t) d.modelo.grupos.size()));
    SEXP pev = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t) d.modelo.grupos.size()));
    SEXP ebv_nomes = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) d.modelo.grupos.size()));
    {
      // O BLOCO DOS GRUPOS COMECA EM x.ncol * t, como em ar1.cpp (M.n_fixo). Com x.ncol so,
      // o cbind() devolvia EBV e PEV deslocados de x.ncol * (t - 1) posicoes: os primeiros
      // valores eram solucoes FIXAS e os ultimos niveis se perdiam (medido contra a MME
      // densa: 0,68 de diferenca maxima). So estava certo com t = 1.
      std::size_t off = d.x.ncol * d.t;
      for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
        const std::size_t larg = d.largura(g);
        SEXP v = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) larg));
        SEXP pv = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) larg));
        for (std::size_t c = 0; c < larg; c++) {
          REAL(v)[c] = r.solucao.empty() ? NA_REAL : r.solucao[off + c];
          REAL(pv)[c] = r.pev.empty() ? NA_REAL : r.pev[off + c];
        }
        {
          SEXP nv = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) larg));
          std::size_t c0 = 0;
          for (std::size_t tm : d.modelo.grupos[g].termos)
            for (const br::DesenhoTermo& a : d.aleatorios)
              if (a.termo == tm) {
                for (std::size_t tau = 0; tau < d.t; tau++)
                  for (std::size_t ct = 0; ct < a.n_coef; ct++)
                    for (std::size_t l = 0; l < a.n_niveis; l++) {
                      std::string nm = a.niveis[l];
                      if (d.t > 1) nm += "|" + d.alvos[tau];
                      if (a.n_coef > 1) nm += "[" + std::to_string(ct) + "]";
                      SET_STRING_ELT(nv,
                          (R_xlen_t) (c0 + (tau * a.n_coef + ct) * a.n_niveis + l),
                          Rf_mkChar(nm.c_str()));
                    }
                c0 += a.z.ncol * d.t;
              }
          Rf_setAttrib(v, R_NamesSymbol, nv);
          Rf_setAttrib(pv, R_NamesSymbol, Rf_duplicate(nv));
          UNPROTECT(1);
        }
        SET_VECTOR_ELT(ebv, (R_xlen_t) g, v);
        SET_VECTOR_ELT(pev, (R_xlen_t) g, pv);
        SET_STRING_ELT(ebv_nomes, (R_xlen_t) g, Rf_mkChar(d.modelo.grupos[g].nome.c_str()));
        UNPROTECT(2);
        off += larg;
      }
    }
    Rf_setAttrib(ebv, R_NamesSymbol, ebv_nomes);
    Rf_setAttrib(pev, R_NamesSymbol, Rf_duplicate(ebv_nomes));

    // Solucoes dos efeitos fixos, na convencao do bloco fixo (identica a da multi): a
    // coluna j da caracteristica tau em j*t + tau. Com t = 1 o nome e so termo=nivel.
    const std::size_t nb = d.x.ncol * d.t;
    SEXP bfix = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) nb));
    SEXP bnms = PROTECT(Rf_allocVector(STRSXP,  (R_xlen_t) nb));
    for (std::size_t j = 0; j < d.x.ncol; j++)
      for (std::size_t tau = 0; tau < d.t; tau++) {
        const std::size_t k = j * d.t + tau;
        REAL(bfix)[k] = r.solucao.empty() ? NA_REAL : r.solucao[k];
        std::string nm = d.nomes_x[j];
        if (d.t > 1) nm += "|" + d.alvos[tau];
        SET_STRING_ELT(bnms, (R_xlen_t) k, Rf_mkChar(nm.c_str()));
      }
    for (std::size_t tau = 0; tau < d.t; tau++)
      br::desfaz_escala_b(d.volta_x[tau],
                          [&](std::size_t j) { return (int) (j * d.t + tau); },
                          REAL(bfix), true);
    Rf_setAttrib(bfix, R_NamesSymbol, bnms);
    SEXP saiu = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) d.saiu_x.size()));
    for (std::size_t k = 0; k < d.saiu_x.size(); k++)
      SET_STRING_ELT(saiu, (R_xlen_t) k, Rf_mkChar(d.saiu_x[k].c_str()));

    // residuos no theta final (uma montagem e um fator a mais), e o sujeito e o tempo de cada
    // linha usada, na ordem em que o motor os encadeia: o que as inovacoes precisam
    br::ResiduosMT res;
    if (!r.theta.empty()) {
      try {
        res = br::residuos_ar1(d, r.theta);
      } catch (const br::Interrompido&) {
        throw;
      } catch (const std::exception& e) {
        res = br::ResiduosMT();
        r.mensagem += std::string(r.mensagem.empty() ? "" : "; ") + "residuals not computed (" +
            e.what() + "): the fit is kept, residual and resid_var are NA";
      }
    }
    if (res.fora > 0)
      r.mensagem += std::string(r.mensagem.empty() ? "" : "; ") + std::to_string(res.fora) +
          " residual-variance pair(s) read outside the sparsity pattern of the factor: the "
          "residual variances are NOT reliable";
    const char* campos[] = {"theta", "se", "neg2logl", "converged", "iters", "reldelta",
                            "message", "n_used", "n_columns", "n_subjects", "ebv", "pev",
                            "b", "dropped_x", "newton_dec", "dense_block", "h_prior",
                            "h_prior_row", "residual", "fitted", "fitted_fixed",
                            "resid_var", "used", "subject_index", "time", "start_from",
                            "start", "dropped_x_reason", "time_step", "vcov"};
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 30));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 30));
    for (int q = 0; q < 30; q++) SET_STRING_ELT(nms, q, Rf_mkChar(campos[q]));
    SET_VECTOR_ELT(out, 29, sexp_vcov(r.vcov, nt2, nms_t));
    poe_residuos(out, 18, res, d.usa, d.alvos);
    {
      SEXP sj = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) d.nlin));
      SEXP tp = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) d.nlin));
      for (std::size_t i = 0; i < d.nlin; i++) {
        INTEGER(sj)[i] = NA_INTEGER;
        REAL(tp)[i] = NA_REAL;
      }
      // o tempo como veio no dado; o motor anda em passos (time_step) quando ha grade
      for (std::size_t si = 0; si < d.sujeitos.size(); si++)
        for (std::size_t i : d.sujeitos[si]) {
          INTEGER(sj)[i] = (int) si + 1;
          REAL(tp)[i] = d.tempo_dado[i];
        }
      SET_VECTOR_ELT(out, 23, sj);
      SET_VECTOR_ELT(out, 24, tp);
      UNPROTECT(2);
    }
    SET_VECTOR_ELT(out, 28, Rf_ScalarReal(std::isfinite(d.passo_tempo) ? d.passo_tempo
                                                                        : NA_REAL));
    SET_VECTOR_ELT(out, 25, Rf_mkString(partida_de.c_str()));
    SET_VECTOR_ELT(out, 26, sexp_partida(r.partida, nms_t));
    {
      SEXP mot = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) d.saiu_motivo.size()));
      for (std::size_t k = 0; k < d.saiu_motivo.size(); k++)
        SET_STRING_ELT(mot, (R_xlen_t) k, Rf_mkChar(d.saiu_motivo[k].c_str()));
      SET_VECTOR_ELT(out, 27, mot);
      UNPROTECT(1);
    }
    SET_VECTOR_ELT(out, 16, sexp_priori_h(diag_h_geno, linha_h_geno, ped));
    SET_VECTOR_ELT(out, 17, sexp_linha_h(linha_h_geno));
    SET_VECTOR_ELT(out, 15, sexp_bloco_denso(r.bloco_denso, r.colunas_fator));
    SET_VECTOR_ELT(out, 0, theta);
    SET_VECTOR_ELT(out, 1, se);
    SET_VECTOR_ELT(out, 2, Rf_ScalarReal(r.neg2logl));
    SET_VECTOR_ELT(out, 3, Rf_ScalarLogical(r.convergiu));
    SET_VECTOR_ELT(out, 4, Rf_ScalarInteger((int) r.iters));
    SET_VECTOR_ELT(out, 5, Rf_ScalarReal(r.reldelta));
    SET_VECTOR_ELT(out, 6, Rf_mkString(r.mensagem.c_str()));
    SET_VECTOR_ELT(out, 7, Rf_ScalarInteger((int) d.n_usadas()));
    SET_VECTOR_ELT(out, 8, Rf_ScalarInteger((int) d.total_colunas()));
    SET_VECTOR_ELT(out, 9, Rf_ScalarInteger((int) d.sujeitos.size()));
    SET_VECTOR_ELT(out, 10, ebv);
    SET_VECTOR_ELT(out, 11, pev);
    SET_VECTOR_ELT(out, 12, bfix);
    SET_VECTOR_ELT(out, 13, saiu);
    SET_VECTOR_ELT(out, 14, Rf_ScalarReal(r.decremento));
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(11);
    return out;
  )
}

// O bootstrap parametrico preparado vive num ponteiro externo do R entre os lotes. O
// finalizador o apaga quando o R o recolhe, e R_residuos_libera o apaga antes, no on.exit de
// quem chamou: o desenho e o fator nao ficam no pico do R alem do necessario.
extern "C++" {
static void finaliza_preparo(SEXP p) {
  br::PreparoNulos* P = static_cast<br::PreparoNulos*>(R_ExternalPtrAddr(p));
  if (P) delete P;
  R_ClearExternalPtr(p);
}

static br::PreparoNulos* preparo_de(SEXP p) {
  if (TYPEOF(p) != EXTPTRSXP || !R_ExternalPtrAddr(p))
    throw br::Erro("the prepared bootstrap was already released");
  return static_cast<br::PreparoNulos*>(R_ExternalPtrAddr(p));
}

// sqrt(w) das linhas usadas (vazio quando nao ha pesos)
static std::vector<double> raizes(const std::vector<char>& usa, const std::vector<double>& pw) {
  std::vector<double> raiz;
  if (pw.empty()) return raiz;
  for (std::size_t i = 0; i < usa.size(); i++)
    if (usa[i]) raiz.push_back(std::sqrt(pw[i]));
  return raiz;
}

// As replicas no formato do R: residual e fixed, matrizes n_usadas x k (as linhas usadas na
// ordem da tabela, de volta a escala do dado), e u (colunas aleatorias x k, ou NULL).
static SEXP sexp_lote(const br::NulosSaida& S, const std::vector<double>& raiz) {
  const std::size_t nu = S.n_usadas, ns = S.nsim;
  SEXP er = PROTECT(Rf_allocMatrix(REALSXP, (int) nu, (int) ns));
  SEXP fx = PROTECT(Rf_allocMatrix(REALSXP, (int) nu, (int) ns));
  for (std::size_t r = 0; r < ns; r++)
    for (std::size_t k = 0; k < nu; k++) {
      const double q = raiz.empty() ? 1.0 : raiz[k];
      REAL(er)[r * nu + k] = S.residuo[r * nu + k] / q;
      REAL(fx)[r * nu + k] = S.fixo[r * nu + k] / q;
    }
  SEXP uu = R_NilValue;
  if (!S.u.empty()) {
    const std::size_t nq = S.u.size() / std::max<std::size_t>(ns, 1);
    uu = PROTECT(Rf_allocMatrix(REALSXP, (int) nq, (int) ns));
    std::copy(S.u.begin(), S.u.end(), REAL(uu));
  } else {
    PROTECT(uu);
  }
  const char* campos[] = {"residual", "fixed", "u"};
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, 3));
  for (int q = 0; q < 3; q++) SET_STRING_ELT(nms, q, Rf_mkChar(campos[q]));
  SET_VECTOR_ELT(out, 0, er);
  SET_VECTOR_ELT(out, 1, fx);
  SET_VECTOR_ELT(out, 2, uu);
  Rf_setAttrib(out, R_NamesSymbol, nms);
  UNPROTECT(5);
  return out;
}

// O que o preparo devolve ao R: o ponteiro (ou NULL), used, o residuo dos dados no fator
// (ou NULL), o ajustado W loc (ou NULL), e por grupo os niveis e a diagonal de K (NULL no
// grupo iid), nomeados pelos grupos.
static SEXP sexp_preparo(SEXP ptr, const std::vector<char>& usa, const std::vector<double>* res,
                         const std::vector<double>* ajust, const std::vector<double>& raiz,
                         const std::vector<std::vector<double>>& diag_k,
                         const std::vector<std::vector<std::string>>& niveis,
                         const std::vector<std::string>& grupos) {
  SEXP us = PROTECT(Rf_allocVector(LGLSXP, (R_xlen_t) usa.size()));
  for (std::size_t i = 0; i < usa.size(); i++) LOGICAL(us)[i] = usa[i] ? TRUE : FALSE;
  SEXP vs[2] = {R_NilValue, R_NilValue};
  const std::vector<double>* fontes[2] = {res, ajust};
  for (int q = 0; q < 2; q++) {
    if (fontes[q]) {
      vs[q] = Rf_allocVector(REALSXP, (R_xlen_t) fontes[q]->size());
      for (std::size_t k = 0; k < fontes[q]->size(); k++)
        REAL(vs[q])[k] = (*fontes[q])[k] / (raiz.empty() ? 1.0 : raiz[k]);
    }
    PROTECT(vs[q]);
  }
  SEXP kd = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t) grupos.size()));
  SEXP nv = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t) grupos.size()));
  SEXP gn = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) grupos.size()));
  for (std::size_t g = 0; g < grupos.size(); g++) {
    SET_STRING_ELT(gn, (R_xlen_t) g, Rf_mkChar(grupos[g].c_str()));
    SEXP rot = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) niveis[g].size()));
    for (std::size_t k = 0; k < niveis[g].size(); k++)
      SET_STRING_ELT(rot, (R_xlen_t) k, Rf_mkChar(niveis[g][k].c_str()));
    SET_VECTOR_ELT(nv, (R_xlen_t) g, rot);
    if (g < diag_k.size() && !diag_k[g].empty()) {
      SEXP v = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) diag_k[g].size()));
      std::copy(diag_k[g].begin(), diag_k[g].end(), REAL(v));
      if (diag_k[g].size() == niveis[g].size()) Rf_setAttrib(v, R_NamesSymbol, rot);
      SET_VECTOR_ELT(kd, (R_xlen_t) g, v);
      UNPROTECT(1);
    }
    UNPROTECT(1);
  }
  Rf_setAttrib(kd, R_NamesSymbol, gn);
  Rf_setAttrib(nv, R_NamesSymbol, gn);
  const char* campos[] = {"ptr", "used", "residual", "fitted", "kdiag", "levels"};
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 6));
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, 6));
  for (int q = 0; q < 6; q++) SET_STRING_ELT(nms, q, Rf_mkChar(campos[q]));
  SET_VECTOR_ELT(out, 0, ptr);
  SET_VECTOR_ELT(out, 1, us);
  SET_VECTOR_ELT(out, 2, vs[0]);
  SET_VECTOR_ELT(out, 3, vs[1]);
  SET_VECTOR_ELT(out, 4, kd);
  SET_VECTOR_ELT(out, 5, nv);
  Rf_setAttrib(out, R_NamesSymbol, nms);
  UNPROTECT(8);
  return out;
}

// b do ajuste contra as colunas fixas do desenho remontado: o mesmo nome, na mesma ordem
static std::vector<double> b_conferido(SEXP b, SEXP bnomes,
                                       const std::vector<std::string>& nomes_x) {
  const std::vector<std::string> bn = textos(bnomes, "names of b");
  if (TYPEOF(b) != REALSXP || (std::size_t) XLENGTH(b) != nomes_x.size() ||
      bn.size() != nomes_x.size())
    throw br::Erro("the fixed-effect solutions of the fit have " +
                   std::to_string(XLENGTH(b)) + " entries and the rebuilt design " +
                   std::to_string(nomes_x.size()) + " columns: the data, the missing code or "
                   "the weights are not the ones of the fit");
  for (std::size_t k = 0; k < bn.size(); k++)
    if (bn[k] != nomes_x[k])
      throw br::Erro("fixed column " + std::to_string(k + 1) + " is '" + nomes_x[k] +
                     "' in the rebuilt design and '" + bn[k] + "' in the fit: the data are "
                     "not the ones of the fit");
  return std::vector<double>(REAL(b), REAL(b) + XLENGTH(b));
}

template <class DES>
static std::vector<std::string> nomes_grupos(const DES& d) {
  std::vector<std::string> out;
  for (const br::Grupo& g : d.modelo.grupos) out.push_back(g.nome);
  return out;
}

// o preparo pronto vira um ponteiro externo com finalizador, e a lista do R
static SEXP entrega_preparo(std::unique_ptr<br::PreparoNulos> P,
                            const std::vector<std::string>& grupos) {
  SEXP ptr = PROTECT(R_MakeExternalPtr(nullptr, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(ptr, finaliza_preparo, TRUE);
  br::PreparoNulos* bruto = P.release();
  R_SetExternalPtrAddr(ptr, bruto);
  SEXP out = sexp_preparo(ptr, bruto->usa, &bruto->residuo_dados, nullptr, bruto->raiz,
                          bruto->diag_k, bruto->niveis, grupos);
  UNPROTECT(1);
  return out;
}
}  // extern "C++"

// O preparo do bootstrap parametrico do model() com theta fixo (os mesmos argumentos de
// desenho do R_ajustar, mais theta, b com os nomes e loc). Com theta: monta e fatora as MME,
// devolve o ponteiro e o residuo dos dados no fator. Sem theta (a conferencia de uma cadeia
// do gibbs()): so o desenho, devolvendo W loc nas linhas usadas, sem fator das MME. Nos dois,
// a diagonal de K_g de cada grupo. loc_padronizado (logico) diz que as x.ncol primeiras
// posicoes de loc ja estao na escala padronizada do desenho (o b_standardized que o gibbs()
// guarda) e entram sem a volta de desfaz_escala_b; sem ele, loc traz b na unidade do dado.
SEXP R_residuos_prepara(SEXP dados, SEXP nomes, SEXP alvo, SEXP tnome, SEXP tcol, SEXP tcov,
                        SEXP test, SEXP tgrp, SEXP tnest, SEXP tbase, SEXP tsoc, SEXP pid,
                        SEXP ppai, SEXP pmae, SEXP ausente, SEXP usa_ausente, SEXP gid, SEXP gm,
                        SEXP mistura, SEXP anucleo, SEXP vk, SEXP mfx, SEXP gmx, SEXP pesos,
                        SEXP kern, SEXP tdil, SEXP tfix, SEXP theta, SEXP b, SEXP bnomes,
                        SEXP loc, SEXP loc_padronizado) {
  GUARDA(
    br::Modelo m = modelo_do_R(alvo, tnome, tcol, tcov, test, tgrp, tnest, tbase, tsoc,
                               ausente, usa_ausente, tdil, tfix);
    br::Tabela t = tabela_do_R(dados, nomes);
    br::Pedigree ped;
    const br::Pedigree* pp = nullptr;
    if (XLENGTH(pid) > 0) {
      auto i = textos(pid, "id");
      auto p = textos(ppai, "sire");
      auto ma = textos(pmae, "dam");
      ped = br::constroi_pedigree(i, p, ma, textos(mfx, "metafounders"), std::vector<double>(REAL(gmx), REAL(gmx) + XLENGTH(gmx)), mgs_de(pmae));
      pp = &ped;
    }
    std::vector<double> pw;
    if (XLENGTH(pesos) > 0) pw.assign(REAL(pesos), REAL(pesos) + XLENGTH(pesos));
    std::vector<br::KernelDecl> ks = kernels_do_R(kern);
    br::Desenho d = br::monta_desenho(m, t, pp, pw.empty() ? nullptr : &pw,
                                      ks.empty() ? nullptr : &ks);
    std::vector<br::KernelDecl>().swap(ks);
    genomica_no_desenho(d, pp, ped, gid, gm, mistura, anucleo, vk);
    std::vector<double> bv = b_conferido(b, bnomes, d.nomes_x);
    br::desfaz_escala_b(d.volta_x[0], [](std::size_t j) { return (int) j; },
                        bv.data(), false);
    const std::vector<std::string> grupos = nomes_grupos(d);
    if (XLENGTH(theta) == 0) {
      // a conferencia da cadeia: W loc com o loc dela, e a diagonal de cada K
      if (TYPEOF(loc) != REALSXP || (std::size_t) XLENGTH(loc) != d.total_colunas())
        throw br::Erro("the posterior means of the chain have " + std::to_string(XLENGTH(loc)) +
                       " entries and the rebuilt design " + std::to_string(d.total_colunas()) +
                       " columns: the data are not the ones of the chain");
      std::vector<std::size_t> off(m.grupos.size());
      std::size_t acc = d.x.ncol;
      for (std::size_t g = 0; g < m.grupos.size(); g++) { off[g] = acc; acc += d.largura(g); }
      const br::LinhasW W = br::linhas_w(d, off);
      std::vector<double> lv(REAL(loc), REAL(loc) + XLENGTH(loc));
      // b na unidade do dado volta a escala padronizada aqui, e com uma covariavel muito
      // deslocada b0 + m b_x se cancela e perde eps |b0| (idade + 1,7e12: 3,27e-6 no W loc).
      // O b padronizado da cadeia, o mesmo E[loc] do W E[loc] que ela guardou, entra como
      // esta, e W loc so difere do ajustado guardado no arredondamento da soma
      if (Rf_asLogical(loc_padronizado) != TRUE)
        br::desfaz_escala_b(d.volta_x[0], [](std::size_t j) { return (int) j; },
                            lv.data(), false);
      std::vector<double> ajust;
      for (std::size_t i = 0; i < d.nlin; i++) {
        if (!d.usa[i]) continue;
        double v = 0.0;
        for (std::size_t p = W.ptr[i]; p < W.ptr[i + 1]; p++) v += W.val[p] * lv[W.col[p]];
        ajust.push_back(v);
      }
      const std::vector<std::vector<double>> dk = br::diagonais_k(d.kinv);
      return sexp_preparo(R_NilValue, d.usa, nullptr, &ajust, raizes(d.usa, pw), dk,
                          br::niveis_grupos(d), grupos);
    }
    std::vector<double> th(REAL(theta), REAL(theta) + XLENGTH(theta));
    if (th.size() != m.ntheta)
      throw br::Erro("theta with " + std::to_string(th.size()) + " entries; the layout asks for " +
                     std::to_string(m.ntheta));
    std::unique_ptr<br::PreparoNulos> P = br::prepara_nulos(std::move(d), th, bv);
    P->raiz = raizes(P->usa, pw);
    return entrega_preparo(std::move(P), grupos);
  )
}

// O mesmo para o model_ar1() univariado (unidades absolutas, sem pesos).
SEXP R_residuos_prepara_ar1(SEXP dados, SEXP nomes, SEXP alvo, SEXP tnome, SEXP tcol,
                            SEXP tcov, SEXP test, SEXP tgrp, SEXP tnest, SEXP tbase, SEXP tsoc,
                            SEXP pid, SEXP ppai, SEXP pmae, SEXP ausente, SEXP usa_ausente,
                            SEXP sujeito, SEXP tempo, SEXP gid, SEXP gm, SEXP mistura,
                            SEXP anucleo, SEXP vk, SEXP mfx, SEXP gmx, SEXP kern, SEXP tdil,
                            SEXP theta, SEXP b, SEXP bnomes) {
  GUARDA(
    SEXP alvo1 = PROTECT(Rf_mkString(CHAR(STRING_ELT(alvo, 0))));
    br::Modelo m = modelo_do_R(alvo1, tnome, tcol, tcov, test, tgrp, tnest, tbase, tsoc,
                               ausente, usa_ausente, tdil);
    UNPROTECT(1);
    br::Tabela t = tabela_do_R(dados, nomes);
    br::Pedigree ped;
    const br::Pedigree* pp = nullptr;
    if (XLENGTH(pid) > 0) {
      auto i = textos(pid, "id");
      auto p2 = textos(ppai, "sire");
      auto ma = textos(pmae, "dam");
      ped = br::constroi_pedigree(i, p2, ma, textos(mfx, "metafounders"), std::vector<double>(REAL(gmx), REAL(gmx) + XLENGTH(gmx)), mgs_de(pmae));
      pp = &ped;
    }
    std::vector<std::string> alv = textos(alvo, "traits");
    std::vector<br::KernelDecl> kd = kernels_do_R(kern);
    br::DesenhoAR d = br::monta_desenho_ar1(m, alv, t, pp, CHAR(STRING_ELT(sujeito, 0)),
                                            CHAR(STRING_ELT(tempo, 0)), &kd);
    std::vector<br::KernelDecl>().swap(kd);
    genomica_no_desenho(d, pp, ped, gid, gm, mistura, anucleo, vk);
    std::vector<double> th(REAL(theta), REAL(theta) + XLENGTH(theta));
    if (th.size() != d.modelo.ntheta)
      throw br::Erro("theta with " + std::to_string(th.size()) + " entries; the layout asks for " +
                     std::to_string(d.modelo.ntheta));
    std::vector<double> bv = b_conferido(b, bnomes, d.nomes_x);
    br::desfaz_escala_b(d.volta_x[0], [](std::size_t j) { return (int) j; },
                        bv.data(), false);
    const std::vector<std::string> grupos = nomes_grupos(d);
    return entrega_preparo(br::prepara_nulos_ar1(std::move(d), th, bv), grupos);
  )
}

// Um lote de k replicas do bootstrap preparado: residual e fixed (n_usadas x k) e, com
// com_u, os u* (colunas aleatorias x k). O RNG do R corre na mesma ordem de sempre.
SEXP R_residuos_lote(SEXP ptr, SEXP k, SEXP com_u, SEXP familia) {
  GUARDA(
    br::PreparoNulos* P = preparo_de(ptr);
    const br::NulosSaida S = P->lote((std::size_t) Rf_asInteger(k), Rf_asLogical(com_u) == TRUE,
                                     Rf_asReal(familia));
    return sexp_lote(S, P->raiz);
  )
}

SEXP R_residuos_libera(SEXP ptr) {
  if (TYPEOF(ptr) == EXTPTRSXP) finaliza_preparo(ptr);
  return R_NilValue;
}

SEXP R_gibbs(SEXP dados, SEXP nomes, SEXP alvo, SEXP tnome, SEXP tcol, SEXP tcov,
             SEXP test, SEXP tgrp, SEXP tnest, SEXP tbase, SEXP tsoc, SEXP pid, SEXP ppai,
             SEXP pmae, SEXP ausente, SEXP usa_ausente, SEXP gid, SEXP gm, SEXP mistura,
             SEXP anucleo, SEXP n_iter, SEXP burnin, SEXP thin, SEXP loc_fixa,
             SEXP theta_fixo, SEXP vk, SEXP verb, SEXP mfx, SEXP gmx, SEXP kern,
             SEXP tfix, SEXP ptipo, SEXP pdf, SEXP pesc, SEXP fam, SEXP tdil,
             SEXP npdm, SEXP pnomes, SEXP pmult, SEXP rtipo, SEXP rdf, SEXP resc,
             SEXP dpfix, SEXP pxda, SEXP part0, SEXP part0_nomes, SEXP asis_r,
             SEXP etipo, SEXP edf, SEXP eesc, SEXP asis_grupos, SEXP recuo) {
  GUARDA(
    // kernel() e kernel(fixed =) chegam como no model(): o amostrador ja era generico em
    // K^-1 (o bloco de localizacao e a InvWishart leem d.kinv[g]); faltava so o transporte
    br::Modelo m = modelo_do_R(alvo, tnome, tcol, tcov, test, tgrp, tnest, tbase, tsoc, ausente,
                               usa_ausente, tdil, tfix);
    br::Tabela t = tabela_do_R(dados, nomes);
    br::Pedigree ped;
    const br::Pedigree* pp = nullptr;
    if (XLENGTH(pid) > 0) {
      auto i = textos(pid, "id");
      auto p = textos(ppai, "sire");
      auto ma = textos(pmae, "dam");
      ped = br::constroi_pedigree(i, p, ma, textos(mfx, "metafounders"), std::vector<double>(REAL(gmx), REAL(gmx) + XLENGTH(gmx)), mgs_de(pmae));
      pp = &ped;
    }
    std::vector<br::KernelDecl> kd = kernels_do_R(kern);
    br::Desenho d = br::monta_desenho(m, t, pp, nullptr, kd.empty() ? nullptr : &kd);
    // os kernels ja estao no desenho (a K^-1 de cada grupo vai em d.kinv): a copia de
    // entrada sai aqui, antes das MME e do fator, e nao fica no pico do ajuste
    std::vector<br::KernelDecl>().swap(kd);
    // a cadeia amostra os componentes, salvo com theta_fixed =
    if (XLENGTH(theta_fixo) == 0) br::confere_covariancias_iid(d.modelo, d.aleatorios, d.usa);
    std::string nota = genomica_no_desenho(d, pp, ped, gid, gm, mistura, anucleo, vk);

    br::PrioriGibbs priori;
    priori.tipo = Rf_asInteger(ptipo);
    priori.gl0 = Rf_asReal(pdf);
    priori.escala = Rf_asReal(pesc);
    // meia-t: o multiplicador de cada variancia pelo nome (vazio = o comum, em pesc)
    priori.nomes_mult = textos(pnomes, "half_t scale names");
    priori.mult.assign(REAL(pmult), REAL(pmult) + XLENGTH(pmult));
    if (priori.mult.size() != priori.nomes_mult.size())
      Rf_error("half_t(scale =): one multiplier per name");
    // a priori do residuo a parte (tipo < 0: a mesma dos grupos, como sempre)
    br::PrioriGibbs priori_res;
    priori_res.tipo = Rf_asInteger(rtipo);
    priori_res.gl0 = Rf_asReal(rdf);
    priori_res.escala = Rf_asReal(resc);

    std::vector<double> thf;
    const std::vector<double>* pthf = nullptr;
    if (XLENGTH(theta_fixo) > 0) {
      thf.assign(REAL(theta_fixo), REAL(theta_fixo) + XLENGTH(theta_fixo));
      if (thf.size() != m.ntheta)
        Rf_error("theta with %d entries; the layout asks for %d", (int) thf.size(),
                 (int) m.ntheta);
      pthf = &thf;
    }
    // a priori normal dos fixos, o passo PX-DA e a partida pedida (start =)
    br::ExtrasGibbs extras;
    extras.dp_fixos = Rf_asReal(dpfix);
    extras.px = Rf_asLogical(pxda) == TRUE;
    extras.asis = Rf_asInteger(asis_r);
    // o padrao do gibbs() no probit: o passo ASIS so nos grupos de covariancia
    extras.asis_so_grupos = Rf_asLogical(asis_grupos) == TRUE;
    extras.partida.assign(REAL(part0), REAL(part0) + XLENGTH(part0));
    extras.nomes_partida = textos(part0_nomes, "start names");
    // a priori dos termos escalares a parte da dos grupos (tipo < 0: a mesma `priori`); e o
    // padrao gaussiano do gibbs(), a meia-t nos escalares e a plana nos grupos
    extras.escalar_a_parte = Rf_asInteger(etipo) >= 0;
    extras.priori_escalar.tipo = Rf_asInteger(etipo);
    extras.priori_escalar.gl0 = Rf_asReal(edf);
    extras.priori_escalar.escala = Rf_asReal(eesc);
    // o padrao gaussiano: o grupo com niveis de menos para a plana recua para a dos escalares
    extras.recuo_poucos_niveis = Rf_asLogical(recuo) == TRUE;
    br::GibbsSaida S = br::gibbs(d, (std::size_t) Rf_asInteger(n_iter),
                                 (std::size_t) Rf_asInteger(burnin),
                                 (std::size_t) Rf_asInteger(thin),
                                 Rf_asLogical(loc_fixa) == TRUE, pthf,
                                 Rf_asLogical(verb) == TRUE, priori,
                                 Rf_asInteger(fam) == 1,
                                 (std::size_t) std::max(0, Rf_asInteger(npdm)),
                                 priori_res.tipo < 0 ? nullptr : &priori_res, &extras);
    if (!nota.empty()) S.mensagem = S.mensagem.empty() ? nota : nota + "; " + S.mensagem;

    std::vector<std::string> nomes_th = m.nomes_theta();
    SEXP amostras = PROTECT(Rf_allocMatrix(REALSXP, (int) S.n_amostras, (int) S.ntheta));
    for (std::size_t i = 0; i < S.n_amostras; i++)
      for (std::size_t k = 0; k < S.ntheta; k++)
        REAL(amostras)[k * S.n_amostras + i] = S.amostras[i * S.ntheta + k];
    SEXP nms_t = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) S.ntheta));
    for (std::size_t k = 0; k < S.ntheta; k++)
      SET_STRING_ELT(nms_t, (R_xlen_t) k, Rf_mkChar(nomes_th[k].c_str()));

    // medias e desvios das localizacoes por grupo, com os niveis como nomes
    br::Montado M = br::monta_mme(d, br::partida(d));
    SEXP ebv = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t) m.grupos.size()));
    SEXP esd = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t) m.grupos.size()));
    SEXP ebv_nomes = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) m.grupos.size()));
    for (std::size_t g = 0; g < m.grupos.size(); g++) {
      const std::size_t larg = d.largura(g);
      SEXP v = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) larg));
      SEXP w = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) larg));
      for (std::size_t c = 0; c < larg; c++) {
        const std::size_t k = M.offset_grupo[g] + c;
        REAL(v)[c] = S.media_loc.empty() ? NA_REAL : S.media_loc[k];
        REAL(w)[c] = S.var_loc.empty() ? NA_REAL : std::sqrt(S.var_loc[k]);
      }
      for (const br::DesenhoTermo& a : d.aleatorios)
        if (a.termo == m.grupos[g].termos[0]) {
          SEXP nv = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) larg));
          for (std::size_t c = 0; c < larg; c++)
            SET_STRING_ELT(nv, (R_xlen_t) c, Rf_mkChar(a.niveis[c % a.n_niveis].c_str()));
          Rf_setAttrib(v, R_NamesSymbol, nv);
          Rf_setAttrib(w, R_NamesSymbol, Rf_duplicate(nv));
          UNPROTECT(1);
          break;
        }
      SET_VECTOR_ELT(ebv, (R_xlen_t) g, v);
      SET_VECTOR_ELT(esd, (R_xlen_t) g, w);
      SET_STRING_ELT(ebv_nomes, (R_xlen_t) g, Rf_mkChar(m.grupos[g].nome.c_str()));
      UNPROTECT(2);
    }
    Rf_setAttrib(ebv, R_NamesSymbol, ebv_nomes);
    Rf_setAttrib(esd, R_NamesSymbol, Rf_duplicate(ebv_nomes));

    // Media e dp a posteriori dos efeitos fixos: as primeiras x.ncol posicoes do vetor de
    // localizacao, com os mesmos nomes termo=nivel que dropped_x usa.
    SEXP bfix = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) d.x.ncol));
    SEXP bsd  = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) d.x.ncol));
    SEXP bnms = PROTECT(Rf_allocVector(STRSXP,  (R_xlen_t) d.x.ncol));
    for (std::size_t k = 0; k < d.x.ncol; k++) {
      REAL(bfix)[k] = S.media_b.size() == d.x.ncol ? S.media_b[k] : NA_REAL;
      REAL(bsd)[k]  = S.var_b.size() == d.x.ncol ? std::sqrt(S.var_b[k]) : NA_REAL;
      SET_STRING_ELT(bnms, (R_xlen_t) k, Rf_mkChar(d.nomes_x[k].c_str()));
    }
    Rf_setAttrib(bfix, R_NamesSymbol, bnms);
    Rf_setAttrib(bsd, R_NamesSymbol, Rf_duplicate(bnms));
    // a mesma media na escala padronizada do motor (as x.ncol primeiras posicoes de E[loc],
    // sem a volta a unidade do dado), interna: residual_diagnostics() confere a cadeia com o
    // W E[loc] dela, e o b na unidade do dado, levado de volta a essa escala, perde eps |b0|
    // no cancelamento b0 + m b_x de uma covariavel muito deslocada
    SEXP bpad = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) d.x.ncol));
    for (std::size_t k = 0; k < d.x.ncol; k++)
      REAL(bpad)[k] = S.media_b.size() == d.x.ncol && S.media_loc.size() >= d.x.ncol
                        ? S.media_loc[k] : NA_REAL;
    Rf_setAttrib(bpad, R_NamesSymbol, Rf_duplicate(bnms));
    SEXP saiu = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) d.saiu_x.size()));
    for (std::size_t k = 0; k < d.saiu_x.size(); k++)
      SET_STRING_ELT(saiu, (R_xlen_t) k, Rf_mkChar(d.saiu_x[k].c_str()));

    SEXP part = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) S.partida.size()));
    for (std::size_t k = 0; k < S.partida.size(); k++) REAL(part)[k] = S.partida[k];
    if (S.partida.size() == S.ntheta) Rf_setAttrib(part, R_NamesSymbol, Rf_duplicate(nms_t));
    const char* campos[] = {"samples", "names", "ebv", "ebv_sd", "message", "n_used",
                            "b", "b_sd", "dropped_x", "dense_block", "start",
                            "residual", "pdm_residuals", "fitted", "fitted_fixed",
                            "b_standardized", "prior_scale", "asis", "prior_fallback"};
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 19));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 19));
    for (int q = 0; q < 19; q++) SET_STRING_ELT(nms, q, Rf_mkChar(campos[q]));
    // os grupos que recuaram da priori deles para a dos escalares (ExtrasGibbs)
    {
      SEXP rc = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) S.grupos_recuo.size()));
      for (std::size_t q = 0; q < S.grupos_recuo.size(); q++)
        SET_STRING_ELT(rc, (R_xlen_t) q, Rf_mkChar(S.grupos_recuo[q].c_str()));
      // os niveis e a dimensao de cada um, em atributos que o R le para a mensagem e tira
      SEXP rnl = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) S.grupos_recuo_nl.size()));
      SEXP rd = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) S.grupos_recuo_d.size()));
      for (std::size_t q = 0; q < S.grupos_recuo_nl.size(); q++) {
        INTEGER(rnl)[q] = S.grupos_recuo_nl[q];
        INTEGER(rd)[q] = S.grupos_recuo_d[q];
      }
      Rf_setAttrib(rc, Rf_install("nl"), rnl);
      Rf_setAttrib(rc, Rf_install("d"), rd);
      SET_VECTOR_ELT(out, 18, rc);
      UNPROTECT(3);
    }
    // o passo ASIS por grupo: tentativas, aceitas e fora (grupos x 3, com os nomes dos grupos;
    // 0 linhas sem o passo). O grupo preso por kernel(fixed =) fica com zeros
    {
      const std::size_t ng = S.asis_tentativas.size();
      SEXP am = PROTECT(Rf_allocMatrix(REALSXP, (int) ng, 3));
      for (std::size_t g = 0; g < ng; g++) {
        REAL(am)[g] = S.asis_tentativas[g];
        REAL(am)[ng + g] = S.asis_aceites[g];
        REAL(am)[2 * ng + g] = S.asis_fora[g];
      }
      SEXP dn = PROTECT(Rf_allocVector(VECSXP, 2));
      SEXP rn = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) ng));
      for (std::size_t g = 0; g < ng; g++)
        SET_STRING_ELT(rn, (R_xlen_t) g, Rf_mkChar(m.grupos[g].nome.c_str()));
      SEXP cn = PROTECT(Rf_allocVector(STRSXP, 3));
      SET_STRING_ELT(cn, 0, Rf_mkChar("attempts"));
      SET_STRING_ELT(cn, 1, Rf_mkChar("accepted"));
      SET_STRING_ELT(cn, 2, Rf_mkChar("outside"));
      SET_VECTOR_ELT(dn, 0, rn);
      SET_VECTOR_ELT(dn, 1, cn);
      Rf_setAttrib(am, R_DimNamesSymbol, dn);
      SET_VECTOR_ELT(out, 17, am);
      UNPROTECT(4);
    }
    // a escala A da meia-t de cada variancia que a usa, nomeada (vazio nas outras prioris)
    {
      std::vector<std::size_t> com;
      for (std::size_t k = 0; k < S.escala_meia_t.size() && k < nomes_th.size(); k++)
        if (std::isfinite(S.escala_meia_t[k])) com.push_back(k);
      SEXP ps = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) com.size()));
      SEXP pn = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) com.size()));
      for (std::size_t q = 0; q < com.size(); q++) {
        REAL(ps)[q] = S.escala_meia_t[com[q]];
        SET_STRING_ELT(pn, (R_xlen_t) q, Rf_mkChar(nomes_th[com[q]].c_str()));
      }
      Rf_setAttrib(ps, R_NamesSymbol, pn);
      SET_VECTOR_ELT(out, 16, ps);
      UNPROTECT(2);
    }
    // o residuo medio a posteriori e os ajustados W E[loc] e X E[b] por linha da tabela (NA
    // fora das usadas), e os sorteios pivotais e_s / sqrt(s2e_s), n_usadas x n_pdm
    {
      const std::vector<double>* fontes[3] = {&S.media_residuo, &S.ajustado, &S.ajustado_fixo};
      const int lugar[3] = {11, 13, 14};
      for (int q = 0; q < 3; q++) {
        SEXP mr = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) d.nlin));
        for (std::size_t i = 0; i < d.nlin; i++) {
          const double v = i < fontes[q]->size() ? (*fontes[q])[i] : NA_REAL;
          REAL(mr)[i] = std::isfinite(v) ? v : NA_REAL;
        }
        SET_VECTOR_ELT(out, lugar[q], mr);
        UNPROTECT(1);
      }
      const std::size_t nu = d.n_usadas();
      SEXP pd = PROTECT(Rf_allocMatrix(REALSXP, (int) nu, (int) S.n_pdm));
      for (std::size_t k = 0; k < S.n_pdm; k++)
        for (std::size_t i = 0; i < nu; i++)
          REAL(pd)[k * nu + i] = S.residuos_pdm[k * nu + i];
      SET_VECTOR_ELT(out, 12, pd);
      UNPROTECT(1);
    }
    SET_VECTOR_ELT(out, 9, sexp_bloco_denso(S.bloco_denso, S.colunas_fator));
    SET_VECTOR_ELT(out, 0, amostras);
    SET_VECTOR_ELT(out, 1, nms_t);
    SET_VECTOR_ELT(out, 2, ebv);
    SET_VECTOR_ELT(out, 3, esd);
    SET_VECTOR_ELT(out, 4, Rf_mkString(S.mensagem.c_str()));
    SET_VECTOR_ELT(out, 5, Rf_ScalarInteger((int) d.n_usadas()));
    SET_VECTOR_ELT(out, 6, bfix);
    SET_VECTOR_ELT(out, 7, bsd);
    SET_VECTOR_ELT(out, 8, saiu);
    SET_VECTOR_ELT(out, 10, part);
    SET_VECTOR_ELT(out, 15, bpad);
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(13);
    return out;
  )
}

// O portao V1 da meia-t (interno): as condicionais da priori sem dado, n x d(d+1)/2 amostras
// do vech de C. variante 1 e a forma errada (nu + nl)/2 da condicional de a.
SEXP R_gibbs_priori_ht(SEXP dim, SEXP nu, SEXP A, SEXP n, SEXP thin, SEXP burnin,
                       SEXP variante) {
  GUARDA(
    const std::size_t d = (std::size_t) std::max(0, Rf_asInteger(dim));
    const std::size_t nn = (std::size_t) std::max(0, Rf_asInteger(n));
    std::vector<double> va(REAL(A), REAL(A) + XLENGTH(A));
    std::vector<double> v = br::gibbs_priori_meia_t(
        d, Rf_asReal(nu), va, nn, (std::size_t) std::max(0, Rf_asInteger(thin)),
        (std::size_t) std::max(0, Rf_asInteger(burnin)), Rf_asInteger(variante));
    const std::size_t np = d * (d + 1) / 2;
    SEXP out = PROTECT(Rf_allocMatrix(REALSXP, (int) nn, (int) np));
    for (std::size_t i = 0; i < nn; i++)
      for (std::size_t k = 0; k < np; k++) REAL(out)[k * nn + i] = v[i * np + k];
    UNPROTECT(1);
    return out;
  )
}

// ---------------------------------------------------------------- curva de crescimento

extern "C++" {

// As opcoes da curva de crescimento chegam numa lista NOMEADA (R/growth.R): os campos sao
// lidos pelo nome, e um campo que falta e um erro que diz qual
static SEXP campo_curva(SEXP l, const char* nome) {
  if (TYPEOF(l) != VECSXP) throw br::Erro("growth settings: expected a named list");
  SEXP nm = Rf_getAttrib(l, R_NamesSymbol);
  for (R_xlen_t k = 0; k < XLENGTH(l); k++)
    if (nm != R_NilValue && std::strcmp(CHAR(STRING_ELT(nm, k)), nome) == 0)
      return VECTOR_ELT(l, k);
  throw br::Erro(std::string("growth settings: '") + nome + "' is missing");
}
static std::vector<double> reais_curva(SEXP l, const char* nome) {
  SEXP v = campo_curva(l, nome);
  if (TYPEOF(v) == REALSXP) return std::vector<double>(REAL(v), REAL(v) + XLENGTH(v));
  if (TYPEOF(v) == INTSXP) return std::vector<double>(INTEGER(v), INTEGER(v) + XLENGTH(v));
  if (v == R_NilValue) return {};
  throw br::Erro(std::string("growth settings: '") + nome + "' must be numeric");
}
static std::vector<int> inteiros_curva(SEXP l, const char* nome) {
  SEXP v = campo_curva(l, nome);
  if (TYPEOF(v) == INTSXP) return std::vector<int>(INTEGER(v), INTEGER(v) + XLENGTH(v));
  if (TYPEOF(v) == LGLSXP) return std::vector<int>(LOGICAL(v), LOGICAL(v) + XLENGTH(v));
  if (v == R_NilValue) return {};
  throw br::Erro(std::string("growth settings: '") + nome + "' must be integer");
}
static double real_curva(SEXP l, const char* nome) {
  const std::vector<double> v = reais_curva(l, nome);
  if (v.size() != 1) throw br::Erro(std::string("growth settings: '") + nome + "' must be one number");
  return v[0];
}
static int inteiro_curva(SEXP l, const char* nome) {
  const std::vector<int> v = inteiros_curva(l, nome);
  if (v.size() != 1) throw br::Erro(std::string("growth settings: '") + nome + "' must be one integer");
  return v[0];
}
static br::CurvaConf curva_do_R(SEXP l) {
  br::CurvaConf c;
  c.nu = real_curva(l, "nu");
  c.log_escala = inteiro_curva(l, "log_scale") != 0;
  c.linear = inteiro_curva(l, "linear") != 0;
  if (!std::isfinite(c.nu) || !(c.nu > -1.0)) throw br::Erro("nu must be finite and above -1");
  return c;
}
// as pesagens: o animal (linha da tabela, base 0), a idade, o y na escala de analise, o log do
// peso e a classe do residuo (base 0)
static br::PesagensCurva pesagens_do_R(SEXP l) {
  return br::agrupa_pesagens(inteiros_curva(l, "animal"), reais_curva(l, "time"),
                             reais_curva(l, "y"), reais_curva(l, "log_weight"),
                             inteiros_curva(l, "class"),
                             (std::size_t) inteiro_curva(l, "n_animals"),
                             (std::size_t) inteiro_curva(l, "n_classes"));
}
static br::Densa matriz3_curva(SEXP l, const char* nome) {
  const std::vector<double> v = reais_curva(l, nome);
  if (v.size() != 9) throw br::Erro(std::string("growth settings: '") + nome + "' must be 3 x 3");
  br::Densa m(3, 3);
  for (std::size_t i = 0; i < 3; i++)
    for (std::size_t j = 0; j < 3; j++) m.at(i, j) = v[j * 3 + i];
  return m;
}
// uma matriz do R (por colunas) a partir de dados por linhas
static SEXP matriz_de_linhas(const double* v, std::size_t nl, std::size_t nc) {
  SEXP m = PROTECT(Rf_allocMatrix(REALSXP, (int) nl, (int) nc));
  for (std::size_t i = 0; i < nl; i++)
    for (std::size_t j = 0; j < nc; j++) REAL(m)[j * nl + i] = v[i * nc + j];
  UNPROTECT(1);
  return m;
}
static SEXP vetor_real(const std::vector<double>& v) {
  SEXP r = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) v.size()));
  for (std::size_t k = 0; k < v.size(); k++) REAL(r)[k] = v[k];
  UNPROTECT(1);
  return r;
}
static br::PrioriCurva priori_curva(SEXP p) {
  br::PrioriCurva q;
  q.tipo = inteiro_curva(p, "type");
  q.gl0 = real_curva(p, "df");
  q.escala = reais_curva(p, "scale");
  if (q.tipo < 1 || q.tipo > 3) throw br::Erro("growth prior: unknown type");
  if (q.tipo == 3 && (!(q.gl0 > 0.0) || q.escala.empty()))
    throw br::Erro("growth prior: a proper prior needs df > 0 and a scale");
  return q;
}
// uma lista nomeada de SEXP ja protegidos pelo chamador
static SEXP lista_nomeada(const std::vector<std::pair<std::string, SEXP>>& itens) {
  SEXP out = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t) itens.size()));
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) itens.size()));
  for (std::size_t k = 0; k < itens.size(); k++) {
    SET_STRING_ELT(nms, (R_xlen_t) k, Rf_mkChar(itens[k].first.c_str()));
    SET_VECTOR_ELT(out, (R_xlen_t) k, itens[k].second);
  }
  Rf_setAttrib(out, R_NamesSymbol, nms);
  UNPROTECT(2);
  return out;
}

}  // extern "C++"

// A curva (escala de analise), o jacobiano dela e o peso em cada idade, num phi dado: o
// portao E1 (interno)
SEXP R_curva_f(SEXP phi, SEXP cfg, SEXP tempos) {
  GUARDA(
    const br::CurvaConf c = curva_do_R(cfg);
    if (XLENGTH(phi) != 3 || TYPEOF(phi) != REALSXP) throw br::Erro("phi must be 3 numbers");
    const std::size_t n = (std::size_t) XLENGTH(tempos);
    std::vector<double> f(n), jac(3 * n), peso(n);
    SEXP ok = PROTECT(Rf_allocVector(LGLSXP, (R_xlen_t) n));
    for (std::size_t i = 0; i < n; i++) {
      LOGICAL(ok)[i] = br::curva_analise(c, REAL(tempos)[i], REAL(phi), f[i], jac.data() + 3 * i);
      peso[i] = br::curva_peso(c, REAL(tempos)[i], REAL(phi));
    }
    SEXP vf = PROTECT(vetor_real(f));
    SEXP vj = PROTECT(matriz_de_linhas(jac.data(), n, 3));
    SEXP vp = PROTECT(vetor_real(peso));
    SEXP out = PROTECT(lista_nomeada({{"f", vf}, {"jacobian", vj}, {"ok", ok}, {"weight", vp}}));
    UNPROTECT(5);
    return out;
  )
}

// O modo condicional de cada animal (as partidas S1 e S2 de gibbs_growth())
SEXP R_curva_modo(SEXP cfg) {
  GUARDA(
    const br::CurvaConf c = curva_do_R(cfg);
    const br::PesagensCurva p = pesagens_do_R(cfg);
    const br::ModosCurva M = br::modos_curva(c, p, reais_curva(cfg, "mu"),
                                             matriz3_curva(cfg, "sigma"), reais_curva(cfg, "s2e"));
    const std::size_t n = p.n_anim;
    SEXP m = PROTECT(matriz_de_linhas(M.m.data(), n, 3));
    SEXP h = PROTECT(matriz_de_linhas(M.h.data(), n, 9));
    SEXP ps = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) n));
    SEXP ok = PROTECT(Rf_allocVector(LGLSXP, (R_xlen_t) n));
    for (std::size_t i = 0; i < n; i++) {
      INTEGER(ps)[i] = M.passos[i];
      LOGICAL(ok)[i] = M.ok[i];
    }
    // a curva no modo em cada pesagem, na ordem de entrada (a partida S4 le o residuo)
    SEXP aj = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) p.t.size()));
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = p.ini[i]; j < p.ini[i + 1]; j++) {
        double f;
        br::curva_analise(c, p.t[j], M.m.data() + 3 * i, f, nullptr);
        REAL(aj)[p.ordem[j]] = f;
      }
    SEXP out = PROTECT(lista_nomeada({{"mode", m}, {"h", h}, {"steps", ps}, {"ok", ok},
                                      {"fitted", aj}}));
    UNPROTECT(6);
    return out;
  )
}

// O passo 1 sozinho para um animal (o portao E2, interno)
SEXP R_curva_passo(SEXP cfg) {
  GUARDA(
    const br::CurvaConf c = curva_do_R(cfg);
    const br::PesagensCurva p = pesagens_do_R(cfg);
    const std::vector<double> mu = reais_curva(cfg, "mu");
    const std::vector<double> phi0 = reais_curva(cfg, "phi0");
    if (mu.size() != 3) throw br::Erro("mu must be 3 numbers");
    if (!phi0.empty() && phi0.size() != 3) throw br::Erro("phi0 must be 3 numbers or empty");
    const br::PassoCurvaSaida S = br::passo_curva_isolado(
        c, p, mu.data(), matriz3_curva(cfg, "sigma"), reais_curva(cfg, "s2e"),
        phi0.empty() ? nullptr : phi0.data(), (std::size_t) inteiro_curva(cfg, "burnin"),
        (std::size_t) inteiro_curva(cfg, "n"), inteiro_curva(cfg, "variant"),
        inteiro_curva(cfg, "moves"));
    SEXP cad = PROTECT(matriz_de_linhas(S.cadeia.data(), S.cadeia.size() / 3, 3));
    SEXP m = PROTECT(vetor_real(std::vector<double>(S.m, S.m + 3)));
    SEXP h = PROTECT(matriz_de_linhas(S.h, 3, 3));
    SEXP ac = PROTECT(vetor_real({S.ac_ind, S.ac_rw}));
    SEXP es = PROTECT(Rf_ScalarReal(S.escala));
    SEXP ps = PROTECT(Rf_ScalarInteger(S.passos));
    SEXP ok = PROTECT(Rf_ScalarLogical(S.ok));
    SEXP out = PROTECT(lista_nomeada({{"chain", cad}, {"mode", m}, {"h", h},
                                      {"acceptance", ac}, {"scale", es}, {"steps", ps},
                                      {"ok", ok}}));
    UNPROTECT(8);
    return out;
  )
}

// A cadeia de gibbs_growth(): o desenho do nivel 2 como no model_mt() (a tabela de um animal
// por linha, com phi nas tres caracteristicas) e as pesagens e as opcoes em `cfg`
SEXP R_gibbs_curva(SEXP dados, SEXP nomes, SEXP alvos, SEXP tnome, SEXP tcol, SEXP tcov,
                   SEXP test, SEXP tgrp, SEXP tnest, SEXP tbase, SEXP tsoc, SEXP pid,
                   SEXP ppai, SEXP pmae, SEXP ausente, SEXP usa_ausente, SEXP gid, SEXP gm,
                   SEXP mistura, SEXP anucleo, SEXP vk, SEXP mfx, SEXP gmx, SEXP kern,
                   SEXP tdil, SEXP cfg) {
  GUARDA(
    SEXP alvo1 = PROTECT(Rf_mkString(CHAR(STRING_ELT(alvos, 0))));
    br::Modelo m = modelo_do_R(alvo1, tnome, tcol, tcov, test, tgrp, tnest, tbase, tsoc,
                               ausente, usa_ausente, tdil);
    UNPROTECT(1);
    br::Tabela t = tabela_do_R(dados, nomes);
    br::Pedigree ped;
    const br::Pedigree* pp = nullptr;
    if (XLENGTH(pid) > 0) {
      auto i = textos(pid, "id");
      auto p2 = textos(ppai, "sire");
      auto ma = textos(pmae, "dam");
      ped = br::constroi_pedigree(i, p2, ma, textos(mfx, "metafounders"), std::vector<double>(REAL(gmx), REAL(gmx) + XLENGTH(gmx)), mgs_de(pmae));
      pp = &ped;
    }
    std::vector<std::string> alv = textos(alvos, "traits");
    std::vector<br::KernelDecl> kd = kernels_do_R(kern);
    br::DesenhoMT d = br::monta_desenho_mt(m, alv, t, pp, &kd);
    std::vector<br::KernelDecl>().swap(kd);
    br::confere_covariancias_iid(d.modelo, d.aleatorios, d.usa);
    std::string nota = genomica_no_desenho(d, pp, ped, gid, gm, mistura, anucleo, vk);

    const br::PesagensCurva p = pesagens_do_R(cfg);
    br::OpcoesCurva o;
    o.curva = curva_do_R(cfg);
    o.n_iter = (std::size_t) inteiro_curva(cfg, "n_iter");
    o.burnin = (std::size_t) inteiro_curva(cfg, "burnin");
    o.thin = (std::size_t) inteiro_curva(cfg, "thin");
    // as prioris dos grupos pelo nome do grupo
    const std::size_t ng = d.modelo.grupos.size();
    {
      SEXP pg = campo_curva(cfg, "prior_groups");
      SEXP pgn = Rf_getAttrib(pg, R_NamesSymbol);
      if (TYPEOF(pg) != VECSXP || pgn == R_NilValue)
        throw br::Erro("growth settings: prior_groups must be a named list");
      for (std::size_t g = 0; g < ng; g++) {
        bool achou = false;
        for (R_xlen_t k = 0; k < XLENGTH(pg); k++)
          if (d.modelo.grupos[g].nome == CHAR(STRING_ELT(pgn, k))) {
            o.priori_grupos.push_back(priori_curva(VECTOR_ELT(pg, k)));
            achou = true;
            break;
          }
        if (!achou) throw br::Erro("growth settings: no prior for group '" +
                                   d.modelo.grupos[g].nome + "'");
      }
    }
    o.priori_pe = priori_curva(campo_curva(cfg, "prior_pe"));
    o.priori_res = priori_curva(campo_curva(cfg, "prior_residual"));
    o.dp_fixos = reais_curva(cfg, "prior_fixed");
    if (o.dp_fixos.size() != 3)
      throw br::Erro("growth settings: 'prior_fixed' must be three numbers, one per parameter");
    for (double v : o.dp_fixos)
      if (!(v >= 0.0) || !std::isfinite(v))
        throw br::Erro("growth settings: 'prior_fixed' must be finite and >= 0 (0 = flat)");
    // a partida dos componentes do nivel 2 pelo nome, na ordem do motor
    const std::vector<std::string> nomes_th = br::nomes_theta_do_mt(d);
    {
      const std::vector<double> th = reais_curva(cfg, "theta0");
      SEXP thn = Rf_getAttrib(campo_curva(cfg, "theta0"), R_NamesSymbol);
      if (th.size() != nomes_th.size() || thn == R_NilValue)
        throw br::Erro("growth settings: theta0 must name the " +
                       std::to_string(nomes_th.size()) + " level-2 components");
      o.theta0.assign(nomes_th.size(), std::nan(""));
      for (std::size_t k = 0; k < nomes_th.size(); k++)
        for (std::size_t q = 0; q < th.size(); q++)
          if (nomes_th[k] == CHAR(STRING_ELT(thn, (R_xlen_t) q))) o.theta0[k] = th[q];
      for (std::size_t k = 0; k < nomes_th.size(); k++)
        if (!std::isfinite(o.theta0[k]))
          throw br::Erro("growth settings: no start for '" + nomes_th[k] + "'");
    }
    o.s2e0 = reais_curva(cfg, "s2e0");
    {
      // phi0 chega como matriz n x 3 do R (por colunas)
      const std::vector<double> ph = reais_curva(cfg, "phi0");
      const std::size_t n = p.n_anim;
      if (ph.size() != 3 * n) throw br::Erro("growth settings: phi0 must be n x 3");
      o.phi0.assign(3 * n, 0.0);
      for (std::size_t i = 0; i < n; i++)
        for (std::size_t q = 0; q < 3; q++) o.phi0[i * 3 + q] = ph[q * n + i];
    }
    o.variante = inteiro_curva(cfg, "variant");
    o.movimentos = inteiro_curva(cfg, "moves");
    o.segura = inteiro_curva(cfg, "hold");
    o.n_pdm = (std::size_t) std::max(0, inteiro_curva(cfg, "residual_draws"));
    o.idades = reais_curva(cfg, "ages");
    o.grade = reais_curva(cfg, "grid");
    o.classe_idade = inteiros_curva(cfg, "age_class");
    o.verboso = inteiro_curva(cfg, "verbose") != 0;
    {
      SEXP gn = campo_curva(cfg, "genetic");
      const std::string nome = TYPEOF(gn) == STRSXP && XLENGTH(gn) == 1 ? CHAR(STRING_ELT(gn, 0)) : "";
      for (std::size_t g = 0; g < ng; g++)
        if (!nome.empty() && d.modelo.grupos[g].nome == nome) o.grupo_genetico = (int) g;
      if (!nome.empty() && o.grupo_genetico < 0)
        throw br::Erro("growth settings: no level-2 group '" + nome + "'");
    }
    br::SaidaCurva S = br::gibbs_curva(d, p, o);
    if (!nota.empty()) S.mensagem = S.mensagem.empty() ? nota : nota + "; " + S.mensagem;

    const std::size_t n = p.n_anim, nw = p.t.size(), T = 3;
    int prot = 0;
    // os componentes: os do motor (o pe ainda chamado res) e os residuos por classe
    SEXP rn = campo_curva(cfg, "residual_names");
    if (TYPEOF(rn) != STRSXP || (std::size_t) XLENGTH(rn) != p.n_classes)
      throw br::Erro("growth settings: one residual name per class");
    SEXP amostras = PROTECT(matriz_de_linhas(S.amostras.data(), S.n_amostras, S.ncomp)); prot++;
    SEXP nms_t = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) S.ncomp)); prot++;
    for (std::size_t k = 0; k < nomes_th.size(); k++)
      SET_STRING_ELT(nms_t, (R_xlen_t) k, Rf_mkChar(nomes_th[k].c_str()));
    for (std::size_t c = 0; c < p.n_classes; c++)
      SET_STRING_ELT(nms_t, (R_xlen_t) (nomes_th.size() + c), STRING_ELT(rn, (R_xlen_t) c));
    // por grupo: medias e desvios das localizacoes, niveis x 3, e os niveis
    SEXP ebv = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t) ng)); prot++;
    SEXP esd = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t) ng)); prot++;
    SEXP niv = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t) ng)); prot++;
    SEXP prim = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t) ng)); prot++;
    SEXP gnm = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) ng)); prot++;
    {
      std::size_t off = d.n_fixa;
      for (std::size_t g = 0; g < ng; g++) {
        const std::size_t larg = d.largura(g);
        const br::DesenhoTermo* a = nullptr;
        for (const br::DesenhoTermo& x : d.aleatorios)
          if (x.termo == d.modelo.grupos[g].termos[0]) { a = &x; break; }
        const std::size_t nl = a ? a->n_niveis : larg / T;
        SEXP v = PROTECT(Rf_allocMatrix(REALSXP, (int) nl, (int) T));
        SEXP w = PROTECT(Rf_allocMatrix(REALSXP, (int) nl, (int) T));
        SEXP f1 = PROTECT(Rf_allocMatrix(REALSXP, (int) nl, (int) T));
        SEXP lv = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) nl));
        for (std::size_t c = 0; c < larg && c < nl * T; c++) {
          REAL(v)[c] = S.media_loc.empty() ? NA_REAL : S.media_loc[off + c];
          REAL(w)[c] = S.var_loc.empty() ? NA_REAL : std::sqrt(S.var_loc[off + c]);
          REAL(f1)[c] = S.sol_primeira.empty() ? NA_REAL : S.sol_primeira[off + c];
        }
        for (std::size_t l = 0; l < nl; l++)
          SET_STRING_ELT(lv, (R_xlen_t) l, Rf_mkChar(a ? a->niveis[l].c_str() : ""));
        SET_VECTOR_ELT(ebv, (R_xlen_t) g, v);
        SET_VECTOR_ELT(esd, (R_xlen_t) g, w);
        SET_VECTOR_ELT(prim, (R_xlen_t) g, f1);
        SET_VECTOR_ELT(niv, (R_xlen_t) g, lv);
        SET_STRING_ELT(gnm, (R_xlen_t) g, Rf_mkChar(d.modelo.grupos[g].nome.c_str()));
        UNPROTECT(4);
        off += larg;
      }
    }
    Rf_setAttrib(ebv, R_NamesSymbol, gnm);
    Rf_setAttrib(esd, R_NamesSymbol, gnm);
    Rf_setAttrib(niv, R_NamesSymbol, gnm);
    Rf_setAttrib(prim, R_NamesSymbol, gnm);
    // os fixos, nomeados "coluna|parametro" na ordem de eq_fixa
    const std::size_t nb = d.n_fixa;
    SEXP bm = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) nb)); prot++;
    SEXP bs = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) nb)); prot++;
    SEXP b1 = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) nb)); prot++;
    SEXP bn = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) nb)); prot++;
    for (std::size_t j = 0; j < d.x.ncol; j++)
      for (std::size_t tau = 0; tau < T; tau++) {
        const int e = d.eq_fixa[j * T + tau];
        if (e < 0) continue;
        REAL(bm)[e] = S.media_b.empty() ? NA_REAL : S.media_b[e];
        REAL(bs)[e] = S.var_b.empty() ? NA_REAL : std::sqrt(S.var_b[e]);
        REAL(b1)[e] = S.b_primeira.empty() ? NA_REAL : S.b_primeira[e];
        SET_STRING_ELT(bn, (R_xlen_t) e, Rf_mkChar((d.nomes_x[j] + "|" + d.alvos[tau]).c_str()));
      }
    Rf_setAttrib(bm, R_NamesSymbol, bn);
    Rf_setAttrib(bs, R_NamesSymbol, bn);
    Rf_setAttrib(b1, R_NamesSymbol, bn);
    SEXP saiu = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) d.saiu_x.size())); prot++;
    for (std::size_t k = 0; k < d.saiu_x.size(); k++)
      SET_STRING_ELT(saiu, (R_xlen_t) k, Rf_mkChar(d.saiu_x[k].c_str()));
    // phi e o desvio dele, n x 3
    std::vector<double> sd_phi(S.var_phi.size());
    for (std::size_t k = 0; k < sd_phi.size(); k++) sd_phi[k] = std::sqrt(S.var_phi[k]);
    SEXP mphi = PROTECT(matriz_de_linhas(S.media_phi.data(), n, T)); prot++;
    SEXP sphi = PROTECT(matriz_de_linhas(sd_phi.data(), n, T)); prot++;
    // a aceitacao por animal
    std::vector<double> ac(4 * n);
    for (std::size_t i = 0; i < n; i++) {
      ac[4 * i] = S.ac_ind[i];
      ac[4 * i + 1] = S.ac_rw[i];
      ac[4 * i + 2] = S.escala_rw[i];
      ac[4 * i + 3] = (double) (p.ini[i + 1] - p.ini[i]);
    }
    SEXP mac = PROTECT(matriz_de_linhas(ac.data(), n, 4)); prot++;
    // por pesagem, na ORDEM DE ENTRADA: WAIC, a curva media e os residuos pivotais
    std::vector<double> lse(nw), mld(nw), m2ld(nw), aj(nw), pdm(S.n_pdm * nw);
    for (std::size_t q = 0; q < nw; q++) {
      const std::size_t k = p.ordem[q];
      lse[k] = S.lse[q];
      mld[k] = S.media_ld[q];
      m2ld[k] = S.m2_ld[q];
      aj[k] = S.media_ajuste[q];
      for (std::size_t s = 0; s < S.n_pdm; s++) pdm[s * nw + k] = S.pdm[s * nw + q];
    }
    SEXP vlse = PROTECT(vetor_real(lse)); prot++;
    SEXP vmld = PROTECT(vetor_real(mld)); prot++;
    SEXP vm2 = PROTECT(vetor_real(m2ld)); prot++;
    SEXP vaj = PROTECT(vetor_real(aj)); prot++;
    SEXP vpdm = PROTECT(Rf_allocMatrix(REALSXP, (int) nw, (int) S.n_pdm)); prot++;
    for (std::size_t k = 0; k < pdm.size(); k++) REAL(vpdm)[k] = pdm[k];
    // os derivados
    const std::size_t nd = o.idades.size() + 2, ngr = o.grade.size(), na_id = o.idades.size();
    const std::size_t nlg = S.media_der.size() / (nd ? nd : 1);
    SEXP dm = PROTECT(matriz_de_linhas(S.media_der.data(), nlg, nd)); prot++;
    SEXP dm2 = PROTECT(matriz_de_linhas(S.m2_der.data(), nlg, nd)); prot++;
    SEXP gcm = PROTECT(matriz_de_linhas(S.media_gc.data(), ngr ? S.media_gc.size() / ngr : 0, ngr)); prot++;
    SEXP gcm2 = PROTECT(matriz_de_linhas(S.m2_gc.data(), ngr ? S.m2_gc.size() / ngr : 0, ngr)); prot++;
    SEXP acm = PROTECT(matriz_de_linhas(S.media_ac.data(), n, ngr)); prot++;
    SEXP popm = PROTECT(matriz_de_linhas(S.pop.data(), ngr ? S.pop.size() / ngr : 0, ngr)); prot++;
    SEXP h2m = PROTECT(matriz_de_linhas(S.h2_idade.data(), na_id ? S.h2_idade.size() / na_id : 0,
                                        na_id)); prot++;
    SEXP ll = PROTECT(vetor_real(S.loglik)); prot++;
    SEXP nk = PROTECT(Rf_ScalarInteger((int) S.n_amostras)); prot++;
    SEXP mf = PROTECT(Rf_ScalarReal((double) S.modos_fora)); prot++;
    SEXP msg = PROTECT(Rf_mkString(S.mensagem.c_str())); prot++;
    SEXP out = PROTECT(lista_nomeada({
        {"samples", amostras}, {"names", nms_t}, {"ebv", ebv}, {"ebv_sd", esd},
        {"levels", niv}, {"first_mean", prim}, {"b", bm}, {"b_sd", bs}, {"first_b", b1},
        {"dropped_x", saiu}, {"phi", mphi}, {"phi_sd", sphi}, {"acceptance", mac},
        {"waic_lse", vlse}, {"waic_mean", vmld}, {"waic_m2", vm2}, {"fitted", vaj},
        {"pdm_residuals", vpdm}, {"derived_mean", dm}, {"derived_m2", dm2},
        {"genetic_mean", gcm}, {"genetic_m2", gcm2}, {"animal_curve", acm},
        {"population", popm}, {"h2_age", h2m}, {"loglik", ll},
        {"n_kept", nk}, {"modes_failed", mf}, {"message", msg}})); prot++;
    UNPROTECT(prot);
    return out;
  )
}

SEXP R_snp_blup(SEXP dados, SEXP nomes, SEXP alvo, SEXP tnome, SEXP tcol, SEXP tcov,
                SEXP test, SEXP tgrp, SEXP tnest, SEXP tbase, SEXP tsoc, SEXP pid,
                SEXP ppai, SEXP pmae, SEXP ausente, SEXP usa_ausente, SEXP gid, SEXP gm,
                SEXP rpg, SEXP theta, SEXP tol, SEXP maxiter, SEXP verb, SEXP mfx, SEXP gmx,
                SEXP tdil) {
  GUARDA(
    br::Modelo m = modelo_do_R(alvo, tnome, tcol, tcov, test, tgrp, tnest, tbase, tsoc, ausente, usa_ausente, tdil);
    br::Tabela t = tabela_do_R(dados, nomes);
    if (XLENGTH(pid) == 0) Rf_error("snp_blup needs a pedigree: the model is single step");
    auto i = textos(pid, "id");
    auto p = textos(ppai, "sire");
    auto ma = textos(pmae, "dam");
    br::Pedigree ped = br::constroi_pedigree(i, p, ma, textos(mfx, "metafounders"),
        std::vector<double>(REAL(gmx), REAL(gmx) + XLENGTH(gmx)), mgs_de(pmae));
    br::Desenho d = br::monta_desenho(m, t, &ped);

    std::vector<double> th(REAL(theta), REAL(theta) + XLENGTH(theta));
    if (th.size() != m.ntheta)
      Rf_error("theta with %d entries; the layout asks for %d", (int) th.size(),
               (int) m.ntheta);
    auto gids = textos(gid, "genotypes");
    if (gids.empty()) Rf_error("snp_blup without genotypes has nothing to solve");
    const br::Genotipos mg = genotipos_do_R(gm, gids.size());
    const int nm2 = (int) mg.m;

    const bool com_mf = std::any_of(ped.eh_mf.begin(), ped.eh_mf.end(),
                                    [](char c) { return c != 0; });
    br::SnpBlup S = br::snp_blup(d, mg, gids, th, Rf_asReal(rpg),
                                 (std::size_t) Rf_asInteger(maxiter), Rf_asReal(tol),
                                 Rf_asLogical(verb) == TRUE, com_mf);

    // fatias como no ajuste exato: b nomeado pelas colunas de X, ebv por grupo
    SEXP bfix = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) S.n_fixo));
    SEXP bnms = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) S.n_fixo));
    for (std::size_t k = 0; k < S.n_fixo; k++) {
      REAL(bfix)[k] = S.solucao[k];
      SET_STRING_ELT(bnms, (R_xlen_t) k, Rf_mkChar(d.nomes_x[k].c_str()));
    }
    if (S.n_fixo == d.x.ncol)
      br::desfaz_escala_b(d.volta_x[0], [](std::size_t j) { return (int) j; },
                          REAL(bfix), true);
    Rf_setAttrib(bfix, R_NamesSymbol, bnms);

    SEXP ebv = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t) m.grupos.size()));
    SEXP ebv_nomes = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) m.grupos.size()));
    for (std::size_t g = 0; g < m.grupos.size(); g++) {
      const std::size_t larg = d.largura(g);
      SEXP v = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) larg));
      for (std::size_t c = 0; c < larg; c++)
        REAL(v)[c] = S.solucao[S.offset_grupo[g] + c];
      for (const br::DesenhoTermo& a : d.aleatorios)
        if (a.termo == m.grupos[g].termos[0]) {
          SEXP nv = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) larg));
          for (std::size_t c = 0; c < larg; c++)
            SET_STRING_ELT(nv, (R_xlen_t) c, Rf_mkChar(a.niveis[c % a.n_niveis].c_str()));
          Rf_setAttrib(v, R_NamesSymbol, nv);
          UNPROTECT(1);
          break;
        }
      SET_VECTOR_ELT(ebv, (R_xlen_t) g, v);
      SET_STRING_ELT(ebv_nomes, (R_xlen_t) g, Rf_mkChar(m.grupos[g].nome.c_str()));
      UNPROTECT(1);
    }
    Rf_setAttrib(ebv, R_NamesSymbol, ebv_nomes);

    // uma fatia: o vetor de sempre; mais de uma: matriz marcador x componente
    SEXP ef = PROTECT(S.fatias.size() == 1
        ? Rf_allocVector(REALSXP, (R_xlen_t) S.efeitos.size())
        : Rf_allocMatrix(REALSXP, nm2, (int) S.fatias.size()));
    for (std::size_t k = 0; k < S.efeitos.size(); k++)
      REAL(ef)[k] = std::isnan(S.efeitos[k]) ? NA_REAL : S.efeitos[k];
    if (S.fatias.size() > 1) {
      SEXP dn = PROTECT(Rf_allocVector(VECSXP, 2));
      SEXP cn = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) S.fatias.size()));
      for (std::size_t k = 0; k < S.fatias.size(); k++)
        SET_STRING_ELT(cn, (R_xlen_t) k, Rf_mkChar(S.fatias[k].c_str()));
      SET_VECTOR_ELT(dn, 1, cn);
      Rf_setAttrib(ef, R_DimNamesSymbol, dn);
      UNPROTECT(2);
    }

    std::string nota = "ssSNPBLUP: " + std::to_string(gids.size()) + " genotyped x " +
        std::to_string(nm2) + " markers, " + std::to_string(S.n_imputados) +
        " missing value(s) imputed by the mean, " + std::to_string(S.n_monomorficos) +
        " monomorphic marker(s) left out";

    const char* campos[] = {"b", "ebv", "g", "converged", "iters", "resnorm", "message",
                            "n_used", "n_columns", "dropped_x"};
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 10));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 10));
    for (int q = 0; q < 10; q++) SET_STRING_ELT(nms, q, Rf_mkChar(campos[q]));
    SEXP saiu = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) d.saiu_x.size()));
    for (std::size_t k = 0; k < d.saiu_x.size(); k++)
      SET_STRING_ELT(saiu, (R_xlen_t) k, Rf_mkChar(d.saiu_x[k].c_str()));
    SET_VECTOR_ELT(out, 0, bfix);
    SET_VECTOR_ELT(out, 1, ebv);
    SET_VECTOR_ELT(out, 2, ef);
    SET_VECTOR_ELT(out, 3, Rf_ScalarLogical(S.convergiu));
    SET_VECTOR_ELT(out, 4, Rf_ScalarInteger((int) S.iters));
    SET_VECTOR_ELT(out, 5, Rf_ScalarReal(S.residuo));
    SET_VECTOR_ELT(out, 6, Rf_mkString(nota.c_str()));
    SET_VECTOR_ELT(out, 7, Rf_ScalarInteger((int) d.n_usadas()));
    SET_VECTOR_ELT(out, 8, Rf_ScalarInteger((int) d.total_colunas()));
    SET_VECTOR_ELT(out, 9, saiu);
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(8);
    return out;
  )
}

// PEGS. y (n x k, NA ausente) e x (n x p, sem NA) como matrizes do R; vb0/ve0 vazios = partida
// padrao; atualiza = FALSE com vb0/ve0 e a resolucao pura do ridge multivariado.
static br::Densa densa_do_R(SEXP m) {
  SEXP dim = Rf_getAttrib(m, R_DimSymbol);
  if (Rf_length(dim) != 2) Rf_error("expected a matrix");
  const std::size_t nl = (std::size_t) INTEGER(dim)[0], nc = (std::size_t) INTEGER(dim)[1];
  br::Densa d(nl, nc);
  const double* v = REAL(m);
  for (std::size_t j = 0; j < nc; j++)
    for (std::size_t i = 0; i < nl; i++) d.at(i, j) = v[j * nl + i];
  return d;
}

static SEXP densa_para_R(const br::Densa& d) {
  SEXP m = PROTECT(Rf_allocMatrix(REALSXP, (int) d.nlin, (int) d.ncol));
  for (std::size_t j = 0; j < d.ncol; j++)
    for (std::size_t i = 0; i < d.nlin; i++) REAL(m)[j * d.nlin + i] = d.at(i, j);
  UNPROTECT(1);
  return m;
}

SEXP R_lanczos_g(SEXP gm, SEXP sondas, SEXP passos) {
  GUARDA(
    const br::Genotipos m = genotipos_do_R(gm, (std::size_t) Rf_nrows(gm));
    br::Densa v = densa_do_R(sondas);
    br::RelatorioG rel;
    br::LanczosG r = br::lanczos_g(m, v, (std::size_t) Rf_asInteger(passos), rel);
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 5));
    const char* campos[] = {"alpha", "beta", "steps", "trace", "dim"};
    for (int k = 0; k < 5; k++) SET_STRING_ELT(nms, k, Rf_mkChar(campos[k]));
    SEXP st = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) r.passos.size()));
    for (std::size_t k = 0; k < r.passos.size(); k++) INTEGER(st)[k] = (int) r.passos[k];
    SET_VECTOR_ELT(out, 0, densa_para_R(r.alfa));
    SET_VECTOR_ELT(out, 1, densa_para_R(r.beta));
    SET_VECTOR_ELT(out, 2, st);
    SET_VECTOR_ELT(out, 3, Rf_ScalarReal(r.traco));
    SET_VECTOR_ELT(out, 4, Rf_ScalarReal((double) r.dim));
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(3);
    return out;
  )
}

// Quantos valores fora de 0, 1, 2 e ausente: a checagem de valida_genotipos() sem os vetores
// logicos do tamanho da matriz que o %in% do R criaria (12 bilhoes de entradas num arquivo de
// rebanho grande).
SEXP R_confere_genotipos(SEXP gm) {
  GUARDA(
    double fora = 0.0;
    const R_xlen_t nx = XLENGTH(gm);
    switch (TYPEOF(gm)) {
      case REALSXP: {
        const double* x = REAL(gm);
        for (R_xlen_t k = 0; k < nx; k++)
          if (!std::isnan(x[k]) && x[k] != 0.0 && x[k] != 1.0 && x[k] != 2.0) fora++;
        break;
      }
      case INTSXP: {
        const int* x = INTEGER(gm);
        for (R_xlen_t k = 0; k < nx; k++)
          if (x[k] != NA_INTEGER && (x[k] < 0 || x[k] > 2)) fora++;
        break;
      }
      case RAWSXP: {
        const Rbyte* x = RAW(gm);
        for (R_xlen_t k = 0; k < nx; k++)
          if (x[k] > 2 && x[k] != 5) fora++;
        break;
      }
      default: Rf_error("the genotype matrix must be double, integer or raw");
    }
    return Rf_ScalarReal(fora);
  )
}

// Um hash de 64 bits do genotipo de cada animal (linha de gm), numa passada pela matriz
// (coluna por coluna, a ordem da memoria), para achar animais com o MESMO genotipo em todo
// marcador sem comparar linhas inteiras no R. Devolve os 52 bits de cima como double: os
// candidatos com o mesmo valor sao conferidos de verdade no R.
SEXP R_hash_linhas_genotipos(SEXP gm) {
  GUARDA(
    SEXP dim = Rf_getAttrib(gm, R_DimSymbol);
    if (TYPEOF(dim) != INTSXP || XLENGTH(dim) != 2) Rf_error("the genotypes must be a matrix");
    const std::size_t n = (std::size_t) INTEGER(dim)[0], m = (std::size_t) INTEGER(dim)[1];
    std::vector<std::uint64_t> h(n, 1469598103934665603ULL);
    const std::uint64_t primo = 1099511628211ULL;
    for (std::size_t j = 0; j < m; j++)
      for (std::size_t i = 0; i < n; i++) {
        const std::size_t k = j * n + i;
        unsigned v = 3;   // ausente
        switch (TYPEOF(gm)) {
          case REALSXP: { const double x = REAL(gm)[k]; if (!std::isnan(x)) v = (unsigned) x; break; }
          case INTSXP: { const int x = INTEGER(gm)[k]; if (x != NA_INTEGER) v = (unsigned) x; break; }
          case RAWSXP: { const Rbyte x = RAW(gm)[k]; if (x != 5) v = (unsigned) x; break; }
          default: Rf_error("the genotype matrix must be double, integer or raw");
        }
        h[i] = (h[i] ^ (std::uint64_t) (v + 1)) * primo;
      }
    SEXP out = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) n));
    for (std::size_t i = 0; i < n; i++) REAL(out)[i] = (double) (h[i] >> 12);
    UNPROTECT(1);
    return out;
  )
}

// A media de cada marcador sobre os valores presentes (NaN sem nenhum), lida da vista.
SEXP R_freq_genotipos(SEXP gm) {
  GUARDA(
    const br::Genotipos g = genotipos_do_R(gm, (std::size_t) Rf_nrows(gm));
    SEXP out = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) g.m));
    medias_genotipos(g, REAL(out));
    UNPROTECT(1);
    return out;
  )
}

// Uma matriz n x n do R por allocVector + dim: allocMatrix recusa mais de 2^31 - 1 entradas
// (n > 46340), e um vetor longo com dim e a mesma matriz. Volta desprotegida.
static SEXP matriz_n_n(std::size_t n) {
  SEXP x = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) n * (R_xlen_t) n));
  SEXP d = PROTECT(Rf_allocVector(INTSXP, 2));
  INTEGER(d)[0] = (int) n;
  INTEGER(d)[1] = (int) n;
  Rf_setAttrib(x, R_DimSymbol, d);
  UNPROTECT(2);
  return x;
}

// A G de VanRaden pelo mesmo nucleo do passo unico (vanraden_g), direto da matriz do R e
// escrita na matriz de saida, sem a copia densa intermediaria.
SEXP R_g_matrix(SEXP gm) {
  GUARDA(
    const br::Genotipos g = genotipos_do_R(gm, (std::size_t) Rf_nrows(gm));
    br::RelatorioG rel;
    SEXP out = PROTECT(matriz_n_n(g.n));
    br::vanraden_g_em(g, rel, false, REAL(out));
    UNPROTECT(1);
    return out;
  )
}

// A D de Vitezica, escrita direto na matriz de saida
SEXP R_g_dominancia(SEXP gm) {
  GUARDA(
    const br::Genotipos g = genotipos_do_R(gm, (std::size_t) Rf_nrows(gm));
    SEXP out = PROTECT(matriz_n_n(g.n));
    br::g_dominancia(g, REAL(out));
    UNPROTECT(1);
    return out;
  )
}

// Heterozigotos e observados por animal, list(het, obs)
SEXP R_homozigose(SEXP gm) {
  GUARDA(
    const br::Genotipos g = genotipos_do_R(gm, (std::size_t) Rf_nrows(gm));
    SEXP het = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) g.n));
    SEXP obs = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) g.n));
    br::homozigose(g, REAL(het), REAL(obs));
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(nms, 0, Rf_mkChar("het"));
    SET_STRING_ELT(nms, 1, Rf_mkChar("obs"));
    SET_VECTOR_ELT(out, 0, het);
    SET_VECTOR_ELT(out, 1, obs);
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(4);
    return out;
  )
}

// Uma mascara logica do R (NULL e "todos") num vetor de char de tamanho n; NA e recusado
static std::vector<char> mascara_do_R(SEXP v, std::size_t n, const char* quem) {
  std::vector<char> out;
  if (v == R_NilValue) return out;
  if (TYPEOF(v) != LGLSXP || (std::size_t) XLENGTH(v) != n)
    throw br::Erro(std::string(quem) + ": expected a logical vector of length " +
                   std::to_string(n));
  out.resize(n);
  const int* x = LOGICAL(v);
  for (std::size_t k = 0; k < n; k++) {
    if (x[k] == NA_LOGICAL) throw br::Erro(std::string(quem) + ": NA in the logical vector");
    out[k] = x[k] ? 1 : 0;
  }
  return out;
}

// Contagens de 0, 1 e 2 por marcador nas linhas marcadas (NULL: todas), list(n0, n1, n2)
SEXP R_contagens_marcadores(SEXP gm, SEXP linhas) {
  GUARDA(
    const br::Genotipos g = genotipos_do_R(gm, (std::size_t) Rf_nrows(gm));
    const std::vector<char> lin = mascara_do_R(linhas, g.n, "rows");
    SEXP n0 = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) g.m));
    SEXP n1 = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) g.m));
    SEXP n2 = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) g.m));
    br::contagens_marcadores(g, lin.empty() ? nullptr : lin.data(), REAL(n0), REAL(n1),
                             REAL(n2));
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 3));
    SET_STRING_ELT(nms, 0, Rf_mkChar("n0"));
    SET_STRING_ELT(nms, 1, Rf_mkChar("n1"));
    SET_STRING_ELT(nms, 2, Rf_mkChar("n2"));
    SET_VECTOR_ELT(out, 0, n0);
    SET_VECTOR_ELT(out, 1, n1);
    SET_VECTOR_ELT(out, 2, n2);
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(5);
    return out;
  )
}

// Contagens de 0, 1 e 2 por marcador so nas linhas dadas (inteiras, base 1), list(n0, n1, n2)
SEXP R_contagens_linhas(SEXP gm, SEXP linhas) {
  GUARDA(
    const br::Genotipos g = genotipos_do_R(gm, (std::size_t) Rf_nrows(gm));
    if (TYPEOF(linhas) != INTSXP) throw br::Erro("marker counts: expected integer rows");
    const R_xlen_t k = XLENGTH(linhas);
    std::vector<std::size_t> lin((std::size_t) k);
    for (R_xlen_t a = 0; a < k; a++) {
      const int x = INTEGER(linhas)[a];
      if (x == NA_INTEGER || x < 1) throw br::Erro("marker counts: rows must be positive integers");
      lin[(std::size_t) a] = (std::size_t) x - 1;
    }
    SEXP n0 = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) g.m));
    SEXP n1 = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) g.m));
    SEXP n2 = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) g.m));
    br::contagens_nas_linhas(g, lin, REAL(n0), REAL(n1), REAL(n2));
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 3));
    SET_STRING_ELT(nms, 0, Rf_mkChar("n0"));
    SET_STRING_ELT(nms, 1, Rf_mkChar("n1"));
    SET_STRING_ELT(nms, 2, Rf_mkChar("n2"));
    SET_VECTOR_ELT(out, 0, n0);
    SET_VECTOR_ELT(out, 1, n1);
    SET_VECTOR_ELT(out, 2, n2);
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(5);
    return out;
  )
}

// Chamadas por animal nos marcadores marcados (NULL: todos)
SEXP R_chamadas_animal(SEXP gm, SEXP colunas) {
  GUARDA(
    const br::Genotipos g = genotipos_do_R(gm, (std::size_t) Rf_nrows(gm));
    const std::vector<char> usa = mascara_do_R(colunas, g.m, "markers");
    SEXP cham = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) g.n));
    br::chamadas_por_animal(g, usa.empty() ? nullptr : usa.data(), REAL(cham));
    UNPROTECT(1);
    return cham;
  )
}

// Pares pai-filho: linhas em base 1 do filho e do pai nos genotipos, nos marcadores marcados
// (NULL: todos), list(both, opposite)
SEXP R_conflitos_parentesco(SEXP gm, SEXP filho, SEXP pai, SEXP colunas) {
  GUARDA(
    const br::Genotipos g = genotipos_do_R(gm, (std::size_t) Rf_nrows(gm));
    if (TYPEOF(filho) != INTSXP || TYPEOF(pai) != INTSXP || XLENGTH(filho) != XLENGTH(pai))
      throw br::Erro("parentage: expected two integer vectors of rows of the same length");
    const std::vector<char> usa = mascara_do_R(colunas, g.m, "markers");
    const R_xlen_t np = XLENGTH(filho);
    std::vector<std::size_t> f((std::size_t) np), p((std::size_t) np);
    for (R_xlen_t k = 0; k < np; k++) {
      const int a = INTEGER(filho)[k], b = INTEGER(pai)[k];
      if (a == NA_INTEGER || b == NA_INTEGER || a < 1 || b < 1)
        throw br::Erro("parentage: rows must be positive integers");
      f[(std::size_t) k] = (std::size_t) a - 1;
      p[(std::size_t) k] = (std::size_t) b - 1;
    }
    SEXP ambos = PROTECT(Rf_allocVector(REALSXP, np));
    SEXP opostos = PROTECT(Rf_allocVector(REALSXP, np));
    br::conflitos_parentesco(g, f, p, usa.empty() ? nullptr : usa.data(), REAL(ambos),
                             REAL(opostos));
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(nms, 0, Rf_mkChar("both"));
    SET_STRING_ELT(nms, 1, Rf_mkChar("opposite"));
    SET_VECTOR_ELT(out, 0, ambos);
    SET_VECTOR_ELT(out, 1, opostos);
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(4);
    return out;
  )
}

// (M - c)' V, m x k, sem formar M - c: V e uma matriz double n x k e c um escalar (o ausente
// vale c, isto e, nao contribui)
SEXP R_genotipos_tv(SEXP gm, SEXP v, SEXP centro) {
  GUARDA(
    const br::Genotipos g = genotipos_do_R(gm, (std::size_t) Rf_nrows(gm));
    if (TYPEOF(v) != REALSXP || (std::size_t) Rf_nrows(v) != g.n)
      throw br::Erro("V must be a double matrix with one row per genotyped animal");
    const std::size_t k = (std::size_t) Rf_ncols(v);
    const std::vector<double> c(g.m, Rf_asReal(centro));
    SEXP out = PROTECT(Rf_allocMatrix(REALSXP, (int) g.m, (int) k));
    br::genotipos_tv(g, c.data(), c.data(), REAL(v), k, REAL(out));
    UNPROTECT(1);
    return out;
  )
}

// Quantos genotipos ausentes (NA, NaN ou 5 no raw)
SEXP R_conta_ausentes(SEXP gm) {
  GUARDA(
    const br::Genotipos g = genotipos_do_R(gm, (std::size_t) Rf_nrows(gm));
    return Rf_ScalarReal(br::conta_ausentes(g));
  )
}

// Fracoes de genes dos metafundadores nos genotipados: pai e mae em base 1 (NA desconhecido),
// do pedigree topologico de pedigree(); mf as linhas dos metafundadores na ordem das colunas
// de Gamma; o as linhas dos genotipados. Devolve length(o) x q.
SEXP R_fracoes_genes(SEXP pai, SEXP mae, SEXP mf, SEXP o) {
  GUARDA(
    if (TYPEOF(pai) != INTSXP || TYPEOF(mae) != INTSXP || TYPEOF(mf) != INTSXP ||
        TYPEOF(o) != INTSXP)
      throw br::Erro("gene fractions: expected integer row numbers");
    const R_xlen_t n = XLENGTH(pai);
    if (XLENGTH(mae) != n) throw br::Erro("sire and dam of different lengths");
    std::vector<std::int64_t> s(n), d(n), col(n, -1);
    for (R_xlen_t i = 0; i < n; i++) {
      s[i] = INTEGER(pai)[i] == NA_INTEGER ? -1 : (std::int64_t) INTEGER(pai)[i] - 1;
      d[i] = INTEGER(mae)[i] == NA_INTEGER ? -1 : (std::int64_t) INTEGER(mae)[i] - 1;
    }
    const std::size_t q = (std::size_t) XLENGTH(mf);
    for (std::size_t k = 0; k < q; k++) {
      const int r = INTEGER(mf)[k];
      if (r == NA_INTEGER || r < 1 || r > n)
        throw br::Erro("metafounder without a line in the pedigree");
      col[r - 1] = (std::int64_t) k;
    }
    const R_xlen_t no = XLENGTH(o);
    for (R_xlen_t t = 0; t < no; t++) {
      const int r = INTEGER(o)[t];
      if (r == NA_INTEGER || r < 1 || r > n)
        throw br::Erro("genotyped animal without a line in the pedigree");
    }
    const std::vector<double> Q = br::fracoes_genes(s, d, col, q);
    SEXP out = PROTECT(Rf_allocMatrix(REALSXP, (int) no, (int) q));
    for (R_xlen_t t = 0; t < no; t++) {
      const int r = INTEGER(o)[t];
      for (std::size_t k = 0; k < q; k++)
        REAL(out)[t + (R_xlen_t) k * no] = Q[(std::size_t) (r - 1) * q + k];
    }
    UNPROTECT(1);
    return out;
  )
}

// A22^-1 X sem formar A22^-1 (o V de estimate_gamma()): o pedigree como em a_inverse()
// (metafundadores e Gamma opcionais), gid os genotipados, X n_geno x q double
SEXP R_a22inv_vezes(SEXP id, SEXP pai, SEXP mae, SEXP mfx, SEXP gmx, SEXP gid, SEXP x) {
  GUARDA(
    br::Pedigree ped = br::constroi_pedigree(textos(id, "id"), textos(pai, "sire"),
        textos(mae, "dam"), textos(mfx, "metafounders"),
        std::vector<double>(REAL(gmx), REAL(gmx) + XLENGTH(gmx)), mgs_de(mae));
    const br::Csc ainv = br::a_inversa(ped, br::endogamia(ped));
    const std::vector<std::size_t> idx = br::linhas_genotipados(ped, textos(gid, "genotypes"));
    if (TYPEOF(x) != REALSXP || (std::size_t) Rf_nrows(x) != idx.size())
      throw br::Erro("X must be a double matrix with one row per genotyped animal");
    const std::size_t q = (std::size_t) Rf_ncols(x);
    const std::vector<double> r = br::a22_inversa_vezes(
        ainv, idx, std::vector<double>(REAL(x), REAL(x) + XLENGTH(x)), q, simbolica_do_padrao);
    SEXP out = PROTECT(Rf_allocMatrix(REALSXP, (int) idx.size(), (int) q));
    std::copy(r.begin(), r.end(), REAL(out));
    UNPROTECT(1);
    return out;
  )
}

// Os tres resumos do ML de um metafundador, no pedigree classico (sem metafundadores)
SEXP R_gamma_ml_resumos(SEXP id, SEXP pai, SEXP mae, SEXP gid, SEXP gm) {
  GUARDA(
    br::Pedigree ped = br::constroi_pedigree(textos(id, "id"), textos(pai, "sire"),
                                             textos(mae, "dam"), {}, {}, mgs_de(mae));
    const br::Csc ainv = br::a_inversa(ped, br::endogamia(ped));
    const auto gids = textos(gid, "genotypes");
    const std::vector<std::size_t> idx = br::linhas_genotipados(ped, gids);
    const br::Genotipos g = genotipos_do_R(gm, gids.size());
    const br::ResumosMl r = br::gamma_ml_resumos(ainv, idx, g);
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 3));
    const char* campos[] = {"ss", "tt", "rr"};
    for (int k = 0; k < 3; k++) SET_STRING_ELT(nms, k, Rf_mkChar(campos[k]));
    SET_VECTOR_ELT(out, 0, Rf_ScalarReal(r.s));
    SET_VECTOR_ELT(out, 1, Rf_ScalarReal(r.t));
    SET_VECTOR_ELT(out, 2, Rf_ScalarReal(r.r));
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(2);
    return out;
  )
}

// Efeitos de SNP com a G* do ajuste: o pedigree como em h_inverse() (metafundadores e Gamma
// do ajuste), os genotipados, u (o valor genetico de cada genotipado, na ordem de gid), a
// mistura, o nucleo APY e o k de Vecchia que o ajuste usou, e hp, a diagonal de G* que o
// ajuste guardou (na ordem de gid; vazia = sem conferencia)
SEXP R_efeitos_snp(SEXP id, SEXP pai, SEXP mae, SEXP mfx, SEXP gmx, SEXP gid, SEXP gm,
                   SEXP u, SEXP mistura, SEXP anucleo, SEXP vk, SEXP hp) {
  GUARDA(
    br::Pedigree ped = br::constroi_pedigree(textos(id, "id"), textos(pai, "sire"),
        textos(mae, "dam"), textos(mfx, "metafounders"),
        std::vector<double>(REAL(gmx), REAL(gmx) + XLENGTH(gmx)), mgs_de(mae));
    const auto gids = textos(gid, "genotypes");
    const br::Genotipos g = genotipos_do_R(gm, gids.size());
    if (TYPEOF(u) != REALSXP || TYPEOF(hp) != REALSXP)
      throw br::Erro("u and the stored diagonal of G* must be double");
    std::vector<std::string> nuc;
    if (XLENGTH(anucleo) > 0) nuc = textos(anucleo, "APY core");
    const br::EfeitosSnp r = br::efeitos_snp(
        ped, gids, g, std::vector<double>(REAL(u), REAL(u) + XLENGTH(u)), Rf_asReal(mistura),
        nuc, (std::size_t) std::max(0, Rf_asInteger(vk)),
        std::vector<double>(REAL(hp), REAL(hp) + XLENGTH(hp)));
    const R_xlen_t nu = (R_xlen_t) r.usados.size();
    SEXP ef = PROTECT(Rf_allocVector(REALSXP, nu));
    SEXP fr = PROTECT(Rf_allocVector(REALSXP, nu));
    SEXP us = PROTECT(Rf_allocVector(INTSXP, nu));
    for (R_xlen_t k = 0; k < nu; k++) {
      REAL(ef)[k] = r.efeito[(std::size_t) k];
      REAL(fr)[k] = r.freq[(std::size_t) k];
      INTEGER(us)[k] = (int) r.usados[(std::size_t) k] + 1;
    }
    SEXP dg = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) r.diag_gstar.size()));
    std::copy(r.diag_gstar.begin(), r.diag_gstar.end(), REAL(dg));
    const char* campos[] = {"effect", "freq", "used", "denom", "a", "b", "diag_gstar", "route"};
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 8));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 8));
    for (int k = 0; k < 8; k++) SET_STRING_ELT(nms, k, Rf_mkChar(campos[k]));
    SET_VECTOR_ELT(out, 0, ef);
    SET_VECTOR_ELT(out, 1, fr);
    SET_VECTOR_ELT(out, 2, us);
    SET_VECTOR_ELT(out, 3, Rf_ScalarReal(r.denom));
    SET_VECTOR_ELT(out, 4, Rf_ScalarReal(r.a));
    SET_VECTOR_ELT(out, 5, Rf_ScalarReal(r.b));
    SET_VECTOR_ELT(out, 6, dg);
    SET_VECTOR_ELT(out, 7, Rf_mkString(r.rota == 2 ? "apy" : r.rota == 3 ? "vecchia" : "exact"));
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(6);
    return out;
  )
}

// A Gram do lado menor de Z para o eigen exato de apy_core_select(), e a escala 2 sum p(1-p).
SEXP R_gram_genotipos(SEXP gm) {
  GUARDA(
    const br::Genotipos g = genotipos_do_R(gm, (std::size_t) Rf_nrows(gm));
    br::RelatorioG rel;
    double escala = 0.0;
    br::Densa c = br::gram_menor(g, rel, escala);
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(nms, 0, Rf_mkChar("gram"));
    SET_STRING_ELT(nms, 1, Rf_mkChar("scale"));
    SET_VECTOR_ELT(out, 0, densa_para_R(c));
    SET_VECTOR_ELT(out, 1, Rf_ScalarReal(escala));
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(2);
    return out;
  )
}

// O arquivo de genotipos do BLUPF90 (SNP_FILE: id, espacos, uma cadeia de 0/1/2 com 5 para o
// ausente, largura fixa), numa matriz RAW, 1 byte por genotipo. Duas passadas: a primeira conta
// as linhas que ficam e confere a largura, a segunda preenche, 256 linhas por vez transpostas
// para a ordem por coluna do R (escrever linha a linha espalharia cada byte numa linha de cache).
SEXP R_le_snp_blupf90(SEXP arq, SEXP quais) {
  GUARDA(
    const std::string f = CHAR(STRING_ELT(arq, 0));
    std::ifstream in(f, std::ios::binary);
    if (!in) Rf_error("cannot open '%s'", f.c_str());
    std::unordered_set<std::string> manter;
    for (R_xlen_t k = 0; k < XLENGTH(quais); k++) manter.insert(CHAR(STRING_ELT(quais, k)));
    auto corta = [](const std::string& l, std::string& id, std::size_t& g0, std::size_t& g1) {
      std::size_t a = 0, e = l.size();
      while (e > 0 && (l[e - 1] == '\r' || l[e - 1] == ' ' || l[e - 1] == '\t')) e--;
      while (a < e && (l[a] == ' ' || l[a] == '\t')) a++;
      std::size_t b = a;
      while (b < e && l[b] != ' ' && l[b] != '\t') b++;
      id = l.substr(a, b - a);
      while (b < e && (l[b] == ' ' || l[b] == '\t')) b++;
      g0 = b; g1 = e;
    };
    std::vector<std::string> ids;
    std::vector<char> fica;
    std::string linha, id;
    std::size_t m = 0, g0 = 0, g1 = 0, nl0 = 0;
    while (std::getline(in, linha)) {
      nl0++;
      corta(linha, id, g0, g1);
      if (id.empty()) { fica.push_back(0); continue; }
      if (g1 <= g0) Rf_error("line %d has an id and no genotypes", (int) nl0);
      if (m == 0) m = g1 - g0;
      else if (g1 - g0 != m)
        Rf_error("line %d has %d genotypes and the first had %d", (int) nl0, (int) (g1 - g0), (int) m);
      const bool ok = manter.empty() || manter.count(id) > 0;
      fica.push_back(ok ? 1 : 0);
      if (ok) ids.push_back(id);
    }
    const std::size_t n = ids.size();
    in.clear();
    in.seekg(0);
    // vetor longo com dim: allocMatrix recusa mais de 2^31 entradas, e um arquivo de rebanho
    // passa disso
    SEXP mat = PROTECT(Rf_allocVector(RAWSXP, (R_xlen_t) n * (R_xlen_t) m));
    SEXP dm = PROTECT(Rf_allocVector(INTSXP, 2));
    INTEGER(dm)[0] = (int) n;
    INTEGER(dm)[1] = (int) m;
    Rf_setAttrib(mat, R_DimSymbol, dm);
    Rbyte* x = RAW(mat);
    const std::size_t L = 256;
    std::vector<unsigned char> buf(L * m);
    std::size_t lido = 0, r0 = 0, nb = 0;
    auto despeja = [&]() {
      for (std::size_t j = 0; j < m; j++) {
        Rbyte* col = x + j * n + r0;
        for (std::size_t q = 0; q < nb; q++) col[q] = buf[q * m + j];
      }
      r0 += nb;
      nb = 0;
    };
    while (std::getline(in, linha)) {
      if (!fica[lido++]) continue;
      corta(linha, id, g0, g1);
      unsigned char* d = &buf[nb * m];
      for (std::size_t j = 0; j < m; j++) {
        const char c = linha[g0 + j];
        if (c == '0' || c == '1' || c == '2' || c == '5') d[j] = (unsigned char) (c - '0');
        else Rf_error("unexpected genotype code '%c' at line %d, marker %d (expected 0, 1, 2 or 5)",
                      c, (int) lido, (int) j + 1);
      }
      if (++nb == L) despeja();
    }
    if (nb > 0) despeja();
    SEXP vid = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) n));
    for (std::size_t k = 0; k < n; k++) SET_STRING_ELT(vid, (R_xlen_t) k, Rf_mkChar(ids[k].c_str()));
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(nms, 0, Rf_mkChar("ids"));
    SET_STRING_ELT(nms, 1, Rf_mkChar("m"));
    SET_VECTOR_ELT(out, 0, vid);
    SET_VECTOR_ELT(out, 1, mat);
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(5);
    return out;
  )
}

SEXP R_pegs(SEXP y, SEXP x, SEXP maxit, SEXP tol, SEXP defl, SEXP atualiza, SEXP vb0,
            SEXP ve0, SEXP estrutura, SEXP nfat) {
  GUARDA(
    br::Densa yd = densa_do_R(y);
    const br::Genotipos xd = genotipos_do_R(x, (std::size_t) Rf_nrows(x));
    br::Densa vb;
    std::vector<double> ve;
    const bool tem_vb = XLENGTH(vb0) > 0, tem_ve = XLENGTH(ve0) > 0;
    if (tem_vb) vb = densa_do_R(vb0);
    if (tem_ve) ve.assign(REAL(ve0), REAL(ve0) + XLENGTH(ve0));
    br::ResultadoPegs r = br::pegs(yd, xd, (std::size_t) Rf_asInteger(maxit), Rf_asReal(tol),
                                   Rf_asReal(defl), Rf_asLogical(atualiza) == TRUE,
                                   tem_vb ? &vb : nullptr, tem_ve ? &ve : nullptr,
                                   Rf_asInteger(estrutura), (std::size_t) Rf_asInteger(nfat));
    const char* campos[] = {"mu", "marker_effects", "gebv", "h2", "Vb", "Ve", "bend",
                            "iters", "converged", "n_imputed"};
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 10));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 10));
    for (int q = 0; q < 10; q++) SET_STRING_ELT(nms, q, Rf_mkChar(campos[q]));
    SEXP mu = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) r.mu.size()));
    std::copy(r.mu.begin(), r.mu.end(), REAL(mu));
    SEXP h2 = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) r.h2.size()));
    std::copy(r.h2.begin(), r.h2.end(), REAL(h2));
    SEXP vev = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) r.ve.size()));
    std::copy(r.ve.begin(), r.ve.end(), REAL(vev));
    SET_VECTOR_ELT(out, 0, mu);
    SET_VECTOR_ELT(out, 1, densa_para_R(r.b));
    SET_VECTOR_ELT(out, 2, densa_para_R(r.gebv));
    SET_VECTOR_ELT(out, 3, h2);
    SET_VECTOR_ELT(out, 4, densa_para_R(r.vb));
    SET_VECTOR_ELT(out, 5, vev);
    SET_VECTOR_ELT(out, 6, Rf_ScalarReal(r.deflate));
    SET_VECTOR_ELT(out, 7, Rf_ScalarInteger((int) r.iters));
    SET_VECTOR_ELT(out, 8, Rf_ScalarLogical(r.convergiu));
    SET_VECTOR_ELT(out, 9, Rf_ScalarReal((double) r.n_imputados));
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(5);
    return out;
  )
}

// a estrutura de covariancia sozinha, para os portoes de identidade
SEXP R_pegs_estrutura(SEXP v, SEXP tipo, SEXP nfat) {
  GUARDA(
    br::Densa d = densa_do_R(v);
    br::estrutura_pegs(d, Rf_asInteger(tipo), (std::size_t) Rf_asInteger(nfat));
    return densa_para_R(d);
  )
}

SEXP R_versao(void) { return Rf_mkString("0.5.0"); }

// Os rotulos que o motor da a numeros (br::rotulo_numero, o de Coluna::rotulo), para o R
// escrever um id numerico do pedigree, dos genotipos ou de uma chave EXATAMENTE como a
// coluna numerica dos dados sai daqui. rotulo_motor() em R/pedigree.R e a unica porta. NA
// (e NaN) segue NA, para o R tratar o ausente como sempre tratou.
//
// Aqui havia uma conferencia de COLISAO, dois numeros distintos com o mesmo rotulo (o %g de
// um nao inteiro guardava 6 digitos significativos), que rodava por chamada e portanto por
// coluna: o id 1.1234567 numa coluna e o pai 1.1234568 noutra davam os dois "1.12346" e
// viravam o mesmo animal, e a coluna de dados do motor nunca passava por ela. O rotulo agora
// volta ao mesmo double (rotulo_numero), entao numeros distintos tem rotulos distintos em
// qualquer uniao de colunas, e a conferencia nao tinha mais o que achar.
SEXP R_rotulos(SEXP x) {
  GUARDA(
    if (TYPEOF(x) != REALSXP) throw br::Erro("rotulos: expected a double vector");
    const R_xlen_t n = XLENGTH(x);
    const double* v = REAL(x);
    SEXP out = PROTECT(Rf_allocVector(STRSXP, n));
    for (R_xlen_t i = 0; i < n; i++) {
      if (std::isnan(v[i])) { SET_STRING_ELT(out, i, NA_STRING); continue; }
      SET_STRING_ELT(out, i, Rf_mkChar(br::rotulo_numero(v[i]).c_str()));
    }
    UNPROTECT(1);
    return out;
  )
}

// As COMPONENTES CONEXAS de um grafo de n nos com as arestas (a[k], b[k]), base 1, por
// uniao-busca. A raiz de cada componente e o MENOR no dela (a uniao pendura a raiz maior na
// menor), entao o resultado nao depende da ordem das arestas nem de thread nenhum: e o rotulo
// que blocos_geneticos() (R/estrutura.R) le para achar os blocos de registros sem parentesco
// entre si. Devolve, por no, o indice (base 1) da raiz.
SEXP R_componentes(SEXP n, SEXP a, SEXP b) {
  GUARDA(
    const int nn = Rf_asInteger(n);
    if (nn == NA_INTEGER || nn < 0) throw br::Erro("componentes: invalid number of nodes");
    if (TYPEOF(a) != INTSXP || TYPEOF(b) != INTSXP || XLENGTH(a) != XLENGTH(b))
      throw br::Erro("componentes: expected two integer vectors of the same length");
    std::vector<int> pai(static_cast<std::size_t>(nn));
    for (int i = 0; i < nn; i++) pai[static_cast<std::size_t>(i)] = i;
    auto raiz = [&](int x) {
      while (pai[static_cast<std::size_t>(x)] != x) {
        pai[static_cast<std::size_t>(x)] = pai[static_cast<std::size_t>(pai[static_cast<std::size_t>(x)])];
        x = pai[static_cast<std::size_t>(x)];
      }
      return x;
    };
    const R_xlen_t m = XLENGTH(a);
    const int* pa = INTEGER(a);
    const int* pb = INTEGER(b);
    for (R_xlen_t k = 0; k < m; k++) {
      const int u = pa[k], v = pb[k];
      if (u == NA_INTEGER || v == NA_INTEGER || u < 1 || v < 1 || u > nn || v > nn)
        throw br::Erro("componentes: edge outside the nodes");
      const int ru = raiz(u - 1), rv = raiz(v - 1);
      if (ru == rv) continue;
      if (ru < rv) pai[static_cast<std::size_t>(rv)] = ru;
      else pai[static_cast<std::size_t>(ru)] = rv;
    }
    SEXP out = PROTECT(Rf_allocVector(INTSXP, nn));
    for (int i = 0; i < nn; i++) INTEGER(out)[i] = raiz(i) + 1;
    UNPROTECT(1);
    return out;
  )
}

// Threads das regioes paralelas: n >= 1 define, qualquer outro valor so consulta; lapack
// TRUE/FALSE escolhe a rota das inversas densas grandes, NA so consulta.
SEXP R_threads(SEXP n, SEXP lapack) {
  const int k = Rf_asInteger(n);
  if (k != NA_INTEGER && k >= 1) br::define_threads(k);
  const int lp = Rf_asLogical(lapack);
  if (lp != NA_LOGICAL) br::define_denso_lapack(lp == TRUE);
  SEXP out = PROTECT(Rf_allocVector(INTSXP, 4));
  INTEGER(out)[0] = br::threads();
  INTEGER(out)[1] = br::threads_disponiveis();
#ifdef _OPENMP
  INTEGER(out)[2] = 1;
#else
  INTEGER(out)[2] = 0;
#endif
  INTEGER(out)[3] = br::denso_lapack() ? 1 : 0;
  UNPROTECT(1);
  return out;
}

// A rota de um posto_esparso para os portoes: colunas pela Cholesky, colunas pelo residuo
// explicito, correcoes semi-normais, 1 se o posto saiu da rota densa de reserva (0 se nao) e a
// coluna (base 1) que a rota esparsa nao certificou (0 se nenhuma)
static SEXP rota_do_posto(const br::PostoInfo& info) {
  SEXP ro = PROTECT(Rf_allocVector(INTSXP, 5));
  INTEGER(ro)[0] = (int) info.rapidas;
  INTEGER(ro)[1] = (int) info.explicitas;
  INTEGER(ro)[2] = (int) info.correcoes;
  INTEGER(ro)[3] = info.densa ? 1 : 0;
  INTEGER(ro)[4] = info.densa ? (int) info.coluna_densa + 1 : 0;
  SEXP rn = PROTECT(Rf_allocVector(STRSXP, 5));
  SET_STRING_ELT(rn, 0, Rf_mkChar("fast"));
  SET_STRING_ELT(rn, 1, Rf_mkChar("explicit"));
  SET_STRING_ELT(rn, 2, Rf_mkChar("corrections"));
  SET_STRING_ELT(rn, 3, Rf_mkChar("dense"));
  SET_STRING_ELT(rn, 4, Rf_mkChar("column"));
  Rf_setAttrib(ro, R_NamesSymbol, rn);
  UNPROTECT(2);
  return ro;
}

// O portao do posto esparso dos efeitos fixos: as colunas que ficam de uma X densa do R (n x
// p, por colunas) sobre as linhas ativas, pela rota esparsa dos desenhos (posto_esparso) e,
// com com_denso, pelo Gram-Schmidt denso de antes (posto_completo). Indices base 1. route diz
// por onde a rota esparsa decidiu (rota_do_posto).
SEXP R_posto_fixo(SEXP x, SEXP ativa, SEXP tol, SEXP com_denso) {
  GUARDA(
    if (!Rf_isMatrix(x) || TYPEOF(x) != REALSXP) throw br::Erro("x must be a double matrix");
    const std::size_t n = (std::size_t) Rf_nrows(x), p = (std::size_t) Rf_ncols(x);
    if (TYPEOF(ativa) != LGLSXP || (std::size_t) XLENGTH(ativa) != n)
      throw br::Erro("ativa must be a logical vector with one entry per row");
    const double* xv = REAL(x);
    const double tl = Rf_asReal(tol);
    std::vector<char> at(n);
    for (std::size_t i = 0; i < n; i++) at[i] = LOGICAL(ativa)[i] == TRUE;
    // a CSC com os nao-zeros da matriz (o NaN entra, como no desenho)
    br::Csc cx(n, p);
    for (std::size_t j = 0; j < p; j++) {
      for (std::size_t i = 0; i < n; i++) {
        const double v = xv[j * n + i];
        if (v == 0.0) continue;
        cx.linha.push_back(static_cast<std::uint32_t>(i));
        cx.valor.push_back(v);
      }
      cx.colptr[j + 1] = cx.linha.size();
    }
    std::vector<std::size_t> fica, sai;
    br::PostoInfo info;
    br::posto_esparso(cx, at, tl, fica, sai, false, &info);
    std::vector<std::size_t> fica_d;
    if (Rf_asLogical(com_denso) == TRUE) {
      std::size_t na = 0;
      for (char a : at) na += a;
      br::Densa xd(na, p);
      std::size_t r = 0;
      for (std::size_t i = 0; i < n; i++) {
        if (!at[i]) continue;
        for (std::size_t j = 0; j < p; j++) xd.at(r, j) = xv[j * n + i];
        r++;
      }
      std::vector<std::size_t> sai_d;
      br::posto_completo(xd, tl, fica_d, sai_d);
    }
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 3));
    SEXP fe = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) fica.size()));
    for (std::size_t k = 0; k < fica.size(); k++) INTEGER(fe)[k] = (int) fica[k] + 1;
    SEXP fd = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) fica_d.size()));
    for (std::size_t k = 0; k < fica_d.size(); k++) INTEGER(fd)[k] = (int) fica_d[k] + 1;
    SEXP ro = PROTECT(rota_do_posto(info));
    SET_VECTOR_ELT(out, 0, fe);
    SET_VECTOR_ELT(out, 1, fd);
    SET_VECTOR_ELT(out, 2, ro);
    SET_STRING_ELT(nms, 0, Rf_mkChar("sparse"));
    SET_STRING_ELT(nms, 1, Rf_mkChar("dense"));
    SET_STRING_ELT(nms, 2, Rf_mkChar("route"));
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(5);
    return out;
  )
}

// O posto do X esparso dos ajustadores do lado R (limiar e sobrevivencia, R/threshold.R): um
// bloco por termo, cols[[b]] a coluna de cada linha (base 1; 0 = sem entrada) e vals[[b]] o
// valor. Devolve as colunas que ficam (base 1), pela rota de posto_esparso, com a rota no
// atributo route (rota_do_posto).
SEXP R_posto_blocos(SEXP cols, SEXP vals, SEXP n, SEXP p, SEXP tol, SEXP relativo) {
  GUARDA(
    const std::size_t nl = (std::size_t) Rf_asInteger(n), nc = (std::size_t) Rf_asInteger(p);
    if (TYPEOF(cols) != VECSXP || TYPEOF(vals) != VECSXP || XLENGTH(vals) != XLENGTH(cols))
      throw br::Erro("cols and vals must be lists of the same length");
    const R_xlen_t nb = XLENGTH(cols);
    br::Csc x(nl, nc);
    for (R_xlen_t b = 0; b < nb; b++) {
      SEXP c = VECTOR_ELT(cols, b), v = VECTOR_ELT(vals, b);
      if (TYPEOF(c) != INTSXP || TYPEOF(v) != REALSXP || (std::size_t) XLENGTH(c) != nl ||
          (std::size_t) XLENGTH(v) != nl)
        throw br::Erro("each block needs an integer column and a double value per row");
      for (std::size_t i = 0; i < nl; i++) {
        const int j = INTEGER(c)[i];
        if (j == NA_INTEGER || j < 0 || (std::size_t) j > nc)
          throw br::Erro("block column outside 0..p");
        if (j > 0 && REAL(v)[i] != 0.0) x.colptr[(std::size_t) j]++;
      }
    }
    for (std::size_t j = 0; j < nc; j++) x.colptr[j + 1] += x.colptr[j];
    x.linha.assign(x.colptr[nc], 0);
    x.valor.assign(x.colptr[nc], 0.0);
    // uma coluna pertence a um bloco so, e cada bloco anda nas linhas em ordem: as linhas de
    // cada coluna saem crescentes
    std::vector<std::size_t> prox(x.colptr.begin(), x.colptr.end() - 1);
    for (R_xlen_t b = 0; b < nb; b++) {
      const int* c = INTEGER(VECTOR_ELT(cols, b));
      const double* v = REAL(VECTOR_ELT(vals, b));
      for (std::size_t i = 0; i < nl; i++) {
        if (c[i] <= 0 || v[i] == 0.0) continue;
        const std::size_t j = (std::size_t) c[i] - 1;
        x.linha[prox[j]] = static_cast<std::uint32_t>(i);
        x.valor[prox[j]++] = v[i];
      }
    }
    std::vector<std::size_t> fica, sai;
    br::PostoInfo info;
    br::posto_esparso(x, std::vector<char>(nl, 1), Rf_asReal(tol), fica, sai,
                      Rf_asLogical(relativo) == TRUE, &info);
    SEXP out = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) fica.size()));
    for (std::size_t k = 0; k < fica.size(); k++) INTEGER(out)[k] = (int) fica[k] + 1;
    // a rota, como no R_posto_fixo: so os portoes a leem
    SEXP ro = PROTECT(rota_do_posto(info));
    Rf_setAttrib(out, Rf_install("route"), ro);
    UNPROTECT(2);
    return out;
  )
}

static const R_CallMethodDef metodos[] = {
  {"R_pedigree",   (DL_FUNC) &R_pedigree,   5},
  {"R_a_inversa",  (DL_FUNC) &R_a_inversa,  5},
  {"R_dominancia_inversa", (DL_FUNC) &R_dominancia_inversa, 6},
  {"R_h_inversa",  (DL_FUNC) &R_h_inversa, 10},
  {"R_inv_pd",       (DL_FUNC) &R_inv_pd,       1},
  {"R_chol_esparsa", (DL_FUNC) &R_chol_esparsa, 5},
  {"R_inv_seletiva", (DL_FUNC) &R_inv_seletiva, 5},
  {"R_resolve",      (DL_FUNC) &R_resolve,      5},
  {"R_avaliar",      (DL_FUNC) &R_avaliar,     29},
  {"R_ajustar",      (DL_FUNC) &R_ajustar,     33},
  {"R_a22_inversa",  (DL_FUNC) &R_a22_inversa,  4},
  {"R_avaliar_mt",   (DL_FUNC) &R_avaliar_mt,  23},
  {"R_ajustar_mt",   (DL_FUNC) &R_ajustar_mt,  31},
  {"R_avaliar_ar1",  (DL_FUNC) &R_avaliar_ar1, 24},
  {"R_ajustar_ar1",  (DL_FUNC) &R_ajustar_ar1, 32},
  {"R_residuos_prepara",     (DL_FUNC) &R_residuos_prepara,     32},
  {"R_residuos_prepara_ar1", (DL_FUNC) &R_residuos_prepara_ar1, 30},
  {"R_residuos_lote",        (DL_FUNC) &R_residuos_lote,         4},
  {"R_residuos_libera",      (DL_FUNC) &R_residuos_libera,       1},
    {"R_gibbs", (DL_FUNC) &R_gibbs, 52},
    {"R_gibbs_priori_ht", (DL_FUNC) &R_gibbs_priori_ht, 7},
  {"R_curva_f",     (DL_FUNC) &R_curva_f,      3},
  {"R_curva_modo",  (DL_FUNC) &R_curva_modo,   1},
  {"R_curva_passo", (DL_FUNC) &R_curva_passo,  1},
  {"R_gibbs_curva", (DL_FUNC) &R_gibbs_curva, 26},
  {"R_snp_blup",   (DL_FUNC) &R_snp_blup,  26},
  {"R_lanczos_g",  (DL_FUNC) &R_lanczos_g,  3},
  {"R_confere_genotipos", (DL_FUNC) &R_confere_genotipos, 1},
  {"R_hash_linhas_genotipos", (DL_FUNC) &R_hash_linhas_genotipos, 1},
  {"R_freq_genotipos", (DL_FUNC) &R_freq_genotipos, 1},
  {"R_g_matrix",   (DL_FUNC) &R_g_matrix,   1},
  {"R_g_dominancia", (DL_FUNC) &R_g_dominancia, 1},
  {"R_homozigose", (DL_FUNC) &R_homozigose, 1},
  {"R_contagens_marcadores", (DL_FUNC) &R_contagens_marcadores, 2},
  {"R_contagens_linhas", (DL_FUNC) &R_contagens_linhas, 2},
  {"R_chamadas_animal", (DL_FUNC) &R_chamadas_animal, 2},
  {"R_conflitos_parentesco", (DL_FUNC) &R_conflitos_parentesco, 4},
  {"R_genotipos_tv", (DL_FUNC) &R_genotipos_tv, 3},
  {"R_conta_ausentes", (DL_FUNC) &R_conta_ausentes, 1},
  {"R_fracoes_genes", (DL_FUNC) &R_fracoes_genes, 4},
  {"R_a22inv_vezes", (DL_FUNC) &R_a22inv_vezes, 7},
  {"R_gamma_ml_resumos", (DL_FUNC) &R_gamma_ml_resumos, 5},
  {"R_efeitos_snp", (DL_FUNC) &R_efeitos_snp, 12},
  {"R_gram_genotipos", (DL_FUNC) &R_gram_genotipos, 1},
  {"R_le_snp_blupf90", (DL_FUNC) &R_le_snp_blupf90, 2},
{"R_versao",     (DL_FUNC) &R_versao,     0},
  {"R_rotulos",    (DL_FUNC) &R_rotulos,    1},
  {"R_componentes", (DL_FUNC) &R_componentes, 3},
  {"R_threads",    (DL_FUNC) &R_threads,    2},
  {"R_posto_fixo", (DL_FUNC) &R_posto_fixo, 4},
  {"R_posto_blocos", (DL_FUNC) &R_posto_blocos, 6},
  {"R_pegs",       (DL_FUNC) &R_pegs,      10},
  {"R_pegs_estrutura", (DL_FUNC) &R_pegs_estrutura, 3},
  {NULL, NULL, 0}
};

void R_init_uaitrem(DllInfo* dll) {
  R_registerRoutines(dll, NULL, metodos, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

}  // extern "C"
