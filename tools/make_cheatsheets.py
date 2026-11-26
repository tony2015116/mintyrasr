#!/usr/bin/env python3
"""Build pkgdown/assets/cheatsheet.html and cheatsheet-zh.html.

The CSS and the small R highlighter are shared with the mintyr cheatsheet
(same author, MIT). Run from the package root:  python3 tools/make_cheatsheets.py
"""
import pathlib

VERSION = "0.5.2"
ROOT = pathlib.Path(__file__).resolve().parents[1]
CSS = (ROOT / "tools" / "cheatsheet.css").read_text(encoding="utf-8")

FNS = ["asr_gblup", "asr_scenarios", "asr_compare", "ainverse", "H.inverse", "cv_scenarios", "nest_cv", "w2l_nest", "c2p_nest",
       "asreml", "vpredict", "lrt", "G.inverse", "G.matrix", "import_xlsx", "vm", "ide"]

SCRIPT = r"""
<script>
(function () {
  var fns = %s;
  var fnRe = new RegExp("^(" + fns.map(function (f) { return f.replace(".", "\\."); }).join("|") + ")\\b");
  function esc(s) { return s.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;"); }
  function hl(src) {
    var out = "", i = 0, m;
    while (i < src.length) {
      var rest = src.slice(i), ch = src[i];
      if (ch === "#") {
        var end = src.indexOf("\n", i); if (end < 0) end = src.length;
        out += '<span class="com">' + esc(src.slice(i, end)) + "</span>"; i = end; continue;
      }
      if (ch === '"') {
        var j = i + 1;
        while (j < src.length && src[j] !== '"') { if (src[j] === "\\") j++; j++; }
        out += '<span class="str">' + esc(src.slice(i, j + 1)) + "</span>"; i = j + 1; continue;
      }
      if ((m = rest.match(fnRe)) && !/[\w.]/.test(src[i - 1] || "")) {
        out += '<span class="fnm">' + m[0] + "</span>"; i += m[0].length; continue;
      }
      if ((m = rest.match(/^([A-Za-z_.][\w.]*)(\s*=)(?!=)/)) && !/[\w.$]/.test(src[i - 1] || "")) {
        out += '<span class="arg">' + m[1] + "</span>" + m[2]; i += m[0].length; continue;
      }
      if ((m = rest.match(/^\d+(\.\d+)?/)) && !/[\w.]/.test(src[i - 1] || "")) {
        out += '<span class="num">' + m[0] + "</span>"; i += m[0].length; continue;
      }
      if ((m = rest.match(/^[A-Za-z_.][\w.]*/))) { out += esc(m[0]); i += m[0].length; continue; }
      out += esc(ch); i++;
    }
    return out;
  }
  var blocks = document.querySelectorAll("pre code");
  for (var k = 0; k < blocks.length; k++) blocks[k].innerHTML = hl(blocks[k].textContent);
})();
</script>
""" % ("[" + ",".join('"%s"' % f for f in FNS) + "]")


def page(lang, t):
    return f"""<!doctype html>
<html lang="{lang}">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{t['title']}</title>
<meta name="description" content="{t['desc']}">
<style>
{CSS}
</style>
</head>
<body>
<div class="wrap">

<header class="top">
  <img src="reference/figures/logo.svg" alt="mintyrasr logo" onerror="this.style.display='none'">
  <div>
    <h1>{t['h1']} <span class="ver">v{VERSION}</span></h1>
    <p>{t['lead']}</p>
  </div>
  <nav class="links">
{t['links']}
  </nav>
</header>

<section class="flow" aria-label="workflow">
{t['flow']}
</section>

<main class="grid">
{t['panels']}
</main>

<footer>
  <span>mintyrasr {VERSION} &middot; MIT &middot; Guo Meng</span>
  <span>{t['footer']}</span>
</footer>
</div>
{SCRIPT}
</body>
</html>
"""


# ----------------------------------------------------------------------------- English
EN = dict(
    title="mintyrasr cheatsheet",
    desc="mintyrasr cheatsheet: GBLUP parameter estimation and cross-validation with ASReml-R on mintyr scenarios.",
    h1="mintyrasr cheatsheet",
    lead="GBLUP parameter estimation and cross-validation with <code>ASReml-R</code>, directly on the nested data built by <code>mintyr</code>: genetic parameters that follow the model, correct LRTs, PA / accuracy / slope per fold.",
    links="""    <a href="index.html">Home</a>
    <a href="reference/index.html">Reference</a>
    <a href="articles/mintyrasr.html">Get started</a>
    <a href="articles/models.html">Models</a>
    <a href="https://github.com/tony2015116/mintyrasr">GitHub</a>
    <a href="cheatsheet-zh.html" hreflang="zh" lang="zh">中文</a>
    <button type="button" onclick="window.print()">Print / PDF</button>""",
    flow="""  <div class="files"><div class="step">Inputs</div><div class="fns">pheno · ginv · F</div></div>
  <div class="io"><div class="step">1 · Scenarios</div><div class="fns">cv_scenarios()<br>nest_cv()</div></div>
  <div class="rs"><div class="step">2 · Full-data fit</div><div class="fns">varcomp · vc_index</div></div>
  <div class="st"><div class="step">3 · Parameters</div><div class="fns">h2 · rg · rp · LRT</div></div>
  <div class="cv"><div class="step">4 · CV fits</div><div class="fns">GEBV · rel · PA<br>accuracy · slope</div></div>
  <div class="io"><div class="step">5 · Output</div><div class="fns">summary · xlsx</div></div>""",
    footer='<a href="https://tony2015116.github.io/mintyrasr/">tony2015116.github.io/mintyrasr</a> &middot; <a href="https://tony2015116.github.io/mintyr/">mintyr</a> &middot; <a href="https://github.com/tony2015116/mintyrasr/issues">Issues</a>',
    panels="""
<section class="panel mt">
  <h2>Setup <small>inputs</small></h2>
  <div class="fn">
    <pre><code>pak::pak("tony2015116/mintyrasr")   # needs licensed ASReml-R 4
library(asreml); library(mintyr); library(mintyrasr)

ginv &lt;- G.inverse(G, sparseform = TRUE)$Ginv.sparse
inbreeding &lt;- data.frame(id = rownames(G), F = diag(G) - 1)</code></pre>
    <table>
      <tr><th>input</th><th>requirement</th></tr>
      <tr><td><code>pheno</code></td><td>one row per individual: id, group, fixed effects, traits</td></tr>
      <tr><td><code>ginv</code></td><td>inverse G (or A) with <code>attr(, "INVERSE") = TRUE</code>; otherwise <code>vm()</code> takes it as G (warning)</td></tr>
      <tr><td><code>inbreeding</code></td><td>data.frame <code>(id, F)</code> or named vector; <code>NULL</code> &rarr; F = 0</td></tr>
    </table>
  </div>
</section>

<section class="panel io">
  <h2><a href="articles/mintyrasr.html">Scenarios</a> <small>mintyr</small></h2>
  <div class="fn">
    <pre><code>sc &lt;- cv_scenarios(pheno, cols = c("adg", "bf"), id = "id",
                   group = "test_herd", v = 5, repeats = 5,
                   seed = 2026, out_type = "df")</code></pre>
    <table>
      <tr><th>scenario</th><th>model</th><th>train</th></tr>
      <tr><td><code>st_sp</code></td><td>1 trait, 1 group</td><td>target group</td></tr>
      <tr><td><code>st_mp</code></td><td>1 trait, groups combined</td><td>+ other groups</td></tr>
      <tr><td><code>mt_sp</code></td><td>2 traits, 1 group</td><td>target group</td></tr>
      <tr><td><code>mt_mp</code></td><td>2 traits, groups combined</td><td>+ other groups</td></tr>
      <tr><td><code>mt_bt</code></td><td>1 trait per group as 2 traits (G&times;E)</td><td>+ other group</td></tr>
    </table>
  </div>
</section>

<section class="panel rs">
  <h2><a href="reference/asr_scenarios.html">Run</a> <small>all scenarios or one</small></h2>
  <div class="fn">
    <div class="sig"><b><a href="reference/asr_scenarios.html">asr_scenarios</a></b>(sc, ginv, scenarios = c("st_sp", "st_mp", "mt_sp", "mt_mp"), dir = ".", tag = deparse(substitute(ginv)), ...)</div>
    <p>One workbook per scenario: <code>&lt;dir&gt;/&lt;scenario&gt;_&lt;tag&gt;.xlsx</code>. <code>"mt_bt"</code>: diag residual, no rp.</p>
  </div>
  <div class="fn">
    <div class="sig"><b><a href="reference/asr_gblup.html">asr_gblup</a></b>(param, cv = NULL, ginv, fixed = ~ 1, random = NULL, id = "id", traits = NULL, genetic_cov = "us", random_cov = "us", residual_cov = "us", cv_vc = "fixed", vpredict = NULL, inbreeding = NULL, cv_min_rel = -Inf, min_rel = 0.05, min_n = 10, neg_to_na = FALSE, test_rg = TRUE, test_rp = FALSE, rp = TRUE, na_x = "include", dense = "auto", maxit = 100, workspace = "auto", max_update = 2, nonconv = "skip", asreml_args = list(), asreml_trace = FALSE, on_error = "warn", keep_pred = TRUE, keep_gebv = TRUE, path = NULL, verbose = TRUE)</div>
    <pre><code>res &lt;- asr_scenarios(sc, ginv, fixed = ~ hhys,
                     inbreeding = inbreeding, dir = "results")
st &lt;- asr_gblup(sc$st_sp$param, sc$st_sp$cv, ginv,
                fixed = ~ hhys, path = "st_sp.xlsx")
# parameters only, per group / trait pair (cv = NULL)
asr_gblup(w2l_nest(pheno, cols = c("adg", "bf"), by = "test_herd"),
          ginv = ginv, fixed = ~ hhys)
# pedigree, genomic, single-step on the same folds
asr_compare(sc, list(ABLUP = ainv, GBLUP = ginv, ssGBLUP = hinv),
            fixed = ~ hhys)</code></pre>
    <div class="note"><code>param</code> and <code>cv</code> rows are matched through shared outer columns (<code>name</code>, <code>pairs</code>, <code>test_herd</code>&hellip;). Single / two-trait from <code>value</code> / <code>value1</code>+<code>value2</code> or <code>traits</code>.</div>
  </div>
</section>

<section class="panel st">
  <h2><a href="articles/models.html">Models</a> <small>random effects &amp; structures</small></h2>
  <div class="fn">
    <pre><code># 1 trait:  value ~ fixed,  ~ vm(id, relmat) + random,  ~ idv(units)
# 2 traits: cbind(a, b) ~ trait + trait:fixed,
#   ~ genetic_cov(trait):vm(id, relmat) + random_cov(trait):random,
#   ~ units:residual_cov(trait)
asr_gblup(..., random = ~ ide(id))            # repeated records: rep
asr_gblup(..., random = ~ litter, random_cov = "idh")
asr_gblup(..., random = ~ vm(dam, ginv) + ide(dam))  # maternal
asr_gblup(..., random = ~ vm(id, ginv):test_herd)    # GxE: rB
asr_gblup(..., genetic_cov = "corgh")         # us too heavy
asr_gblup(..., residual_cov = "diag", rp = FALSE)</code></pre>
    <table>
      <tr><th>argument</th><th>values</th></tr>
      <tr><td><code>genetic_cov</code>, <code>random_cov</code></td><td><code>us</code> <code>idh</code>/<code>diag</code> <code>corv</code> <code>corh</code> <code>corgh</code></td></tr>
      <tr><td><code>residual_cov</code></td><td><code>us</code> <code>idh</code>/<code>diag</code></td></tr>
      <tr><td><code>ginv</code></td><td>A&#8315;&sup1; <code>ainverse(ped)</code> &middot; G&#8315;&sup1; &middot; H&#8315;&sup1;; <code>dense = "auto"</code>; F from <code>ainverse()</code></td></tr>
      <tr><td><code>na_x</code></td><td><code>include</code> (NA covariate = 0) <code>fail</code> <code>omit</code></td></tr>
    </table>
  </div>
</section>

<section class="panel st">
  <h2><a href="articles/models.html">Genetic parameters</a> <small>by component name</small></h2>
  <div class="fn">
    <p>Each varcomp row is matched to <code>names(fit$vparameters)</code> (= <code>V</code> of <code>vpredict</code>) and classified as genetic / random / residual and var / cov / cor. Scale rows <code>units!R</code> (after <code>units!units</code>) and <code>units:trait!R</code> are skipped.</p>
    <table>
      <tr><th>parameter</th><th>formula</th></tr>
      <tr><td>h2<sub>i</sub></td><td>Va<sub>i</sub> / (all variances of trait i)</td></tr>
      <tr><td>rg</td><td>Cov<sub>a</sub> / &radic;(Va<sub>1</sub>Va<sub>2</sub>); cor* &rarr; the correlation</td></tr>
      <tr><td>rp</td><td>&Sigma; covariances / &radic;(Vp<sub>1</sub>Vp<sub>2</sub>); cor term &rarr; r&radic;(V<sub>1</sub>V<sub>2</sub>)</td></tr>
      <tr><td><code>prop_&lt;term&gt;</code></td><td>V<sub>term</sub> / Vp (c&sup2;, pe&sup2;, m&sup2;)</td></tr>
      <tr><td><code>rep</code> / <code>rB_X</code></td><td>(Va + Vpe) / Vp &middot; Va / (Va + V<sub>G&times;X</sub>)</td></tr>
    </table>
    <pre><code># us G + us R
h2_value1 = V1/(V1 + V5); rg = V2/sqrt(V1 * V3)
# corgh G + diag R
h2_nls = V2/(V2 + V5); rg = V1
# own formulas (check res$vc_index first)
asr_gblup(..., random = ~ ide(id),
  vpredict = list(rep = rep ~ (V1 + V2)/(V1 + V2 + V3)))</code></pre>
    <div class="note">fa() / rr() and several genetic terms (maternal <code>str()</code>): write formulas in <code>vpredict</code>.</div>
  </div>
</section>

<section class="panel ut">
  <h2>Tests <small>lrt(boundary = FALSE)</small></h2>
  <div class="fn">
    <table>
      <tr><th>argument</th><th>reduced model</th><th>output</th></tr>
      <tr><td><code>test_rg = TRUE</code></td><td>genetic <code>idh(trait)</code></td><td><code>rg_lrt</code> <code>rg_p</code></td></tr>
      <tr><td><code>test_rp = TRUE</code></td><td>all terms <code>idh</code></td><td><code>cov_df</code> <code>cov_lrt</code> <code>cov_p</code></td></tr>
    </table>
    <div class="note"><code>test_rp</code>: are the traits correlated at all (all covariances = 0), not rp = 0 alone.</div>
  </div>
</section>

<section class="panel cv">
  <h2><a href="articles/mintyrasr.html">Cross-validation</a> <small>per fold</small></h2>
  <div class="fn">
    <table>
      <tr><th><code>cv_vc</code></th><th>variance components in each fold</th></tr>
      <tr><td><code>"fixed"</code></td><td>full-data values, <code>Constraint = "F"</code>, <code>maxit = 1</code> &mdash; fast, slight leakage</td></tr>
      <tr><td><code>"warm"</code></td><td>re-estimated, starting from full-data values</td></tr>
      <tr><td><code>"estimate"</code></td><td>re-estimated from defaults &mdash; no leakage, slowest</td></tr>
    </table>
    <table>
      <tr><th>measure</th><th>definition</th></tr>
      <tr><td><code>rel</code></td><td>1 &minus; PEV / ((1 + F) Va), Va of the fold; <b>not filtered</b> in CV (<code>cv_min_rel = -Inf</code>); fold needs &ge; <code>min_n</code> = 10</td></tr>
      <tr><td><code>pa</code></td><td>r(y, GEBV)</td></tr>
      <tr><td><code>accuracy</code></td><td>pa / &radic;h2, full-data h2</td></tr>
      <tr><td><code>slope</code></td><td>b(y ~ GEBV), ideally 1</td></tr>
    </table>
  </div>
</section>

<section class="panel io">
  <h2>Output <small>one sheet each</small></h2>
  <div class="fn">
    <table>
      <tr><th>element</th><th>content</th></tr>
      <tr><td><code>varcomp</code></td><td>summary(fit)$varcomp + <code>V</code></td></tr>
      <tr><td><code>vc_index</code></td><td>V &rarr; term, role, kind, traits, used</td></tr>
      <tr><td><code>param</code></td><td>h2, rg, rp (+ SE), <code>formulas</code>, converge, loglik, boundary, tests</td></tr>
      <tr><td><code>accuracy</code></td><td>rep &times; fold &times; trait: status, n, pa, accuracy, slope, va, cv_iter</td></tr>
      <tr><td><code>summary</code></td><td>pa / accuracy / slope: n, mean, sd, se</td></tr>
      <tr><td><code>pred</code></td><td>validation: obs, gebv, pev, inb, rel</td></tr>
      <tr><td><code>gebv</code></td><td>all individuals (full data): gebv, pev, rel, <code>reliable</code> = rel &ge; <code>min_rel</code></td></tr>
      <tr><td><code>log</code></td><td>errors, non-convergence, boundary, invalid LRT</td></tr>
    </table>
  </div>
</section>

<section class="panel mt">
  <h2>Checks &amp; tips</h2>
  <div class="fn">
    <ul>
      <li>warnings: no <code>INVERSE</code> attribute; IDs not in <code>ginv</code>; components at the boundary (<code>B</code>, also <code>param$boundary</code>); <code>gammaPar = TRUE</code>; no convergence</li>
      <li><code>workspace = "auto"</code>: sized from individuals &times; traits; a fixed large value (e.g. 30Gb) slows every fit</li>
      <li>inside the model the matrix is <code>relmat</code>: components read <code>vm(id, relmat)</code>; in <code>random</code> write <code>vm(dam, ginv)</code></li>
      <li>one progress line per fit: status, iterations, time, ~time left (ASReml monitor off: <code>asreml_trace = FALSE</code>)</li>
      <li><code>max_update = 2</code> &rarr; <code>update()</code>; still not converged: <code>nonconv = "skip"</code> (CV rows skipped) or <code>"keep"</code></li>
      <li>errors &amp; non-convergence &rarr; <code>status</code> + <code>$log</code>, the loop goes on; rows with status &ne; ok are not in <code>summary</code></li>
      <li>keep the computer awake: sleep shows up as a long pause</li>
      <li>models are not kept &rarr; low memory</li>
      <li><code>asreml_args = list(ai.sing = TRUE)</code> for other <code>asreml()</code> arguments</li>
    </ul>
  </div>
</section>
""")

# ----------------------------------------------------------------------------- 中文
ZH = dict(
    title="mintyrasr 速查表",
    desc="mintyrasr 速查表：在 mintyr 场景数据上用 ASReml-R 完成 GBLUP 参数估计与交叉验证。",
    h1="mintyrasr 速查表",
    lead="在 <code>mintyr</code> 生成的嵌套数据上直接用 <code>ASReml-R</code> 完成 GBLUP 参数估计和交叉验证：遗传参数公式跟着模型走、检验正确、每折输出 PA / 准确性 / 回归系数。",
    links="""    <a href="index.html">主页</a>
    <a href="reference/index.html">函数参考</a>
    <a href="articles/mintyrasr-zh.html">入门</a>
    <a href="articles/models-zh.html">模型</a>
    <a href="https://github.com/tony2015116/mintyrasr">GitHub</a>
    <a href="cheatsheet.html" hreflang="en" lang="en">English</a>
    <button type="button" onclick="window.print()">打印 / PDF</button>""",
    flow="""  <div class="files"><div class="step">输入</div><div class="fns">表型 · ginv · F</div></div>
  <div class="io"><div class="step">1 · 场景</div><div class="fns">cv_scenarios()<br>nest_cv()</div></div>
  <div class="rs"><div class="step">2 · 全数据拟合</div><div class="fns">varcomp · vc_index</div></div>
  <div class="st"><div class="step">3 · 遗传参数</div><div class="fns">h2 · rg · rp · LRT</div></div>
  <div class="cv"><div class="step">4 · 交叉验证</div><div class="fns">GEBV · rel · PA<br>accuracy · slope</div></div>
  <div class="io"><div class="step">5 · 输出</div><div class="fns">summary · xlsx</div></div>""",
    footer='<a href="https://tony2015116.github.io/mintyrasr/">tony2015116.github.io/mintyrasr</a> &middot; <a href="https://tony2015116.github.io/mintyr/">mintyr</a> &middot; <a href="https://github.com/tony2015116/mintyrasr/issues">反馈问题</a>',
    panels="""
<section class="panel mt">
  <h2>准备 <small>输入</small></h2>
  <div class="fn">
    <pre><code>pak::pak("tony2015116/mintyrasr")   # 需要已授权的 ASReml-R 4
library(asreml); library(mintyr); library(mintyrasr)

ginv &lt;- G.inverse(G, sparseform = TRUE)$Ginv.sparse
inbreeding &lt;- data.frame(id = rownames(G), F = diag(G) - 1)</code></pre>
    <table>
      <tr><th>输入</th><th>要求</th></tr>
      <tr><td><code>pheno</code></td><td>每个个体一行：id、分组、固定效应、性状</td></tr>
      <tr><td><code>ginv</code></td><td>G（或 A）的逆，带 <code>attr(, "INVERSE") = TRUE</code>；否则 <code>vm()</code> 当作 G 本身（会警告）</td></tr>
      <tr><td><code>inbreeding</code></td><td>数据框 <code>(id, F)</code> 或命名向量；<code>NULL</code> &rarr; F = 0</td></tr>
    </table>
  </div>
</section>

<section class="panel io">
  <h2><a href="articles/mintyrasr-zh.html">场景</a> <small>mintyr</small></h2>
  <div class="fn">
    <pre><code>sc &lt;- cv_scenarios(pheno, cols = c("adg", "bf"), id = "id",
                   group = "test_herd", v = 5, repeats = 5,
                   seed = 2026, out_type = "df")</code></pre>
    <table>
      <tr><th>场景</th><th>模型</th><th>训练集</th></tr>
      <tr><td><code>st_sp</code></td><td>单性状、单群体</td><td>目标群体</td></tr>
      <tr><td><code>st_mp</code></td><td>单性状、多群体合并</td><td>+ 其他群体</td></tr>
      <tr><td><code>mt_sp</code></td><td>两性状、单群体</td><td>目标群体</td></tr>
      <tr><td><code>mt_mp</code></td><td>两性状、多群体合并</td><td>+ 其他群体</td></tr>
      <tr><td><code>mt_bt</code></td><td>同一性状在两个群体当作两个性状（G&times;E）</td><td>+ 另一群体</td></tr>
    </table>
  </div>
</section>

<section class="panel rs">
  <h2><a href="reference/asr_scenarios.html">运行</a> <small>全部场景或单个场景</small></h2>
  <div class="fn">
    <div class="sig"><b><a href="reference/asr_scenarios.html">asr_scenarios</a></b>(sc, ginv, scenarios = c("st_sp", "st_mp", "mt_sp", "mt_mp"), dir = ".", tag = deparse(substitute(ginv)), ...)</div>
    <p>每个场景一个工作簿：<code>&lt;dir&gt;/&lt;场景&gt;_&lt;tag&gt;.xlsx</code>。<code>"mt_bt"</code>：diag 残差，不算 rp。</p>
  </div>
  <div class="fn">
    <div class="sig"><b><a href="reference/asr_gblup.html">asr_gblup</a></b>(param, cv = NULL, ginv, fixed = ~ 1, random = NULL, id = "id", traits = NULL, genetic_cov = "us", random_cov = "us", residual_cov = "us", cv_vc = "fixed", vpredict = NULL, inbreeding = NULL, cv_min_rel = -Inf, min_rel = 0.05, min_n = 10, neg_to_na = FALSE, test_rg = TRUE, test_rp = FALSE, rp = TRUE, na_x = "include", dense = "auto", maxit = 100, workspace = "auto", max_update = 2, nonconv = "skip", asreml_args = list(), asreml_trace = FALSE, on_error = "warn", keep_pred = TRUE, keep_gebv = TRUE, path = NULL, verbose = TRUE)</div>
    <pre><code>res &lt;- asr_scenarios(sc, ginv, fixed = ~ hhys,
                     inbreeding = inbreeding, dir = "results")
st &lt;- asr_gblup(sc$st_sp$param, sc$st_sp$cv, ginv,
                fixed = ~ hhys, path = "st_sp.xlsx")
# 只估计参数：每个群体 / 每对性状（cv = NULL）
asr_gblup(w2l_nest(pheno, cols = c("adg", "bf"), by = "test_herd"),
          ginv = ginv, fixed = ~ hhys)
# 系谱、基因组、一步法在同一套分折上比较
asr_compare(sc, list(ABLUP = ainv, GBLUP = ginv, ssGBLUP = hinv),
            fixed = ~ hhys)</code></pre>
    <div class="note"><code>param</code> 和 <code>cv</code> 通过共同的外层列（<code>name</code>、<code>pairs</code>、<code>test_herd</code>&hellip;）对应。单/两性状由 <code>value</code> / <code>value1</code>+<code>value2</code> 或 <code>traits</code> 判断。</div>
  </div>
</section>

<section class="panel st">
  <h2><a href="articles/models-zh.html">模型</a> <small>随机效应与方差结构</small></h2>
  <div class="fn">
    <pre><code># 单性状：value ~ fixed，~ vm(id, relmat) + random，~ idv(units)
# 两性状：cbind(a, b) ~ trait + trait:fixed，
#   ~ genetic_cov(trait):vm(id, relmat) + random_cov(trait):random，
#   ~ units:residual_cov(trait)
asr_gblup(..., random = ~ ide(id))            # 重复记录：rep
asr_gblup(..., random = ~ litter, random_cov = "idh")
asr_gblup(..., random = ~ vm(dam, ginv) + ide(dam))  # 母体效应
asr_gblup(..., random = ~ vm(id, ginv):test_herd)    # G×E：rB
asr_gblup(..., genetic_cov = "corgh")         # us 难收敛时
asr_gblup(..., residual_cov = "diag", rp = FALSE)</code></pre>
    <table>
      <tr><th>参数</th><th>取值</th></tr>
      <tr><td><code>genetic_cov</code>、<code>random_cov</code></td><td><code>us</code> <code>idh</code>/<code>diag</code> <code>corv</code> <code>corh</code> <code>corgh</code></td></tr>
      <tr><td><code>residual_cov</code></td><td><code>us</code> <code>idh</code>/<code>diag</code></td></tr>
      <tr><td><code>ginv</code></td><td>A&#8315;&sup1; <code>ainverse(ped)</code> &middot; G&#8315;&sup1; &middot; H&#8315;&sup1;；<code>dense = "auto"</code>；F 取自 <code>ainverse()</code></td></tr>
      <tr><td><code>na_x</code></td><td><code>include</code>（缺失协变量 = 0）<code>fail</code> <code>omit</code></td></tr>
    </table>
  </div>
</section>

<section class="panel st">
  <h2><a href="articles/models-zh.html">遗传参数</a> <small>按组分名字定位</small></h2>
  <div class="fn">
    <p>varcomp 每一行与 <code>names(fit$vparameters)</code>（即 <code>vpredict</code> 的 <code>V</code>）对应，并分为 遗传 / 随机 / 残差 和 方差 / 协方差 / 相关。尺度参数 <code>units!R</code>（在 <code>units!units</code> 之后）和 <code>units:trait!R</code> 不使用。</p>
    <table>
      <tr><th>参数</th><th>公式</th></tr>
      <tr><td>h2<sub>i</sub></td><td>Va<sub>i</sub> / 性状 i 的所有方差之和</td></tr>
      <tr><td>rg</td><td>Cov<sub>a</sub> / &radic;(Va<sub>1</sub>Va<sub>2</sub>)；cor* 结构直接取相关</td></tr>
      <tr><td>rp</td><td>所有协方差之和 / &radic;(Vp<sub>1</sub>Vp<sub>2</sub>)；相关项按 r&radic;(V<sub>1</sub>V<sub>2</sub>) 换算</td></tr>
      <tr><td><code>prop_&lt;项&gt;</code></td><td>V<sub>项</sub> / Vp（c&sup2;、pe&sup2;、m&sup2;）</td></tr>
      <tr><td><code>rep</code> / <code>rB_X</code></td><td>(Va + Vpe) / Vp &middot; Va / (Va + V<sub>G&times;X</sub>)</td></tr>
    </table>
    <pre><code># us G + us R
h2_value1 = V1/(V1 + V5); rg = V2/sqrt(V1 * V3)
# corgh G + diag R
h2_nls = V2/(V2 + V5); rg = V1
# 自定义公式（先看 res$vc_index）
asr_gblup(..., random = ~ ide(id),
  vpredict = list(rep = rep ~ (V1 + V2)/(V1 + V2 + V3)))</code></pre>
    <div class="note">fa() / rr() 以及多个遗传项（母体效应 <code>str()</code>）：在 <code>vpredict</code> 中自己写公式。</div>
  </div>
</section>

<section class="panel ut">
  <h2>检验 <small>lrt(boundary = FALSE)</small></h2>
  <div class="fn">
    <table>
      <tr><th>参数</th><th>简化模型</th><th>输出</th></tr>
      <tr><td><code>test_rg = TRUE</code></td><td>遗传项 <code>idh(trait)</code></td><td><code>rg_lrt</code> <code>rg_p</code></td></tr>
      <tr><td><code>test_rp = TRUE</code></td><td>所有项 <code>idh</code></td><td><code>cov_df</code> <code>cov_lrt</code> <code>cov_p</code></td></tr>
    </table>
    <div class="note"><code>test_rp</code>：两性状是否有任何相关（所有协方差 = 0），不是单独检验 rp = 0。</div>
  </div>
</section>

<section class="panel cv">
  <h2><a href="articles/mintyrasr-zh.html">交叉验证</a> <small>每折</small></h2>
  <div class="fn">
    <table>
      <tr><th><code>cv_vc</code></th><th>每折的方差组分</th></tr>
      <tr><td><code>"fixed"</code></td><td>全数据估计值，<code>Constraint = "F"</code>，<code>maxit = 1</code> &mdash; 快，有轻微泄漏</td></tr>
      <tr><td><code>"warm"</code></td><td>以全数据值为初值重新估计</td></tr>
      <tr><td><code>"estimate"</code></td><td>从默认初值重新估计 &mdash; 无泄漏，最慢</td></tr>
    </table>
    <table>
      <tr><th>指标</th><th>定义</th></tr>
      <tr><td><code>rel</code></td><td>1 &minus; PEV / ((1 + F) Va)，Va 取当折；交叉验证中<b>不过滤</b>（<code>cv_min_rel = -Inf</code>）；每折至少 <code>min_n</code> = 10 个</td></tr>
      <tr><td><code>pa</code></td><td>r(y, GEBV)，预测能力</td></tr>
      <tr><td><code>accuracy</code></td><td>pa / &radic;h2，h2 用全数据估计值</td></tr>
      <tr><td><code>slope</code></td><td>b(y ~ GEBV)，理想值 1</td></tr>
    </table>
  </div>
</section>

<section class="panel io">
  <h2>输出 <small>每个元素一张表</small></h2>
  <div class="fn">
    <table>
      <tr><th>元素</th><th>内容</th></tr>
      <tr><td><code>varcomp</code></td><td>summary(fit)$varcomp + <code>V</code></td></tr>
      <tr><td><code>vc_index</code></td><td>V &rarr; 所属项、角色、类型、性状、是否使用</td></tr>
      <tr><td><code>param</code></td><td>h2、rg、rp（含 SE）、<code>formulas</code>、converge、loglik、boundary、检验</td></tr>
      <tr><td><code>accuracy</code></td><td>重复 &times; 折 &times; 性状：status、n、pa、accuracy、slope、va、cv_iter</td></tr>
      <tr><td><code>summary</code></td><td>pa / accuracy / slope：n、mean、sd、se</td></tr>
      <tr><td><code>pred</code></td><td>验证个体：obs、gebv、pev、inb、rel</td></tr>
      <tr><td><code>gebv</code></td><td>全数据所有个体：gebv、pev、rel、<code>reliable</code> = rel &ge; <code>min_rel</code></td></tr>
      <tr><td><code>log</code></td><td>出错、未收敛、边界、LRT 无效</td></tr>
    </table>
  </div>
</section>

<section class="panel mt">
  <h2>检查与提示</h2>
  <div class="fn">
    <ul>
      <li>警告：ginv 无 <code>INVERSE</code> 属性；个体不在 <code>ginv</code> 中；组分在边界上（<code>B</code>，也记在 <code>param$boundary</code>）；<code>gammaPar = TRUE</code>；未收敛</li>
      <li><code>workspace = "auto"</code>：按个体数 &times; 性状数估算；固定的大值（如 30Gb）会拖慢每次拟合</li>
      <li>模型内部矩阵叫 <code>relmat</code>：组分名为 <code>vm(id, relmat)</code>；<code>random</code> 中写 <code>vm(dam, ginv)</code></li>
      <li>每个模型一行进度：状态、迭代、用时、预计剩余（关闭 ASReml 迭代输出：<code>asreml_trace = FALSE</code>）</li>
      <li><code>max_update = 2</code> &rarr; <code>update()</code>；仍不收敛：<code>nonconv = "skip"</code>（跳过其交叉验证）或 <code>"keep"</code></li>
      <li>出错、未收敛 &rarr; 记 <code>status</code> 和 <code>$log</code>，循环继续；status &ne; ok 的行不进 <code>summary</code></li>
      <li>长时间运行请保持电脑唤醒：休眠表现为长时间停顿</li>
      <li>不保存模型对象 &rarr; 内存占用低</li>
      <li>其他 <code>asreml()</code> 参数：<code>asreml_args = list(ai.sing = TRUE)</code></li>
    </ul>
  </div>
</section>
""")

out = ROOT / "pkgdown" / "assets"
out.mkdir(parents=True, exist_ok=True)
(out / "cheatsheet.html").write_text(page("en", EN), encoding="utf-8")
(out / "cheatsheet-zh.html").write_text(page("zh-CN", ZH), encoding="utf-8")
print("written:", out / "cheatsheet.html", out / "cheatsheet-zh.html")
