# mintyrasr 入门

[English
version](https://tony2015116.github.io/mintyrasr/articles/mintyrasr.md)

ASReml-R 是商业软件，网站构建时不运行本页代码；在装有 ASReml-R 4
的电脑上，每段代码都可以照原样运行。

## 这个包做什么

比较基因组预测模型时，mintyr
为每个场景准备两张表：**参数估计数据**（全部记录）和**训练集 /
验证集**（每个“重复 × 折”一行），所有场景在同一批个体上验证。mintyrasr
接过这两张表，对每一行：

1.  用全部数据拟合模型：方差组分、h2；两性状再加 rg、rp 和 rg = 0
    的似然比检验；
2.  用每个训练集重新拟合（方差组分固定、以全数据为初值重估，或完全重估）；
3.  预测验证个体的 GEBV，计算可靠性、预测能力、准确性和回归系数；
4.  按折汇总，每个场景写一个 Excel 工作簿。

| 函数 | 用途 |
|----|----|
| `asr_gblup(param, cv, ginv, ...)` | 跑一个场景（或任意嵌套表 + 交叉验证表） |
| `asr_scenarios(sc, ginv, ...)` | 一次跑完 [`mintyr::cv_scenarios()`](https://tony2015116.github.io/mintyr/reference/cv_scenarios.html) 的所有场景 |

## 安装

``` r

# 先安装并激活 ASReml-R 4（https://vsni.co.uk/）
install.packages("mintyr")
pak::pak("tony2015116/mintyrasr")
```

## 1. 准备输入

需要三样东西：

- **表型**：每个个体一行，包括
  ID、分组列（场、品系……）、固定效应和性状；
- **`ginv`**：`vm()` 能用的基因组关系逆矩阵。必须带
  `attr(ginv, "INVERSE") <- TRUE`，否则 `vm()` 会把它当成 G
  本身（mintyrasr 会警告）；
- **近交系数**（可选）：数据框 `(id, F)`，用于可靠性
  `1 - PEV / ((1 + F) Va)`。GBLUP 常用 `F = diag(G) - 1`；不提供时 F =
  0。

``` r

library(asreml)
library(ASRgenomics)
library(mintyr)
library(mintyrasr)

pheno <- import_xlsx("pheno.xlsx")             # id, test_herd, hhys, adg, bf, fcr

G    <- G.matrix(M = snp, method = "VanRaden")$G
ginv <- G.inverse(G = G, sparseform = TRUE)$Ginv.sparse   # 自带 INVERSE = TRUE
inbreeding <- data.frame(id = rownames(G), F = diag(G) - 1)
```

## 2. 用 mintyr 生成场景

``` r

sc <- cv_scenarios(pheno, cols = c("adg", "bf", "fcr"), id = "id",
                   group = "test_herd", v = 5, repeats = 5, seed = 2026,
                   out_type = "df")      # ASReml-R 需要 data.frame
sc
```

| 场景 | 模型 | 参数估计数据 | 训练集 |
|----|----|----|----|
| `st_sp` | 单性状、单群体 | 性状 × 群体 | 目标群体 |
| `st_mp` | 单性状、多群体合并 | 性状 | 目标群体 + 其他所有群体 |
| `mt_sp` | 两性状、单群体 | 性状对 × 群体 | 目标群体 |
| `mt_mp` | 两性状、多群体合并 | 性状对 | 目标群体 + 其他所有群体 |
| `mt_bt` | 同一性状在不同群体当作不同性状 | 性状 | 目标群体 + 其他所有群体 |

## 3. 一次跑完所有场景

``` r

res <- asr_scenarios(sc, ginv, fixed = ~ hhys, inbreeding = inbreeding,
                     dir = "results")
# results/st_sp_ginv.xlsx、st_mp_ginv.xlsx、mt_sp_ginv.xlsx、mt_mp_ginv.xlsx

res$st_sp$summary
res$mt_sp$param
```

在 `scenarios` 中加入
`"mt_bt"`，就是把同一性状在两个群体当作两个性状（G×E）。没有个体同时有两个群体的记录，所以自动用
`diag` 残差，不算 rp：

``` r

res_bt <- asr_scenarios(sc, ginv, scenarios = "mt_bt", fixed = ~ hhys,
                        genetic_cov = "corgh", inbreeding = inbreeding)
res_bt$mt_bt$param        # 各群体 h2、群体间 rg（含 SE）、rg_p
```

## 4. 单独跑一个场景

[`asr_gblup()`](https://tony2015116.github.io/mintyrasr/reference/asr_gblup.md)
接收一张嵌套参数表和一张交叉验证表，二者通过共同的外层列（`name`、`pairs`、`test_herd`……）对应。

``` r

st <- asr_gblup(sc$st_sp$param, sc$st_sp$cv, ginv = ginv,
                fixed = ~ hhys, inbreeding = inbreeding,
                path = "results/st_sp_ginv.xlsx")
st
```

用
[`w2l_nest()`](https://tony2015116.github.io/mintyr/reference/w2l_nest.html)
/
[`c2p_nest()`](https://tony2015116.github.io/mintyr/reference/c2p_nest.html) +
[`nest_cv()`](https://tony2015116.github.io/mintyr/reference/nest_cv.html)
自己生成的表也可以。响应列不叫 `value` / `value1` + `value2` 时，用
`traits` 指定：

``` r

nested <- w2l_nest(pheno, cols = c("adg", "bf"), by = "test_herd", out_type = "df")
cvtab  <- nest_cv(nested, v = 5, repeats = 2, seed = 1, out_type = "df")
asr_gblup(nested, cvtab, ginv, fixed = ~ hhys)

asr_gblup(all_nest5, all_cv5, ginv, fixed = ~ hhys, traits = c("bfs", "nls"),
          genetic_cov = "corgh", residual_cov = "diag", rp = FALSE)
```

## 5. 输出内容

每个结果是一个列表，每个元素是工作簿中的一张表。

| 元素 | 内容 |
|----|----|
| `varcomp` | 每个参数估计模型的 `summary(fit)$varcomp`，附 `V` 编号 |
| `vc_index` | 每个 `V` 是哪个组分：所属项、角色（genetic / random / residual）、类型（var / cov / cor）、性状、是否使用 |
| `param` | h2（两性状分性状）、rg、rp 及 SE；实际用的 `formulas`；`status`、`converge`、`updates`、`iterations`、`loglik`、`boundary`；`rg_lrt`、`rg_p`（`test_rp = TRUE` 时还有 `cov_*`） |
| `accuracy` | 每个重复 × 折 × 性状：`status`、`n`、`pa`、`accuracy`、`slope`、`va`、`h2_full`、`cv_converge`、`cv_updates`、`cv_iter` |
| `summary` | `status == "ok"` 的折上 `pa`、`accuracy`、`slope` 的 n、mean、sd、se |
| `pred` | 验证个体的 `obs`、`gebv`、`pev`、`inb`（F）、`rel` |
| `gebv` | 全数据模型中所有个体的 `gebv`、`pev`、`inb`、`rel`、`reliable` |
| `log` | 每一次出错、未收敛、组分在边界上、LRT 无效的记录，含时间和标签 |

- `pa`（预测能力）= 验证个体的 r(y, GEBV)；
- `accuracy`（准确性）= pa / sqrt(h2)，h2 用全数据估计值（Dekkers
  2007）；
- `slope` = b(y ~ GEBV)，等于 1 表示 GEBV 既没有膨胀也没有压缩。

**交叉验证中不按可靠性过滤**（`cv_min_rel = -Inf`）。验证个体在训练集中没有自己的表型，可靠性偏低是正常的；只留下可靠性高的个体会让
PA 偏向亲缘联系强的个体，样本也会变少。阈值 `min_rel`（默认
0.05）改为用于全数据模型的最终 GEBV，即 `gebv$reliable`
标记。每折至少要有 `min_n`（默认 10）个验证个体，否则该折的
PA、准确性和回归系数为
`NA`。负值保留（`neg_to_na = FALSE`）：去掉负值会让均值偏高。

## 6. 交叉验证中的方差组分

| `cv_vc` | 每折的处理 | 说明 |
|----|----|----|
| `"fixed"`（默认） | 固定为全数据估计值（`Constraint = "F"`，`maxit = 1`） | 最快；验证个体的表型参与了方差估计，有轻微信息泄漏 |
| `"warm"` | 用训练集重新估计，以全数据值为初值 |  |
| `"estimate"` | 从默认初值重新估计 | 无泄漏，最慢 |

计算可靠性用的 `va` 总是取当折模型的值。

## 7. 长时间运行：进度、update() 与安全运行

每个模型一行进度：状态、迭代次数、用时和预计剩余时间；mintyrasr
运行期间关闭 ASReml 自己的迭代输出（`asreml_trace = TRUE` 可打开）：

    asr_gblup: two-trait model | 4 parameter fit(s) + 60 CV fit(s) | cv_vc = "fixed"
    [param 1/4] pairs=adg-bf test_herd=nls | ok | 9 it, 21s
    [param 2/4] pairs=adg-bf test_herd=bfs | ok after 1 update() | 31 it, 48s
    [cv 12/60  20%] pairs=adg-bf test_herd=nls rep=Repeat2 fold=Fold2 | ok | 2 it, 4s | 6.1 min elapsed, ~19.5 min left
    Done in 31.2 min | param: 4 ok, 0 not converged, 0 failed | cv: ok 60 | 1 log entry

- **[`update()`](https://rdrr.io/r/stats/update.html)**：未收敛的模型从当前估计值继续迭代，最多
  `max_update = 2` 次（记在 `param$updates`、`accuracy$cv_updates`）。
- **仍不收敛**（`nonconv = "skip"`，默认）：参数估计模型标为
  `status = "not_converged"`，其交叉验证行跳过（`status = "param_not_converged"`）；交叉验证模型则该折为
  `NA`。`nonconv = "keep"` 照常使用这些模型。
- **出错**：该次拟合标为
  `status = "failed"`，循环继续（`on_error = "stop"` 立即停止）。
- **`$log`** 记录每一次出错、未收敛、组分在边界上、LRT
  无效，含时间、阶段和标签；结束时给出一条警告，说明有多少模型失败或未收敛。
- `status` 不是 `"ok"` 的行不进入 `summary`。

其他提示：

- 不保存模型对象，折数很多时内存也不会爆。
- 两性状 GBLUP 需要比 ASReml
  默认（`maxit = 13`、`workspace = "128mb"`）更大的设置；mintyrasr 默认
  `maxit = 100`、`workspace = "30Gb"`。
- `na_x = "include"`（默认，与常见脚本一致）把缺失协变量当作
  0，数值协变量请先中心化；ASReml 自身默认是 `"fail"`。
- 其他 `asreml()` 参数：`asreml_args = list(ai.sing = TRUE)`。
- 长时间运行时让电脑保持唤醒：电脑休眠（或在线许可证等待）表现为长时间停顿，而不是拟合慢。

下一步：[模型与遗传参数](https://tony2015116.github.io/mintyrasr/articles/models-zh.md)——随机效应、方差结构、V
编号如何确定、检验以及自定义 `vpredict()` 公式。
