## 准备 输入

    pak::pak("tony2015116/mintyrasr")   # 需要已授权的 ASReml-R 4
    library(asreml); library(mintyr); library(mintyrasr)

    ginv <- G.inverse(G, sparseform = TRUE)$Ginv.sparse
    inbreeding <- data.frame(id = rownames(G), F = diag(G) - 1)

| 输入 | 要求 |
|----|----|
| `pheno` | 每个个体一行：id、分组、固定效应、性状 |
| `ginv` | G（或 A）的逆，带 `attr(, "INVERSE") = TRUE`；否则 `vm()` 当作 G 本身（会警告） |
| `inbreeding` | 数据框 `(id, F)` 或命名向量；`NULL` → F = 0 |

## [场景](https://tony2015116.github.io/mintyrasr/articles/mintyrasr-zh.md) mintyr

    sc <- cv_scenarios(pheno, cols = c("adg", "bf"), id = "id",
                       group = "test_herd", v = 5, repeats = 5,
                       seed = 2026, out_type = "df")

| 场景    | 模型                                  | 训练集      |
|---------|---------------------------------------|-------------|
| `st_sp` | 单性状、单群体                        | 目标群体    |
| `st_mp` | 单性状、多群体合并                    | \+ 其他群体 |
| `mt_sp` | 两性状、单群体                        | 目标群体    |
| `mt_mp` | 两性状、多群体合并                    | \+ 其他群体 |
| `mt_bt` | 同一性状在两个群体当作两个性状（G×E） | \+ 另一群体 |

## [运行](https://tony2015116.github.io/mintyrasr/reference/asr_scenarios.md) 全部场景或单个场景

**[asr_scenarios](https://tony2015116.github.io/mintyrasr/reference/asr_scenarios.md)**(sc,
ginv, scenarios = c("st_sp", "st_mp", "mt_sp", "mt_mp"), dir = ".", tag
= deparse(substitute(ginv)), ...)

每个场景一个工作簿：`<dir>/<场景>_<tag>.xlsx`。`"mt_bt"`：diag
残差，不算 rp。

**[asr_gblup](https://tony2015116.github.io/mintyrasr/reference/asr_gblup.md)**(param,
cv, ginv, fixed = ~ 1, random = NULL, id = "id", traits = NULL,
genetic_cov = "us", random_cov = "us", residual_cov = "us", cv_vc =
"fixed", vpredict = NULL, inbreeding = NULL, cv_min_rel = -Inf, min_rel
= 0.05, min_n = 10, neg_to_na = FALSE, test_rg = TRUE, test_rp = FALSE,
rp = TRUE, na_x = "include", dense = TRUE, maxit = 100, workspace =
"30Gb", max_update = 2, nonconv = "skip", asreml_args = list(),
asreml_trace = FALSE, on_error = "warn", keep_pred = TRUE, keep_gebv =
TRUE, path = NULL, verbose = TRUE)

    res <- asr_scenarios(sc, ginv, fixed = ~ hhys,
                         inbreeding = inbreeding, dir = "results")
    st <- asr_gblup(sc$st_sp$param, sc$st_sp$cv, ginv,
                    fixed = ~ hhys, path = "st_sp.xlsx")
    # 自己的嵌套表、其他列名
    asr_gblup(all_nest, all_cv, ginv, fixed = ~ hhys,
              traits = c("bfs", "nls"))

`param` 和 `cv`
通过共同的外层列（`name`、`pairs`、`test_herd`…）对应。单/两性状由
`value` / `value1`+`value2` 或 `traits` 判断。

## [模型](https://tony2015116.github.io/mintyrasr/articles/models-zh.md) 随机效应与方差结构

    # 单性状：value ~ fixed，~ vm(id, ginv) + random，~ idv(units)
    # 两性状：cbind(a, b) ~ trait + trait:fixed，
    #   ~ genetic_cov(trait):vm(id, ginv) + random_cov(trait):random，
    #   ~ units:residual_cov(trait)
    asr_gblup(..., random = ~ ide(id))            # 重复记录
    asr_gblup(..., random = ~ litter, random_cov = "idh")
    asr_gblup(..., genetic_cov = "corgh")         # us 难收敛时
    asr_gblup(..., residual_cov = "diag", rp = FALSE)

| 参数                        | 取值                                     |
|-----------------------------|------------------------------------------|
| `genetic_cov`、`random_cov` | `us` `idh`/`diag` `corv` `corh` `corgh`  |
| `residual_cov`              | `us` `idh`/`diag`                        |
| `na_x`                      | `include`（缺失协变量 = 0）`fail` `omit` |

## [遗传参数](https://tony2015116.github.io/mintyrasr/articles/models-zh.md) 按组分名字定位

varcomp 每一行与 `names(fit$vparameters)`（即 `vpredict` 的
`V`）对应，并分为 遗传 / 随机 / 残差 和 方差 / 协方差 / 相关。尺度参数
`units!R`（在 `units!units` 之后）和 `units:trait!R` 不使用。

| 参数   | 公式                                               |
|--------|----------------------------------------------------|
| h2_(i) | Va_(i) / 性状 i 的所有方差之和                     |
| rg     | Cov_(a) / √(Va₁Va₂)；cor\* 结构直接取相关          |
| rp     | 所有协方差之和 / √(Vp₁Vp₂)；相关项按 r√(V₁V₂) 换算 |

    # us G + us R
    h2_value1 = V1/(V1 + V5); rg = V2/sqrt(V1 * V3)
    # corgh G + diag R
    h2_nls = V2/(V2 + V5); rg = V1
    # 自定义公式（先看 res$vc_index）
    asr_gblup(..., random = ~ ide(id),
      vpredict = list(rep = rep ~ (V1 + V2)/(V1 + V2 + V3)))

fa() / rr() 以及多个遗传项（母体效应 `str()`）：在 `vpredict`
中自己写公式。

## 检验 lrt.asreml(boundary = FALSE)

| 参数             | 简化模型            | 输出                       |
|------------------|---------------------|----------------------------|
| `test_rg = TRUE` | 遗传项 `idh(trait)` | `rg_lrt` `rg_p`            |
| `test_rp = TRUE` | 所有项 `idh`        | `cov_df` `cov_lrt` `cov_p` |

`test_rp`：两性状是否有任何相关（所有协方差 = 0），不是单独检验 rp = 0。

## [交叉验证](https://tony2015116.github.io/mintyrasr/articles/mintyrasr-zh.md) 每折

| `cv_vc` | 每折的方差组分 |
|----|----|
| `"fixed"` | 全数据估计值，`Constraint = "F"`，`maxit = 1` — 快，有轻微泄漏 |
| `"warm"` | 以全数据值为初值重新估计 |
| `"estimate"` | 从默认初值重新估计 — 无泄漏，最慢 |

| 指标 | 定义 |
|----|----|
| `rel` | 1 − PEV / ((1 + F) Va)，Va 取当折；交叉验证中**不过滤**（`cv_min_rel = -Inf`）；每折至少 `min_n` = 10 个 |
| `pa` | r(y, GEBV)，预测能力 |
| `accuracy` | pa / √h2，h2 用全数据估计值 |
| `slope` | b(y ~ GEBV)，理想值 1 |

## 输出 每个元素一张表

| 元素 | 内容 |
|----|----|
| `varcomp` | summary(fit)\$varcomp + `V` |
| `vc_index` | V → 所属项、角色、类型、性状、是否使用 |
| `param` | h2、rg、rp（含 SE）、`formulas`、converge、loglik、boundary、检验 |
| `accuracy` | 重复 × 折 × 性状：status、n、pa、accuracy、slope、va、cv_iter |
| `summary` | pa / accuracy / slope：n、mean、sd、se |
| `pred` | 验证个体：obs、gebv、pev、inb、rel |
| `gebv` | 全数据所有个体：gebv、pev、rel、`reliable` = rel ≥ `min_rel` |
| `log` | 出错、未收敛、边界、LRT 无效 |

## 检查与提示

- 警告：ginv 无 `INVERSE` 属性；个体不在 `ginv`
  中；组分在边界上（`B`，也记在
  `param$boundary`）；`gammaPar = TRUE`；未收敛
- ASReml 默认 `maxit = 13`、`workspace = "128mb"`，两性状 GBLUP 不够用
- 每个模型一行进度：状态、迭代、用时、预计剩余（关闭 ASReml
  迭代输出：`asreml_trace = FALSE`）
- `max_update = 2` →
  `update()`；仍不收敛：`nonconv = "skip"`（跳过其交叉验证）或 `"keep"`
- 出错、未收敛 → 记 `status` 和 `$log`，循环继续；status ≠ ok 的行不进
  `summary`
- 长时间运行请保持电脑唤醒：休眠表现为长时间停顿
- 不保存模型对象 → 内存占用低
- 其他 `asreml()` 参数：`asreml_args = list(ai.sing = TRUE)`
