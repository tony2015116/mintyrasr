# 模型与遗传参数

[English
version](https://tony2015116.github.io/mintyrasr/articles/models.md)

## 模型

**单性状**（嵌套数据有 `value` 列，或 `traits = "adg"`）：

``` text
asreml(value ~ <fixed>,
       random   = ~ vm(id, ginv) + <random>,
       residual = ~ idv(units), dense = ~ vm(id, ginv), ...)
```

**两性状**（`value1` + `value2`，或 `traits = c("a", "b")`）：

``` text
asreml(cbind(value1, value2) ~ trait + trait:<fixed>,
       random   = ~ <genetic_cov>(trait):vm(id, ginv) + <random_cov>(trait):<random>,
       residual = ~ units:<residual_cov>(trait), dense = ~ <genetic_cov>(trait):vm(id, ginv), ...)
```

| 参数 | 取值 | 默认 |
|----|----|----|
| `fixed` | 单边公式，如 `~ hhys + sex` | `~ 1` |
| `random` | 其他随机效应，如 `~ ide(id)`、`~ litter + pen` | 无 |
| `genetic_cov` | `"us"`、`"idh"` / `"diag"`、`"corv"`、`"corh"`、`"corgh"` | `"us"` |
| `random_cov` | 同上，作用于 `random` 中的每一项 | `"us"` |
| `residual_cov` | `"us"`、`"idh"` / `"diag"` | `"us"` |
| `na_x` | `na.method(x = )`：`"include"`、`"fail"`、`"omit"` | `"include"` |
| `asreml_args` | 其他 `asreml()` 参数 | [`list()`](https://rdrr.io/r/base/list.html) |

常见写法：

``` r

# 重复记录：加永久环境效应，否则高估 h2
asr_gblup(param, cv, ginv, fixed = ~ hhys, random = ~ ide(id))

# 窝效应，性状间不相关
asr_gblup(param, cv, ginv, fixed = ~ hhys, random = ~ litter, random_cov = "idh")

# us() 难收敛时用相关参数化
asr_gblup(param, cv, ginv, fixed = ~ hhys, genetic_cov = "corgh")

# 同一性状在两个场：没有个体同时有两个记录 -> diag 残差，不算 rp
asr_gblup(param, cv, ginv, fixed = ~ hhys, traits = c("bfs", "nls"),
          genetic_cov = "corgh", residual_cov = "diag", rp = FALSE)
```

## 为什么不能写死 `V1 / (V1 + V2)`？

`vpredict()` 中的 `V1`、`V2`……按模型 `vparameters`
的顺序编号，这个顺序随模型变化：

| 模型                 | V1  | V2    | V3  | V4                  | V5  | V6              | V7  |
|----------------------|-----|-------|-----|---------------------|-----|-----------------|-----|
| us G + us R          | Va1 | Cov_a | Va2 | `units:trait!R` = 1 | Ve1 | Cov_e           | Ve2 |
| corgh G + diag R     | r_a | Va1   | Va2 | `units:trait!R` = 1 | Ve1 | Ve2             | \-  |
| us G + idh 窝 + us R | Va1 | Cov_a | Va2 | Vc1                 | Vc2 | `units:trait!R` | ……  |

`rg ~ V2 / sqrt(V1 * V3)` 对第一个模型是对的；对第二个模型是错的（V1
已经是相关系数），带 `V7` 的 rp 公式还会报错。所以 mintyrasr
读取方差组分的**名字**，自动生成公式。

## 如何识别方差组分

`summary(fit)$varcomp` 的每一行先和 `names(fit$vparameters)` 对应（得到
`V` 编号），再分类：

| 名字（示例） | 所属项 | 角色 | 类型 |
|----|----|----|----|
| `vm(id, ginv)` | `vm(id, ginv)` | genetic | var |
| `ide(id)`、`litter` | 本身 | random | var |
| `units!units` | `units` | residual | var |
| idv 中 `units!units` 之后的 `units!R`、`units:trait!R` |  |  | 固定为 1 的尺度参数，不使用 |
| `trait:vm(id, ginv)!trait_adg:adg` | `trait:vm(id, ginv)` | genetic | var |
| `trait:vm(id, ginv)!trait_bf:adg` |  | genetic | cov |
| `trait:litter!trait_adg` | `trait:litter` | random | var（idh） |
| `...!trait!bf:!trait!adg.cor`、`...!trait!cor` |  |  | cor（corgh、corh、corv） |
| `...!trait!var` |  |  | 两性状共用的方差（corv） |
| `fa(...)`、`rr(...)` 项 |  |  | 不解析（给出警告） |

性状名在组分名中按完整单词匹配，任何列名都可以。结果中的 `vc_index`
表列出了分类，新模型第一次运行后请看一眼。

## 公式

对每个性状 *i*，Vp_i = 该性状**所有**方差之和（遗传 + 其他随机效应 +
残差）：

- **h2_i** = Va_i / Vp_i
- **rg** = Cov_a / sqrt(Va1 Va2)；`corv` / `corh` / `corgh`
  时直接取相关参数
- **rp** = 所有协方差之和 / sqrt(Vp1 Vp2)；相关参数化的项按 r sqrt(V1
  V2) 换算成协方差

自动生成的例子（`param$formulas`）：

    # 单性状，G + e
    h2 = V1/(V1 + V2)
    # 单性状，G + pe + e
    h2 = V1/(V1 + V2 + V3)
    # us G + us R
    h2_value1 = V1/(V1 + V5); rg = V2/sqrt(V1 * V3); rp = (V2 + V6)/sqrt((V1 + V5) * (V3 + V7))
    # corgh G + diag R，同一性状在两个场
    h2_nls = V2/(V2 + V5); h2_bfs = V3/(V3 + V6); rg = V1

## 自定义公式

其他参数——重复力、窝效应比例、母体遗传力或总遗传力——用 `vpredict`
以命名列表给出。先在 `vc_index` 中确认 `V` 编号，它随模型变化。

``` r

res <- asr_gblup(param, cv, ginv, fixed = ~ hhys, random = ~ ide(id),
                 vpredict = list(rep = rep ~ (V1 + V2) / (V1 + V2 + V3)))
res$vc_index
res$param[, c("h2", "h2_se", "rep", "rep_se")]
```

有多个遗传项的模型（如母体效应的
[`str()`](https://rdrr.io/r/utils/str.html) 模型）以及 `fa()` / `rr()`
结构需要自己写公式：自动公式只使用 `vm(<id>, ginv)` 这一项。

## 检验

| 参数 | 完整模型 | 简化模型 | `lrt.asreml()` |
|----|----|----|----|
| `test_rg = TRUE`（默认） | 按设定 | 遗传结构换成 `idh(trait)`，其余不变 | `boundary = FALSE`（相关为 0 不在参数边界上） |
| `test_rp = TRUE` | 按设定 | 遗传、其他随机效应和残差都换成 `idh` | `boundary = FALSE`；`cov_df` = 协方差个数 |

`test_rp` 检验的是两个性状是否有任何相关，不是单独对 rp = 0
的检验。结果在 `param` 的 `rg_lrt`、`rg_p`、`cov_df`、`cov_lrt`、`cov_p`
列。

## 检查与警告

- `ginv` 没有 `attr(, "INVERSE") = TRUE`：`vm()` 会把它当作 G。
- 有表型的个体不在 `ginv` 中。
- 方差组分被固定在边界上（bound 为 `B`）：给出警告，并记在
  `param$boundary`。
- `asreml.options(gammaPar = TRUE)`：自动公式基于默认的方差（sigma）参数化。
- 模型未收敛（`max_update` 可调用
  [`update()`](https://rdrr.io/r/stats/update.html)）。
