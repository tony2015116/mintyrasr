## 用 mintyr + mintyrasr 直接在 R 里跑 ASReml GBLUP（无需导出数据文件）
library(asreml)
library(mintyr)
library(mintyrasr)

# 1. 表型：每行一个个体；ginv：G 逆矩阵；inbreeding：两列 (id, F)
# pheno <- import_xlsx("pheno.xlsx")
# ginv  <- ...
# inbreeding <- data.frame(id = rownames(G), F = diag(G) - 1)   # GBLUP 常用

# 2. 一次生成 4 个场景的参数估计数据 + 交叉验证数据（同一批验证个体）
sc <- cv_scenarios(pheno, cols = c("adg", "bf", "fcr"), id = "id",
                   group = "test_herd", v = 5, repeats = 5, seed = 2026,
                   scenarios = c("st_sp", "st_mp", "mt_sp", "mt_mp"),
                   out_type = "df")

# 3a. 一次跑完 4 个场景，输出 st_sp_ginv.xlsx / st_mp_ginv.xlsx / mt_sp_ginv.xlsx / mt_mp_ginv.xlsx
res <- asr_scenarios(sc, ginv, fixed = ~ hhys, inbreeding = inbreeding, dir = "results")

# 3b. 或者单独跑某一个场景（等价于原来的 st_blup(..., sp = TRUE)）
res_st_sp <- asr_gblup(sc$st_sp$param, sc$st_sp$cv, ginv = ginv,
                       fixed = ~ hhys, inbreeding = inbreeding,
                       path = "results/st_sp_ginv.xlsx")
res_st_sp$summary

# 3c. 原来的 all_nest / all_cv（w2l_nest / c2p_nest + nest_cv）也能直接用
# asr_gblup(all_nest1, all_cv1, ginv = ginv, fixed = ~ hhys, inbreeding = inbreeding)

# 4. 模型变了（加随机效应、换结构），h2 / rg / rp 公式自动按组分名字生成
res_pe <- asr_gblup(sc$st_sp$param, sc$st_sp$cv, ginv = ginv, fixed = ~ hhys,
                    random = ~ ide(id),                                   # 重复记录：永久环境
                    vpredict = list(rep = rep ~ (V1 + V2) / (V1 + V2 + V3)),  # 编号先看 vc_index
                    inbreeding = inbreeding)
res_pe$vc_index        # 每个 V 对应哪个方差组分
res_pe$param$formulas  # 实际使用的 vpredict 公式

res_corgh <- asr_gblup(sc$mt_sp$param, sc$mt_sp$cv, ginv = ginv, fixed = ~ hhys,
                       genetic_cov = "corgh",                  # us 太重时
                       random = ~ litter, random_cov = "idh",
                       inbreeding = inbreeding)

# 5. CV 中方差组分：固定（默认，快）/ 以全数据为初值重估 / 完全重估（无泄漏）
res_noleak <- asr_gblup(sc$st_sp$param, sc$st_sp$cv, ginv = ginv, fixed = ~ hhys,
                        cv_vc = "estimate", inbreeding = inbreeding)
res_noleak$summary     # pa = r(y, GEBV)；accuracy = pa / sqrt(h2)；slope

# 6. 同一性状在两个场（test_herd = nls / bfs）当作两个性状：每个性状一组结果，corgh 直接给出 rg 及其 SE
sc_bt <- cv_scenarios(pheno, cols = c("bf", "adg"), id = "id", group = "test_herd",
                      v = 5, repeats = 5, seed = 2026, scenarios = "mt_bt", out_type = "df")
res_bt <- asr_scenarios(sc_bt, ginv, scenarios = "mt_bt", fixed = ~ hhys,
                        genetic_cov = "corgh", inbreeding = inbreeding)
res_bt$mt_bt$param     # h2_nls、h2_bfs、rg（= corgh 相关参数）、rg_p（corgh vs idh，boundary = FALSE）

# 7. 已有嵌套数据、响应列不叫 value1/value2 时
# asr_gblup(all_nest5, all_cv5, ginv, fixed = ~ hhys, traits = c("bfs", "nls"),
#           genetic_cov = "corgh", residual_cov = "diag", rp = FALSE)
