# BCIntrog 0.1.0

## Initial release

- 实现 `read_vcf()`：读取 VCF/VCF.GZ，提取 CHROM、POS、REF、ALT、GT、GQ。
- 实现 `select_site()`：筛选双亲相反纯合的双等位标记，按 GQ 过滤，判定后代 D/R/H。
- 实现 `emission_probability()` 和 `transition_probability()`。
- 实现 `run_hmm()`：三状态 HMM 的 Viterbi 和 Forward-Backward，缺失观测作为无信息处理。
- 实现 `segment_detection()`：按染色体识别片段，使用 1-based 闭区间，处理小片段，最多迭代 2 轮。`segments` 表只包含 `CHROM`、`segment_id`、`Start`、`End`、`Length_bp`、`SNP_number`、`Observed_SNP_number`、`Viterbi` 和 `State`。
- 实现 `summary.donor_segments()`：供体基因组当量、轮回亲本背景恢复率和按状态/染色体汇总。
- 实现 `visualization()`：复现示例风格的全基因组标记和片段图。
- 实现 `run_analysis()`：一键运行完整流程，可选写出结果文件。
- 添加 5 个 vignettes、单元测试和完整 Demo 集成测试。
- `R CMD check`：0 errors、0 warnings、0 notes。
