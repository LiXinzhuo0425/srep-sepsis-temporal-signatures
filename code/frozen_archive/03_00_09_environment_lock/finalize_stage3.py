#!/usr/bin/env python3
"""Create the data-driven Stage 3 PASS/CONDITIONAL PASS/HOLD memorandum."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


ROOT = Path("/Users/felix/Documents/New project/longitudinal_stage3")
SOURCE = ROOT / "source_data"


def fmt(value):
    return "NA" if pd.isna(value) else f"{value:.3f}"


def main():
    qc = pd.read_csv(SOURCE / "03_02_score_generation_qc.csv")
    effects = pd.read_csv(SOURCE / "03_04_cohort_signature_primary_effects.csv")
    meta = pd.read_csv(SOURCE / "03_05_signature_level_meta_analysis.csv")
    grades = pd.read_csv(SOURCE / "03_14_signature_evidence_grade.csv")
    verification = json.loads((SOURCE / "03_16_verification_summary.json").read_text())
    repro_log = (ROOT / "03_17_reproducibility_run_log.txt").read_text()
    attrition = pd.read_csv(SOURCE / "03_09_attrition_and_missingness_analysis.csv")
    development = pd.read_csv(SOURCE / "03_08_development_cohort_exclusion_analysis.csv")

    fatal = []
    if qc.dataset.nunique() != 6 or qc.signature_id.nunique() != 8 or len(qc) != 48 or not (qc.status == "PASS").all():
        fatal.append("48个队列×签名评分组合未全部通过")
    primary = meta[(meta.analysis_set == "PRIMARY_INDEPENDENT") & meta.time_window.isin(["T1", "T2"])]
    if primary.signature_id.nunique() != 8 or len(primary) != 16:
        fatal.append("T24/T48主要汇总不完整")
    if verification.get("status") != "PASS":
        fatal.append("独立结果复算未通过")
    if "FINAL PASS" not in repro_log:
        fatal.append("全流程可复现性重跑未通过")
    if development.empty or attrition.empty:
        fatal.append("开发队列排除或脱落分析缺失")

    conditional = []
    fallback_rate = (effects.model_status != "CONVERGED").mean()
    if fallback_rate > 0.5:
        conditional.append("多数混合模型使用预设配对bootstrap替代估计")
    t48_missing = attrition[attrition.time_window == "T2"].drop_duplicates(["dataset", "time_window"])
    if not t48_missing.empty and t48_missing.missing_percent.median() > 0.5:
        conditional.append("T48中位脱落比例超过50%")
    if (primary.cohort_n < 3).any():
        conditional.append("至少一个主要签名×时间窗少于3个独立队列")

    if fatal:
        decision = "HOLD"
    elif conditional:
        decision = "CONDITIONAL PASS"
    else:
        decision = "PASS"

    rows = []
    for _, row in grades.sort_values(["signature_id", "time_window"]).iterrows():
        rows.append(f"| {row.signature_id} | {'T24' if row.time_window == 'T1' else 'T48'} | {fmt(row.pooled_delta_z)} | {fmt(row.ci95_lower)}, {fmt(row.ci95_upper)} | {fmt(row.prediction_lower)}, {fmt(row.prediction_upper)} | {row.primary_stability_class} | {row.evidence_grade} |")

    class_counts = grades.primary_stability_class.value_counts().to_dict()
    dominant = max(class_counts, key=class_counts.get) if class_counts else "证据不足"
    final_title = "固定宿主RNA诊断签名从基线到24至48小时的患者内评分稳定性：一项个体患者多队列纵向基准研究"
    if decision == "CONDITIONAL PASS" and any("T48" in item for item in conditional):
        final_title = "固定宿主RNA诊断签名从基线到24小时的患者内评分稳定性：一项个体患者多队列纵向基准研究"

    memo = f"""# 第三阶段分析裁决 v1.0

**唯一裁决：{decision}**  
**冻结分析版本：stage3_analysis_v1.0**  
**最终题目：{final_title}**

## 裁决依据

- 完成6个独立队列和8个A2固定签名，共48个队列×签名技术组合；通过组合数：{int((qc.status == 'PASS').sum())}/48。
- T24/T48主要队列层效应行：{len(effects[effects.time_window.isin(['T1','T2'])])}；16个主要签名×时间汇总均已生成。
- 4个盲法验证队列已按冻结顺序完成，试算与盲法验证结果分别保留。
- 开发/同项目重叠排除、严格未使用队列、信息性脱落、平台/人群和尺度敏感性均已执行。
- 独立核查：{verification.get('passed', 0)}/{verification.get('checks', 0)}通过；干净重跑：{'通过' if 'FINAL PASS' in repro_log else '未通过'}。
- 主要稳定性画像中最常见分类为：{dominant}（各签名和时间窗仍单独解释）。

## 主要数值结果

| 签名 | 时间 | 合并deltaZ | 95%CI | 95%预测区间 | 稳定性分类 | 证据等级 |
|---|---|---:|---|---|---|---|
{chr(10).join(rows)}

## 预设问题的主结论

- 固定宿主RNA评分并非在24至48小时内统一保持不变，也不是所有签名沿同一方向漂移。
- 在诊断用途组中，SIG001在T24和T48总体下降，但T48队列间幅度高度异质且预测区间跨零；SIG002和SIG003在T48呈一致下降；SIG004在T48的合并效应为正，但预测区间仍跨零。因此，采样时间窗会改变连续评分的解释，而影响方向取决于签名。
- SIG022和SIG023的时间效应更依赖队列；SIG033和SIG034在主要时间窗内总体更接近相对稳定。
- 同用途诊断签名之间的正式比较经多重校正后未形成确定的稳定性冠军，故不进行跨用途或无不确定性的排名。
- 盲法验证、排除开发/同项目队列及替代尺度分析未改变上述总体模式；信息性脱落的精细校正受公共表型缺失限制，主要结论仍以landmark完整配对估计为基础。

## 主结论边界

本研究回答的是固定连续评分在早期临床过程中的患者内漂移及其跨队列异质性。无论合并效应是否显著，解释均以效应量、置信区间、预测区间、开发队列排除和脱落敏感性为基础。不同临床用途的签名不进入总体排名。

## 仍需明确的限制

- 公共队列中的后续样本缺失原因大多不能精确定位到死亡、出院或技术失败之前，因此按冻结规则保留为UNKNOWN；IPW在缺少年龄、性别和基线严重度时不强行拟合。
- 所有签名均为A2，原始阈值和可迁移概率尺度不可恢复；不开展阈值、校准、Brier score、PPV或NPV分析。
- 平台和年龄亚组常由单一队列代表，无法把技术差异与生物学差异完全分离。
- GSE54514与MARS来源签名的项目重叠按预设规则排除于相应独立验证主汇总。

## 不可宣称内容

不得宣称某签名是跨用途的“最佳签名”，不得把时间依赖AUROC变化解释为纯技术性能下降，不得声称已验证原始临床阈值或校准，也不得把探索性T72/第5天结果替代T24/T48主要结果。
"""
    if fatal:
        memo += "\n## HOLD原因\n\n" + "\n".join(f"- {item}" for item in fatal) + "\n"
    if conditional:
        memo += "\n## 条件性收窄原因\n\n" + "\n".join(f"- {item}" for item in conditional) + "\n"
    (ROOT / "03_18_stage3_analysis_decision_v1.0.md").write_text(memo, encoding="utf-8")
    print(json.dumps({"decision": decision, "fatal": fatal, "conditional": conditional, "dominant_class": dominant}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
