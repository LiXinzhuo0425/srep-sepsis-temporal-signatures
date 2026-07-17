import fs from "node:fs/promises";
import path from "node:path";
const artifactModule = process.env.ARTIFACT_TOOL_MJS;
if (!artifactModule) throw new Error("Set ARTIFACT_TOOL_MJS to artifact_tool.mjs");
const { SpreadsheetFile, Workbook } = await import(artifactModule);

const ROOT = process.env.SEPSIS_SIGNATURE_ANALYSIS_ROOT || process.cwd();
const SOURCE = path.join(ROOT, "04_source_data");
const QA = path.join(ROOT, "04_qa", "workbooks");
const NAVY = "#17365D";
const BLUE = "#DCE6F1";
const PALE = "#F4F7FA";
const BORDER = "#C8D1DC";
const TEAL = "#0F6B6D";

const specs = [
  {
    output: "04_02_evidence_level_and_claim_boundary.xlsx",
    title: "Stage 4 evidence levels and claim boundaries",
    status: "FROZEN",
    note: "Defines the interpretation permissions for E1-E5 evidence. E2-E4 findings cannot be promoted to causal or clinical-utility claims.",
    sheets: [["EvidenceLevels", "04_02_evidence_level_and_claim_boundary.csv"]],
  },
  {
    output: "04_03_signature_decomposition_specification.xlsx",
    title: "Exact signature-decomposition specification",
    status: "FROZEN BEFORE DECOMPOSITION",
    note: "Contains the formula structure, exact attribution method and frozen signature-gene coefficients. Nonlinear signatures use exact Shapley or Owen attribution, not linear approximation.",
    sheets: [["DecompositionSpec", "04_03_signature_decomposition_specification.csv"], ["GeneCoefficients", "04_03_signature_gene_coefficients.csv"]],
  },
  {
    output: "04_04_gene_level_data_qc.xlsx",
    title: "Gene-level longitudinal dataset quality control",
    status: "PASS",
    note: "Cross-file and reconstruction checks for the gene-level dataset. The complete patient-level reconstruction ledger is included for traceability.",
    sheets: [["GeneLevelQC", "04_04_gene_level_data_qc.csv"], ["ReconstructionLedger", "04_05_decomposition_reconstruction_qc.csv"]],
  },
  {
    output: "04_06_cohort_level_gene_contributions.xlsx",
    title: "Cohort-level exact gene contributions",
    status: "COMPLETE",
    note: "All component genes are retained. Contribution estimates are mathematical sources of score change and are not interpreted as causal biological drivers.",
    sheets: [["GeneContributions", "04_06_cohort_level_gene_contributions.csv"], ["ArchitectureMetrics", "04_06_cohort_level_signature_architecture_metrics.csv"]],
  },
  {
    output: "04_07_gene_contribution_meta_analysis.xlsx",
    title: "Cross-cohort gene-contribution synthesis",
    status: "COMPLETE",
    note: "Random-effects results include uncertainty, heterogeneity, prediction intervals, development-overlap exclusions and leave-one-dataset-out diagnostics.",
    sheets: [["GeneMetaAnalysis", "04_07_gene_contribution_meta_analysis.csv"], ["LeaveOneOut", "04_07_gene_contribution_leave_one_out.csv"]],
  },
  {
    output: "04_08_signature_drift_architecture.xlsx",
    title: "Signature drift architecture",
    status: "FROZEN-RULE CLASSIFICATION",
    note: "Classifications were assigned from prespecified dominance, cancellation, heterogeneity and direction-consistency rules. This is not a best-to-worst ranking.",
    sheets: [["DriftArchitecture", "04_08_signature_drift_architecture.csv"]],
  },
  {
    output: "04_09_pathway_feasibility_audit.xlsx",
    title: "Pathway-module feasibility audit",
    status: "START GATE PASSED",
    note: "The pathway module was activated because at least four independent cohorts supported the prespecified gene sets and time windows without systematic platform failure.",
    sheets: [["ExpressionAudit", "04_09_full_transcriptome_expression_audit.csv"], ["GeneCoverage", "04_09_pathway_gene_coverage.csv"], ["WindowSupport", "04_09_pathway_time_window_support.csv"]],
  },
  {
    output: "04_11_longitudinal_pathway_changes.xlsx",
    title: "Longitudinal pathway changes",
    status: "COMPLETE",
    note: "All prespecified pathways are reported. Primary singscore results are accompanied by ssGSEA sensitivity results and leave-one-dataset-out diagnostics.",
    sheets: [["CohortChanges", "04_11_cohort_longitudinal_pathway_changes.csv"], ["MetaSingscore", "04_11_meta_longitudinal_pathway_changes.csv"], ["MetaSsGSEA", "04_11_meta_longitudinal_pathway_changes_ssgsea.csv"], ["LeaveOneOut", "04_11_pathway_leave_one_out.csv"]],
  },
  {
    output: "04_12_signature_pathway_coupling.xlsx",
    title: "Signature-pathway longitudinal coupling",
    status: "COMPLETE WITH BOUNDED E2 CLAIMS",
    note: "Associations describe concordant within-patient change. They do not establish that a pathway causes signature drift or that the score monitors treatment response.",
    sheets: [["CohortSingscore", "04_12_cohort_signature_pathway_coupling.csv"], ["MetaSingscore", "04_12_meta_signature_pathway_coupling.csv"], ["CohortSsGSEA", "04_12_cohort_signature_pathway_coupling_ssgsea.csv"], ["MetaSsGSEA", "04_12_meta_signature_pathway_coupling_ssgsea.csv"]],
  },
  {
    output: "04_13_cell_source_annotation.xlsx",
    title: "Potential blood-cell source annotation",
    status: "ANNOTATION COMPLETE",
    note: "Annotations are based on external reference expression and indicate possible source lineages only. Bulk RNA cannot separate abundance shifts from within-cell transcriptional changes.",
    sheets: [["GeneAnnotation", "04_13_cell_source_annotation.csv"], ["LeadingSources", "04_14_leading_cell_source_by_signature.csv"], ["CellSourceMap", "04_14_gene_contribution_cell_source_map.csv"]],
  },
  {
    output: "04_15_clinical_anchor_availability_audit.xlsx",
    title: "Clinical-anchor availability audit",
    status: "FORMAL CLINICAL GATE FAILED",
    note: "No common time-matched Level 1 clinical anchor was available in at least two cohorts and approximately 80 patients. This decision precedes and prevents association testing.",
    sheets: [["AvailabilityAudit", "04_15_clinical_anchor_availability_audit.csv"]],
  },
  {
    output: "04_17_signature_clinical_anchor_analysis.xlsx",
    title: "Clinical-anchor analysis status",
    status: "STOPPED BY PREDEFINED GATE",
    note: "Zero clinical association models and zero P values were generated. Mortality and static severity labels were not substituted for concurrent severity change.",
    sheets: [["AnalysisStatus", "04_17_signature_clinical_anchor_analysis_status.csv"]],
  },
  {
    output: "04_18_context_of_use_evidence_matrix.xlsx",
    title: "Context-of-use evidence matrix",
    status: "COMPLETE",
    note: "Separates supported temporal descriptions from unsupported threshold, prognostic, treatment-monitoring and clinical-decision claims.",
    sheets: [["ContextOfUse", "04_18_context_of_use_evidence_matrix.csv"]],
  },
  {
    output: "04_20_integrated_signature_interpretation_matrix.xlsx",
    title: "Integrated signature interpretation matrix",
    status: "COMPLETE",
    note: "Integrates Stage 3 temporal effects with exact gene architecture, limited pathway evidence, potential cell sources, development-overlap sensitivity and claim boundaries.",
    sheets: [["IntegratedMatrix", "04_20_integrated_signature_interpretation_matrix.csv"]],
  },
  {
    output: "04_21_multiplicity_and_robustness_audit.xlsx",
    title: "Multiplicity and robustness audit",
    status: "COMPLETE",
    note: "Documents each inferential family, multiplicity rule and required robustness check. Unadjusted findings are not used as main claims.",
    sheets: [["Audit", "04_21_multiplicity_and_robustness_audit.csv"], ["FigureExportQC", "04_figure_export_qc.csv"]],
  },
  {
    output: "04_22_independent_extension_verification.xlsx",
    title: "Independent Stage 4 extension verification",
    status: "38 OF 38 CHECKS PASS",
    note: "Independent checks cover exact reconstruction, formula rosters, protected source hashes, pathway gates, meta-analysis inputs, multiplicity, cell annotation, clinical stopping rules, figure exports and workbook integrity.",
    sheets: [["VerificationChecks", "04_22_independent_extension_verification.csv"]],
  },
];

function parseCsv(text) {
  const rows = [];
  let row = [], field = "", quoted = false;
  for (let i = 0; i < text.length; i++) {
    const ch = text[i];
    if (quoted) {
      if (ch === '"' && text[i + 1] === '"') { field += '"'; i++; }
      else if (ch === '"') quoted = false;
      else field += ch;
    } else {
      if (ch === '"') quoted = true;
      else if (ch === ',') { row.push(field); field = ""; }
      else if (ch === '\n') { row.push(field); rows.push(row); row = []; field = ""; }
      else if (ch !== '\r') field += ch;
    }
  }
  if (field.length || row.length) { row.push(field); rows.push(row); }
  return rows;
}

function coerce(value) {
  const v = value.trim();
  if (v === "") return "";
  if (/^(true|false)$/i.test(v)) return v.toLowerCase() === "true";
  if (/^(nan|na|null|none)$/i.test(v)) return "";
  if (/^[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?$/.test(v) && !/^0\d+$/.test(v)) {
    const n = Number(v);
    if (Number.isFinite(n)) return n;
  }
  return value;
}

async function loadCsv(filename) {
  const text = await fs.readFile(path.join(SOURCE, filename), "utf8");
  const rows = parseCsv(text);
  return rows.map((r, i) => i === 0 ? r : r.map(coerce));
}

function styleHeader(range) {
  range.format.fill = NAVY;
  range.format.font = { bold: true, color: "#FFFFFF", size: 10 };
  range.format.wrapText = true;
  range.format.rowHeightPx = 34;
  range.format.borders = { preset: "all", style: "thin", color: BORDER };
}

function setColumnWidths(sheet, matrix) {
  const rows = matrix.slice(0, Math.min(matrix.length, 151));
  for (let c = 0; c < matrix[0].length; c++) {
    let maxLen = 8;
    for (const row of rows) maxLen = Math.max(maxLen, String(row[c] ?? "").length);
    const width = Math.max(72, Math.min(260, 12 + maxLen * 6.2));
    sheet.getRangeByIndexes(0, c, matrix.length, 1).format.columnWidthPx = width;
  }
}

function styleDataSheet(sheet, matrix) {
  const rows = matrix.length, cols = matrix[0].length;
  sheet.showGridLines = false;
  sheet.freezePanes.freezeRows(1);
  const used = sheet.getRangeByIndexes(0, 0, rows, cols);
  used.format.font = { name: "Aptos", size: 9, color: "#1F2937" };
  used.format.verticalAlignment = "center";
  used.format.borders = { preset: "all", style: "thin", color: BORDER };
  if (rows > 1) {
    const body = sheet.getRangeByIndexes(1, 0, rows - 1, cols);
    body.format.wrapText = true;
    body.format.fill = "#FFFFFF";
  }
  styleHeader(sheet.getRangeByIndexes(0, 0, 1, cols));
  setColumnWidths(sheet, matrix);
  const headers = matrix[0].map(String);
  for (let c = 0; c < cols; c++) {
    const h = headers[c].toLowerCase();
    const col = sheet.getRangeByIndexes(1, c, Math.max(1, rows - 1), 1);
    if (/(^n$|count|number|patients|cohorts|genes|rows|windows)/.test(h)) col.setNumberFormat("0");
    else if (/(p_value|fdr|rho|effect|estimate|contribution|ratio|index|mean|median|sd|se$|ci_|pi_|tau|i2|coverage|proportion|error)/.test(h)) col.setNumberFormat("0.0000");
  }
}

function addReadme(workbook, spec) {
  const sheet = workbook.worksheets.add("README");
  sheet.showGridLines = false;
  sheet.getRange("A1:F1").merge();
  sheet.getRange("A1").values = [[spec.title]];
  sheet.getRange("A1:F1").format.fill = NAVY;
  sheet.getRange("A1:F1").format.font = { bold: true, color: "#FFFFFF", size: 16 };
  sheet.getRange("A1:F1").format.rowHeightPx = 34;
  const rows = [
    ["Stage", "4 - Biological interpretation and clinical anchoring"],
    ["Version", "v1.0"],
    ["Date", "2026-07-16"],
    ["Status", spec.status],
    ["Purpose", spec.note],
    ["Scope lock", "Does not change Stage 3 signatures, patients, time windows, score directions or principal conclusions."],
    ["Claim boundary", "Mathematical decomposition is E1. Pathway coupling is E2. Cell annotation is potential-source evidence. No causal, prognostic, threshold or treatment-decision claim is permitted."],
    ["Source data", spec.sheets.map((x) => `04_source_data/${x[1]}`).join("\n")],
    ["Data sheets", spec.sheets.map((x) => x[0]).join(", ")],
  ];
  sheet.getRangeByIndexes(2, 0, rows.length, 2).values = rows;
  sheet.getRangeByIndexes(2, 0, rows.length, 1).format.fill = BLUE;
  sheet.getRangeByIndexes(2, 0, rows.length, 1).format.font = { bold: true, color: NAVY };
  sheet.getRangeByIndexes(2, 0, rows.length, 2).format.borders = { preset: "all", style: "thin", color: BORDER };
  sheet.getRangeByIndexes(2, 0, rows.length, 2).format.wrapText = true;
  sheet.getRangeByIndexes(2, 0, rows.length, 2).format.verticalAlignment = "top";
  sheet.getRange("A3:A11").format.columnWidthPx = 125;
  sheet.getRange("B3:B11").format.columnWidthPx = 650;
  sheet.getRange("A3:B11").format.rowHeightPx = 34;
  sheet.getRange("A7:B10").format.rowHeightPx = 62;
  sheet.freezePanes.freezeRows(1);
}

async function build(spec) {
  const workbook = Workbook.create();
  addReadme(workbook, spec);
  let previewRange = "A1:F25";
  for (const [sheetName, filename] of spec.sheets) {
    const matrix = await loadCsv(filename);
    const sheet = workbook.worksheets.add(sheetName.slice(0, 31));
    sheet.getRangeByIndexes(0, 0, matrix.length, matrix[0].length).values = matrix;
    styleDataSheet(sheet, matrix);
    if (filename === spec.sheets[0][1]) {
      let col = matrix[0].length, letters = "";
      while (col > 0) { col--; letters = String.fromCharCode(65 + (col % 26)) + letters; col = Math.floor(col / 26); }
      previewRange = `A1:${letters}${Math.min(matrix.length, 26)}`;
    }
  }
  const inspect = await workbook.inspect({ kind: "sheet", include: "id,name" });
  await fs.writeFile(path.join(QA, "inspect", `${spec.output}.json`), JSON.stringify(inspect, null, 2));
  const preview = await workbook.render({ sheetName: spec.sheets[0][0].slice(0, 31), range: previewRange, scale: 0.9, format: "png" });
  await fs.writeFile(path.join(QA, "previews", `${spec.output}.png`), new Uint8Array(await preview.arrayBuffer()));
  const output = await SpreadsheetFile.exportXlsx(workbook);
  await output.save(path.join(ROOT, spec.output));
}

await fs.mkdir(path.join(QA, "inspect"), { recursive: true });
await fs.mkdir(path.join(QA, "previews"), { recursive: true });
for (const spec of specs) await build(spec);
console.log(JSON.stringify({ workbooks: specs.length, outputs: specs.map((s) => s.output) }, null, 2));
