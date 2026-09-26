"""Reproducible reports and an explicit handoff to later, real CFD observations."""
from __future__ import annotations

import csv
import json
import platform
from pathlib import Path

import numpy as np

from . import __version__
from .core import (Baseline, CONTEXT, REFERENCE, DATA_DIR, INPUTS, TARGETS, TERMS,
                   QuadraticRSM, best_on_slice, endpoint_contrasts, in_support,
                   load_k1, loo_diagnostics, metrics)

OBS_FIELDS = ["case_id", "allocation_id", "context_id", "phase_id", "mesh_id", "reference_id",
              *INPUTS, "delta_eta_pp", "delta_torque_nm", "sealing_effectiveness",
              "ingress_kg_s", "status", "origin", "role"]


def write_csv(path, rows, fields=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows and fields is None:
        raise ValueError("Empty CSV requires explicit column names.")
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields or list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf-8")


def get_model():
    ids, x, y, digest = load_k1()
    return ids, x, y, digest, QuadraticRSM.fit(x, y)


def prepare_output(path):
    out = Path(path)
    out.mkdir(parents=True, exist_ok=True)
    # Refuse accidental overwrite of a previous run or any existing user files.
    if any(out.iterdir()):
        raise ValueError(f"Output directory is not empty: {out}. Choose a new run directory.")
    return out


def reproduce(out, plots=True):
    out = prepare_output(out)
    ids, x, y, digest, model = get_model()
    predicted = model.predict(x)
    folds = loo_diagnostics(ids, x, y)
    valid = [r for r in folds if r["status"] == "ok"]
    report = {
        "package_version": __version__, "python_version": platform.python_version(),
        "numpy_version": np.__version__, "source_sha256": digest,
        "n_samples": len(ids), "n_coefficients_per_target": len(TERMS), "rank": model.rank,
        "condition_number": float(model.singular_values[0]/model.singular_values[-1]),
        "in_sample": {t: metrics(y[:, j], predicted[:, j]) for j, t in enumerate(TARGETS)},
        "loo": {"valid_folds": len(valid), "total_folds": len(folds),
                "undefined_cases": [r["case_id"] for r in folds if r["status"] != "ok"],
                "complete_loo_score_available": len(valid) == len(folds),
                "valid_folds_only": {t: metrics([r[f"observed_{t}"] for r in valid],
                                                [r[f"predicted_{t}"] for r in valid]) for t in TARGETS}},
        "independent_cfd_validation": "not_performed",
        "warning": "Fit residuals and LOO are not CFD numerical uncertainty or out-of-condition validation."
    }
    write_json(out/"metrics.json", report)
    write_json(out/"model.json", {"terms": TERMS, "targets": TARGETS,
                                 "coefficients": model.coefficients.tolist(), "source_sha256": digest})
    write_csv(out/"coefficients.csv", [{"term": term, **{t: float(model.coefficients[i,j])
               for j,t in enumerate(TARGETS)}} for i, term in enumerate(TERMS)])
    write_csv(out/"reproduction.csv", [{"case_id": case, **dict(zip(INPUTS, x[i])),
              **{f"observed_{t}": y[i,j] for j,t in enumerate(TARGETS)},
              **{f"predicted_{t}": predicted[i,j] for j,t in enumerate(TARGETS)}}
              for i, case in enumerate(ids)])
    write_csv(out/"loo_diagnostics.csv", folds)
    write_csv(out/"endpoint_contrasts.csv", endpoint_contrasts(x,y))
    profile = []
    for k, variable in enumerate(INPUTS):
        points = np.zeros((101,4)); points[:,k] = np.linspace(-1,1,101)
        result = model.predict(points)
        for p, v in zip(points, result):
            profile.append({"variable": variable, "coded_value": p[k],
                            **dict(zip(TARGETS,v)), "kind": "RSM prediction, others at design"})
    write_csv(out/"one_factor_profiles.csv", profile)
    eta = report["in_sample"][TARGETS[0]]
    text = f"""# K1 공개 표 재현 보고서

- 출처: DOI 10.5293/kfma.2025.28.2.007, Table 2, p.11.
- 25개 공개 CFD 결과, 각 출력마다 15계수 2차식. 새 CFD 계산이 아닙니다.
- 효율 변화의 단위는 %가 아니라 **percentage point (%p)**입니다. 덜 음수일수록 좋습니다.
- 효율 학습 R²: {eta['r2']:.9f}; RMSE: {eta['rmse']:.9f} %p.
- 독립 CFD 검증: **아직 수행하지 않았습니다**.

## 중요한 교차검증 한계
일반적인 leave-one-out 25회 중 {len(valid)}회만 15계수 모델을 고유하게 정할 수 있습니다.
유일한 중앙점(case 00)을 빼면 나머지 점에서 x3²+x4²+x5²+x6²=2이므로,
상수항과 제곱항의 합이 구분되지 않습니다. 따라서 해당 fold에는 임의의 예측값을 채우지 않습니다.
metrics.json의 valid_folds_only 수치는 전체 LOO 점수도, 독립 물리 검증도 아닙니다.

## 민감도 해석
endpoint_contrasts.csv는 다른 조건이 같은 행들의 +1 minus -1 차이입니다.
동일한 kg/s 변화가 아니라 각각 자기 기준유량의 50%→150% 변화입니다.
④의 U자형 곡률은 양 끝점 차이 하나에 드러나지 않으므로 one_factor_profiles.csv도 보세요.
⑤·⑥을 작은 영향이라는 이유만으로 항상 무시할 수 있다고 결론 내리지 않습니다.

## 다음 단계
실제 m3,0와 m4,0가 제공되면 고정 총유량 단면을 계산합니다.
지금 데이터에는 ε, ingress, 금속온도, 가스별 결과가 없습니다.
이 값을 RSM에서 만들어 내지 않으며, 새 CFD/실험으로 별도 확보해야 합니다.
"""
    (out/"report.md").write_text(text, encoding="utf-8")
    if plots:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        for j,t in enumerate(TARGETS):
            fig, ax = plt.subplots(figsize=(6,5))
            ax.scatter(y[:,j], predicted[:,j], label="Published points, in-sample fit")
            extent = [float(y[:,j].min()),float(y[:,j].max())]
            ax.plot(extent, extent, linestyle="--", label="Identity")
            ax.set(xlabel=f"Published {t}", ylabel=f"Fitted {t}", title="Reproduction, NOT independent validation")
            ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(out/f"parity_{t}.png",dpi=150); plt.close(fig)
        fig, ax = plt.subplots(figsize=(7,5))
        for k,v in enumerate(INPUTS):
            p = np.zeros((101,4)); p[:,k] = np.linspace(-1,1,101)
            ax.plot(p[:,k],model.predict(p)[:,0],label=v)
        ax.set(xlabel="Coded flow (-1 = 50%, +1 = 150% of own design flow)",
               ylabel="Predicted delta efficiency [percentage points]",
               title="One-factor RSM profiles; other flows at design")
        ax.legend(); fig.tight_layout(); fig.savefig(out/"one_factor_eta.png",dpi=150); plt.close(fig)
    return report


def sweep(config, out, budget_scale=1., points=201, plots=True):
    if not isinstance(points,int) or points < 3 or points > 100001:
        raise ValueError("points must be an integer from 3 to 100001.")
    baseline = Baseline.load(config)
    lo, hi = baseline.alpha_bounds(budget_scale)
    _, _, _, digest, model = get_model()
    a = np.linspace(lo, hi, points)
    m3, m4, x = baseline.allocation(a, budget_scale)
    predicted = model.predict(x)
    nominal = model.predict([0,0,0,0])[0]
    optimum, score = best_on_slice(model, baseline, budget_scale)
    out = prepare_output(out)
    total6 = baseline.total6
    if total6 is not None:
        total6 += baseline.m34*(budget_scale-1)
    summary = {"source_kind": baseline.source_kind, "provenance": baseline.provenance,
               "source_sha256": digest, "baseline_config": json.loads(Path(config).read_text(encoding="utf-8")),
               "budget_scale": budget_scale, "m34_kg_s": baseline.m34*budget_scale,
               "k1_six_passage_total_kg_s": total6, "alpha_bounds": [lo,hi],
               "design_alpha": baseline.alpha0, "rsm_eta_maximizer_alpha": optimum,
               "predicted_delta_eta_pp": score, "predicted_gain_vs_design_pp": score-float(nominal[0]),
               "sealing_constraint": "NOT_EVALUATED_NO_DATA",
               "meaning": "Unconstrained RSM candidate, not a validated or thermally safe turbine design"}
    rows = []
    for i, alpha in enumerate(a):
        rows.append({"alpha": alpha, "m3_kg_s": m3[i], "m4_kg_s": m4[i],
                     **dict(zip(INPUTS,x[i])), **dict(zip(TARGETS,predicted[i])),
                     "gain_vs_design_eta_pp": predicted[i,0]-nominal[0],
                     "inside_design_hull": bool(in_support(x[i])[0]), "source_kind": baseline.source_kind})
    write_json(out/"allocation_summary.json",summary)
    write_csv(out/"allocation_sweep.csv",rows)
    # Five-point screen plus baseline and exact surrogate optimum, merged if coincident.
    candidates = [("screen",float(v)) for v in np.linspace(lo,hi,5)] + [("rsm_eta_best",optimum)]
    if lo <= baseline.alpha0 <= hi:
        candidates.append(("baseline" if abs(budget_scale-1)<1e-12 else "design_split", baseline.alpha0))
    selected = []
    for label, alpha in candidates:
        match = next((r for r in selected if abs(r["alpha"]-alpha)<1e-10),None)
        if match is not None:
            match["purpose"] += "+"+label
            continue
        f3,f4,coords = baseline.allocation(alpha,budget_scale)
        value = model.predict(coords)[0]
        selected.append({"case_id": f"A{len(selected):02d}", "purpose": label, "alpha": alpha,
                         "m3_kg_s": float(f3[0]),"m4_kg_s":float(f4[0]),
                         **dict(zip(INPUTS,coords[0])),
                         "predicted_delta_eta_pp":float(value[0]),
                         "predicted_delta_torque_nm":float(value[1]),
                         "source_kind":baseline.source_kind,"sealing_status":"unknown"})
    write_csv(out/"suggested_cases.csv",selected)
    pending = []
    for row in selected:
        pending.append({"case_id":row["case_id"],"allocation_id":row["case_id"],"context_id":CONTEXT,
                        "phase_id":"paper_zero_offset", "mesh_id":"paper_medium", "reference_id":REFERENCE,
                        **{v:row[v] for v in INPUTS}, **{v:None for v in TARGETS},
                        "sealing_effectiveness":None,"ingress_kg_s":None,
                        "status":"pending","origin":"not_run","role":"unassigned"})
    write_csv(out/"cfd_results_template.csv",pending,OBS_FIELDS)
    (out/"report.md").write_text(f"""# ③·④ 배분 탐색

자료 상태: **{baseline.source_kind}**. {baseline.provenance}

m3+m4={baseline.m34*budget_scale:.8g} kg/s인 단면에서, 각각 자기 기준유량의 50~150% 안만 탐색했습니다.
허용 alpha=[{lo:.8g}, {hi:.8g}], RSM의 공력 후보 alpha={optimum:.8g}.
예측 효율 변화={score:.8g} %p. 이는 CFD로 검증한 개선량이 아닙니다.

①②⑤⑥은 고정입니다. 실제 여섯 유량이 모두 제공되지 않으면 여섯 유로 총량은 null입니다.
budget_scale=1일 때 기준 총량도 유지됩니다. 다른 값이면 총량 변화 실험입니다.

## 실제 연구로 넘길 파일
suggested_cases.csv의 predicted 열은 계산 예약에 참고할 추정치입니다.
cfd_results_template.csv의 관측 결과 칸은 의도적으로 비어 있습니다.
실제 계산 후 status=completed, origin=cfd, role=holdout 또는 verification으로 지정하고,
물리조건/격자/위상/기준 정의를 확인한 뒤 실제 관측값을 넣으세요.
이 파일은 CFX 실행 입력도, 자동 CFX 실행기도 아닙니다.

ε 제약을 확인하지 않았으므로 최종 설계/안전 배분이라고 부르지 마세요.
""",encoding="utf-8")
    if plots:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7,5))
        ax.plot(a,predicted[:,0],label="RSM, no sealing constraint")
        ax.axvline(optimum,linestyle="--",label="Surrogate candidate")
        ax.set(xlabel="alpha = m3 / (m3 + m4)",ylabel="Predicted delta efficiency [percentage points]",
               title=f"{baseline.source_kind.upper()} baseline; not CFD validation")
        ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(out/"allocation_eta.png",dpi=150); plt.close(fig)
    return summary


def read_observations(path):
    rows = list(csv.DictReader(Path(path).open(encoding="utf-8-sig")))
    if not rows:
        raise ValueError("No observation rows.")
    completed = []
    seen = set()
    for row in rows:
        if row.get("status") != "completed":
            continue
        if row.get("origin") != "cfd" or row.get("role") not in ("holdout","verification"):
            raise ValueError("Completed observations need origin=cfd and role=holdout/verification.")
        for field in ("case_id","allocation_id","context_id","phase_id","mesh_id","reference_id"):
            if not row.get(field,"").strip():
                raise ValueError(f"Missing observation metadata: {field}.")
        key = tuple(row[k] for k in ("case_id","context_id","phase_id","mesh_id","reference_id"))
        if key in seen:
            raise ValueError("Duplicate observation key.")
        seen.add(key)
        for field in (*INPUTS,*TARGETS,"sealing_effectiveness","ingress_kg_s"):
            value = row.get(field,"")
            if field in INPUTS or field == "delta_eta_pp" or value not in (None,""):
                value = float(value)
                if not np.isfinite(value):
                    raise ValueError(f"Non-finite observation {field}.")
            else:
                value = None
            row[field] = value
        if row["sealing_effectiveness"] is not None and not 0<=row["sealing_effectiveness"]<=1:
            raise ValueError("Sealing effectiveness must be in [0,1].")
        if row["ingress_kg_s"] is not None and row["ingress_kg_s"]<0:
            raise ValueError("Gross ingress must be non-negative.")
        completed.append(row)
    if not completed:
        raise ValueError("No completed real CFD observations; template predictions are not observations.")
    return completed


def validate_cfd(results, out):
    observations = read_observations(results)
    _, training_x, _, _, model = get_model()
    errors = []
    for row in observations:
        if row["context_id"] != CONTEXT or row["reference_id"] != REFERENCE or row["phase_id"] != "paper_zero_offset" or row["mesh_id"] != "paper_medium":
            raise ValueError("Gas/phase/mesh/reference differs from K1. Use paired comparison; do not silently pool conditions.")
        x = np.array([row[k] for k in INPUTS])
        training_location = bool(np.any(np.all(np.isclose(training_x,x,atol=1e-10,rtol=0),axis=1)))
        if row["role"] == "holdout" and training_location:
            raise ValueError("A published training location cannot be labeled a new-location holdout.")
        prediction = model.predict(x)[0]
        item = {"case_id":row["case_id"],"role":row["role"],"training_location":training_location}
        for j,t in enumerate(TARGETS):
            item[f"observed_{t}"] = row[t]
            item[f"predicted_{t}"] = float(prediction[j])
            item[f"error_{t}"] = float(prediction[j]-row[t]) if row[t] is not None else None
        errors.append(item)
    out = prepare_output(out)
    write_csv(out/"cfd_errors.csv",errors)
    summary = {"n_completed":len(errors),"scores":{}}
    for role in ("holdout","verification"):
        subset = [r for r in errors if r["role"]==role]
        summary["scores"][role]={"n":len(subset)}
        for t in TARGETS:
            measured=[r for r in subset if r[f"observed_{t}"] is not None]
            summary["scores"][role][t]=metrics([r[f"observed_{t}"] for r in measured],
                                              [r[f"predicted_{t}"] for r in measured]) if measured else None
    write_json(out/"validation.json",summary)
    return summary


def paired_comparison(results, reference_allocation, out):
    """Subtract a matched baseline within each gas/phase/mesh/reference group.

    This avoids treating absolute clocking spread as numerical random noise.
    A few frozen phases do NOT constitute transient validation.
    """
    observations = read_observations(results)
    groups = {}
    group_fields=("context_id","phase_id","mesh_id","reference_id")
    for row in observations:
        key=tuple(row[k] for k in group_fields)
        group=groups.setdefault(key,{})
        if row["allocation_id"] in group:
            raise ValueError("Multiple runs for one allocation/group; resolve replicates explicitly.")
        group[row["allocation_id"]]=row
    contrasts=[]
    for key,group in groups.items():
        if reference_allocation not in group:
            raise ValueError(f"Missing matched baseline {reference_allocation} for {key}.")
        reference=group[reference_allocation]
        for name,row in group.items():
            if name==reference_allocation:
                continue
            item={**dict(zip(group_fields,key)),"allocation_id":name,"baseline_allocation":reference_allocation}
            for t in (*TARGETS,"sealing_effectiveness","ingress_kg_s"):
                item[f"difference_{t}"]=row[t]-reference[t] if row[t] is not None and reference[t] is not None else None
            contrasts.append(item)
    if not contrasts:
        raise ValueError("No candidate/baseline pairs found.")
    out=prepare_output(out)
    write_csv(out/"paired_differences.csv",contrasts)
    write_json(out/"paired_summary.json",{"n_pairs":len(contrasts),"warning":
               "Matched deterministic contrasts, not confidence intervals; no assertion of URANS or experimental validation."})
    return {"n_pairs":len(contrasts)}
