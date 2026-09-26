# Turbine — 고압터빈 이차공기 배분 연구의 재현 가능한 출발점

**공개 논문 데이터 재현 + 고정 총유량 배분 탐색 + 후속 CFD 검증 인터페이스.**
이 저장소는 CFD solver도, 밀봉 성능 예측기도 아닙니다. 확인되지 않은 ε/수소 효과를 만들어 내지 않습니다.

## 지금 가능한 것

1. K1 논문의 Table 2(25개 case)를 원래 단위와 부호로 읽습니다.
2. 효율 변화/동익 한 개의 토크 변화에 각각 15계수 2차 RSM을 맞춥니다.
3. 학습 재현 오차와 leave-one-out의 식별 가능성을 구분합니다.
4. 실제 기준유량이 주어지면 ③+④ 총량을 고정한 배분 단면을 탐색합니다.
5. 새 CFD에서 채워 넣을 결과 양식을 만들고, 이후 독립 관측값 검증/동일 위상 paired comparison에 재사용합니다.

## 설치 및 첫 실행

Python 3.10 이상. CFX, GPU, Cantera는 이 단계에 필요하지 않습니다.

```bash
git clone https://github.com/slwjackie/Turbine.git
cd Turbine
python3 -m venv .venv
source .venv/bin/activate
# Windows: .venv\Scripts\activate
python -m pip install -e '.[dev]'
python -m pytest -q
python -m turbine_poc demo --baseline configs/illustrative.json --out outputs/day1
```

이미 존재하는 결과를 덮어쓰지 않습니다. 재실행할 때는 `outputs/day1b`처럼 새 경로를 주세요.
그래프가 필요 없으면 `--no-plots`를 붙입니다.

**주의: illustrative.json의 0.4/0.6 kg/s는 설명용 가정값입니다. 논문/연구실의 실제 유량이 아닙니다.**
따라서 demo가 제시하는 배분비는 실제 터빈 권장 설계가 아닙니다.

## 생성되는 결과

| 파일 | 의미 |
|---|---|
| `reproduction/report.md` | 한국어 재현 보고서 및 해석 주의사항 |
| `reproduction/metrics.json` | 학습 오차, 데이터 해시, 소프트웨어 버전, LOO 유효 fold 수 |
| `reproduction/coefficients.csv` | 두 RSM의 15개 계수 |
| `reproduction/endpoint_contrasts.csv` | ③~⑥의 서로 짝지어진 -1→+1 변화량 |
| `reproduction/one_factor_profiles.csv` | 곡률까지 확인할 수 있는 일변수 예측 곡선 |
| `allocation/allocation_summary.json` | 사용한 기준유량의 출처, 허용 범위, 미검증 RSM 후보 |
| `allocation/allocation_sweep.csv` | alpha별 m3, m4, coded variables와 RSM 예측 |
| `allocation/suggested_cases.csv` | 후속 CFD 후보. `predicted_*`는 관측값이 아님 |
| `allocation/cfd_results_template.csv` | 실제 CFD 결과를 기록할 빈 관측 칸 |

각 그래프는 독립적인 PNG로 생성됩니다.

## 단위·총량·제약을 구분하세요

- `delta_eta_pp = eta_with_leakage - eta_no_leakage`의 **percentage point** 값입니다.
  예: -0.3091 %p = 효율 fraction 차이 -0.003091. 덜 음수일수록 좋으므로 **최대화**합니다.
- `delta_torque_nm`는 **동익 한 개**의 토크 변화입니다. 엔진 총출력이 아닙니다.
- `x_i = 2*(m_i/m_i0 - 1)`. -1/0/+1은 각 유량의 자체 기준 대비 50/100/150%입니다.
- `M34=m3+m4`, `alpha=m3/M34`; ①②⑤⑥은 고정입니다.
- 총량 고정은 `m3_0*x3 + m4_0*x4 = 0`이지 일반적으로 `x3+x4=0`이 아닙니다.
- 여섯 유로 총량은 `m1+...+m6`. ①②⑤⑥을 고정하고 M34도 고정하면 여섯 유로 총량도 고정됩니다.
  여섯 값이 모두 주어지지 않으면 실제 총량은 `null`로 남깁니다. 엔진 전체 bleed budget과도 구분합니다.
- 50~150%의 각 유로 범위와 K1 Box–Behnken convex hull 밖 예측을 막습니다.
  표본 영역 안이라고 모델 정확도가 보증되는 것은 아닙니다.
- ε 제약, 열적 안전, 실제 공급압력/네트워크 실행 가능성은 **미평가**입니다.

## 재현 결과와 숨기지 않은 한계

공개된 반올림 수치에 대한 학습 R² ≈ 0.999917, 효율 RMSE ≈ 0.000647662 %p.
이는 새로운 CFD의 예측 오차나 격자 불확실성이 아닙니다.

**Case 00을 제외하는 LOO fold는 정의되지 않습니다.** 나머지 24개 점에서는
`x3²+x4²+x5²+x6²=2`이므로 상수항과 제곱항 합을 구분할 수 없고 설계행렬 rank가 15→14가 됩니다.
프로그램은 이를 감추기 위해 최소노름 해를 임의로 선택하지 않습니다.
`valid_folds_only`는 나머지 24개 fold의 참고값이며, 전체 LOO 점수나 독립 검증으로 보고하지 않습니다.

⑤⑥의 endpoint effect가 작다고 해서 모든 조건에서 중요하지 않다고 단정하지 않습니다.
특히 끝점 차이는 ④의 U자 곡률을 요약하지 못하고, 각 유로의 실제 kg/s 변화량도 다릅니다.

## 실제 연구자료를 받은 다음

```bash
cp configs/laboratory.template.json configs/lab.local.json
# ③④의 실제 기준유량과 출처를 채우세요. 다른 유량도 알려진 경우만 입력하세요.
python -m turbine_poc sweep --baseline configs/lab.local.json --out outputs/lab_slice
```

lab.local.json에 넣는 유량은 전주(360°) 환산 또는 공통 계산영역 기준 등 **서로 같은 유량 기준**이어야 합니다.
서로 다른 pitch/sector의 유량을 그대로 더하지 마세요. `lab_provided` 표시는 사용자가 출처를 선언한 것이며 자동 검증이 아닙니다.

`--budget-scale 0.9` 등으로 총량 변화 단면도 계산할 수 있습니다. 이 경우는 기준 총량 고정 연구가 아니라
총량 변화 연구로 명시됩니다. 모든 예측은 여전히 원논문 범위 내로 제한됩니다.

실제 CFD 후 관측값과 메타데이터를 채운 결과에 대해서만:

```bash
python -m turbine_poc validate-cfd --results private_data/new_cfd.csv --out outputs/validation
python -m turbine_poc paired --results private_data/phase_cases.csv --reference-allocation A00 --out outputs/paired
```

`validate-cfd`는 원논문과 일치하는 물리조건·위상·격자·기준 정의에서만 허용합니다.
수소, 새 격자, 다른 상대위치 데이터를 기존 학습집합에 몰래 섞지 않습니다.
`paired`는 같은 조건·위상·격자에서 후보와 기준안의 **차이**를 계산합니다.
절대 위상 산포를 무작위 수치잡음으로 간주하지 않으며, frozen-rotor 비교를 비정상 검증으로 부르지 않습니다.

## 하루 구현 범위와 이후 연결

오늘: PoC 1(RSM), PoC 2(배분), PoC 5(보고서/개념도), 실제 CFD 결과 저장·검증 양식.
후속: 기준 CFD 재현 → cavity/평가면 점검 → tracer와 ingress 계산 → 밀봉 제약 → 상대위치 검증.
Cantera/IWM 및 수소 확장은 지금 구현하지 않았습니다. 필요한 검증자료가 없는데 ε를 만드는 것보다 기반을 먼저 확보합니다.

- [필독 논문과 읽을 질문](docs/reading_guide.md)
- [면담용 개념과 단계별 연결](docs/research_handoff.md)

## 데이터 출처·공개 범위

Choi, M.; Kim, M.; Song, S. J. (2025), *Effects of Leakage Flows from Secondary Flow Passages on the Flow and Efficiency of a High-Pressure Turbine*, KSFM 28(2), 7–13.
DOI: https://doi.org/10.5293/kfma.2025.28.2.007
Table 2 p.11: https://ksfmjournal.org/xml/44747/44747.pdf

원문 그림/전문은 재배포하지 않습니다. 공개 수치 표의 전사와 출처만 포함합니다.
테스트의 합성 데이터는 소프트웨어 검증 전용이며 과학적 결과에 포함하지 않습니다.
연구실의 비공개 형상·CFX 파일·허가되지 않은 결과는 공개 저장소에 올리지 마세요.
`.gitignore`는 보조 장치일 뿐입니다. 업로드 전 자료 공개 허가를 직접 확인해야 합니다.
