# IWM × wheel-space swirl — 면담 전 분석 코드

**실험 전 준비 코드입니다. 새로운 swirler 법칙이나 실험 검증 결과가 아닙니다.**
이 폴더는 상위 저장소의 기존 `turbine_poc` 이차공기 배분 연구와 독립된 Python 프로젝트입니다.
기존 코드는 수정하지 않습니다. GPU, CFD solver, Cantera가 필요하지 않습니다.

## 설치와 3분 실행

Python 3.10 이상. 저장소를 처음 받는 경우:

```bash
git clone --branch feat/iwm-swirl-premeeting https://github.com/slwjackie/Turbine.git
cd Turbine/rimseal
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
python -m pytest -q
python -m rimseal_iwm published-check --out outputs/published
python -m rimseal_iwm demo --scenario swirl --out outputs/swirl --draws 100
```

Windows에서는 가상환경 활성화 명령만 `.venv\Scripts\activate`로 바꿉니다.
재실행은 `outputs/swirl_v2`처럼 새 폴더를 지정합니다. 기존 결과는 덮어쓰지 않습니다.

반드시 확인할 대조 실험:

```bash
# cavity swirl의 추가 설명력이 없는 가상 세계
python -m rimseal_iwm demo --scenario null --out outputs/null
# 한 곡선 안에서 Phi_A가 purge에 따라 변하는 가상 세계: 상수 진폭 모델의 실패 진단
python -m rimseal_iwm demo --scenario purge_dependence --out outputs/purge_dependence
# baseline만 학습해서 swirler 보정계수를 알아낼 수 없다는 식별성 확인
python -m rimseal_iwm demo --split configuration --out outputs/configuration_holdout --no-plots
```

## 무엇이 구현되어 있나

| 기능 | 범위 |
|---|---|
| IWM 핵심식 | Tang Eqs.18,30–35: 무차원 유량, ingress/egress, effectiveness, 역산 및 목표유량 |
| 곡선별 보정 | 오차가 알려진 sealing-effectiveness 곡선에서 상수 Phi_A를 WLS로 추정 |
| M0 | annulus swirl에 대한 선형 진폭: 계수는 해당 rig의 학습자료로 다시 추정 |
| M1 | 기준 wheel-space swirl 편차에 따른 경험적 보정계수 k 추가 |
| 독립 검증 | 동일 flow coefficient의 baseline·swirler 곡선을 함께 제외 |
| 불확실성 | 측정값 perturbation → 학습 전 과정 재적합 → 조건부 민감도 구간 |
| 입력 방어 | 측정 반경·기준 purge·무차원화·밀도비·출처 및 데이터 누출 검사 |
| 공개 수치 확인 | Choi Table4의 점별 역산과 Fig9 임계유량을 서로 다른 반경으로 구분 |

**M1의 좋은 합성자료 성능은 의도한 수식을 코드가 회복하는지 확인하는 것일 뿐입니다.**
부호가 반대인 k도 허용하며, null 시나리오에서는 보정을 추가해도 좋아지지 않을 수 있습니다.
실험자료 없이 최적 k, 실제 절감률, 논문의 novelty를 확정하지 않습니다.

## 생성되는 파일

`report.md`부터 읽으세요. `heldout_effectiveness.png`, `heldout_target_flow.png`는 면담 시연용입니다.
`predictions.csv`와 `thresholds.csv`에는 held-out 예측, 외삽 여부, MC 유효 횟수가 기록됩니다.
`coefficients.csv`에는 fold별 k와 annulus/baseline-swirl 회귀계수가 기록됩니다.
`fold_membership.csv`는 어떤 곡선이 학습·검증에 포함됐는지 보여줍니다.
`curve_fits_DIAGNOSTIC_ONLY.csv`는 전체 자료의 적합 진단일 뿐 독립 검증 성능이 아닙니다.
`manifest.json`에는 데이터 SHA256, seed, 패키지 버전 및 불확실성에서 제외한 요소가 기록됩니다.

## 연구실 자료를 받은 뒤

```bash
python -m rimseal_iwm template --out private_data/measurements.csv
# docs/DATA_AND_MODEL.md에 맞춰 승인받은 실험자료로 빈 CSV를 채웁니다.
python -m rimseal_iwm analyze --data private_data/measurements.csv --out outputs/lab_v1 --draws 500
```

모든 `sigma_*`는 **절대값 기준 1 표준불확실성**입니다. 논문의 95%/퍼센트 수치를 그대로 넣지 마세요.
`beta_c_ref`는 모든 조건에서 미리 정한 동일 purge와 반경의 독립 측정값입니다.
검증조건의 effectiveness는 보정에 쓰지 않지만, 해당 조건의 reference swirl은 입력으로 사용합니다.
따라서 새로운 조건에서 아무것도 측정하지 않는 설계 예측이 아닙니다.

## 중요한 한계

현재 DR=1만 허용합니다. 압력장, 비정상 구조의 주파수, 밀도비 scaling, swirler 형상 최적화,
수소 연소 또는 실제 엔진 효율 계산은 구현하지 않았습니다.
상수 Phi_A가 곡선에 맞지 않으면 k를 늘려 덮기보다 그 가정부터 검토해야 합니다.
`Phi95`는 effectiveness 0.95의 모델상 목표유량이며 완전 밀봉 임계유량과 다릅니다.
모든 불확실성을 포함한 통계적 유의성 검증도 아직 아닙니다.

자세한 식과 입력 정의: [DATA_AND_MODEL.md](docs/DATA_AND_MODEL.md)
면담에서 설명할 내용: [INTERVIEW.md](docs/INTERVIEW.md)
실행 검증 기록: [VALIDATION.md](docs/VALIDATION.md)

연구실 원자료, 형상, 논문 PDF 전문은 공개 저장소에 넣지 마세요.
`.gitignore`는 보조 장치입니다. `outputs`의 입력 스냅샷에도 비공개 실험자료가 포함될 수 있습니다.
