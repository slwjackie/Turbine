# 식·가정·데이터 계약

## 1. 원논문에서 구현한 부분

Tang, Scobie, Wang & Lock (2024), *A theoretical model for ingress through turbine rim seals based on physically-observed unsteadiness*, IJHFF 106,109300. DOI: https://doi.org/10.1016/j.ijheatfluidflow.2024.109300

Eq.18: `Phi0 = mdot / (2*pi*rho0*Omega*r1^2*sc)`.
Eq.30: `xi = asin(min(Phi0/PhiA,1))`.
Eqs.31–33: `Phii = PhiA/pi*cos(xi)+Phi0/pi*(xi-pi/2)`;
`Phie = PhiA/pi*cos(xi)+Phi0/pi*(xi+pi/2)`; `epsilon=Phi0/Phie`.
코드는 Phi0=0에서 나눗셈 특이점을 피하고 질량보존 `Phie-Phii=Phi0`를 검사합니다.

**PhiA가 purge sweep 전체에서 일정할 때만** Eq.35의 `Phi_min=PhiA`를 사용합니다.
동일한 가정에서 `Phi95=0.7448170868845*PhiA`입니다. 이는 이 식의 수치적 성질이지 새 실험 결과가 아닙니다.
`amplitude_from_point`는 Phi0>0, 0<epsilon<1에서만 유일한 역산을 허용합니다.
압력 모델과 밀도 차이의 부록식은 이번 구현에 포함하지 않았습니다.

## 2. 기존 가설과 새 후보식을 구분

Vella, Salvatori, Scobie, Lock, Sangan & Tang (2026), *Scaling Sealing Performance Across Engine Operating Conditions*, J. Eng. Gas Turbines Power148(6),061002. DOI: https://doi.org/10.1115/1.4070053

이 논문의 annulus-swirl 선형 가설을 M0 `PhiA=a+b*beta_a`로 구현했습니다.
다른 시험장치의 수치 계수를 SNU에 옮기지 않습니다. M0와 M1은 동일한 학습 곡선을 사용합니다.
논문의 low-TRL 검증을 실제 수소터빈 엔진 검증으로 해석하지 않습니다.

M1은 본 코드의 **검증 대기 중인 경험식**입니다:

```
beta_c_baseline(beta_a) = d0 + d1*beta_a
Delta_beta_c = beta_c_ref - beta_c_baseline(beta_a)
PhiA = (a+b*beta_a)*exp(k*Delta_beta_c)
```

기준 추세 d0,d1은 학습집합의 baseline reference-swirl만으로 추정합니다.
진폭 응답에 대한 계수는 a,b에서 a,b,k로 하나 늘지만, **보조 회귀계수 d0,d1도 데이터로 추정**합니다.
전체 자유도가 단지 하나 늘었다고 해석하거나 이 기준선 추정의 비용을 숨기지 않습니다.
학습 범위 양 끝에서 양수 진폭을 로그 변수로 매개화해 M0의 선형성을 유지합니다.
k의 부호는 강제하지 않습니다. 범위 밖에서는 외삽을 표시하고 선형 진폭이 음수가 되면 예측을 중단합니다.

`beta_c_ref`는 측정 위치의 설명변수입니다. 서로 다른 반경의 beta 차이를 실제 seal 입구의
속도차·전단응력과 동일시하지 않습니다. swirl과 형상 효과가 함께 변하므로 인과효과를 분리했다는 주장도 불가합니다.

## 3. 입력 열

한 행은 한 purge setting의 한 관측입니다. 한 `curve_id`는 하나의 조건·형상에서 얻은 purge sweep입니다.
반복 측정 행은 가능하나 상관성 구조는 사전에 검토하세요.

| 열 | 정의 |
|---|---|
| curve_id | purge 곡선 식별자 |
| condition_id | 동일 flow coefficient의 baseline·swirler를 묶는 식별자 |
| configuration | `baseline` 또는 `rotor_swirler` |
| flow_coefficient | 동일 기준의 nominal flow coefficient; 조건 묶음의 공통 설정값 |
| phi0, epsilon | 일관된 정의의 purge와 측정 effectiveness |
| sigma_epsilon, sigma_phi0 | 절대 1-sigma 표준불확실성; epsilon에는 양수 필요 |
| beta_a, sigma_beta_a | annulus swirl 및 절대 1-sigma; 같은 조건의 두 형상에 공유된 값 |
| beta_c_ref, sigma_beta_c_ref | 사전 정의 reference purge에서의 cavity swirl 및 절대 1-sigma |
| phi_ref | reference swirl을 측정한 purge; 모든 곡선에서 동일 |
| r_epsilon_over_b, r_swirl_over_b | 각 측정의 반경; 서로 같을 필요는 없지만 정의는 모든 곡선에서 고정 |
| density_ratio | rho_purge/rho_annulus; 현재 1만 허용 |
| rig_id, seal_id, normalization | rig·seal·Phi 및 beta 기준 정의를 식별하는 문자열 |
| source_kind | synthetic / published_digitized / lab_measured |
| source | 자료 출처, 파일/측정 run/처리 버전 등 추적 정보 |

프로그램이 검사하는 것은 선언값의 일관성입니다. normalization 문자열이 같다고 실제 단위를 자동 검증하는 것은 아닙니다.
회전수·형상·센서 위치가 달라지는 캠페인은 이 단일-CF 검증 범위를 넘으므로 먼저 나누어 분석합니다.
동일 CF를 임의로 다른 condition_id로 나누는 것은 금지합니다.
최소 3개 이상의 정보가 있는 비포화 purge 수준이 있어야 곡선 진폭을 적합합니다.

`radial_reference`는 동일 purge의 방사방향 측정 profile만 선형 보간하며 외삽하지 않습니다.
이 함수에 넣기 전에 run·purge·센서 위치가 맞는지 확인해야 합니다. `phi_ref` 선택이나
유량방향 보간·추세모델은 자동화하지 않았습니다. 검증 오차가 최소가 되도록 phi_ref를 고르면 데이터 누출입니다.

## 4. Choi 논문과 무차원화·측정 위치

Choi, Goo, Cho & Song (2024), *Swirl Enhancement Effect on Turbine Rim Seal Performance*, J. Turbomachinery146(5),051002. DOI: https://doi.org/10.1115/1.4064230

Table4는 Phi0=0.074, r/b=0.923의 epsilon·beta 값입니다. Fig9 임계유량은 r/b=0.941입니다.
`published-check`는 이 사실을 각각의 CSV에 보존하며 서로 짝지어 k를 적합하지 않습니다.
annulus beta 약1.15는 직접 측정값이 아니라 가정·상관식을 사용한 추정값입니다.

Choi의 Phi0는 b 및 axial clearance 기준입니다. Tang의 r1, sc와 같은 값인지 확인해야 합니다.
`renormalize_phi`는 물리적 질량유량을 보존하는 변환만 수행합니다. 선택한 반경/간극의 물리적 적절성을 대신 판단하지 않습니다.
stator-surface tracer effectiveness를 IWM의 혼합 effectiveness와 대응시키는 것도 검증해야 할 가정입니다.

## 5. 불확실성과 검증

곡선별 WLS 목적함수는 `sum(((epsilon_model-epsilon_measured)/sigma_epsilon)^2)`입니다.
조건부 local `se_phi_a`는 Jacobian으로 계산하며 model discrepancy까지 포함하지 않습니다.
고정 purge에 대한 평균 진폭으로 적합하므로 purge-dependent amplitude를 발견하는 데는 잔차 구조도 필요합니다.
`reduced_chi2>3` 표시는 휴리스틱 경고이지 사전 유의수준이 설정된 검정이 아닙니다.

MC는 학습 epsilon/Phi0, 조건 공유 beta_a, 곡선 공유 beta_c_ref를 perturb하고 모든 학습단계를 다시 실행합니다.
예측 epsilon은 검증조건의 nominal Phi0에서 계산합니다. 검증 epsilon은 예측을 만드는 데 사용하지 않습니다.
구간은 가정한 측정오차에 대한 2.5/97.5 percentile **민감도 구간**으로, 교정된 신뢰구간/예측구간이 아닙니다.
최소20개 draw 및80% 성공이 안 되면 구간을 주지 않습니다. 실제 분석은 draw 수에 따른 분위수 안정성을 확인하세요.

현재 센서 공통 bias, 형상·반경/간극 오차, 전체 측정 공분산, 모델오차, 관측 Phi95 보간오차는 전파하지 않습니다.
관측 Phi95는 유일하게 bracket된 경우만 선형 보간하며, 이를 오차 없는 참값이라고 부르지 않습니다.
Choi Table3의95% 및% 표기를 절대1-sigma로 그대로 넣지 않습니다. 상대/절대, coverage factor와
반복 측정인지 공통 보정오차인지 확인한 뒤 환산합니다. 근거 없는 sigma는 채우지 마세요.

실제 예측 개선의 통계적 유의성을 주장하려면 paired threshold uncertainty와 shared bias/
model discrepancy 처리를 추가하고, 기존 조건 외의 독립 운전조건에서 확인해야 합니다.
