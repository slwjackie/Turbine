# 실행 검증 기록

2026-10-05, Linux CPU, Python 3.13.5. 테스트 환경의 정확한 버전은 첨부 실행 manifest를 참조합니다.
새 `rimseal_iwm` 모듈만 실행 검증했습니다. 기존 `turbine_poc`는 변경하지 않았으며 여기서 재실행하지 않았습니다.

## 소프트웨어 확인

`PYTHONPATH=src python -m pytest -q`: 39 passed (최종 전체 검사).
- 독립적인 각도 방향 수치적분과 ingress/egress 해석식 비교.
- 질량보존, 극한, 유일한 역함수 조건, 무차원화, tracer covariance, radial interpolation.
- noiseless 진폭·k 회복 및 null 경우 k=0.
- 검증 epsilon을 바꿔도 해당 held-out 예측/계수가 변하지 않는 누출 회귀시험.
- baseline-only로는 swirler correction을 식별할 수 없다는 실패 조건.
- 출처/반경/밀도비 등 입력 오류와 MC 재현성, 기존 출력 보호.

## 고정 seed 20261005의 합성자료 실행

5개 CF × 2개 형상 × 17개 purge =170행. 세 시나리오 모두 실제 연구실 측정자료가 아닙니다.
아래 값은 5개 전체 조건 holdout의 평균 epsilon RMSE입니다. 무차원 epsilon fraction 단위입니다.

| 합성자료 생성 가정 | annulus M0 | swirl M1 | 해석 |
|---|---:|---:|---|
| 후보 swirl 보정식이 참인 가상 세계 | 0.02593047 | 0.00559046 | 생성식을 회복하는 소프트웨어 확인 |
| swirl 추가효과가 없는 null | 0.00558726 | 0.00561980 | 추가 보정이 자동으로 유리하지 않음 |
| purge에 따라 진폭이 변함 | 0.02879766 | 0.01151910 | 10개 곡선 중 6개에서 reduced chi2>3 경고 |

swirl 시나리오만100 MC draws를 실행했고 1,000개 model-fold refit 중 실패0개였습니다.
이 구간은 조건부 measurement-perturbation 민감도이며 실제 모델의 예측구간 coverage 검증이 아닙니다.

## 공개 수치 확인

Choi Table4의 `(Phi0,epsilon)=(.074,.811)`과 `(.074,.908)`을 Tang 식에 역대입하면
유효 진폭은 각각 약0.152022,0.115048입니다. 같은 점을 되찾는 항등적 계산입니다.
다른 반경인 Fig9의0.151/0.117과 가깝다는 사실만으로 독립 검증/동일 물리량을 주장하지 않습니다.

원문 출처와 정확한 식·가정은 DATA_AND_MODEL.md에 기재했습니다.
