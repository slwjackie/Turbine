# 면담에서 보여줄 내용

## 30초 설명

“공개 논문을 바탕으로 IWM의 sealing-effectiveness 부분을 구현하고,
실험 곡선에서 유효 Phi_A를 추정하는 분석 코드를 준비했습니다.
annulus swirl만 사용하는 모델과 wheel-space reference swirl을 추가한 후보를
같은 학습자료로 비교하며, 한 flow coefficient의 baseline과 swirler 자료는 함께 빼고 예측합니다.
지금은 합성자료와 공개된 몇 가지 수치만 확인했기 때문에 실제 swirler 보정식의 타당성을 검증했다고 보지는 않습니다.
연구실 자료의 측정 위치와 불확실성을 먼저 확인하고 적용 범위를 조정하고 싶습니다.”

## 5분 시연 순서

1. `published-check`: 논문 Table4의 역산과 Fig9 threshold의 측정 반경 차이 설명.
2. `swirl`: 한 조건 전체를 제외한 곡선 예측과 k 표시. 합성자료의 생성 가정을 명시.
3. `null`: 정보가 없는 변수를 추가하면 항상 좋아지는 것이 아님을 설명.
4. `purge_dependence`: 상수 Phi_A가 깨지면 좋은 보정 계수만으로 해결되지 않는 잔차를 설명.
5. CSV 양식과 uncertainty 가정을 보여주고 교수님 의견을 구함.

## 자료를 받기 전에 확인할 질문

- 원시 tracer 농도와 보정값, 반복 횟수, 공통 calibration error가 보존되어 있는가?
- baseline과 swirler의 epsilon·swirl이 동일한 CF/회전수/purge에서 짝지어져 있는가?
- radial swirl은 어느 반경·축방향 위치이며, seal 경계의 swirl과 어떻게 대응하는가?
- 이미 확보된 여러 CF별 baseline·swirler purge sweep이 있는가?
- annulus swirl은 측정값인가, vane/deviation 가정으로 추정한 값인가?
- Phi0의 b/간극 및 beta의 Omega*r 기준을 다른 논문과 어떻게 통일할 것인가?
- epsilon95와 완전 밀봉 판정은 무엇으로 정의하며 센서 분해능은 충분한가?
- 연구실의 현재 연구와 겹치지 않는 분석 범위 및 공개 가능한 결과는 무엇인가?

## 주장하지 않을 것

전체 IWM 압력장 재현, 새로운 불안정성 이론, 보편적인 swirler 법칙, 수소터빈 검증,
CFD의 실험 대체, 계수의 인과적 해석, 단 두 점으로의 일반화 검증을 주장하지 않습니다.
검증조건의 cavity swirl이 필요한 현재 모델과 swirl까지 미리 예측하는 설계도구는 구분합니다.

면담 후 첫 단계는 모델 복잡화가 아니라 자료와 측정 정의를 확인하는 것입니다.
