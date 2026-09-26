# 현재 purge-allocation 주제에 맞춘 필독 순서

일반 추천이 아니라 이 프로젝트의 질문을 이해하기 위한 목록이다. 논문 데이터와 리그는 서로 다르므로
스월러 실험의 ε 곡선을 K1 고압터빈의 정답으로 그대로 사용하지 않는다.

## K1을 포함한 면담 전 중심 논문

### 1. Choi, Kim, Song (2025) — 연구 형상·공개 데이터의 출처
Effects of Leakage Flows from Secondary Flow Passages on the Flow and Efficiency of a High-Pressure Turbine.
KSFM 28(2), 7–13. DOI: https://doi.org/10.5293/kfma.2025.28.2.007
원문: https://ksfmjournal.org/xml/44747/44747.pdf

Fig.1에서 ①~⑥ 위치, §2의 계산조건, Table 2의 coded variables/단위, Fig.4의 상대위치 효과,
Fig.9~11에서 허브 누설과 부하·역류의 해석을 읽는다.
답할 질문: ①②가 고정인데도 여섯 유로 총량에 포함되는 이유는? ③④ 자체 기준유량을 왜 받아야 하는가?

### 2. Hong, Choi, Song (2025) — 냉각유량과 손실 메커니즘
Aerodynamic impact of the secondary air mass flow rate and swirl in a turbine blade cascade with tip clearance.
JMST 39, 4017–4026. DOI: https://doi.org/10.1007/s12206-025-0626-2

§2–3의 실험/CFD, §4.3–4.4와 결론을 읽는다. MFR/SR 증가가 팁 누설과 통로와류 손실을
동일한 방향으로 움직이지 않는 점을 이해한다. K1과 다른 저속 선형 캐스케이드임을 기억한다.
답할 질문: 유량을 늘렸는데 한 손실은 감소하고 다른 손실은 증가할 수 있는가?

### 3. Mingeun Choi, Goo, Cho, Song (2024) — 밀봉과 스월
Swirl Enhancement Effect on Turbine Rim Seal Performance.
Journal of Turbomachinery 146(5), 051002. DOI: https://doi.org/10.1115/1.4064230
공식 초록: https://snu.elsevierpure.com/en/publications/swirl-enhancement-effect-on-turbine-rim-seal-performance-2

실험장치, sealing-effectiveness 계측 정의, 주유동/휠스페이스 압력과 스월, 논문이 제시한 밀봉 해석을 읽는다.
로터 측 스월러 형상을 이번 프로젝트에 새로 도입하라는 뜻이 아니다.
답할 질문: ε는 어디에서 어떤 기준으로 측정했는가? 압력만으로 밀봉 변화를 설명할 수 있는가?

### 4. Kang et al. (2025) — 최근 연구실의 밀봉 설계 방향
Impact of the Stator-Side Swirler on Turbine Rim Sealing Effectiveness.
ASME Turbo Expo, GT2025-151746. DOI: https://doi.org/10.1115/GT2025-151746
공식 초록: https://snu.elsevierpure.com/en/publications/impact-of-the-stator-side-swirler-on-turbine-rim-sealing-effectiv/

2024 rotor-side와 2025 stator-side의 차이, 밀봉유량 절감과 동력 비용의 관계를 중심으로 읽는다.
공식 초록은 확인됐지만 전체 본문 접근이 필요하다. 초록만으로 계수·불확실성·검증 조건을 확정하지 않는다.
고정 총유량 재배분에 논문의 엔진출력 환산계수를 그대로 적용하지 않는다.

## 선택적으로 추가할 송 교수님 논문

### Lee, Yoon, Song (2026)
Effect of the outer stator swirler on wheel space sealing performance.
JMST 40, 5133–5139. DOI: https://doi.org/10.1007/s12206-026-0622-1

outer/inner 반경과 스월러 위치에 따른 밀봉 변화. 면담 전에는 초록/서론/결론을 우선 확인한다.
현재 배분 연구에서 스월러 형상 최적화까지 범위를 늘릴 필요는 없다.

### Hong, Choi, Kim, Song (2025)
블레이드 하류에서의 케이싱 측 냉각 유동이 선형 터빈 캐스케이드 유동장 및 손실에 미치는 영향.
대한기계학회논문집 B 49(10), 571–578. DOI: https://doi.org/10.3795/KSME-B.2025.49.10.571

⑤⑥ 또는 후속 익렬 영향으로 확장할 때 읽는다. 후단 주입은 현재 동익 효율 민감도가 작아도
하류 유동장까지 영향이 없다는 뜻이 아님을 이해한다. 원논문은 후속 정익의 실제 성능을 측정한 논문이 아니다.

## 이전 외부 문헌 목록에서 남길 것

- Zerobin et al. (2018), Aerodynamic Performance of Turbine Center Frames With Purge Flows—Part II:
  The Influence of Individual Hub and Tip Purge Flows. DOI https://doi.org/10.1115/1.4039363
  위치별 purge 영향이 다름을 읽되, K1 효율 목적함수와 TCF 손실을 구분한다.
- Wang et al. (2026), Multi-Objective Optimization of Multi-Channel Cooling Flow Distribution for
  Turbine Vanes Under Constant Total Cooling Air Flow. https://www.mdpi.com/2079-6412/16/8/985
  총량 제약·대리모델·독립 CFD 확인 방법을 참고한다. 내부 냉각 온도 최적화이며 purge ingress 연구가 아니다.
  제목과 달리 선택된 최적안 총유량은 기준 대비 약 1.13% 작다. 엄밀한 등식 제약 구현의 표준으로 그대로 복사하지 않는다.
- Lee et al. (2022), Controlling the flow distribution characteristics in the secondary air system of a gas turbine.
  DOI https://doi.org/10.1016/j.ijheatfluidflow.2022.108973
  유량을 지정하는 CFD와 실제 공급 네트워크/출구 면적 설계를 구분하기 위해 읽는다.
- Tang et al. (2025), Ingress Wave Model with Purge-Mainstream Density Ratio.
  DOI https://doi.org/10.1016/j.ijheatmasstransfer.2024.126372
  수소·밀도비 확장에 착수하기 직전에 읽는다. 무보정 IWM을 K1에 바로 적용하지 않는다.

## 면담 전 시간이 짧다면
K1 전체 → Hong 2025 메커니즘 → Choi 2024 밀봉 정의 → Kang 2025 최근 방향.
외부 논문은 Zerobin의 물리 질문과 Wang의 방법/한계부터 읽고, Tang은 수소 확장 전으로 미룬다.
원문을 못 읽은 경우 면담에서 초록/공개자료를 검토했다고 정확히 말한다.
