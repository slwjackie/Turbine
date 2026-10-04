# Rim-seal IWM 연구 준비

새 연구 코드는 [rimseal/README.md](rimseal/README.md)에 있습니다.
기존 `turbine_poc` 이차공기 배분 코드와 독립된 프로젝트이며 기존 파일은 수정하지 않습니다.

```bash
cd rimseal
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
python -m pytest -q
python -m rimseal_iwm demo --out outputs/interview --draws 100
```

합성자료 데모는 실험 검증이 아닙니다. 공개 논문 전체/PDF나 비공개 연구실 자료는 포함하지 않습니다.
