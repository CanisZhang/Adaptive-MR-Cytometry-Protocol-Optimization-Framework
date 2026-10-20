# IMPULSED Protocol Optimization

```text
journal_submission/
|-- run_optimization.m       Search and independent validation
|-- optimization_config.m    Tissue, acquisition, fitting, and search settings
|-- +impulsed/
|   |-- protocol_model.m     Pulse construction and hardware feasibility
|   |-- evaluate_protocol.m  Signal simulation and parameter fitting
|   |-- fit_error.m          Estimation error and fitting failures
|   |-- sample_tissues.m     Truncated tissue perturbations
|   `-- s2_protocols.m       Supplementary Table S2 protocols
|-- tests/
|   `-- test_framework.m     Numerical checks
|-- mati/
|   `-- +mati/              MATI source subset and original license
|-- .gitignore              Generated output exclusions
`-- README.md               Directory structure
```
