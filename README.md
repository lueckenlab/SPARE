# SPARE Benchmark for Sample Representation from Single-Cell Data

As single-cell datasets are growing, it is becoming possible to analyse differences between groups of samples on a cellular and molecular level. The promise of patient stratification, disease classification, and early-stage diagnosis has led to the development of several so-called sample representation methods. However, consistent standards for the evaluation of sample representation methods are lacking. We developed SPARE – a modular and extendable sample representation benchmark, defining 3 application-inspired metrics, and used these to compare 8 sample representation methods on 5 datasets, testing different preprocessing regimes. We find that the density-based method Gloscope outperforms other methods on most datasets and identify general best-practice preprocessing strategies for sample representation methods. We envision that this study will set standards for the development of sample representation methods and facilitate users in selecting an optimal tool, leading to improved outcomes for single-cell applications in precision medicine. 

![Benchmark overview](./figures/1_benchmark_overview.png)

This is an overview of the current pipeline:

```mermaid
graph TD
    A[download_data] --> D1
    A --> D2
    A --> D3
    A --> D4
    A --> D5
    A --> D6
    D1[(Synthetic data)] --> B[clean_data]
    D2[(COPD)] --> B
    D3[(COMBAT)] --> B
    D4[(Stephenson)] --> B
    D5[(HLCA)] --> B
    D6[(onek1k)] --> B
    B --> C[preprocess]

    C --> D[represent]

    D --> R1[[composition]]
    D --> R2[[pseudobulk]]
    D --> R3[[grouped_pseudobulk]]
    D --> R4[[random]]
    D --> R5[[scPoli]]
    D --> R6[[gloscope]]
    D --> R7[[MOFA]]
    D --> R8[[MrVI]]
    D --> R9[[PILOT]]
    R1 --> E[aggregate_representations]
    R2 --> E
    R3 --> E
    R4 --> E
    R5 --> E
    R6 --> E
    R7 --> E
    R8 --> E
    R9 --> E
    E --> G[evaluate]
    G --> M1([signal retention])
    G --> M2([batch removal])
    G --> M3([trajectory preservation])
    G --> M4([replicate robustness])
    G --> M5([scalability])
    
```

The development is still ongoing and this repo is subject to change.