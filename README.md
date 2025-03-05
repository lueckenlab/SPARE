# Benchmark for sample representation from single-cell data

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
