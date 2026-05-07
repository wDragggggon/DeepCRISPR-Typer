# 重复序列分类模型配置
REPEAT_MODEL_CONFIG = {
    "num_subtypes": 33,  # CRISPR-Cas亚型数量
    "max_repeat_length": 50,  # 最大重复序列长度
    "embedding_dim": 4,   # 独热编码维度
    "cnn_channels": 64,
    "transformer_heads": 8,
    "transformer_layers": 2,
    "dropout_rate": 0.1
}

# Cas蛋白分类模型配置（基于ESM-2）
CAS_PROTEIN_MODEL_CONFIG = {
    "esm_model_name": "esm2_t33_650M_UR50D",  # ESM-2模型名称
    "num_cas_families": 6,  # Cas蛋白功能家族数量
    "lora_r": 8,           # LoRA秩
    "lora_alpha": 16,      # LoRA缩放因子
    "lora_dropout": 0.1,
    "max_protein_length": 2000  # 最大蛋白序列长度
}

# HMM扫描配置
HMM_CONFIG = {
    "evalue_threshold": 1e-5,
    "bitscore_threshold": 30.0,
    "domain_evalue_threshold": 1e-3
}

# 规则矩阵配置
RULE_CONFIG = {
    "adaptation_score_threshold": 0.5,
    "interference_score_threshold": 0.5,
    "final_confidence_threshold": 0.7
}